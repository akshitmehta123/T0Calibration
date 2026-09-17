package org.clas.modules.analysis;

import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataBank;

import org.jlab.groot.data.H1F;
import org.jlab.groot.math.F1D;

import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.geom.base.ConstantProvider;
import org.jlab.detector.calib.utils.DatabaseConstantProvider;

import java.util.*;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;

/**
 * ============================================================================
 * AHDC T0 CALIBRATION
 * ============================================================================
 *   SECTION 1: GEOMETRY       - detector wire-line lookups
 *   SECTION 2: HIT EXTRACTION - turns HIPO banks into a simple list of hits
 *   SECTION 3: CLUSTERING     - pure algorithm: hits in -> clusters out
 *   SECTION 4: FITTING        - histograms + Gaussian fit -> T0 per wire
 *   SECTION 5: SUMMARY        - reduces per-wire T0s to mean/RMS/%fit per layer
 *   SECTION 6: RUNNER (main)  - layers 1-5 together, no canvas, no GUI, just
 *                              prints the three output values (mean/RMS/%fit).

 * Histograms (H1F) are used only as an internal tool to do the Gaussian fit - they are never drawn.
 * ============================================================================
 */
public class T0Calibration {

    // ---- quality cuts applied to every hit before it is used -----------------
    private static final int MIN_ADC = 500;
    private static final int MAX_PED = 300;
    private static final double MIN_TOT = 1.0;

    // ---- clustering thresholds -------------------------------------------
    private static final int MIN_HITS_PER_CLUSTER = 7;
    private static final int MIN_SUPERLAYERS = 4;
    private static final double SEED_TIME_WINDOW_NS = 150.0;
    private static final double CANDIDATE_TIME_WINDOW_NS = 150.0;
    private static final double CANDIDATE_PHI_WINDOW_DEG = 25.0;
    private static final double CANDIDATE_MAX_PERP_DIST = 20.0;

    // ---- fit search / quality windows --------------------------------------
    private static final double SEARCH_MIN_NS = 150.0;
    private static final double SEARCH_MAX_NS = 500.0;
    private static final double SEED_WINDOW_NS = 35.0;
    private static final double MIN_ENTRIES = 20;
    private static final double MAX_MEAN_SHIFT_FROM_PEAK = 80.0;
    private static final double MIN_SIGMA = 2.0;
    private static final double MAX_SIGMA = 150.0;
    private static final double BACKGROUND_TOLERANCE = -5.0;


    private static final double VERTEX_Z_OFFSET = 5.0; // vertex-Z correction used throughout

    // ---- input file -----------------------------------------------------
    // Edit this path to point at whichever hipo file you want to run over.
    private static final String HIPO_FILE = "/path/to/your/file.hipo";
            
    private static final int SECTOR = 1;


    private static final PrintStream REAL_STDOUT = System.out; // the real stdout, saved once so it can always be restored

    /**
     * The jlab-io / groot / minuit libraries print/log a lot of internal
     * diagnostic noise ("reader:: ...", "[fit-benchmark] ...", timestamped
     * Minuit "MnPosDef"/"MnHesse" warnings) that has nothing to do with our
     * own output. This turns all of that off up front.
     */
    private static void silenceLibraryLogging() {
        LogManager.getLogManager().reset();
        Logger.getLogger("").setLevel(Level.OFF);
    }

    /** Runs a library call while discarding anything it writes to System.out, then restores stdout. */
    private static void runQuietly(Runnable libraryCall) {
        System.setOut(new PrintStream(new OutputStream() {
            @Override public void write(int b) { /* discard */ }
        }));
        try {
            libraryCall.run();
        } finally {
            System.setOut(REAL_STDOUT);
        }
    }


    // ============================================================================
    // SECTION 1: GEOMETRY
    // All AlertDCDetector wire-line lookups live here. If a geometry call ever
    // needs to change, this is the only place to look.
    // ============================================================================

    private final AlertDCDetector detector;

    private T0Calibration() {
        ConstantProvider cp = new DatabaseConstantProvider(1, "default", "default");
        this.detector = new AlertDCFactory().createDetectorCLAS(cp);
    }

    /** Project a wire onto the plane Z = vz. Returns null if the wire doesn't cross that Z. */
    private Point3D projectWireToVertexZ(int sector, int layerCode, int component, double vz) {
        int superlayer = layerCode / 10;
        int layer = layerCode % 10;
        Line3D wireLine;
        try {
            wireLine = detector.getSector(sector).getSuperlayer(superlayer).getLayer(layer).getComponent(component).getLine();
        } catch (Exception e) {
            return null;
        }
        Point3D p1 = wireLine.origin();
        Point3D p2 = wireLine.end();
        double dz = p2.z() - p1.z();
        if (Math.abs(dz) < 1e-6) return null;

        double q = (vz - p1.z()) / dz;
        if (q < 0 || q > 1) return null;

        return new Point3D(p1.x() + q * (p2.x() - p1.x()), p1.y() + q * (p2.y() - p1.y()), vz);
    }

    private int numWiresInLayer(int sector, int layerCode) {
        try {
            int superlayer = layerCode / 10;
            int layer = layerCode % 10;
            return detector.getSector(sector).getSuperlayer(superlayer).getLayer(layer).getNumComponents();
        } catch (Exception e) {
            return 0;
        }
    }

    // ============================================================================
    // SECTION 2: HIT EXTRACTION
    // Turns one HIPO event into a simple list of hits, grouped by layer.
    // This is the ONLY part of the file that calls event.getBank(...) /
    // DataBank.get*(...) - if bank/column names ever change, fix it here.
    // ============================================================================

    /** One AHDC wire hit, already projected to the event vertex Z. */
    private static class Hit {
        final Point3D pos;      // (x, y, vz)
        final double time;      // leadingEdgeTime - startTime (ns)
        final int layer;        // superlayer*10 + layer
        final int wire;
        final int adc;
        final float ped;
        final float tot;
        final int trackid;      // -1 if no official track used this hit

        Hit(Point3D pos, double time, int layer, int wire, int adc, float ped, float tot, int trackid) {
            this.pos = pos; this.time = time; this.layer = layer; this.wire = wire;
            this.adc = adc; this.ped = ped; this.tot = tot; this.trackid = trackid;
        }

        String key() { return layer + "_" + wire; }

        String hitID() { return layer + "_" + wire + "_" + String.format("%.3f", time); }

        boolean passesQualityCuts() { return adc > MIN_ADC && ped < MAX_PED && tot > MIN_TOT; }
    }

    /** Event vertex Z (cm), offset-corrected. Returns 0 if REC::Particle is missing/empty. */
    private double eventVertexZ(DataEvent event) {
        double vz = 0.0;
        if (event.hasBank("REC::Particle")) {
            DataBank rec = event.getBank("REC::Particle");
            if (rec.rows() > 0) vz = rec.getFloat("vz", 0);
        }
        return vz - VERTEX_Z_OFFSET;
    }

    /** (layerCode_wire) -> trackid, from AHDC::hits. */
    private Map<String, Integer> trackIdMap(DataEvent event) {
        Map<String, Integer> trackMap = new HashMap<>();
        if (!event.hasBank("AHDC::hits")) return trackMap;

        DataBank hitsBank = event.getBank("AHDC::hits");
        for (int i = 0; i < hitsBank.rows(); i++) {
            int superlayer = hitsBank.getInt("superlayer", i);
            int layer = hitsBank.getInt("layer", i);
            int wire = hitsBank.getInt("wire", i);
            int trackid = hitsBank.getInt("trackid", i);
            trackMap.put((superlayer * 10 + layer) + "_" + wire, trackid);
        }
        return trackMap;
    }

    /** Extract all usable AHDC::adc hits for this event, projected to vertex Z, grouped by layer. */
    private Map<Integer, List<Hit>> extractHitsByLayer(DataEvent event, double vertexZ) {

        Map<Integer, List<Hit>> hitsByLayer = new HashMap<>();
        if (!event.hasBank("AHDC::adc")) return hitsByLayer;

        DataBank adcBank = event.getBank("AHDC::adc");
        Map<String, Integer> trackMap = trackIdMap(event);
        float startTime = 0f;
        if (event.hasBank("REC::Event")) {
            startTime = safeGetFloat(event.getBank("REC::Event"), "startTime", 0);
        }

        for (int i = 0; i < adcBank.rows(); i++) {

            short wfType = safeGetShort(adcBank, "wfType", i);
            if (wfType == 4 || wfType == 5) continue; // rejected waveform quality

            int sector = safeGetInt(adcBank, "sector", i);
            int layerCode = safeGetInt(adcBank, "layer", i);
            int component = safeGetInt(adcBank, "component", i);
            float leadTime = safeGetFloat(adcBank, "leadingEdgeTime", i);
            int adc = safeGetInt(adcBank, "ADC", i);
            float ped = safeGetFloat(adcBank, "ped", i);
            float tot = safeGetFloat(adcBank, "timeOverThreshold", i);

            Point3D projected = projectWireToVertexZ(sector, layerCode, component, vertexZ);
            if (projected == null) continue;

            int trackid = trackMap.getOrDefault(layerCode + "_" + component, -1);
            double hitTime = leadTime - startTime;

            hitsByLayer.computeIfAbsent(layerCode, k -> new ArrayList<>())
                    .add(new Hit(projected, hitTime, layerCode, component, adc, ped, tot, trackid));
        }
        return hitsByLayer;
    }

    private static int safeGetInt(DataBank bank, String col, int row) {
        try { return bank.getInt(col, row); } catch (Exception ex) { return 0; }
    }
    private static float safeGetFloat(DataBank bank, String col, int row) {
        try { return bank.getFloat(col, row); } catch (Exception ex) { return 0f; }
    }
    private static short safeGetShort(DataBank bank, String col, int row) {
        try { return bank.getShort(col, row); } catch (Exception ex) { return 0; }
    }


    // ============================================================================
    // SECTION 3: CLUSTERING
    // Pure algorithm. Input: this event's hits, grouped by layer.
    // Output: a list of clusters (each cluster = a List<Hit> that passed the
    // timing/phi/distance cuts and the minimum hit-count / superlayer-count
    // requirement). This method has no idea a HIPO file or a bank exists
    // ============================================================================

    private List<List<Hit>> findClusters(Map<Integer, List<Hit>> hitsByLayer) {

        List<List<Hit>> clusters = new ArrayList<>();
        List<Integer> layers = new ArrayList<>(hitsByLayer.keySet());
        Collections.sort(layers);

        Set<String> used = new HashSet<>();

        for (int i = 0; i < layers.size(); i++) {
            List<Hit> hits1 = hitsByLayer.get(layers.get(i));
            if (hits1 == null) continue;

            for (int j = i + 1; j < layers.size(); j++) {
                List<Hit> hits2 = hitsByLayer.get(layers.get(j));
                if (hits2 == null) continue;

                for (Hit h1 : hits1) {
                    for (Hit h2 : hits2) {

                        if (used.contains(h1.hitID()) || used.contains(h2.hitID())) continue;
                        if (Math.abs(h2.time - h1.time) > SEED_TIME_WINDOW_NS) continue;

                        Vector3D dir = new Vector3D(h2.pos.x() - h1.pos.x(), h2.pos.y() - h1.pos.y(), 0);
                        dir.unit();

                        List<Hit> cluster = new ArrayList<>();
                        cluster.add(h1);
                        cluster.add(h2);

                        Set<String> tempUsed = new HashSet<>();
                        tempUsed.add(h1.hitID());
                        tempUsed.add(h2.hitID());

                        for (int k = j + 1; k < layers.size(); k++) {
                            List<Hit> candidates = hitsByLayer.get(layers.get(k));
                            if (candidates == null) continue;

                            Hit best = findBestCandidate(h1, dir, candidates);
                            if (best != null && !used.contains(best.hitID())) {
                                cluster.add(best);
                                tempUsed.add(best.hitID());
                            }
                        }

                        if (cluster.size() >= MIN_HITS_PER_CLUSTER && countSuperlayers(cluster) >= MIN_SUPERLAYERS) {
                            clusters.add(cluster);
                            used.addAll(tempUsed);
                        }
                    }
                }
            }
        }
        return clusters;
    }

    /** Score every candidate against the seed hit + direction; return the single best match (or null). */
    private Hit findBestCandidate(Hit seed, Vector3D dir, List<Hit> candidates) {

        double bestScore = Double.MAX_VALUE;
        Hit best = null;

        for (Hit c : candidates) {

            double dt = c.time - seed.time;
            if (Math.abs(dt) > CANDIDATE_TIME_WINDOW_NS) continue;

            double phiRef = Math.atan2(seed.pos.y(), seed.pos.x());
            double phiCand = Math.atan2(c.pos.y(), c.pos.x());
            double dphi = Math.atan2(Math.sin(phiRef - phiCand), Math.cos(phiRef - phiCand));
            if (Math.abs(dphi) > Math.toRadians(CANDIDATE_PHI_WINDOW_DEG)) continue;

            Vector3D v = new Vector3D(c.pos.x() - seed.pos.x(), c.pos.y() - seed.pos.y(), 0);
            double perpDist = v.cross(dir).mag();
            if (perpDist > CANDIDATE_MAX_PERP_DIST) continue;

            double score = perpDist + Math.abs(dt) + Math.abs(dphi);
            if (score < bestScore) {
                bestScore = score;
                best = c;
            }
        }
        return best;
    }

    private int countSuperlayers(List<Hit> cluster) {
        Set<Integer> sl = new HashSet<>();
        for (Hit h : cluster) sl.add(h.layer / 10);
        return sl.size();
    }


    // ============================================================================
    // SECTION 4: FITTING
    // Fills a clustered-timing histogram per wire, then fits each one with a
    // Gaussian + flat background, seeded from the histogram's own peak, and
    // applies the accept/reject quality checks. Histograms are used only as
    // an internal fitting tool - they are never drawn to a canvas.
    // ============================================================================

    /** Result of fitting one wire's clustered-timing histogram. */
    private static class WireFitResult {
        final String key;
        final boolean accepted;
        final String rejectReason;
        final double t0;   // peak - 50 ns

        WireFitResult(String key, boolean accepted, String rejectReason, double t0) {
            this.key = key; this.accepted = accepted; this.rejectReason = rejectReason; this.t0 = t0;
        }
    }

    private WireFitResult fitWire(String key, H1F hCl) {

        if (hCl == null || hCl.getEntries() < MIN_ENTRIES) {
            return new WireFitResult(key, false, "insufficient entries", Double.NaN);
        }

        int firstBin = hCl.getAxis().getBin(SEARCH_MIN_NS);
        int lastBin = hCl.getAxis().getBin(SEARCH_MAX_NS);

        int maxBin = firstBin;
        double maxContent = -1.0;
        for (int b = firstBin; b <= lastBin; b++) {
            double c = hCl.getBinContent(b);
            if (c > maxContent) { maxContent = c; maxBin = b; }
        }
        if (maxContent <= 0) return new WireFitResult(key, false, "no peak found", Double.NaN);

        double peak = hCl.getAxis().getBinCenter(maxBin);

        // seed mean/sigma from a weighted window around the peak
        double weightedSum = 0, weightedXSum = 0, weightedX2Sum = 0;
        for (int b = firstBin; b <= lastBin; b++) {
            double x = hCl.getAxis().getBinCenter(b);
            double y = hCl.getBinContent(b);
            if (Math.abs(x - peak) > SEED_WINDOW_NS || y <= 0) continue;
            weightedSum += y;
            weightedXSum += y * x;
            weightedX2Sum += y * x * x;
        }
        double seedMean = peak;
        double seedSigma = 15.0;
        if (weightedSum > 0) {
            seedMean = weightedXSum / weightedSum;
            double variance = weightedX2Sum / weightedSum - seedMean * seedMean;
            if (variance > 0) seedSigma = Math.sqrt(variance);
        }
        seedSigma = Math.max(5.0, Math.min(50.0, seedSigma));

        double fitMin = Math.max(SEARCH_MIN_NS, peak - 100.0);
        double fitMax = Math.min(SEARCH_MAX_NS, peak + 50.0);

        double bgSum = 0; int bgCount = 0;
        for (int b = firstBin; b <= lastBin; b++) {
            double x = hCl.getAxis().getBinCenter(b);
            if ((x >= fitMin && x < fitMin + 5.0) || (x <= fitMax && x > fitMax - 5.0)) {
                bgSum += hCl.getBinContent(b);
                bgCount++;
            }
        }
        double background = bgCount > 0 ? bgSum / bgCount : 0.0;
        double amplitude = Math.max(1.0, maxContent - background);

        F1D gaus = new F1D("gaus_" + key, "[bg] + [amp]*gaus(x,[mean],[sigma])", fitMin, fitMax);
        gaus.setParameter(0, background);
        gaus.setParameter(1, amplitude);
        gaus.setParameter(2, seedMean);
        gaus.setParameter(3, seedSigma);

        try {
            runQuietly(() -> hCl.fit(gaus));
        } catch (Exception fitException) {
            return new WireFitResult(key, false, "fit failed: " + fitException.getMessage(), Double.NaN);
        }

        double fitBackground = gaus.getParameter(0);
        double fitAmplitude = gaus.getParameter(1);
        double mean = gaus.getParameter(2);
        double sigma = Math.abs(gaus.getParameter(3));
        double t0 = peak - 50.0;

        String rejectReason = null;
        if (Math.abs(mean - peak) > MAX_MEAN_SHIFT_FROM_PEAK) rejectReason = "mean moved too far from peak";
        else if (sigma < MIN_SIGMA || sigma > MAX_SIGMA) rejectReason = "sigma out of range: " + sigma;
        else if (fitAmplitude <= 0) rejectReason = "non-positive amplitude";
        else if (fitBackground < BACKGROUND_TOLERANCE) rejectReason = "negative background: " + fitBackground;

        return new WireFitResult(key, rejectReason == null, rejectReason, t0);
    }


    // ============================================================================
    // SECTION 5: SUMMARY
    // Reduces all the per-wire fit results down to the three output values: mean T0, RMS, and percent-of-wires-fit, PER LAYER.
    // ============================================================================

    private void printSummary(Map<Integer, List<Double>> t0ValuesPerLayer,
                              Map<Integer, Integer> acceptedPerLayer,
                              int sector, int totalWires, int acceptedWires, int rejectedWires) {

        System.out.println("==================================================");
        System.out.println("              T0 FIT SUMMARY");
        System.out.println("==================================================");
        System.out.println("Total wires tested = " + totalWires);
        System.out.println("Accepted fits      = " + acceptedWires);
        System.out.println("Rejected fits      = " + rejectedWires);
        System.out.println();
        System.out.println("Layer | MeanT0(ns) | RMS(ns) | %WiresFit");

        for (int layerCode : t0ValuesPerLayer.keySet()) {

            List<Double> values = t0ValuesPerLayer.get(layerCode);
            if (values.isEmpty()) continue;

            double sum = 0;
            for (double v : values) sum += v;
            double mean = sum / values.size();

            double sq = 0;
            for (double v : values) sq += (v - mean) * (v - mean);
            double rms = Math.sqrt(sq / values.size());

            int accepted = acceptedPerLayer.getOrDefault(layerCode, 0);
            int total = numWiresInLayer(sector, layerCode);
            double percentFit = total > 0 ? 100.0 * accepted / total : 0.0;

            System.out.printf("%5d | %10.3f | %7.3f | %9.1f%n", layerCode, mean, rms, percentFit);
        }
    }

    // ============================================================================
    // SECTION 6: RUNNER
    // Reads the hipo file, calls extraction -> clustering -> fitting -> summary.
    // ============================================================================

    private void run(String hipoFile, int sector) {

        HipoDataSource reader = new HipoDataSource();
        runQuietly(() -> reader.open(hipoFile));

        // per-wire clustered-timing histograms, used only for fitting (never drawn)
        Map<String, H1F> clustMap = new LinkedHashMap<>();

        int eventCount = 0;
        while (reader.hasEvent()) {
            DataEvent event = reader.getNextEvent();
            eventCount++;

            if (!event.hasBank("AHDC::adc") || !event.hasBank("REC::Event")) continue;

            double vertexZ = eventVertexZ(event);
            Map<Integer, List<Hit>> hitsByLayer = extractHitsByLayer(event, vertexZ);
            if (hitsByLayer.isEmpty()) continue;

            // ---- SECTION 3 call: clustering ----
            List<List<Hit>> clusters = findClusters(hitsByLayer);

            // ---- fill per-wire clustered histograms ----
            for (List<Hit> cluster : clusters) {
                for (Hit h : cluster) {
                    if (!h.passesQualityCuts()) continue;
                    String key = h.key();
                    clustMap.computeIfAbsent(key, k -> {
                        H1F h1 = new H1F("hClust_" + k, "Cluster Timing " + k, 200, 0, 800);
                        h1.setTitleX("LeadTime - StartTime (ns)");
                        h1.setTitleY("Counts");
                        return h1;
                    });
                    clustMap.get(key).fill(h.time);
                }
            }
        }
        reader.close();

        // ---- SECTION 4 call: fit every wire ----
        Map<Integer, List<Double>> t0ValuesPerLayer = new TreeMap<>();
        Map<Integer, Integer> acceptedPerLayer = new TreeMap<>();
        int totalWires = 0, acceptedWires = 0, rejectedWires = 0;

        for (String key : clustMap.keySet()) {
            totalWires++;
            WireFitResult result = fitWire(key, clustMap.get(key));

            if (!result.accepted) {
                rejectedWires++;
                continue;
            }
            acceptedWires++;
            int layerCode = Integer.parseInt(key.split("_")[0]);
            t0ValuesPerLayer.computeIfAbsent(layerCode, k -> new ArrayList<>()).add(result.t0);
            acceptedPerLayer.merge(layerCode, 1, Integer::sum);
        }

        // ---- SECTION 5 call: print mean/RMS/%fit per layer ----
        printSummary(t0ValuesPerLayer, acceptedPerLayer, sector, totalWires, acceptedWires, rejectedWires);
    }

    public static void main(String[] args) {

        silenceLibraryLogging(); // turn off the timestamped Minuit/library log spam
        String hipoFile = args.length > 0 ? args[0] : HIPO_FILE;
        int sector = args.length > 1 ? Integer.parseInt(args[1]) : SECTOR;

        T0Calibration calib = new T0Calibration();
        calib.run(hipoFile, sector);
    }
}
