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
import org.jlab.geom.base.ConstantProvider;
import org.jlab.detector.calib.utils.DatabaseConstantProvider;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;


/**
 * ============================================================================
 * AHDC T0 CALIBRATION
 * ============================================================================
 * SECTION 1: GEOMETRY       - wire geometry
 * SECTION 2: HIT EXTRACTION - HIPO -> hits
 * SECTION 3: CLUSTERING     - hits -> clusters
 * SECTION 4: FITTING        - timing histogram -> T0
 * SECTION 5: SUMMARY        - mean / RMS / % fit
 * SECTION 6: RUNNER         - complete analysis
 * H1F histograms are used only internally for Gaussian fitting.
 * ============================================================================
 */
public class T0Calibration {

    // ------------------------------------------------------------------------
    // QUALITY CUTS
    // ------------------------------------------------------------------------
    private static final int MIN_ADC = 500;
    private static final int MAX_PED = 300;
    private static final double MIN_TOT = 1.0;

    // ------------------------------------------------------------------------
    // CLUSTERING
    // ------------------------------------------------------------------------

    private static final int MIN_HITS_PER_CLUSTER = 7;
    private static final int MIN_SUPERLAYERS = 4;
    private static final double SEED_TIME_WINDOW_NS = 150.0;
    private static final double CANDIDATE_TIME_WINDOW_NS = 150.0;
    private static final double CANDIDATE_PHI_WINDOW = Math.toRadians(25.0);
    private static final double CANDIDATE_MAX_PERP_DIST = 20.0;

    // ------------------------------------------------------------------------
    // FITTING
    // ------------------------------------------------------------------------

    private static final double SEARCH_MIN_NS = 150.0;
    private static final double SEARCH_MAX_NS = 500.0;
    private static final double SEED_WINDOW_NS = 35.0;
    private static final double MIN_ENTRIES = 20;
    private static final double MAX_MEAN_SHIFT_FROM_PEAK = 80.0;
    private static final double MIN_SIGMA = 2.0;
    private static final double MAX_SIGMA = 150.0;
    private static final double BACKGROUND_TOLERANCE = -5.0;

    // ------------------------------------------------------------------------
    // INPUT
    // ------------------------------------------------------------------------

    private static final double VERTEX_Z_OFFSET = 5.0;
    private static final String HIPO_FILE = "/path/to/hipo/file.hipo";
    private static final int SECTOR = 1;

    // ------------------------------------------------------------------------
    // OUTPUT CONTROL
    // ------------------------------------------------------------------------

    private static final PrintStream QUIET_STDOUT = new PrintStream(new OutputStream() {
        @Override
        public void write(int b) {}
        @Override
        public void write(byte[] b, int off, int len) {}
    });

    private static void silenceLibraryLogging() {
        LogManager.getLogManager().reset();
        Logger.getLogger("").setLevel(Level.OFF);
    }

    private static void runQuietly(Runnable action) {
        PrintStream old = System.out;
        System.setOut(QUIET_STDOUT);
        try {
            action.run();
        } finally {
            System.setOut(old);
        }
    }
    // =========================================================================
    // SECTION 1: GEOMETRY
    // =========================================================================

    private final AlertDCDetector detector;
    /*
     * Wire geometry never changes during the run.
     */
    private final Map<Long, Line3D> wireCache = new HashMap<>();

    private T0Calibration() {
        ConstantProvider cp = new DatabaseConstantProvider(1, "default", "default");
        detector = new AlertDCFactory().createDetectorCLAS(cp);
    }
    private long wireKey(int sector, int layer, int wire) {
        return ((long) sector << 32) | ((long) layer << 16) | (wire & 0xffffL);
    }
    private Line3D getWire(int sector, int layerCode, int wire) {

        long key = wireKey(sector, layerCode, wire);
        Line3D line = wireCache.get(key);

        if (line != null) return line;

        try {
            int superlayer = layerCode / 10;
            int layer = layerCode % 10;

            line = detector.getSector(sector)
                    .getSuperlayer(superlayer)
                    .getLayer(layer)
                    .getComponent(wire)
                    .getLine();

            wireCache.put(key, line);
            return line;

        } catch (Exception e) {
            return null;
        }
    }
    private Point3D projectWireToVertexZ(int sector, int layerCode, int wire, double vz) {

        Line3D line = getWire(sector, layerCode, wire);
        if (line == null) return null;

        Point3D p1 = line.origin();
        Point3D p2 = line.end();

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

            return detector.getSector(sector)
                    .getSuperlayer(superlayer)
                    .getLayer(layer)
                    .getNumComponents();

        } catch (Exception e) {
            return 0;
        }
    }
    // =========================================================================
    // SECTION 2: HIT EXTRACTION
    // =========================================================================
    private static class Hit {
        final int id;
        final double x;
        final double y;
        final double phi;
        final double time;
        final int layer;
        final int wire;
        final int adc;
        final float ped;
        final float tot;
        Hit(int id, double x, double y, double time, int layer, int wire, int adc, float ped, float tot) {

            this.id = id;
            this.x = x;
            this.y = y;
            this.phi = Math.atan2(y, x);
            this.time = time;
            this.layer = layer;
            this.wire = wire;
            this.adc = adc;
            this.ped = ped;
            this.tot = tot;
        }
        int wireKey() {
            return (layer << 16) | (wire & 0xffff);
        }
        boolean passesQualityCuts() {
            return adc > MIN_ADC && ped < MAX_PED && tot > MIN_TOT;
        }
    }
    private double eventVertexZ(DataEvent event) {
        double vz = 0.0;
        if (event.hasBank("REC::Particle")) {
            DataBank bank = event.getBank("REC::Particle");
            if (bank.rows() > 0)
                vz = bank.getFloat("vz", 0);
        }
        return vz - VERTEX_Z_OFFSET;
    }

    private Map<Integer, List<Hit>> extractHitsByLayer(DataEvent event, double vertexZ) {
        Map<Integer, List<Hit>> hitsByLayer = new HashMap<>();
        DataBank adcBank = event.getBank("AHDC::adc");
        DataBank eventBank = event.getBank("REC::Event");
        float startTime = eventBank.getFloat("startTime", 0);
        int hitID = 0;
        for (int i = 0; i < adcBank.rows(); i++) {
            short wfType = adcBank.getShort("wfType", i);
            if (wfType == 4 || wfType == 5) continue;
            int sector = adcBank.getInt("sector", i);
            int layer = adcBank.getInt("layer", i);
            int wire = adcBank.getInt("component", i);
            float leadTime = adcBank.getFloat("leadingEdgeTime", i);
            int adc = adcBank.getInt("ADC", i);
            float ped = adcBank.getFloat("ped", i);
            float tot = adcBank.getFloat("timeOverThreshold", i);
            Point3D p = projectWireToVertexZ(sector, layer, wire, vertexZ);

            if (p == null) continue;

            Hit hit = new Hit(hitID++, p.x(), p.y(), leadTime - startTime, layer, wire, adc, ped, tot);

            hitsByLayer.computeIfAbsent(layer, k -> new ArrayList<>()).add(hit);
        }
        return hitsByLayer;
    }
    // =========================================================================
    // SECTION 3: CLUSTERING
    // =========================================================================
    private List<List<Hit>> findClusters(Map<Integer, List<Hit>> hitsByLayer) {
        List<List<Hit>> clusters = new ArrayList<>();
        List<Integer> layers = new ArrayList<>(hitsByLayer.keySet());
        Collections.sort(layers);
        int totalHits = 0;
        for (List<Hit> hits : hitsByLayer.values())
            totalHits += hits.size();
        boolean[] used = new boolean[totalHits];

        for (int i = 0; i < layers.size(); i++) {
            List<Hit> hits1 = hitsByLayer.get(layers.get(i));

            for (int j = i + 1; j < layers.size(); j++) {

                List<Hit> hits2 = hitsByLayer.get(layers.get(j));

                for (Hit h1 : hits1) {

                    if (used[h1.id]) continue;

                    for (Hit h2 : hits2) {

                        if (used[h1.id] || used[h2.id]) continue;

                        if (Math.abs(h2.time - h1.time) > SEED_TIME_WINDOW_NS) continue;

                        double dx = h2.x - h1.x;
                        double dy = h2.y - h1.y;

                        double length = Math.hypot(dx, dy);
                        if (length == 0) continue;

                        double dirX = dx / length;
                        double dirY = dy / length;

                        List<Hit> cluster = new ArrayList<>();
                        cluster.add(h1);
                        cluster.add(h2);

                        for (int k = j + 1; k < layers.size(); k++) {

                            List<Hit> candidates = hitsByLayer.get(layers.get(k));

                            Hit best = findBestCandidate(h1, dirX, dirY, candidates);

                            if (best != null && !used[best.id]) cluster.add(best);
                        }

                        if (cluster.size() >= MIN_HITS_PER_CLUSTER && countSuperlayers(cluster) >= MIN_SUPERLAYERS) {
                            clusters.add(cluster);
                            for (Hit h : cluster) used[h.id] = true;
                        }
                    }
                }
            }
        }
        return clusters;
    }
    private Hit findBestCandidate(Hit seed, double dirX, double dirY, List<Hit> candidates) {
        Hit best = null;
        double bestScore = Double.MAX_VALUE;

        for (Hit c : candidates) {

            double dt = c.time - seed.time;

            if (Math.abs(dt) > CANDIDATE_TIME_WINDOW_NS) continue;

            // Both phi values are already calculated in Hit.
            double dphi = seed.phi - c.phi;

            if (dphi > Math.PI) dphi -= 2 * Math.PI;
            else if (dphi < -Math.PI) dphi += 2 * Math.PI;

            if (Math.abs(dphi) > CANDIDATE_PHI_WINDOW) continue;

            double dx = c.x - seed.x;
            double dy = c.y - seed.y;

            // Since dir is a unit vector, this is the perpendicular distance to the seed direction.
            double perpDist = Math.abs(dx * dirY - dy * dirX);

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

        int mask = 0;

        for (Hit h : cluster) mask |= 1 << (h.layer / 10);

        return Integer.bitCount(mask);
    }
    // =========================================================================
    // SECTION 4: FITTING
    // =========================================================================
    private static class WireFitResult {
        final int key;
        final boolean accepted;
        final String rejectReason;
        final double t0;
        WireFitResult(int key, boolean accepted, String rejectReason, double t0) {
            this.key = key;
            this.accepted = accepted;
            this.rejectReason = rejectReason;
            this.t0 = t0;
        }
    }
    private WireFitResult fitWire(int key, H1F hCl) {

        if (hCl.getEntries() < MIN_ENTRIES)
            return new WireFitResult(key, false, "insufficient entries", Double.NaN);

        int firstBin = hCl.getAxis().getBin(SEARCH_MIN_NS);

        int lastBin = hCl.getAxis().getBin(SEARCH_MAX_NS);

        // Find timing peak.
        int maxBin = firstBin;
        double maxContent = -1;

        for (int b = firstBin; b <= lastBin; b++) {

            double content = hCl.getBinContent(b);

            if (content > maxContent) {
                maxContent = content;
                maxBin = b;
            }
        }

        if (maxContent <= 0)
            return new WireFitResult(key, false, "no peak found", Double.NaN);

        double peak = hCl.getAxis().getBinCenter(maxBin);

        // ---------------------------------------------------------------------
        // Gaussian seed
        // ---------------------------------------------------------------------

        double sum = 0;
        double xSum = 0;
        double x2Sum = 0;

        for (int b = firstBin; b <= lastBin; b++) {

            double x = hCl.getAxis().getBinCenter(b);

            double y = hCl.getBinContent(b);

            if (y <= 0 || Math.abs(x - peak) > SEED_WINDOW_NS) continue;

            sum += y;
            xSum += y * x;
            x2Sum += y * x * x;
        }

        double seedMean = peak;
        double seedSigma = 15.0;

        if (sum > 0) {

            seedMean = xSum / sum;

            double variance = x2Sum / sum - seedMean * seedMean;

            if (variance > 0) seedSigma = Math.sqrt(variance);
        }

        seedSigma = Math.max(5.0, Math.min(50.0, seedSigma));
        // ---------------------------------------------------------------------
        // Fit range
        // ---------------------------------------------------------------------
        double fitMin = Math.max(SEARCH_MIN_NS, peak - 100.0);

        double fitMax = Math.min(SEARCH_MAX_NS, peak + 50.0);
        // ---------------------------------------------------------------------
        // Background
        // ---------------------------------------------------------------------

        double bgSum = 0;
        int bgCount = 0;

        for (int b = firstBin; b <= lastBin; b++) {

            double x = hCl.getAxis().getBinCenter(b);

            if ((x >= fitMin && x < fitMin + 5.0) || (x <= fitMax && x > fitMax - 5.0)) {

                bgSum += hCl.getBinContent(b);
                bgCount++;
            }
        }

        double background = bgCount > 0 ? bgSum / bgCount : 0.0;

        double amplitude = Math.max(1.0, maxContent - background);

        // ---------------------------------------------------------------------
        // Gaussian + flat background
        // ---------------------------------------------------------------------

        F1D gaus = new F1D("gaus_" + key, "[bg] + [amp]*gaus(x,[mean],[sigma])", fitMin, fitMax);
        gaus.setParameter(0, background);
        gaus.setParameter(1, amplitude);
        gaus.setParameter(2, seedMean);
        gaus.setParameter(3, seedSigma);

        try {
            hCl.fit(gaus);
        } catch (Exception e) {
            return new WireFitResult(key, false, "fit failed", Double.NaN);
        }

        double fitBackground = gaus.getParameter(0);

        double fitAmplitude = gaus.getParameter(1);

        double mean = gaus.getParameter(2);

        double sigma = Math.abs(gaus.getParameter(3));

        // Keep your original T0 definition.
        double t0 = peak - 50.0;

        if (Math.abs(mean - peak) > MAX_MEAN_SHIFT_FROM_PEAK)

            return new WireFitResult(key, false, "mean moved too far from peak", Double.NaN);

        if (sigma < MIN_SIGMA || sigma > MAX_SIGMA)

            return new WireFitResult(key, false, "sigma out of range", Double.NaN);

        if (fitAmplitude <= 0)

            return new WireFitResult(key, false, "non-positive amplitude", Double.NaN);

        if (fitBackground < BACKGROUND_TOLERANCE)

            return new WireFitResult(key, false, "negative background", Double.NaN);

        return new WireFitResult(key, true, null, t0);
    }
    // =========================================================================
    // SECTION 5: SUMMARY
    // =========================================================================
    private void printSummary(Map<Integer, List<Double>> t0ValuesPerLayer, Map<Integer, Integer> acceptedPerLayer, int sector, int totalWires, int acceptedWires, int rejectedWires) {

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
    // =========================================================================
    // SECTION 6: RUNNER
    // =========================================================================
    private void run(String hipoFile, int sector) {

        HipoDataSource reader = new HipoDataSource();

        runQuietly(() -> reader.open(hipoFile));

        // Integer key = layer + wire.
        Map<Integer, H1F> clustMap = new HashMap<>();

        while (reader.hasEvent()) {

            DataEvent event = reader.getNextEvent();

            if (!event.hasBank("AHDC::adc") || !event.hasBank("REC::Event")) continue;

            double vertexZ = eventVertexZ(event);

            Map<Integer, List<Hit>> hitsByLayer = extractHitsByLayer(event, vertexZ);

            if (hitsByLayer.isEmpty()) continue;

            List<List<Hit>> clusters = findClusters(hitsByLayer);

            // ---------------------------------------------------------------
            // Fill timing histograms
            // ---------------------------------------------------------------

            for (List<Hit> cluster : clusters) {

                for (Hit h : cluster) {

                    if (!h.passesQualityCuts()) continue;

                    int key = h.wireKey();

                    H1F hist = clustMap.get(key);

                    if (hist == null) {

                        hist = new H1F("hClust_" + key, "Cluster Timing", 200, 0, 800);
                        clustMap.put(key, hist);
                    }

                    hist.fill(h.time);
                }
            }
        }

        reader.close();

        // =========================================================================
        // FIT ALL WIRES
        //
        // Suppress Minuit/GROOT output once for the entire loop instead of
        // creating/replacing System.out for every individual wire.
        // =========================================================================

        Map<Integer, List<Double>> t0ValuesPerLayer = new TreeMap<>();

        Map<Integer, Integer> acceptedPerLayer = new TreeMap<>();

        int totalWires = 0;
        int acceptedWires = 0;
        int rejectedWires = 0;

        PrintStream oldOut = System.out;
        System.setOut(QUIET_STDOUT);

        try {

            for (Map.Entry<Integer, H1F> entry : clustMap.entrySet()) {

                totalWires++;

                WireFitResult result = fitWire(entry.getKey(), entry.getValue());

                if (!result.accepted) {
                    rejectedWires++;
                    continue;
                }
                acceptedWires++;

                int layerCode = result.key >>> 16;

                t0ValuesPerLayer.computeIfAbsent(layerCode,k -> new ArrayList<>()).add(result.t0);

                acceptedPerLayer.merge(layerCode, 1, Integer::sum);
            }

        } finally {
            System.setOut(oldOut);
        }

        // =========================================================================
        // PRINT FINAL RESULTS
        // =========================================================================

        printSummary(t0ValuesPerLayer, acceptedPerLayer, sector, totalWires, acceptedWires, rejectedWires);
    }
    // =========================================================================
    // MAIN
    // =========================================================================

    public static void main(String[] args) {

        silenceLibraryLogging();

        String hipoFile = args.length > 0 ? args[0] : HIPO_FILE;

        int sector = args.length > 1 ? Integer.parseInt(args[1]) : SECTOR;

        T0Calibration calibration = new T0Calibration();

        calibration.run(hipoFile, sector);
    }
}
