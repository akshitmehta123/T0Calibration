# AHDC T0 Calibration

This code computes the T0 timing calibration for the AHDC (ALERT Drift
Chamber) wires from a HIPO data file. It reads the raw ADC hits, clusters
them into tracks, fits the timing distribution of each wire, and reports
the calibration quality per layer.

## What the code does

The pipeline is: **read banks → find hits → cluster hits → fit each wire →
print results.** Everything runs headless — no windows, no canvases, no
histograms are ever drawn to the screen. Histograms are used only
internally as a tool to do the Gaussian fit.

Step by step, for every event in the HIPO file:

1. **Hit extraction** — Reads the `AHDC::adc` bank for raw hits and the
   `AHDC::hits` bank for which hits the official tracking used. Each hit's
   wire position is projected onto the plane at the event's vertex Z
   (`REC::Particle`), and hits from bad waveforms are dropped.
2. **Clustering** — Groups hits into candidate tracks based on timing,
   angular (phi), and spatial-distance consistency. Only clusters that hit
   enough wires across enough superlayers are kept.
3. **Histogram filling** — For every wire, the timing (`leadingEdgeTime -
   startTime`) of every hit belonging to an accepted cluster is filled into
   a per-wire histogram, after quality cuts on ADC, pedestal, and
   time-over-threshold.
4. **Fitting** — Each wire's histogram is fit with a Gaussian on top of a
   flat background, seeded from the histogram's own peak. Fits that don't
   meet quality checks (peak shift, sigma range, amplitude, background)
   are rejected and excluded from the results.
5. **Summary** — For each detector layer, the code computes the **mean
   T0**, the **RMS of T0** across that layer's wires, and the **percentage
   of that layer's wires with an accepted fit.** These three numbers per
   layer are the only thing printed.

The code is organized into six clearly labeled sections inside a single
file (`T0Calibration.java`) so that any one piece — the clustering
algorithm, the fitting logic, etc. — can be read and pulled out on its own
without needing the rest of the file.

## Output

Running the code prints a table like this:


==================================================
              T0 FIT SUMMARY
==================================================
Total wires tested = 377
Accepted fits      = 327
Rejected fits      = 50

Layer | MeanT0(ns) | RMS(ns) | %WiresFit
   11 |    164.000 |  16.452 |      51.1
   21 |    165.692 |  11.645 |      46.4
   ...


All internal library log/console noise (from the HIPO reader and the
Minuit fitter) is suppressed, so this table is the only substantial output.

## How to use it

1. Open `T0Calibration.java` and edit the `HIPO_FILE` constant near the
   top of the class to point at the `.hipo` file you want to run over.
   

  
   private static final String HIPO_FILE = "/path/to/your/file.hipo";
   
   

2. Compile and run it like any other class in this analysis package, e.g.:

  
   java org.clas.modules.analysis.T0Calibration
  

3. Read the printed summary table for the mean T0, RMS, and percent of
   wires fit per layer.

## Tuning

All of the cuts and thresholds used by the clustering and fitting logic
are declared as named constants at the top of the class (quality cuts,
clustering windows, fit search ranges, quality-check tolerances). Change
these directly rather than hunting for magic numbers elsewhere in the file.

