# Line-noise and taper simulation

Run `demo_line_noise_taper_leakage.m` in MATLAB with FieldTrip and the Signal Processing Toolbox on the path. Inputs are at the top. The script generates the plots, CSV tables and `simulation_results.mat` in this folder. Re-running replaces these demo outputs.

## Settings and calibration

The pipeline defaults are a 1-second epoch, `method='mtmfft'`, `taper='dpss'` and `tapsmofrq=1` Hz. This gives NW=1 and **one retained DPSS taper**, following FieldTrip's construction. The simulation assumes a 1000-Hz sample rate; set `fs` to your recording's `data.fsample` if different.

The simulated clean signal has colored and white noise plus rhythms at 10, 40 and 55 Hz. Its RMS is calibrated to 5 microvolts. `lineToCleanPowerDb=20` gives line power 100 times the total clean signal power, or line RMS of 50 microvolts. This is a controlled simulation, not an amplitude fitted to a recording. Change this setting to test weaker or stronger interference. The strength sweep runs from -20 to +40 dB.

The same clean trials are reused in both scenarios:

- Exactly 50 Hz: an integer number of cycles in a 1-second epoch.
- 50.3 Hz: a stationary tone slightly offset from the FFT grid and the removal frequency.

This does not simulate fluctuating line frequency or amplitude within an epoch.

## What is compared

The figures compare clean/noisy PSDs with a rectangular window and the pipeline DPSS window. They also show fitting and subtracting 50/100-Hz sinusoids **before** DPSS, using the same `ft_preproc_dftfilter(...,'dftreplace','zero')` call as the DFT connectivity code. A filtered-clean reference reveals the clean signal power removed by this operation.

Line-only PSDs measure leakage directly. Subtracting clean PSD from noisy PSD is not an exact measure of line leakage because finite-sample signal/noise cross terms remain. The separate band-error table does compare the complete contaminated or filtered PSD with the corresponding estimator's unfiltered clean PSD.

The calculation uses one-sided PSD in microvolt squared/Hz. Summing bins times their width gives mean-square power; multiplying by epoch duration gives epoch energy. The kernel figure uses zero padding only to draw the spectral window densely; the simulation spectra and numerical summaries retain the native 1-Hz frequency grid.

## Default results

| Line frequency | Window | Line power outside 46–54 Hz, as % of total line spectral power |
|---|---|---:|
| 50 Hz | Rectangular | Approximately zero |
| 50 Hz | DPSS | 0.01128% |
| 50.3 Hz | Rectangular | 2.93943% |
| 50.3 Hz | DPSS | 0.31822% |

For exactly 50 Hz, DPSS places 10.5877% of the line's spectral power outside the 50-Hz bin, mostly nearby. Rectangular-window cancellation at the other native bins depends on the ideal stationary tone being exactly on the FFT grid. The off-grid example shows why this result does not mean that omitting tapers generally prevents leakage: DPSS substantially reduces distant leakage for 50.3 Hz.

At the chosen +20-dB strength, the DPSS line contribution remaining in 0–200 Hz after excluding 46–54 Hz is 1.22% of the clean power in those retained bins for 50 Hz, and 33.62% for 50.3 Hz. A small leakage fraction of a strong line can therefore still matter. The fixed 50-Hz sinusoid fit removes the exact tone to numerical precision, but leaves substantial residual contamination for 50.3 Hz.

These are PSD results; they do not quantify Granger causality bias. No connectivity is computed.

## Files and validation

- `scenario_50Hz.png` and `scenario_50.3Hz.png`: time signals, clean/noisy spectra, isolated line leakage, filtering and band-power comparisons.
- `taper_kernels_and_strength_sweep.png`: taper shapes, spectral windows and retained contamination versus noise strength. Values near -250 dB or lower represent floating-point residuals, not meaningful physical leakage estimates.
- `leakage_summary.csv`: isolated line power and leakage fractions. The masks are nested and are not meant to sum to 100%. Filtered-row percentages use the original unfiltered DPSS line spectral power as denominator.
- `band_power_error.csv`: signed band-power errors and power ratios relative to clean data, using half-open bands [low, high).
- `simulation_results.mat`: settings, taper vectors, spectra, sweep results, tables and the FieldTrip validation result.

The explicit PSD estimator was checked against `ft_freqanalysis` with relative error approximately 1.02e-15. Runtime assertions also check injected power, rectangular-window Parseval consistency and removal of the known exact stationary tone.

References: [FieldTrip taper and smoothing documentation](https://www.fieldtriptoolbox.org/tutorial/sensor/timefrequencyanalysis/) and [MATLAB DPSS documentation](https://www.mathworks.com/help/signal/ref/dpss.html).
