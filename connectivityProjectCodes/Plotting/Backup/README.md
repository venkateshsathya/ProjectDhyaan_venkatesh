# Granger causality versus electrode distance

Run `run_gc_distance_example.m` in MATLAB with `plot_gc_distance.m` and `actiCap64_UOL_angles.csv` in the same folder. Edit `subject` and `protocols` at the top. The example defaults to 019CKa / EO1 and produces separate outgoing and incoming figures. Uncomment the eight-protocol list for a 4-by-2 grid in each direction. No FieldTrip or EEGLAB installation is needed to plot these saved results.

## Data actually inspected

All eight requested subject files exist. They contain numeric `connPre` and `connPost` arrays, each 64-by-64-by-201, and `freqPre`, `freqPost`, `numGoodTrials`. Frequency samples are 0:200 Hz. Labels and geometry are absent from the GC files; the function loads the matching `BK1/data/ftData/<subject>/<protocol>_ep_v8.mat`, using `data.label` and label-matched `data.elec.chanpos`. The seed indices for this subject are Oz=17, O2=18, POz=48, O1=16. These indices are discovered by label, never hardcoded.

The inspected `saveConnData.m` passes the companion data to FieldTrip, stores `abs(temp.grangerspctrm)`, and masks `data.badElecs` rows and columns with NaN. The new plotting code uses those saved values without an additional absolute value, normalization or symmetrization. The old absolute-value operation cannot be undone from saved data. The files do not retain the FieldTrip version or output labels; channel order is inferred from the inspected saving pipeline and its companion input.

FieldTrip's inspected full-matrix GC implementation stores the driver in the first index: `C(source,target,f)`. Therefore outgoing uses `C(seed,:,f)` and incoming uses `C(:,seed,f)`. `MatrixOrder='target-source'` supports explicitly transposed input, if needed. Reference: https://github.com/fieldtrip/fieldtrip/blob/master/connectivity/ft_connectivity_granger.m

## Matching Figure 1(d)

The attached paper plots PPC; this code adopts its display and binning for **GC**, whose numerical values and scale are different. Defaults: four colored seed curves, a black equal-seed mean, 20–32 Hz inclusive (13 saved samples), six equal cosine bins with edges `linspace(-1,1,7)`, and a reversed x-axis from near (+1) to far (-1). Internal result arrays run from far to near. Equal cosine intervals are **not** equal 30-degree angular intervals.

`Distance='project'` follows the local `getElectrodeGroupsConn.m`: cosine of the Euclidean norm of wrapped azimuth and elevation differences, with Cz azimuth difference set to zero. The bundled CSV contains exact `sph_theta`/`sph_phi` values from the local actiCap64_UOL montage, mapped using its companion label table. These differ slightly from angles reconstructed from the saved Cartesian coordinates (and by 90 degrees in the arbitrary azimuth at Cz), so project mode uses the original montage angles. For another cap, supply its label/angle CSV through `AngleFile`. `Distance='spherical'` instead computes the dot product of unit electrode-position vectors, which is the true cosine of their 3-D angular separation. Both assume coordinates relative to the cap sphere origin, as in these companion files. Do not apply the spherical option to translated head coordinates without centering them first.

`Epoch='auto'` uses the mean of Pre/Post for EO1, EC1, M1, EO2, EC2 and Post for G1, G2, M2, following the paper's stimulus versus non-stimulus treatment. Explicit `pre`, `post`, or `mean` overrides this. This averages saved GC estimates; it does not recompute GC from pooled data.

Self-connections are excluded. Other occipital seeds remain eligible target/source electrodes. A pair contributes only if every selected frequency and requested epoch is finite. Within each seed/bin, electrodes are equally weighted. Empty bins remain NaN; the black curve is shown only where all four seeds have valid bin means. No interpolation, zero filling, standard-error bars or subject-level inference is applied.

## Returned values and missing data

`results(q)` includes `pairOutgoing`, `pairIncoming`, `cosDistance`, `bin`, `outgoing`, `incoming`, `countOutgoing`, `countIncoming`, `frequencies`, and source file paths. Each binned matrix is 4-by-6 in seed order Oz/O2/POz/O1. Counts are electrode-pair counts, not trials or independent observations. The example includes optional PNG and MAT export lines.

M1 has only 10,050 finite elements per saved epoch, versus roughly 500,000 in most other protocols. Its occipital off-diagonal connections are unusable in the requested band; expect an empty M1 panel, not a zero-valued curve. G2 also has no valid outgoing O1 estimates in this band. These are properties of the saved data requiring upstream investigation.

## Validation limits

The files, source code, montage and paper figure were inspected directly. An independent numerical calculation checked the real-data bin memberships and outgoing/incoming summaries. MATLAB is installed here but its batch process exits before executing commands, so the MATLAB files and figure rendering could not be run in this environment. `test_gc_distance.m` provides runnable synthetic checks for directed orientation, transposed storage, reordered electrode geometry, self exclusion, cosine endpoints and missing frequency samples. Run it in your working MATLAB session, then run the example.
