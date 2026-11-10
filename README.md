# GRASP DCE Viewer (MATLAB)

3-plane viewer and ROI analysis platform for dynamic GRASP DCE-MRI reconstructions
stored as MATLAB `.mat` files (deep learning subspace recon: `vk_all` + `basis_all`).

- No toolboxes required (own NIfTI writer/reader, own drawing tools). MATLAB R2019b or
  newer recommended; `.xlsx` export needs R2019a+ (`writecell`/`writematrix`), `.csv` works anywhere.
- Runs fully offline. Nothing is written except the files you export.

![Frame view with SSS control and enhancing-tumor ROIs and their time-intensity curves](docs/screenshot_main.png)

![Subtraction map (last - first frame) with the enhancing rim grown from the crosshair](docs/screenshot_subtraction.png)

![Auto SSS: ROI in the posterior superior sagittal sinus just above the torcula](docs/screenshot_sss.png)

![Peak enhancement normalized to the SSS as a colour map](docs/screenshot_params.png)

![Parameters window: settings, colour scale and tumor ROI statistics](docs/screenshot_param_window.png)

*Synthetic phantoms, rendered with GNU Octave; controls look slightly different in MATLAB.
Top: frame view after 2x temporal binning with the SSS control and the tumor ROI. Middle:
subtraction map (frame 60 - frame 1), threshold preview in green and the enhancing rim grown
from the crosshair (necrotic core excluded). Bottom: Auto SSS on a venous phantom - the ROI
(red) sits in the posterior SSS just above the torcula; the straight sinus (anterior to it on
the axial view), the transverse sinuses and a scalp vein are excluded. Last two: peak
enhancement / SSS as a jet colour map (0-0.8) with the ROI outlines, and the Parameters
window with the statistics of the tumor ROI.*

## Quick start

```matlab
addpath('/path/to/GRASP_DCE_Viewer')
GraspDCEViewer                          % then File > Open .mat (or the Open button)
GraspDCEViewer('recon.mat')             % opens the load dialog for this file
GraspDCEViewer(vk_all, basis_all)       % variables already in the workspace
GraspDCEViewer(data1)                   % output of ReadData.m ([nx ny t nz])
GraspDCEViewer('recon.mat', 'Voxel', [1 1 1], 'dt', 3.3, 'NoDialog', true)
GraspDCEViewer('recon.mat', 'BinSize', 2)   % average every 2 frames on load
```

The load dialog lists the variables in the file, detects `vk_all`/`basis_all`, shows the
resulting matrix/K/T and a memory estimate, and asks for the acquisition geometry: voxel
orientation (default LPS, see below), pixel spacing, slice spacing and temporal resolution.
The geometry is remembered between sessions and can be changed later under View > Geometry.

## How frames are computed (identical to ReadData.m)

For slice z and frame t:

    data1(:,:,t,z) = abs( reshape(vk_all(z,:,:,:), nx*ny, K) * Phi_z' )     Phi_z = basis_all(z,:,:) as [T x K]

i.e. the `TempPCASub` adjoint (`PCA' * Recon`), followed by the `CropImg` central crop
(factor 2, removes the 2x readout oversampling) and division by the global maximum over
all slices and frames (taken before the crop, as in ReadData.m). The viewer keeps only
the cropped coefficients and synthesizes frames, slices, voxel time courses and maps on
demand. A ~2 GB file therefore needs roughly 0.5 GB after loading (peak ~2.5 GB while
reading), and the uncropped 4D array that ReadData.m builds is never created.

`tests/gdce_validate_readdata('recon.mat')` re-runs the original ReadData.m arithmetic
with your `@TempPCASub` class on a real file and reports the difference (expected
~1e-7, single-precision rounding). Basis arrays stored as `[slices K T]` or a single
shared `[T K]` basis are also accepted.

## Orientation

The stored array `data1(:,:,t,z)` has rows, columns and slices (partitions). Its orientation
is described by three **axis codes**: the anatomical direction in which the row, column and
slice index increases (L/R, P/A, S/I), as in NIfTI.

| Code | Rows | Columns | Slices | Use |
|---|---|---|---|---|
| **LPS** (default) | R -> L | A -> P | I -> S | scanner x, y, z order; the deep-learning GRASP reconstructions are stored this way: `imagesc(data1(:,:,t,z))` shows the axial slice with anterior on the **left** |
| PLS | A -> P | R -> L | I -> S | arrays that `imagesc` already shows radiologically |
| ILP | S -> I | R -> L | A -> P | coronal slices |
| IPL | S -> I | A -> P | R -> L | sagittal slices |

Earlier versions assumed PLS, which showed the GRASP brains rotated 90 degrees (anterior to
the left on the axial view) with coronal and sagittal anatomy swapped between the two
reformat panels. LPS is now the default for new files.

**Left-right cannot be judged from symmetric anatomy.** LPS follows from the anterior-left,
superior-up appearance of the stored data (scanner axis order); confirm left-right once per
reconstruction pipeline with a lateralized finding (lesion side in the report) or a marker.
If it is mirrored, press Flip L-R.

View > Geometry shows the three codes as menus, preset buttons (LPS, PLS, ILP, IPL) and
**Rotate axial 90 CW**, **Flip L-R / A-P / S-I** buttons that apply immediately. The
setting is remembered for the next files, stored in ROI sessions and in the NIfTI JSON
sidecar. ROIs follow orientation changes. After loading, a quick check warns in the status
bar if the head looks wider left-right than front-back (view probably rotated 90 degrees).

## Coordinates

All exports use one canonical voxel order: i = Right->Left, j = Anterior->Posterior,
k = Inferior->Superior (LPS), with a NIfTI sform/qform that maps voxels to RAS mm with the
volume centred at 0. Image, mask and map NIfTIs therefore overlay in ITK-SNAP, 3D Slicer,
FSLeyes, nibabel etc. The geometry is approximate (voxel sizes and orientation from the
dialog, not scanner coordinates); masks are meant to be used with the image NIfTI exported
from this tool.

## Window layout

- Left: data info, segmentation tools, threshold, display, export.
- Top: axial / coronal / sagittal views with crosshair and orientation labels.
- Middle strip: frame slider (1-frame steps), first/prev/next/last buttons, typed frame
  number, temporal resolution (s/frame), time of the frame, play/pause with fps.
- Bottom: time-intensity curves (all visible ROIs overlaid) and the ROI manager.

Mouse: left = active tool, right-drag = window/level, shift- or middle-drag = pan,
wheel = slice, Ctrl+wheel = zoom, click on the TIC plot = jump to that frame.
Keys: Left/Right = frame -/+1 (Shift: 10), Home/End, PgUp/PgDn, Up/Down = slice,
Space = play, N/C/B/E = tools, [ ] = brush radius, D = 3D toggle, P = threshold preview,
S = subtraction view, G = grow from crosshair, M = parameter colour map on/off,
Ctrl+Z / Ctrl+Y = undo / redo. Help > Mouse and keyboard lists everything.

## Segmentation tools

- Navigate: click/drag moves the crosshair; the crosshair voxel's curve is shown dashed in
  the TIC plot (Cursor TIC), useful for locating the superior sagittal sinus.
- Circle: press at the centre, drag to the radius (radius shown in mm, anisotropic voxels
  handled). With 3D on it fills a sphere.
- Brush / Eraser: radius in mm; 2D (disk in the plane of the view) or 3D (sphere).
- Exclusive ROIs (default on): painting removes the voxels from other ROIs, so the label
  map has no overlap.
- Undo / redo: 100 steps, stored as voxel diffs (Ctrl+Z, Ctrl+Y).
- Multiple ROIs: `Auto SSS` (below), `+ SSS` (empty SSS control ROI to draw by hand),
  `+ Tumor` (auto-numbered Tumor 1..n), `+ Custom`. Rename, recolor, show/hide, delete;
  the table shows voxels and mL.
- Clear slice / Clear ROI, Fill holes (slice by slice in the plane of the last-used view).

## Automatic superior sagittal sinus ROI (Auto SSS)

**Auto SSS** (ROI panel, or Edit > Detect superior sagittal sinus) places the control ROI
`SSS (auto)` in the thickest part of the posterior superior sagittal sinus just above the
torcula, excluding the transverse sinuses and the straight sinus. Requires correct
orientation (LPS etc., see Orientation).

1. Peak enhancement map (maximum over time minus baseline) and a head mask (Otsu on the
   baseline image).
2. Venous structure: voxels enhancing >= 35% of the brightest voxel near the posterior
   midline; the largest 3D connected region grown from the brightest posterior-midline
   voxels (the SSS-torcula-transverse sinus complex; arteries and scalp veins form small
   separate regions).
3. Torcula level: the most superior axial slice where this structure is wider than 30 mm
   left-right (transverse sinuses leaving the midline).
4. On slices 5-40 mm above the torcula where the structure stays within 20 mm of width
   (midline SSS only), the SSS cross-section is its most posterior part against the
   occipital inner table, no deeper front-to-back than it is wide (excludes the straight
   sinus, which joins from anterior). The thickest cross-section (3-slice average) is chosen.
5. ROI: 5 mm of slices centred there; in each, the cross-section voxels enhancing >= 50% of
   their maximum (partial-volume edges excluded).

The crosshair moves to the ROI and the status bar reports the torcula slice, ROI slices,
voxels and volume. Running it again replaces `SSS (auto)`; Undo removes it; the eraser and
brush can edit it. If the transverse sinuses are not in the slab, the lowest posterior SSS
segment is used and the status bar says so. Always check the ROI on the three views.
Scripting: `[lin, info] = gdce_detect_sss(D)` (options for every threshold in `help
gdce_detect_sss`).

## Enhancing tumor from the subtraction map

The subtraction map is `frame(post) - frame(pre)`, by default **last GRASP volume minus
first GRASP volume** (post = last, pre = first; both editable). It is the default threshold
map, and **Show** (or key S) displays it in all three views instead of the current frame,
with its own window/level.

1. `+ Tumor`, Map = *Subtraction: post - pre frame*, press **Show**.
2. Set the window: type Low (High = Inf), drag the slider, or draw a rough sphere around the
   lesion and press **Otsu**. Tick **Preview** to see in-window voxels in green.
3. Turn the subtracted enhancement into the ROI mask, either:
   - **Grow from crosshair** (key G): put the crosshair on the enhancing rim (Navigate) and
     add the 3D, 6-connected region of in-window voxels touching it, up to *max* mm from the
     crosshair. Vessels, choroid plexus and dura that do not touch the lesion stay out, and a
     non-enhancing necrotic core is excluded.
   - **Keep: 3D / Keep: slice**: draw a generous circle/sphere first, then keep only its
     in-window voxels.
   - **All -> ROI**: add every in-window voxel in the volume (includes vessels; erase or
     intersect afterwards).
4. Undo/redo, eraser and Fill holes work as usual on the result.

## Threshold maps

1. Signal at current frame
2. Subtraction: post - pre frame (default last - first)
3. Current frame - baseline
4. Current frame relative enhancement (%)
5. Peak - baseline (maximum after the baseline frames)
6. Peak relative enhancement (%)

Baseline = mean of the chosen pre-contrast frames (default 1-5). Maps 3-6 are computed on
first use (one pass over the data, a few seconds). Every map works with Preview, Otsu,
Keep: slice/3D, Grow from crosshair, All -> ROI and "Constrain brush/circle" (painting only
adds in-window voxels). The line under the buttons shows how many voxels of the active ROI
the current window keeps.

## Temporal resolution ("time unknown")

The temporal resolution is the time per reconstructed frame in seconds (spokes per frame x
TR). It is not stored in the `.mat` file, so it is unknown until you enter it; the viewer
then shows "(time unknown)" / "s/frame unknown" and plots curves against frame number.
Everything else (frames, maps, ROIs, exports) works without it.

Enter it in any of: the **s/frame** box in the frame strip (next to the frame number), the
load dialog, or View > Geometry. It is remembered for the next files. With it set, the
frame strip shows the time of the current frame, the TIC X axis can show seconds, time to
peak is exported in seconds and the NIfTI/JSON/spreadsheet exports carry frame times.
Enter the resolution of the original frames; after binning the viewer multiplies it by n.

## Temporal binning

Binning lowers the temporal resolution by averaging the signal intensity (magnitude) of n
consecutive reconstructed frames into one frame, e.g. to compare an 8-spoke reconstruction
binned x2 or x4 with the 16- and 32-spoke reconstructions of the same scan.

How to:

1. Open the reconstruction and (optionally) enter the s/frame of the original frames.
2. In the DATA section type n in **Average every [n] frames** (e.g. 2).
3. Press **Apply binning** (or Enter in the box). The status bar reports the result, e.g.
   `131 raw frames -> 65 frames of 2 averaged, 6.8 s/frame. Raw frame 131 dropped`.
4. Work as usual: views, slider, curves, enhancement maps and the subtraction map
   (last and first binned frames) all use the binned series. The current frame and the
   baseline frames are moved to the same time points.
5. Export. Image NIfTI: n x dt in pixdim(4), frame times and bin layout in the JSON sidecar.
   Spreadsheet: the file name gets the frame count, e.g.
   `case_roi_timeseries_65frames.xlsx` (CSV: `..._65frames_voxels.csv`), and the Info sheet
   records the bin size and dropped frames.
6. To return to the original series, set n = 1 and press Apply binning.

Details: frame i of the binned series is the mean of raw frames (i-1)n+1 ... in, and its
time is the mean of their times (bin centre). If the number of frames is not a multiple
of n, the incomplete last bin is dropped (reported); View > Temporal binning can keep it as
a shorter last bin instead. ROIs are not affected (they are spatial). Batch:
`gdce_convert_mat2nifti(file, '', 'BinSize', n)`; scripting: `D.setBinning(n)`.
Binning averages reconstructed magnitude frames; it is not a reconstruction with more
spokes per frame.

## Time-intensity curves

Mean +/- SD band, median [IQR] band, or mean only; Y as signal, relative enhancement
(%), S/S0 or S-S0 (S0 = mean of the ROI curve over the baseline frames); X in seconds
(if temporal resolution is set) or frames. File > Save TIC plot as PNG.

## Semi-quantitative parameters (wash-in, washout, AUC, peak enhancement)

Analysis > Semi-quantitative parameters (or the **Parameters...** button above the curves) opens the
parameter window. Press **Compute parameters** (one pass over the data, a few seconds); the
voxel maps, the colour map and the statistics are then available. Voxel values and ROI
values are computed from the same time-intensity curve S(t) on the current (binned)
series.

Definitions (S0 = mean of the baseline frames; times are measured from the last baseline
frame, i.e. the injection; seconds if s/frame is set, otherwise frames):

| Parameter | Definition | Unit |
|---|---|---|
| Peak enhancement (PE) | Smax - S0, Smax = maximum after the baseline | a.u. |
| Relative peak enhancement | 100 x (Smax - S0) / S0 | % |
| Wash-in slope (WIS) | PE / (TP - To): line from the onset To to the peak TP (Khalifa 2014). Option: maximum up-slope between consecutive frames before the peak | a.u./s |
| Washout slope (WOS) | least-squares slope from the peak to the last frame (at least the last 3 frames, so a curve still rising at the end gets a positive slope). Option: 2-point (S_last - Smax)/(t_last - TP). Negative = washout | a.u./s |
| Area under the curve (AUC) | trapezoidal area of S - S0 from the end of the baseline to the end, or over the first N s (AUC window 60 = IAUC60, 90 = IAUC90) | a.u. x s |
| Time to peak (TTP) | TP measured from the end of the baseline | s |

Onset To = last frame before the peak whose enhancement is below 10% of PE.

**Normalization to the SSS.** "/SSS" parameters are the voxel value divided by the same
parameter of the reference ROI mean curve, in line with the AUC ratio to a venous/vascular
reference used for subtle blood-brain barrier leakage (Heye et al. 2016). The reference is
the first ROI whose name contains "SSS" (Auto SSS creates "SSS (auto)"), or the ROI chosen
under SSS reference. Normalization cancels the arbitrary scaling of the reconstruction
(global-maximum normalization, coil and receiver gain) and makes values comparable between
scans and between the 8/16/32-spoke reconstructions of one scan. The washout ratio is a
ratio of slopes: the SSS normally washes out (negative slope), so a positive WOS/SSS means
a tumor washout, a negative one a persistent rise.

**Colour maps.** Choose the parameter, Show on images (or key M), colour map (Jet or Turbo
rainbow, Hot, Gray), Min/Max (type values, Auto head = 1st-99.5th percentile in the head,
Auto ROIs = 2nd-98th percentile in the tumor ROIs), opacity, Hide below min and Inside head
only. The colour bar appears in the sagittal panel; the status bar shows the value under
the cursor. ROIs are drawn as outlines while a map is shown.

**Statistics.** The table lists, for every tumor ROI (all ROIs except the reference; tick
Include the reference to add it) and every parameter: N voxels, mean, SD, median, Q1, Q3,
10th and 90th percentile, min, max (voxel values; NaN voxels are ignored) and the
parameter of the ROI mean curve ("ROI curve"). Voxel-wise statistics describe
heterogeneity; the ROI-curve value is less noisy.

**Export.** Export statistics + voxel values writes `.xlsx` (sheet Stats, sheet Settings
with the definitions, options and the reference values, one sheet per ROI with one row
per voxel: Voxel, i, j, k, x_mm, y_mm, z_mm and one column per parameter) or `.csv`
(`_param_stats.csv`, `_param_voxels.csv`, `_param_settings.txt`); after binning the file
name ends in `_<N>frames`. Export parameter maps writes every parameter as a NIfTI map
(`<base>_PE.nii.gz`, `<base>_PE_SSS.nii.gz`, ...) in the geometry of the image export.
ROI sessions (Save ROIs) keep the parameter options and the reference ROI.

Scripting:

```matlab
app.computeParams('WashIn', 'onset', 'WashOut', 'regression', 'AUCWindow', 60);
[hdr, rows] = app.paramStatsTable();       % cell table as shown in the window
app.setParamKey('WIS_SSS'); app.setParamDisplay(true); app.setParamCLim([0 2]);
app.exportParams('tumor_parameters.xlsx');
P = gdce_curve_params(TC, D.paramTimes(), 1:5);   % any curves [N x T]
```

Other semi-quantitative descriptors used in DCE-MRI that are not (yet) included: bolus
arrival / onset time, time to 90% of the peak (T90), maximum intensity-time ratio (MITR =
PE / TTP), signal enhancement ratio (SER = (S1 - S0)/(S2 - S0), early vs late post-contrast
frame), relative signal intensity at fixed times, plateau level and washout percentage
((Smax - S_end)/Smax), and the curve-type classification (I persistent, II plateau, III
washout). Several are one-line additions to `gdce_curve_params.m`.

## Exports

| Button / menu | Output |
|---|---|
| Image -> NIfTI | 4D float32 `.nii.gz`/`.nii` (all frames, streamed frame by frame) or the current frame (3D), plus a JSON sidecar (frame times, source, geometry, normalization) |
| Mask -> NIfTI | uint8 label map (value = ROI number), `<name>_labels.txt` (label, name, RGB, voxels, mL), `<name>_itksnap_labels.txt`; optionally one binary mask per ROI |
| ROI voxel time series -> spreadsheet | `.xlsx`: sheets Info, Summary (per-ROI mean, SD, median, relative enhancement, n per frame) and one sheet per ROI with one row per voxel: Voxel, i, j, k, x_mm, y_mm, z_mm, F001..Fnnn. `.csv`: `_voxels.csv` (all ROIs, ROI_ID/ROI_Name columns), `_summary.csv`, `_info.txt`. After binning the file name ends in `_<N>frames` |
| File > Export enhancement maps / Parameters > Export parameter maps | Baseline, Max and every parameter (PE, RE, WIS, WOS, AUC, TTP; PE, WIS, WOS, AUC also /SSS) as NIfTI |
| Parameters > Export statistics + voxel values | parameter statistics per tumor ROI and voxel values (see above) |
| (binning) | binned series everywhere; JSON sidecar and Info sheet record bin size, raw frames per frame and dropped frames |
| Save / Load ROIs | `.mat` session (ROI names, colors, voxels in native indices, geometry); reloads onto any reconstruction with the same matrix, e.g. the 8-, 16- and 32-spoke versions of one scan |
| Import mask | label NIfTI from this tool or edited elsewhere; axis permutations/flips are corrected from the header; names/colors are read from `_labels.txt` if present |

When a new file with the same matrix size is opened, the viewer offers to keep the
current ROIs, so the same ROIs can be applied across reconstructions.

i, j, k in the spreadsheet are 1-based indices into the exported NIfTI arrays (as returned
by `niftiread` in MATLAB); x, y, z are RAS mm from the affine. Excel limits (1,048,575
voxels per ROI, 16,384 columns) apply to `.xlsx`; use `.csv` for very large ROIs.

## Scripting and batch use

```matlab
% batch conversion of every reconstruction in a folder to 4D NIfTI
f = dir('*.mat');
for k = 1:numel(f)
    gdce_convert_mat2nifti(f(k).name, '', 'Voxel', [1 1 1], 'dt', 3.3);   % AxisCodes 'LPS' by default
end
gdce_convert_mat2nifti('recon_8spoke.mat', '', 'dt', 1.6, 'BinSize', 2);   % binned x2

% data model only
D = GraspData.fromFile('recon.mat');
D.setAxisCodes('LPS'); D.VoxNative = [1 1 1]; D.dt = 3.3;
D.setBinning(2);                   % optional: average frame pairs
V  = D.canonVolume(40);            % 3D frame 40 (canonical order)
S  = D.canonSubtraction(D.NT, 1);  % last - first
TC = D.timecourses(find(mask));    % voxel x frame matrix for a canonical mask
gdce_export_timeseries(D, {mask}, {'Tumor 1'}, 'tumor.csv', 'BaseFrames', 1:5);
[sss, info] = gdce_detect_sss(D);  % automatic SSS ROI (canonical indices)

% drive the GUI
app = GraspDCEViewer(D);
app.setFrame(45); app.addTumor();
app.setThrSource('sub'); app.setThreshold(0.2, Inf);
app.setPos([90 50 24]); app.growFromCrosshair();
app.autoSSS();
app.exportMaskNifti('mask.nii.gz');
```

Standalone build (MATLAB Compiler): `mcc -m gdce_app.m -a GraspDCEViewer.m -a GraspData.m -a gdce_*.m`.

## Files

| File | Purpose |
|---|---|
| `GraspDCEViewer.m` | GUI (classic figure, no toolboxes) |
| `GraspData.m` | data model: subspace/dense loading, frame synthesis, orientation, maps, affine |
| `gdce_load_dialog.m` | variable/format/geometry dialog |
| `gdce_nifti_write.m`, `gdce_nifti_read.m` | NIfTI-1 I/O (.nii/.nii.gz), streaming 4D writes |
| `gdce_write_dce_nifti.m` | series -> NIfTI + JSON sidecar |
| `gdce_export_timeseries.m` | voxel time series -> .xlsx/.csv |
| `gdce_curve_params.m` | semi-quantitative curve parameters (vectorized over curves) |
| `gdce_export_params.m` | parameter statistics and voxel values -> .xlsx/.csv |
| `gdce_reorient_to_ref.m` | mask import reorientation |
| `gdce_detect_sss.m` | automatic superior sagittal sinus ROI |
| `gdce_convert_mat2nifti.m` | batch .mat -> NIfTI |
| `gdce_app.m` | launcher / compiler entry point |
| `docs/` | README screenshots |
| `legacy/` | the original ReadData.m, CropImg.m and @TempPCASub (reference for validation) |
| `tests/gdce_selftest.m` | synthetic end-to-end test (`gdce_selftest(true)` also drives the GUI) |
| `tests/gdce_validate_readdata.m` | compares the viewer with ReadData.m on a real file |
| `tests/gdce_sss_phantom.m` | venous phantom (SSS, torcula, transverse/straight sinuses, scalp vein, artery) |

## Notes and limitations

- Geometry comes from the dialog; there is no DICOM header import yet. Voxel sizes,
  orientation (especially left-right) and frame timing should be checked once per protocol.
- Auto SSS was tested on synthetic phantoms; check its ROI on real data, in particular
  near the torcula and with a hypoplastic or asymmetric confluence.
- Signal is the magnitude of the reconstruction (normalized to the global maximum when
  ReadData.m normalization is on); no T1 mapping or concentration conversion is done.
  Parameters are therefore semi-quantitative (signal-based): compare absolute values only
  between scans with the same protocol and scaling, or use the /SSS values.
- Voxel-wise slopes are noisy at high temporal resolution; binning or the ROI-curve values
  are more robust.
- Fill holes is 2D (slice by slice). Grow from crosshair is the connected-component tool
  (6-connectivity, bounded by the max distance).
- Binning averages magnitudes of reconstructed frames; it does not re-reconstruct with
  more spokes per frame.
