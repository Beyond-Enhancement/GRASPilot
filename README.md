# GRASP DCE Viewer (MATLAB)

3-plane viewer and ROI analysis platform for dynamic GRASP DCE-MRI reconstructions
stored as MATLAB `.mat` files (deep learning subspace recon: `vk_all` + `basis_all`).

- No toolboxes required (own NIfTI writer/reader, own drawing tools). MATLAB R2019b or
  newer recommended; `.xlsx` export needs R2019a+ (`writecell`/`writematrix`), `.csv` works anywhere.
- Runs fully offline. Nothing is written except the files you export.

![Frame view with SSS control and enhancing-tumor ROIs and their time-intensity curves](docs/screenshot_main.png)

![Subtraction map (last - first frame) with the enhancing rim grown from the crosshair](docs/screenshot_subtraction.png)

*Synthetic brain phantom (ring-enhancing lesion, superior sagittal sinus), rendered with GNU
Octave; controls look slightly different in MATLAB. Top: frame view after 2x temporal binning
with the SSS control and the tumor ROI. Bottom: subtraction map (frame 60 - frame 1),
threshold preview in green and the enhancing rim grown from the crosshair (necrotic core
excluded).*

## Quick start

```matlab
addpath('/path/to/GRASP_DCE_Viewer')
GraspDCEViewer                          % then File > Open .mat (or the Open button)
GraspDCEViewer('recon.mat')             % opens the load dialog for this file
GraspDCEViewer(vk_all, basis_all)       % variables already in the workspace
GraspDCEViewer(data1)                   % output of ReadData.m ([nx ny t nz])
GraspDCEViewer('recon.mat', 'Voxel', [0.94 0.94 2], 'dt', 3.3, 'NoDialog', true)
GraspDCEViewer('recon.mat', 'BinSize', 2)   % average every 2 frames on load
```

The load dialog lists the variables in the file, detects `vk_all`/`basis_all`, shows the
resulting matrix/K/T and a memory estimate, and asks for the acquisition geometry
(native slice plane, pixel spacing, slice spacing, temporal resolution). The geometry is
remembered between sessions and can be changed later under View > Geometry.

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

## Orientation and coordinates

The "native" array is `data1(:,:,t,z)` (rows, cols, partitions). Choose the native slice
plane (axial / coronal / sagittal); the native plane is then displayed exactly as
`imagesc(data1(:,:,t,z))` shows it, and the two reformats are built from it. If the
anatomy is not radiological (anterior up, patient right on the left; superior up on
coronal/sagittal), use the flips / row-col swap under View > Geometry. ROIs follow
orientation changes.

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
  number, time of the frame, play/pause with fps.
- Bottom: time-intensity curves (all visible ROIs overlaid) and the ROI manager.

Mouse: left = active tool, right-drag = window/level, shift- or middle-drag = pan,
wheel = slice, Ctrl+wheel = zoom, click on the TIC plot = jump to that frame.
Keys: Left/Right = frame -/+1 (Shift: 10), Home/End, PgUp/PgDn, Up/Down = slice,
Space = play, N/C/B/E = tools, [ ] = brush radius, D = 3D toggle, P = threshold preview,
S = subtraction view, G = grow from crosshair, Ctrl+Z / Ctrl+Y = undo / redo. Help > Mouse and keyboard lists everything.

## Segmentation tools

- Navigate: click/drag moves the crosshair; the crosshair voxel's curve is shown dashed in
  the TIC plot (Cursor TIC), useful for locating the superior sagittal sinus.
- Circle: press at the centre, drag to the radius (radius shown in mm, anisotropic voxels
  handled). With 3D on it fills a sphere.
- Brush / Eraser: radius in mm; 2D (disk in the plane of the view) or 3D (sphere).
- Exclusive ROIs (default on): painting removes the voxels from other ROIs, so the label
  map has no overlap.
- Undo / redo: 100 steps, stored as voxel diffs (Ctrl+Z, Ctrl+Y).
- Multiple ROIs: `+ SSS control`, `+ Tumor` (auto-numbered Tumor 1..n), `+ Custom`.
  Rename, recolor, show/hide, delete; the table shows voxels and mL.
- Clear slice / Clear ROI, Fill holes (slice by slice in the plane of the last-used view).

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
5. Peak - baseline (max over time)
6. Peak relative enhancement (%)

Baseline = mean of the chosen pre-contrast frames (default 1-5). Maps 3-6 are computed on
first use (one pass over the data, a few seconds). Every map works with Preview, Otsu,
Keep: slice/3D, Grow from crosshair, All -> ROI and "Constrain brush/circle" (painting only
adds in-window voxels). The line under the buttons shows how many voxels of the active ROI
the current window keeps.

## Temporal binning

**Average every n frames -> Apply binning** (or View > Temporal binning) lowers the temporal
resolution: each frame becomes the mean signal intensity (magnitude) of n consecutive
reconstructed frames, e.g. to compare an 8-spoke reconstruction binned x2 or x4 with the
16- and 32-spoke reconstructions. Everything follows the binned series: views, slider,
curves, maps, the subtraction map (last and first binned frames), NIfTI export (pixdim(4) =
n x dt) and spreadsheets (frame times = bin centres). An incomplete last bin is dropped by
default (reported in the status bar, the JSON sidecar and the spreadsheet Info sheet); the
View menu dialog can keep it instead. n = 1 restores the original series. The current frame
and the baseline frames are mapped to the same time points.

## Time-intensity curves

Mean +/- SD band, median [IQR] band, or mean only; Y as signal, relative enhancement
(%), S/S0 or S-S0 (S0 = mean of the ROI curve over the baseline frames); X in seconds
(if temporal resolution is set) or frames. File > Save TIC plot as PNG.

## Exports

| Button / menu | Output |
|---|---|
| Image -> NIfTI | 4D float32 `.nii.gz`/`.nii` (all frames, streamed frame by frame) or the current frame (3D), plus a JSON sidecar (frame times, source, geometry, normalization) |
| Mask -> NIfTI | uint8 label map (value = ROI number), `<name>_labels.txt` (label, name, RGB, voxels, mL), `<name>_itksnap_labels.txt`; optionally one binary mask per ROI |
| ROI voxel time series -> spreadsheet | `.xlsx`: sheets Info, Summary (per-ROI mean, SD, median, relative enhancement, n per frame) and one sheet per ROI with one row per voxel: Voxel, i, j, k, x_mm, y_mm, z_mm, F001..Fnnn. `.csv`: `_voxels.csv` (all ROIs, ROI_ID/ROI_Name columns), `_summary.csv`, `_info.txt` |
| File > Export enhancement maps | Baseline, Max, PeakEnh, RelEnh, TTP (s), AUC as NIfTI |
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
    gdce_convert_mat2nifti(f(k).name, '', 'Voxel', [0.94 0.94 2], 'dt', 3.3);
end
gdce_convert_mat2nifti('recon_8spoke.mat', '', 'dt', 1.6, 'BinSize', 2);   % binned x2

% data model only
D = GraspData.fromFile('recon.mat');
D.setOrientation('axial', false, [0 0 0]); D.VoxNative = [0.94 0.94 2]; D.dt = 3.3;
D.setBinning(2);                   % optional: average frame pairs
V  = D.canonVolume(40);            % 3D frame 40 (canonical order)
S  = D.canonSubtraction(D.NT, 1);  % last - first
TC = D.timecourses(find(mask));    % voxel x frame matrix for a canonical mask
gdce_export_timeseries(D, {mask}, {'Tumor 1'}, 'tumor.csv', 'BaseFrames', 1:5);

% drive the GUI
app = GraspDCEViewer(D);
app.setFrame(45); app.addTumor();
app.setThrSource('sub'); app.setThreshold(0.2, Inf);
app.setPos([90 50 24]); app.growFromCrosshair();
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
| `gdce_reorient_to_ref.m` | mask import reorientation |
| `gdce_convert_mat2nifti.m` | batch .mat -> NIfTI |
| `gdce_app.m` | launcher / compiler entry point |
| `docs/` | README screenshots |
| `legacy/` | the original ReadData.m, CropImg.m and @TempPCASub (reference for validation) |
| `tests/gdce_selftest.m` | synthetic end-to-end test (`gdce_selftest(true)` also drives the GUI) |
| `tests/gdce_validate_readdata.m` | compares the viewer with ReadData.m on a real file |

## Notes and limitations

- Geometry comes from the dialog; there is no DICOM header import yet. Voxel sizes,
  orientation and frame timing should be checked once per protocol.
- Signal is the magnitude of the reconstruction (normalized to the global maximum when
  ReadData.m normalization is on); no T1 mapping or concentration conversion is done.
- Fill holes is 2D (slice by slice). Grow from crosshair is the connected-component tool
  (6-connectivity, bounded by the max distance).
- Binning averages magnitudes of reconstructed frames; it does not re-reconstruct with
  more spokes per frame.
