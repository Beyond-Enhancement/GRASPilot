# GRASP DCE Viewer (MATLAB)

3-plane viewer and ROI analysis platform for dynamic GRASP DCE-MRI reconstructions
stored as MATLAB `.mat` files (deep learning subspace recon: `vk_all` + `basis_all`).

- No toolboxes required (own NIfTI writer/reader, own drawing tools). MATLAB R2019b or
  newer recommended; `.xlsx` export needs R2019a+ (`writecell`/`writematrix`), `.csv` works anywhere.
- Runs fully offline. Nothing is written except the files you export.

## Quick start

```matlab
addpath('/path/to/GRASP_DCE_Viewer')
GraspDCEViewer                          % then File > Open .mat (or the Open button)
GraspDCEViewer('recon.mat')             % opens the load dialog for this file
GraspDCEViewer(vk_all, basis_all)       % variables already in the workspace
GraspDCEViewer(data1)                   % output of ReadData.m ([nx ny t nz])
GraspDCEViewer('recon.mat', 'Voxel', [0.94 0.94 2], 'dt', 3.3, 'NoDialog', true)
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
- Multiple ROIs: `+ SSS control`, `+ Tumor` (auto-numbered Tumor 1..n), `+ Custom`.
  Rename, recolor, show/hide, delete; the table shows voxels and mL.
- Clear slice / Clear ROI, Fill holes (slice by slice in the plane of the last-used view).

## Threshold (enhancing tumor only)

Threshold maps:
1. Signal at current frame
2. Current frame - baseline
3. Current frame relative enhancement (%)
4. Peak - baseline (max over time)
5. Peak relative enhancement (%)

Baseline = mean of the chosen pre-contrast frames (default 1-5). Maps 2-5 are computed on
first use (one pass over the data, a few seconds).

Typical workflow: `+ Tumor`, draw a generous circle/sphere around the lesion, choose map
3 or 5, press Otsu (or type Low/High), tick Preview (green = inside the window), then
Keep: 3D (or Keep: slice). Alternatively tick "Constrain brush/circle" so painting only
adds voxels inside the window. The line under the buttons shows how many ROI voxels the
current window keeps.

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

% data model only
D = GraspData.fromFile('recon.mat');
D.setOrientation('axial', false, [0 0 0]); D.VoxNative = [0.94 0.94 2]; D.dt = 3.3;
V  = D.canonVolume(40);            % 3D frame 40 (canonical order)
TC = D.timecourses(find(mask));    % voxel x frame matrix for a canonical mask
gdce_export_timeseries(D, {mask}, {'Tumor 1'}, 'tumor.csv', 'BaseFrames', 1:5);

% drive the GUI
app = GraspDCEViewer(D);
app.setFrame(45); app.addTumor(); app.exportMaskNifti('mask.nii.gz');
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
| `legacy/` | the original ReadData.m, CropImg.m and @TempPCASub (reference for validation) |
| `tests/gdce_selftest.m` | synthetic end-to-end test (`gdce_selftest(true)` also drives the GUI) |
| `tests/gdce_validate_readdata.m` | compares the viewer with ReadData.m on a real file |

## Notes and limitations

- Geometry comes from the dialog; there is no DICOM header import yet. Voxel sizes,
  orientation and frame timing should be checked once per protocol.
- Signal is the magnitude of the reconstruction (normalized to the global maximum when
  ReadData.m normalization is on); no T1 mapping or concentration conversion is done.
- Fill holes is 2D (slice by slice); there is no connected-component tool yet.
