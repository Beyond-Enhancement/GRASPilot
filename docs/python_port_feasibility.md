# Porting GRASPilot to Python: feasibility

Assessment only; nothing has been ported. Based on the code at v0.1.0+ (November 2026).

## Verdict

**Feasible, with no technical blockers.** GRASPilot uses no MATLAB toolboxes. It has its own
NIfTI and DICOM readers and writers, and its numerical code is plain array arithmetic
(vectorized curve parameters, batched Levenberg-Marquardt, exact exponential convolution,
SVD, Nelder-Mead). All of it has direct numpy/scipy equivalents. The cost is effort and
validation, not missing capability. The GUI is the largest single part and the one that needs
real redesign rather than translation.

## What there is to port

| Part | MATLAB lines | Python counterpart |
|---|---|---|
| Viewer and windows (`GRASPilot.m`: 3 linked views, painting tools, ROI manager, curve plot, parameter, PK and motion windows, ~390 UI controls) | ~6,400, of which roughly 4,000 are UI | PySide6 (Qt) + pyqtgraph (fast image views, overlays, ROI painting); Qt table views; native dark theme |
| Data model (`GRASPilotData.m`: subspace synthesis, orientation/affine, binning, maps, PK storage, fast-open copy) | ~1,700 | numpy class; `scipy.io.loadmat` (v5-v7.2) and `h5py`/`mat73` (v7.3) for GRASP `.mat` files |
| Curve parameters, SSS window, baseline detection | ~280 | numpy, near line-for-line |
| PK: Tofts, extended Tofts, 2CXM, RATIO, VCAI, population AIF, signal-ratio conversion, varpro + LM | ~800 | numpy; `scipy.optimize` where `fminsearch` is used |
| Auto SSS, Auto AIF | ~530 | numpy + `scipy.ndimage` |
| Motion estimation/correction | ~380 | `scipy.ndimage.map_coordinates` (also gives the cubic interpolation Octave lacked) |
| NIfTI and DICOM I/O | ~1,140 | `nibabel` and `pydicom` replace most of it (pydicom is mature and fast); keep the folder-sorting logic |
| Exports (xlsx/csv, label tables, case export) | ~600 | `openpyxl` / `pandas` |
| Self-tests and phantoms | ~2,300 | `pytest`; port the phantoms first, since they are the reference data |

Total: about 14,700 lines of MATLAB, of which about 4,500 are the numerical core.

## Main risks and how to handle them

1. **Indexing and memory order.** MATLAB is 1-based and column-major, numpy is 0-based and
   row-major. This touches every index computation: canonical LPS order, linear indices in ROI
   masks, `ind2sub`, exported i/j/k columns (currently 1-based) and NIfTI affines. Keep one
   explicit convention (e.g. arrays in `(i, j, k, t)` order, `order='F'` only at file
   boundaries) and test it on asymmetric phantoms.
2. **Numerical parity.** Quantile and percentile definitions, Otsu, `interp1` edge handling,
   `filter` initial conditions, LM stopping rules and grid choices can all shift results
   slightly. Mitigation: produce golden outputs with the MATLAB version (phantoms plus one or
   two real anonymized cases: maps, ROI statistics, PK fits, Auto SSS/AIF voxels) and require
   agreement to set tolerances before switching.
3. **GUI behaviour.** Painting with brush/circle/sphere, undo/redo, crosshair-linked views,
   threshold preview, keyboard shortcuts and live curve updates have to be rebuilt, not
   translated. pyqtgraph handles large 4D data interactively; matplotlib would be too slow for
   painting.
4. **File compatibility.** ROI sessions (`*_rois.mat`) and fast-open copies (`-v7.3`) should
   stay readable in both versions. `h5py` reads v7.3. Writing v7.3 that MATLAB accepts needs
   care (the `hdf5storage` package), or switch the session format to NIfTI + JSON.
5. **Performance.** The heavy loops are already vectorized and map to numpy at similar speed.
   Batched PK fitting (thousands of voxels/s in Octave) should match or exceed that. `numba`
   is optional for the per-slice loops.

## What Python would gain

- No MATLAB licence needed. Could be shipped as a pip/conda package or a standalone app
  (PyInstaller).
- A headless library and command-line tool for batch processing of many cases, with no GUI
  needed.
- Mature libraries (pydicom, nibabel), and easier integration with other pipelines (e.g. deep
  learning segmentation, BIDS, XNAT, PACS tools).

## Suggested plan (if and when it is done)

1. **Golden data (about 1 week):** script the MATLAB version to export everything (case export)
   for the phantoms and real cases; freeze it as reference outputs.
2. **Core library, headless (about 4-6 weeks):** data model, I/O, curve parameters, Auto SSS/AIF,
   baseline, PK/RATIO/VCAI, motion, exports; pytest against the golden data; a CLI
   (`graspilot export-case input.mat --rois rois.mat`).
3. **GUI (about 6-10 weeks):** PySide6 + pyqtgraph, reusing the core. Port window by window
   (main viewer, parameters, PK, motion, load and geometry dialogs).
4. **Parallel running (a few weeks of use):** both versions side by side on real cases, then
   switch.

Total for one experienced developer: roughly 3-4 months, with the core library usable after
about 1.5 months. AI-assisted translation shortens the coding but not the validation, which is
the critical path. Porting the core first carries little risk and pays off on its own: batch
processing in Python while the MATLAB GUI keeps working.
