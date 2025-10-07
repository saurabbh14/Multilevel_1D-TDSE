# Multi-level 1D TDSE solver — Documentation (updated)

This document summarizes the repository layout, how to build and run the code (CLI and GUI), where to find important modules, and the main input / output artifacts.

## Quick start — command line

- Build (Nix):
  ```
  nix build
  ```
  Executable: `result/bin/ML-TDSE`

- Run (CLI):
  ```
  ./result/bin/ML-TDSE input.ini
  ```
  Use the provided [`input.ini`](../input.ini) at project root as a template.

## Quick start — GUI

A development PySide6 GUI is provided to edit `input.ini`, build, run simulations, view run logs and plot outputs interactively.

- Install Python deps (recommended inside a venv):
  ```
  python3 -m pip install -r tools/requirements_pyqt.txt
  ```
  See: [`tools/requirements_pyqt.txt`](../tools/requirements_pyqt.txt)

- Launch GUI:
  ```
  python3 tools/pyqt_gui.py
  ```
  GUI script: [`tools/pyqt_gui.py`](../tools/pyqt_gui.py)

- Features:
  - Raw editor + quick fields for common parameters
  - Build (nix build), Run, Stop controls
  - Runtime log viewer; logs are saved to `run_logs/run_log_N.txt` after each run (see [`run_logs/`](../run_logs/))
  - Interactive matplotlib plot (zoom/pan) with live updates (plots a selected output file periodically once it is created)
  - File chooser to select any output file to plot

See [`documentation/GUI_DOC.md`](GUI_DOC.md) for a short GUI guide.

## Project layout (important files)

- Top-level:
  - [`input.ini`](../input.ini) — example input
  - [`README.md`](../README.md) — project overview
  - [`documentation/USAGE.md`](USAGE.md) — input parameter reference and quick examples
  - [`documentation/GUI_DOC.md`](GUI_DOC.md) — GUI quick guide
  - [`tools/pyqt_gui.py`](../tools/pyqt_gui.py) — PySide6 GUI
  - [`tools/requirements_pyqt.txt`](../tools/requirements_pyqt.txt) — GUI Python requirements

- Source:
  - [`src/main.f90`](../src/main.f90) — program entry and initialization
  - [`src/initializer.f90`](../src/initializer.f90)
  - [`src/nuclear_wv.f90`](../src/nuclear_wv.f90)
  - [`src/propagation.f90`](../src/propagation.f90)
  - [`src/pulse_gen/pulse.f90`](../src/pulse_gen/pulse.f90)
  - IO & helpers:
    - [`src/IO_modules/`](../src/IO_modules/) — input parsing, variables, FFTW wrapper, etc.
    - [`src/fft/fftw3.f90`](../src/fft/fftw3.f90)
    - [`src/blas_module/`](../src/blas_module/)

## Inputs

The code uses Fortran namelist-style `input.ini`. See [`documentation/USAGE.md`](USAGE.md) for a detailed description of sections and parameters (grid, masses, time grid, electronic states, vibrational state options, laser pulses, absorber, output directory, parallelization flags).

## Outputs

- All outputs are written to the directory specified by `output_data_dir` in `input.ini`. Typical files:
  - `density_1d.out`, `Pdensity_1d.out`, `avgR_1d.out`, `norm_1d.out`
  - KER & momentum spectra: `KER_spectra_...`, `momt_spectra_...`
  - Nuclear/vibrational data in `output_data/nuclear_wavepacket_data/...`
- GUI live plotting watches the selected output file (relative to `output_data_dir`) and updates the figure periodically after the file appears.

## Parallelization / performance

- FFTW threading can be controlled via `prop_par_FFTW` and `ITP_par_FFTW` namelists.
- OpenMP / threaded BLAS behavior depends on your build and environment variables (see [`src/meson.build`](../src/meson.build) and meson configuration).

## Tests / verification

- Use provided `input.ini` for the H2+ test case.
- After a run, compare produced files in the `output_data_dir` and check saved run logs in [`run_logs/`](../run_logs/).
- See [`documentation/TEST_DOC.md`](TEST_DOC.md) for test notes (if present).

## Where to look to change behavior

- Pulse shapes and IO: [`src/pulse_gen/pulse.f90`](../src/pulse_gen/pulse.f90)
- Potentials / dipoles read logic: [`src/main.f90`](../src/main.f90) (`pot_read`, `trans_dipole_read`)
- Propagation / observables: [`src/propagation.f90`](../src/propagation.f90)
- Input parsing and shared variables: [`src/IO_modules/`](../src/IO_modules/)

If you want, the GUI can be extended to validate parameters, present tooltips, or plot multiple files simultaneously — open an issue or request specific enhancements.
```// filepath: /home/saurabh/Saurabh/TDSE_1D/documentation/DOCUMENTATION.md
# Multi-level 1D TDSE solver — Documentation (updated)

This document summarizes the repository layout, how to build and run the code (CLI and GUI), where to find important modules, and the main input / output artifacts.

## Quick start — command line

- Build (Nix):
  ```
  nix build
  ```
  Executable: `result/bin/ML-TDSE`

- Run (CLI):
  ```
  ./result/bin/ML-TDSE input.ini
  ```
  Use the provided [`input.ini`](../input.ini) at project root as a template.

## Quick start — GUI

A development PySide6 GUI is provided to edit `input.ini`, build, run simulations, view run logs and plot outputs interactively.

- Install Python deps (recommended inside a venv):
  ```
  python3 -m pip install -r tools/requirements_pyqt.txt
  ```
  See: [`tools/requirements_pyqt.txt`](../tools/requirements_pyqt.txt)

- Launch GUI:
  ```
  python3 tools/pyqt_gui.py
  ```
  GUI script: [`tools/pyqt_gui.py`](../tools/pyqt_gui.py)

- Features:
  - Raw editor + quick fields for common parameters
  - Build (nix build), Run, Stop controls
  - Runtime log viewer; logs are saved to `run_logs/run_log_N.txt` after each run (see [`run_logs/`](../run_logs/))
  - Interactive matplotlib plot (zoom/pan) with live updates (plots a selected output file periodically once it is created)
  - File chooser to select any output file to plot

See [`documentation/GUI_DOC.md`](GUI_DOC.md) for a short GUI guide.

## Project layout (important files)

- Top-level:
  - [`input.ini`](../input.ini) — example input
  - [`README.md`](../README.md) — project overview
  - [`documentation/USAGE.md`](USAGE.md) — input parameter reference and quick examples
  - [`documentation/GUI_DOC.md`](GUI_DOC.md) — GUI quick guide
  - [`tools/pyqt_gui.py`](../tools/pyqt_gui.py) — PySide6 GUI
  - [`tools/requirements_pyqt.txt`](../tools/requirements_pyqt.txt) — GUI Python requirements

- Source:
  - [`src/main.f90`](../src/main.f90) — program entry and initialization
  - [`src/initializer.f90`](../src/initializer.f90)
  - [`src/nuclear_wv.f90`](../src/nuclear_wv.f90)
  - [`src/propagation.f90`](../src/propagation.f90)
  - [`src/pulse_gen/pulse.f90`](../src/pulse_gen/pulse.f90)
  - IO & helpers:
    - [`src/IO_modules/`](../src/IO_modules/) — input parsing, variables, FFTW wrapper, etc.
    - [`src/fft/fftw3.f90`](../src/fft/fftw3.f90)
    - [`src/blas_module/`](../src/blas_module/)

## Inputs

The code uses Fortran namelist-style `input.ini`. See [`documentation/USAGE.md`](USAGE.md) for a detailed description of sections and parameters (grid, masses, time grid, electronic states, vibrational state options, laser pulses, absorber, output directory, parallelization flags).

## Outputs

- All outputs are written to the directory specified by `output_data_dir` in `input.ini`. Typical files:
  - `density_1d.out`, `Pdensity_1d.out`, `avgR_1d.out`, `norm_1d.out`
  - KER & momentum spectra: `KER_spectra_...`, `momt_spectra_...`
  - Nuclear/vibrational data in `output_data/nuclear_wavepacket_data/...`
- GUI live plotting watches the selected output file (relative to `output_data_dir`) and updates the figure periodically after the file appears.

## Parallelization / performance

- FFTW threading can be controlled via `prop_par_FFTW` and `ITP_par_FFTW` namelists.
- OpenMP / threaded BLAS behavior depends on your build and environment variables (see [`src/meson.build`](../src/meson.build) and meson configuration).

## Tests / verification

- Use provided `input.ini` for the H2+ test case.
- After a run, compare produced files in the `output_data_dir` and check saved run logs in [`run_logs/`](../run_logs/).
- See [`documentation/TEST_DOC.md`](TEST_DOC.md) for test notes (if present).

## Where to look to change behavior

- Pulse shapes and IO: [`src/pulse_gen/pulse.f90`](../src/pulse_gen/pulse.f90)
- Potentials / dipoles read logic: [`src/main.f90`](../src/main.f90) (`pot_read`, `trans_dipole_read`)
- Propagation / observables: [`src/propagation.f90`](../src/propagation.f90)
- Input parsing and shared variables: [`src/IO_modules/`](../src/IO_modules/)

If you want, the GUI can be extended to validate parameters, present tooltips, or plot multiple files simultaneously — open an issue or request specific enhancements.