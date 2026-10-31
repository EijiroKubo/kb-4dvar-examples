# KB-4DVar: Shallow-water equation example

This directory contains research code for an example of **Koopman-based four-dimensional variational data assimilation (KB-4DVar)** using a shallow-water-equation model. The code combines a Fortran physical simulation with Python code for extended dynamic mode decomposition (EDMD), an approximate adjoint, and iterative initial-state correction.

The repository is provided as **research/example code**, not as a production flood-forecasting system. The original numerical algorithms, configuration, and observations are retained in this release package; they have not been refactored or modified as part of the GitHub packaging work.

## Included files

| File | Role |
| --- | --- |
| `main.py` | Driver for physical simulations and iterative data assimilation |
| `calc_edmd_adj.py` | EDMD approximation, approximate adjoint, and result visualization |
| `calc_sea_field.py` | Python interface to the Fortran shared library |
| `calc_dct.py` | Two-dimensional DCT-II filtering routines |
| `obsdata_controller.py` | Observation data and observation-operator handling |
| `config.yml` | Physical-model and data-assimilation configuration |
| `observe.zip` | Supplied water-depth observation series and observation locations |
| `seafield/CalcSeaFieldMod_Alter.f90` | Fortran shallow-water simulation source |
| `seafield/mesh.csv` | Optional mesh input when `read_initial_file: True` |
| `build.sh` / `build_win.bat` | Linux / Windows build scripts for the shared library |
| `requirements.txt` | Python dependencies |

## Requirements

- Python 3 with NumPy, SciPy, pandas, PyYAML, Matplotlib, and Pillow. See `requirements.txt`.
- A Fortran compiler that supports shared libraries and OpenMP. The supplied scripts use `gfortran`.
- Sufficient memory and CPU time for the configured mesh and assimilation loops.

The examples are intended to be run **from this directory**, because configuration and observation paths are relative to the working directory. Running `main.py` from another directory is supported by its internal working-directory change, but running imports or functions directly may require setting the working directory yourself.

### Linux

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
bash build.sh
python main.py
```

The build produces `seafield/CalcSeaFieldMod.so`.

### Windows

Install Python, the packages in `requirements.txt`, and a compatible `gfortran` toolchain, then run from this directory:

```bat
python -m pip install -r requirements.txt
build_win.bat
python main.py
```

The build produces `seafield/CalcSeaFieldMod.dll`. **Windows-specific note:** `calc_sea_field.py` currently uses a hard-coded dependent-DLL search directory, `C:\TDM-GCC-64\bin`. If your compiler runtime is installed elsewhere, that path must be adapted locally or an equivalent environment made available before running the original code. This was deliberately left unchanged in this packaging-only release.

## Default configuration and DCT filtering

The provided `config.yml` uses a 200 x 200 grid, 80 retained EDMD modes, 10 local iterations and 100 global iterations. These settings may make a full execution computationally demanding.

The observation archive is read automatically from `observe.zip` by `obsdata_controller.py`. Do **not** unpack or rename it for the default run.

The setting below is preserved from the supplied code:

```yaml
DMD:
  dct_threshold: 0.8
```

In `main.py`, this retains the low-frequency 80% range of DCT coefficients **along each spatial axis** and zeroes coefficients outside that range (approximately 20% at the high-frequency end of each axis). The filter is applied to the corrected initial **water-depth** field after each global assimilation loop, before the next physical-model run. It is not a guarantee that all numerical instabilities or NaN values will be eliminated.

The script also has `DIMENSION_STANDARDIZE = False` and `SAVE_ANIM = False` by default. Despite the latter, the original script still saves two final GIF visualizations to `result_anim/`. That behavior has not been changed.

> **Reproducibility note:** This is the supplied research implementation. The README documents what the code actually executes. It does not assert that this package reproduces every figure or numerical value in a published paper. For a publication, cross-check the filtering stage, model settings, and reported experiments against the manuscript.

## Generated files and rebuilding

Generated GIF images are written to `result_anim/`. Compiled `.so`/`.dll` files and Fortran intermediate files are excluded from this source distribution; regenerate them using the supplied build scripts. `.gitignore` prevents typical compiled artifacts, Python caches, and visualizations from being committed.

## Attribution and licensing

The DCT source file acknowledges an external Japanese-language tutorial as a reference. Details are recorded in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md). If any source was copied or substantially adapted, confirm permission and the required notices before public redistribution.

**No license has been assigned here to the original project code.** The copyright holder should select and include a suitable repository license after confirming the ownership and release permissions for the research software and observation data. Publishing a GitHub repository without a license does not automatically authorize others to reuse or modify it.

## Scope of validation

The source files have been checked for Python syntax, the Fortran code can be compiled on Linux with `gfortran`, and a standalone DCT round-trip/filtering check can be performed. A full 100-global-loop assimilation run and Windows execution have **not** been verified as part of the packaging process.
