# SAMap Environment Installation

This document describes how to install and validate the Python environment used for SAMap analyses in the *Macrostomum lignano* single-cell RNA-seq atlas project.

The environment is maintained separately from the R/Seurat environment used by the main atlas pipeline. It is intended primarily for the SAMap Jupyter notebooks under `notebooks/`.

## Overview

The SAMap environment uses:

- Python 3.12
- SAMap 3 (`sc-samap`)
- `samap-extension`
- Scanpy and the scientific Python stack
- Jupyter/IPython kernel support

The environment has been validated by creating a clean Conda environment and running the project SAMap notebooks, including Sankey plotting with `samap-extension`.

The validated core versions are:

| Package | Version |
|---|---:|
| Python | 3.12 |
| sc-samap | 3.0.1 |
| samap-extension | 0.1.1 |

SAMap 3 is distributed on PyPI as `sc-samap`, but it is imported in Python as `samap`:

```python
from samap import SAMAP
```

In SAMap 3.0.1, the `SAMAP` class is implemented as `samap.core.mapping.SAMAP`.

## Repository Structure

```text
environments/samap/
├── environment.yml
├── environment.lock.yml
├── conda-linux-64.lock.txt
└── setup.sh
```

### `environment.yml`

Human-maintained specification containing the direct dependencies:

```yaml
name: mlig_samap

channels:
  - conda-forge

dependencies:
  - python=3.12
  - pip
  - ipykernel

  - pip:
      - sc-samap[viz]==3.0.1
      - samap-extension==0.1.1
```

Transitive dependencies such as NumPy, SciPy, Pandas, Scanpy, AnnData, hnswlib, UMAP, matplotlib, Plotly, and Seaborn are resolved through package dependency metadata.

### `environment.lock.yml`

Records the complete resolved Conda and pip environment validated with the project notebooks. Use it when reproducibility is more important than obtaining the newest compatible transitive dependencies.

### `conda-linux-64.lock.txt`

Records the exact Conda Linux-64 artifacts in the validated environment. Pip-installed packages are not represented by `conda list --explicit`, so this file does not replace `environment.lock.yml`.

### `setup.sh`

Automated installation and validation script. It creates the SAMap Conda environment from the validated environment specification and performs package and import checks.

## Recommended Installation

Clone the repository and move to its root directory:

```bash
git clone <repository-url>
cd Mlig_scRNASeq_atlas
```

Run:

```bash
./environments/samap/setup.sh
```

By default, the installer creates `mlig_samap`. Activate it with:

```bash
conda activate mlig_samap
```

Confirm the interpreter:

```bash
which python
python --version
```

The interpreter should belong to the new environment and Python should be version 3.12.

## Installing Under a Different Environment Name

```bash
ENV_NAME=mlig_samap_test \
  ./environments/samap/setup.sh

conda activate mlig_samap_test
```

This is useful for clean-room reproducibility testing without modifying an existing environment.

## Manual Installation

For the validated locked environment:

```bash
conda env create \
  -n mlig_samap \
  -f environments/samap/environment.lock.yml
```

Alternatively, use the smaller human-maintained specification:

```bash
conda env create \
  -f environments/samap/environment.yml
```

For routine reproduction of project analyses, the validated lock file is preferred.

## Validate the Installation

```bash
conda activate mlig_samap
python -m pip check
```

A successful environment should report `No broken requirements found.`

Check key package versions:

```bash
python - <<'PY'
from importlib.metadata import version

for package in [
    "sc-samap",
    "samap-extension",
    "scanpy",
    "anndata",
    "numpy",
    "pandas",
    "scipy",
]:
    print(f"{package:20s} {version(package)}")
PY
```

Test project imports:

```bash
python - <<'PY'
from samap import SAMAP
from samap.analysis import (
    get_mapping_scores,
    GenePairFinder,
    CellTypeTriangles,
    sankey_plot,
)
from samap_extension import plotting as splt

print("SAMAP:", SAMAP)
print("SAMAP module:", SAMAP.__module__)
print("samap-extension:", splt.__file__)
print("SAMap imports successful.")
PY
```

For SAMap 3.0.1, `SAMAP module` should be `samap.core.mapping`.

## Configure a Jupyter Kernel

Activate the environment first:

```bash
conda activate mlig_samap
```

Register it:

```bash
python -m ipykernel install \
  --user \
  --name mlig_samap \
  --display-name "M. lignano SAMap"
```

## Verify the Jupyter Interpreter

The Jupyter kernel display name is only a label. It does not guarantee which Python executable the kernel uses.

Inside a notebook:

```python
import sys
print(sys.executable)
```

It should point to the intended environment, for example:

```text
.../conda/envs/mlig_samap/bin/python
```

Check SAMap directly:

```python
import samap
from importlib.metadata import version

print(samap.__file__)
print(version("sc-samap"))
```

The expected distribution is `sc-samap 3.0.1`.

## Important: SAMap 3 vs Older SAMap Releases

This project uses `sc-samap 3.0.1` and expects:

```python
from samap import SAMAP
```

The resulting class is `samap.core.mapping.SAMAP`.

An older environment may instead contain the older `samap` PyPI package. A symptom of accidentally using it is:

```text
ImportError: cannot import name 'SAMAP' from 'samap'
```

If this occurs:

```python
import sys
import samap

print(sys.executable)
print(samap.__file__)
```

If the path points to an old Conda environment, the Jupyter kernel is using the wrong interpreter. Do not change the project notebooks to work around this error; select or recreate the correct kernel.

## Diagnosing Jupyter Kernel Problems

List kernels:

```bash
jupyter kernelspec list
```

Inspect a kernel:

```bash
cat ~/.local/share/jupyter/kernels/<kernel-name>/kernel.json
```

The first `argv` entry identifies the Python executable actually used by the kernel.

If a kernel points to the wrong environment:

```bash
jupyter kernelspec remove <kernel-name>
conda activate mlig_samap
python -m ipykernel install \
  --user \
  --name mlig_samap \
  --display-name "M. lignano SAMap"
```

Restart the notebook kernel afterward.

## `samap-extension`

The validated version is `samap-extension 0.1.1`.

Version 0.1.1 includes a compatibility fix for Sankey plotting with modern Pandas/NumPy arrays. The implementation explicitly creates a writable NumPy array:

```python
A = mat_df.to_numpy(dtype=float, copy=True)
```

This prevents:

```text
ValueError: assignment destination is read-only
```

when `sankey_plot_3()` modifies the array during preprocessing.

For reproducible atlas analyses, use the published package rather than an editable development checkout.

Verify its location with:

```python
import samap_extension
print(samap_extension.__file__)
```

It should resolve inside the Conda environment's `site-packages` directory.

## Clean-Room Validation

The environment was validated by:

1. Creating a new Conda environment independently of the original development environment.
2. Installing SAMap 3 and `samap-extension` from published packages.
3. Running `pip check`.
4. Confirming that `SAMAP` resolves to `samap.core.mapping.SAMAP`.
5. Running the actual project SAMap notebook.
6. Testing `samap-extension` Sankey plotting.
7. Recreating the environment from the project installation infrastructure and rerunning the notebook.

This validates the actual analysis workflow, not only package importability.

## Updating the Environment

Do not manually edit `environment.lock.yml` to update individual packages.

For an intentional update:

1. Modify direct requirements in `environment.yml`.
2. Create a fresh test environment.
3. Run `python -m pip check`.
4. Run the project SAMap notebooks.
5. Verify all analyses and plots, including Sankey plots.
6. Regenerate lock files only after validation.

Regenerate the full environment lock:

```bash
conda env export \
  > environments/samap/environment.lock.yml

sed -i '/^prefix:/d' \
  environments/samap/environment.lock.yml
```

Regenerate the exact Conda Linux-64 artifact list:

```bash
conda list --explicit \
  > environments/samap/conda-linux-64.lock.txt
```

Commit the updated human-maintained specification and lock files together.

## Troubleshooting

### `ImportError: cannot import name 'SAMAP' from 'samap'`

Check:

```python
import sys
import samap
print(sys.executable)
print(samap.__file__)
```

The most likely cause is a Jupyter kernel pointing to an older SAMap environment. The project expects `sc-samap==3.0.1`.

### Jupyter shows the correct kernel name but SAMap comes from another environment

Check:

```python
import sys
print(sys.executable)
```

Then inspect the corresponding `kernel.json`. Kernel display names do not determine the interpreter path.

### `ValueError: assignment destination is read-only`

Check:

```bash
python -m pip show samap-extension
```

The project requires version 0.1.1 or later containing the Sankey writable-array fix.

### `pip check` reports errors involving `labcore`

`labcore` is not required by the SAMap notebooks in this project and is not part of the reproducible SAMap environment. If these errors appear, verify that the notebook is running the intended clean environment rather than an older development environment containing an editable `labcore` installation.

## Reproducibility Summary

For routine use:

```bash
./environments/samap/setup.sh
conda activate mlig_samap
```

For development, modify:

```text
environments/samap/environment.yml
```

For exact records of the validated environment, retain:

```text
environments/samap/environment.lock.yml
environments/samap/conda-linux-64.lock.txt
```

The lock files should only be regenerated after the updated environment has successfully run the project SAMap notebooks.

