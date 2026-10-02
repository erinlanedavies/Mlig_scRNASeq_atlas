# Biowulf installation and reproducibility

## Purpose

This document describes the reproducible software environment used by
the **Mlig_scRNASeq_atlas** project on NIH Biowulf. It is intended both
as an installation guide and as a troubleshooting record.

The project uses three distinct layers of software management:

1.  **Biowulf environment modules** provide R, the compiler toolchain,
    native system libraries, Pandoc, and TeX.
2.  **renv** restores the R package environment recorded in `renv.lock`.
3.  **Conda** provides only the Python runtime and Python packages
    required by R/Python interoperability, principally through
    `reticulate` and `zellkonverter`.

Keeping these responsibilities separate is important. Several
installation failures encountered while making the project reproducible
were caused by allowing settings from one layer---especially Conda or
global compiler flags---to interfere with another.

------------------------------------------------------------------------

## 1. Repository files involved in environment management

The relevant files are under `biowulf/` plus the standard renv files in
the project root:

``` text
Mlig_scRNASeq_atlas/
├── biowulf/
│   ├── setup.sh
│   ├── setup_conda.sh
│   ├── setup_renv.R
│   ├── run_r.sh
│   ├── render_rmd.R
│   ├── environment.yml
│   ├── environment.lock.yml
│   └── conda-linux-64.lock.txt
├── renv/
│   └── activate.R
├── renv.lock
└── renv/settings.json
```

### `biowulf/setup.sh`

This is the top-level installer. A new user should normally run this
script once after cloning the repository.

Its responsibilities are to:

-   determine the project root;
-   establish the default or user-overridden environment locations;
-   create the Python environment if it does not already exist;
-   load the required Biowulf modules;
-   configure the renv library and cache locations;
-   restore the R environment from `renv.lock`;
-   validate the resulting R and Python environments.

It must **not** update or snapshot the R environment.

### `biowulf/setup_conda.sh`

Creates the Python environment from the explicit Linux lock file:

``` text
biowulf/conda-linux-64.lock.txt
```

The normal environment location is:

``` text
/data/$USER/conda_envs/mlig_scrna_runtime
```

The location can be overridden with `ENV_PREFIX`, which is useful for
clean-room testing.

### `biowulf/setup_renv.R`

Restores the R environment from the committed `renv.lock` using
`renv::restore()` and then checks `renv::status()`.

This script deliberately does **not** call:

``` r
renv::update()
renv::snapshot()
```

The lockfile is treated as the authoritative description of the R
environment during reproduction.

### `biowulf/run_r.sh`

Runtime launcher for `.R` and `.Rmd` work after installation is
complete.

It loads the required modules, establishes the renv and Python runtime
paths, and launches the requested analysis. It is intentionally **not an
installer** and does not configure package-compilation flags.

### `biowulf/render_rmd.R`

R Markdown rendering helper. It explicitly activates renv because
`run_r.sh` invokes R with `--vanilla`, which bypasses the project's
`.Rprofile`.

### Conda environment files

The three Conda files have different roles:

-   `environment.yml` --- small, human-maintained list of direct Python
    requirements.
-   `environment.lock.yml` --- fully resolved Conda environment, useful
    for auditing exact versions/builds.
-   `conda-linux-64.lock.txt` --- explicit Linux-64 package artifact
    list used for the strongest reproduction on Biowulf.

The explicit lock file is the installation source used by
`setup_conda.sh`.

------------------------------------------------------------------------

## 2. Validated software stack

The environment was validated with the following Biowulf modules:

``` text
hdf5/1.12.2
netcdf/4.9.0_gcc-11.3.0
pandoc/2.18
tex/2024
pcre2/10.40_gcc-11.3.0
gcc/11.3.0
openmpi/5.0.5/gcc-11.3.0
libtiff/4.6.0_gcc-11.3.0
R/4.5.2
libwebp/1.6.0-gcc-11.3.0
```

The validated compiler was:

``` text
/usr/local/GCC/11.3.0/bin/gcc
/usr/local/GCC/11.3.0/bin/g++
```

with GCC/G++ 11.3.0.

The expected `pkg-config` executable is:

``` text
/usr/bin/pkg-config
```

These exact versions matter because Biowulf modules can change over
time. If a previously reproducible installation starts failing after an
HPC module update, compare the currently loaded module versions and
compiler/library paths against this list before changing the project
lockfiles.

------------------------------------------------------------------------

## 3. Python environment

The Python environment is deliberately minimal. Conda is used for Python
and Python packages only; it is **not** used as the native build
environment for R packages.

The human-maintained requirements are conceptually:

``` yaml
name: mlig_scrna_runtime

channels:
  - conda-forge

dependencies:
  - python=3.10
  - anndata
  - h5py
  - numpy
  - pandas
  - scipy
```

A validated environment contained:

``` text
Python   3.10.21
anndata  0.11.4
h5py     3.16.0
numpy    2.2.6
pandas   2.3.3
scipy    1.15.2
```

For exact reproduction, `setup_conda.sh` creates the environment from
`conda-linux-64.lock.txt`, rather than solving `environment.yml` again.

### Important: do not activate Conda for the R runtime

The launcher does not run `conda activate` before starting R. Instead it
sets:

``` bash
RETICULATE_PYTHON=/data/$USER/conda_envs/mlig_scrna_runtime/bin/python
```

and lets `reticulate` use that interpreter directly.

This is intentional. Activating Conda can alter `PATH`,
`LD_LIBRARY_PATH`, `PKG_CONFIG_PATH`, and related variables. That can
cause R packages to compile or load against Conda-provided native
libraries instead of the Biowulf/system libraries with which the R
environment was validated.

------------------------------------------------------------------------

## 4. R environment

The project uses R 4.5.2 and renv.

The project uses a complete, committed `renv.lock`. Examples of
important package versions in the validated environment include:

``` text
Seurat               5.5.0
EnhancedVolcano      1.28.2
DoubletFinder        2.0.2
harmony              2.0.2
leidenbase           0.1.36
scater               1.38.1
sctransform           0.4.3
igraph                2.3.0
ragg                  1.5.2
haven                 2.5.5
selectr                0.6-0
presto                 1.1.0
```

Relevant renv settings include:

``` text
snapshot.type = all
lockfile.sanitize = true
use.cache = true
```

The validated lockfile includes `selectr 0.6-0` and `presto 1.1.0`.
`presto` is recorded as a URL source pinned to the immutable GitHub
commit `b5df6ee6097eb62522f2e557aa93ac21bab2d05f`, so restoration does
not need to resolve the repository's moving HEAD through the GitHub API.

`snapshot.type = all` is intentional. The objective is to preserve the
complete known-good project environment rather than rely on source-code
dependency discovery to decide which packages belong in the lockfile.

### renv storage

The normal persistent locations are:

``` text
RENV_PATHS_LIBRARY_ROOT=/vf/users/$USER/renv_libs
RENV_PATHS_CACHE=/vf/users/$USER/renv_cache
```

These can be overridden for testing.

------------------------------------------------------------------------

## 5. Normal installation

After cloning the repository, installation should normally be a single
command:

``` bash
git clone <repository-url>
cd Mlig_scRNASeq_atlas
./biowulf/setup.sh
```

The installer should:

1.  create the locked Python environment;
2.  load the validated Biowulf modules;
3.  configure renv paths;
4.  restore the R environment from `renv.lock`;
5.  validate the R packages;
6.  validate Python imports;
7.  validate Python access through `reticulate`.

After installation, analyses can be run with:

``` bash
./biowulf/run_r.sh scripts/KnitReports.R
```

or, for an individual R Markdown report:

``` bash
./biowulf/run_r.sh Rmarkdown/example.Rmd
```

### Validated clean-room result

The current installer was validated on **October 2, 2026** from a fresh
Git clone using a new Conda prefix, a new renv library, and a new empty
renv cache. The first restore completed successfully:

``` text
Successfully installed 315 packages in 1700 seconds.

=== Checking renv status ===
No issues found -- the project is in a consistent state.
```

The installer then validated the expected R packages, Python packages,
and Python access through `reticulate`. Finally,

``` bash
./biowulf/run_r.sh scripts/KnitReports.R
```

completed as expected using the same clean-room environment.

This end-to-end test is the reference validation for the current
installation workflow.

### Installation logs

`setup.sh` writes its output to a timestamped log under `logs/`.
Preserve the log when diagnosing a failed installation because the
earliest compiler or linker error is often more informative than the
final list of packages that failed.

------------------------------------------------------------------------

## 6. Why R is launched with `--vanilla`

The launcher uses `Rscript --vanilla` to reduce dependence on
user-specific R startup files.

A consequence is that `.Rprofile` is not read, so renv is not activated
automatically. For `.R` files, `run_r.sh` therefore explicitly performs
the equivalent of:

``` r
renv::load(project = project_root)
```

before sourcing the target script.

For `.Rmd` files, `render_rmd.R` explicitly sources:

``` text
renv/activate.R
```

before calling `rmarkdown::render()`.

This avoids a situation where a report silently uses packages from the
user's global R library instead of the project renv library.

------------------------------------------------------------------------

## 7. Compiler handling: GCC version versus C++ standard

One of the most important installation issues found during
reproducibility testing was the distinction between the **compiler
version** and the **C++ language standard**.

The project uses one compiler toolchain:

``` text
GCC/G++ 11.3.0
```

However, individual R packages can require different C++ standards. For
example, one package may compile as C++14 while another requires C++17.

This is normal. GCC 11.3.0 supports these standards.

The project therefore does **not** globally set variables such as:

``` bash
CXX11
CXX11STD
CXX11FLAGS
PKG_CXXFLAGS
```

and in particular does not globally force:

``` text
-std=gnu++14
```

Each R package and R's own build system must be allowed to request the
language standard it requires.

### Failures encountered during clean-room testing

An earlier launcher globally exported settings similar to:

``` bash
export CXX11="g++"
export CXX11STD="-std=gnu++14"
export CXX11FLAGS="-O2 -march=haswell -mtune=generic"
export PKG_CXXFLAGS="-std=gnu++14"
```

These settings were introduced while investigating an older `presto`
build that forced C++11 even though the locked `RcppArmadillo`/Armadillo
headers required C++14.

The global workaround created a broader reproducibility problem.
Packages requiring C++17 then failed. In particular:

-   `BiocNeighbors` uses C++17 features such as `if constexpr`;
-   `RcppTOML` 0.2.3 explicitly requests C++17;
-   the global `PKG_CXXFLAGS=-std=gnu++14` was appended after
    package-specific C++17 flags.

The resulting compiler command could contain:

``` text
g++ -std=gnu++17 ... -std=gnu++14 ...
```

Because the last `-std` option takes precedence, the package was
actually compiled as C++14.

The global C++ overrides were therefore removed permanently. Two
remaining clean-room failures were then traced to outdated package
versions and fixed at the package level:

-   **`selectr 0.5-1` → `0.6-0`**. The old archived version depended on
    `stringr`, and renv's dependency planning installed `selectr` before
    `stringr` during the first clean restore. `selectr 0.6-0` no longer
    has that dependency, eliminating the failure without a two-pass
    restore.
-   **`presto 1.0.0` → `1.1.0`**. The old pinned source commit
    (`7636b3d0465c468c35853f82f1717d3a64b3c8f6`) explicitly contained
    `CXX_STD = CXX11` in both `src/Makevars` and `src/Makevars.win`.
    With current `RcppArmadillo`, compilation failed because Armadillo
    requires C++14. The lockfile now pins `presto 1.1.0` to immutable
    commit `b5df6ee6097eb62522f2e557aa93ac21bab2d05f`, which removes the
    obsolete C++11 build constraint.

After these package updates, a completely fresh first-pass restore
succeeded with no failed packages. No global or package-specific
compiler-standard workaround is required.

### Troubleshooting implication

If compilation failures appear after a future Biowulf update, **do not
immediately introduce a global C++ standard**.

First inspect the actual compiler command. A package requesting C++17
should be allowed to compile with something similar to:

``` text
g++ -std=gnu++17 ...
```

A package requiring C++14 may use:

``` text
g++ -std=gnu++14 ...
```

Both can use the same GCC 11.3.0 compiler.

If one specific package requires special handling, solve it at the
package level rather than changing the C++ standard globally for the
entire renv restore.

------------------------------------------------------------------------

## 8. Native-library issue encountered with `ragg`

Another reproducibility problem involved `ragg` 1.5.2.

An older cached build of `ragg` had been linked against libraries that
were no longer available in the active runtime, including SONAMEs such
as:

``` text
libjpeg.so.8
libbz2.so.1.0
```

The correct solution was **not** to create manual compatibility symlinks
and not to add Conda library directories to `LD_LIBRARY_PATH`.

Instead, `ragg` was rebuilt under the correct Biowulf module
environment:

``` r
renv::install("ragg@1.5.2", rebuild = TRUE)
```

The rebuilt `ragg.so` linked against available Biowulf/system libraries,
including:

``` text
libpng16.so.16
libtiff.so.5
libjpeg.so.62
libwebp.so.7
libwebpmux.so.3
libbz2.so.1
```

with `libwebp` and `libwebpmux` supplied by the Biowulf
`libwebp/1.6.0-gcc-11.3.0` module.

### Troubleshooting implication

If an R shared object starts failing with messages such as:

``` text
error while loading shared libraries
...
not found
```

inspect its linkage before modifying environment variables. For example:

``` bash
ldd /path/to/package/libs/package.so
```

If the binary was built against obsolete libraries, rebuilding the
package under the current validated module stack is safer than
fabricating SONAME symlinks.

------------------------------------------------------------------------

## 9. `pkg-config`, `PKG_CONFIG_PATH`, and native libraries

The validated environment uses:

``` text
/usr/bin/pkg-config
```

The Biowulf modules establish the paths required for native dependencies
such as `libwebp`.

The setup/runtime scripts verify that `pkg-config` can locate:

``` text
libwebp
libwebpmux
```

For example:

``` bash
pkg-config --exists libwebp
pkg-config --exists libwebpmux
```

Do not replace `/usr/bin/pkg-config` with a Conda `pkg-config` during R
package compilation unless there is a specific, understood reason to do
so.

Likewise, preserve the `PKG_CONFIG_PATH` and `LD_LIBRARY_PATH`
established by the Biowulf module system rather than prepending Conda
library directories.

------------------------------------------------------------------------

## 10. Do not use Conda as the R native-library environment

A central design decision is:

> Conda provides Python; Biowulf modules provide the native toolchain
> and libraries used by R.

This separation avoids accidentally compiling an R package against one
set of native libraries and then loading it under another.

Therefore the setup and launcher should not generally do things such as:

``` bash
export LD_LIBRARY_PATH="$ENV_PREFIX/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$ENV_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
export CPATH="$ENV_PREFIX/include:$CPATH"
export LIBRARY_PATH="$ENV_PREFIX/lib:$LIBRARY_PATH"
```

Such changes can make an installation appear to work while creating
binaries that depend on Conda libraries not present in the intended R
runtime.

------------------------------------------------------------------------

## 11. renv staging directories

An earlier version of `run_r.sh` attempted to manage `renv/staging` by
replacing it with a symlink to a `/data` location.

This caused the launcher to fail after a fresh `renv::restore()` because
renv had created a real `renv/staging` directory. The launcher then
reported that the path existed and was not a symlink.

This exposed another important separation-of-responsibilities issue:

-   package installation/staging belongs to the installation process;
-   running already installed analyses should not manipulate renv
    installation staging directories.

The runtime launcher therefore no longer creates, deletes, redirects, or
validates `renv/staging`.

If future installation performance or filesystem limits require staging
redirection, implement that in the installation workflow, not in
`run_r.sh`.

------------------------------------------------------------------------

## 12. Clean-room reproducibility testing

The environment was tested using a fresh clone, a new Conda prefix, a
new renv library, and a new empty renv cache. This is substantially
stronger than reinstalling into an existing project environment because
it detects dependencies that were previously being satisfied
accidentally by cached packages or user libraries.

The scripts support environment-variable overrides for this purpose.

Example:

``` bash
export ENV_PREFIX="/data/$USER/conda_envs/mlig_scrna_runtime_fresh_test"
export RENV_PATHS_LIBRARY_ROOT="/vf/users/$USER/mlig_fresh_test/renv_libs"
export RENV_PATHS_CACHE="/vf/users/$USER/mlig_fresh_test/renv_cache"

./biowulf/setup.sh
```

Then test the actual runtime using the same exported variables:

``` bash
./biowulf/run_r.sh scripts/KnitReports.R
```

The launcher prints its selected paths. For a clean-room test, verify
that the output points to the fresh-test environments rather than the
production locations.

A rigorous clean-room test should satisfy all of the following:

-   repository was freshly cloned from Git;
-   Conda prefix did not previously exist;
-   renv project library did not previously exist;
-   renv cache was empty/new;
-   no production renv library is being used;
-   no production Conda environment is being used;
-   `renv.lock` is not modified;
-   `renv::snapshot()` is not run;
-   missing packages are not manually installed during the test;
-   `renv::restore()` is allowed to finish before diagnosing failures.

For the current validated environment, the clean-room restore should
finish with zero failed packages, `renv::status()` should report no
issues, and `run_r.sh scripts/KnitReports.R` should execute successfully
using the fresh-test paths.

------------------------------------------------------------------------

## 13. Important lesson from interrupted `renv::restore()`

During one clean-room test, `renv::restore()` was interrupted after some
package compilation failures appeared. A subsequent `renv::status()`
naturally reported a very large number of packages as missing.

That did **not** mean that all of those packages had independently
failed. They had simply never been reached because the restore was
interrupted.

Therefore:

> Do not diagnose the completeness of an environment from
> `renv::status()` after intentionally interrupting `renv::restore()`.

Allow restore to finish. Then diagnose the packages reported in the
final installation error. Downstream failures often originate from only
one or a few root compilation failures.

This distinction was important when diagnosing the C++14/C++17 issue:
many packages appeared as failed, but the useful root errors were
concentrated in a small number of packages such as `BiocNeighbors`,
`RcppTOML`, and `beachmat`, with additional dependency failures
downstream.

------------------------------------------------------------------------

## 14. Expected runtime validation

A validated runtime should show approximately:

``` text
R version: R 4.5.2
Compiler: GCC/G++ 11.3.0
pkg-config: /usr/bin/pkg-config
Python: 3.10.21
```

Important R packages should include the locked versions listed above.

Important Python imports should succeed:

``` python
import anndata
import h5py
import numpy
import pandas
import scipy
```

`reticulate::py_config()` should resolve to the project's Conda Python,
not a system Python or another Conda environment.

The production environment previously validated successfully with:

``` text
Seurat 5.5.0
EnhancedVolcano 1.28.2
DoubletFinder 2.0.2
harmony 2.0.2
leidenbase 0.1.36
scater 1.38.1
sctransform 0.4.3
igraph 2.3.0
ragg 1.5.2
haven 2.5.5
```

A namespace warning involving `S4Arrays::makeNindexFromArrayViewport`
and `DelayedArray::makeNindexFromArrayViewport` was observed during
validation. It did not prevent successful runtime validation. Treat it
as a warning unless future package changes make it associated with an
actual failure.

------------------------------------------------------------------------

## 15. Diagnosing failures after Biowulf module changes

Because the HPC environment is externally maintained, future module
updates are a likely source of reproducibility failures. When an
installation that previously worked begins failing, investigate the HPC
layer before changing `renv.lock` or the Conda lock.

### Step 1: inspect loaded modules

``` bash
module list
```

Compare against the validated stack documented above.

### Step 2: inspect R

``` bash
which R
which Rscript
R --version
```

Expected major setup:

``` text
R 4.5.2
```

### Step 3: inspect compiler selection

``` bash
which gcc
which g++
gcc --version
g++ --version
```

Validated toolchain:

``` text
GCC/G++ 11.3.0
```

### Step 4: inspect `pkg-config`

``` bash
which pkg-config
pkg-config --version
pkg-config --modversion libwebp
pkg-config --modversion libwebpmux
```

The validated executable is:

``` text
/usr/bin/pkg-config
```

### Step 5: inspect relevant environment variables

``` bash
printf '%s\n' "$PATH"
printf '%s\n' "${LD_LIBRARY_PATH:-}"
printf '%s\n' "${PKG_CONFIG_PATH:-}"
printf '%s\n' "${CPATH:-}"
printf '%s\n' "${LIBRARY_PATH:-}"
```

Look especially for unexpected Conda paths contaminating the R native
build environment.

### Step 6: check for global compiler overrides

``` bash
env | grep -E '^(CC|CXX|CXX11|CXX14|CXX17|CXX20|PKG_CXXFLAGS|PKG_CPPFLAGS|PKG_LIBS)='
```

Unexpected global `-std=` flags are suspicious. A package's own C++
requirement should normally control its language standard.

### Step 7: inspect the first/root package compilation errors

Do not focus only on the final list of failed packages. Search the build
log for the earliest real compiler/linker errors, for example:

``` bash
grep -nE 'ERROR:|error:|not found|requires C\\+\\+' install.log
```

A long list of failed packages may be downstream of one native
dependency.

### Step 8: inspect shared-library linkage

For packages that install but fail to load:

``` bash
ldd /path/to/package/libs/package.so
```

Look for:

``` text
not found
```

If the package was built against obsolete module libraries, rebuild it
under the intended module stack rather than immediately adding
compatibility symlinks.

### Step 9: only then consider changing lockfiles

A changed compiler, system library, or module does not automatically
mean that the R or Python package versions should change.

Modify `renv.lock` or the Conda locks only when the project
intentionally adopts a new software environment and that new environment
has been validated.

------------------------------------------------------------------------

## 16. Useful diagnostic commands

### renv status

``` bash
Rscript --vanilla -e '
renv::load()
renv::status()
'
```

A correctly restored environment should report that the project is
synchronized / has no issues.

### R library paths

``` bash
Rscript --vanilla -e '
renv::load()
print(.libPaths())
'
```

### reticulate Python

``` bash
Rscript --vanilla -e '
renv::load()
print(reticulate::py_config())
'
```

### Python package versions

``` bash
"${ENV_PREFIX:-/data/$USER/conda_envs/mlig_scrna_runtime}/bin/python" - <<'PY'
import anndata, h5py, numpy, pandas, scipy
print("anndata", anndata.__version__)
print("h5py", h5py.__version__)
print("numpy", numpy.__version__)
print("pandas", pandas.__version__)
print("scipy", scipy.__version__)
PY
```

### Shell-script syntax

``` bash
bash -n biowulf/setup.sh
bash -n biowulf/setup_conda.sh
bash -n biowulf/run_r.sh
```

### R setup-script syntax

``` bash
Rscript --vanilla -e '
parse(file = "biowulf/setup_renv.R")
cat("setup_renv.R syntax OK\n")
'
```

------------------------------------------------------------------------

## 17. Updating the environment intentionally

Reproduction and environment maintenance are different operations.

The installation scripts are designed for **reproduction** and therefore
must not automatically update package versions or rewrite lockfiles.

When intentionally updating the project environment:

1.  work in a controlled development environment;
2.  change packages deliberately;
3.  run the project's analyses/tests;
4.  verify native library linkage where appropriate;
5.  validate R/Python interoperability;
6.  update the appropriate lockfile only after validation;
7.  perform a new clean-room installation from the updated locks;
8.  document any changes to required Biowulf modules in this file.

For R, the committed `renv.lock` should represent a known-good complete
environment. For Python, regenerate the resolved/explicit lock only
after the desired environment has been validated.

Do not use the reproduction installer itself to update the environment.

------------------------------------------------------------------------

## 18. Summary of problems already encountered

The following table provides a quick reference for future
troubleshooting.

  -----------------------------------------------------------------------------
  Symptom                 Root cause found              Correct response
  ----------------------- ----------------------------- -----------------------
  `RcppTOML` says C++17   Global                        Remove global C++
  is required             `PKG_CXXFLAGS=-std=gnu++14`   standard override
                          overrode package C++17        
                          request                       

  `BiocNeighbors`         Package requires C++17 but    Allow package/R build
  compilation fails       build was forced to C++14     system to select C++17
  around `if constexpr`                                 

  Many                    Downstream dependency         Diagnose earliest/root
  Seurat/Bioconductor     failures                      compilation errors
  packages fail after a                                 first
  few compiler errors                                   

  `selectr 0.5-1` fails   Historical dependency         Use the validated
  because `stringr` is    metadata / installation       `selectr 0.6-0`
  unavailable during      ordering for the archived     lockfile entry; do not
  first restore           package                       add a two-pass restore
                                                        workaround

  `presto 1.0.0` compiles Old pinned `presto` source    Use the validated
  with C++11 and fails in explicitly forced             `presto 1.1.0` commit
  `RcppArmadillo`         `CXX_STD = CXX11`             pinned in `renv.lock`;
                                                        do not globally force
                                                        C++14

  `ragg.so` cannot find   Cached binary linked against  Rebuild `ragg` under
  old JPEG/bzip2 SONAMEs  obsolete/incompatible         validated Biowulf
                          libraries                     modules

  Temptation to add Conda Mixing Python environment     Keep Conda limited to
  `lib/` to               with R native build/runtime   Python; use
  `LD_LIBRARY_PATH`                                     Biowulf/system native
                                                        libraries

  `run_r.sh` refuses      Runtime launcher was trying   Do not manage
  because `renv/staging`  to manage installation        `renv/staging` in
  is a real directory     staging                       runtime launcher

  Huge number of packages `renv::restore()` was stopped Let restore finish
  appear missing after    before installing them        before interpreting
  interrupted restore                                   final failures

  `Rscript --vanilla`     `--vanilla` skips `.Rprofile` Explicitly call
  does not use project                                  `renv::load()` / source
  renv automatically                                    `renv/activate.R`

  Wrong Python selected   User/project startup          Set `RETICULATE_PYTHON`
  by reticulate           configuration can override    explicitly in
                          interpreter                   launcher/setup

  Conda environment works Conda activation changes      Do not activate Conda
  but R                   native-library/toolchain      before running/building
  compilation/linking     paths                         R packages
  becomes unstable                                      
  -----------------------------------------------------------------------------

------------------------------------------------------------------------

## 19. Design principles to preserve

Future modifications to the installation should preserve the following
rules unless there is a specific, documented reason to change them:

1.  **Biowulf modules own the R compiler/native-library environment.**
2.  **Conda owns only the Python runtime and Python packages.**
3.  **renv owns R package versions.**
4.  **Do not activate Conda before running R.**
5.  **Do not prepend Conda native-library paths to the R environment.**
6.  **Do not globally force a C++ standard.**
7.  **Do not modify `renv.lock` during reproduction.**
8.  **Do not call `renv::snapshot()` during a clean install.**
9.  **Keep installation (`setup.sh`) separate from execution
    (`run_r.sh`).**
10. **Test important environment changes with a fresh clone, fresh
    library, and fresh cache.**
11. **When HPC modules change, diagnose the
    module/compiler/native-library layer before changing package
    locks.**

These principles are the result of actual failures encountered while
reproducing the project and are therefore part of the project's
reproducibility strategy, not merely stylistic preferences.

------------------------------------------------------------------------

## 20. Recommended record when the environment is revalidated

Whenever the project is successfully validated after a significant
Biowulf/module change, record at least:

``` bash
module list
R --version
gcc --version
g++ --version
which pkg-config
pkg-config --modversion libwebp
pkg-config --modversion libwebpmux
```

and from R:

``` r
renv::status()
.libPaths()
reticulate::py_config()
```

Also record the versions of the key R and Python packages. Updating this
document with the date and validated module stack will make it
substantially easier to distinguish a project dependency regression from
an HPC infrastructure change in the future.

