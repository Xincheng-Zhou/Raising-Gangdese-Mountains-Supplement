# Raising-Gangdese-Mountains-Supplement
This is the supporting information:
Zhou, X., Cao, W., Yang, J., Kaus, B. J. P., Ji, W.-Q., Gordon, S. M., & Zuza, A. V. (2026). Raising the Gangdese Mountains via Subduction of Indian Continental Crust, Not Slab Breakoff, During Subduction to Collision Transition. Tectonics, 45, e2026TC009363. https://doi.org/10.1029/2026TC009363

# User Guide for LaMEM Geodynamic Modeling Scripts

This repository contains Julia scripts and Jupyter notebooks for setting up and visualizing LaMEM geodynamic simulations of India–Eurasia convergence and Gangdese arc elevation evolution.

## File Overview

| File | Purpose |
|------|---------|
| `Ref.jl` | Julia script that sets up the reference model geometry, material phases, and boundary conditions, appends a two-stage pushing-block velocity condition, and launches the LaMEM simulation. Can be run directly on a personal computer or submitted to an HPC cluster. |
| `Ref.ipynb` | Jupyter notebook version of `Ref.jl` with identical model setup; intended for personal computer use only. Runs the setup cell by cell and generates preview plots of material phases and temperatures before launching the simulation. |
| `submit_ref.slm` | SLURM batch script for submitting `Ref.jl` to an HPC cluster |
| `_plot_isostatic_topo_julia_thickening_v12_2.ipynb` | Reads simulation output; generates animations of lithospheric phase field and three elevation indicators; produces topography–time evolution diagrams for the whole domain and the arc region |
| `_plot_field_julia_v2.ipynb` | Reads simulation output; plots phase, temperature, density, viscosity, strain rate, x/z velocities, and x-component of deviatoric stress as functions of time |
| `_plot_summary_v1_no_grid.ipynb` | Compares temporal evolution of arc surface topography (ST) across different model runs using `.csv` files exported by the plotting notebooks |

---

## 1. Install Julia

1. Go to <https://julialang.org/downloads/> and download **Julia v1.10.3** (the version used here; listed under "Older releases").
2. Run the installer and follow the on-screen instructions.
3. Verify the installation by opening a terminal and running:
   ```
   julia --version
   ```
   Expected output: `julia version 1.10.3`

---

## 2. Install Required Julia Packages

Open a Julia terminal and run the following. Installing all packages at once is recommended:

```julia
using Pkg

Pkg.add([
    PackageSpec(name="LaMEM",                     version="0.3.8"),
    PackageSpec(name="GeophysicalModelGenerator", version="0.7.7"),
    PackageSpec(name="Plots",                     version="1.40.5"),
    PackageSpec(name="SavitzkyGolay"),
    PackageSpec(name="CSV"),
    PackageSpec(name="DataFrames"),
    PackageSpec(name="Colors"),
    PackageSpec(name="ColorSchemes"),
])
```

`Statistics` is part of Julia's standard library and does not need to be installed separately.

> **Note on versions:** The version numbers above match the tested configuration. Omitting `version=` will install the latest available version, which may introduce breaking changes.

> **(Optional) Suppressing ReadVTK warnings:** The steps above are sufficient to run all scripts without errors. However, warnings may appear due to the ReadVTK package version not being pinned. To eliminate these warnings, manually downgrade ReadVTK to version 0.1.6:
> 
> ```julia
> using Pkg
> Pkg.add(PackageSpec(name="ReadVTK", version="0.1.6"))
> ```
---

## 3. Install Jupyter and IJulia (only required for the plotting notebooks)

The model setup script `Ref.jl` is a plain Julia script and does **not** require Jupyter. Jupyter is only needed to run the three `.ipynb` plotting notebooks.

If you do not already have Jupyter installed:

```bash
# Install JupyterLab via pip (requires Python)
pip install jupyterlab
```

Then, inside a Julia session, install the Julia kernel for Jupyter:

```julia
using Pkg
Pkg.add("IJulia")
```

Launch JupyterLab from the terminal:

```bash
jupyter lab
```

---

## 4. (Optional) Install ffmpeg for video export

The plotting notebooks call `ffmpeg` to assemble PNG frames into `.mp4` animations.


---

## 5. Workflow

### Step 1 — Set up and run the model (`Ref.jl`)

`Ref.jl` is a self-contained Julia script. It builds the full model geometry, assigns material phases and temperatures, defines rheological parameters, sets up phase transitions (eclogite transition and slab break-off), configures the multigrid solver, generates LaMEM input files via `prepare_lamem`, appends a two-stage pushing-block boundary condition to `output.dat`, and finally calls `run_lamem` to start the simulation.

Edit the parameter block at the top of `Ref.jl` to adjust model settings (e.g., `breakoff_depth`, `TarcMoho`, crustal thicknesses) before running.

**Option A — Run locally on a personal computer:**

```bash
julia Ref.jl
```
**Alternatively**, run `Ref.ipynb` cell by cell in JupyterLab, which provides preview plots of material phases and temperatures.

**Option B — Submit to an HPC cluster via SLURM:**

`submit_ref.slm` is provided as an **example** only. Before using it, you must edit the file to match your own cluster environment. At minimum update the working directory, output/error log paths, the absolute path to `Ref.jl`, and the account/partition names. The relevant lines to modify are:

```bash
#!/bin/bash
#SBATCH --chdir=/data/gpfs/assoc/lithospheric_dynamics/Zenodo
#SBATCH --output=/data/gpfs/assoc/lithospheric_dynamics/Zenodo/slurm-%j.out
#SBATCH --error=/data/gpfs/assoc/lithospheric_dynamics/Zenodo/slurm-%j.err
#SBATCH --job-name=Ref
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=compute_bound
#SBATCH --mem-per-cpu=3500
#SBATCH --account=cpu-s2-lithospheric_dynamics-0
#SBATCH --partition=cpu-s2-core-0
#SBATCH --mail-type=END,FAIL

julia /data/gpfs/assoc/lithospheric_dynamics/Zenodo/Ref.jl
```

Also adjust `--ntasks`, `--cpus-per-task`, `--time`, and `--mem-per-cpu` as appropriate for your cluster. The `prepare_lamem(model, 8)` call inside `Ref.jl` must match the total number of cores requested. Once edited, submit with:

```bash
sbatch submit_ref.slm
```

Simulation output (`.pvd`, `.vts`, surface files) will be written to the `output/` subdirectory of the working directory.

---

### Step 2 — Plot phase fields and physical fields

Place `_plot_field_julia_v2.ipynb` in the **same directory** as the simulation results (the folder containing the `output/` subdirectory) and run all cells in JupyterLab. It will:

- Loop over all saved time steps.
- Generate 8-panel PNG plots (phase, temperature, viscosity, strain rate, density, deviatoric stress, Vx, Vz) in a `plots_phase_vico_density_v2/` subfolder.
- Assemble the PNGs into `output_video.mp4` using ffmpeg.

---

### Step 3 — Plot isostatic topography and arc elevation

Place `_plot_isostatic_topo_julia_v12_2.ipynb` in the same simulation directory and run all cells. It will:

- Read passive tracer data to track arc crustal thickness over time.
- Compute isostatic elevation from the density field at each time step.
- Generate combined phase + topography PNG frames in a `plots_isostatic_arc_smooth_no_ecl_*` subfolder.
- Assemble frames into an `.mp4` animation.
- Save a topography–time heatmap and an arc elevation evolution (accounting for uplift due to magmatic thickening) versus time plot as standalone PNGs.
- Export arc simulated elevation (ST) to a `.csv` file (`data_<folder_name>.csv`) for downstream comparison.

---

### Step 4 — Compare multiple model runs

1. Collect the `.csv` files produced by Step 3 from all model runs into one folder.
2. Place `_plot_summary_v1_no_grid.ipynb` in that folder and run all cells.
3. The notebook reads each `.csv` and overlays arc elevation time series for different parameter choices (Indian crustal thickness, slab break-off depth, density, convergence rate, etc.).

---

## 6. Software Versions

| Software | Version |
|----------|---------|
| Julia | 1.10.3 |
| LaMEM.jl | 0.3.8 |
| GeophysicalModelGenerator.jl | 0.7.7 |
| Plots.jl | 1.40.5 |

More information about LaMEM.jl: <https://github.com/JuliaGeodynamics/LaMEM.jl>
