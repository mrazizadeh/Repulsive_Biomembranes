# Repulsive Biomembranes - HPC version with hard exclusion

This branch contains the batch (non-GUI) version of Repulsive Biomembranes that runs on HPC clusters. The source is the source that is inside the Singularity image `biorsurfaces.sif` used for the production simulations. The `main` branch contains the interactive GUI version.

The solver minimizes the bending energy of a closed membrane with the bilayer-couple (BC) and area-difference elasticity (ADE) models. A tangent-point energy (TPE) prevents self-intersection, also at low reduced volumes. The method extends [Repulsive Surfaces](https://dl.acm.org/doi/abs/10.1145/3478513.3480521) from the [Keenan Crane Lab](https://www.cs.cmu.edu/~kmcrane/). The discrete bending energy follows [Mem3DG](https://github.com/RangamaniLabUCSD/Mem3DG).

## Differences from `main`

- No GUI. The program reads a scene file, runs a fixed number of optimization steps and writes a NetCDF trajectory.
- The scene file sets the model parameters directly: `reduced_volume`, `reduced_area`, `helfrich_weight`, `alpha_ade`, `write_index` and `iteration_limit`.
- Optional linear ramp of the reduced volume (`gradual`).

## Provenance

| Item | Value |
|---|---|
| Source history | `BioRepulsiveSurfaces`, branch `hpc-netcdf`, commit `0ed5d0f` (2024-02-15) |
| February 2024 image | Exactly commit `0ed5d0f` |
| May 2024 image (production runs) | Commit `0ed5d0f` plus the `gradual` reduced-volume ramp (second commit of this branch) |
| Executable in the May 2024 image | `/app/simulation/biorsurfaces`, SHA-256 `ccaa1b129bb4948ebca96f9661c3f95da4644ad8fbd7137f1d7531c1301c4c58` |
| Image software | Ubuntu 22.04.3, Intel oneAPI 2024.0 (MKL 2024.0, TBB 2021.11), Clang 14, CMake 3.22, Eigen 3.3.7 |

Editor settings, macOS metadata and old run outputs from the image are not included.

## Build

The image file is not part of this repository. Build an image from the recipe (route B). Then use the executable inside the image, or compile this source inside the image (route A).

### A. Inside the existing image (tested route)

On a Linux x86-64 node with Apptainer or Singularity:

```sh
bash singularity/build_in_image.sh /path/to/biorsurfaces.sif build 16
```

The script mounts this source read-only and uses the compiler, MKL, TBB, NetCDF and Eigen inside the image. The build does not download anything. The Release flags include `-march=native`. Build on the node type that runs the solver.

### B. New image from the recipe

`singularity/Dockerfile` is the recipe of the image. Build from the repository root:

```sh
docker build -f singularity/Dockerfile -t biorsurfaces .
docker save biorsurfaces > biorsurfaces.tar
singularity build biorsurfaces.sif docker-archive://biorsurfaces.tar
```

The original recipe used `intel/oneapi-hpckit:latest`. The recipe here pins `intel/oneapi-hpckit:2024.0.1-devel-ubuntu22.04`, which matches the software in the image. A new image from this pinned recipe has not been built yet. During the build, CMake downloads Eigen.

### C. Native build

Dependencies: Intel oneAPI MKL (`source /opt/intel/oneapi/setvars.sh`), TBB, OpenMP, NetCDF-C, NetCDF-C++4 and Clang (`export CC=clang CXX=clang++`).

```sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 8
```

The solver uses MKL sparse routines and PARDISO directly. A build without MKL needs code changes.

## Run

Always use a new, empty directory for each run. The solver replaces `mesh_data.nc` and `log_Energy.txt` in the working directory.

```sh
mkdir -p run1/mesh && cp sample_input/oblate.obj sample_input/sceneADE.txt run1/ && cd run1
singularity exec --cleanenv --bind "$PWD:/work" --pwd /work \
    /path/to/biorsurfaces.sif /app/simulation/biorsurfaces sceneADE.txt --threads 8
```

For a binary from route A, mount the build directory and call that binary:

```sh
singularity exec --cleanenv --bind "$PWD:/work" --pwd /work \
    --bind /path/to/build:/biorsurfaces/build:ro \
    /path/to/biorsurfaces.sif /biorsurfaces/build/bin/biorsurfaces sceneADE.txt --threads 8
```

Options: `--threads N` sets the number of threads. `--theta T` sets the Barnes-Hut parameter (default 0.3; 0 is exact).

### Scene file

```text
repel_mesh oblate.obj 6 12      # input mesh and TPE exponents
reduced_volume 0.9
reduced_area 1.05               # preferred reduced area difference
helfrich_weight 5000
alpha_ade 1000                  # ADE stiffness
write_index 10                  # mesh frame interval
constrain barycenter
constrain volume 1
constrain reduced_volume 1
minimize helfrichADE 1
iteration_limit 1000
```

Optional reduced-volume ramp: `gradual <initial_rv> <final_rv> <steps>`. The target changes by `(final_rv - initial_rv)/steps` after each step, for `steps` steps.

## Output

`mesh_data.nc` (NetCDF-4):

| Name | Content |
|---|---|
| `vertices`, `faces` | Mesh frames, one frame each `write_index` steps |
| `HelfrichEnergy` | Bending energy divided by `helfrich_weight`, one record per step |
| `TPEEnergy` | Tangent-point energy, one record per step |
| `dela`, `area`, `volume`, `ReducedVolume` | One record per step |
| Attributes | `reduced_volume`, `reduced_area`, `weight_Helfrich`, `ADE_alpha`, `NetCDF_Step` |

Step indexing: the loop counter `i` starts at 0, and each record is written after step `i + 1`. Energy record `k` is the state after step `k + 1`. Mesh frame `k` is the state after step `k * write_index + 1`.

`log_Energy.txt` contains one line per step (bending energy / weight, TPE / weight). Each line is written before its step.

`mesh/frame0001.obj` is the final mesh (6 significant digits). Create the `mesh/` directory before the run; otherwise this file is not written.

## License and acknowledgments

MIT license (`LICENSE`). This project is derived from [Repulsive Surfaces](https://github.com/chrisyu-cs/repulsive-surfaces) (`LICENSE_Repulsive_Surfaces`) and [Mem3DG](https://github.com/RangamaniLabUCSD/Mem3DG) (`LICENSE_MEM3DG.txt`). Vendored libraries in `deps/` keep their own licenses.

## TPE-off version (this branch)

This branch adds one scene line to the HPC version:

```text
repulsion off
```

With `repulsion off`, the tangent-point energy (TPE) is removed from the objective and from the gradient. The TPE object stays in the flow, so the Hs (fractional Sobolev) metric is built exactly as with the repulsion on. TPE off therefore means no TPE energy and no TPE force. It does not remove the proximity terms of the Hs metric. `repulsion on` (or no `repulsion` line) gives the original behavior.

The option needs the Barnes-Hut TPE (`--theta` greater than 0, default 0.3). At startup the program prints `*** REPULSION OFF ... ***`.

Output changes with `repulsion off`:

- `TPEEnergy` is 0, because the TPE is not in the objective.
- `TPEEnergyDiagnostic` contains the TPE value of the current shape. It is a proximity diagnostic only.
- The NetCDF attribute `repulsion` is `on` or `off`.

Example: `sample_input/sceneADE_repulsion_off.txt`.

Changed files relative to `hpc-version`: `include/energy/tpe_barnes_hut_0.h`, `src/energy/tpe_barnes_hut_0.cpp`, `include/scene_file.h`, `src/scene_file.cpp`, `src/main.cpp`. This source was compiled inside the May 2024 image with `-DGC_EIGEN_LOCATION` set to the Eigen copy in the image (route A above).

## Hard-exclusion version (this branch)

This branch adds a geometric exclusion guard to the TPE-off version. The guard prevents self-intersection independently of the tangent-point energy. Thus `repulsion on` and `repulsion off` can be compared with the same admissible-motion rule. The guard uses continuous collision detection (CCD) from the [Tight-Inclusion](https://github.com/Continuous-Collision-Detection/Tight-Inclusion) library, pinned at commit `6f84001a790a9d3362b7c13d410037140d83fa2c`, in double precision.

### How the guard works

Each optimization step has three phases, and the guard checks each phase separately:

1. **Descent.** The complete linear motion of a line-search trial is checked before its energy is evaluated. Unsafe or uncertain trials are shortened and checked again.
2. **Constraint correction.** Each correction is checked and, if necessary, shortened. The cached Schur correction is repeated until the relative area and volume errors are at most `hard_constraint_tolerance`, with a maximum of 20 corrections per step.
3. **Barycenter projection.** The candidate motion is checked before it is applied.

A step is accepted only if all three phases are safe, the connectivity is unchanged, area and volume are finite and positive, the relative area and volume errors are at most `hard_constraint_tolerance`, and the objective does not increase (relative 1e-8). Otherwise the program restores the state at the start of the step, writes the reason to `hard_exclusion_status.json` and stops with a nonzero exit code. A rejected step is never written as a trajectory frame.

The CCD tests nonincident vertex-face and edge-edge pairs, with candidates from swept bounding boxes. A quadratic normal test rejects triangle collapse or inversion during a motion. With `hard_exclusion on`, the GMRES solve of the Hs direction is also checked for convergence (status, residual, finite values). A failed solve stops the run and restores the step start.

### Scene lines

```text
method hs
remesh off
hard_exclusion on
repulsion off                         # or on
hard_min_separation_relative 1e-06    # clearance / initial mean edge length
hard_ccd_tolerance_relative 1e-8      # CCD tolerance / initial mean edge length
hard_constraint_tolerance 0.001       # maximum relative area and volume error
```

`remesh on|off` is also available without the guard. The guard supports only this restricted scene set:

- One closed, connected, manifold membrane; no obstacles, barriers, pins or ramps (`gradual`).
- `method hs`, `remesh off` and `write_index 1`.
- Exactly `minimize helfrichADE 1`.
- Exactly the constraints `constrain barycenter`, `constrain volume 1`, `constrain reduced_volume 1`, in this order.
- `hard_ccd_tolerance_relative` smaller than `hard_min_separation_relative`; `hard_constraint_tolerance` at most 0.001.

Other scenes stop with an error before the optimization starts. The clearance and CCD tolerance are numerical contact thresholds, not a membrane thickness.

### Starting mesh

**The starting mesh must be free of self-intersections.** The static vertex-face and edge-edge tests cannot detect every intersection that already exists in the input. Check the input mesh independently first, for example with CGAL (`PMP::does_self_intersect`) or the MeshLab self-intersection filter. Then confirm the check to the solver with an environment variable:

```sh
export HARD_EXCLUSION_INPUT_AUDITED=1
```

Without this variable, `hard_exclusion on` stops with an error that asks for an audited input. The message refers to the audited runner of the original study; the variable is the same gate.

The guard computes the area target from the start volume and `reduced_volume`. If the start mesh differs from the area or volume target by more than `hard_constraint_tolerance`, the run stops with `invalid_input`. Set `reduced_volume` to the reduced volume of the start mesh.

Example: `sample_input/scene_hard_exclusion.txt`, 20 steps from `sample_input/oblate.obj`. This mesh is closed, genus 0, has no self-intersecting faces (MeshLab check) and has reduced volume 0.96757.

### Output

In addition to `mesh_data.nc` (every accepted step, full precision):

| File | Content |
|---|---|
| `hard_exclusion_steps.csv` | Accepted steps: area and volume residuals, objective change, displacement |
| `hard_exclusion_motions.csv` | Every checked motion, including rejected line-search trials: fraction, backtracks, query counts |
| `hard_exclusion_linear_solves.csv` | GMRES status and residual of each solve |
| `hard_exclusion_config.json` | Thresholds and targets |
| `hard_exclusion_status.json` | `iteration_limit_reached`, `stopped_inconclusive` or `invalid_input`, with the reason |
| `hard_exclusion_failure_state.obj`, `hard_exclusion_restored_state.obj` | Geometry at a caught failure (before rollback) and the restored geometry |
| `mesh/frame0001.obj` | Final accepted mesh, 17 significant digits (create `mesh/` before the run) |

`iteration_limit_reached` means that the run completed its steps. `stopped_inconclusive` can be caused by contact, constraints, no descent or numerical stagnation. Neither status shows equilibrium.

### Tests

The build makes three test programs in `build/bin/`:

- `hard_guard_test`: synthetic CCD cases (vertex-face tunneling, edge-edge crossing, coplanar overlap, legitimate adjacency, folding, intermediate collapse, safe rigid motion). `hard_guard_test --inspect mesh.obj [relative_clearance]` checks a mesh for topology and feature clearance. This is not a replacement for an independent self-intersection check of the input.
- `hard_gmres_test`: rejection of an unconverged GMRES solve.
- `hard_constraint_test`: constraint correction, budget exhaustion and invalid geometry.

### Implementation notes

- New code: `include/hard_exclusion/`, `src/hard_exclusion_bridge.cpp`, `src/hard_exclusion_guard.cpp`, `tests/`. Changed files: `CMakeLists.txt`, `include/scene_file.h`, `include/sobolev/hs_iterative.h`, `include/sobolev/hs_schur.h`, `src/line_search.cpp`, `src/main.cpp`, `src/obj_writer.cpp`, `src/scene_file.cpp`, `src/sobolev/hs.cpp`.
- `third_party/` contains Tight-Inclusion and spdlog 1.11.0 with their licenses. The build does not download them.
- The solver stays in C++11; the CCD library and tests use C++17. Floating-point contraction and fast-math are disabled for the solver and the CCD code.
- The Helfrich/ADE energy, the constraint derivatives and the TPE-based Hs metric are unchanged. With `hard_exclusion off`, the guard and the GMRES check are inactive.
- This source was compiled inside the May 2024 image (route A), and the three test programs and two-step guarded runs passed there. In this branch, only the `HARD_PACKAGE` path in `CMakeLists.txt` differs, because `third_party/` and `tests/` sit at the repository root.
