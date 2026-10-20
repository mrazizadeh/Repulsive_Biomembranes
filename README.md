# Repulsive Biomembranes - HPC version with TPE off

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
