#!/usr/bin/env bash
# Compile this source tree with the compiler and libraries inside an existing
# biorsurfaces.sif image. The image and the source are mounted read-only; all
# output goes to BUILD_DIR. Eigen is taken from the copy inside the image, so
# the build does not download anything.
#
# Usage: bash singularity/build_in_image.sh IMAGE.sif BUILD_DIR [JOBS]
# Run it on the node type that will run the solver: the Release flags include
# -march=native.
set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo 'Usage: bash singularity/build_in_image.sh IMAGE.sif BUILD_DIR [JOBS]' >&2
    exit 2
fi
image=$(cd -- "$(dirname -- "$1")" && pwd -P)/$(basename -- "$1")
[[ -r "$image" ]] || { echo "Cannot read image: $image" >&2; exit 2; }
source_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
mkdir -p "$2"
build_dir=$(cd -- "$2" && pwd -P)
jobs=${3:-$(nproc)}
runtime=$(command -v apptainer || command -v singularity || true)
[[ -n "$runtime" ]] || { echo 'Needs apptainer or singularity in PATH.' >&2; exit 2; }

run=("$runtime" exec --cleanenv
     --bind "$source_dir:/biorsurfaces/src:ro"
     --bind "$build_dir:/biorsurfaces/build"
     "$image")
"${run[@]}" cmake -S /biorsurfaces/src -B /biorsurfaces/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=/usr/bin/clang \
    -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
    -DGC_EIGEN_LOCATION=/app/biorsurfaces/build/deps/geometry-central/deps/eigen-src \
    2>&1 | tee "$build_dir/configure.log"
"${run[@]}" cmake --build /biorsurfaces/build --parallel "$jobs" \
    2>&1 | tee "$build_dir/compile.log"
[[ -x "$build_dir/bin/biorsurfaces" ]] || { echo "No executable at $build_dir/bin/biorsurfaces" >&2; exit 1; }
"${run[@]}" ldd /biorsurfaces/build/bin/biorsurfaces > "$build_dir/ldd.txt"
if grep -q 'not found' "$build_dir/ldd.txt"; then
    echo "Unresolved libraries; see $build_dir/ldd.txt" >&2
    exit 1
fi
echo "BUILD OK: $build_dir/bin/biorsurfaces"
echo "Run inside the same image, for example:"
echo "  $runtime exec --cleanenv --bind \"\$PWD:/work\" --pwd /work --bind \"$build_dir:/biorsurfaces/build:ro\" \\"
echo "      $image /biorsurfaces/build/bin/biorsurfaces scene.txt --threads 8"
