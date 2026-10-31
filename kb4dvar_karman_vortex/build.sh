#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="${SCRIPT_DIR}/ns_model"
FC="${FC:-gfortran}"
FFLAGS=(-O3 -fPIC -cpp -fno-range-check -I "${SRC_DIR}" -J "${SRC_DIR}")

SOURCES=(
    mod_variables
    sub_initial
    sub_utils
    sub_bc_outer
    sub_bc_wall
    sub_rhs3rd_fwd
    sub_hsmac_fwd
    calc_forward_ns
)

cleanup() {
    rm -f "${SRC_DIR}"/*.o "${SRC_DIR}"/*.mod
}
trap cleanup EXIT

echo "===== Fortran shared library build start ====="

OBJECTS=()
for source in "${SOURCES[@]}"; do
    object="${SRC_DIR}/${source}.o"
    "${FC}" "${FFLAGS[@]}" -c "${SRC_DIR}/${source}.f90" -o "${object}"
    OBJECTS+=("${object}")
done

"${FC}" -fopenmp -shared -o "${SRC_DIR}/calc_forward_ns.so" "${OBJECTS[@]}"

echo
echo "Fortran shared library compiled successfully."
