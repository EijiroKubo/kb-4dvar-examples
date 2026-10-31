#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="${SCRIPT_DIR}/seafield"
FC="${FC:-gfortran}"
SOURCE="${SRC_DIR}/CalcSeaFieldMod_Alter.f90"

cleanup() {
    rm -f "${SRC_DIR}"/*.o "${SRC_DIR}"/*.obj "${SRC_DIR}"/*.mod
    rm -f "${SCRIPT_DIR}"/*.mod
}
trap cleanup EXIT

echo "===== Fortran shared library build start ====="

"${FC}" -O3 -fopenmp -shared -fPIC -J "${SRC_DIR}" \
    -o "${SRC_DIR}/CalcSeaFieldMod.so" "${SOURCE}"

echo
echo "Fortran shared library compiled successfully."
