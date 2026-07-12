#!/usr/bin/env bash
set -eo pipefail

CP2K=${CP2K:-cp2k.psmp}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
CP2K_PARALLEL_JOBS=${CP2K_PARALLEL_JOBS:-20}

phases=(Ih II III IV VI VII VIII IX XI XIII XIV XV XVII)
methods=(GFN1 GFN2)

cd "$(dirname "$0")/.."

run_one() {
  local method=$1
  local phase=$2
  local run_dir="runs/${method}/${phase}"
  local input="ice_${phase}_${method}.inp"
  local output="ice_${phase}_${method}.out"
  mkdir -p "${run_dir}"
  cp "inputs/${input}" "${run_dir}/"
  (
    cd "${run_dir}"
    if [[ ${FORCE_RERUN:-0} != 1 ]] && grep -q "ENERGY| Total FORCE_EVAL" "${output}" 2>/dev/null; then
      echo "SKIP ${method} ${phase}"
    else
      echo "RUN  ${method} ${phase}"
      rm -f "${output}"
      "${CP2K}" -i "${input}" -o "${output}"
    fi
  )
}
export -f run_one
export CP2K

job_file=$(mktemp)
trap 'rm -f "${job_file}"' EXIT
for method in "${methods[@]}"; do
  for phase in "${phases[@]}"; do
    printf '%s %s\n' "${method}" "${phase}" >> "${job_file}"
  done
done

xargs -n 2 -P "${CP2K_PARALLEL_JOBS}" bash -c 'run_one "$@"' _ < "${job_file}"

python3 scripts/dmc_ice13_pipeline.py > results.out
