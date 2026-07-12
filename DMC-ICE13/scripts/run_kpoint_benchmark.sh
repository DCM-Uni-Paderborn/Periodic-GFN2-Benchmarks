#!/usr/bin/env bash
set -eo pipefail

CP2K=${CP2K:-cp2k.psmp}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
CP2K_PARALLEL_JOBS=${CP2K_PARALLEL_JOBS:-20}

phases=(Ih II III IV VI VII VIII IX XI XIII XIV XV XVII)
methods=(GFN1 GFN2)
meshes=(k111 k222 k333 k444 k555)

cd "$(dirname "$0")/.."

python3 scripts/dmc_ice13_kpoint_benchmark.py prepare

run_one() {
  local mesh=$1
  local method=$2
  local phase=$3
  local run_dir="runs_kpoints/${mesh}/${method}/${phase}"
  local input="ice_${phase}_${method}_${mesh}.inp"
  local output="ice_${phase}_${method}_${mesh}.out"
  mkdir -p "${run_dir}"
  cp "kpoint_inputs/${mesh}/${input}" "${run_dir}/"
  (
    cd "${run_dir}"
    if [[ ${FORCE_RERUN:-0} != 1 ]] && grep -q "ENERGY| Total FORCE_EVAL" "${output}" 2>/dev/null; then
      echo "SKIP ${mesh} ${method} ${phase}"
    else
      echo "RUN  ${mesh} ${method} ${phase}"
      rm -f "${output}"
      "${CP2K}" -i "${input}" -o "${output}"
    fi
  )
}
export -f run_one
export CP2K

job_file=$(mktemp)
trap 'rm -f "${job_file}"' EXIT
for mesh in "${meshes[@]}"; do
  for method in "${methods[@]}"; do
    for phase in "${phases[@]}"; do
      printf '%s %s %s\n' "${mesh}" "${method}" "${phase}" >> "${job_file}"
    done
  done
done

xargs -n 3 -P "${CP2K_PARALLEL_JOBS}" bash -c 'run_one "$@"' _ < "${job_file}"

python3 scripts/dmc_ice13_kpoint_benchmark.py analyse > kpoint_results.out
