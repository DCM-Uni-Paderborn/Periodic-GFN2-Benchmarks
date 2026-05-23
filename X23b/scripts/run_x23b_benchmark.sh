#!/usr/bin/env bash
set -eo pipefail

CP2K=${CP2K:-cp2k.psmp}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
CP2K_PARALLEL_JOBS=${CP2K_PARALLEL_JOBS:-20}
X23B_TASKS=${X23B_TASKS:-molecule_geoopt crystal_sp cellopt_gamma}

cd "$(dirname "$0")/.."

run_input() {
  local input=$1
  local rel=${input#inputs/}
  local run_dir="runs/${rel%.inp}"
  local output="$(basename "${input%.inp}").out"
  mkdir -p "${run_dir}"
  cp "${input}" "${run_dir}/"
  (
    cd "${run_dir}"
    if grep -q "ENERGY| Total FORCE_EVAL" "${output}" 2>/dev/null; then
      echo "SKIP ${rel}"
    else
      echo "RUN  ${rel}"
      rm -f "${output}"
      "${CP2K}" -i "$(basename "${input}")" -o "${output}"
    fi
  )
}
export -f run_input
export CP2K

job_file=$(mktemp)
trap 'rm -f "${job_file}"' EXIT
if [ "${X23B_PREPARE:-0}" = "1" ] || [ ! -d inputs ]; then
  python3 scripts/x23b_pipeline.py prepare
fi
for task in ${X23B_TASKS}; do
  find "inputs/${task}" -type f -name '*.inp' | sort >> "${job_file}"
done

xargs -n 1 -P "${CP2K_PARALLEL_JOBS}" bash -c 'run_input "$@"' _ < "${job_file}"

python3 scripts/x23b_pipeline.py analyse > x23b_results.out
