#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
o2_tag="${NA6P_CI_O2_TAG:-latest}"

set +u
source /cvmfs/alice.cern.ch/etc/login.sh
set -u
export ALIBUILD_WORK_DIR="${ALIBUILD_WORK_DIR:-/cvmfs/alice.cern.ch/sw}"
architecture="el9-$(uname -m)"
source "${ALIBUILD_WORK_DIR}/${architecture}/O2/${o2_tag}/etc/profile.d/init.sh"

build_dir="${repo_root}/build"
install_dir="${repo_root}/install"
cmake -S "${repo_root}" -B "${build_dir}" -DCMAKE_INSTALL_PREFIX="${install_dir}"
cmake --build "${build_dir}" --parallel "$(nproc)"
cmake --install "${build_dir}"

source "${install_dir}/init.sh"
params_test_dir="${build_dir}/Testing/params"
mkdir -p "${params_test_dir}"
(
  cd "${params_test_dir}"
  root -l -b -q "${repo_root}/test/ci/checkRecoParams.C"
  grep -Fq '[reco]' na6pRecoParam.ini
  grep -Fq 'mtZMatching=41.25' na6pRecoParam.ini
)
