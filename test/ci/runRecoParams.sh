#!/usr/bin/env bash
# Called inside AlmaLinux 9 by the GitHub-hosted workflow.
set -eo pipefail

if [[ "${1:-}" == --test ]]; then
  cmake -S . -B build -DBUILD_TESTING=ON -DWITH_GEANT4_UIVIS=OFF \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
  cmake --build build --parallel 2
  cmake --install build
  source install/init.sh
  ctest --test-dir build --output-on-failure -R '^na6prec_params$'
  exit 0
fi

dnf install -y cmake make gcc gcc-c++ gcc-gfortran python3 which findutils tcl \
  git tar gzip bzip2 unzip file procps-ng perl zlib-devel \
  libX11 libXpm libXft libXext libXrender libSM libICE libXi libXmu \
  libXrandr libXinerama libXcursor mesa-libGL mesa-libGLU \
  libxml2 libnsl libtirpc openssl-libs

module_dir=/cvmfs/alice.cern.ch/el9-x86_64/Modules/modulefiles/O2PDPSuite
test -x /cvmfs/alice.cern.ch/bin/alienv
test -d "$module_dir"
release_tag=${NA6P_CI_O2_TAG:-}
if [[ -z "$release_tag" ]]; then
  release_tag=$(find "$module_dir" -maxdepth 1 -name 'daily-*' -printf '%f\n' |
    sort -V | tail -n 1)
fi
if [[ -z "$release_tag" || "$release_tag" == *[!A-Za-z0-9._-]* ]]; then
  echo "Invalid or missing O2PDPSuite release: $release_tag" >&2
  exit 1
fi
if [[ ! -e "$module_dir/$release_tag" ]]; then
  echo "O2PDPSuite release not available: $release_tag" >&2
  exit 1
fi
echo "Using O2PDPSuite/$release_tag on AlmaLinux 9"
exec /cvmfs/alice.cern.ch/bin/alienv setenv "O2PDPSuite/$release_tag" \
  -c bash test/ci/runRecoParams.sh --test
