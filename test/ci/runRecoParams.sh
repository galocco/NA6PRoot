#!/usr/bin/env bash
# Static parameter validation: Python standard library only.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$script_dir/checkRecoParams.py" "$@"
