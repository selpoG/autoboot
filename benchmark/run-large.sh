#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
exec timeout -s KILL "${BOOTSTRAP_TIMEOUT_SECONDS:-180}s" \
  wolframscript -f "$repo/benchmark/ValidateLarge.wls" "$@"
