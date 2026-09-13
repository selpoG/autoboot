#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
# A hard external limit also stops non-interruptible symbolic kernel work.
exec timeout -s KILL "${BOOTSTRAP_TIMEOUT_SECONDS:-150}s" \
  wolframscript -f "$repo/benchmark/Bootstrap.wls" "$@"
