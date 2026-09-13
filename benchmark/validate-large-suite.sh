#!/bin/sh
# Opt-in: each case uses a fresh kernel, and each has an external hard limit.
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
results=${1:-$(mktemp -d /tmp/autoboot-large.XXXXXX)}
mkdir -p "$results"
results=$(CDPATH= cd -- "$results" && pwd)
printf 'Results: %s\n' "$results"
if [ "$#" -gt 0 ]; then shift; fi
if [ "$#" -eq 0 ]; then set -- SO7 SO10 SO12 SO5Adj SO5Sym2 A3Adj Sp4 F4; fi
for example do
  for mode in exact numeric; do
    printf 'Validating %s %s\n' "$mode" "$example"
    if "$repo/benchmark/run-large.sh" "$mode" "$example" "$results/$example.m" > "$results/$example-$mode.log" 2>&1; then
      cat "$results/$example-$mode.log"
    else
      cat "$results/$example-$mode.log"
      exit 1
    fi
  done
done
