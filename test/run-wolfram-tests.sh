#!/bin/sh

set -eu

repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)

if ! command -v wolframscript >/dev/null 2>&1; then
	echo "wolframscript is required to run Wolfram Language tests." >&2
	exit 127
fi

for test_file in "$repo"/test/*.wls; do
	echo "Running ${test_file#"$repo"/}"
	wolframscript -f "$test_file"
done
