#!/bin/sh

set -eu

repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)

for command in dotnet gap make python3; do
	if ! command -v "$command" >/dev/null 2>&1; then
		echo "$command is required to run license-free tests." >&2
		exit 127
	fi
done

python3 "$repo/test/check-pseudoreal-source.py"
make -C "$repo/GAPToMathematica"

test_tmp=$(mktemp -d "${TMPDIR:-/tmp}/autoboot-gap-test.XXXXXX")
trap 'rm -rf "$test_tmp"' EXIT HUP INT TERM

cd "$test_tmp"
printf 'Y\n8\n3\n' |
	dotnet "$repo/GAPToMathematica/bin/GAPToMathematica.dll"

if [ ! -s sg.8.3.m ]; then
	echo "GAPToMathematica did not generate sg.8.3.m" >&2
	exit 1
fi

echo "GAP_TO_MATHEMATICA_SMOKE=PASS"
