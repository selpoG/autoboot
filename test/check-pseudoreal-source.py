#!/usr/bin/env python3

"""Check that exact and numerical pseudoreal sF definitions stay in sync."""

from pathlib import Path
import re
import sys


REPO = Path(__file__).resolve().parent.parent
PSEUDOREAL_SF = re.compile(
    r"^x:sF\[[^\n]*/; isPseaudo\[s\] := x =\n(?P<body>.*?)(?=\n\n)",
    re.DOTALL | re.MULTILINE,
)
SIGNED_PAIR = re.compile(
    r"\\\[Alpha\]\[o1, o2, o3, o4\]\[o, n, m\]\s*-\s*"
    r"\\\[Alpha\]\[o1, o2, o3, o4\]\[do, n, m\]"
)


def extract(path: Path) -> str:
    match = PSEUDOREAL_SF.search(path.read_text(encoding="utf-8"))
    if match is None:
        raise ValueError(f"pseudoreal sF definition not found in {path.name}")
    return match.group("body")


def main() -> int:
    try:
        exact = extract(REPO / "inv.m")
        numerical = extract(REPO / "ninv.m")
    except ValueError as error:
        print(error, file=sys.stderr)
        return 1

    if exact != numerical:
        print(
            "pseudoreal sF definitions differ between inv.m and ninv.m",
            file=sys.stderr,
        )
        return 1

    if SIGNED_PAIR.search(exact) is None:
        print(
            "pseudoreal sF does not apply the negative dual-operator sign",
            file=sys.stderr,
        )
        return 1

    print("PSEUDOREAL_SOURCE_SYNC=PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
