#!/usr/bin/env python3
"""Generate blocks_3d tables for one identical 3d primary, without shell evaluation.

This is an explicit finite-spin, finite-derivative truncation. Increase all
truncations and precision to test convergence before interpreting a bound.
"""
import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path
import subprocess


def fraction(value):
    return Fraction(value)


def jobs(spin):
    qs = [Fraction(k, 2) for k in range(-int(2 * spin), int(2 * spin) + 1, 2)]
    for q in itertools.product(qs, repeat=4):
        if (4 * spin - sum(q)) % 2:
            continue
        nonzero = [v for v in q if v]
        if nonzero and nonzero[0] < 0:
            continue
        for sign, left, right in itertools.product((1, -1), range(int(2 * spin) + 1), range(int(2 * spin) + 1)):
            yield q, sign, left, right


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--blocks-3d', required=True)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--spin', type=fraction, default=Fraction(1, 2))
    parser.add_argument('--dimension', type=fraction, default=Fraction(6, 5))
    parser.add_argument('--max-spin', type=int, default=4)
    parser.add_argument('--lambda', dest='derivatives', type=int, default=3)
    parser.add_argument('--order', type=int, default=12)
    parser.add_argument('--kept-pole-order', type=int, default=6)
    parser.add_argument('--precision', type=int, default=300)
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    if args.spin < 0 or (2 * args.spin).denominator != 1 or args.dimension <= 0:
        parser.error('spin must be a nonnegative half integer and dimension positive')
    if min(args.max_spin, args.derivatives) < 0 or min(args.order, args.kept_pole_order, args.threads) <= 0 or args.precision < 128:
        parser.error('invalid truncation, thread count, or precision')
    # Decimal external dimensions must be exact (no binary-float rounding).
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = args.precision + 20
        dimension = format(Decimal((2 * args.dimension).numerator) / Decimal((2 * args.dimension).denominator), 'f')
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = []
    for number, (q, sign, left, right) in enumerate(jobs(args.spin)):
        directory = args.output / str(number)
        cmd = [args.blocks_3d,
               '--j-external', ','.join([str(float(args.spin))] * 4),
               '--j-internal', ','.join(map(str, range(args.max_spin + 1))),
               '--j-12', str(left), '--j-43', str(right),
               '--delta-12', '0', '--delta-43', '0', '--delta-1-plus-2', dimension,
               '--four-pt-struct', ','.join(str(float(v)) for v in q),
               '--four-pt-sign', str(sign), '--order', str(args.order),
               '--lambda', str(args.derivatives), '--kept-pole-order', str(args.kept_pole_order),
               '--precision', str(args.precision), '--num-threads', str(args.threads),
               '--coordinates', 'xt', '--profile', '0', '-o', str(directory / 'derivs_{}.json')]
        print(f'blocks_3d job {number + 1}: q={q}, sign={sign}, j12={left}, j43={right}', flush=True)
        subprocess.run(cmd, check=True)
        manifest.extend(str(Path(str(number)) / f'derivs_{j}.json') for j in range(args.max_spin + 1))
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
