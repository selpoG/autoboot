#!/usr/bin/env python3
"""Regenerate independent GAP Lie-algebra fixtures (no autoboot imports)."""
import ast
from pathlib import Path
import subprocess
cases = [
    ('A', 2, [2,2], [1,1]), ('A', 3, [1,0,1], [1,0,1]),
    ('A', 4, [1,0,0,1], [1,0,0,1]),
    ('B', 2, [2,0], [2,0]), ('B', 3, [0,1,0], [0,1,0]),
    ('B', 5, [1,0,0,0,0], [1,0,0,0,0]),
    ('C', 4, [0,1,0,0], [0,1,0,0]),
    ('C', 5, [1,0,0,0,0], [1,0,0,0,0]),
    ('D', 5, [0,0,0,0,1], [0,0,0,0,1]),
    ('D', 6, [1,0,0,0,0,0], [1,0,0,0,0,0]),
    ('E', 6, [1,0,0,0,0,0], [1,0,0,0,0,0]),
    ('E', 7, [0,0,0,0,0,0,1], [0,0,0,0,0,0,1]),
    ('E', 8, [0,0,0,0,0,0,0,1], [0,0,0,0,0,0,0,1]),
    ('F', 4, [0,0,0,1], [0,0,0,1]), ('G', 2, [0,1], [0,1]),
]
script = 'SizeScreen([100000,100000]);;\n'
for t, r, a, b in cases:
    # GAP F4 nodes in Bourbaki order are [2, 4, 3, 1] (one-based).
    perm = [1, 3, 2, 0] if t == 'F' else list(range(r))
    inverse = [perm.index(i) for i in range(r)]
    ga, gb = [a[i] for i in inverse], [b[i] for i in inverse]
    script += f'''L:=SimpleLieAlgebra("{t}",{r},Rationals);;
a:={ga};; b:={gb};;
Print("DATA ", [CartanMatrix(RootSystem(L)), DimensionOfHighestWeightModule(L,a), DominantCharacter(L,a), DecomposeTensorProduct(L,a,b)],"\\n");
'''
script += 'Print("VERSION ", GAPInfo.Version,"\\n");\nQUIT;\n'
p = subprocess.run(['gap','-q','-b'], input=script, text=True, capture_output=True, check=True, timeout=120)
rows = [ast.literal_eval(line[5:]) for line in p.stdout.splitlines() if line.startswith('DATA ')]
if len(rows) != len(cases) or p.stderr:
    raise RuntimeError(p.stdout + p.stderr)
# Convert GAP F4 numbering to Bourbaki on both matrix axes and all weights.
for case, row in zip(cases, rows):
    if case[0] == 'F':
        perm = [1, 3, 2, 0]
        row[0] = [[row[0][i][j] for j in perm] for i in perm]
        for data in (row[2], row[3]):
            data[0] = [[w[i] for i in perm] for w in data[0]]
version = next(line[8:] for line in p.stdout.splitlines() if line.startswith('VERSION '))
def wl(x):
    if isinstance(x, str): return '"'+x+'"'
    if isinstance(x, (list,tuple)): return '{'+','.join(map(wl,x))+'}'
    return str(x)
path = Path(__file__).resolve().parents[1] / 'fixtures/LieGAP.m'
path.write_text('(* Generated independently with GAP '+version+'; see test/reference/generate-gap.py. *)\n'+wl([list(c)+row for c,row in zip(cases,rows)])+'\n')
print(f'Wrote {len(rows)} cases from GAP {version}: {path}')
