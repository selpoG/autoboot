# Lie groups

Load `group.m` for exact arithmetic, or `ngroup.m` followed by
`setPrecision[50]` for numerical arithmetic. All constructors below work in
both modes. Returned objects implement `id`, `dim`, `dual`, `prod`, `isrep`,
`minrep`, `gA`, and `gG`, and can be passed to `setGroup` and `product`.

| Constructor | Group | Range |
| --- | --- | --- |
| `getSU[n]` | SU(n) | integer n >= 2; existing labels |
| `getO[n]` | O(n), including reflections | integer n >= 2 |
| `getSO[n]` | SO(n) | integer n >= 2 |
| `getSpin[n]` | Spin(n) | integer n >= 4; includes spinors |
| `getSp[n]` | compact Sp(n), also called USp(2n) | integer n >= 1; defining dimension 2n |
| `getLie["A",r]` | simply connected type A | integer r >= 1 |
| `getLie["B",r]` | Spin(2r+1) | integer r >= 2 |
| `getLie["C",r]` | Sp(r) | integer r >= 1 |
| `getLie["D",r]` | Spin(2r) | integer r >= 2 |
| `getLie["E",r]` | E6, E7, E8 | r = 6, 7, 8 |
| `getLie["F",4]` | F4 | rank 4 |
| `getLie["G",2]` | G2 | rank 2 |

D2 denotes Spin(4), the semisimple product SU(2) x SU(2).
Spin(3) can be accessed through `getSU[2]` or `getLie["A",1]`.
`getO[n]` includes the disconnected component for all n >= 2. See
[Orthogonal groups](OrthogonalGroups.md) for extension and induced-pair labels.
The Lie algebra generators span the complexified algebra; the implementation
constructs an orthonormal basis with each raising matrix adjoint to its lowering
matrix. The `gG` list is empty for these connected groups.

## Representation labels

`getLie`, `getSpin`, `getSp`, and `getSO[n]` for n >= 4 use
`v[a1,...,ar]`, where each entry is a nonnegative integer **Dynkin label**.
The trivial irrep is `v[0,...,0]`.

`getSU` retains its existing spin / Young diagram labels, and SO(2) and SO(3)
retain their existing labels. Thus `getSU[3]` and `getLie["A",2]` describe
isomorphic groups with different label conventions: the SU(3) adjoint is
`v[2,1]` in the first and `v[1,1]` in the second.

The Cartan convention is `A[i,j] = <alpha_j, alpha_i coroot>`, with Bourbaki
node numbering. See the [SageMath Cartan matrix reference](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/cartan_matrix.html)
and [type E diagrams](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/type_E.html).
For clarity, the matrices used here are available as
``LieRepresentations`cartan[type,rank]``.

- A, B, C and F have chain nodes 1,2,...,r. The last B node is short;
  the last C node is long. F has long nodes 1,2 and short nodes 3,4.
- D has chain 1,...,r-1 and a second end node r attached to r-2.
  In D2 the two nodes are disconnected.
- E has chain 1-3-4-5-...-r, with node 2 attached to node 4.
- G has short node 1 and long node 2.

SO(n) accepts only Spin(n) representations on which the covering kernel acts
trivially: the last label must be even for odd n; the sum of the last two
labels must be even for even n. `getSpin` accepts all nonnegative Dynkin labels.
These restrictions are checked by `isrep` and by all representation operations.
For example, `getSO[5][isrep[v[0,1]]]` is `False`, but
`getSpin[5][isrep[v[0,1]]]` is `True` (the 4-dimensional spinor).

## Examples

```Mathematica
Get["group.m"];
G = getSp[2];
G[dim[v[1,0]]]          (* 4 *)
G[prod[v[1,0],v[1,0]]]  (* irreps of dimensions 1, 5, 10 *)
setGroup[G];
isPseaudo[v[1,0]]       (* True *)

G = getLie["G",2];
G[dim[v[1,0]]]          (* 7 *)
G[dim[v[0,1]]]          (* 14 *)
G[prod[v[1,0],v[1,0]]]  (* irreps of dimensions 1, 7, 14, 27 *)

G = getLie["E",6];
G[dim[v[1,0,0,0,0,0]]]  (* 27 *)
G[dual[v[1,0,0,0,0,0]]] (* v[0,0,0,0,0,1] *)

G = getSO[6];
G[isrep[v[0,0,1]]]      (* False: use getSpin[6] for this spinor *)
```

## Construction and practical limits

Dimensions use the Weyl formula; characters use Freudenthal recursion.
Tensor products use Klimyk's shifted Weyl action on the weights of the
smaller factor, avoiding the characters of large output irreps.
Duals are obtained by Weyl reflection of the negative highest weight.

Generator matrices are built lazily from the highest weight. At each height,
the Chevalley relations determine the Gram matrix of lowered vectors. Its
nullspace is removed, and the remaining basis is orthonormalized. This works
uniformly across the supported Cartan types, without constructing a large
tensor product of fundamental representations.

The numerical mode converts these shared exact matrices to the configured
precision. It does not avoid the cost of the initial exact construction.
Large representations or tensor products can still require substantial time
and memory. The new regression suite exercises all seven series, including
E8's 248-dimensional representation, and checks dimensions, characters,
Chevalley and Serre relations, duals, tensor products and CG singlets.

## End-to-end bootstrap checks and timing

`test/BootstrapExact.wls` and `test/BootstrapNumeric.wls` run `setOps`,
`bootAll[]`, and `makeSDP` for SO(5), Sp(2), and G2. They check CG isometries
and completeness, the full four-index crossing tensors, crossing involution,
and the analytic singlet / symmetric-traceless / antisymmetric projectors for
SO(5). Final equations are compared as linear systems, allowing different row
bases and numerical rounding. The G2 case includes a singlet external scalar.

Representative `bootAll[]` times on an Intel Core i5-13500H, Wolfram Engine 15,
using separate exact / numerical test kernels and 50-digit numerical precision
(the three examples run in the listed order within each kernel):

| External operators | Exact | Numerical |
| --- | ---: | ---: |
| SO(5) vector | 0.5 s | 0.4 s |
| Sp(2) defining representation, with its automatically registered dual | 1.0 s | 0.8 s |
| G2 defining representation + singlet | 1.4 s | 1.1 s |

F4's 26-dimensional defining representation now completes the exact bootstrap
in about 37 seconds (about 34 seconds at 50-digit numerical precision), after
previously exceeding 120 seconds. Both modes produce five equations and pass
full tensor checks; see [validation and performance](LargeValidation.md).
The completed representation-matrix tests for E6/E7/E8 do not establish that
full bootstrap generation for those groups is practical; it has not been
benchmarked here. Group support therefore does not imply uniformly fast
bootstrap generation for all irreps.

These timings measure `bootAll[]` after ordinary `setOps` registration,
without manual CG precomputation. A separate standalone run of the numerical
G2 mixed benchmark took 3.3 seconds, including a cold runtime for that example. Registration initializes the
external singlet metric and takes about 0.01–0.2 seconds in these examples.
Kernel startup is excluded. `makeSDP` takes about 0.01 seconds or less. Times vary by machine and
load. Before optimization, the numerical G2 mixed case did not complete
within 120 seconds. Numerical equation reduction also needed a correctness
fix: terms with an implicit coefficient of 1 must be retained.

CG construction now solves only for highest-weight vectors of the tensor
product, then applies the exact basis changes recorded during target
representation construction. Numerical Lie groups reuse that exact solution
in the same generator basis. Singlet
multiplicities use Schur orthogonality, without decomposing a large tensor
square. CG and invariant tensors are stored as sparse arrays; three- and
four-point tensor contractions use sparse matrix products.

Reproduce a timed run and optionally save the equations with:

```sh
./benchmark/run-bootstrap.sh exact SO5 single /tmp/so5-equations.m
./benchmark/run-bootstrap.sh numeric G2 mixed /tmp/g2-equations.m
```

The wrapper applies a hard process time limit (150 seconds by default,
configurable with `BOOTSTRAP_TIMEOUT_SECONDS`) in addition to the script's
120-second bootstrap and 1.5 GB incremental-memory limits. Its failure exit
status distinguishes incomplete runs from successful equation generation.

For higher-rank classical groups, non-defining representations, and independent
GAP / analytic / built-in CG comparisons, see [LargeValidation.md](LargeValidation.md).
Run `make test-large` for the extended end-to-end suite.
