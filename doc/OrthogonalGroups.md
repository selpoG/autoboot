# Orthogonal groups O(n)

`getO[n]` supports integer degrees `n >= 2`, in both `group.m` and `ngroup.m`.
O(2) and O(3) keep their original representation labels. For `n >= 4`, put
`r = Floor[n/2]` and use `v[a1,...,ar,p]`: the first r entries are SO(n)
Dynkin labels and the final entry describes the disconnected component.
Spinor labels that do not descend to SO(n) are rejected.

| Case | Valid label | Meaning |
| --- | --- | --- |
| Odd n | `v[a1,...,ar,+1]` or `v[a1,...,ar,-1]` | Eigenvalue of central inversion `-I` |
| Even n, `a[r-1] == a[r]` | `v[a1,...,ar,+1]` or `v[a1,...,ar,-1]` | Two extensions of a reflection-invariant SO(n) irrep |
| Even n, `a[r-1] > a[r]` | `v[a1,...,ar,0]` | Induced irrep containing both SO(n) labels related by exchanging the last two nodes |

For even n, the canonical representative of a two-element orbit has the
penultimate Dynkin label greater than the last one. The reverse ordering is
rejected, as are `p=+/-1` on a non-fixed orbit and `p=0` on a fixed orbit.
The dimension of an induced irrep is twice that of its SO(n) component.
For a fixed label, the diagram-reflection intertwiner is normalized to fix
the highest-weight vector; `p` multiplies that intertwiner.

Examples after importing the group API:

```wl
G = getO[4];
G[dim[v[1,1,1]]]       (* 4: defining vector *)
G[dim[v[2,0,0]]]       (* 6: the two SO(4) chiral two-forms together *)
G[prod[v[1,1,1],v[1,1,1]]]  (* channels of dimensions 1, 6, 9 *)

G = getO[5];
G[dim[v[1,0,-1]]]      (* 5: defining vector, inversion-odd *)

G = getO[6];
G[dim[v[1,0,0,1]]]     (* 6: defining vector *)
G[dim[v[0,2,0,0]]]     (* 20: both SO(6) chiral three-forms *)
G[prod[v[0,2,0,0],v[0,0,0,-1]]]  (* {v[0,2,0,0]} *)
G[prod[v[0,1,1,1],v[0,1,1,1]]]
(* Both v[0,1,1,1] and v[0,1,1,-1] occur once. *)
```

The determinant representation is `v[0,...,0,-1]`, and the identity is
`v[0,...,0,+1]`. The irreps accepted here are tensor irreps of the compact
orthogonal group and are self-dual. O(2) and O(3) retain their existing
identity labels instead.

For odd n, the implementation uses `O(n) = SO(n) x {1,-I}`. For even n,
reflection exchanges the last two D-series simple roots. Its matrix is
transported from the highest vector using the same basis changes as the
SO(n) representation. In an induced irrep it exchanges the two SO(n) blocks.
It obeys `R^2=I` and `R F_i R=F_tau(i)` (also for raising generators).
Tensor multiplicities are split by reflection eigenvalues in the appropriate
highest-weight spaces. CG solutions impose that reflection equation as well
as every Lie-algebra intertwining equation.

`test/OrthogonalExact.wls` and `test/OrthogonalNumeric.wls` cover O(4) through
O(8) vector products, CG orthonormality and both types of intertwining equation,
three-equation vector bootstraps and SDP construction. Additional O(6) cases
check induced three-forms, determinant twists, invalid orbit labels, and the
two reflection parities in the adjoint square. Induced representations on input
legs are also checked for O(4) and O(6), including all CG channels. They also check compatibility
with the O(2)/O(3) dimension APIs and vector bootstraps. Arbitrary highest weights and degrees are
supported by the algorithm, but their running time is not uniformly bounded.
