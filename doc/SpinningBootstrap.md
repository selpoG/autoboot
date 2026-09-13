# Spinning bootstrap in three dimensions

The spinning packages generate spacetime tensor structures separately from the
internal-symmetry representations used by `group.m`. They are an initial,
explicit API for issue #2; the scalar `bootAll[]` interface is unchanged.

## Scope and research basis

The numerical target is **3d, parity-invariant, identical Hermitian primaries in
long conformal representations**. Integer and half-integer external spins share
the same q-basis and SO(3)-basis implementation. The worked numerical example is
four Majorana fermions. Mixed external spins are supported at the tensor-structure
and component-crossing level, but automatic mixed-correlator/flavor SDP assembly
is not implemented.

This choice follows the availability of general-spin 3d blocks in
[blocks_3d](https://gitlab.com/bootstrapcollaboration/blocks_3d), whose
[paper](https://arxiv.org/abs/2011.01959) specifies the bases, normalizations,
identity contribution, and a complete Majorana example. Fermionic mixed
correlators already produce numerical islands, as illustrated by the
[Gross–Neveu–Yukawa archipelago](https://arxiv.org/abs/2210.02492).

Conserved currents and stress tensors are important next targets, but require
conservation equations, Ward identities, and treatment of shortened exchanged
representations. The
[3d Ising stress-tensor bootstrap](https://arxiv.org/abs/2411.15300) and
[3d U(1) current atlas](https://arxiv.org/abs/2412.01608) demonstrate their physical
value. A spin label and a dimension at the unitarity bound do not implement these
constraints. The numerical adapter therefore rejects external unitarity-bound
operators and exchanged gaps at or below a long-block pole. In particular, the
example does **not** include the conserved stress-tensor contribution needed for
a local-CFT bound.

There are also numerical spinning results in four dimensions, including
[bounds on abelian currents](https://arxiv.org/abs/2512.20710). They do not make
3d blocks a backend for arbitrary dimensions; 4d tensor conventions and block
integration are outside this implementation.

The SDPB target checked on 2026-09-14 is
[SDPB 3.1.0](https://github.com/davidsd/sdpb/releases/tag/3.1.0).
SDPB consumes polynomial matrix programs; it does not supply spinning conformal
blocks or conservation constraints.

## Structure and crossing API

Run from the repository root:

```Mathematica
Get["spinning.m"];
basis = spinThreePointBasis[{1/2, 1/2, 2}, 1, -1];
basis["QCoefficients"]
basis["SO3Coefficients"]

(* Scalar/fermion mixed component crossing; s and t are ordered channel functions. *)
spinCrossingEquations[{0, 1/2, 0, 1/2}, 0, s, t, z, zb]
```

Spins must be exact nonnegative half integers, with an even number of fermions.
Three-point q labels obey `Total[q]==0`. `spinSO3ToQ` expresses SO(3) structures
as rows in the q basis; `spinQToSO3` is its inverse, computed using CG
orthogonality rather than generic algebraic matrix inversion.

`spinThreePointBasis[spins, parity, exchange, constraints]` intersects parity,
(12) exchange, and optional homogeneous linear constraints on q coefficients.
`exchange` is `0` for distinct operators or the required eigenvalue `+1`/`-1`.
For identical fermions without flavor it is `-1`. With internal representations,
the spacetime eigenvalue must also include the CG exchange eigenvalue.
`constraints` must be supplied explicitly; the function does not derive
conservation equations.

The returned `RealityPhase` is `I` when two three-point operators are fermionic,
and `1` otherwise. Thus four identical fermions contribute **minus** a real OPE
quadratic form before contracting with the block. The right three-point basis is
ordered **43**, following blocks_3d. Crossing already includes the fermionic
permutation sign; adding another statistics sign changes the equations.

`spinCrossingMatrix` returns both ordered q-label lists. Its columns use
`CrossedStructures`, not an independently sorted canonical basis. This matters
for mixed external spins.

## Numerical Majorana example

Install an activated Wolfram Engine and build blocks_3d using its upstream
instructions. autoboot neither downloads nor builds external tools at runtime.
For a small **input-generation smoke test**:

```sh
python3 examples/spinning/generate-blocks.py \
  --blocks-3d /path/to/blocks_3d --output /tmp/majorana-blocks \
  --max-spin 4 --lambda 3 --order 12 --kept-pole-order 6 --precision 300
wolframscript -f examples/spinning/majorana.wls \
  /tmp/majorana-blocks /tmp/majorana-pmp.json
pmp2sdp --precision=256 --input=/tmp/majorana-pmp.json --output=/tmp/majorana-sdp
sdpb --sdpDir=/tmp/majorana-sdp
```

The generator defaults to external spin `1/2`, dimension `6/5`. It also accepts
`--spin`, `--dimension`, and explicit truncations. Output includes a manifest;
use a new directory for a different calculation. Generated tables and SDP output
belong outside the source tree.

The Mathematica steps are available separately:

- `readSpinningBlocks[file, digits]`: parses JSON decimal strings without
  evaluating code and checks working precision and required metadata.
- `spinBlockPolynomial[data,{j120,j430},{m,n},x]`: extracts the `xt` derivative
  numerator. It restores the imaginary phase for half-integer exchanged spin.
- `spinIdenticalSDP[tables,sectors,x,lambda]`: contracts allowed three-point
  structures and forms real OPE matrices. Each sector specifies `Spin`, `Parity`,
  and `Gap`; the polynomial variable is `x = Delta - Gap >= 0`.
- `spinMajoranaRows[sdp]`: selects the independent Majorana derivative components
  of [Counting Conformal Correlators, appendix A.1](https://arxiv.org/abs/1612.08987),
  including removal of dependent `++++` derivatives on `z=zb`. At `lambda=3`
  there are 13 components.
- `spinWritePMP[file,sdp,x,digits]`: emits SDPB 3.x PMP JSON, rejecting symbolic
  blocks, insufficient numerical precision, asymmetric matrices, and prefactors
  with a pole on the positive half line.

Different SO(3) structures can have different pole sets. The adapter takes a
common denominator with the maximum multiplicity of each pole and multiplies
individual numerators by the missing factors **before** adding them. The gap
shift applies to both polynomial numerators and prefactors.

blocks_3d uses `x=(z+zb-1)/2`, `t=((z-zb)/2)^2` for derivative coordinates
(distinct from the dimension polynomial variable). Its odd component is divided
by `(z-zb)/2` before differentiation. Crossing therefore has an extra minus sign
for this regularized odd component.

`spinMakeSDP[identity,sectors]` is also usable independently with user-supplied
homogeneous real OPE quadratic forms. Off-diagonal entries are half the mixed
monomial coefficient. The exported objective is zero and the identity vector is
the functional normalization, for a feasibility/exclusion problem.

## Validation and limits

For an external integration check, generate two Majorana table directories
with `--order 12 --kept-pole-order 6` and
`--order 24 --kept-pole-order 12` (both `--precision 300 --max-spin 4 --lambda 3`),
then run:

```sh
wolframscript -f benchmark/ValidateSpinning.wls \
  /tmp/majorana-low /tmp/majorana-high /path/to/sdpb/mathematica/SDPB.m
```


`make test-wolfram` includes exact q/SO(3) inversion and reconstruction, the
published Majorana selection rules and basis conversion, crossing phases,
identity normalization, quadratic-form assembly, JSON interfaces, and rejection
of invalid inputs. The synthetic block interface tests are not physical block
checks.

External integration was checked with blocks_3d commit `58c97173` and the SDPB
3.1.0 Mathematica writer. Actual generated scalar, Majorana, and nonconserved-vector
tables were assembled into PMP files. Majorana's two-dimensional parity-even
OPE sector appears at exchanged spin 2, so this exercises matrix positivity,
not only scalar inequalities. Increasing recursion/pole orders from `12/6` to
`24/12` also checks that omitted-component residuals decrease. These are tests
of conventions and numerical integration, not a reproduction of a published
exclusion bound. An SDPB optimization run has not been validated here.

`spinIndependentRows[sdp,x,tolerance]` provides an explicitly approximate row
selection for exploration with other identical spins. Its tolerance is part of
the approximation; near-dependence from truncated blocks can otherwise become a
spurious constraint. Use analytic independent structures where available. For
Majorana, use `spinMajoranaRows`, not numerical rank detection.

Before interpreting any numerical exclusion one must include the required
protected sectors and verify convergence in exchanged spin, derivatives,
recursion depth, retained poles, and precision. This package does not infer an
infinite-spin tail, supply Ward identities, or certify a bound from one finite
truncation. It does not automatically combine the new spacetime structures with
the existing global-symmetry `bootAll[]` pipeline.
