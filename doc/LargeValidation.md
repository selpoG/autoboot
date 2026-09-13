# Larger representations and independent validation

The extended suite runs a fresh Wolfram kernel for each group and arithmetic
mode. It generates the complete bootstrap system and constructs the symbolic
SDP object. It does not run a numerical SDP optimization or find CFT bounds.

```sh
make test                       # Includes the external reference checks
make test-large                 # Extended exact + 50-digit numerical checks
# Or retain equations and logs in a chosen directory:
./benchmark/validate-large-suite.sh /tmp/autoboot-large-results
python3 test/reference/generate-gap.py  # Regenerate the independent GAP data
```

Each extended case has a 180-second external hard limit, configurable with
`BOOTSTRAP_TIMEOUT_SECONDS`. Failure, missing references, and timeouts produce
nonzero exit codes. Exact equations are saved first; the numerical run compares
coefficient row spaces against that exact result, with tolerance `10^-25`.
Startup and `setOps` are excluded from `BOOT_SECONDS`. `TOTAL_SECONDS` also
includes registration and the subsequent component tests, but excludes startup.

## What is checked independently

* **Dimensions, weight multiplicities, tensor decomposition multiplicities:**
  15 fixtures generated with GAP 4.11.1, using its Lie-algebra implementation.
  Examples include SU(3) `27 x 8`, SU(4) `15 x 15`, SU(5) `24 x 24`,
  SO(5) `14 x 14`, SO(7) `21 x 21`, SO(11) `11 x 11`,
  Sp(4) `27 x 27`, Sp(5) `10 x 10`, Spin(10) `16 x 16`,
  SO(12) `12 x 12`, E6 `27 x 27`, E7 `56 x 56`, E8 `248 x 248`,
  F4 `26 x 26`, and G2 `14 x 14`.
  Sp rank conventions are those of `getSp`: Sp(4) has an 8-dimensional
  defining representation. Spin(10)'s chiral spinor product distinguishes
  the two conjugate 126-dimensional representations by Dynkin label.
  See [GAP section 64.13](https://gap-system.github.io/gap/doc/ref/chap64.html).
* **Individual CG components:** `test/CGExternal.wls` compares the general
  highest-weight backend for SO(3) spins `10 x 8 -> 2, 9, 18` against
  [Wolfram's built-in ClebschGordan](https://reference.wolfram.com/language/ref/ClebschGordan.html).
  Only one overall phase per embedding may differ. The full arrays include
  1,785, 6,783, and 13,209 components respectively, including all required zeros.
* **Higher-rank CG subspaces:** SO(7), SO(10), SO(12) vector CGs are checked
  against the analytic projectors
  `P_S = |J><J|/N`, `P_A = (I-Swap)/2`,
  `P_T = (I+Swap)/2-P_S`, where `J` is the invariant metric in the
  representation's basis. The trace metric necessarily comes from the
  representation basis; the symmetric/antisymmetric projectors are independent
  of the CG solver. These are the normalized singlet, antisymmetric and
  symmetric-traceless tensor structures of
  [Kos, Poland, Simmons-Duffin, arXiv:1307.6856, eq. (2.2)](https://arxiv.org/abs/1307.6856),
  with index orientation and normalization adjusted for Hermitian projectors.

The GAP fixtures do **not** contain CG coefficients. We do not claim direct
component comparisons against an external CG library for every Lie type.
These checks cover the stated representations, not arbitrary rank or labels.
F4's defining-26 bootstrap is also included in the extended integration tests.

## Extended integration checks

For each case, the CG matrices for **all** tensor-product channels and
multiplicity copies are joined into a square change-of-basis matrix. The
suite verifies its unitarity and every raising/lowering intertwining equation.
It then reconstructs every component of all four-index invariant tensors,
including both multiplicity indices, and checks the crossing matrix and its
involution. SU(4)'s adjoint square includes two copies of the adjoint and
conjugate complex channels; Sp(4) exercises the pseudoreal external operator
and its automatically registered dual.

These structural checks and the exact/numerical equation comparison supplement
the independent reference checks above. They are not described as independent
external calculations of the whole bootstrap system.

## Performance

Wolfram Engine 15.0, Intel Core i5-13500H, fresh kernels and 50-digit
numerical precision. The table reports `bootAll[]` only; startup,
registration and validation are excluded. Timings depend on machine load.


Fresh kernels on the same machine, with all component checks and exact/numerical
row-space comparisons passing (14 classical runs):

| External representation | Exact | Numerical | Equations |
| --- | ---: | ---: | ---: |
| SO(7) vector 7 | 0.19 s | 0.11 s | 3 |
| SO(10) vector 10 | 0.50 s | 0.29 s | 3 |
| SO(12) vector 12 | 1.07 s | 0.63 s | 3 |
| SO(5) adjoint 10 | 1.25 s | 0.75 s | 6 |
| SO(5) symmetric traceless 14 | 3.08 s | 1.91 s | 6 |
| SU(4) adjoint 15 | 2.57 s | 1.32 s | 6 |
| Sp(4) defining 8, with pseudoreal dual | 1.47 s | 0.98 s | 21 |
| F4 defining 26 | 37.1 s | 33.9 s | 5 |

F4 full-component validation took 116.5 seconds in exact mode on this machine.
The SU(4) example contains 9 four-point invariant tensors and 6 equations;
Sp(4) includes the automatically registered pseudoreal dual. These measurements
do not establish performance for arbitrary highest weights or E6/E7/E8 bootstrap.

## Validation limits

The saved bootstrap systems are regression data, not independently published
answers. Both arithmetic modes share code, so their agreement alone is not
an independent proof. Tests cover CG covariance/completeness, tensor crossing
and construction of the SDP object. Direct comparisons with published final
sum rules, an independent end-to-end equation-reduction implementation, and
independent PSD-normalization checks remain future work.
