# autoboot

Automatical Generator of Conformal Bootstrap Equation

For more information, see
[autoboot: A generator of bootstrap equations with global symmetry](https://arxiv.org/abs/1903.10522),
or [An Automated Generation of Bootstrap Equations for Numerical Study of Critical Phenomena](https://arxiv.org/abs/2006.04173).

Some usages also can be checked by typing `?someSymbolName` (for example, `?getGroup`) in [Mathematica](http://reference.wolfram.com/language/tutorial/GettingInformationAboutWolframLanguageObjects.html).
``?somePackageName`*`` (for example, ``?ClebschGordan`*``) will give package-lebel information.
These usages are also documented in [inv.md](doc/inv.md), [group.md](doc/group.md).

- [Setup](#setup)
- [Testing](#testing)
- [Usage](#usage)
- [Example](#example)

## Setup

```sh
tar xvf sgd.tar.xz
```

## Testing

Run the complete test suite with:

```sh
make test
```

This runs both the license-free checks and the Wolfram Language regression
tests. The latter require an installed and activated `wolframscript`.

GitHub Actions runs only the license-free part so that CI does not require a
Wolfram license or entitlement:

```sh
make test-free
```

To run only the tests that exercise autoboot with the Wolfram Engine:

```sh
make test-wolfram
```

## Usage

To use autoboot, you have to load either `group.m` or `ngroup.m` (**not both**).
`ngroup.m` requires **just one call** of `setPrecision`, which controls the precision of calculation in autoboot.

If you load `group.m`, all numerical values are rigorous and it takes much time in calculation (in some cases, Mathematica will freeze).

If you load `ngroup.m`, all numerical values (except for signs, multiplicities and so on) are approximated and it takes much less time.

### Groups

autoboot supports many **finite groups**, some **Lie groups** and **product groups** of them.
More formally, groups which we can treat as a global symmetry of CFT are defined by:

```EBNF
$group = $finite_group | $lie_group | pGroup[$group,$group]
$finite_group = group[$n,$n] | dih[$n] | dic[$n]
$lie_group = su[$degree] | so[$degree] | spin[$spin_degree] | sp[$n] | lie[$type,$rank] | o[$degree]
$spin_degree = integer_at_least_4
$type,$rank = supported_finite_Cartan_type_and_rank
$degree = integer_at_least_2
$n = positive_integer
```

Once you get a group `G`:

1. Set `G` as a global symmetry by `setGroup[G]` (if you call this more than once, all values calculated by autoboot previously will be cleared).

1. Register fundamental operators by `setOps[...]`. 'Fundamental' means that the operator in summation of all primary operator will be treated independently.

1. Get bootstrap equations by `bootAll[]`. If you need human-readable format, use `format[...]`.

1. *(Optionally)* You can get a Python code for [cboot](https://github.com/tohtsky/cboot) by `toCboot[makeSDP[eq]]`.

For `getSU[n]`, `n` is the degree (the Lie algebra rank is `n-1`). Both exact
and numerical modes support all integer `n >= 2`. SU(2) retains spin labels;
SU(n) for `n >= 3` uses Young diagram row lengths `v[l1,...,l(n-1)]`.
Generator matrices are built lazily from tensor products of exterior powers.
Large highest weights can require substantial memory and time in either mode;
the numerical mode converts the shared exact matrices to the configured precision.

SO(n), Spin(n), compact Sp(n), and the exceptional groups are also supported
in both modes. For example, `getSO[5]`, `getSpin[6]`, `getSp[2]`, and
`getLie["G",2]` return group objects. `getLie[type,rank]` supports all finite
Cartan types A–G, including E6, E7 and E8. These new constructors use Dynkin
labels; SO rejects spinor labels that do not descend from Spin.
See [Lie groups](doc/LieGroups.md) for conventions, supported ranks and examples.

Reproducible bootstrap benchmarks and end-to-end checks are described in
[Lie group timing and validation](doc/LieGroups.md#end-to-end-bootstrap-checks-and-timing).

### Irreps

Please see [IrrepLabels.md](/doc/IrrepLabels.md).

### Custom Group-Object

Please see [CustomGroup.md](/doc/CustomGroup.md).

## Example

This example generates a bootstrap equation of D8-symmetric CFT.
More example codes can be found in [`sample` folder](/sample).

```Mathematica
(* change the path of autoboot properly *)
SetDirectory["~/autoboot/"];
<< "group.m"
(* getGroup[8, 3] is isomorphic to getDihedral[4] *)
d8 = getGroup[8, 3];
setGroup[d8];
(* set e and v as fundamental operators. rep[5] is the unique 2-dim irrep of d8 (you can check this by d8[ct]). *)
setOps[{op[e, d8[id], 1, 1], op[v, rep[5], 1, 1]}]
format[eq = bootAll[]]
ans = makeSDP[eq];
WriteString["d8.py", toCboot[ans]]
(* specify how to convert operators to latex code *)
opToTeX[e] := "\\epsilon"
opToTeX[v] := "v"
(* specify how to convert irreps to latex code *)
repToTeX[rep[n_]] := TemplateApply["\\mathbf{`n`}", <|"n" -> n|>]
(* you can paste printed string to your latex file *)
Print[toTeX[eq]]
```

Larger classical-group bootstrap examples and independent GAP / CG reference
checks are documented in [doc/LargeValidation.md](doc/LargeValidation.md).
Run `make test-large` for the extended exact/numerical integration suite.

General O(n), including the disconnected reflection component, is described
in [doc/OrthogonalGroups.md](doc/OrthogonalGroups.md).
