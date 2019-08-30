# autoboot

Automatical Generator of Conformal Bootstrap Equation

For more information, see
[autoboot: A generator of bootstrap equations with global symmetry](https://arxiv.org/abs/1903.10522).

Some usages also can be checked by typing `?someSymbolName` (for example, `?getGroup`) in [Mathematica](http://reference.wolfram.com/language/tutorial/GettingInformationAboutWolframLanguageObjects.html).
``?somePackageName`*`` (for example, ``?ClebschGordan`*``) will give package-lebel information.

## Setup

```sh
tar xvf sgd.tar.xz
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
$lie_group = su[2] | so[2] | so[3] | o[2] | o[3] | su[4]
$n = positive_integer
```

Once you get a group `G`:

1. Set `G` as a global symmetry by `setGroup[G]` (if you call this more than once, all values calculated by autoboot previously will be cleared).

1. Register fundamental operators by `setOps[...]`. 'Fundamental' means that the operator in summation of all primary operator will be treated independently.

1. Get bootstrap equations by `bootAll[]`. If you need human-readable format, use `format[...]`.

1. *(Optionally)* You can get a Python code for [cboot](https://github.com/tohtsky/cboot) by `toCboot[makeSDP[eq]]`.

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
