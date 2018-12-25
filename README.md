# autoboot

Automatical Generator of Conformal Bootstrap Equation

## Usage

To use autoboot, you have to load either `group.m` or `ngroup.m` (not both).
`ngroup.m` requires just one call of `setPrecision`, which controls the precision of calculation in autoboot.

If you load `group.m`, all numerical values are rigorous and it takes much time in calculation (in some cases, Mathematica will freeze).

If you load `ngroup.m`, all numerical values (except for signs, multiplicities and so on) are approximated and it takes much less time.

### Groups

autoboot supports many finite groups, some Lie groups and product groups of them.
More formally, groups which we can treat as a global symmetry of CFT are defined by:

```EBNF
group = finite_group | lie_group | "pGroup[", group , "," , group , "]";
finite_group = "group[", n, ",", n, "]" | "dih[", n, "]" | "dic[", n, "]";
lie_group = "su[2]" | "so[2]" | "so[3]" | "o[2]" | "o[3]";
n = positive_integer;
```

## Example

This example generates a bootstrap equation of D8-symmetric CFT.
You can execute this in `sample.nb`.

```Mathematica
SetDirectory["~/autoboot/"];
Block[{$CharacterEncoding = "UTF-8"}, << "group.m"]
d8 = getGroup[8, 3]; (* getGrou[8, 3] == getDihedral[4] *)
setGroup[d8];
setOps[{op[e, d8[id], 1, 1], op[v, rep[5], 1, 1]}]
format[eq = bootAll[]]
ans = makeSDP[eq];
buf = OpenWrite["d8.py"];
WriteString[buf, toString[ans]]
Close[buf];
```
