# autoboot

Automatical Generator of Conformal Bootstrap Equation

## group.m

### Example

Example of generation of bootstrap equation of D8-symmetric CFT.
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
