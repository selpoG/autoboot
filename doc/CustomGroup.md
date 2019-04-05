# How to Make a Custom Group-Object

Any custom group-object `g` must have attributes: `id, dim, prod, dual, isrep, gG, gA, minrep`.
You also attach labels of irrep-objects for all (unique) irreps.

[`getO[3]`](https://github.com/selpoG/autoboot/tree/master/grouplie.m) is
a good example for creating custom group-objects.

## `isrep`

`g[isrep[r]]` must give `True` or `False` (`True` if and only if `r` is a valid irrep-object for `g`).

## `minrep`

`g[minrep[r1,r2]]` must give `r1` or `r2` for irrep-objects `r1` and `r2`.
The set of all irrep-objects are assumed to be totally ordered by `minrep`.

## `dim`

`g[dim[r]]` must give the dimension of irrep-object `r`.
`Head[g[dim[r]]]===Integer` is assumed.

## `dual`

`g[dual[r]]` must give irrep-object for the dual irrep of irrep-object `r`.
`g[dual[g[dual[r]]]===r` is assumed.
If `r` is pseudo-real or strictly real, `g[dual[r]]===r` is also assumed.

## `id`

`g[id]` must give irrep-object for trivial irrep.
`g[dim[g[id]]]===1` and `g[dual[g[id]]]===g[id]` are assumed.

## `prod`

`g[prod[r1,r2]]` must give a list of all irrep-objects with multiplicity
arising in irreducible decomposition of direct product representation of
`r1` and `r2` for irrep-objects `r1` and `r2`.
`Head[g[prod[r,s]]]===List`, `And[g[minrep[#]]&/@g[prod[r,s]]]===True` and
`g[prod[g[id],r]]===g[prod[r,g[id]]]==={r}` are assumed.

## `gG` and `gA`

`g[gG]` must give a list of labels for generators of discrete part of `g`.
If `g[gG]={a,b,c,...}` and `e` is the identity of `g`,
discrete part of `g` must be `{e,a,b,c,aa,ab,ac,ba,bb,bc,ca,cb,cc,...}`.
`Head[g[gG]]===List` is assumed.

`g[gA]` must give a list of labels for generators of lie-algebra part of `g`.
If `g[gA]={X,Y,Z,...}`, lie-algebra of `g` over complex numbers must be spanned by
`{X,Y,Z,[X,Y],[X,Z],[Y,Z],...}`.
`Head[g[gA]]===List` is assumed.

These labels must give explicit matrix representation for irrep-object of `g`.
All irreps must be unitary.

If `a` is an element of `g[gG]`, `a[r]` must be a unitary matrix which represents `a` in `r`.
If `X` is an element of `g[gA]`, `X[r]` must be a matrix which represents `X` in `r`
(unitary condition for `X` is not obvious; `X` may not be an element of lie-algebra of `g` over real numbers).
If `p` is an element of `g[gG]` or `g[gA]`,
`Head[p[r]]===List` and `And[Head[#]===List&&Length[p[r]]===Length[#]&/@p[r]]` are assumed.
Each elements `v` of the matrix must satisfy `ExactNumberQ[v]===True` for `inv.m`
or `NumericQ[v]===True` for `ninv.m`.

These matrices are used for write down the `g`-invariance conditions for `ope[r,s,t][n][a,b,c]`.
Let `f[a,b,c]` is a tensor in which `a` transforms as an indice of `r`, `b` of `s` and `c` of `g[dual[t]]`.
We need `g`-invariant subspace, so we have:

- `Sum[f[a,b,c2] p[t][[c2,c]],{c2,l}] == Sum[p[r][[a,a2]] p[s][[b,b2]] f[a2,b2,c],{a2,n},{b2,m}]`
  for `p` in `g[gG]`,

- `Sum[f[a,b,c2] p[t][[c2,c]],{c2,l}] == Sum[p[r][[a,a2]] f[a2,b,c],{a2,n}] + Sum[p[s][[b,b2]] f[a,b2,c],{b2,m}]`
  for `p` in `g[gA]`.

Here we used `n:=g[dim[r]]`, `m:=g[dim[s]]` and `l:=g[dim[t]]`.
