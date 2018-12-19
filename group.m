Needs["CommonFunctions`", "common.m"]

BeginPackage["GroupInfo`"]

(* Supported groups are direct products of finite groups and some lie groups. *)
(* Supported finite groups are those whose irreps was calculated by GAP (in sgd folder), and any dihedral and quartenion groups. *)
(* Supported lie groups are su[2], o[2], so[2]. We support only compact groups, so we can assume any finite dim. irrep can be unitarized. *)
getGroup::usage = "getGroup[g,i] loads data from sg.g.i.m and returns group-object group[g,i]. g is the order of the finite group, i is the number of the group assined by GAP."
product::usage = "product[g1,g2] returns group-object pGroup[g1,g2] which represents direct product of two group-object g1, g2."
group::usage = "group[g,i] is a group-object whose order is g and whose number assigned by GAP is i. Before using this value, you have to call getGroup[g,i] to get proper group-object."
pGroup::usage = "pGroup[g1,g2] is a group-object which is a direct product of g1, g2. Before using this value, you have to call product[g1,g2] to get proper group-object."
(* group-object g has attributes ncg, ct, id, dim, prod, dual, isrep, gG, gA, minrep. *)
(* You can evaluate attributes in putting it in g[...]. *)
(* For example, g[dim[r]] gives the dimension of irrep r. *)
ncg::usage = "ncg is the number of conjugacy classes, which is also the number of inequivalent irreps. This is not defined for lie groups."
ct::usage = "ct is the character table. This is not defined for lie groups."
id::usage = "id is the trivial representation."
dim::usage = "dim[r] is the dimension of irrep r."
prod::usage = "prod[r,s] gives a list of all irreps arising in irreducible decomposition of direct product representation of r and s. prod[r,s] may not be duplicate-free."
dual::usage = "dual[r] gives dual representation of irrep r."
isrep::usage = "isrep[r] gives whether r is recognised as a irrep-object of the group-object or not."
gG::usage = "gG is a list of all generator-objects of finite group part of the group-object."
gA::usage = "gA is a list of all generator-objects of lie algebra part of the group-object."
rep::usage = "rep[n] is n-th irrep-object (n is assined by GAP and corresponds to the index of ct). This is recognised only by group[g,i].
rep[r1,r2] is natural irrep-object of pGroup[g1,g2] where r1 is irrep-object of g1, r2 is irrep-object of g2. This is recognised only by pGroup[g1,g2]."
v::usage = "v[n] is spin-n irrep-object. This is recognised only by dih[n], dic[n], su[2], so[3], o[2] and so[2].
v[n,s] is spin-n irrep-object with sign s. This is recognised only by o[3]."
i::usage = "i[a] is one-dimensional irrep-object with sign a. This is recognised only by dih[n] (n: odd) and o[2].
i[a,b] is one-dimensional irrep-object with sign a,b. This is recognised only by dih[n] (n:even), dic[n]."
(* We need all irreps to be sorted in some linear order. *)
minrep::usage = "minrep[r,s] gives r if r < s else s. r and s are irrep-objects."
setGroup::usage = "setGroup[G] loads inv.m with global symmetry G. This action clears all values calculated by inv.m previously."
available::usage = "available[g,i] gives whether group[g,i] are supported or not."

Begin["`Private`"]

e = CommonFunctions`e
MyReap = CommonFunctions`MyReap
dual /: G_[dual[x_List]] := G[dual[#]] & /@ x;

setGroup[G_] := (Block[{$CharacterEncoding = "UTF-8"}, <<"inv.m"]; ClebschGordan`Private`setGroup[G];)

getGroup[1, 1] := getGroup[1, 1] = AbortProtect @ Module[{G},
	G = group[1, 1];
	G[id] = rep[1];
	G[ncg] = 1;
	G[ct] = {{1}};
	G[dim[rep[1]]] = 1;
	G[prod[rep[1], rep[1]]] = {rep[1]};
	G[dual[rep[1]]] = rep[1];
	G[isrep[_]] := False;
	G[isrep[rep[1]]] = True;
	G[gG] = {};
	G[gA] = {};
	G[minrep[rep[1], rep[1]]] = rep[1];
	G
]

available[1, 1] = True
x:available[g_, i_] := x = FileExistsQ[TemplateApply["sgd/sg.`g`.`i`.m", <|"g" -> g, "i" -> i|>]]

getGroup[g_, i_] := getGroup[g, i] = Module[{dat, ncg$, ct$, rep$, mul$, ip, file, G, r, s, a, b, mats, j, n, m, myprod, z},
	file = TemplateApply["sgd/sg.`g`.`i`.m", <|"g" -> g, "i" -> i|>];
	If[! FileExistsQ[file], Return[Null]];
	dat = FullSimplify[Import[file]];
	If[dat[[1]] != {g, i}, Return[Null]];
	ncg$ = Length[mul$ = dat[[2]]];
	ct$ = Expand[dat[[3]]];
	rep$ = FullSimplify @ Expand[dat[[5]]];
	ip[r_, s_] := Simplify[Total[Conjugate[r] s mul$] / g];
	z:myprod[rep[a_]] := z = Module[{x}, MyReap[Do[Do[Sow[rep[x]], N[FullSimplify[ip[ct$[[x]], ct$[[a]]^2]], 200]], {x, ncg$}]]];
	AbortProtect[
		G = group[g, i];
		G[id] = rep[1];
		G[ncg] = ncg$;
		G[ct] = ct$;
		G[dim[rep[a_]]] := G[dim[rep[a]]] = ct$[[a, 1]];
		G[prod[rep[a_], rep[a_]]] := myprod[rep[a]];
		G[prod[rep[a_], rep[b_]]] := G[prod[rep[a], rep[b]]] =
			Module[{x}, MyReap[Do[Do[Sow[rep[x]], N[FullSimplify[ip[ct$[[x]], ct$[[a]] ct$[[b]]]], 200]], {x, ncg$}]]];
		G[dual[rep[a_]]] := G[dual[rep[a]]] =
			Module[{r, x}, r = Conjugate[ct$[[a]]]; Catch[Do[If[ip[ct$[[x]], r] != 0, Throw[rep[x]]], {x, ncg$}]]];
		G[isrep[_]] := False;
		G[isrep[rep[a_]]] := 1 <= a <= G[ncg];
		G[gA] = {};
		G[gG] = ToExpression[TemplateApply["GroupInfo`a`group[`g`, `i`][GroupInfo`a`Private`a``x`]", <|"a" -> "`", "g" -> g, "i" -> i, "x" -> #|>]] & /@ dat[[4]];
		set[r:rep[_], mats_] := Module[{j}, Do[Evaluate[G[gG][[j]][r]] = mats[[j]], {j, Length[G[gG]]}]];
		set[mats_] := Do[set[rep[j], mats[[j]]], {j, Length[mats]}];
		set[rep$];
		G[minrep[rep[n_], rep[m_]]] := rep[Min[n, m]];
		G
	]
]

product[g1_, g2_] := AbortProtect @ Module[{G, r1, r2, s1, s2, g},
	G = pGroup[g1, g2];
	G[id] = rep[g1[id], g2[id]];
	G[dim[rep[r1_, r2_]]] := g1[dim[r1]] g2[dim[r2]];
	G[prod[rep[r1_, r2_], rep[s1_, s2_]]] := G[prod[rep[r1, r2], rep[s1, s2]]] = Module[{T1, T2, t1, t2},
		T1 = g1[prod[r1, s1]];
		T2 = g2[prod[r2, s2]];
		MyReap @ Do[Sow[rep[t1, t2]], {t1, T1}, {t2, T2}]
	];
	G[dual[rep[r1_, r2_]]] := rep[g1[dual[r1]], g2[dual[r2]]];, su[2], so[3], o[2], so[2]
	G[isrep[_]] := False;, su[2], so[3], o[2], so[2]
	G[isrep[rep[r1_, r2_]]] := g1[isrep[r1]] && g2[isrep[r2]];, su[2], so[3], o[2], so[2]
	G[gG] = Join[G[1, #]& /@ g1[gG], G[2, #]& /@ g2[gG]];, su[2], so[3], o[2], so[2]
	G[gA] = Join[G[1, #]& /@ g1[gA], G[2, #]& /@ g2[gA]];, su[2], so[3], o[2], so[2]
	G[1, g_][rep[r1_, r2_]] := G[1, g][rep[r1, r2]] = Kronecke, su[2], so[3], o[2], so[2]rProduct[g[r1], IdentityMatrix[g2[dim[r2]]]];
	G[2, g_][rep[r1_, r2_]] := G[2, g][rep[r1, r2]] = Kronecke, su[2], so[3], o[2], so[2]rProduct[IdentityMatrix[g1[dim[r1]]], g[r2]];
	G[minrep[rep[r1_, r2_], rep[s1_, s2_]]] := rep[g1[minrep[r, su[2], so[3], o[2], so[2]1, s1]], g2[minrep[r2, s2]]];
	G
]

End[ ]

EndPackage[ ]

<< "groupd.m"
<< "grouplie.m"
