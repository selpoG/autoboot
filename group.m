Needs["CommonFunctions`", "common.m"]

BeginPackage["GroupInfo`"]

getGroup::usage = "getGroup[g,i]"
product::usage = "product[g1,g2]"
group::usage = "group[g,i]"
pGroup::usage = "pGroup[g1,g2]"
ncg::usage = "ncg"
ct::usage = "ct[g]"
id::usage = "id[g]"
dim::usage = "dim[r]"
prod::usage = "prod[r,s]"
dual::usage = "dual[r]"
isrep::usage = "isrep[r]"
gG::usage = "gG"
gA::usage = "gA"
rep::usage = "rep[n]"
v::usage = "v[n]
v[n,s]"
i::usage = "i[a,b]"
minrep::usage = "minrep[r,s]"
setGroup::usage = "setGroup[G]"
available::usage = "available[g,i]"

Begin["`Private`"]

e = CommonFunctions`e
MyReap = CommonFunctions`MyReap
dual /: G_[dual[x_List]] := G[dual[#]] & /@ x;

(*setGroup[G_] := (<<"cg.m"; ClebschGordan`Private`setGroup[G];)*)
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
	G[dual[rep[r1_, r2_]]] := rep[g1[dual[r1]], g2[dual[r2]]];
	G[isrep[_]] := False;
	G[isrep[rep[r1_, r2_]]] := g1[isrep[r1]] && g2[isrep[r2]];
	G[gG] = Join[G[1, #]& /@ g1[gG], G[2, #]& /@ g2[gG]];
	G[gA] = Join[G[1, #]& /@ g1[gA], G[2, #]& /@ g2[gA]];
	G[1, g_][rep[r1_, r2_]] := G[1, g][rep[r1, r2]] = KroneckerProduct[g[r1], IdentityMatrix[g2[dim[r2]]]];
	G[2, g_][rep[r1_, r2_]] := G[2, g][rep[r1, r2]] = KroneckerProduct[IdentityMatrix[g1[dim[r1]]], g[r2]];
	G[minrep[rep[r1_, r2_], rep[s1_, s2_]]] := rep[g1[minrep[r1, s1]], g2[minrep[r2, s2]]];
	G
]

End[ ]

EndPackage[ ]

<< "groupd.m"
<< "grouplie.m"
