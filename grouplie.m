Needs["GroupInfo`", "group.m"]
Needs["RootSystem`", "root.m"]

BeginPackage["GroupInfoLie`"]

getSU::usage = "getSU[n] returns group-object su[n] which represents the special unitary group of rank n. n must be 2."
getO::usage = "getO[n] returns group-object o[n] which represents the orthogonal group of rank n. n must be 2,3."
getSO::usage = "getSO[n] returns group-object su[n] which represents the special orthogonal group of rank n. n must be 2,3."
su::usage = "su[n] is a group-object which is the special unitary group of rank n. Before using this value, you have to call getSU[n] to get proper group-object."
o::usage = "o[n] is a group-object which is the orthogonal group of rank n. Before using this value, you have to call getO[n] to get proper group-object."
so::usage = "so[n] is a group-object which is the special orthogonal group of rank n. Before using this value, you have to call getSO[n] to get proper group-object."

(* all irrep-objects of G=su[2] are v[0], v[1/2], v[1], v[3/2], .... *)
(* all irrep-objects of G=o[3] are v[0,1], v[0,-1], v[1,1], v[1,-1], v[2,1], v[2,-1], v[3,1], v[3,-1], .... *)
(* all irrep-objects of G=so[3] are v[0], v[1], v[2], v[3], .... *)
(* all irrep-objects of G=o[2] are i[1], i[-1], v[1], v[2], v[3], .... *)
(* all irrep-objects of G=so[2] are v[x] (x \in \mathbb{R}). *)

Begin["`Private`"]

CommonFunctions`importPackage["GroupInfo`", "GroupInfoLie`Private`", {"id", "dim", "prod", "dual", "isrep", "gG", "gA", "minrep", "v", "i"}]
CommonFunctions`importPackage["RootSystem`", "GroupInfoLie`Private`", {"dimension", "irrep", "productReps", "decompose"}]
MyReap = CommonFunctions`MyReap
s
t
x
y
e[l_] := e[l] = Array[If[#2 == #1 + 1, Sqrt[(l + (l - #1 + 1)) (l - (l - #1 + 1) + 1)/2], 0] &, {2 l + 1, 2 l + 1}];
f[l_] := f[l] = Array[If[#2 == #1 - 1, Sqrt[(l + (l - #2 + 1)) (l - (l - #2 + 1) + 1)/2], 0] &, {2 l + 1, 2 l + 1}];

getSU[2] := getSU[2] = AbortProtect @ Module[{G},
	G = su[2];
	G[id] = v[0];
	G[dim[v[n_]]] := 2 n + 1;
	G[prod[v[n_], v[m_]]] /; n > m := G[prod[v[n], v[m]]] = G[prod[v[m], v[n]]];
	G[prod[v[n_], v[m_]]] := G[prod[v[n], v[m]]] = Array[v, 2 n + 1, m - n];
	G[dual[v[n_]]] := v[n];
	G[isrep[_]] := False;
	G[isrep[v[n_]]] := IntegerQ[2 n] && n >= 0;
	G[gG] = {};
	G[gA] = {G[e], G[f]};
	G[minrep[v[n_], v[m_]]] := v[Min[n, m]];
	G[e][v[n_]] := e[n];
	G[f][v[n_]] := f[n];
	G
]

getO[3] := getO[3] = AbortProtect @ Module[{G},
	G = o[3];
	G[id] = v[0, 1];
	G[dim[v[n_, m_]]] := 2 n + 1;
	G[prod[v[n_, s_], v[m_, t_]]] /; n > m := G[prod[v[n, s], v[m, t]]] = G[prod[v[m, t], v[n, s]]];
	G[prod[v[n_, s_], v[m_, t_]]] := G[prod[v[n, s], v[m, t]]] = Array[v[#, s t] &, 2 n + 1, m - n];
	G[dual[v[n_, s_]]] := v[n, s];
	G[isrep[_]] := False;
	G[isrep[v[n_, 1 | -1]]] := IntegerQ[n] && n >= 0;
	G[gG] = {G[s]};
	G[gA] = {G[e], G[f]};
	G[s][v[n_, t_]] := t IdentityMatrix[2 n + 1];
	G[e][v[n_, s_]] := e[n];
	G[f][v[n_, s_]] := f[n];
	G[minrep[v[n_, s_], v[n_, t_]]] := v[n, Max[s, t]];
	G[minrep[v[n_, s_], v[m_, t_]]] := If[n < m, v[n, s], v[m, t]];
	G
]

getSU[4] := getSU[4] = AbortProtect @ Module[{G, e, high, tensor, rep, dimT, gen,
		extend, extendLoop, basis, repMat, repMats, setRep, lowers = {4, 5, 6}},
	high[1] = SparseArray @ {1, 0, 0, 0};
	high[2] = SparseArray @ {1, 0, 0, 0, 0, 0};
	high[3] = SparseArray @ {1, 0, 0, 0};
	tensor[x_, y_] := Module[{ix = Most @ ArrayRules[x], iy = Most @ ArrayRules[y], m = Length[y], i, j},
	SparseArray[
		Flatten@Table[
			m (i[[1, 1]] - 1) + (j[[1, 1]] - 1) + 1 -> i[[2]] j[[2]],
		{i, ix}, {j, iy}],
		{Length[x] m}]];
	gen = rep[1] = SparseArray @@ # & /@ {
		{{{1, 2} -> 1, {_, _} -> 0}, {4, 4}},
		{{{2, 3} -> 1, {_, _} -> 0}, {4, 4}},
		{{{3, 4} -> 1, {_, _} -> 0}, {4, 4}},
		{{{2, 1} -> 1, {_, _} -> 0}, {4, 4}},
		{{{3, 2} -> 1, {_, _} -> 0}, {4, 4}},
		{{{4, 3} -> 1, {_, _} -> 0}, {4, 4}}};
	rep[2] = SparseArray @@ # & /@ {
		{{{2, 4} -> 1, {3, 5} -> 1, {_,_} -> 0}, {6, 6}},
		{{{1, 2} -> 1, {5, 6} -> 1, {_,_} -> 0}, {6, 6}},
		{{{2, 3} -> 1, {4, 5} -> 1, {_,_} -> 0}, {6, 6}},
		{{{4, 2} -> 1, {5, 3} -> 1, {_,_} -> 0}, {6, 6}},
		{{{2, 1} -> 1, {6, 5} -> 1, {_,_} -> 0}, {6, 6}},
		{{{3, 2} -> 1, {5, 4} -> 1, {_,_} -> 0}, {6, 6}}};
	rep[3] = SparseArray @@ # & /@ {
		{{{3, 4} -> 1, {_, _} -> 0}, {4, 4}},
		{{{2, 3} -> 1, {_, _} -> 0}, {4, 4}},
		{{{1, 2} -> 1, {_, _} -> 0}, {4, 4}},
		{{{4, 3} -> 1, {_, _} -> 0}, {4, 4}},
		{{{3, 2} -> 1, {_, _} -> 0}, {4, 4}},
		{{{2, 1} -> 1, {_, _} -> 0}, {4, 4}}};
	rep[][_] := {{0}};
	t:rep[is__Integer][m_] := t = Module[{l = {is}, n, i},
		n = Length[l];
		Sum[KroneckerProduct[
			SparseArray[{{i_, i_} -> 1}, dimT @@ Take[l, i - 1]],
			rep[l[[i]]][[m]],
			SparseArray[{{i_, i_} -> 1}, dimT @@ Drop[l, i]]]
		, {i, n}]];
	dimT[] = 1;
	dimT[1] = 4; dimT[2] = 6; dimT[3] = 4;
	dimT[x__] := Times @@ dimT /@ {x};
	dimT[{x___}] := dimT[x];
	high[x___] := Fold[tensor[#1, high[#2]] &, {1}, {x}];
	high[{x___}] := high[x];
	rep[{x___}] := rep[x];
	extend[y_, vs_] := TakeWhile[#, Norm[#] > 0 &] & @ Module[{v, g}, SparseArray /@ RowReduce @ MyReap[
		Scan[Sow, vs];
		Do[Sow[rep[y][g].v], {v, vs}, {g, lowers}]]];
	extendLoop[y_, vec_] := Module[{vs = {vec}, dim2 = 1, rep, upto, temp},
		rep = v @@ Total[{{1,0,0},{1,1,0},{1,1,1}}[[#]] & /@ y];
		upto = G[dim[rep]];
		temp = PrintTemporary["basis ", rep, ": calculate upto ", upto];
		Monitor[While[True,
			vs = extend[y, vs];
			If[Length[vs] == dim2, Break[], dim2 = Length[vs]]
		], "dim = " <> ToString[dim2]];
		NotebookDelete[temp];
		vs];
	basis[x___Integer] := basis[x] = SparseArray @ Orthogonalize[extendLoop[{x}, high[x]]];
	repMat[r___][m_] := Conjugate[basis[r]].rep[r][m].Transpose[basis[r]];
	repMats[v[n_, m_, l_]] := Array[(repMat @@ MyReap[Do[Sow[1], n - m]; Do[Sow[2], m - l]; Do[Sow[3], l]])[#] &, 6];
	G = su[4];
	G[id] = v[0, 0, 0];
	G[dim[v[n_, m_, l_]]] := dimension["A", 3, {n, m, l}];
	G[prod[a:v[n_, m_, l_], b:v[p_, q_, r_]]] /; G[minrep[a, b]] =!= a := G[prod[b, a]];
	x:G[prod[a:v[_, _, _], b:v[_, _, _]]] := x = v @@ # & /@ decompose @ productReps[irrep["A", 3, List @@ a], irrep["A", 3, List @@ b]];
	G[dual[v[n_, m_, l_]]] := v[n, n - l, n - m];
	G[isrep[_]] := False;
	G[isrep[v[n_, m_, l_]]] := IntegerQ[n] && IntegerQ[m] && IntegerQ[l] && n >= m >= l >= 0;
	G[gG] = {};
	G[gA] = {G[x[1]], G[x[2]], G[x[3]], G[y[1]], G[y[2]], G[y[3]]};
	e[i_, j_] := Array[If[#1 == i && #2 == j, 1, 0] &, {4, 4}];
	G[_][v[0, 0, 0]] := {{0}};
	t:G[x[i_]][a:v[n_, m_, l_]] := (setRep[a]; t);
	t:G[y[i_]][a:v[n_, m_, l_]] := (setRep[a]; t);
	s:setRep[a:v[n_, m_, l_]] := s = Module[{temp},
		temp = PrintTemporary["generating irrep ", a, "..."];
		setRep[a, repMats[a]];
		NotebookDelete[temp];
	];
	s:setRep[a:v[n_, m_, l_], r_] := s = Module[{i},
		Do[
			G[x[i]][a] = r[[i]];
			G[y[i]][a] = r[[i + 3]];
		, {i, 3}];
	];
	setRep[v[1, 0, 0], rep[1]];
	setRep[v[1, 1, 0], rep[2]];
	setRep[v[1, 1, 1], rep[3]];
	G[minrep[a:v[n_, m_, l_], b:v[p_, q_, r_]]] := Which[
		n + m + l < p + q + r, a,
		n + m + l > p + q + r, b,
		n > p, a,
		n < p, b,
		m > q, a,
		True, b];
	G
]

getSO[3] := getSO[3] = AbortProtect @ Module[{G},
	G = so[3];
	G[id] = v[0];
	G[dim[v[n_]]] := 2 n + 1;
	G[prod[v[n_], v[m_]]] /; n > m := G[prod[v[n], v[m]]] = G[prod[v[m], v[n]]];
	G[prod[v[n_], v[m_]]] := G[prod[v[n], v[m]]] = Array[v, 2 n + 1, m - n];
	G[dual[v[n_]]] := v[n];
	G[isrep[_]] := False;
	G[isrep[v[n_]]] := IntegerQ[n] && n >= 0;
	G[gG] = {};
	G[gA] = {G[e], G[f]};
	G[minrep[v[n_], v[m_]]] := v[Min[n, m]];
	G[e][v[n_]] := e[n];
	G[f][v[n_]] := f[n];
	G
]

getO[2] := getO[2] = AbortProtect @ Module[{G},
	G = o[2];
	G[id] = i[1];
	G[dim[i[_]]] := 1;
	G[dim[v[_]]] := 2;
	G[prod[i[p_], i[q_]]] := {i[p q]};
	G[prod[i[_], v[n_]]] := {v[n]};
	G[prod[v[n_], i[_]]] := {v[n]};
	G[prod[v[n_], v[n_]]] := {i[1], i[-1], v[2 n]};
	G[prod[v[n_], v[m_]]] := {v[n + m], v[Abs[n - m]]};
	G[dual[v[n_]]] := v[n];
	G[dual[i[n_]]] := i[n];
	G[dual[x_List]] := G[dual[#]] & /@ x;
	G[isrep[_]] := False;
	G[isrep[i[1 | -1]]] := True;
	G[isrep[v[n_Integer]]] := n > 0;
	G[gG] = {G[t]};
	G[gA] = {G[s]};
	G[s][v[n_]] := G[s][v[n]] = {{0, -n}, {n, 0}};
	G[t][v[n_]] := G[t][v[n]] = {{1, 0}, {0, -1}};
	G[s][i[p_]] := G[s][i[p]] = {{0}};
	G[t][i[p_]] := G[t][i[p]] = {{p}};
	G[minrep[v[n_], v[m_]]] := v[Min[n, m]];
	G[minrep[i[p_], v[_]]] := i[p];
	G[minrep[v[_], i[p_]]] := i[p];
	G[minrep[i[p_], i[q_]]] := i[Max[p, q]];
	G
]

getSO[2] := getSO[2] = AbortProtect @ Module[{G},
	G = so[2];
	G[id] = v[0];
	G[dim[v[_]]] := 1;
	G[prod[v[n_], v[m_]]] := {v[n + m]};
	G[dual[v[n_]]] := v[-n];
	G[dual[x_List]] := G[dual] /@ x;
	G[isrep[_]] := False;
	G[isrep[v[x_]]] := NumericQ[x];
	G[gG] = {};
	G[gA] = {G[s]};
	G[s][v[n_]] := {{n}};
	G[minrep[v[n_], v[m_]]] := v[Min[n, m]];
	G
]

End[ ]

EndPackage[ ]
