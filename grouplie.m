Needs["GroupInfo`", "group.m"]

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
s
t
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
