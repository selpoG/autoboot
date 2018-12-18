Needs["NGroupInfo`", "ngroup.m"]

BeginPackage["NGroupInfoD`"]

getDihedral::usage = "getDihedral[n]"
getDicyclic::usage = "getDicyclic[n]"
dih::usage = "dih[n]"
dic::usage = "dic[n]"

Begin["`Private`"]

CommonFunctions`importPackage["NGroupInfo`", "NGroupInfoD`Private`", {"id", "dim", "prod", "dual", "isrep", "gG", "gA", "minrep", "v", "i"}]
e = CommonFunctions`e
MyReap = CommonFunctions`MyReap
s
t

getDihedral[n_] /; Mod[n, 2] == 0 := getDihedral[n] = AbortProtect @ Module[{G, tmp},
	G = dih[n];
	G[id] = i[1, 1];
	G[dim[i[_, _]]] := 1;
	G[dim[v[_]]] := 2;
	G[prod[i[a_, b_], i[c_, d_]]] := {i[a c, b d]};
	G[prod[i[_, 1], v[a_]]] := {v[a]};
	G[prod[i[_, -1], v[a_]]] := {v[n / 2 - a]};
	G[prod[v[a_], i[_, 1]]] := {v[a]};
	G[prod[v[a_], i[_, -1]]] := {v[n / 2 - a]};
	tmp[0] := (Sow[i[1, 1]]; Sow[i[-1, 1]]);
	tmp[n / 2] := (Sow[i[1, -1]]; Sow[i[-1, -1]]);
	tmp[a_] /; a < 0 := tmp[-a];
	tmp[a_] /; a >= n := tmp[Mod[a, n]];
	tmp[a_] := If[a > n / 2, Sow[v[n - a]], Sow[v[a]]];
	G[prod[v[a_], v[b_]]] := G[prod[v[a], v[b]]] = MyReap[tmp[a + b]; tmp[a - b]];
	G[dual[i[a_, b_]]] := i[a, b];
	G[dual[v[a_]]] := v[a];
	G[isrep[_]] := False;
	G[isrep[i[1|-1, 1|-1]]] := True;
	G[isrep[v[p_Integer]]] := 1 <= p < n / 2;
	G[gG] = {G[s], G[t]};
	G[gA] = {};
	G[s][v[a_]] := G[s][v[a]] = {{e[n, a], 0}, {0, e[n, -a]}};
	G[t][v[a_]] := G[t][v[a]] = {{0, 1}, {1, 0}};
	G[s][i[a_, b_]] := G[s][i[a, b]] = {{b}};
	G[t][i[a_, b_]] := G[t][i[a, b]] = {{a}};
	G[minrep[i[a_, b_], i[c_, d_]]] := If[a > c, i[a, b], If[c > a, i[c, d], i[a, Max[b, d]]]];
	G[minrep[i[a_, b_], v[_]]] := i[a, b];
	G[minrep[v[_], i[a_, b_]]] := i[a, b];
	G[minrep[v[a_], v[b_]]] := v[Min[a, b]];
	G
]

getDihedral[n_] := getDihedral[n] = AbortProtect @ Module[{G, tmp},
	G = dih[n];
	G[id] = i[1];
	G[dim[i[_]]] := 1;
	G[dim[v[_]]] := 2;
	G[prod[i[a_], i[b_]]] := {i[a b]};
	G[prod[i[_], v[a_]]] := {v[a]};
	G[prod[v[a_], i[_]]] := {v[a]};
	tmp[0] := (Sow[i[1]]; Sow[i[-1]]);
	tmp[a_] /; a < 0 := tmp[-a];
	tmp[a_] /; a >= n := tmp[Mod[a, n]];
	tmp[a_] := If[a > n / 2, Sow[v[n - a]], Sow[v[a]]];
	G[prod[v[a_], v[b_]]] := G[prod[v[a], v[b]]] = MyReap[tmp[a + b]; tmp[a - b]];
	G[dual[i[a_]]] := i[a];
	G[dual[v[a_]]] := v[a];
	G[isrep[_]] := False;
	G[isrep[i[1|-1]]] := True;
	G[isrep[v[p_Integer]]] := 1 <= p < n / 2;
	G[gG] = {G[s], G[t]};
	G[gA] = {};
	G[s][v[a_]] := G[s][v[a]] = {{e[n, a], 0}, {0, e[n, -a]}};
	G[t][v[a_]] := G[t][v[a]] = {{0, 1}, {1, 0}};
	G[s][i[a_]] := G[s][i[a, b]] = {{1}};
	G[t][i[a_]] := G[t][i[a, b]] = {{a}};
	G[minrep[i[a_], i[b_]]] := i[Max[a, b]];
	G[minrep[i[a_], v[_]]] := i[a];
	G[minrep[v[_], i[a_]]] := i[a];
	G[minrep[v[a_], v[b_]]] := v[Min[a, b]];
	G
]

getDicyclic[n_] := getDicyclic[n] = AbortProtect @ Module[{G, tmp},
	G = dic[n];
	G[id] = i[1, 1];
	G[dim[i[_, _]]] := 1;
	G[dim[v[_]]] := 2;
	If[Mod[n, 2] == 0,
		G[prod[i[a_, b_], i[c_, d_]]] := {i[a c, b d]},
		G[prod[i[1, a_], i[1, b_]]] := {i[1, a b]};
		G[prod[i[-1, a_], i[-1, b_]]] := {i[1, -a b]};
		G[prod[i[-1, a_], i[1, b_]]] := {i[-1, a b]};
		G[prod[i[1, a_], i[-1, b_]]] := {i[-1, a b]}
	];
	G[prod[i[1, _], v[a_]]] := {v[a]};
	G[prod[i[-1, _], v[a_]]] := {v[n - a]};
	G[prod[v[a_], i[1, _]]] := {v[a]};
	G[prod[v[a_], i[-1, _]]] := {v[n - a]};
	tmp[0] := (Sow[i[1, 1]]; Sow[i[1, -1]]);
	tmp[n] := (Sow[i[-1, 1]]; Sow[i[-1, -1]]);
	tmp[a_] /; a < 0 := tmp[-a];
	tmp[a_] /; a >= 2 * n := tmp[Mod[a, 2 * n]];
	tmp[a_] := If[a > n, Sow[v[2 * n - a]], Sow[v[a]]];
	G[prod[v[a_], v[b_]]] := G[prod[v[a], v[b]]] = MyReap[tmp[a + b]; tmp[a - b]];
	If[Mod[n, 2] == 0,
		G[dual[i[a_, b_]]] := i[a, b],
		G[dual[i[1, b_]]] := i[1, b];
		G[dual[i[-1, b_]]] := i[-1, -b]
	];
	G[dual[v[a_]]] := v[a];
	G[isrep[_]] := False;
	G[isrep[i[1|-1, 1|-1]]] := True;
	G[isrep[v[p_Integer]]] := 1 <= p < n;
	G[gG] = {G[s], G[t]};
	G[gA] = {};
	G[s][v[a_]] := G[s][v[a]] = {{e[2 n, a], 0}, {0, e[2 n, -a]}};
	G[t][v[a_]] := G[t][v[a]] = {{0, (-1)^a}, {1, 0}};
	G[s][i[a_, b_]] := G[s][i[a, b]] = {{a}};
	If[Mod[n, 2] == 0,
		G[t][i[a_, b_]] := G[t][i[a, b]] = {{b}},
		G[t][i[1, b_]] := G[t][i[1, b]] = {{b}};
		G[t][i[-1, b_]] := G[t][i[-1, b]] = {{b I}}
	];
	G[minrep[i[a_, b_], i[c_, d_]]] := If[a > c, i[a, b], If[c > a, i[c, d], i[a, Max[b, d]]]];
	G[minrep[i[a_, b_], v[_]]] := i[a, b];
	G[minrep[v[_], i[a_, b_]]] := i[a, b];
	G[minrep[v[a_], v[b_]]] := v[Min[a, b]];
	G
]

End[ ]

EndPackage[ ]
