Needs["NGroupInfo`", "ngroup.m"]
Needs["RootSystem`", "root.m"]
Needs["SURepresentations`", "su-representations.m"]
Needs["LieRepresentations`", "lie-representations.m"]
Needs["OrthogonalRepresentations`", "orthogonal-representations.m"]

BeginPackage["NGroupInfoLie`"]

getSU::usage = "getSU[n] returns group-object su[n] which represents the special unitary group of degree n (rank n-1). n must be an integer >= 2."
getSU::degree = "SU degree `1` must be an integer >= 2."
getO::usage = "getO[n] returns O(n), n >= 2. O(2,3) retain their labels; n >= 4 uses v[a1,...,ar,p] with Dynkin labels and extension parity p=+/-1, or p=0 for an induced pair."
getSO::usage = "getSO[n] returns the special orthogonal group so[n], for integer n >= 2. SO(2,3) retain their existing labels; n >= 4 uses Dynkin labels that descend from Spin(n)."
su::usage = "su[n] is a group-object which is the special unitary group of degree n (rank n-1). Before using this value, you have to call getSU[n] to get proper group-object."
o::usage = "o[n] is a group-object which is the orthogonal group of degree n. Before using this value, you have to call getO[n] to get proper group-object."
so::usage = "so[n] is a group-object which is the special orthogonal group of degree n. Before using this value, you have to call getSO[n] to get proper group-object."

getLie::usage = "getLie[type,rank] constructs a compact simply connected Lie group of Cartan type A, B, C, D, E, F or G. Irreps v[a1,...,ar] use Dynkin labels. D2 denotes Spin(4)."
getLie::type = "Unsupported finite Cartan type and rank: `1`."
lie::usage = "lie[type,rank] is initialized by getLie[type,rank]. Its irreps use Dynkin labels."
getSp::usage = "getSp[n] constructs compact Sp(n), of rank n with defining dimension 2n. Irreps use Dynkin labels. n must be an integer >= 1."
getSp::degree = "Sp requires one integer rank >= 1; got `1`."
sp::usage = "sp[n] is the compact symplectic group initialized by getSp[n]."
getSpin::usage = "getSpin[n] constructs Spin(n), including spinor irreps with Dynkin labels. n must be an integer >= 4."
getSpin::degree = "Spin requires one integer degree >= 4; got `1`. Use getSU[2] for Spin(3)."
spin::usage = "spin[n] is initialized by getSpin[n]."
getO::degree = "O requires one integer degree >= 2; got `1`."
getSO::degree = "SO requires one integer degree >= 2; got `1`."

(* all irrep-objects of G=su[2] are v[0], v[1/2], v[1], v[3/2], .... *)
(* all irrep-objects of G=o[3] are v[0,1], v[0,-1], v[1,1], v[1,-1], v[2,1], v[2,-1], v[3,1], v[3,-1], .... *)
(* all irrep-objects of G=so[3] are v[0], v[1], v[2], v[3], .... *)
(* all irrep-objects of G=o[2] are i[1], i[-1], v[1], v[2], v[3], .... *)
(* all irrep-objects of G=so[2] are v[x] (x \in \mathbb{R}). *)

Begin["`Private`"]

CommonFunctions`importPackage["NGroupInfo`", "NGroupInfoLie`Private`", {"id", "dim", "prod", "dual", "isrep", "gG", "gA", "minrep", "v", "i"}]
CommonFunctions`importPackage["RootSystem`", "NGroupInfoLie`Private`", {"dimension", "irrep", "productReps", "decompose"}]
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

getSU[n_Integer] /; n >= 3 := NGroupInfo`Private`checkPrec[getSU[n] = AbortProtect @ Module[{G = su[n], rank = n - 1, gen},
	G[id] = v @@ ConstantArray[0, rank];
	G[isrep[_]] := False;
	G[isrep[r_v]] := SURepresentations`validWeight[n, List @@ r];
	G[dim[r_v]] /; G[isrep[r]] := dimension["A", rank, List @@ r];
	G[dual[r_v]] /; G[isrep[r]] := v @@ Prepend[r[[1]] - Reverse[Rest[List @@ r]], r[[1]]];
	G[minrep[a_v, b_v]] /; G[isrep[a]] && G[isrep[b]] :=
		First @ SortBy[{a, b}, {Total[List @@ #] &, (-(List @@ #)) &}];
	G[prod[a_v, b_v]] /; G[isrep[a]] && G[isrep[b]] :=
		G[prod[a, b]] = G[prod[b, a]] = (v @@ # & /@ decompose @
			productReps[irrep["A", rank, List @@ a], irrep["A", rank, List @@ b]]);
	G[gG] = {};
	G[gA] = Array[G[gen[#]] &, 2 rank];
	G[gen[j_Integer]][r_v] /; 1 <= j <= 2 rank && G[isrep[r]] := Module[{m},
		m = SURepresentations`matrices[n, List @@ r];
		If[m === $Failed, $Failed, G[gen[j]][r] = NGroupInfo`Private`num[m[[j]]]]
	];
	G
]]

getSU[n_] := (Message[getSU::degree, n]; $Failed)

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

lieCheck = NGroupInfo`Private`checkPrec
lieNumber = NGroupInfo`Private`num
Get["lie-groups.m"]

End[ ]

EndPackage[ ]
