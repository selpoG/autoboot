(* Loaded by separate exact and numeric kernels. *)
Begin["SUChecks`"];
failed = False;
assert[name_, condition_] := If[TrueQ[condition], Print[name, "=PASS"],
	Print[name, "=FAIL"]; failed = True];
nearZero[m_] := Max[Abs[Flatten[N[Simplify[Normal[m]], 30]]]] < 10^-20;
comm[a_, b_] := a.b - b.a;

assert["invalid degrees", And @@ (Quiet[get[#]] === $Failed & /@ {1, 0, -1, 3/2, 3., q})];
g2 = get[2];
assert["SU2 spin labels", g2[dim[v[1/2]]] == 2 && g2[dual[v[1/2]]] === v[1/2]];

Do[
	g = get[n]; rank = n - 1;
	fund = v @@ PadRight[{1}, rank];
	anti = v @@ ConstantArray[1, rank];
	adj = v @@ Prepend[ConstantArray[1, rank - 1], 2];
	assert["SU" <> ToString[n] <> " labels", g[id] === v @@ ConstantArray[0, rank] &&
		g[dim[fund]] == n && g[dim[adj]] == n^2 - 1 && g[dual[fund]] === anti &&
		g[dual[anti]] === fund];
	assert["invalid weights", !Or @@ (g[isrep[#]] & /@
		{v[], v[1], v @@ ConstantArray[0, n], v @@ PadRight[{1/2}, rank],
		v @@ PadRight[{-1}, rank], v @@ PadRight[{0, 1}, rank], v @@ PadRight[{q}, rank]})];
	assert["tensor decomposition", Sort[g[prod[fund, anti]]] === Sort[{g[id], adj}] &&
		g[prod[fund, anti]] === g[prod[anti, fund]]];
	assert["ordering", g[minrep[fund, anti]] === g[minrep[anti, fund]] &&
		g[minrep[fund, fund]] === fund];
	(* Interleave different group objects before asking for uncached matrices. *)
	get[n + 1];
	Do[
		m = (#[r] & /@ g[gA]); d = g[dim[r]];
		assert["matrix dimensions", Dimensions[m] === {2 rank, d, d}];
		lower = Take[m, rank]; upper = Drop[m, rank];
		h = MapThread[comm, {upper, lower}];
		assert["adjoints", And @@ Table[nearZero[upper[[i]] - ConjugateTranspose[lower[[i]]]], {i, rank}]];
		assert["Chevalley relations", And @@ Flatten[Table[
			{nearZero[comm[upper[[i]], lower[[j]]] - If[i == j, h[[i]], 0 lower[[j]]]],
			nearZero[comm[h[[i]], lower[[j]]] + (2 Boole[i == j] - Boole[Abs[i-j] == 1]) lower[[j]]],
			nearZero[comm[h[[i]], h[[j]]]]}, {i, rank}, {j, rank}]]];
		assert["Serre relations", And @@ Flatten[Table[If[i == j, True,
			If[Abs[i-j] == 1, nearZero[comm[lower[[i]], comm[lower[[i]], lower[[j]]]]],
				nearZero[comm[lower[[i]], lower[[j]]]]]], {i, rank}, {j, rank}]]];
	, {r, {g[id], fund, anti, v @@ PadRight[{2}, rank], adj}}];
, {n, {3, 4, 5}}];

g = get[3];
assert["SU3 repeated multiplicity", Count[g[prod[v[2,1], v[2,1]]], v[2,1]] == 2 &&
	Total[g[dim[#]] & /@ g[prod[v[2,1], v[2,1]]]] == 64];
setGroup[g];
Do[
	{a,b,c} = triple;
	equations = eq[a, b, c];
	assert["sparse equations", Head[equations] === SparseArray];
	assert["nullity", Length[NullSpace[equations]] == 1];
	cg = Flatten[Table[ope[a,b,c][1][i,j,k], {i,g[dim[a]]}, {j,g[dim[b]]}, {k,g[dim[c]]}]];
	assert["CG intertwiner", VectorQ[cg, NumericQ] && nearZero[equations.cg] && nearZero[{Conjugate[cg].cg - g[dim[c]]}]];
, {triple, {{v[1,0],v[1,1],v[0,0]}, {v[1,0],v[1,0],v[1,1]}, {v[1,0],v[1,0],v[2,0]}}}];
Print["SU_REGRESSION=", If[failed, "FAIL", "PASS"]];
Exit[If[failed, 1, 0]];
End[];
