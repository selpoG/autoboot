(* Shared suite, run in isolated exact and numerical kernels. *)
Begin["LieChecks`"];
failed = False;
assert[name_, condition_] := If[TrueQ[condition], Print[name,"=PASS"], Print[name,"=FAIL"]; failed = True];
nearZero[m_] := Max[Append[Abs[N[RootReduce[Last /@ Most[ArrayRules[SparseArray[m]]]],30]],0]] < 10^-20;
comm[a_, b_] := a.b-b.a;

assert["invalid Cartan types", And @@ (Quiet[getLie @@ #] === $Failed & /@
	{{"Z",2},{"E",5},{"F",3},{"G",3},{"B",1},{"D",1},{"A",0},{"C",1.5},{"A",x}})];
assert["invalid classical degrees", Quiet[getSp[0]] === $Failed && Quiet[getSpin[3]] === $Failed && Quiet[getSO[1]] === $Failed];

samples = {
	{"A",1,{2},3}, {"B",2,{1,0},5}, {"B",2,{0,1},4}, {"B",3,{0,0,1},8},
	{"C",1,{1},2}, {"C",2,{1,0},4}, {"C",3,{0,1,0},14},
	{"D",2,{1,1},4}, {"D",3,{0,0,1},4}, {"D",4,{0,0,0,1},8},
	{"G",2,{1,0},7}, {"G",2,{0,1},14}, {"G",2,{1,1},64},
	{"F",4,{0,0,0,1},26}, {"E",6,{1,0,0,0,0,0},27},
	{"E",7,{0,0,0,0,0,0,1},56}, {"E",8,{0,0,0,0,0,0,0,1},248}};
Do[
	{t,rank,labels,d} = sample; g = getLie[t,rank]; r = v @@ labels;
	assert[ToString[{t,rank,labels}] <> " dimension", g[dim[r]] === d];
	assert["labels and duality", g[isrep[r]] && g[dual[g[dual[r]]]] === r &&
		!g[isrep[v[]]] && !g[isrep[v @@ ConstantArray[-1,rank]]] && !g[isrep[v @@ ConstantArray[1/2,rank]]]];
	(* Interleave group initialization to detect leaking generator definitions. *)
	getLie["A",rank+1];
	m = (#[r] & /@ g[gA]);
	assert["matrix dimensions", Dimensions[m] === {2 rank,d,d}];
	lower = Take[m,rank]; upper = Drop[m,rank]; h = MapThread[comm,{upper,lower}];
	c = LieRepresentations`cartan[t,rank];
	assert["adjoints", And @@ Table[nearZero[upper[[i]]-ConjugateTranspose[lower[[i]]]],{i,rank}]];
	assert["Chevalley relations", And @@ Flatten[Table[
		{nearZero[comm[upper[[i]],lower[[j]]] - If[i==j,h[[i]],0 lower[[j]]]],
		 nearZero[comm[h[[i]],lower[[j]]] + c[[i,j]] lower[[j]]], nearZero[comm[h[[i]],h[[j]]]]},
		{i,rank},{j,rank}]]];
	assert["Serre relations", And @@ Flatten[Table[If[i==j,True,
		nearZero[Nest[comm[lower[[i]],#] &,lower[[j]],1-c[[i,j]]]]],{i,rank},{j,rank}]]];
	computedCharacter = LieRepresentations`repCharacter[t,rank,labels];
	observed = Association[Rule @@@ Tally[Round[Transpose[Diagonal[Normal[#]] & /@ h]]]];
	assert["character agrees with matrices", Sort[Normal[computedCharacter]] === Sort[Normal[observed]] && Total[Values[computedCharacter]] === d];
, {sample,samples}];

assert["SO excludes spinors", !getSO[5][isrep[v[0,1]]] && getSpin[5][isrep[v[0,1]]] &&
	!getSO[6][isrep[v[0,0,1]]] && getSpin[6][isrep[v[0,0,1]]] && getSO[6][isrep[v[0,1,1]]]];
assert["classical aliases", getSO[4][dim[v[1,1]]] === 4 && getSp[1][dim[v[1]]] === 2 &&
	getSp[3][dim[v[1,0,0]]] === 6 && getSO[7][dim[v[1,0,0]]] === 7];
assert["complex duals", getSpin[6][dual[v[0,0,1]]] === v[0,1,0] &&
	getLie["E",6][dual[v[1,0,0,0,0,0]]] === v[0,0,0,0,0,1]];

Do[
	{g,a,expected} = sample;
	decomp = g[prod[a,a]];
	assert["known tensor product", Sort[g[dim[#]] & /@ decomp] === Sort[expected] && And @@ (g[isrep[#]] & /@ decomp)];
	assert["tensor unit", g[prod[a,g[id]]] === {a} && g[prod[g[id],a]] === {a}];
, {sample, {{getSO[5],v[1,0],{1,10,14}}, {getSp[2],v[1,0],{1,5,10}},
	{getLie["G",2],v[1,0],{1,7,14,27}}, {getSpin[4],v[1,0],{1,3}}}}];
g = getLie["A",2];
assert["tensor multiplicity", Count[g[prod[v[1,1],v[1,1]]],v[1,1]] == 2];

Do[
	{g,a,pseudo} = sample;
	setGroup[g]; b = g[dual[a]]; one = g[id];
	equations = eq[a,b,one];
	cg = Flatten[Table[ope[a,b,one][1][i,j,1],{i,g[dim[a]]},{j,g[dim[b]]}]];
	assert["CG singlet", VectorQ[cg,NumericQ] && nearZero[equations.cg] && nearZero[{Conjugate[cg].cg-1}]];
	assert["reality type", isPseaudo[a] === pseudo];
, {sample, {{getSO[5],v[1,0],False},{getSp[2],v[1,0],True},{getLie["G",2],v[1,0],False}}}];

Print["LIE_REGRESSION=",If[failed,"FAIL","PASS"]];
Exit[If[failed,1,0]];
End[];
