(* Do not import this package directly. *)
Needs["CommonFunctions`", "common.m"]
Needs["ToPython`", "topy.m"]
ClebschGordan`clearCG[];

BeginPackage["ClebschGordan`"]

symmetryGroup::usage = "symmetryGroup represents the global symmetry of CFT."
clearCG::usage = "clearCG[] clears all values calculated by this package."

(* r,s,t,r1,r2 must be irrep-objects, i.e. such as symmetryGroup[isrep[r]] must be True *)
(* id means the trivial representation of symmetryGroup. *)
inv::usage = "inv[r,s,t] gives multiplicity of id in decomposition of r\[TensorProduct]s\[TensorProduct]t.
inv[{r,s},{t}] gives multiplicity of t in decomposition of r\[TensorProduct]s.
inv[r1,r2,r3,r4] gives an association whose key is a irrep s s.t. inv[{r1,r2},{s}]>0 and inv[r3,r4,s]>0, whose value is a list {inv[{r1,r2},{s}],inv[r3,r4,s]}."
invs::usage = "invs[r1,r2,r3,r4] gives a list of {s,n,m} s.t. cor[r1,r2,r3,r4][s,n,m] is defined."

(* a,b,c,a1,a2,... are indices of irreps. An index a of irrep r must sutisfy 1<=a<=symmetryGroup[dim[r]]. *)
(* n,m represents multiplicity of Clebsch-Gordan coefficients. For more information, please refer to autoboot.pdf. *)
ope::usage = "ope[r,s,t][n][a,b,c] gives Clebsch-Gordan coefficient defined by |r,a\[RightAngleBracket]\[TensorProduct]|s,b\[RightAngleBracket]=\[Sum]ope[r,s,t][n][a,b,c]|t,c,n\[RightAngleBracket].
ope[r][a,b] gives Sqrt[dim[r]]ope[r,dual[r],id][1][a,b,1]."
cor::usage = "cor[r,s,t][n][a,b,c] gives \[Sum]ope[r,s,dual[t]][n][a,b,c2] ope[dual[t]][c2,c]/Sqrt[dim[t]].
cor[r1,r2,r3,r4][s,n,m][a1,a2,a3,a4] gives \[Sum]cor[r1,r2,dual[s]][n][a1,a2,b] ope[r3,r4,dual[s]][a3,a4,b]."

\[Sigma]::usage = "\[Sigma][r,s,t][n] gives the sign change in swap between r and s in cor, i.e. cor[r,s,t][n][a,b,c]=\[Sigma][r,s,t][n]cor[s,r,t][n][b,a,c].
\[Sigma][r] gives \[Sigma][r,dual[r],id].
\[Sigma][op[x,r,p,q]] gives p."
\[Tau]::usage = "\[Tau][r,s,t][n,m] describes the behavior of cor under cyclic permutation of r,s,t, i.e. cor[s,t,r][m][b,c,a]=\[Sum]cor[r,s,t][n][a,b,c]\[Tau][r,s,t][n,m]."
\[Omega]::usage = "\[Omega][r,s,t][n,m] describes the behavior of cor under complex conjugate, i.e. \[Sum]ope[dual[r]][a2,a] ope[dual[s]][b2,b] ope[dual[t]][c2,c] cor[dual[r],dual[s],dual[t]][m][a2,b2,c2]\[Conjugate]=\[Sum]cor[r,s,t][n][a,b,c]\[Omega][r,s,t][n,m]."
six::usage = "six[r1,r2,r3,r4][s,n,m,t,k,l] describes the behavior of cor under swap between r2 and r4, i.e. cor[r1,r4,r3,r2][t,k,l][a1,a4,a3,a2]=\[Sum]cor[r1,r2,r3,r4][s,n,m][a1,a2,a3,a4] six[r1,r2,r3,r4][s,n,m,t,k,l]."

isReal::usage = "isReal[r] gives whether r is a real representation or not."
isComplex::usage = "isComplex[r] gives whether r is a complex representation or not."
isPseaudo::usage = "isPseaudo[r] gives whether r is a pseaudo real representation or not."

op::usage = "op[x,r,p,q] represents a primary operator object O whose name is x, which belongs to irrep r, which sign (=\[Sigma][O]) is p and whose parity of spin is q (=±1). If O appears in summation, we use 'op' as its name."
dualOp::usage = "dualOp[op[x,r,p,q]] gives dual operator object of op[x,r,p,q]."

format::usage = "format[eq] gives readable format of eq with little loss of information. We need much redundancy to calculate properly, so formatted value cannot be used for any argument of our function. format[eq] is assumed to be used only for human-readability of last output."

sum::usage = "sum[x,op[o,r,1,q]] represents sum of x over all intermediate primary operator O=op[op,r,1,q] which belongs to irrep r and whose parity of spin is q. sum[...] is automatically expanded so as x to be (F or H)*\[Beta]^2."
single::usage = "single[x] represents x. We need to wrap x for redundancy. single[...] is automatically expanded so as x to be (Fp or Hp)*\[Beta]^2, or (Fp or Hp)."

F::usage = "F[o1,o2,o3,o4,s,p] represents generalized conformal block "<>ToString[Subsuperscript["F", {"p", "s"}, {1, 2, 3, 4}], StandardForm]<>" , where o1,...,o4 are primary scalars.
F[a,b,c,d] represents normal conformal block of type-F."
H::usage = "H[a,b,c,d] represents normal conformal block of type-H."
Fp::usage = "Fp[o1,o2,o3,o4,o] represents generalized conformal block "<>ToString[Subsuperscript["F", {"o"}, {1, 2, 3, 4}], StandardForm]<>" , where o1,...,o4 and intermediate o are primary scalars.
Fp[a,b,c,d,o] represents normal conformal block of type-F with intermediate o."
Hp::usage = "Hp[a,b,c,d,o] represents normal conformal block of type-H with intermediate o."

\[Lambda]::usage = "\[Lambda][o1,o2,o3][n] gives n-th OPE coefficient of o1\[Times]o2\[Rule]"<>ToString[OverBar["o3"], StandardForm]<>"."
\[Alpha]::usage = "\[Alpha][o1,o2,o3][n] gives n-th coefficient of three-point function \[LeftAngleBracket]o1 o2 o3\[RightAngleBracket]."
\[Mu]::usage = "\[Mu][o1,o2,o3][n] gives n-th OPE coefficient of o1\[Times]o2\[Rule]"<>ToString[OverBar["o3"], StandardForm]<>", where o3 is registered as a fundamental scalar."
\[Nu]::usage = "\[Nu][o1,o2,o3][n] gives n-th coefficient of three-point function \[LeftAngleBracket]o1 o2 o3\[RightAngleBracket], where o3 is registered as a fundamental scalar."
(* \[Lambda], \[Alpha], \[Mu] and \[Nu] are linear combination of \[Beta]. *)
\[Beta]::usage = "\[Beta][o1,o2,o3][n] gives minimal basis to describe \[Lambda], \[Alpha], etc.."

eqn::usage = "eqn[{a,b,...}] represents bootstrap equation that claims a,b,... must equals to 0, where a,b,... are real-linear combination of sum and single.
eqn[sec,{a,b,...}] represents extracted part of eqn[{...}] which contains only sum or single related to sec."
bootAll::usage = "bootAll[ops] generates bootstrap equation from all four-point function of ops.
bootAll[] generates bootstrap equation from all four-point function of fundamental scalars."

(* fundamental scalars are used to seperate sum of conformal blocks over scalars to single[...] and sum[...]. *)
setOps::usage = "setOps[ops] registers ops and duals of ops as fundamental scalars."
one::usage = "one represents unit operator. This is implicitly registered as a fundamental scalar."

extract::usage = "extract[x,op[op,r,1,p]]
extract[x,scalar]
extract[x,unit]"
unit::usage = "unit"
scalar::usage = "scalar"
sector::usage = "sector[eq]"

makeG::usage = "makeG[eqn[sec,{a,b,...}]]"
makeMat::usage = "makeMat[eqn[sec,{a,b,...}]]"
makeSDP::usage = "makeSDP[eqn[{a,b,...}]]"
sdp::usage = "sdp[secs,scalarnum,vals,mat]"
toString::usage = "toString[sdp]"

Begin["`Private`"]

setOPE::usage = "setOPE[r,s,t]
setOPE[r,s]"
setAllOPE::usage = "setAllOPE[reps]"
boot::usage = "boot[o1,o2,o3,o4,s,p]"
connectedComponents::usage = "connectedComponents[z]"

CommonFunctions`importPackage["CommonFunctions`", "ClebschGordan`Private`", {"importPackage", "lhs", "rhs", "MyReap", "newUF", "newSet", "classify", "keys", "add", "reverseIndex", "zero", "root", "unite", "has"}]
setGroup[G_] := AbortProtect @ Module[{x, name = "ClebschGordan`Private`"},
	symmetryGroup = G;
	Scan[(Evaluate[ToExpression[name <> #]][x___] := G[ToExpression[#][x]]) &,
		{"isrep", "dim", "dual", "minrep", "inv"}];
	Evaluate[ToExpression[name <> "gprod"]][x___] := G[ToExpression["prod"][x]];
	Evaluate[ToExpression[name <> "prod"]][r_, s_] := ToExpression[name <> "prod"][r, s] =
		Module[{y = <||>}, Scan[If[KeyExistsQ[y, #], y[#]++, y[#] = 1]&, gprod[r, s]]; y];
	Scan[(Evaluate[ToExpression[name <> #]] = G[ToExpression[#]]) &, {"id", "gG", "gA"}];
	one = op[0, id, 1, 1];
]
clearCG[] := ClearAll["ClebschGordan`*", "ClebschGordan`Private`*"]

x:inv[{r_, s_}, {t_}] := x = If[KeyExistsQ[prod[r, s], t], prod[r, s][t], 0]
x:inv[r_, s_, t_] := x = inv[{r, s}, {dual[t]}]
x:inv[r1_, r2_, r3_, r4_] := x = Module[{s}, Association @ Table[
	If[inv[s, r3, r4] > 0, s -> {prod[r1, r2][s], inv[s, r3, r4]}, Nothing], {s, Keys @ prod[r1, r2]}]]
x:invs[r1_, r2_, r3_, r4_] := x = Module[{t = inv[r1, r2, r3, r4], s, n, m},
	MyReap @ Do[Sow[{s, n, m}], {s, Keys[t]}, {n, t[s][[1]]}, {m, t[s][[2]]}]]

eq[r_, s_, t_] /; inv[{r, s}, {t}] > 0 :=
Module[{a, b, c, x, y, z, X, g, e, d1 = dim[r], d2 = dim[s], d3 = dim[t], d},
	d = d1 d2 d3;
	MyReap[
		Do[
			e = ConstantArray[0, d];
			Do[e[[d2 d3 (x - 1) + d3 (b - 1) + c]] -= X[r][[a, x]], {x, d1}];
			Do[e[[d2 d3 (a - 1) + d3 (y - 1) + c]] -= X[s][[b, y]], {y, d2}];
			Do[e[[d2 d3 (a - 1) + d3 (b - 1) + z]] += X[t][[z, c]], {z, d3}];
			Sow[e],
		{a, d1}, {b, d2}, {c, d3}, {X, gA}];
		Do[
			e = ConstantArray[0, d];
			Do[e[[d2 d3 (x - 1) + d3 (y - 1) + c]] -= g[r][[a, x]] g[s][[b, y]], {x, d1}, {y, d2}];
			Do[e[[d2 d3 (a - 1) + d3 (b - 1) + z]] += g[t][[z, c]], {z, d3}];
			Sow[e],
		{a, d1}, {b, d2}, {c, d3}, {g, gG}]
	]
]

x:setOPE[r_, r_, t_] /; inv[{r, r}, {t}] > 0 := x =
Module[{sol, a, b, c, d = dim[r], d3 = dim[t], l, n, sym, tmp, v, p, even, odd},
	sym[v_, p_] := (
		tmp = v;
		Do[tmp[[d d3 (a - 1) + d3 (b - 1) + c]] += p v[[d d3 (b - 1) + d3 (a - 1) + c]], {a, d}, {b, d}, {c, d3}];
		Expand[tmp]);
	l = Length[sol = NullSpace @ Expand @ eq[r, r, t]];
	If[l == 0, Message[setOPE::imcmpt, r, r, t]; Return[]];
	If[l != inv[{r, r}, {t}], Message[setOPE::diff, r, r, t, l, inv[{r, r}, {t}]]; Return[]];
	even = Select[simp @ Orthogonalize[simp[sym[#, +1] & /@ sol], simp[Conjugate[#1].#2] &], AnyTrue[#, simp @ # != 0 &] &];
	odd  = Select[simp @ Orthogonalize[simp[sym[#, -1] & /@ sol], simp[Conjugate[#1].#2] &], AnyTrue[#, simp @ # != 0 &] &];
	l = Length[even];
	AbortProtect[
		Do[ope[r, r, t][n][a, b, c] = even[[n, d d3 (a - 1) + d3 (b - 1) + c]] Sqrt[d3]
		, {n, l}, {a, d}, {b, d}, {c, d3}];
		Do[ope[r, r, t][n + l][a, b, c] = odd[[n, d d3 (a - 1) + d3 (b - 1) + c]] Sqrt[d3]
		, {n, Length[odd]}, {a, d}, {b, d}, {c, d3}];
		\[Sigma][r, r, dual[t]][n_] := If[n <= l, 1, -1]]
]

x:setOPE[r_, s_, t_] /; inv[{r, s}, {t}] > 0 := x = setOPE[s, r, t] =
Module[{sol, a, b, c, d1 = dim[r], d2 = dim[s], d3 = dim[t], l, n},
	l = Length[sol = NullSpace @ Expand @ eq[r, s, t]];
	If[l == 0, Message[setOPE::imcmpt, r, s, t]; Return[]];
	If[l != inv[{r, s}, {t}], Message[setOPE::diff, r, s, t, l, inv[{r, s}, {t}]]; Return[]];
	sol = simp @ Orthogonalize[sol, simp[Conjugate[#1].#2] &];
	AbortProtect[
		Do[ope[r, s, t][n][a, b, c] = ope[s, r, t][n][b, a, c] = sol[[n, d2 d3 (a - 1) + d3 (b - 1) + c]] Sqrt[d3]
		, {n, l}, {a, d1}, {b, d2}, {c, d3}];
		\[Sigma][r, s, dual[t]][_] := 1;
		\[Sigma][s, r, dual[t]][_] := 1]
]
setOPE[r_, s_] := Scan[setOPE[r, s, #] &, Keys[prod[r, s]]]
setAllOPE[reps_List] := Module[{r, s}, Do[setOPE[r, s]; Scan[setOPE[#] &, Keys @ prod[r, s]], {r, reps}, {s, reps}]]
setOPE::imcmpt = "setOPE[`1`,`2`,`3`]: There is no compatible CG coefficient.";
setOPE::diff = "setOPE[`1`,`2`,`3`]: The dimension of CG coefficients is `4`, which must be `5`.";
ope[r_, i_, r_][1][a_, 1, b_] /; i === id := If[a == b, 1, 0]
ope[i_, r_, r_][1][1, a_, b_] /; i === id := If[a == b, 1, 0]
x:ope[r_, s_, t_][n_][a_, b_, c_] /; 1 <= n <= inv[{r, s}, {t}] := (setOPE[r, s, t]; x)
ope[_, _, _][_][_, _, _] := 0

x:setOPE[r_] := x = Module[{r2 = dual[r], a, b, d = dim[r], sd},
	sd = Sqrt[d]; AbortProtect @ Do[ope[r][a, b] = simp[sd ope[r, r2, id][1][a, b, 1]], {a, d}, {b, d}]]
x:ope[r_][a_, b_] := (setOPE[r]; x)

x:setCor[r_, s_, t_] := x = Module[{n, a, b, c, c2, d = dim[t], sd, t2 = dual[t], res},
	sd = Sqrt[d];
	Do[res[n, a, b, c] = simp[Sum[ope[r, s, t2][n][a, b, c2] ope[t2][c2, c], {c2, d}] / sd]
	, {n, inv[r, s, t]}, {a, dim[r]}, {b, dim[s]}, {c, d}];
	AbortProtect @ Do[cor[r, s, t][n][a, b, c] = res[n, a, b, c]
	, {n, inv[r, s, t]}, {a, dim[r]}, {b, dim[s]}, {c, d}];
]
x:cor[r_, s_, t_][n_][a_, b_, c_] /; 1 <= n <= inv[r, s, t] := (setCor[r, s, t]; x)
cor[_, _, _][_][_, _, _] := 0

x:setCor[r1_, r2_, r3_, r4_, s_] := x = Module[{n, m, a1, a2, a3, a4, s2 = dual[s], res},
	Do[res[n, m, a1, a2, a3, a4] = simp @ Sum[cor[r1, r2, s2][n][a1, a2, b] ope[r3, r4, s2][m][a3, a4, b], {b, dim[s]}]
	, {n, inv[{r1, r2}, {s}]}, {m, inv[r3, r4, s]}, {a1, dim[r1]}, {a2, dim[r2]}, {a3, dim[r3]}, {a4, dim[r4]}];
	AbortProtect @ Do[cor[r1, r2, r3, r4][s, n, m][a1, a2, a3, a4] = res[n, m, a1, a2, a3, a4],
		{n, inv[{r1, r2}, {s}]}, {m, inv[r3, r4, s]}, {a1, dim[r1]}, {a2, dim[r2]}, {a3, dim[r3]}, {a4, dim[r4]}];
]
x:cor[r1_, r2_, r3_, r4_][s_, n_, m_][a1_, a2_, a3_, a4_] /;
	1 <= n <= inv[{r1, r2}, {s}] && 1 <= m <= inv[r3, r4, s] := (setCor[r1, r2, r3, r4, s]; x)
cor[_, _, _, _][_, _, _][_, _, _, _] := 0

dec[r_, s_, t_][f_] := Module[{vec},
	decPrep[r, s, t]; vec = Array[{f @@ bas[r, s, t][#]} &, inv[r, s, t]]; simp @ Flatten[invmat[r, s, t].vec]]

x:decPrep[r_, s_, t_] /; inv[r, s, t] > 0 := x =
Module[{n, f = cor[r, s, t], max = inv[r, s, t], mat, rev, new, old, sc, a, b, c, tmp, tmp2, res},
	Do[If[f[1][a, b, c] != 0, res[1] = {a, b, c}; mat = {{f[1][a, b, c]}}; rev = Inverse[mat]; Break[]]
	, {a, dim[r]}, {b, dim[s]}, {c, dim[t]}];
	Do[Do[
		new = Array[{f[n] @@ res[#]} &, n - 1];
		old = {Array[f[#][a, b, c] &, n - 1]};
		sc = f[n][a, b, c];
		tmp2 = rev.new;
		If[(tmp = simp[sc - (old.tmp2)[[1, 1]]]) != 0,
			res[n] = {a, b, c};
			mat = ArrayFlatten[{{mat, new}, {old, {{sc}}}}];
			rev = ArrayFlatten[{{rev + tmp2.old.rev / tmp, -tmp2 / tmp}, {-old.rev / tmp, {{1 / tmp}}}}];
			Break[]]
	, {a, dim[r]}, {b, dim[s]}, {c, dim[t]}], {n, 2, max}];
	Assert[simp[rev.mat] == IdentityMatrix[max]];
	rev = simp @ rev;
	AbortProtect[invmat[r, s, t] = rev; Do[bas[r, s, t][n] = res[n], {n, max}]]
]

dec[r1_, r2_, r3_, r4_][f_] := Module[{vec, t = {r1, r2, r3, r4}},
	decPrep @@ t;
	vec = Array[{f @@ ((bas @@ t)[#])} &, Length[invs @@ t]];
	simp @ Flatten[(invmat @@ t).vec]]
x:decPrep[r1_, r2_, r3_, r4_] := x =
Module[{f = cor[r1, r2, r3, r4], mat, rev, new, old, sc, a1, a2, a3, a4, t,
	ks = invs[r1, r2, r3, r4], k, res, tmp, tmp2},
	Do[t = {a1, a2, a3, a4};
		If[(new = (f @@ First[ks]) @@ t) != 0, res[1] = t; mat = {{new}}; rev = Inverse[mat]; Break[]]
	, {a1, dim[r1]}, {a2, dim[r2]}, {a3, dim[r3]}, {a4, dim[r4]}];
	Do[Do[
		t = {a1, a2, a3, a4}; k = ks[[n]];
		new = Array[{(f @@ k) @@ res[#]} &, n - 1];
		old = {Array[(f @@ ks[[#]]) @@ t &, n - 1]};
		sc = (f @@ k) @@ t;
		tmp2 = rev.new;
		If[(tmp = simp[sc - (old.tmp2)[[1, 1]]]) != 0,
			res[n] = t;
			mat = ArrayFlatten[{{mat, new}, {old, {{sc}}}}];
			rev = ArrayFlatten[{{rev + tmp2.old.rev / tmp, -tmp2 / tmp}, {-old.rev / tmp, {{1 / tmp}}}}];
			Break[]]
	, {a1, dim[r1]}, {a2, dim[r2]}, {a3, dim[r3]}, {a4, dim[r4]}], {n, 2, Length[ks]}];
	rev = simp @ rev;
	AbortProtect[invmat[r1, r2, r3, r4] = rev; Do[bas[r1, r2, r3, r4][n] = res[n], {n, Length[ks]}]]
]

\[Sigma][r_, i_, r_][1] /; i === id := 1
\[Sigma][i_, r_, r_][1] /; i === id := 1
x:\[Sigma][r_, s_, t_][n_] := (setOPE[r, s, dual[t]]; x)
x:\[Sigma][r_] := x = \[Sigma][r, dual[r], id][1]

\[Sigma][op[_, _, p : 1 | -1, 1 | -1]] := p
\[Sigma][x_Times] := \[Sigma] /@ x
\[Sigma][x_^n_] := \[Sigma][x]^n
\[Sigma][1] = 1

x:setTau[r_, s_, t_] := x = Module[{n, m, l = inv[r, s, t], v, f, res},
	Do[
		f[a_, b_, c_] := cor[s, t, r][m][b, c, a];
		v = dec[r, s, t][f];
		Do[res[n, m] = v[[n]], {n, l}]
	, {m, l}];
	AbortProtect @ Do[\[Tau][r, s, t][n, m] = res[n, m], {m, l}, {n, l}]]
x:\[Tau][r_, s_, t_][n_, m_] /; 1 <= n <= inv[r, s, t] && 1 <= m <= inv[r, s, t] := (setTau[r, s, t]; x)
\[Tau][_, _, _][n_, n_] := 1
\[Tau][_, _, _][_, _] := 0

x:setOmega[r_, s_, t_] := x = Module[{n, m, l = inv[r, s, t], v, f, res},
	Do[
		f[a_, b_, c_] := simp @
			Sum[ope[dual[r]][a2, a] ope[dual[s]][b2, b] ope[dual[t]][c2, c] Conjugate[cor[dual[r], dual[s], dual[t]][m][a2, b2, c2]]
			, {a2, dim[r]}, {b2, dim[s]}, {c2, dim[t]}];
		v = dec[r, s, t][f];
		Do[res[n, m] = v[[n]], {n, l}]
	, {m, l}];
	AbortProtect @ Do[\[Omega][r, s, t][n, m] = res[n, m], {m, l}, {n, l}]]
x:\[Omega][r_, s_, t_][n_, m_] /; 1 <= n <= inv[r, s, t] && 1 <= m <= inv[r, s, t] := (setOmega[r, s, t]; x)
\[Omega][_, _, _][n_, n_] := 1
\[Omega][_, _, _][_, _] := 0

x:setSix[r1_, r2_, r3_, r4_] := x =
Module[{f, ss = invs[r1, r2, r3, r4], ts = invs[r1, r4, r3, r2], t, k, l, tkl, s, n, m, v, i, res},
	Do[
		{t, k, l} = tkl;
		f[a1_, a2_, a3_, a4_] := cor[r1, r4, r3, r2][t, k, l][a1, a4, a3, a2];
		v = dec[r1, r2, r3, r4][f];
		Do[res[tkl, i] = v[[i]], {i, Length[ss]}]
	, {tkl, ts}];
	AbortProtect @ Do[{t, k, l} = tkl; {s, n, m} = ss[[i]];
		six[r1, r2, r3, r4][s, n, m, t, k, l] = res[tkl, i], {tkl, ts}, {i, Length[ss]}]]
x:six[r1_, r2_, r3_, r4_][s_, n_, m_, t_, k_, l_] /;
1 <= n <= inv[{r1, r2}, {s}] && 1 <= m <= inv[r3, r4, s] && 1 <= k <= inv[{r1, r4}, {t}] && 1 <= l <= inv[r3, r2, t] :=
	(setSix[r1, r2, r3, r4]; x)
six[_, _, _, _][_, _, _, _, _, _] := 0

simp[x_] := FullSimplify @ Expand @ x

isReal[r_] := dual[r] === r && \[Sigma][r] == 1
isComplex[r_] := dual[r] =!= r && \[Sigma][r] == 1
isPseaudo[r_] := dual[r] === r && \[Sigma][r] == -1
dualOp[op[x_, r_, p: 1 | -1, q: 1 | -1]] := If[isPseaudo[r], op[x, r, -p, q], op[x, dual[r], p, q]]

reim[x_ == y_] := ComplexExpand[Re[x] == Re[y] && Im[x] == Im[y]]
reim[x_List] := And @@ reim /@ x
reim[a_And] := reim /@ a
reim[True] = True
eqList[x_ == y_] := {x == y}
eqList[True] = {}
eqList[a_And] := List @@ a
sum[a_?NumericQ b_, x_op] := a sum[b, x]
sum[a_Plus, x_op] := Plus @@ (sum[#, x] &) /@ List @@ a
sum[0, x_op] := 0
single[a_?NumericQ b_] := a single[b]
single[a_Plus] := Plus @@ single /@ List @@ a
single[0] = 0

flip12[L[L[o1 : op[_, r_, 1 | -1, 1], o2 : op[_, s_, 1 | -1, 1], o3 : op[_, t_, 1 | -1, q : 1 | -1], n_], val_]] := L[L[o2, o1, o3, n], q \[Sigma][r, s, t][n] val]
SetAttributes[flip12, Listable]
rotate[Ls_List] := Module[{N = Length[Ls], n, m, r, s, t, x = List @@ Take[Ls[[1, 1]], 3], tau, o1, o2, o3, sg},
	{r, s, t} = dual[#[[2]]] & /@ x;
	{d1, d2, d3} = dualOp /@ ({o1, o2, o3} = x);
	sg = \[Sigma][dualOp[o2] dualOp[o3]];
	tau = \[Tau][t, r, s];
	Table[L[L[o3, o1, o2, n], simp[sg Sum[tau[n, m] Ls[[m, 2]], {m, N}]]], {n, N}]
]
flipim[Ls_List] := Module[{N = Length[Ls], n, m, r, s, t, x = List @@ Take[Ls[[1, 1]], 3], omega, o1, o2, o3, d1, d2, d3, sg},
	{r, s, t} = dual[#[[2]]] & /@ x;
	{d1, d2, d3} = dualOp /@ ({o1, o2, o3} = x);
	sg = \[Sigma][o1 o2 d3];
	omega = \[Omega][r, s, t];
	Table[L[L[d1, d2, d3, n], simp @ ComplexExpand[sg Sum[omega[n, m] Conjugate[Ls[[m, 2]]], {m, N}]]], {n, N}]
]

setOPE[o1 : op[_, _, 1 | -1, 1], o2 : op[_, _, 1 | -1, 1], o3 : op[_, _, -1, 1 | -1]] := setOPE[dualOp[o1], dualOp[o2], dualOp[o3]]
setOPE[o1 : op[_, _, -1, 1], o2 : op[_, _, 1 | -1, 1], o3 : op[_, t_, 1, 1 | -1]] /; ! isPseaudo[t] := setOPE[dualOp[o1], dualOp[o2], dualOp[o3]]
setOPE[o1 : op[_, r_, 1, 1], o2 : op[_, _, -1, 1], o3 : op[_, t_, 1, 1 | -1]] /; ! isPseaudo[t] && ! isPseaudo[r] := setOPE[dualOp[o1], dualOp[o2], dualOp[o3]]
setOPE[o1 : op[_, r_, 1 | -1, 1], o2 : op[_, s_, 1 | -1, 1], o3 : op[_, t_, 1 | -1, 1 | -1]] :=
	Module[{re, im, a, b, c, d, e, v, cnt, d1 = dualOp[o1], d2 = dualOp[o2], d3 = dualOp[o3], n, N = inv[r, s, t], sol, add, eq, ev0, res},
		add[x_] := Sow[Join[Array[D[x, re[#]] &, N], Array[D[x, im[#]] &, N]]];
		a = Array[L[L[o1, o2, o3, #], re[#] + I im[#]] &, N];
		b = flip12[a];
		c = flipim[a];
		d = flip12[c];
		cnt = GroupBy[Join[a, b, c, d], First];
		eq = MyReap[Scan[Scan[add, zero @ reim[# == 0 & /@ Differences[#[[2]] & /@ #]]] &, cnt]];
		sol = NullSpace @ If[eq === {}, {ConstantArray[0, 2 N]}, eq];
		sol = Array[\[Beta][o1, o2, o3][#] &, Length[sol]] . sol;
		If[sol === 0, sol = ConstantArray[0, 2 N]];
		Do[re[n] = sol[[n]], {n, N}];
		Do[im[n] = sol[[n + N]], {n, N}];
		ev0[L[a_, b_, c_, n_]] := res[a, b, c, n] = simp @ cnt[L[a, b, c, n]][[1, 2]];
		ev[L[a_, b_, c_, n_]] := \[Alpha][a, b, c][n] = res[a, b, c, n];
		Scan[ev0, Keys[cnt]];
		AbortProtect @ Scan[ev, Keys[cnt]]]
setOPES[o1 : op[_, r_, 1 | -1, 1], o2 : op[_, s_, 1 | -1, 1], o3 : op[_, t_, 1 | -1, 1]] :=
	Module[{re, im, eq, a123, b123, e, v, cnt, x, N = inv[r, s, t], a213, a132, a231, a312, a321, b213, b132, b231, b312, b321, sol, add, n, ev0, res},
		add[x_] := Sow[Join[Array[D[x, re[#]] &, N], Array[D[x, im[#]] &, N]]];
		a123 = Array[L[L[o1, o2, o3, #], re[#] + I im[#]] &, N];
		a312 = rotate[a123];
		a231 = rotate[a312];
		a213 = flip12[a123];
		a132 = flip12[a312];
		a321 = flip12[a231];
		b123 = flipim[a123];
		b312 = rotate[b123];
		b231 = rotate[b312];
		b213 = flip12[b123];
		b132 = flip12[b312];
		b321 = flip12[b231];
		cnt = GroupBy[Join[a123, a213, a132, a231, a312, a321, b123, b213, b132, b231, b312, b321], First];
		eq = MyReap[Scan[Scan[add, zero @ reim[# == 0 & /@ Differences[#[[2]] & /@ #]]] &, cnt]];
		sol = NullSpace @ If[eq === {}, {ConstantArray[0, 2 N]}, eq];
		sol = Array[\[Beta][o1, o2, o3][#] &, Length[sol]] . sol;
		If[sol === 0, sol = ConstantArray[0, 2 N]];
		Do[re[n] = sol[[n]], {n, N}];
		Do[im[n] = sol[[n + N]], {n, N}];
		ev0[L[a_, b_, c_, n_]] := res[a, b, c, n] = simp @ cnt[L[a, b, c, n]][[1, 2]];
		ev[L[a_, b_, c_, n_]] := \[Nu][a, b, c][n] = res[a, b, c, n];
		Scan[ev0, Keys[cnt]];
		AbortProtect @ Scan[ev, Keys[cnt]]]
\[Lambda][o1_, o2_, o3_][n_] := \[Alpha][o1, o2, dualOp[o3]][n] / Sqrt[dim[o3[[2]]]]
x:\[Alpha][o1 : op[_, r_, 1 | -1, 1], o2 : op[_, s_, 1 | -1, 1], o3 : op[_, t_, 1 | -1, 1 | -1]][n_] /; 1 <= n <= inv[r, s, t] := (setOPE[o1, o2, o3]; x)
\[Alpha][_, _, _][_] := 0
\[Mu][o1_, o2_, o3_][n_] := \[Nu][o1, o2, dualOp[o3]][n] / Sqrt[dim[o3[[2]]]]
x:\[Nu][o1 : op[_, r_, 1 | -1, 1], o2 : op[_, _, 1 | -1, 1], o3 : op[0, _, 1, 1]][1] /; o3 === one := x = If[o1 === dualOp[o2], Sqrt[dim[r]]\[Sigma][o1], 0]
x:\[Nu][o1 : op[_, r_, 1 | -1, 1], o2 : op[0, _, 1, 1], o3 : op[_, _, 1 | -1, 1]][1] /; o2 === one := x = If[o1 === dualOp[o3], Sqrt[dim[r]], 0]
x:\[Nu][o1 : op[0, _, 1, 1], o2 : op[_, r_, 1 | -1, 1], o3 : op[_, _, 1 | -1, 1]][1] /; o1 === one := x = If[o2 === dualOp[o3], Sqrt[dim[r]], 0]
x:\[Nu][o1 : op[_, r_, 1 | -1, 1], o2 : op[_, s_, 1 | -1, 1], o3 : op[_, t_, 1 | -1, 1]][n_] /; 1 <= n <= inv[r, s, t] && KeyExistsQ[allops, o3] := (setOPES[o1, o2, o3]; x)
\[Nu][_, _, _][_] := 0

x:\[Alpha][o1_, o2_, o3_, o4_][o_, n_, m_] := x = \[Lambda][o1, o2, o][n] \[Alpha][o3, o4, o][m]
x:\[Nu][o1_, o2_, o3_, o4_][o_, n_, m_] := x = \[Mu][o1, o2, o][n] \[Nu][o3, o4, o][m]

x:sF[o1_, o2_, o3_, o4_][s_, n_, m_, p : 1 | -1] /; ! isPseaudo[s] && s === minrep[s, dual[s]] := x =
	Simplify[Module[{sp}, Sum[With[{o = op[op, s, 1, sp]}, sum[ComplexExpand[\[Alpha][o1, o2, o3, o4][o, n, m] F[{o1, o2, o3, o4}, p]], o]], {sp, {1, -1}}]]]
x:sF[o1_, o2_, o3_, o4_][s_, n_, m_, p : 1 | -1] /; s =!= minrep[s, dual[s]] := x =
	Simplify[Module[{sp},
		Sum[With[{o = op[op, dual[s], 1, sp], do = op[op, s, 1, sp]}, sum[ComplexExpand[\[Alpha][o1, o2, o3, o4][do, n, m] F[{o1, o2, o3, o4}, p]], o]], {sp, {1, -1}}]]]
x:sF[o1_, o2_, o3_, o4_][s_, n_, m_, p : 1 | -1] /; isPseaudo[s] := x =
	Simplify[Module[{sp}, Sum[With[{o = op[op, s, 1, sp], do = op[op, s, -1, sp]},
			sum[ComplexExpand[(\[Alpha][o1, o2, o3, o4][o, n, m] + \[Alpha][o1, o2, o3, o4][do, n, m]) F[{o1, o2, o3, o4}, p]], o]], {sp, {1, -1}}]]]

x:sFp[o1_, o2_, o3_, o4_][s_, n_, m_, p : 1 | -1] /; KeyExistsQ[allreps, s] := x =
	Module[{o}, Sum[single[ComplexExpand[\[Nu][o1, o2, o3, o4][o, n, m] Fp[{o1, o2, o3, o4}, p, o]]], {o, Keys[allreps[s]]}]]
sFp[_, _, _, _][_, _, _, 1 | -1] := 0

sF2[o1_, o2_, o3_, o4_][s_, n_, m_, p : 1 | -1] := sF[o1, o2, o3, o4][s, n, m, p] + sFp[o1, o2, o3, o4][s, n, m, p]

boot[o1 : op[_, r1_, 1 | -1, 1], o2 : op[_, r2_, 1 | -1, 1], o3 : op[_, r3_, 1 | -1, 1], o4 : op[_, r4_, 1 | -1, 1]][s_, n_, m_, p : 1 | -1] :=
	Module[{t, k, l}, 0 == sF2[o1, o2, o3, o4][s, n, m, p]
		- p Sum[six[r1, r2, r3, r4][s, n, m, t, k, l] sF2[o1, o4, o3, o2][t, k, l, p], {t, Keys[prod[r1, r4]]}, {k, inv[{r1, r4}, {t}]}, {l, inv[r3, r2, t]}]]
boot[o1 : op[_, r1_, 1 | -1, 1], o2 : op[_, r2_, 1 | -1, 1], o3 : op[_, r3_, 1 | -1, 1], o4 : op[_, r4_, 1 | -1, 1]][s_] :=
	Reduce @ reim @ MyReap @ Module[{n, m, p}, Do[Sow[boot[o1, o2, o3, o4][s, n, m, p]], {n, inv[{r1, r2}, {s}]}, {m, inv[r3, r4, s]}, {p, {1, -1}}]]
bootAll[ops_List] := eqn @ zero @ Reduce @ reim @ MyReap @ Module[{a, b, c, d, s}, Do[Sow[boot[a, b, c, d][s]], {a, ops}, {b, ops}, {s, prodOp[a, b]}, {c, ops}, {d, ops}]]
bootAll[] := bootAll @ keys @ allopsForBoot

format[F[x_, y_, z_, w_]] := Superscript[F, {x, y, z, w}]
format[H[x_, y_, z_, w_]] := Superscript[H, {x, y, z, w}]
format[Fp[x_, y_, z_, w_, o_]] := Subsuperscript[F, o, {x, y, z, w}]
format[Hp[x_, y_, z_, w_, o_]] := Subsuperscript[H, o, {x, y, z, w}]
format[\[Beta][o1_, o2_, o3_][n_]] := Subscript[\[Beta], format[o1], format[o2], format[o3], n]
format[sum[x_, op[op, r_, 1, p : 1 | -1]]] := Underscript["\[CapitalSigma]", Element[op, sign[p, r]]][format[x]]
format[single[x_]] := format[x]
format[x_Association] := Association[format[#] -> format[x[#]] & /@ Keys[x]]
format[x_[y___]] := x @@ format /@ List @@ x[y]
format[x_] := x
format[op[x_, r_, p : 1 | -1, 1 | -1]] /; !isComplex[r] := tmp[p x]
format[op[x_, r_?isComplex, 1, 1 | -1]] := tmp[If[r === minrep[r, dual[r]], 1, -1] x]
tmp[-x_] := OverBar[x]
tmp[x_] := x
sign[1, x_] := SuperPlus[x]
sign[-1, x_] := SuperMinus[x]

extract2[a_Plus, o_] := extract2[#, o] & /@ a
extract2[x_?NumericQ y_, o_] := x extract2[y, o]
extract2[sum[x_, o_op], o_] := sum[x, o]
extract2[sum[_, _op], _] := 0
extract2[x_single, unit | scalar] := x
extract2[_single, _op] := 0
extract2[single[(Fp | Hp)[_, _, _, _, 0]], scalar] := 0
extract2[single[Times[__, (Fp | Hp)[_, _, _, _, _]]], unit] := 0
extract2[0, _] := 0
extract[eqn[eq_List], o_] := eqn[o, extract2[#, o] & /@ eq]

F[{x__}, -1] := F[x]
F[{x__}, 1] := H[x]
F[op[x1_, _, 1 | -1, 1], op[x2_, _, 1 | -1, 1], op[x3_, _, 1 | -1, 1], op[x4_, _, 1 | -1, 1]] := F[x1, x2, x3, x4]
H[op[x1_, _, 1 | -1, 1], op[x2_, _, 1 | -1, 1], op[x3_, _, 1 | -1, 1], op[x4_, _, 1 | -1, 1]] := H[x1, x2, x3, x4]
Fp[{x__}, -1, o: op[_, _, 1 | -1, 1]] := Fp[x, o]
Fp[{x__}, 1, o: op[_, _, 1 | -1, 1]] := Hp[x, o]
Fp[op[x1_, _, 1 | -1, 1], op[x2_, _, 1 | -1, 1], op[x3_, _, 1 | -1, 1], op[x4_, _, 1 | -1, 1], op[o_, _, 1 | -1, 1]] := Fp[x1, x2, x3, x4, o]
Hp[op[x1_, _, 1 | -1, 1], op[x2_, _, 1 | -1, 1], op[x3_, _, 1 | -1, 1], op[x4_, _, 1 | -1, 1], op[o_, _, 1 | -1, 1]] := Hp[x1, x2, x3, x4, o]
prodOp[op[_, r_, 1 | -1, 1], op[_, s_, 1 | -1, 1]] := Keys @ prod[r, s]
prodOp[op[_, r_, 1 | -1, 1], s_] := Keys @ prod[r, s]
prodOp[r_, op[_, s_, 1 | -1, 1]] := Keys @ prod[r, s]

tmpsec = newSet[]
tmpf[a_ + b_] := (tmpf[a]; tmpf[b];)
tmpf[a_?NumericQ b_] := (tmpf[b];)
tmpf[sum[x_, o_op]] := add[tmpsec, o]
sector[eqn[eq_List]] := Module[{pr},
	clear[tmpsec];
	add[tmpsec, unit];
	add[tmpsec, scalar];
	Scan[tmpf, eq];
	pr[unit, _] := True;
	pr[_, unit] := False;
	pr[scalar, _] := True;
	pr[_, scalar] := False;
	pr[op[op, r_, 1, p_], op[op, s_, 1, q_]] := (r === s && p > q) || minrep[r, s] =!= s;
	Sort[keys[tmpsec], pr]]

F[a_, b_, a_, b_] /; ord[b] < ord[a] := F[b, a, b, a]
H[a_, b_, a_, b_] /; ord[b] < ord[a] := H[b, a, b, a]
F[a_, b_, c_, d_] /; ord[c] < ord[a] || (ord[c] == ord[a] && ord[d] < ord[b]) := F[c, d, a, b]
H[a_, b_, c_, d_] /; ord[c] < ord[a] || (ord[c] == ord[a] && ord[d] < ord[b]) := H[c, d, a, b]
Fp[a_, b_, a_, b_, o_] /; ord[b] < ord[a] := Fp[b, a, b, a, o]
Hp[a_, b_, a_, b_, o_] /; ord[b] < ord[a] := Hp[b, a, b, a, o]
Fp[a_, b_, c_, d_, o_] /; ord[c] < ord[a] || (ord[c] == ord[a] && ord[d] < ord[b]) := Fp[c, d, a, b, o]
Hp[a_, b_, c_, d_, o_] /; ord[c] < ord[a] || (ord[c] == ord[a] && ord[d] < ord[b]) := Hp[c, d, a, b, o]
addOp[o : op[x_, r_, 1 | -1, 1]] := AbortProtect[allops[o] = 1; add[allopsForBoot, o];
	If[!KeyExistsQ[allreps, r], allreps[r] = <|o->1|>, allreps[r][o] = 1];
	If[!KeyExistsQ[ord, x], ord[x] = Length[ord]];]
setOps[ops_List] := AbortProtect[allops = <|one->1|>; allreps = <|id-><|one->1|>|>; ord = <||>; allopsForBoot = newSet[]; Scan[(addOp[#]; addOp[dualOp[#]]) &, ops]]

connectedComponents[eqn[_, z_]] :=
	Module[{f, uf = newUF[]},
		f[a_Plus] := Scan[f, a];
		f[x_?NumericQ y_] := f[y];
		f[sum[(_F | _H) a_ b_, _]] := unite[uf, a, b];
		f[sum[(_F | _H) a_^2, _]] := add[uf, a];
		f[single[(_Fp | _Hp) a_ b_]] := unite[uf, a, b];
		f[single[(_Fp | _Hp) a_^2]] := add[uf, a];
		Scan[f, z];
		classify[uf]]

makeMatHelper[z_, con_] :=
	Module[{ind, size = Length[con], v, mat, f, add},
		Do[ind[con[[v]]] = v, {v, size}];
		ind[_] := 0;
		mat = Array[0&, {size, size}];
		f[a_Plus] := Scan[f, a];
		f[x_?NumericQ y_] := f[x, y];
		f[x : _sum | _single] := f[1, x];
		add[_, _, ___, 0, ___] := Null;
		add[c_, x_, n_, m_] := (mat[[n, m]] += c x / 2; mat[[m, n]] += c x / 2;);
		add[c_, x_, n_] := (mat[[n, n]] += c x;);
		f[c_, sum[(x:_F | _H) a_ b_, _]] := add[c, x, ind[a], ind[b]];
		f[c_, sum[(x:_F | _H) a_^2, _]] := add[c, x, ind[a]];
		f[c_, single[(x:_Fp | _Hp) a_ b_]] := add[c, x, ind[a], ind[b]];
		f[c_, single[(x:_Fp | _Hp) a_^2]] := add[c, x, ind[a]];
		f[z];
		mat
	]

makeMat[z:eqn[_, eq_List]] := Module[{con = connectedComponents[z], c},
	If[Length[con] == 0, {{{# /. single[x_] :> x}} & /@ eq}, Table[makeMatHelper[#, c]& /@ eq, {c, con}]]]

makeSDP[z:eqn[eq_List]] := Module[{mat, tmp, sec = sector[z], secs = <||>, scalnum, val = newSet[], f, s, i},
		f[x_Plus] := f /@ x;
		f[0] = 0;
		f[x_?NumericQ y_] := If[x > 0, add[val, x]; block[1, x, y], If[x != -1, add[val, -x]]; block[-1, -x, y]];
		f[y_] := block[1, 1, y];
		SetAttributes[f, Listable];
		mat[unit, 1] = f[makeMat[extract[z, unit]][[1]]];
		tmp = makeMat[extract[z, scalar]];
		scalnum = Length[tmp];
		Do[mat[scalar, i] = f[tmp[[i]]], {i, scalnum}];
		Do[tmp = makeMat[extract[z, s]]; secs[s] = Length[tmp]; Do[mat[s, i] = f[tmp[[i]]], {i, secs[s]}];, {s, Drop[sec, 2]}];
		sdp[secs, scalnum, keys[val], mat]
	]

toString[F[x_, y_, z_, w_]] := TemplateApply["get(F, \"`x`\", \"`y`\", \"`z`\", \"`w`\")", <|"x"->x, "y"->y, "z"->z, "w"->w|>];
toString[H[x_, y_, z_, w_]] := TemplateApply["get(H, \"`x`\", \"`y`\", \"`z`\", \"`w`\")", <|"x"->x, "y"->y, "z"->z, "w"->w|>];
toString[Fp[x_, y_, z_, w_, 0]] := TemplateApply["get(F, \"`x`\", \"`y`\", \"`z`\", \"`w`\")", <|"x"->x, "y"->y, "z"->z, "w"->w|>];
toString[Hp[x_, y_, z_, w_, 0]] := TemplateApply["get(H, \"`x`\", \"`y`\", \"`z`\", \"`w`\")", <|"x"->x, "y"->y, "z"->z, "w"->w|>];
toString[Fp[x_, y_, z_, w_, o_]] := TemplateApply["get(F, \"`x`\", \"`y`\", \"`z`\", \"`w`\", \"`o`\")", <|"x"->x, "y"->y, "z"->z, "w"->w, "o"->o|>];
toString[Hp[x_, y_, z_, w_, o_]] := TemplateApply["get(H, \"`x`\", \"`y`\", \"`z`\", \"`w`\", \"`o`\")", <|"x"->x, "y"->y, "z"->z, "w"->w, "o"->o|>];

(*
secs = <|op[op, rep[1], 1, 1] -> 3, ...|>
scalarnum = 2
vals = {1/2, Sqrt[2], Pi, ...}
mat[sec, num] = {{{block[-1,Sqrt[2],F[e,v,e,v]]+block[1,Sqrt[3]/2,F[e,e,v,v]], ...}, ...}, ...}
*)
toString[sdp[secs_, scalarnum_, vals_, mat_]] :=
	Module[{s, v, o, revval = reverseIndex[vals], blk, f, convert, make, tmp, secsStr, valsStr, rmats, smats, umats, filename},
		secsStr = StringRiffle[TemplateApply["(\"`r`\", `p`): `n`", <|"r" -> #[[2]], "p" -> (1 - #[[4]])/2, "n" -> secs[#]|>] & /@ Keys[secs], ", "];
		valsStr = tensorToString[vals, ToPython`pyeval];
		rmats = newSet[];
		smats = newSet[];
		filename = StringRiffle[TemplateApply["{0[`k`]}", <|"k" -> #|>] & /@ DeleteDuplicates[#[[1]] & /@ keys[allopsForBoot]], "-"];
		convert[0] = "0";
		convert[block[1, 1, bl_]] := TemplateApply["bl[`i`]", <|"i" -> blk[bl]|>];
		convert[block[-1, 1, bl_]] := TemplateApply["-bl[`i`]", <|"i" -> blk[bl]|>];
		convert[block[1, v_, bl_]] := TemplateApply["val[`j`] * bl[`i`]", <|"i" -> blk[bl], "j" -> revval[v]|>];
		convert[block[-1, v_, bl_]] := TemplateApply["-val[`j`] * bl[`i`]", <|"i" -> blk[bl], "j" -> revval[v]|>];
		convert[x_Plus] := convert[x[[1]]] <> StringJoin[Table[If[y[[1]] > 0, " + " <> convert[y], " - " <> convert[block[1, y[[2]], y[[3]]]]], {y, Rest[x]}]];
		f[x_List] := Scan[f, x];
		f[x_Plus] := Scan[f, List @@ x];
		f[block[_, _, bl_]] := If[!KeyExistsQ[blk, bl], blk[bl] = Length[blk];];
		make[m_, init_] := Module[{tmp = <||>},
			Scan[(tmp[#] = init[#]) &, Keys[init]];
			blk = <||>; f[m];
			tmp["bl"] = tensorToString[Keys[blk], toString];
			tmp["mats"] = If[Length[m[[1]]] != 1, tensorToString[m, convert], tensorToString[#[[1, 1]] & /@ m , convert]];
			tmp];
		Do[add[rmats, make[mat[s, n], <|"r" -> s[[2]], "p" -> (1 - s[[4]])/2, "n" -> n - 1|>]], {s, Keys[secs]}, {n, secs[s]}];
		Do[add[smats, make[mat[scalar, n], <|"n" -> n - 1|>]], {n, scalarnum}];
		rmats = StringRiffle[TemplateApply[rmatstemplate, #] & /@ keys[rmats], "\n"];
		smats = StringRiffle[TemplateApply[smatstemplate, #] & /@ keys[smats], "\n"];
		umats = TemplateApply[umatstemplate, make[mat[unit, 1], <||>]];
		ToPython`createPython[secsStr, scalarnum, valsStr, rmats, smats, umats, filename]
	]

tensorToString[tensor_List, conv_] := "[" <> StringRiffle[tensorToString[#, conv] & /@ tensor, ", "] <> "]"
tensorToString[tensor_, conv_] := conv[tensor]

rmatstemplate = "	if sector == \"`r`\" and spin % 2 == `p` and num == `n`:
		bl = `bl`
		return `mats`"

smatstemplate = "	if num == `n`:
		bl = `bl`
		return `mats`"

umatstemplate = "bl = `bl`
	return `mats`"

makeG[eqn[sec_, eq_List]] :=
	Module[{vs = <||>, es = <||>, f, add, label},
		If[sec === unit, Return[ToString[unit, StandardForm]]];
		If[sec === scalar, label = scalar, label = sign[sec[[4]], sec[[2]]]];
		add[a_] := vs[a] = 1;
		add[a_, a_] := add[a];
		add[a_, b_] := (vs[a] = 1; vs[b] = 1; es[Sort[{a, b}]] = 1);
		f[a_Plus] := Scan[f, List @@ a];
		f[_?NumericQ y_] := f[y];
		f[sum[(_F | _H) a_ b_, _]] := add[a, b];
		f[single[(_Fp | _Hp) a_ b_]] := add[a, b];
		f[sum[(_F | _H) a_^2, _]] := add[a, a];
		f[single[(_Fp | _Hp) a_^2]] := add[a, a];
		Scan[f, eq];
		Graph[Labeled[#, format[#]] & /@ Keys[vs], #[[1]] <-> #[[2]] & /@ Keys[es], PlotLabel -> ToString[label, StandardForm]]]

End[ ]

EndPackage[ ]
