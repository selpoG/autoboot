Needs["CommonFunctions`", "common.m"]
BeginPackage["RootSystem`"]

dimension::usage = "dimension[t,r,l] gives the dimension of the irrep with the highest weight l. The Lie algebra is of type t and rank r."
irrep::usage = "irrep[t,r,l] gives the character object of the irrep with the highest weight l. The Lie algebra is of type t and rank r."
product::usage = "product[ch1,ch2] gives the character object of the product representation of two representation ch1 and ch2."
decompose::usage = "decompose[c] gives the irreducible decomposition of the character c."
ch::usage = "ch[t,r,x] represents a character object in the Lie algebra of type t and rank r. x is a list of w->d in which w is a weight and d is the dimension of the weight space with weight w."

Begin["`Private`"]

bilinear::usage = "bilinear[t,r][x,y] is a bilinear form of weight x and y."
lexComp::usage = "lexComp[x,y] compares two weights x and y. Let e be an infinitesimal, the weyl chamber is selected by the vector (e,e^2,e^3,...). If lexComp[x,y] > 0, x is higher than y."

allPublicSymbol = {dimension, irrep, product, decompose, ch}

protectCounter = 0
myUnprotect[expr_] := Module[{ans},
	Unprotect[Evaluate[allPublicSymbol]];
	protectCounter++;
	ans = expr;
	protectCounter--;
	If[protectCounter == 0, Protect[Evaluate[allPublicSymbol]]];
	ans]
myAbortProtect[expr_] := AbortProtect[myUnprotect[expr]]
SetAttributes[myUnprotect, HoldFirst]
SetAttributes[myAbortProtect, HoldFirst]

CommonFunctions`importPackage["CommonFunctions`", "RootSystem`Private`", {"importPackage", "MyReap", "newPQueue", "empty", "enqueue", "dequeue", "size", "top"}]

vec[dim_, k_] /; k <= dim := vec[dim, k] = Table[If[x == k, 1, 0], {x, dim}]
vec[dim_, k_] /; k == dim + 1 := vec[dim, k] = Table[-1, {x, dim}]

bilinear["A", rank_][x_, y_] := (rank + 1) x.y - Total[x] Total[y]
bilinear[_, _][x_, y_] := x.y
norm[type_, rank_][x_] := bilinear[type, rank][x, x]

lexComp[x_List, y_List] := Module[{i, k = 0}, Do[If[(k = Sign[x[[i]] - y[[i]]]) != 0, Break[]], {i, Length[x]}]; k]

x:simpleRoots["A", rank_] /; IntegerQ[rank] && rank >= 1 := x =
	Table[vec[rank, i] - vec[rank, i + 1], {i, rank}]
x:simpleRoots["B", rank_] /; IntegerQ[rank] && rank >= 2 := x =
	Join[Table[vec[rank, i] - vec[rank, i + 1], {i, rank - 1}], {vec[rank, rank]}]
x:simpleRoots["C", rank_] /; IntegerQ[rank] && rank >= 2 := x =
	Join[Table[vec[rank, i] - vec[rank, i + 1], {i, rank - 1}], {2 vec[rank, rank]}]
x:simpleRoots["D", rank_] /; IntegerQ[rank] && rank >= 2 := x =
	Join[Table[vec[rank, i] - vec[rank, i + 1], {i, rank - 1}], {vec[rank, rank - 1] + vec[rank, rank]}]

x:positiveRoots["A", rank_] /; IntegerQ[rank] && rank >= 1 := x =
	Flatten[Table[vec[rank, i] - vec[rank, j], {i, rank}, {j, i + 1, rank + 1}], 1]
x:positiveRoots["B", rank_] /; IntegerQ[rank] && rank >= 2 := x =
	Join[positiveRoots["D", rank], Table[vec[rank, i], {i, rank}]]
x:positiveRoots["C", rank_] /; IntegerQ[rank] && rank >= 2 := x =
	Join[positiveRoots["D", rank], Table[2 vec[rank, i], {i, rank}]]
x:positiveRoots["D", rank_] /; IntegerQ[rank] && rank >= 2 := x =
	Flatten[Table[vec[rank, i] + s vec[rank, j], {i, rank - 1}, {j, i + 1, rank}, {s, {1, -1}}], 2]

simpleRoots["E", 8] = {
	{0, 0, 0, 0, 0, 0, 1, -1}, {0, 0, 0, 0, 0, 1, -1, 0},
	{0, 0, 0, 0, 0, 0, 1, 1}, {0, 0, 0, 0, 1, -1, 0, 0},
	{0, 0, 0, 1, -1, 0, 0, 0}, {0, 0, 1, -1, 0, 0, 0, 0},
	{0, 1, -1, 0, 0, 0, 0, 0}, {1, -1, -1, -1, -1, -1, -1, 1}/2};

positiveRoots["E", 8] = Join[positiveRoots["D", 8],
	MyReap @ Do[
		If[Mod[DigitCount[i, 2, 1], 2] == 0, Sow[# - 1/2 & /@ IntegerDigits[i, 2, 8]]]
	, {i, 2^7, 2^8 - 1}]];

x:\[Rho][type_, rank_] := myAbortProtect[x = Total @ positiveRoots[type, rank] / 2]

dimension[type_, rank_, \[Lambda]_] :=
	Product[bilinear[type, rank][\[Lambda] + \[Rho][type, rank], a], {a, positiveRoots[type, rank]}] /
	Product[bilinear[type, rank][\[Rho][type, rank], a], {a, positiveRoots[type, rank]}]

x:irrep[type_, rank_, \[Lambda]_] := Module[{ans},
	ans = MyReap @ Module[{
		queue = newPQueue[lexComp], w, set = <||>, sum, \[Alpha], k,
		bi = bilinear[type, rank], nm = norm[type, rank], \[Rho] = \[Rho][type, rank], denom, add},
		add[k_] := If[!KeyMemberQ[set, k], enqueue[queue, k]; True, False];
		add[\[Lambda]];
		Monitor[While[size[queue] > 0,
			w = dequeue[queue];
			If[KeyMemberQ[set, w], Continue[]];
			If[w === \[Lambda],
				set[w] = 1; Sow[w -> 1],
				denom = nm[\[Lambda] + \[Rho]] - nm[w + \[Rho]];
				If[denom == 0, Continue[]];
				sum = 0;
				Do[
					k = \[Alpha] + w;
					While[KeyMemberQ[set, k], sum += set[k] bi[\[Alpha], k]; k += \[Alpha]]
				, {\[Alpha], positiveRoots[type, rank]}];
				If[sum <= 0, Continue[]];
				set[w] = 2 sum/denom;
				Sow[w -> set[w]]];
			Do[add[w - \[Alpha]], {\[Alpha], simpleRoots[type, rank]}]
			], TemplateApply["q=`q`, len=`l`, set=`s`, w=`w`", <|"q" -> queue, "l" -> size[queue], "s" -> set, "w" -> w|>]
		];];
	myAbortProtect[x = ch[type, rank, ans]]]

product[ch[type_, rank_, p_List], ch[type_, rank_, q_List]] :=
	ch[type, rank, Module[{set = <||>, x, y, z},
		Do[
			z = x[[1]] + y[[1]];
			If[!KeyMemberQ[set, z], set[z] = 0];
			set[z] += x[[2]] y[[2]]
		, {x, p}, {y, q}];
		# -> set[#] & /@ Sort[Keys[set], lexComp]]]

decompose[ch[type_, rank_, p_List]] :=
	MyReap @ Module[{set = <||>, x, y, c},
		Do[set[x[[1]]] = x[[2]], {x, p}];
		Do[
			If[set[x] == 0, Continue[]];
			c = set[x];
			Do[Sow[x], c];
			Do[set[y[[1]]] -= c y[[2]], {y, irrep[type, rank, x][[3]]}]
		, {x, First /@ p}]]

Protect[Evaluate[allPublicSymbol]]

End[ ]

EndPackage[ ]
