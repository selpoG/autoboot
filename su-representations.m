(* Shared exact construction; no GroupInfo / NGroupInfo symbols are imported. *)
Needs["RootSystem`", "root.m"]
BeginPackage["SURepresentations`"]
validWeight::usage = "validWeight[n, rows] tests SU(n) Young diagram row lengths."
matrices::usage = "matrices[n, rows] constructs lowering then raising generators in an orthonormal basis."
matrices::dimension = "SU(`1`) representation `2` generated `3` basis vectors; expected `4`."
Begin["`Private`"]

validWeight[n_Integer, rows_List] := Length[rows] == n - 1 &&
	AllTrue[rows, IntegerQ[#] && # >= 0 &] && Reverse[Sort[rows]] === rows;

(* Exterior powers of the defining representation. Adjacent replacements
   require no permutation sign in the increasing-subset basis. *)
fundamental[n_, k_] := fundamental[n, k] = Module[{subsets, lower, rules, target},
	subsets = Subsets[Range[n], {k}];
	lower = Table[
		rules = CommonFunctions`MyReap @ Do[
			If[MemberQ[subsets[[col]], i] && !MemberQ[subsets[[col]], i + 1],
				target = subsets[[col]] /. i -> i + 1;
				Sow[{First @ FirstPosition[subsets, target], col} -> 1]]
		, {col, Length[subsets]}];
		SparseArray[rules, {Length[subsets], Length[subsets]}]
	, {i, n - 1}];
	Join[lower, Transpose /@ lower]
];

matrices[n_Integer, rows_List] /; n >= 3 && validWeight[n, rows] := Module[{result},
	result = construct[n, rows];
	If[result === $Failed, $Failed, matrices[n, rows] = result]
];

construct[n_, rows_] := Module[
	{rank = n - 1, counts, factors, dims, size, expected, generators,
	 layer, basis, candidates, reduced, total = 0, result},
	counts = rows - Append[Rest[rows], 0];
	factors = Flatten[MapIndexed[ConstantArray[First[#2], #1] &, counts]];
	If[factors === {}, Return[ConstantArray[SparseArray[{}, {1, 1}], 2 rank]]];
	dims = Binomial[n, #] & /@ factors;
	size = Times @@ dims;
	expected = RootSystem`dimension["A", rank, rows];
	generators = Table[
		Sum[KroneckerProduct[
			IdentityMatrix[Times @@ Take[dims, j - 1], SparseArray],
			fundamental[n, factors[[j]]][[i]],
			IdentityMatrix[Times @@ Drop[dims, j], SparseArray]]
		, {j, Length[factors]}]
	, {i, 2 rank}];
	layer = {SparseArray[{1 -> 1}, size]};
	(* Each lowering step changes the height by one, so different layers are
	   orthogonal. Row reduction removes dependencies within each layer. *)
	basis = CommonFunctions`MyReap[
		While[layer =!= {} && total < expected,
			Scan[Sow, SparseArray /@ Orthogonalize[Normal /@ layer]];
			total += Length[layer];
			If[total >= expected, Break[]];
			candidates = Flatten[Table[generators[[i]].vec,
				{vec, layer}, {i, rank}], 1];
			reduced = Select[Normal @ RowReduce[SparseArray[candidates]],
				AnyTrue[#, !TrueQ[# == 0] &] &];
			(* Keep the frontier rational; normalize only the saved basis. *)
			layer = SparseArray /@ reduced;
		]
	];
	If[total != expected,
		Message[matrices::dimension, n, rows, total, expected]; Return[$Failed]];
	basis = SparseArray[basis];
	result = Simplify[Conjugate[basis].#.Transpose[basis]] & /@ generators;
	result
];
End[]
EndPackage[]
