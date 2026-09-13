(* Finite-type highest-weight representations, in Dynkin coordinates. *)
Needs["CommonFunctions`", "common.m"]
BeginPackage["LieRepresentations`"]
supported::usage = "supported[type, rank] tests supported finite Cartan types."
validLabel::usage = "validLabel[rank, labels] tests nonnegative integral Dynkin labels."
cartan::usage = "cartan[type, rank] gives A_ij = <alpha_j, alpha_i coroot>."
repDimension::usage = "repDimension[type, rank, labels] evaluates the Weyl dimension formula."
repCharacter::usage = "repCharacter[type, rank, labels] gives weight multiplicities."
repProduct::usage = "repProduct[type, rank, a, b] gives highest weights with multiplicity."
repProduct::multiplicity = "Negative tensor multiplicities in (`1`, `2`, `3`, `4`)."
repDual::usage = "repDual[type, rank, labels] gives the dual highest weight."
repMatrices::usage = "repMatrices[type, rank, labels] gives lowering then raising matrices."
repMatrices::dimension = "Representation (`1`, `2`, `3`) generated `4` basis vectors; expected `5`."
intertwiners::usage = "Internal group capability for exact CG intertwiners in the group's generator basis."
repIntertwiners::usage = "repIntertwiners[type,rank,a,b,c] computes exact intertwiners using weight conservation."
Begin["`Private`"]
(* Large tensor layers repeat the same algebraic entries many times. Canonicalize
   each distinct scalar once, keeping zero entries implicit in sparse arrays. *)
scalarReduce[x_] := scalarReduce[x] = RootReduce[x];
arrayReduce[x_SparseArray] := SparseArray[(First[#] -> scalarReduce[Last[#]] &) /@ ArrayRules[x],Dimensions[x]];
arrayReduce[x_List] := arrayReduce /@ x;
arrayReduce[x_] := scalarReduce[x];

supported[t_, r_] := IntegerQ[r] && Switch[t,
	"A" | "C", r >= 1, "B" | "D", r >= 2,
	"E", MemberQ[{6,7,8}, r], "F", r == 4, "G", r == 2, _, False];
validLabel[r_, a_List] := Length[a] == r && AllTrue[a, IntegerQ[#] && # >= 0 &];
validLabel[_, _] := False;

(* Bourbaki numbering: E has chain 1-3-4-5-... with 2 attached to 4.
   B's final node is short, C's is long; G's first node is short. *)
cartan[t_, r_] /; supported[t,r] := cartan[t,r] = Module[{a = 2 IdentityMatrix[r], edges},
	edges = Switch[t,
		"D", Join[Table[{i,i+1}, {i,r-2}], {{r-2,r}}] /. {0,_} -> Nothing,
		"E", Join[{{1,3},{2,4}}, Table[{i,i+1}, {i,3,r-1}]],
		_, Table[{i,i+1}, {i,r-1}]];
	Do[a[[e[[1]],e[[2]]]] = a[[e[[2]],e[[1]]]] = -1, {e,edges}];
	Switch[t, "B", a[[r,r-1]] = -2, "C", If[r > 1, a[[r-1,r]] = -2],
		"F", a[[3,2]] = -2, "G", a[[1,2]] = -3];
	a
];

(* Reflection closure, stored as coefficients of simple roots. *)
positive[a_List] := positive[a] = Module[{roots = IdentityMatrix[Length[a]], next, old = {}},
	While[roots =!= old,
		old = roots;
		next = Flatten[Table[c - (a[[i]].c) UnitVector[Length[a],i], {c,roots}, {i,Length[a]}],1];
		roots = Union[roots,next]];
	Select[roots, AllTrue[#, # >= 0 &] &]
];

metric[t_, r_] := metric[t,r] = Module[{d},
	d = Switch[t, "B", Append[ConstantArray[2,r-1],1],
		"C", Append[ConstantArray[1,r-1],2], "F", {2,2,1,1},
		"G", {1,3}, _, ConstantArray[1,r]];
	DiagonalMatrix[d].Inverse[cartan[t,r]]
];

repDimension[t_, r_, a_List] /; supported[t,r] && validLabel[r,a] :=
	Times @@ (((a + 1).#/Total[#]) & /@ positive[Transpose[cartan[t,r]]]);

repDual[t_, r_, a_List] /; supported[t,r] && validLabel[r,a] := Module[{w = -a, i, c = cartan[t,r]},
	While[AnyTrue[w, # < 0 &], i = First @ FirstPosition[w, _?Negative]; w -= w[[i]] c[[All,i]]]; w
];

(* Freudenthal recursion, processing one root-height layer at a time. *)
repCharacter[t_, r_, a_List] /; supported[t,r] && validLabel[r,a] :=
 repCharacter[t,r,a] = Module[{c = cartan[t,r], m = metric[t,r], roots, rho = ConstantArray[1,r],
	weights = <||>, layer = {a}, next, w, denom, mult, sum, u, depth = 0},
	roots = ({c.#, Total[#]} &) /@ positive[c];
	While[layer =!= {},
		next = {};
		Do[
			denom = (a+rho).m.(a+rho) - (w+rho).m.(w+rho);
			mult = If[w === a, 1, If[denom == 0, 0,
				sum = Sum[Sum[u = w + k root[[1]];
					If[KeyMemberQ[weights,u], weights[u] (u.m.root[[1]]), 0],
					{k,Floor[depth/root[[2]]]}], {root,roots}]; 2 sum/denom]];
			If[mult > 0,
				weights[w] = mult;
				next = Join[next, Table[w - c[[All,i]], {i,r}]]]
		, {w,layer}];
		layer = DeleteDuplicates[next]; depth++];
	weights
];

(* Klimyk's shifted Weyl action: only the weights of the smaller factor
   are needed. Singular shifted weights contribute zero. This avoids building
   characters of much larger output representations (e.g. E8's 30380). *)
repProduct[t_, r_, a_List, b_List] /; supported[t,r] && validLabel[r,a] && validLabel[r,b] :=
 Module[{small=a, other=b, ch, c=cartan[t,r], rho=ConstantArray[1,r], out=<||>, w, sign, i, target},
	If[repDimension[t,r,b] < repDimension[t,r,a], {small,other}={b,a}];
	ch=repCharacter[t,r,small];
	KeyValueMap[Function[{weight,multiplicity},
		w=weight+other+rho; sign=1;
		While[Min[w]<0 && !MemberQ[w,0],
			i=First @ FirstPosition[w,_?Negative];
			w-=w[[i]] c[[All,i]]; sign=-sign];
		If[!MemberQ[w,0], target=w-rho;
			out[target]=If[KeyMemberQ[out,target],out[target],0]+sign multiplicity]
	],ch];
	If[AnyTrue[Values[out],#<0&], Message[repProduct::multiplicity,t,r,a,b]; Return[$Failed]];
	Flatten[KeyValueMap[ConstantArray[#1,#2] &,Select[out,#>0&]],1]
];

(* Retained as an independent algorithmic cross-check for moderate examples. *)
productByCharacters[t_, r_, a_List, b_List] /; supported[t,r] && validLabel[r,a] && validLabel[r,b] :=
 Module[{wa = repCharacter[t,r,a], wb = repCharacter[t,r,b], weights = <||>, w, n, inv = Inverse[cartan[t,r]], ch},
	Do[w = u+v; weights[w] = If[KeyMemberQ[weights,w], weights[w], 0] + wa[u] wb[v],
		{u,Keys[wa]}, {v,Keys[wb]}];
	CommonFunctions`MyReap[
		While[Length[weights] > 0,
			w = First @ MaximalBy[Keys[weights], Total[inv.#] &];
			n = weights[w];
			Do[Sow[w], {n}];
			ch = repCharacter[t,r,w];
			Do[weights[u] -= n ch[u], {u,Keys[ch]}];
			weights = Select[weights, # != 0 &]
		]
	]
];

(* The contravariant form obeys <F_i u,v> = <u,E_i v>.
   For an orthonormal height layer, the Gram matrix of its lowered vectors is
   <F_i a,F_j b> = <E_j a,E_i b> + delta_ij <a,H_i b>.
   Quotient its nullspace before proceeding, avoiding an exponential word basis. *)
repMatrices[t_, r_, a_List] /; supported[t,r] && validLabel[r,a] := Module[{result},
	result = construct[t,r,a];
	If[result === $Failed, $Failed, repMatrices[t,r,a] = result]
];
construct[t_, r_, a_] := Module[{c = cartan[t,r], expected = repDimension[t,r,a],
	weights = {a}, allWeights = {a}, steps = {}, oldLower = ConstantArray[{{0}},r], gram, pivots, transform,
	candidateWeights, reduced, block, rules, offsets = 0, total = 1, lower, width},
	rules = Table[{}, {r}];
	While[total < expected && weights =!= {},
		width = Length[weights];
		candidateWeights = Flatten[Table[w - c[[All,i]], {i,r}, {w,weights}],1];
		gram = arrayReduce @ ArrayFlatten @ Table[
			oldLower[[j]].Transpose[oldLower[[i]]] +
			If[i == j, DiagonalMatrix[weights[[All,i]]], ConstantArray[0,{width,width}]],
			{i,r}, {j,r}];
		reduced = Select[RowReduce[gram], AnyTrue[#, # != 0 &] &];
		pivots = (First @ FirstPosition[#, _?(# != 0 &)] &) /@ reduced;
		If[pivots === {}, weights = {}; Break[]];
		weights = candidateWeights[[pivots]];
		With[{form = gram[[pivots,pivots]]},
			transform = Orthogonalize[IdentityMatrix[Length[pivots]], (#1.form.#2) &]];
		AppendTo[steps,{pivots,transform}];
		Do[
			block = arrayReduce[transform.gram[[pivots,Range[(i-1) width+1,i width]]]];
			oldLower[[i]] = block;
			rules[[i]] = Join[rules[[i]],
				(#[[1]] + {total,offsets} -> #[[2]] & /@ Most[ArrayRules[SparseArray[block]]])]
		, {i,r}];
		offsets = total; total += Length[weights]; allWeights = Join[allWeights,weights]
	];
	If[total != expected, Message[repMatrices::dimension,t,r,a,total,expected]; Return[$Failed]];
	basisWeights[t,r,a] = allWeights;
	basisSteps[t,r,a] = steps;
	lower = SparseArray[#, {expected,expected}] & /@ rules;
	Join[lower, Transpose /@ lower]
];

(* Find highest-weight vectors in the tensor product, then follow precisely
   the basis changes used to construct the target representation. Only the
   highest-weight space needs a nullspace solve, not all tensor coefficients. *)
(* The dual metric is a transported highest-weight vector for the
   contragredient generators. No tensor-square nullspace is necessary. *)
repMetric[t_,rank_,a_] := repMetric[t,rank,a] = Module[{b=repDual[t,rank,a],ma,mb,d,pos,current,images,candidates,next,ops},
 {ma,mb}=repMatrices[t,rank,#]&/@{a,b};d=Length[ma[[1]]];
 pos=First@FirstPosition[basisWeights[t,rank,a],-b];
 current=ArrayReshape[SparseArray[{pos->1},d],{d,1}];images=current;
 ops=-Transpose[#]&/@Take[ma,rank];
 Do[candidates=Join[Sequence@@(#.current&/@ops),2];
  next=arrayReduce[candidates[[All,step[[1]]]].Transpose[step[[2]]]];
  images=Join[images,next,2];current=next,
 {step,basisSteps[t,rank,b]}];images];

(* Frobenius reciprocity keeps the largest representation on the output leg.
   Reusing an already constructed invariant tensor avoids solving a much larger
   highest-weight system after cyclic permutations of its three legs. *)
repIntertwiners[t_,rank_,a_List,b_List,c_List] /;
 supported[t,rank] && And@@(validLabel[rank,#]&/@{a,b,c}) :=
 repIntertwiners[t,rank,a,b,c] = Module[{da,db,dc,prime,tri,mb,mc},
 {da,db,dc}=repDimension[t,rank,#]&/@{a,b,c};
 If[Total[c]==0,Return[If[b===repDual[t,rank,a],{Flatten[Normal[repMetric[t,rank,a]]]},{}]]];
 If[da>db && da>dc,Return[Flatten[Transpose[ArrayReshape[#,{db,da,dc}],{2,1,3}]]&/@repIntertwiners[t,rank,b,a,c]]];
 If[db>dc,
  prime=repIntertwiners[t,rank,a,repDual[t,rank,c],repDual[t,rank,b]];
  mb=repMetric[t,rank,repDual[t,rank,b]];mc=Transpose[repMetric[t,rank,c]];
  Return[Table[
   tri=ArrayReshape[SparseArray[ArrayReshape[v,{da dc,db}]].mb,{da,dc,db}];
   Flatten[Normal[arrayReduce[ArrayReshape[Transpose[tri,{1,3,2}],{da db,dc}].mc]]],{v,prime}]]];
 constructIntertwiners[t,rank,a,b,c]
];

constructIntertwiners[t_, rank_, a_List, b_List, target_List] /;
 supported[t,rank] && And @@ (validLabel[rank,#] & /@ {a,b,target}) :=
 Module[{ma,mb,mc,da,db,dc,pairs,columns,operators,blocks,block,active,highest,
	maps,current,next,candidates,images},
	{ma,mb,mc} = repMatrices[t,rank,#] & /@ {a,b,target};
	If[MemberQ[{ma,mb,mc},$Failed], Return[$Failed]];
	{da,db,dc} = Length[#[[1]]] & /@ {ma,mb,mc};
	pairs = Flatten[Table[u+v,{u,basisWeights[t,rank,a]},{v,basisWeights[t,rank,b]}],1];
	columns = Flatten[Position[pairs,target,{1}]];
	If[columns === {}, Return[{}]];
	operators = Table[KroneckerProduct[ma[[i]],IdentityMatrix[db,SparseArray]] +
		KroneckerProduct[IdentityMatrix[da,SparseArray],mb[[i]]],{i,2 rank}];
	blocks = Table[
		block = operators[[rank+i]][[All,columns]];
		active = Union[First /@ (First /@ Most[ArrayRules[block]])];
		If[active === {}, Nothing, block[[active]]]
	, {i,rank}];
	block = If[blocks === {}, SparseArray[{}, {1,Length[columns]}], Join @@ blocks];
	highest = NullSpace[block];
	maps = Table[
		current = ArrayReshape[SparseArray[Thread[columns -> h],da db],{da db,1}];
		images = current;
		Do[
			candidates = Join[Sequence @@ Table[operators[[i]].current,{i,rank}],2];
			next = arrayReduce[candidates[[All,step[[1]]]].Transpose[step[[2]]]];
			images = Join[images,next,2]; current = next
		, {step,basisSteps[t,rank,target]}];
		Flatten[Normal[images]]
	, {h,highest}];
	maps
];
End[]
EndPackage[]
