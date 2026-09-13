(* Loaded inside each Lie package's Private context, after importing group
   symbols. lieCheck and lieNumber select exact or configured numerical mode. *)
getLie[t_, r_] /; LieRepresentations`supported[t,r] := lieCheck[
	getLie[t,r] = makeLieGroup[lie[t,r], t, r, (True &)]
];
getLie[args___] := (Message[getLie::type, HoldForm[{args}]]; $Failed);

getSp[n_Integer] /; n >= 1 := lieCheck[
	getSp[n] = makeLieGroup[sp[n], "C", n, (True &)]
];
getSp[args___] := (Message[getSp::degree, HoldForm[{args}]]; $Failed);

getSpin[n_Integer] /; n >= 4 := lieCheck[
	getSpin[n] = makeLieGroup[spin[n], If[OddQ[n], "B", "D"], Floor[n/2], (True &)]
];
getSpin[args___] := (Message[getSpin::degree, HoldForm[{args}]]; $Failed);

getSO[n_Integer] /; n >= 4 := lieCheck[
	getSO[n] = makeLieGroup[so[n], If[OddQ[n], "B", "D"], Floor[n/2],
		If[OddQ[n], (EvenQ[Last[#]] &), (EvenQ[Total[Take[#,-2]]] &)]]
];
getSO[args___] := (Message[getSO::degree, HoldForm[{args}]]; $Failed);

makeLieGroup[group_, type_, rank_, allowed_] := AbortProtect @ Module[{G = group, gen},
	G[id] = v @@ ConstantArray[0,rank];
	G[isrep[_]] := False;
	G[isrep[r_v]] := LieRepresentations`validLabel[rank,List @@ r] && allowed[List @@ r];
	G[dim[r_v]] /; G[isrep[r]] := LieRepresentations`repDimension[type,rank,List @@ r];
	G[dual[r_v]] /; G[isrep[r]] := v @@ LieRepresentations`repDual[type,rank,List @@ r];
	G[minrep[a_v,b_v]] /; G[isrep[a]] && G[isrep[b]] :=
		First @ SortBy[{a,b}, {Total[List @@ #] &, (List @@ #) &}];
	G[prod[a_v,b_v]] /; G[isrep[a]] && G[isrep[b]] :=
		G[prod[a,b]] = G[prod[b,a]] = (v @@ # & /@ LieRepresentations`repProduct[type,rank,List @@ a,List @@ b]);
	G[LieRepresentations`intertwiners[a_v,b_v,c_v]] /; G[isrep[a]] && G[isrep[b]] && G[isrep[c]] :=
		LieRepresentations`repIntertwiners[type,rank,List @@ a,List @@ b,List @@ c];
	G[gG] = {};
	G[gA] = Array[G[gen[#]] &, 2 rank];
	G[gen[j_Integer]][r_v] /; 1 <= j <= 2 rank && G[isrep[r]] := Module[{m},
		m = LieRepresentations`repMatrices[type,rank,List @@ r];
		If[m === $Failed, $Failed, G[gen[j]][r] = lieNumber[m[[j]]]]
	];
	G
];

getO[n_Integer] /; n >= 4 := lieCheck[getO[n] = Module[{G=o[n],gen,reflect,rank=Floor[n/2]},
 G[id] = v @@ Append[ConstantArray[0,rank],1];
 G[isrep[_]] := False;
 G[isrep[q_v]] := OrthogonalRepresentations`oValid[n,List@@q];
 G[dim[q_v]] /; G[isrep[q]] := OrthogonalRepresentations`oDimension[n,List@@q];
 G[dual[q_v]] /; G[isrep[q]] := q;
 G[minrep[a_v,b_v]] /; G[isrep[a]] && G[isrep[b]] := First@Sort[{a,b}];
 G[prod[a_v,b_v]] /; G[isrep[a]] && G[isrep[b]] := G[prod[a,b]]=G[prod[b,a]]=
  (v@@#&/@OrthogonalRepresentations`oProduct[n,List@@a,List@@b]);
 G[gA] = Array[G[gen[#]]&,2 rank]; G[gG] = {G[reflect]};
 G[gen[j_Integer]][q_v] /; 1<=j<=2 rank && G[isrep[q]] := G[gen[j]][q]=lieNumber[OrthogonalRepresentations`oMatrices[n,List@@q][[j]]];
 G[reflect][q_v] /; G[isrep[q]] := G[reflect][q]=lieNumber[OrthogonalRepresentations`oReflection[n,List@@q]];
 G[LieRepresentations`intertwiners[a_v,b_v,c_v]] /; G[isrep[a]] && G[isrep[b]] && G[isrep[c]] :=
  OrthogonalRepresentations`oIntertwiners[n,List@@a,List@@b,List@@c];
 G]];
getO[args___] := (Message[getO::degree,HoldForm[{args}]];$Failed);
