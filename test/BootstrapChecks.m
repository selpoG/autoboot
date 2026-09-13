(* End-to-end equations, independent tensor identities, and saved exact systems. *)
Begin["BootstrapChecks`"];
failed = False;
assert[name_, condition_] := If[TrueQ[condition], Print[name,"=PASS"], Print[name,"=FAIL"]; failed=True];
nearZero[x_] := Max[Append[Abs[Flatten[N[RootReduce[Normal[x]],40]]],0]] < 10^-25;
normalize[x_] := FixedPoint[Expand[# /. {
 sum[y_Plus,o_] :> Total[sum[#,o] & /@ List @@ y], single[y_Plus] :> Total[single /@ List @@ y],
 sum[c_?NumericQ y_,o_] :> c sum[y,o], single[c_?NumericQ y_] :> c single[y]}] &, x];
equivalent[a_,b_] := Module[{pa=normalize[a[[1]]],pb=normalize[b[[1]]],atoms,ma,mb},
 atoms = Union[Cases[{pa,pb},_sum|_single,Infinity]];
 ma = RowReduce[Table[Coefficient[p,t],{p,pa},{t,atoms}]];
 mb = RowReduce[Table[Coefficient[p,t],{p,pb},{t,atoms}]];
 Dimensions[ma] === Dimensions[mb] && nearZero[ma-mb]
];

Do[
 {name,g,r,mixed,expectedCount} = sample;
 setGroup[g];
 setOps[If[mixed,{op[s,g[id]],op[phi,r]},{op[phi,r]}]];
 {seconds,result} = AbsoluteTiming[bootAll[]];
 assert[name<>" bootstrap", Head[result]===eqn && Length[result[[1]]]===expectedCount];
 assert[name<>" SDP", Head[makeSDP[result]]===sdpobj];
 If[name=="SO5",assert["SO5 reference equations",equivalent[result,Get["test/fixtures/BootstrapSO5.m"]]]];
 If[name=="Sp2",assert["Sp2 reference equations",equivalent[result,Get["test/fixtures/BootstrapSp2.m"]]]];
 If[name=="G2 mixed",assert["G2 exact/numeric equations",equivalent[result,Get["test/fixtures/BootstrapG2Mixed.m"]]]];
 Print[name," BOOT_SECONDS=",seconds];
 d = g[dim[r]]; channels = g[prod[r,r]];
 matrices = Table[ArrayReshape[Table[ope[r,r,t][1][a,b,c],{a,d},{b,d},{c,g[dim[t]]}],{d^2,g[dim[t]]}],{t,channels}];
 assert[name<>" CG isometries",And @@ MapThread[nearZero[ConjugateTranspose[#1].#1-IdentityMatrix[g[dim[#2]]]] &,{matrices,channels}]];
 projectors = (#.ConjugateTranspose[#] &) /@ matrices;
 assert[name<>" CG completeness",nearZero[Total[projectors]-IdentityMatrix[d^2]]];
 assert[name<>" intertwining equations",And @@ Flatten[Table[
  With[{x=g[gA][[i]][r],y=g[gA][[i]][channels[[k]]],c=matrices[[k]]},
   nearZero[(KroneckerProduct[x,IdentityMatrix[d]]+KroneckerProduct[IdentityMatrix[d],x]).c-c.y]],
  {i,Length[g[gA]]},{k,Length[channels]}]]];
 (* All index components, not the interpolation samples used by setSix. *)
 tensors = Table[Table[cor[r,r,r,r][t,1,1][a,b,c,e],{a,d},{b,d},{c,d},{e,d}],{t,channels}];
 basis = Flatten /@ tensors;
 crossed = Flatten[Transpose[#,{1,4,3,2}]] & /@ tensors;
 crossing = Table[six[r,r,r,r][t,1,1,u,1,1],{t,channels},{u,channels}];
 assert[name<>" tensor normalization",nearZero[Conjugate[basis].Transpose[basis]-IdentityMatrix[Length[channels]]]];
 assert[name<>" crossing reconstruction",nearZero[Conjugate[basis].Transpose[crossed]-crossing]];
 assert[name<>" crossing involution",nearZero[crossing.crossing-IdentityMatrix[Length[channels]]]];
 If[name=="SO5",
  metric = Flatten[Table[ope[r][a,b],{a,d},{b,d}]];
  singletProjector = Outer[Times,metric,Conjugate[metric]]/d;
  swap = SparseArray[Flatten[Table[{d(a-1)+b,d(b-1)+a}->1,{a,d},{b,d}],1],{d^2,d^2}];
  analytic = Table[Switch[g[dim[t]],1,singletProjector,10,(IdentityMatrix[d^2]-swap)/2,
   14,(IdentityMatrix[d^2]+swap)/2-singletProjector],{t,channels}];
  assert["SO5 analytic projectors",And @@ MapThread[nearZero[#1-#2]&,{projectors,analytic}]];
 ];
,{sample,{{"SO5",getSO[5],v[1,0],False,3},{"Sp2",getSp[2],v[1,0],False,21},
 {"G2 mixed",getLie["G",2],v[1,0],True,9}}}];
(* The sparse tensor cache must retain multiplicity axes and complex duals. *)
g = getLie["A",2]; setGroup[g]; r = v[1,1];
cgCopies = Table[ArrayReshape[Table[ope[r,r,r][n][a,b,c],{a,8},{b,8},{c,8}],{64,8}],{n,2}];
assert["CG multiplicity orthogonality",And @@ Flatten[Table[
 nearZero[ConjugateTranspose[cgCopies[[i]]].cgCopies[[j]]-If[i==j,IdentityMatrix[8],ConstantArray[0,{8,8}]]],{i,2},{j,2}]]];
triples = Table[Flatten[Table[cor[r,r,r][n][a,b,c],{a,8},{b,8},{c,8}]],{n,2}];
assert["three-point multiplicity normalization",nearZero[Conjugate[triples].Transpose[triples]-IdentityMatrix[2]]];
r = v[1,0]; rb = v[0,1]; channels = g[prod[r,rb]];
tensors = Table[Table[cor[r,rb,r,rb][t,1,1][a,b,c,e],{a,3},{b,3},{c,3},{e,3}],{t,channels}];
basis = Flatten /@ tensors; crossed = Flatten[Transpose[#,{1,4,3,2}]] & /@ tensors;
crossing = Table[six[r,rb,r,rb][t,1,1,u,1,1],{t,channels},{u,channels}];
assert["complex dual crossing",nearZero[Conjugate[basis].Transpose[crossed]-crossing] && nearZero[crossing.crossing-IdentityMatrix[2]]];

(* Frobenius reciprocity must also work for complex and pseudoreal legs.
   Compare with the direct highest-weight solve, bypassing reciprocity. *)
Do[
 {type,rank,aa,bb,cc}=sample;g=getLie[type,rank];setGroup[g];{a,b,c}=(v@@#&/@{aa,bb,cc});
 {da,db,dc}=g[dim[#]]&/@{a,b,c};
 cm=ArrayReshape[Table[ope[a,b,c][1][i,j,k],{i,da},{j,db},{k,dc}],{da db,dc}];
 direct=First[LieRepresentations`Private`constructIntertwiners[type,rank,aa,bb,cc]];
 direct=Sqrt[dc] direct/Sqrt[Conjugate[direct].direct];direct=ArrayReshape[direct,{da db,dc}];
 assert[type<>" reciprocity isometry",nearZero[ConjugateTranspose[cm].cm-IdentityMatrix[dc]]];
 assert[type<>" reciprocity direct-solve projector",nearZero[cm.ConjugateTranspose[cm]-direct.ConjugateTranspose[direct]]];
,{sample,{{"A",2,{1,0},{1,1},{1,0}},{"C",2,{1,0},{2,0},{1,0}},{"G",2,{1,0},{0,1},{1,0}}}}];

Print["BOOTSTRAP_REGRESSION=",If[failed,"FAIL","PASS"]];
Exit[If[failed,1,0]];
End[];
