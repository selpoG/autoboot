Begin["OrthogonalChecks`"];
failed=False;
assert[label_,p_]:=If[TrueQ[p],Print[label,"=PASS"],Print[label,"=FAIL"];failed=True];
nearZeroO[x_]:=Max[Append[Abs[Flatten[N[Normal[RootReduce[x]],40]]],0]]<10^-25;
Do[
 {n,r}=sample;g=getO[n];rank=Floor[n/2];d=g[dim[r]];
 Print["O",n," vector dimension=",d];assert["vector dimension",d==n];
 det=v@@Append[ConstantArray[0,rank],-1];
 assert["determinant square",g[prod[det,det]]=={g[id]}];
 assert["determinant twist",g[prod[r,det]]=={ReplacePart[r,-1->-Last[r]]}];
 rr=g[gG][[1]][r];ms=#[r]&/@g[gA];
 perm=If[OddQ[n],Range[rank],Join[Range[rank-2],{rank,rank-1}]];
 assert["reflection involution",nearZeroO[rr.rr-IdentityMatrix[d]]];
 assert["reflection conjugation",And@@Table[nearZeroO[rr.ms[[i]].rr-ms[[Join[perm,rank+perm][[i]]]]],{i,2 rank}]];
 channels=g[prod[r,r]];assert["vector tensor dimensions",Sort[g[dim[#]]&/@channels]==Sort[{1,n(n-1)/2,n(n+1)/2-1}]];
 setGroup[g];
 cs=Table[ArrayReshape[Table[ope[r,r,t][1][a,b,c],{a,d},{b,d},{c,g[dim[t]]}],{d^2,g[dim[t]]}],{t,channels}];
 all=Join[Sequence@@cs,2];assert["CG completeness",nearZeroO[ConjugateTranspose[all].all-IdentityMatrix[d^2]]];
 assert["CG reflection",And@@MapThread[nearZeroO[KroneckerProduct[rr,rr].#1-#1.g[gG][[1]][#2]]&,{cs,channels}]];
 assert["CG Lie equations",And@@Flatten[Table[nearZeroO[(KroneckerProduct[x[r],IdentityMatrix[d]]+KroneckerProduct[IdentityMatrix[d],x[r]]).cs[[k]]-cs[[k]].x[channels[[k]]]],{x,g[gA]},{k,Length[channels]}]]];
 setOps[{op[phi,r]}];{secs,result}=AbsoluteTiming[bootAll[]];
 assert["vector bootstrap",Head[result]===eqn && Length[result[[1]]]==3];
 assert["vector SDP",Head[makeSDP[result]]===sdpobj];Print["O",n," BOOT_SECONDS=",secs];
,{sample,{{4,v[1,1,1]},{5,v[1,0,-1]},{6,v[1,0,0,1]},{7,v[1,0,0,-1]},{8,v[1,0,0,0,1]}}}];
(* A non-fixed D3 pair: the two chiral three-forms form one real O(6) irrep. *)
g=getO[6];r=v[0,2,0,0];det=v[0,0,0,-1];
assert["induced dimension",g[dim[r]]==20];
assert["induced determinant twist",g[prod[r,det]]=={r}];
assert["reject duplicate orbit",!g[isrep[v[0,0,2,0]]]];
assert["reject parity on pair",!g[isrep[v[0,2,0,1]]]];
assert["reject spinor",!g[isrep[v[0,1,0,0]]]];
setGroup[g];j=Table[ope[r][a,b],{a,20},{b,20}];
assert["induced real metric",nearZeroO[j-Transpose[j]] && nearZeroO[j.ConjugateTranspose[j]-IdentityMatrix[20]]];
rr=g[gG][[1]][r];assert["induced metric reflection",nearZeroO[rr.j.Transpose[rr]-j]];
(* Reflection separates the two SO(6) adjoint copies into opposite extensions. *)
g=getO[6];r=v[0,1,1,1];rp=v[0,1,1,-1];channels=g[prod[r,r]];
assert["adjoint reflection multiplicities",Count[channels,r]==1 && Count[channels,rp]==1 && Total[g[dim[#]]&/@channels]==225];
setGroup[g];rr=g[gG][[1]][r];
cs=Table[ArrayReshape[Table[ope[r,r,t][1][a,b,c],{a,15},{b,15},{c,15}],{225,15}],{t,{r,rp}}];
assert["adjoint parity CG orthogonality",nearZeroO[ConjugateTranspose[cs[[1]]].cs[[2]]]];
assert["adjoint parity CG reflection",And@@MapThread[nearZeroO[KroneckerProduct[rr,rr].#1-#1.g[gG][[1]][#2]]&,{cs,{r,rp}}]];
(* Induced representations must also work on input legs, not only as outputs. *)
Do[
 {n,a,b}=sample;g=getO[n];setGroup[g];da=g[dim[a]];db=g[dim[b]];
 channels=Tally[g[prod[a,b]]];
 labels=Flatten[Table[{t[[1]],k},{t,channels},{k,t[[2]]}],1];
 cs=Table[ArrayReshape[Table[ope[a,b,t[[1]]][t[[2]]][i,j,k],
  {i,da},{j,db},{k,g[dim[t[[1]]]]}],{da db,g[dim[t[[1]]]]}],{t,labels}];
 all=Join[Sequence@@cs,2];
 assert["induced input completeness",Dimensions[all]=={da db,da db} && nearZeroO[ConjugateTranspose[all].all-IdentityMatrix[da db]]];
 assert["induced input reflection",And@@MapThread[
  nearZeroO[KroneckerProduct[g[gG][[1]][a],g[gG][[1]][b]].#1-#1.g[gG][[1]][#2[[1]]]]&,{cs,labels}]];
 assert["induced input Lie equations",And@@Flatten[Table[
  nearZeroO[(KroneckerProduct[x[a],IdentityMatrix[db]]+KroneckerProduct[IdentityMatrix[da],x[b]]).cs[[k]]-cs[[k]].x[labels[[k,1]]]],
  {x,g[gA]},{k,Length[labels]}]]];
,{sample,{{4,v[2,0,0],v[1,1,1]},{6,v[0,2,0,0],v[1,0,0,1]}}}];
(* Existing low-degree APIs are retained. *)
assert["O2 compatibility",getO[2][dim[v[2]]]==2];
assert["O3 compatibility",getO[3][dim[v[2,-1]]]==5];
Do[
 {n,r}=sample;g=getO[n];setGroup[g];setOps[{op[phi,r]}];result=bootAll[];
 assert["O"<>ToString[n]<>" bootstrap compatibility",Head[result]===eqn && Length[result[[1]]]==3 && Head[makeSDP[result]]===sdpobj];
,{sample,{{2,v[1]},{3,v[1,-1]}}}];
Print["ORTHOGONAL_REGRESSION=",If[failed,"FAIL","PASS"]];Exit[If[failed,1,0]];
End[];
