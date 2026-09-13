(* Separate process per case: includes cold caches and bounds expensive tests. *)
SetDirectory[DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]];
If[Length[$ScriptCommandLine]<3,Print["Usage: ValidateLarge.wls exact|numeric CASE [equations.m]"];Exit[2]];
mode=$ScriptCommandLine[[2]]; name=$ScriptCommandLine[[3]];
If[!MemberQ[{"exact","numeric"},mode],Exit[2]];
prefix=If[mode=="numeric","N",""];Get[If[mode=="numeric","ngroup.m","group.m"]];
If[mode=="numeric",NGroupInfo`setPrecision[50]];
CommonFunctions`importPackage[prefix<>"GroupInfo`","Large`",{"id","dim","prod","v","gA","setGroup"}];
CommonFunctions`importPackage[prefix<>"GroupInfoLie`","Large`",{"getSO","getSp","getLie"}];
Large`name=name; Large`mode=mode;
Large`output=If[Length[$ScriptCommandLine]>=4,$ScriptCommandLine[[4]],None];
Begin["Large`"];
{g,r,analyticN}=Switch[name,
 "SO7",{getSO[7],v[1,0,0],7},
 "SO10",{getSO[10],v[1,0,0,0,0],10},
 "SO12",{getSO[12],v[1,0,0,0,0,0],12},
 "SO5Adj",{getSO[5],v[0,2],0},
 "SO5Sym2",{getSO[5],v[2,0],0},
 "A3Adj",{getLie["A",3],v[1,0,1],0},
 "Sp4",{getSp[4],v[1,0,0,0],0},
 "F4",{getLie["F",4],v[0,0,0,1],0},
 _,Print["Unknown case"];Exit[2]];
setGroup[g];End[];
CommonFunctions`importPackage[prefix<>"ClebschGordan`","Large`",{"op","setOps","bootAll","eqn","makeSDP","sdpobj","ope","cor","six","sum","single"}];
Begin["Large`"];
failed=False;
assert[label_,p_]:=If[TrueQ[p],Print[label,"=PASS"],Print[label,"=FAIL"];failed=True];
nearZero[x_]:=Max[Append[Abs[Flatten[Normal[N[x,40]]]],0]]<10^-25;
  normalize[x_]:=FixedPoint[Expand[#/.{
   sum[y_Plus,o_]:>Total[sum[#,o]&/@List@@y],single[y_Plus]:>Total[single/@List@@y],
   sum[c_?NumericQ y_,o_]:>c sum[y,o],single[c_?NumericQ y_]:>c single[y]}]&,x];
start=AbsoluteTime[];
Print["START ",mode," ",name];
setOps[{op[phi,r]}];
{seconds,result}=AbsoluteTiming[bootAll[]];
assert["bootstrap",Head[result]===eqn];
If[Head[result]=!=eqn,Exit[1]];
assert["SDP",Head[makeSDP[result]]===sdpobj];
Print["BOOT_SECONDS=",seconds," EQUATIONS=",Length[result[[1]]]];
If[name=="F4",
 fixture=Get["test/fixtures/BootstrapF4.m"];
 pa=normalize[result[[1]]];pb=normalize[fixture[[1]]];atoms=Union[Cases[{pa,pb},_sum|_single,Infinity]];
 ma=RowReduce[Table[Coefficient[p,t],{p,pa},{t,atoms}]];mb=RowReduce[Table[Coefficient[p,t],{p,pb},{t,atoms}]];
 assert["F4 saved equation regression",Dimensions[ma]===Dimensions[mb] && nearZero[ma-mb]];
];
If[StringQ[output],
 If[mode=="exact",Put[result,output],
  If[!FileExistsQ[output],Print["Missing exact reference: ",output];Exit[2]];
  reference=Get[output];pa=normalize[result[[1]]];pb=normalize[reference[[1]]];
  atoms=Union[Cases[{pa,pb},_sum|_single,Infinity]];
  ma=RowReduce[Table[Coefficient[p,t],{p,pa},{t,atoms}]];
  mb=RowReduce[Table[Coefficient[p,t],{p,pb},{t,atoms}]];
  assert["exact/numeric equation row space",Dimensions[ma]===Dimensions[mb] && nearZero[ma-mb]];
 ]];
d=g[dim[r]];channels=Tally[g[prod[r,r]]];
indices=Flatten[Table[{t[[1]],n},{t,channels},{n,t[[2]]}],1];
cgs=Table[SparseArray[N[ArrayReshape[Table[ope[r,r,t[[1]]][t[[2]]][a,b,c],{a,d},{b,d},{c,g[dim[t[[1]]]]}],{d^2,g[dim[t[[1]]]]}],40]],{t,indices}];
(* Each multiplicity is a separate embedding; joining gives a unitary change of basis. *)
whole=Join[Sequence@@cgs,2];
assert["CG orthonormality and completeness",Dimensions[whole]=={d^2,d^2} && nearZero[ConjugateTranspose[whole].whole-IdentityMatrix[d^2]]];
assert["all CG intertwining equations",And@@Flatten[Table[
 With[{x=N[g[gA][[i]][r],40],y=N[g[gA][[i]][indices[[k,1]]],40],c=cgs[[k]]},
 nearZero[(KroneckerProduct[x,IdentityMatrix[d,SparseArray]]+KroneckerProduct[IdentityMatrix[d,SparseArray],x]).c-c.y]],
 {i,Length[g[gA]]},{k,Length[indices]}]]];
If[analyticN>0,
 metric=Flatten[Table[ope[r][a,b],{a,d},{b,d}]];
 p0=N[Outer[Times,metric,Conjugate[metric]]/d,40];
 swap=SparseArray[Flatten[Table[{d(a-1)+b,d(b-1)+a}->1,{a,d},{b,d}],1],{d^2,d^2}];
 projectors=Table[Switch[g[dim[t[[1]]]],1,p0,d(d-1)/2,(IdentityMatrix[d^2,SparseArray]-swap)/2,
 d(d+1)/2-1,(IdentityMatrix[d^2,SparseArray]+swap)/2-p0],{t,indices}];
 assert["SO(N) analytic CG projectors",And@@MapThread[nearZero[#1.ConjugateTranspose[#1]-#2]&,{cgs,projectors}]];
];
(* Reconstruct every four-index component, including multiplicity pairs. *)
quadIndices=Flatten[Table[{t[[1]],n,m},{t,channels},{n,t[[2]]},{m,t[[2]]}],2];
tensors=Table[N[RootReduce[Table[cor[r,r,r,r][t[[1]],t[[2]],t[[3]]][a,b,c,e],{a,d},{b,d},{c,d},{e,d}]],40],{t,quadIndices}];
basis=Flatten/@tensors;crossed=Flatten[Transpose[#,{1,4,3,2}]]&/@tensors;
crossing=N[Table[six[r,r,r,r][Sequence@@Join[t,u]],{t,quadIndices},{u,quadIndices}],40];
assert["four-tensor orthonormality",nearZero[Conjugate[basis].Transpose[basis]-IdentityMatrix[Length[basis]]]];
assert["all-component crossing reconstruction",nearZero[Conjugate[basis].Transpose[crossed]-crossing]];
assert["crossing involution",nearZero[crossing.crossing-IdentityMatrix[Length[basis]]]];
Print["TOTAL_SECONDS=",AbsoluteTime[]-start," RESULT=",If[failed,"FAIL","PASS"]];
Exit[If[failed,1,0]];
