(* Run from the repository root with low- and high-order Majorana table directories.
   Optional fourth argument: SDPB 3.1.0 mathematica/SDPB.m. *)
Get["spinning-identical.m"];
If[!MemberQ[{3,4},Length[$ScriptCommandLine]],Print["Usage: ValidateSpinning.wls LOW_DIR HIGH_DIR [SDPB.m]"];Exit[2]];
failed=False;
assert[label_,ok_]:=If[TrueQ[ok],Print[label,"=PASS"],Print[label,"=FAIL"];failed=True];
models={};residuals={};
Do[
 files=Import[FileNameJoin[{dir,"manifest.json"}],"RawJSON"];
 {seconds,tables}=AbsoluteTiming[readSpinningBlocks[FileNameJoin[{dir,#}],75]&/@files];
 assert["all tables read",AllTrue[tables,AssociationQ]];If[failed,Exit[1]];
 Print["TABLE_IMPORT_SECONDS=",seconds];
 spins=Union[Lookup[tables,"j_internal"]];
 sectors=Flatten[Table[<|"Spin"->j,"Parity"->p,"Gap"->If[j==0,1,j+11/10]|>,{j,spins},{p,{-1,1}}]];
 {seconds,model}=AbsoluteTiming[spinIdenticalSDP[tables,sectors,x,3]];
 assert["Majorana assembly",AssociationQ[model]];If[failed,Exit[1]];
 Print["ASSEMBLY_SECONDS=",seconds];
 reduced=spinMajoranaRows[model];assert["13 independent components",Length[reduced["Identity"]]==13];
 assert["two OPE structures",Dimensions[reduced["Sectors"]["{2, 1}"]["Matrices"]]=={13,2,2}];
 (* Independently project every original equation onto the published selection,
    using all polynomial coefficients and the identity as columns. *)
 columns={model["Identity"]};
 Do[mat=block["Matrices"];degree=Max[0,Sequence@@Exponent[Flatten[mat],x]];
  Do[AppendTo[columns,Coefficient[mat[[All,a,b]],x,k]],{a,Length[mat[[1]]]},{b,Length[mat[[1]]]},{k,0,degree}],{block,Values[model["Sectors"]]}];
 vectors=Transpose[columns];basis=Orthogonalize[vectors[[reduced["RetainedRows"]]],Method->"GramSchmidt"];
 residual=Max[Table[If[Sign[Norm[row]]===0,0,Norm[row-Total[(#.row) #&/@basis]]/Norm[row]],{row,vectors}]];
 AppendTo[residuals,residual];Print["OMITTED_COMPONENT_RESIDUAL=",residual];
 file=FileNameJoin[{$TemporaryDirectory,"autoboot-pmp-"<>CreateUUID[]<>".json"}];
 assert["real blocks PMP export",spinWritePMP[file,reduced,x,40]===file];
 If[FileExistsQ[file],DeleteFile[file]];AppendTo[models,reduced];
,{dir,$ScriptCommandLine[[2;;3]]}];
assert["omitted-component convergence",Last[residuals]<10^-12 && Last[residuals]<First[residuals]/10^4];
If[Length[$ScriptCommandLine]==4,
 Get[$ScriptCommandLine[[4]]];sdp=Last[models];
 official=SDP[ConstantArray[0,Length[sdp["Identity"]]],sdp["Identity"],Map[Function[b,PositiveMatrixWithPrefactor[DampedRational[b["Prefactor"][[1]],b["Prefactor"][[2]],b["Prefactor"][[3]],x],Transpose[b["Matrices"],{3,1,2}]]],Values[sdp["Sectors"]]]];
 f1=FileNameJoin[{$TemporaryDirectory,"autoboot-pmp-"<>CreateUUID[]<>".json"}];f2=FileNameJoin[{$TemporaryDirectory,"official-pmp-"<>CreateUUID[]<>".json"}];
 spinWritePMP[f1,sdp,x,40];WritePmpJson[f2,official,40];
 (* Decimal parser from the block reader accepts no executable Wolfram input. *)
 normalize[value_String]:=SpinningBlocks`Private`numberString[value];
 normalize[value_List]:=normalize/@value;
 normalize[value_Association]:=Map[normalize,KeySort[value]];
 own=normalize[Import[f1,"RawJSON"]];other=normalize[Import[f2,"RawJSON"]];
 same[a_Association,b_Association]:=Keys[a]===Keys[b] && AllTrue[Keys[a],same[a[#],b[#]]&];
 same[a_List,b_List]:=Length[a]==Length[b] && And@@MapThread[same,{a,b}];
 same[a_?NumericQ,b_?NumericQ]:=TrueQ[Abs[a-b]<=10^-37 Max[1,Abs[a],Abs[b]]];
 assert["official SDPB writer",same[own,other]];
 DeleteFile/@{f1,f2}];
Print["SPINNING_INTEGRATION=",If[failed,"FAIL","PASS"]];Exit[If[failed,1,0]];
