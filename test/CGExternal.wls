(* Independent component oracle: Wolfram's built-in Wigner ClebschGordan,
   not autoboot's SU(2) backend. Integer spins also descend to SO(3). *)
SetDirectory[DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]];
Get["group.m"];
CommonFunctions`importPackage["GroupInfo`","CGExternal`",{"v","setGroup"}];
CommonFunctions`importPackage["GroupInfoLie`","CGExternal`",{"getLie"}];
CGExternal`setGroup[CGExternal`getLie["A",1]];
CommonFunctions`importPackage["ClebschGordan`","CGExternal`",{"ope"}];
Begin["CGExternal`"];
failed=False;
j1=10;j2=8;
Do[
 {seconds,cg}=AbsoluteTiming[Flatten[Table[ope[v[2 j1],v[2 j2],v[2 j3]][1][a,b,c],
 {a,2 j1+1},{b,2 j2+1},{c,2 j3+1}]]];
 reference=Flatten[Table[If[m1+m2==m3,
 System`ClebschGordan[{j1,m1},{j2,m2},{j3,m3}],0],
 {m1,j1,-j1,-1},{m2,j2,-j2,-1},{m3,j3,-j3,-1}]];
 (* Global phase is conventional; no component-dependent adjustment is allowed. *)
 pivot=First@FirstPosition[reference,_?(#!=0&)];phase=RootReduce[cg[[pivot]]/reference[[pivot]]];
 ok=TrueQ[Abs[phase]==1] && And@@(TrueQ[#==0]&/@RootReduce[cg-phase reference]);
 Print["SO3 spins ",{j1,j2,j3}," BUILTIN_CG_COMPONENTS=",Length[cg]," SECONDS=",seconds," ",If[ok,"PASS","FAIL"]];
 If[!ok,failed=True];
,{j3,{2,9,18}}];
Exit[If[failed,1,0]];
End[];
