SetDirectory[DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]];
If[Length[$ScriptCommandLine] < 3, Print["Usage: Bootstrap.wls exact|numeric SO5|SO7|SO10|SO12|SO5Adj|SO5Sym2|A3Adj|Sp2|Sp3|Sp4|G2|F4 [mixed|single] [output.m]"]; Exit[2]];
mode = $ScriptCommandLine[[2]]; family = $ScriptCommandLine[[3]];
If[!MemberQ[{"exact","numeric"},mode] || !MemberQ[{"SO5","SO7","SO10","SO12","SO5Adj","SO5Sym2","A3Adj","Sp2","Sp3","Sp4","G2","F4"},family], Print["Unsupported mode or example"]; Exit[2]];
Get[If[mode == "numeric", "ngroup.m", "group.m"]];
prefix = If[mode == "numeric", "N", ""];
If[mode == "numeric", NGroupInfo`setPrecision[50]];
CommonFunctions`importPackage[prefix <> "GroupInfo`", "Bench`", {"id","dim","prod","v","setGroup"}];
CommonFunctions`importPackage[prefix <> "GroupInfoLie`", "Bench`", {"getSO","getSp","getLie"}];
Bench`family = family; Bench`mode = mode;
Bench`output = If[Length[$ScriptCommandLine] >= 5,$ScriptCommandLine[[5]],None];
Bench`mixed = Length[$ScriptCommandLine] >= 4 && $ScriptCommandLine[[4]] == "mixed";
Begin["Bench`"];
Print["START ", mode," ",family];
start = AbsoluteTime[];
{g,r} = Switch[family,
 "SO5",{getSO[5],v[1,0]}, "SO7",{getSO[7],v[1,0,0]},
 "SO10",{getSO[10],v[1,0,0,0,0]}, "SO12",{getSO[12],v[1,0,0,0,0,0]},
 "SO5Adj",{getSO[5],v[0,2]}, "SO5Sym2",{getSO[5],v[2,0]},
 "A3Adj",{getLie["A",3],v[1,0,1]},
 "Sp2",{getSp[2],v[1,0]}, "Sp3",{getSp[3],v[1,0,0]},
 "Sp4",{getSp[4],v[1,0,0,0]}, "G2",{getLie["G",2],v[1,0]},
 "F4",{getLie["F",4],v[0,0,0,1]}];
setGroup[g];
End[];
CommonFunctions`importPackage[prefix <> "ClebschGordan`", "Bench`", {"setOps","op","bootAll","eqn","makeSDP","sdpobj"}];
Begin["Bench`"];
Print["SETUP ",AbsoluteTime[]-start];
Print["SETOPS ", First[AbsoluteTiming[setOps[If[mixed,{op[s,g[id]],op[phi,r]},{op[phi,r]}]]]]];
{elapsed,result} = AbsoluteTiming[TimeConstrained[MemoryConstrained[bootAll[],1500000000,$Failed],120,$Aborted]];
Print["BOOT ",elapsed," HEAD ",Head[result]," LEAFCOUNT ",LeafCount[result]];
If[Head[result]===eqn,
 Print["EQUATIONS ",Length[result[[1]]]];
 Print["SDP ",First[AbsoluteTiming[sdp=makeSDP[result]]]," HEAD ",Head[sdp]];
 If[StringQ[output], Put[result,output]];
];
Print["TOTAL ",AbsoluteTime[]-start];
Exit[If[Head[result]===eqn && Head[sdp]===sdpobj,0,1]];
End[];
