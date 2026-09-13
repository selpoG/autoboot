repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
SetDirectory[repo];
Get["ngroup.m"];
NGroupInfo`setPrecision[50];
CommonFunctions`importPackage["NGroupInfo`", "BootstrapChecks`", {"id","dim","prod","v","gA","setGroup"}];
CommonFunctions`importPackage["NGroupInfoLie`", "BootstrapChecks`", {"getSO","getSp","getLie"}];
(* Load CG before importing its public symbols; each setGroup clears caches. *)
NGroupInfo`setGroup[NGroupInfoLie`getSO[5]];
CommonFunctions`importPackage["NClebschGordan`", "BootstrapChecks`", {"op","setOps","bootAll","eqn","makeSDP","sdpobj","ope","cor","six","sum","single","F","H","Fp","Hp","\[Beta]"}];
Get["test/BootstrapChecks.m"];

(* The shared suite must terminate explicitly; load/syntax failures must fail. *)
Exit[1];
