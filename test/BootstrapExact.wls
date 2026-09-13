repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
SetDirectory[repo];
Get["group.m"];

CommonFunctions`importPackage["GroupInfo`", "BootstrapChecks`", {"id","dim","prod","v","gA","setGroup"}];
CommonFunctions`importPackage["GroupInfoLie`", "BootstrapChecks`", {"getSO","getSp","getLie"}];
(* Load CG before importing its public symbols; each setGroup clears caches. *)
GroupInfo`setGroup[GroupInfoLie`getSO[5]];
CommonFunctions`importPackage["ClebschGordan`", "BootstrapChecks`", {"op","setOps","bootAll","eqn","makeSDP","sdpobj","ope","cor","six","sum","single","F","H","Fp","Hp","\[Beta]"}];
Get["test/BootstrapChecks.m"];

(* The shared suite must terminate explicitly; load/syntax failures must fail. *)
Exit[1];
