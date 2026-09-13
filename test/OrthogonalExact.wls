SetDirectory[DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]];
Get["group.m"];

CommonFunctions`importPackage["GroupInfo`","OrthogonalChecks`",{"v","id","dim","prod","gA","gG","isrep","setGroup"}];
CommonFunctions`importPackage["GroupInfoLie`","OrthogonalChecks`",{"getO"}];
GroupInfo`setGroup[GroupInfoLie`getO[4]];
CommonFunctions`importPackage["ClebschGordan`","OrthogonalChecks`",{"ope","setOps","op","bootAll","eqn","makeSDP","sdpobj"}];
Get["test/OrthogonalChecks.m"];

(* The shared suite must terminate explicitly; load/syntax failures must fail. *)
Exit[1];
