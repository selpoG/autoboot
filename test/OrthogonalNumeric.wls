SetDirectory[DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]];
Get["ngroup.m"];
NGroupInfo`setPrecision[50];
CommonFunctions`importPackage["NGroupInfo`","OrthogonalChecks`",{"v","id","dim","prod","gA","gG","isrep","setGroup"}];
CommonFunctions`importPackage["NGroupInfoLie`","OrthogonalChecks`",{"getO"}];
NGroupInfo`setGroup[NGroupInfoLie`getO[4]];
CommonFunctions`importPackage["NClebschGordan`","OrthogonalChecks`",{"ope","setOps","op","bootAll","eqn","makeSDP","sdpobj"}];
Get["test/OrthogonalChecks.m"];

(* The shared suite must terminate explicitly; load/syntax failures must fail. *)
Exit[1];
