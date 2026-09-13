repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
SetDirectory[repo];
Get["ngroup.m"];
NGroupInfo`setPrecision[50];
CommonFunctions`importPackage["NGroupInfo`", "LieChecks`", {"id", "dim", "dual", "prod", "v", "isrep", "minrep", "gA", "setGroup"}];
CommonFunctions`importPackage["NGroupInfoLie`", "LieChecks`", {"getLie", "getSO", "getSpin", "getSp"}];
LieChecks`eq = NClebschGordan`Private`eq;
LieChecks`ope = NClebschGordan`ope;
LieChecks`isPseaudo = NClebschGordan`isPseaudo;
Get["test/LieChecks.m"];

(* The shared suite must terminate explicitly; load/syntax failures must fail. *)
Exit[1];
