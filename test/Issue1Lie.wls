repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
SetDirectory[repo];
Get["group.m"];

CommonFunctions`importPackage["GroupInfo`", "LieChecks`", {"id", "dim", "dual", "prod", "v", "isrep", "minrep", "gA", "setGroup"}];
CommonFunctions`importPackage["GroupInfoLie`", "LieChecks`", {"getLie", "getSO", "getSpin", "getSp"}];
LieChecks`eq = ClebschGordan`Private`eq;
LieChecks`ope = ClebschGordan`ope;
LieChecks`isPseaudo = ClebschGordan`isPseaudo;
Get["test/LieChecks.m"];

(* The shared suite must terminate explicitly; load/syntax failures must fail. *)
Exit[1];
