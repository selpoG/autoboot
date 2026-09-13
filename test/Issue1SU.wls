repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
SetDirectory[repo];
Get["group.m"];

CommonFunctions`importPackage["GroupInfo`", "SUChecks`", {"id", "dim", "dual", "prod", "v", "isrep", "minrep", "gA", "setGroup"}];
SUChecks`get = GroupInfoLie`getSU;
SUChecks`eq = ClebschGordan`Private`eq;
SUChecks`ope = ClebschGordan`ope;
Get["test/SUChecks.m"];

(* The shared suite must terminate explicitly; load/syntax failures must fail. *)
Exit[1];
