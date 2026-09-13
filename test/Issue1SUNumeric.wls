repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
SetDirectory[repo];
Get["ngroup.m"];
NGroupInfo`setPrecision[50];
CommonFunctions`importPackage["NGroupInfo`", "SUChecks`", {"id", "dim", "dual", "prod", "v", "isrep", "minrep", "gA", "setGroup"}];
SUChecks`get = NGroupInfoLie`getSU;
SUChecks`eq = NClebschGordan`Private`eq;
SUChecks`ope = NClebschGordan`ope;
Get["test/SUChecks.m"];

(* The shared suite must terminate explicitly; load/syntax failures must fail. *)
Exit[1];
