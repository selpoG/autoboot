(* Numerical regression test for https://github.com/selpoG/autoboot/issues/8. *)
repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
SetDirectory[repo];

Get["ngroup.m"];
NGroupInfo`setPrecision[50];
numericGroup = NGroupInfoD`getDicyclic[2];
NGroupInfo`setGroup[numericGroup];
numericA = NClebschGordan`op[a, NGroupInfo`i[1, 1], 1, 1];
numericB = NClebschGordan`op[b, NGroupInfo`v[1], 1, 1];
NClebschGordan`setOps[{numericA, numericB}];

Unprotect[NClebschGordan`\[Alpha]];
NClebschGordan`\[Alpha][numericA, numericB, numericA, numericB][
	NClebschGordan`op[op, NGroupInfo`v[1], 1, _], 1, 1] = 7;
NClebschGordan`\[Alpha][numericA, numericB, numericA, numericB][
	NClebschGordan`op[op, NGroupInfo`v[1], -1, _], 1, 1] = 2;
Protect[NClebschGordan`\[Alpha]];

numericResult = NClebschGordan`Private`sF[
	numericA, numericB, numericA, numericB][NGroupInfo`v[1], 1, 1, 1];
numericCoefficient = numericResult /.
	NClebschGordan`sum[x_, _] :> x /.
	NClebschGordan`H[__] :> 1;
numericIntermediate = NClebschGordan`op[
	op, NGroupInfo`v[1], 1, 1];
numericPassed = NClebschGordan`isPseaudo[NGroupInfo`v[1]] &&
	NClebschGordan`\[Sigma][numericIntermediate] === 1 &&
	NClebschGordan`\[Sigma][
		NClebschGordan`dualOp[numericIntermediate]] === -1 &&
	numericCoefficient === 10;

Print["ISSUE8_PSEUDOREAL_SIGN_NUMERIC=", If[numericPassed, "PASS", "FAIL"]];
Exit[If[numericPassed, 0, 1]];
