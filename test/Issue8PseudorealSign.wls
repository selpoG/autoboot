(* Regression test for https://github.com/selpoG/autoboot/issues/8. *)
repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
SetDirectory[repo];

Get["group.m"];
exactGroup = GroupInfoD`getDicyclic[2];
GroupInfo`setGroup[exactGroup];
exactA = ClebschGordan`op[a, GroupInfo`i[1, 1], 1, 1];
exactB = ClebschGordan`op[b, GroupInfo`v[1], 1, 1];
ClebschGordan`setOps[{exactA, exactB}];

Unprotect[ClebschGordan`\[Alpha]];
ClebschGordan`\[Alpha][exactA, exactB, exactA, exactB][
	ClebschGordan`op[op, GroupInfo`v[1], 1, _], 1, 1] = 7;
ClebschGordan`\[Alpha][exactA, exactB, exactA, exactB][
	ClebschGordan`op[op, GroupInfo`v[1], -1, _], 1, 1] = 2;
Protect[ClebschGordan`\[Alpha]];

exactResult = ClebschGordan`Private`sF[
	exactA, exactB, exactA, exactB][GroupInfo`v[1], 1, 1, 1];
exactCoefficient = exactResult /.
	ClebschGordan`sum[x_, _] :> x /.
	ClebschGordan`H[__] :> 1;
exactIntermediate = ClebschGordan`op[
	op, GroupInfo`v[1], 1, 1];
exactPassed = ClebschGordan`isPseaudo[GroupInfo`v[1]] &&
	ClebschGordan`\[Sigma][exactIntermediate] === 1 &&
	ClebschGordan`\[Sigma][
		ClebschGordan`dualOp[exactIntermediate]] === -1 &&
	exactCoefficient === 10;

Quiet[Get["ngroup.m"], General::shdw];
NGroupInfo`setPrecision[50];
numericGroup = NGroupInfoD`getDicyclic[2];
Quiet[NGroupInfo`setGroup[numericGroup], General::shdw];
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

passed = exactPassed && numericPassed;
Print["ISSUE8_PSEUDOREAL_SIGN=", If[passed, "PASS", "FAIL"]];
Exit[If[passed, 0, 1]];
