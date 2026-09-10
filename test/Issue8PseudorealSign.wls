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

Print["ISSUE8_PSEUDOREAL_SIGN_EXACT=", If[exactPassed, "PASS", "FAIL"]];
Exit[If[exactPassed, 0, 1]];
