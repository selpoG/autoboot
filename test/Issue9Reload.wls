(* Regression test for https://github.com/selpoG/autoboot/issues/9. *)
repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
SetDirectory[repo];

reloadPassed = Check[
	Get["group.m"];
	Get["group.m"];
	True,
	False,
	{Set::wrsym}
];

ImportPackageTestSource`x = 7;
ImportPackageTestTarget`a = 1;
ImportPackageTestTarget`b = 2;
CommonFunctions`importPackage[
	"ImportPackageTestSource`",
	"ImportPackageTestTarget`",
	"x"
];
CommonFunctions`importPackage[
	"ImportPackageTestSource`",
	"ImportPackageTestTarget`",
	"x"
];
aliasPassed = ImportPackageTestTarget`x === 7;

Quiet[Check[
	CommonFunctions`importPackage[
		"ImportPackageTestSource`",
		"ImportPackageTestTarget`",
		"*"
	], Null]];
wildcardPassed = ImportPackageTestTarget`a === 1 &&
	ImportPackageTestTarget`b === 2;

passed = reloadPassed && aliasPassed && wildcardPassed;

Print["ISSUE9_RELOAD=", If[passed, "PASS", "FAIL"]];
Exit[If[passed, 0, 1]];
