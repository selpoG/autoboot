(* Regression test for https://github.com/selpoG/autoboot/issues/9. *)
repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
SetDirectory[repo];

passed = Check[
	Get["group.m"];
	Get["group.m"];
	True,
	False,
	{Set::wrsym}
];

Print["ISSUE9_RELOAD=", If[passed, "PASS", "FAIL"]];
Exit[If[passed, 0, 1]];
