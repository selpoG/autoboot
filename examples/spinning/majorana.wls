(* Run from the repository root:
   wolframscript -f examples/spinning/majorana.wls /path/to/tables /path/to/pmp.json
   This is an input-generation example, not a converged exclusion bound. *)
Get["spinning-identical.m"];
If[Length[$ScriptCommandLine]!=3,Print["Usage: majorana.wls TABLE_DIRECTORY OUTPUT.json"];Exit[2]];
dir=$ScriptCommandLine[[2]];output=$ScriptCommandLine[[3]];
files=Quiet[Check[Import[FileNameJoin[{dir,"manifest.json"}],"RawJSON"],$Failed]];
If[!VectorQ[files,StringQ],Exit[1]];
tables=readSpinningBlocks[FileNameJoin[{dir,#}],60]&/@files;
If[!AllTrue[tables,AssociationQ],Exit[1]];
spins=Union[Lookup[tables,"j_internal"]];lambda=Min[Lookup[tables,"lambda"]];
(* Strictly above the long-block poles. A physical local CFT analysis also
   needs a separate conserved stress-tensor block with Ward constraints. *)
sectors=Flatten[Table[<|"Spin"->j,"Parity"->p,"Gap"->If[j==0,1,j+11/10]|>,{j,spins},{p,{-1,1}}]];
sdp=spinIdenticalSDP[tables,sectors,x,lambda];If[!AssociationQ[sdp],Exit[1]];
sdp=spinMajoranaRows[sdp];If[!AssociationQ[sdp],Exit[1]];
If[spinWritePMP[output,sdp,x,40]===$Failed,Exit[1]];
Print["Wrote ",output," (",Length[sdp["Identity"]]," equations; ",Length[sdp["Sectors"]]," long-operator sectors)"];
