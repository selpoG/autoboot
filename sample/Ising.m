(* ::Package:: *)

(* ::Input:: *)
SetDirectory[ParentDirectory[NotebookDirectory[]]]


(* ::Input:: *)
<<"group.m"


(* ::Input:: *)
z2 = getGroup[2, 1];
setGroup[z2];


(* ::Input:: *)
setOps[{op[e, z2[id]], op[a, rep[2]]}];
format[eq = bootAll[]]
sdp = makeSDP[eq];
py = toCboot[sdp];
buf = OpenWrite["Ising.py"];
WriteString[buf, py];
Close[buf];


(* ::Input:: *)
opToTeX[e] := "\\epsilon"
opToTeX[a] := "\\sigma"
repToTeX[rep[1]] := "I^+"
repToTeX[rep[2]] := "I^-"
Print[toTeX[eq]]
