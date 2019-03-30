(* ::Package:: *)

(* ::Input:: *)
SetDirectory[ParentDirectory[NotebookDirectory[]]]


(* ::Input:: *)
<<"group.m"


(* ::Input:: *)
g = getGroup[8, 3];
setGroup[g];


(* ::Input:: *)
setOps[{op[e, rep[1]], op[v, rep[5]]}];
format[eq = bootAll[]]
sdp = makeSDP[eq];
py = toCboot[sdp];
buf = OpenWrite["D8.py"];
WriteString[buf, py];
Close[buf];


(* ::Input:: *)
opToTeX[e] = "\\epsilon";
opToTeX[v] = "v";
repToTeX[rep[1]] = "\\mathbf{I}^{+,+}";
repToTeX[rep[2]] = "\\mathbf{I}^{-,+}";
repToTeX[rep[3]] = "\\mathbf{I}^{+,-}";
repToTeX[rep[4]] = "\\mathbf{I}^{-,-}";
repToTeX[rep[5]] = "J";
Print[toTeX[eq]]
