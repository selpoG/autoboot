(* ::Package:: *)

(* ::Input:: *)
SetDirectory[ParentDirectory[NotebookDirectory[]]]


(* ::Input:: *)
<<"group.m"


(* ::Input:: *)
o2 = getO[2];
setGroup[o2];


(* ::Input:: *)
setOps[{op[s, o2[id]], op[v, v[1]], op[t, v[2]]}];
format[eq = bootAll[]]
sdp = makeSDP[eq];
py = toCboot[sdp];
buf = OpenWrite["O2.py"];
WriteString[buf, py];
Close[buf];


(* ::Input:: *)
opToTeX[s] = "s";
opToTeX[v] = "\\phi";
opToTeX[t] = "t";
repToTeX[v[1]] = "\\mathbf{V}";
repToTeX[v[2]] = "\\mathbf{T}";
repToTeX[v[n_]] := TemplateApply["\\mathbf{S}^{``}", n]
repToTeX[i[n_]] := TemplateApply["\\mathbf{I}^``", If[n > 0, "+", "-"]]
Print[toTeX[eq]]
