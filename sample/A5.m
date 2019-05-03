SetDirectory[ParentDirectory[NotebookDirectory[]]]


<<"ngroup.m"
setPrecision[100]


a5 = getGroup[60, 5];
setGroup[a5];


setOps[{op[e, rep[1]], op[v, rep[2]], op[t, rep[5]]}];
format[eq = bootAll[]]
sdp = makeSDP[eq];
py = toCboot[sdp];
buf = OpenWrite["A5.py"];
WriteString[buf, py];
Close[buf];


opToTeX[e] = "\\epsilon";
opToTeX[v] = "v";
opToTeX[t] = "t";
repToTeX[rep[1]] = "\\mathbf{I}";
repToTeX[rep[2]] = "\\mathbf{V}";
repToTeX[rep[3]] = "\\tilde{\\mathbf{V}}";
repToTeX[rep[4]] = "\\mathbf{4}";
repToTeX[rep[5]] = "\\mathbf{5}";
Print[toTeX[eq]]
