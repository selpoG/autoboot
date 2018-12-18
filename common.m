BeginPackage["CommonFunctions`"]

e::usage = "e[n]
e[n, m]"
MyReap::usage = "MyReap[x]"
DoMany::usage = "DoMany[x,{y}]"
dupcomb::usage = "dupcomb[x,k]"
removeSGs::usage = "removeSGs[]"
lhs::usage = "lhs[x==y]"
rhs::usage = "rhs[x==y]"
importPackage::usage = "importPackage[from,to,symbol]
importPackage[from,to,symbols]"
zero::usage = "zero[x==y]"

newSet::usage = "newSet[]
newSet[x]"
set::usage = "set[_]"
newUF::usage = "newUF[]"
uf::usage = "uf[_]"
add::usage = "add[s, v]
add[u, v]"
has::usage = "has[s, v]"
remove::usage = "add[s, v]"
clear::usage = "clear[s]"
keys::usage = "keys[s]"
size::usage = "size[s]"
root::usage = "root[u, v]"
unite::usage = "unite[u, u, v]"
classify::usage = "classify[u]"

reverseIndex::usage = "reverseIndex[x]"

Begin["`Private`"]

MyReap[x_] := Module[{t}, t = Reap[x][[2]]; If[Length[t] == 0, {}, t[[1]]]]
SetAttributes[MyReap, HoldFirst]

e[n_] := Cos[2 Pi/n] + I Sin[2 Pi/n]
e[n_, m_] := Cos[2 Pi m/n] + I Sin[2 Pi m/n]

DoMany[x_, {y__}] := Do[x, y]
SetAttributes[DoMany, HoldFirst]

dupcombHelper[n_, k_] := Module[{i, j}, MyReap[DoMany[Sow[Array[i, k]], Join[{{i[1], n}}, Table[{i[j], i[j - 1], n}, {j, 2, k}]]]]]
dupcomb[x_List, k_] := Module[{i, ind}, Table[Table[x[[i]], {i, ind}], {ind, dupcombHelper[Length[x], k]}]]

pat = StringMatchQ["SmallGroup" | "LieGroup" ~~ "*`"]
removeSGs[] := ($ContextPath = MyReap @ Do[If[! pat[x], Sow[x]], {x, $ContextPath}];)

lhs[a_ == _] := a
rhs[_ == b_] := b
SetAttributes[lhs, Listable]
SetAttributes[rhs, Listable]
zero2[x_ == y_] := x - y
zero[x_And] := MyReap[Scan[Sow[zero2[#]] &, List @@ x]]
zero[x_] := {zero2[x]}
zero[True] = {}

importPackage[from_, to_, symbols_List] := Scan[importPackage[from, to, #]&, symbols]
importPackage[from_, to_, symbol_] := (Evaluate[Symbol[to<>symbol]] = Symbol[from<>symbol];)

newSet[] := Module[{s}, vs[s] ^= <||>; sz[s] ^= 0; set[s]]
newSet[x__] := Module[{s, ans, t}, vs[s] ^= <||>; sz[s] ^= 0; ans = set[s]; Do[add[ans, t], {t, List[x]}]; ans]
has[set[s_], v_] ^:= KeyExistsQ[vs[s], v]
add[x : set[s_], v_] ^:= If[! has[x, v], AbortProtect[vs[s] ^= AssociateTo[vs[s], v -> 1]; sz[s] ^= sz[s] + 1;]]
remove[x : set[s_], v_] ^:= If[has[x, v], AbortProtect[vs[s] ^= KeyDropFrom[vs[s], v]; sz[s] ^= sz[s] - 1;]]
clear[set[s_]] ^:= AbortProtect[vs[s] ^= <||>; sz[s] ^= 0;]
keys[set[s_]] ^:= Keys[vs[s]]
size[set[s_]] ^:= sz[s]

newUF[] := Module[{par}, vs[par] ^= newSet[]; uf[par]]
add[uf[par_], v_] ^:= If[! has[vs[par], v], AbortProtect[add[vs[par], v]; par[v] = v;]]
root[x : uf[par_], v_] ^:= (add[x, v]; If[par[v] =!= v, par[v] = root[x, par[v]]]; par[v]);
unite[x : uf[par_], u_, v_] ^:= (par[root[x, u]] = root[x, v];)
classify[x : uf[par_]] ^:= Module[{rt = <||>, v, r}, Do[r = root[x, v]; If[! KeyExistsQ[rt, r], rt[r] = newSet[]]; add[rt[r], v], {v, keys[vs[par]]}]; keys /@ List @@ rt]

reverseIndex[x_List] := Module[{d = <||>}, Scan[(d[#] = Length[d]) &, x]; d]

End[ ]

EndPackage[ ]
