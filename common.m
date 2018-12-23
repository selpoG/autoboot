BeginPackage["CommonFunctions`"]

e::usage = "e[n] gives n-th root of 1 with minimum positive argument.
e[n, m] gives m-th power of e[n]."
MyReap::usage = "MyReap[x] gives a list of values sowed in evaluation of x."
DoMany::usage = "DoMany[x,{y}] evaluates Do[x,y]. For example, DoMany[x,{{a,3},{b,4},{c,2,5}}] is equivalent to Do[x,{a,3},{b,4},{c,2,5}]."
dupcomb::usage = "dupcomb[x,k] gives a list of sorted selection of k elements from x with repetition."
removeSGs::usage = "removeSGs[] removes all SmallGroup and LieGroup packages from ContextPath."
lhs::usage = "lhs[x==y] gives left-hand side x."
rhs::usage = "rhs[x==y] gives right-hand side y."
importPackage::usage = "importPackage[from,to,symbol] imports 'symbols' in package 'from' to package 'to'. 'symbols' can be one symbol or a list of symbols."
zero::usage = "zero[x==y] gives x-y, which equation asserts to be zero."

(* any instance of a set or union-find tree has inner state and methods are destructive (not pure function). *)
newSet::usage = "newSet[] creates an instance of an empty set.
newSet[x] creates an instance of a set with elements in x."
set::usage = "set[_] represents an instance of a set."
newUF::usage = "newUF[] creates an instance of an empty union-find tree."
uf::usage = "uf[_] represents an instance of a union-find tree."
add::usage = "add[s, v] adds a element v to a set s.
add[u, v] adds a vertex v to a union-find tree s."
has::usage = "has[s, v] gives a set s has a element v or not."
remove::usage = "remove[s, v] removes a element v from a set s."
clear::usage = "clear[s] removes all elements from s."
keys::usage = "keys[s] gives a list of all elements in s."
size::usage = "size[s] gives a number of elements in s."
root::usage = "root[u, v] gives a root vertex of a vertex v in a union-find tree u."
unite::usage = "unite[u, v1, v2] puts vertices v1 and v2 in the same category in union-find tree u."
classify::usage = "classify[u] gives a list of connected components in union-find tree u. Each component is a list of vertices in the same category."

reverseIndex::usage = "reverseIndex[x] gives an association d s.t. d[x[[i+1]]]=i for all 0<=i<Length[x]. x is assumed to be a duplicate-free list."

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
