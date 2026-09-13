(* O(n), n>=4: tensor representations of SO(n), including the disconnected
   component. The last label is +/-1 for an extension, or 0 for an induced pair. *)
Needs["LieRepresentations`", "lie-representations.m"];
BeginPackage["OrthogonalRepresentations`"];
oValid::usage="Internal O(n) label validator.";
oDimension::usage="Dimension of an O(n) irrep.";
oMatrices::usage="Lowering and raising matrices of an O(n) irrep.";
oReflection::usage="Disconnected generator: central inversion for odd n, diagram reflection for even n.";
oProduct::usage="O(n) tensor product, including reflection multiplicities.";
oIntertwiners::usage="Exact O(n) CG solutions satisfying Lie and reflection equations.";
Begin["`Private`"];
type[n_]:=If[OddQ[n],"B","D"];
swap[a_List]:=Join[Drop[a,-2],Reverse[Take[a,-2]]];
canonical[a_]:=If[a[[-2]]>=a[[-1]],a,swap[a]];
oValid[n_,q_List]:=IntegerQ[n] && n>=4 && Length[q]==Floor[n/2]+1 &&
 LieRepresentations`validLabel[Floor[n/2],Most[q]] &&
 If[OddQ[n],EvenQ[q[[-2]]] && MemberQ[{-1,1},Last[q]],
  EvenQ[Total[Take[Most[q],-2]]] && If[q[[-3]]==q[[-2]],MemberQ[{-1,1},Last[q]],q[[-3]]>q[[-2]] && Last[q]==0]];
oValid[_,_]:=False;
parts[n_,q_]:=If[EvenQ[n] && Last[q]==0,{Most[q],swap[Most[q]]},{Most[q]}];
oDimension[n_,q_]:=Total[LieRepresentations`repDimension[type[n],Floor[n/2],#]&/@parts[n,q]];
blockDiagonal[blocks_]:=SparseArray[ArrayFlatten[Table[If[i==j,blocks[[i]],ConstantArray[0,{Length[blocks[[i]]],Length[blocks[[j]]]}]],{i,Length[blocks]},{j,Length[blocks]}]]];
oMatrices[n_,q_]:=oMatrices[n,q]=Module[{ms=LieRepresentations`repMatrices[type[n],Floor[n/2],#]&/@parts[n,q]},
 Table[blockDiagonal[#[[i]]&/@ms],{i,2 Floor[n/2]}]];
weights[n_,q_]:=(oMatrices[n,q];Join@@(LieRepresentations`Private`basisWeights[type[n],Floor[n/2],#]&/@parts[n,q]));
(* Apply the same target basis changes used by the Lie representation builder. *)
propagate[ops_,n_,a_,h_]:=Module[{current=ArrayReshape[SparseArray[h],{Length[h],1}],images,next,candidates},
 LieRepresentations`repMatrices[type[n],Floor[n/2],a];images=current;
 Do[candidates=Join[Sequence@@(#.current&/@Take[ops,Floor[n/2]]),2];
  next=RootReduce[candidates[[All,step[[1]]]].Transpose[step[[2]]]];
  images=Join[images,next,2];current=next,
 {step,LieRepresentations`Private`basisSteps[type[n],Floor[n/2],a]}];images];
(* Source a -> tau(a), fixing the highest vector. Consequently J_tau J_a=I. *)
diagramMap[n_,a_]:=diagramMap[n,a]=Module[{rank=n/2,perm,ms},
 perm=Join[Range[rank-2],{rank,rank-1}];
 ms=LieRepresentations`repMatrices["D",rank,swap[a]];
 propagate[Join[ms[[perm]],ms[[rank+perm]]],n,a,UnitVector[Length[ms[[1]]],1]]];
oReflection[n_,q_]:=oReflection[n,q]=If[OddQ[n],Last[q] IdentityMatrix[oDimension[n,q],SparseArray],
 If[Last[q]!=0,Last[q] diagramMap[n,Most[q]],
  With[{j=diagramMap[n,Most[q]],d=oDimension[n,q]/2},
   SparseArray[ArrayFlatten[{{ConstantArray[0,{d,d}],Transpose[j]},{j,ConstantArray[0,{d,d}]}}]]]]];
tensorOps[n_,a_,b_]:=tensorOps[n,a,b]=Module[{ma=oMatrices[n,a],mb=oMatrices[n,b],da=oDimension[n,a],db=oDimension[n,b]},
 Table[KroneckerProduct[ma[[i]],IdentityMatrix[db,SparseArray]]+KroneckerProduct[IdentityMatrix[da,SparseArray],mb[[i]]],{i,Length[ma]}]];
(* Highest vectors, optionally restricted to a reflection eigenvalue. *)
highest[n_,a_,b_,target_,parity_]:=highest[n,a,b,target,parity]=Module[{pairs,cols,ops,blocks,block,active,rr,sol,d=oDimension[n,a] oDimension[n,b]},
 pairs=Flatten[Table[u+v,{u,weights[n,a]},{v,weights[n,b]}],1];cols=Flatten[Position[pairs,target,{1}]];
 If[cols=={},Return[{}]];ops=tensorOps[n,a,b];
 blocks=Table[block=x[[All,cols]];active=Union[First/@(First/@Most[ArrayRules[block]])];
  If[active=={},Nothing,block[[active]]],{x,Drop[ops,Floor[n/2]]}];
 If[parity!=0,rr=KroneckerProduct[oReflection[n,a],oReflection[n,b]][[cols,cols]]-parity IdentityMatrix[Length[cols],SparseArray];AppendTo[blocks,rr]];
 block=If[blocks=={},SparseArray[{}, {1,Length[cols]}],Join@@blocks];sol=NullSpace[RootReduce[block]];
 SparseArray[Thread[cols->#],d]&/@sol];
oProduct[n_,a_,b_]:=oProduct[n,a,b]=Module[{restr,targets,counts,result={},m,c},
 If[OddQ[n],Return[Append[#,Last[a] Last[b]]&/@LieRepresentations`repProduct[type[n],Floor[n/2],Most[a],Most[b]]]];
 restr=Flatten[Table[LieRepresentations`repProduct["D",n/2,u,v],{u,parts[n,a]},{v,parts[n,b]}],2];
 counts=Counts[restr];targets=DeleteDuplicates[canonical/@restr];
 Do[If[c=!=swap[c],result=Join[result,ConstantArray[Append[c,0],counts[c]]],
  Do[m=Length[highest[n,a,b,c,p]];result=Join[result,ConstantArray[Append[c,p],m]],{p,{1,-1}}]],{c,targets}];result];
oIntertwiners[n_,a_,b_,c_]:=Module[{hs,ops,maps,first,second,rr,j},
 If[OddQ[n],Return[If[Last[a] Last[b]==Last[c],LieRepresentations`repIntertwiners[type[n],Floor[n/2],Most[a],Most[b],Most[c]],{}]]];
 hs=highest[n,a,b,Most[c],Last[c]];If[hs=={},Return[{}]];ops=tensorOps[n,a,b];
 If[Last[c]==0,rr=KroneckerProduct[oReflection[n,a],oReflection[n,b]];j=diagramMap[n,Most[c]]];
 maps=Table[first=propagate[ops,n,Most[c],h];
  If[Last[c]==0,second=RootReduce[rr.first.Transpose[j]];first=Join[first,second,2]];
  Flatten[Normal[first]],{h,hs}];maps];
End[];
EndPackage[];
