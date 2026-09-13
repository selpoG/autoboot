(* Spacetime tensor structures are independent of internal-symmetry irreps.
   Conventions: arXiv:2011.01959 sections 2-3 and arXiv:1612.08987 appendix B. *)
BeginPackage["SpinningBootstrap`"];
spinThreePointBasis::usage = "spinThreePointBasis[{j1,j2,j3}, parity, exchange] returns q- and SO(3)-basis coefficients. parity is +/-1; exchange is 0 (distinct operators) or the required (12) eigenvalue +/-1. Optional homogeneous constraints act on q coefficients.";
spinQBasis::usage = "spinQBasis[spins] enumerates three-point (sum q=0) or four-point q labels for exact nonnegative half-integral spins.";
spinSO3Basis::usage = "spinSO3Basis[{j1,j2,j3}] lists {j12,j123} labels.";
spinQToSO3::usage = "spinQToSO3[spins] gives the inverse basis transformation using CG orthogonality.";
spinSO3ToQ::usage = "spinSO3ToQ[spins] gives rows expressing SO(3) structures in the q basis, in blocks_3d conventions.";
spinCrossingMatrix::usage = "spinCrossingMatrix[spins, parity] returns the (13) crossing matrix and q labels. Includes the fermionic permutation sign. parity=0 retains both space parities.";
spinCrossingEquations::usage = "spinCrossingEquations[spins, parity, s, t, z, zb] generates component crossing equations for ordered coefficient functions s[q,z,zb] and t[q,z,zb].";
spinIdentityBlock::usage = "spinIdentityBlock[{j1,j2,j3,j4}, delta1, q, z, zb] gives the unit-normalized (12)(34) identity block when each pair consists of identical operators.";
spinMakeSDP::usage = "spinMakeSDP[identity, sectors] converts real OPE quadratic forms into symmetric matrices. Each sector has Couplings and Equations, and optionally Prefactor={constant,poles,base}.";
spinWritePMP::usage = "spinWritePMP[file, sdp, x, digits] writes polynomial matrices to the SDPB 3.x PMP JSON format. All blocks must already be evaluated as real polynomials in x, on x>=0.";
spinThreePointBasis::input = "Invalid spins, parity, exchange symmetry, or q-basis constraint matrix.";
spinCrossingMatrix::input = "Expected four exact half-integral spins, an even number of fermions, and parity 0 or +/-1.";
spinMakeSDP::input = "Expected a nonempty identity vector and homogeneous quadratic forms in distinct real OPE variables, with matching equation counts.";
spinWritePMP::input = "Cannot export: expected real polynomial matrices, sufficient precision, a nonzero identity normalization, and a positive pole-free prefactor on x>=0.";
Begin["`Private`"];
spinQ[j_] := IntegerQ[2 j] && TrueQ[j>=0];
validSpins[js_,n_] := ListQ[js] && Length[js]==n && AllTrue[js,spinQ] && IntegerQ[Total[js]];
spinQBasis[js_List] /; MemberQ[{3,4},Length[js]] && validSpins[js,Length[js]] :=
 If[Length[js]==3,Select[Tuples[Range[-#, #]& /@ js],Total[#]==0&],Tuples[Range[-#, #]& /@ js]];
spinQBasis[_] := (Message[spinThreePointBasis::input];$Failed);
spinSO3Basis[js_List] /; validSpins[js,3] := Flatten[Table[{k,l},
 {k,Abs[js[[1]]-js[[2]]],js[[1]]+js[[2]]}, {l,Abs[js[[3]]-k],js[[3]]+k}],1];
spinSO3Basis[_] := (Message[spinThreePointBasis::input];$Failed);
spinSO3ToQ[js_List] /; validSpins[js,3] := spinSO3ToQ[js] = Module[{qs=spinQBasis[js],bs=spinSO3Basis[js]},
 Table[If[Abs[q[[3]]]>b[[1]],0,(-1)^(js[[1]]-js[[3]]+q[[2]]) Sqrt[Times@@MapThread[Binomial[2 #1,#1+#2]&,{js,q}]]
  ClebschGordan[{js[[1]],q[[1]]},{js[[2]],q[[2]]},{b[[1]],-q[[3]]}]
  ClebschGordan[{b[[1]],-q[[3]]},{js[[3]],q[[3]]},{b[[2]],0}]],{b,bs},{q,qs}]];
spinSO3ToQ[_] := (Message[spinThreePointBasis::input];$Failed);
spinQToSO3[js_List] /; validSpins[js,3] := spinQToSO3[js] =
 RootReduce[DiagonalMatrix[(1/Times@@MapThread[Binomial[2 #1,#1+#2]&,{js,#}])&/@spinQBasis[js]].Transpose[spinSO3ToQ[js]]];
spinQToSO3[_] := (Message[spinThreePointBasis::input];$Failed);
spinThreePointBasis[js_,p_,exchange_:0,constraints_:{}] := Module[{qs,d,reflect,swap,eqs,rows,so3},
 If[!validSpins[js,3] || !MemberQ[{-1,1},p] || !MemberQ[{-1,0,1},exchange] ||
  (exchange!=0 && js[[1]]!=js[[2]]),Message[spinThreePointBasis::input];Return[$Failed]];
 qs=spinQBasis[js];d=Length[qs];
 If[constraints=!={} && (!MatrixQ[constraints] || Last[Dimensions[constraints]]!=d),
  Message[spinThreePointBasis::input];Return[$Failed]];
 reflect=Table[Boole[q===-r],{q,qs},{r,qs}];eqs=reflect-p IdentityMatrix[d];
 If[exchange!=0,
  swap=Table[(-1)^(js[[1]]+js[[2]]-js[[3]]) Boole[q===-r[[{2,1,3}]]],{q,qs},{r,qs}];
  eqs=Join[eqs,swap-exchange IdentityMatrix[d]]];
 rows=NullSpace[Join[eqs,constraints]];so3=spinSO3ToQ[js];
 <|"Spins"->js,"Parity"->p,"QBasis"->qs,"SO3Basis"->spinSO3Basis[js],
  "QCoefficients"->rows,"SO3Coefficients"->If[rows=={},{},RootReduce[rows.spinQToSO3[js]]],
  "RealityPhase"->If[AllTrue[js,IntegerQ],1,I]|>];
spinCrossingMatrix[js_,p_:0] := Module[{qs,ts,fs,statistics,mat},
 If[!validSpins[js,4] || !MemberQ[{-1,0,1},p],Message[spinCrossingMatrix::input];Return[$Failed]];
 qs=Select[spinQBasis[js],p==0 || (-1)^Total[js-#]==p&];ts=#[[{3,2,1,4}]]&/@qs;
 fs=Mod[2 js,2];statistics=(-1)^(fs[[1]] fs[[2]]+fs[[1]] fs[[3]]+fs[[2]] fs[[3]]);
 mat=SparseArray[DiagonalMatrix[statistics (-1)^(#[[1]]+#[[2]]+#[[3]]-#[[4]])&/@qs]];
 <|"Spins"->js,"CrossedSpins"->js[[{3,2,1,4}]],"Structures"->qs,"CrossedStructures"->ts,"Matrix"->mat|>];
spinCrossingEquations[js_,p_,s_,t_,z_,zb_] := Module[{c=spinCrossingMatrix[js,p]},
 If[c===$Failed,Return[$Failed]];
 (s[#,z,zb]&/@c["Structures"])-c["Matrix"].(t[#,1-z,1-zb]&/@c["CrossedStructures"])];
spinIdentityBlock[js_,delta_,q_,z_,zb_] /; validSpins[js,4] && js[[1]]==js[[2]] && js[[3]]==js[[4]] && MemberQ[spinQBasis[js],q] :=
 If[q[[1]]!=q[[2]] || q[[3]]!=q[[4]],0,
  I^(2 q[[1]]+2 q[[4]]) Binomial[2 js[[1]],js[[1]]+q[[1]]] Binomial[2 js[[4]],js[[4]]+q[[4]]]
  z^(-delta+q[[1]]) zb^(-delta-q[[1]])];
spinIdentityBlock[___] := (Message[spinThreePointBasis::input];$Failed);
(* The symmetric matrix uses half the coefficient of an off-diagonal monomial.
   Its quadratic form is exactly the original expression, including all copies. *)
spinMakeSDP[identity_List,sectors_Association] := Module[{out=<||>,n=Length[identity],vars,eqs,mats,bad=False},
 If[n==0,Message[spinMakeSDP::input];Return[$Failed]];
 KeyValueMap[Function[{label,data},
  If[!AssociationQ[data],bad=True;Return[]];
  vars=Lookup[data,"Couplings",{}];eqs=Lookup[data,"Equations",{}];
  If[!ListQ[vars] || vars=={} || !DuplicateFreeQ[vars] || !AllTrue[vars,MatchQ[#,_Symbol]&] ||
   !ListQ[eqs] || Length[eqs]!=n || !AllTrue[eqs,PolynomialQ[#,vars]&],bad=True;Return[]];
  mats=Table[Table[Expand[D[e,a,b]/2],{a,vars},{b,vars}],{e,eqs}];
  If[!FreeQ[mats,Alternatives@@vars] || !AllTrue[eqs,AllTrue[First/@CoefficientRules[#,vars],Total[#]==2&]&],bad=True;Return[]];
  out[label]=<|"Couplings"->vars,"Matrices"->mats,"Prefactor"->Lookup[data,"Prefactor",{1,{},1}]|>
 ],sectors];
 If[bad,Message[spinMakeSDP::input];$Failed,<|"Identity"->identity,"Sectors"->out|>]];
spinMakeSDP[___] := (Message[spinMakeSDP::input];$Failed);
realNumberQ[v_,digits_] := NumericQ[v] && TrueQ[Im[N[v,digits]]==0] &&
 FreeQ[v,Indeterminate|_DirectedInfinity] && (Precision[v]===Infinity || TrueQ[Precision[v]>=digits] || (Sign[v]===0 && TrueQ[Accuracy[v]>=digits]));
decimal[v_,digits_] := ToString[NumberForm[N[Re[v],digits],digits,NumberPadding->{"",""},ExponentFunction->(Null&)],OutputForm];
spinWritePMP[file_String,sdp_Association,x_Symbol,digits_Integer:50] := Module[{id,blocks={},bad=False,mat,pref,coeffs,polys,json},
 id=Lookup[sdp,"Identity",{}];
 If[digits<20 || !StringEndsQ[file,".json"] || !ListQ[id] || id=={} || !AllTrue[id,realNumberQ[#,digits]&] || AllTrue[id,TrueQ[#==0]&] ||
  !AssociationQ[Lookup[sdp,"Sectors",None]],Message[spinWritePMP::input];Return[$Failed]];
 KeyValueMap[Function[{label,data},
  mat=Lookup[data,"Matrices",{}];pref=Lookup[data,"Prefactor",{}];
  If[Length[Dimensions[mat]]!=3 || First[Dimensions[mat]]!=Length[id] || Dimensions[mat][[2]]!=Dimensions[mat][[3]] ||
   !AllTrue[Flatten[CoefficientList[#,x]&/@Flatten[Expand[#-Transpose[#]]&/@mat]],Sign[#]===0&] ||
   !AllTrue[Flatten[mat],PolynomialQ[#,x]&] || !MatchQ[pref,{_,_List,_}],bad=True;Return[]];
  If[!AllTrue[Flatten[pref],realNumberQ[#,digits]&] || !TrueQ[pref[[1]]>0] || !TrueQ[pref[[3]]>0] || !AllTrue[pref[[2]],TrueQ[#<0]&],bad=True;Return[]];
  polys=Map[If[TrueQ[#==0],{0},CoefficientList[#,x]]&,Transpose[mat,{3,1,2}],{3}];
  If[!AllTrue[Flatten[polys],realNumberQ[#,digits]&],bad=True;Return[]];
  AppendTo[blocks,<|"prefactor"-><|"constant"->decimal[pref[[1]],digits],"poles"->(decimal[#,digits]&/@pref[[2]]),"base"->decimal[pref[[3]],digits]|>,
   "polynomials"->Map[decimal[#,digits]&,polys,{4}]|>]
 ],sdp["Sectors"]];
 If[bad || blocks=={},Message[spinWritePMP::input];Return[$Failed]];
 json=<|"objective"->ConstantArray["0",Length[id]],"normalization"->(decimal[#,digits]&/@id),"PositiveMatrixWithPrefactorArray"->blocks|>;
 Export[file,json,"RawJSON"]];
spinWritePMP[___] := (Message[spinWritePMP::input];$Failed);
End[];
EndPackage[];
