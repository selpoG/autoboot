(* Parity-invariant four-point functions of one Hermitian 3d primary.
   Long external representations; conservation/Ward identities are additional
   constraints and are deliberately not inferred from the scaling dimension. *)
Needs["SpinningBootstrap`","spinning.m"];
Needs["SpinningBlocks`","spinning-blocks.m"];
BeginPackage["SpinningIdentical`",{"SpinningBootstrap`","SpinningBlocks`"}];
spinIdenticalSDP::usage="spinIdenticalSDP[tables, sectors, x, lambda] builds crossing from blocks_3d tables. Each sector specifies Spin, Parity, Gap. Includes all q components, both z/zb parities, and real OPE matrices. Rows may be redundant.";
spinMajoranaRows::usage="spinMajoranaRows[sdp] selects the independent Majorana derivative components of arXiv:1612.08987 appendix A.1, excluding the dependent ++++ line derivatives. Requires ExternalSpin=1/2.";
spinIndependentRows::usage="spinIndependentRows[sdp,x,tolerance] selects numerically independent equation rows from all polynomial coefficients and the identity. Returns the retained row indices and relative residuals; tolerance must be explicitly chosen.";
spinIdenticalSDP::input="Missing, inconsistent, or invalid blocks_3d tables or sector parameters.";
spinIndependentRows::input="Expected finite real numerical polynomial matrices and a positive numerical tolerance.";
Begin["`Private`"];
representative[q_] := If[Select[q,#!=0&]=={} || First[Select[q,#!=0&]]>0,q,-q];
(* Taylor coefficients of (.5+x+y)^a (.5+x-y)^b. t=y^2;
   the odd component is divided by y before taking xt derivatives. *)
powerDerivative[a_,b_,m_,k_] := 2^(-a-b+m+k) Sum[
 Binomial[a,i+h] Binomial[i+h,i] Binomial[b,m-i+k-h] Binomial[m-i+k-h,m-i] (-1)^(k-h),{i,0,m},{h,0,k}];
identityDerivative[j_,delta_,q_,s_,m_,n_] := If[q[[1]]!=q[[2]] || q[[3]]!=q[[4]],0,
 I^(2 q[[1]]+2 q[[4]]) Binomial[2 j,j+q[[1]]] Binomial[2 j,j+q[[4]]]
 m! n! powerDerivative[-delta+q[[1]],-delta-q[[1]],m,2 n+Boole[s==-1]]];
spinIdenticalSDP[tables_List,sectors_List,x_Symbol,lambda_Integer] := Module[
 {ans,first,j,delta,digits,qs,rows,id,out=<||>,index=<||>,key,tag,bas,c,bs,vars,polys,ref,gap,l,p,raw,cross,q,t,s,m,n,phase,eqs,poleSets,commonPoles,missingPoles},
 ans=Catch[
 If[tables=={} || sectors=={} || lambda<0 || !AllTrue[tables,AssociationQ],Throw[$Failed]];
 first=First[tables];j=first["j_external"][[1]];delta=first["delta_1_plus_2"]/2;digits=Min[Lookup[tables,"Digits"]];
 If[!IntegerQ[2 j] || j<0 || !TrueQ[delta>If[j<1,j+1/2,j+1]],Throw[$Failed]];
 Do[If[!And@@(data[#]===first[#]&/@{"j_external","delta_12","delta_43","delta_1_plus_2","order","kept_pole_order"}) ||
   data["j_external"]!=ConstantArray[j,4] || data["delta_12"]!=0 || data["delta_43"]!=0 || data["lambda"]<lambda,Throw[$Failed]];
  key=ToString[{data["j_internal"],data["four_pt_struct"],data["four_pt_sign"],data["j_12"],data["j_43"]},InputForm];
  If[KeyExistsQ[index,key],Throw[$Failed]];index[key]=data,{data,tables}];
 qs=DeleteDuplicates[representative/@Select[spinQBasis[ConstantArray[j,4]],(-1)^Total[ConstantArray[j,4]-#]==1&]];
 rows=Flatten[Table[Table[{q,s,m,n},{m,0,lambda-Boole[s==-1]},{n,0,Floor[(lambda-Boole[s==-1]-m)/2]}],{q,qs},{s,{1,-1}}],3];
 phase[q_]:=(-1)^(q[[1]]+q[[2]]-q[[3]]-q[[4]]);
 cross[q_,s_,m_]:=phase[q] (-1)^m s If[representative[q[[{3,2,1,4}]]]===q[[{3,2,1,4}]],1,s];
 id=N[Map[Function[row,With[{q0=row[[1]],s0=row[[2]],m0=row[[3]],n0=row[[4]]},
 identityDerivative[j,delta,q0,s0,m0,n0]-cross[q0,s0,m0] identityDerivative[j,delta,representative[q0[[{3,2,1,4}]]],s0,m0,n0]]],rows],digits+5];
 Do[
  If[!AssociationQ[sector] || !And@@(KeyExistsQ[sector,#]&/@{"Spin","Parity","Gap"}),Throw[$Failed]];
  l=sector["Spin"];p=sector["Parity"];gap=sector["Gap"];
  If[!IntegerQ[l] || l<0 || !MemberQ[{-1,1},p] || !NumericQ[gap] || !TrueQ[gap>If[l==0,1/2,l+1]],Throw[$Failed]];
  bas=spinThreePointBasis[{j,j,l},p,(-1)^(2 j)];c=bas["SO3Coefficients"];bs=bas["SO3Basis"];
  If[c=={},Continue[]];vars=Table[Unique["ope"],{Length[c]}];
  ref=SelectFirst[tables,#["j_internal"]==l&,Missing[]];If[MissingQ[ref],Throw[$Failed]];
  poleSets=Lookup[Select[tables,#["j_internal"]==l&],"pole_list_x"];
  commonPoles=Flatten[Table[ConstantArray[pole,Max[Count[#,pole]&/@poleSets]],{pole,Union[Flatten[poleSets]]}]];
  ref["pole_list_x"]=commonPoles;
  raw[q_,s_,m_,n_]:=Module[{a,b,data,poly,total=0},
   Do[If[AllTrue[c[[All,a]],TrueQ[#==0]&] || AllTrue[c[[All,b]],TrueQ[#==0]&],Continue[]];
    key=ToString[{l,q,s,bs[[a,1]],bs[[b,1]]},InputForm];If[!KeyExistsQ[index,key],Throw[$Failed]];data=index[key];
    If[data["delta_minus_x"]=!=ref["delta_minus_x"],Throw[$Failed]];
    poly=spinBlockPolynomial[data,{bs[[a,2]],bs[[b,2]]},{m,n},x];If[poly===$Failed,Throw[$Failed]];
    missingPoles=Fold[DeleteCases[#1,#2,1,1]&,commonPoles,data["pole_list_x"]];
    poly*=Times@@(x-#&/@missingPoles);
    total+=bas["RealityPhase"]^2 (vars.c[[All,a]]) (vars.c[[All,b]]) poly,{a,Length[bs]},{b,Length[bs]}];Expand[total]];
  eqs=Map[Function[row,With[{q0=row[[1]],s0=row[[2]],m0=row[[3]],n0=row[[4]]},
    raw[q0,s0,m0,n0]-cross[q0,s0,m0] raw[representative[q0[[{3,2,1,4}]]],s0,m0,n0]]],rows];
  tag=ToString[{l,p},InputForm];If[KeyExistsQ[out,tag],Throw[$Failed]];
  out[tag]=<|"Couplings"->vars,"Equations"->Expand[eqs/.x->x+gap-ref["delta_minus_x"]],"Prefactor"->spinBlockPrefactor[ref,gap]|>,{sector,sectors}];
 If[Length[out]==0,Throw[$Failed]];
 polys=spinMakeSDP[id,out];If[polys===$Failed,Throw[$Failed]];
 Join[polys,<|"Rows"->rows,"ExternalSpin"->j,"ExternalDimension"->delta,"Digits"->digits|>]
 ];If[ans===$Failed,Message[spinIdenticalSDP::input]];ans];
spinIdenticalSDP[___]:=(Message[spinIdenticalSDP::input];$Failed);
spinMajoranaRows[sdp_Association] /; Lookup[sdp,"ExternalSpin",None]===1/2 := Module[{keep,out=sdp},
 keep=Select[Range[Length[sdp["Rows"]]],Function[i,With[{q=sdp["Rows"][[i,1]],s=sdp["Rows"][[i,2]],m=sdp["Rows"][[i,3]],n=sdp["Rows"][[i,4]]},
  (q=={1/2,1/2,1/2,1/2} && ((s==1 && OddQ[m] && n>=1)||(s==-1 && EvenQ[m]))) ||
  (q=={1/2,-1/2,1/2,-1/2} && s==1 && OddQ[m]) || (q=={1/2,1/2,-1/2,-1/2} && s==1)]]];
 out["Identity"]=sdp["Identity"][[keep]];out["Rows"]=sdp["Rows"][[keep]];
 out["Sectors"]=Map[Join[#,<|"Matrices"->#["Matrices"][[keep]]|>]&,sdp["Sectors"]];Join[out,<|"RetainedRows"->keep|>]];
spinMajoranaRows[___]:=(Message[spinIdenticalSDP::input];$Failed);

spinIndependentRows[sdp_Association,x_Symbol,tol_] := Module[{id=sdp["Identity"],v,basis={},keep={},residual={},w,norm,r,out=sdp,blocks,scale},
 If[!NumericQ[tol] || !TrueQ[0<tol<1],Message[spinIndependentRows::input];Return[$Failed]];
 blocks=Values[sdp["Sectors"]];
 (* Pad each polynomial separately so coefficient columns have the same meaning. *)
 v=Transpose[Join[{id},Flatten[Table[With[{mat=b["Matrices"],degree=Max[0,Sequence@@Exponent[Flatten[b["Matrices"]],x]]},
  Flatten[Table[Table[Coefficient[mat[[All,a,c]],x,k],{k,0,degree}],{a,Length[mat[[1]]]},{c,Length[mat[[1]]]}],2]],{b,blocks}],1]]];
 If[!MatrixQ[v,NumericQ] || !FreeQ[v,_Complex|Indeterminate|_DirectedInfinity],Message[spinIndependentRows::input];Return[$Failed]];
 scale=Max[Norm/@v];
 Do[norm=Norm[v[[i]]];w=If[TrueQ[norm<=tol scale],ConstantArray[0,Length[v[[i]]]],v[[i]]/norm];
  Do[Do[w-= (u.w) u,{u,basis}],{2}];r=Norm[w];AppendTo[residual,r];
  If[TrueQ[r>tol],AppendTo[keep,i];AppendTo[basis,w/r]],{i,Length[v]}];
 out["Identity"]=id[[keep]];out["Sectors"]=Map[Function[b,Join[b,<|"Matrices"->b["Matrices"][[keep]]|>]],sdp["Sectors"]];
 If[KeyExistsQ[sdp,"Rows"],out["Rows"]=sdp["Rows"][[keep]]];
 Join[out,<|"RetainedRows"->keep,"RowResiduals"->residual,"RowTolerance"->tol|>]];
End[];
EndPackage[];
