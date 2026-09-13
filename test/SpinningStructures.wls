SetDirectory[DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]];
Get["spinning.m"];
failed=False;
assert[label_,ok_]:=If[TrueQ[ok],Print[label,"=PASS"],Print[label,"=FAIL"];failed=True];
assert["invalid spins rejected",Quiet[spinQBasis[{0,0,1/3}]]===$Failed && Quiet[spinQBasis[{0,0,1/2}]]===$Failed];
Do[
 q=spinQBasis[js];b=spinSO3Basis[js];m=spinSO3ToQ[js];
 assert["q/SO3 dimensions",Length[q]==Length[b] && Dimensions[m]=={Length[q],Length[q]}];
 assert["q/SO3 invertibility",RootReduce[m.spinQToSO3[js]]==IdentityMatrix[Length[q]]];
 Do[r=spinThreePointBasis[js,p];assert["SO3 reconstruction",r["QCoefficients"]=={} || RootReduce[r["SO3Coefficients"].m-r["QCoefficients"]]==ConstantArray[0,Dimensions[r["QCoefficients"]]]],{p,{-1,1}}];
,{js,{{0,0,4},{1/2,1/2,0},{1/2,1/2,3},{1/2,1,3/2},{1,1,2},{2,2,4}}}];
(* arXiv:2011.01959 eqs. (5.6)-(5.8): identical Majorana selection rules. *)
Do[
 counts=Length[spinThreePointBasis[{1/2,1/2,j},#,-1]["QCoefficients"]]&/@{1,-1};
 assert["Majorana allowed structures",counts==If[j==0,{1,1},If[EvenQ[j],{2,1},{0,1}]]];
,{j,0,6}];
(* Published conversion (5.36), including normalization and sign. *)
Do[
 js={1/2,1/2,j};q=spinQBasis[js];v=(Boole[#=={1/2,1/2,-1}]+Boole[#=={-1/2,-1/2,1}]&/@q);
 result=RootReduce[v.spinQToSO3[js]];
 expected=Table[If[b=={1,j},(-1)^(1-j) Sqrt[2/Binomial[2 j,j-1]],0],{b,spinSO3Basis[js]}];
 assert["published Majorana basis conversion",result==expected];
,{j,1,5}];
Do[
 c=spinCrossingMatrix[js];d=spinCrossingMatrix[js[[{3,2,1,4}]]];
 assert["crossing involution",c["Matrix"].c["Matrix"]==IdentityMatrix[Length[c["Structures"]],SparseArray]];
 eq=spinCrossingEquations[js,0,s,t,z,zb];
 assert["all four-point components",Length[eq]==Times@@(2 js+1)];
,{js,{{0,0,0,0},{1/2,1/2,1/2,1/2},{0,1/2,0,1/2},{1,1,1,1},{2,0,2,0}}}];
assert["fermion crossing sign included",spinCrossingEquations[{1/2,1/2,1/2,1/2},1,s,s,z,zb][[-1]]==s[{1/2,1/2,1/2,1/2},z,zb]-s[{1/2,1/2,1/2,1/2},1-z,1-zb]];
assert["scalar identity",spinIdentityBlock[{0,0,0,0},delta,{0,0,0,0},z,zb]==z^-delta zb^-delta];
assert["fermion identity phase",spinIdentityBlock[{1/2,1/2,1/2,1/2},delta,{1/2,1/2,1/2,1/2},z,zb]==-z^(-delta+1/2) zb^(-delta-1/2)];
model=spinMakeSDP[{1,2},<|"spin2"-><|"Couplings"->{lambda1,lambda2},"Equations"->{lambda1^2 f+3 lambda1 lambda2 g+lambda2^2 h,lambda1^2+lambda2^2},"Prefactor"->{1,{-1},1/2}|>|>];
assert["off-diagonal factor",model["Sectors"]["spin2"]["Matrices"][[1,1,2]]==3 g/2];
assert["nonquadratic rejection",Quiet[spinMakeSDP[{1},<|"bad"-><|"Couplings"->{lambda1},"Equations"->{lambda1^3}|>|>]]===$Failed];
file=FileNameJoin[{$TemporaryDirectory,"autoboot-spinning-"<>CreateUUID[]<>".json"}];
sdp=model/.{f->1+x,g->x,h->2+x^2};
assert["SDPB PMP export",spinWritePMP[file,sdp,x,40]===file];
json=Import[file,"RawJSON"];
assert["PMP axes",Length[json["PositiveMatrixWithPrefactorArray"][[1]]["polynomials"]]==2 && Length[json["PositiveMatrixWithPrefactorArray"][[1]]["polynomials"][[1,1]]]==2];
DeleteFile[file];
assert["unevaluated blocks rejected",Quiet[spinWritePMP[file,model,x,40]]===$Failed && !FileExistsQ[file]];
Print["SPINNING_STRUCTURES=",If[failed,"FAIL","PASS"]];Exit[If[failed,1,0]];
