(* Reader for blocks_3d's documented JSON output; no executable expressions
   are accepted in decimal strings. See arXiv:2011.01959 section 3.3. *)
Needs["SpinningBootstrap`","spinning.m"];
BeginPackage["SpinningBlocks`"];
readSpinningBlocks::usage="readSpinningBlocks[file,digits] reads blocks_3d JSON with at least digits decimal digits of working precision.";
spinBlockPolynomial::usage="spinBlockPolynomial[data,{j120,j430},{m,n},x] returns the xt derivative numerator polynomial. Imaginary fermion-exchange coefficients are restored.";
spinBlockPrefactor::usage="spinBlockPrefactor[data,gap] returns {constant,poles,base} after setting Delta=x+gap; uses the blocks_3d normalization a0=1,bj=1.";
readSpinningBlocks::input="Invalid or insufficient-precision blocks_3d data in `1`.";
spinBlockPolynomial::input="The requested basis labels, coordinates, or derivative are unavailable.";
Begin["`Private`"];
integerString[s_] := If[StringStartsQ[s,"-"],-FromDigits[StringDrop[s,1]],FromDigits[StringDelete[s,StartOfString~~"+"]]];
numberString[s_String] := Module[{parts,mant,exponent=0,sign=1,dot},
 If[!StringMatchQ[s,RegularExpression["[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?"]],Throw[$Failed,"blocks"]];
 parts=StringSplit[s,RegularExpression["[eE]"]];mant=First[parts];If[Length[parts]==2,exponent=integerString[Last[parts]]];
 If[StringStartsQ[mant,"-"],sign=-1;mant=StringDrop[mant,1]];mant=StringDelete[mant,StartOfString~~"+"];
 dot=StringPosition[mant,"."];If[dot=!={},exponent-=StringLength[mant]-dot[[1,1]]];
 sign FromDigits[StringDelete[mant,"."]] 10^exponent];
numberString[_] := Throw[$Failed,"blocks"];
readSpinningBlocks[file_String,digits_Integer:60] := Module[{data,ans},
 ans=Catch[
  data=Quiet[Check[Import[file,"RawJSON"],$Failed]];
  If[!AssociationQ[data] || !And@@(KeyExistsQ[data,#]&/@{"j_external","j_internal","delta_minus_x","j_12","j_43","four_pt_struct","four_pt_sign","delta_12","delta_43","delta_1_plus_2","pole_list_x","precision_actual","index_values","index_names","derivs","order","lambda","kept_pole_order"}),Throw[$Failed,"blocks"]];
  If[digits<20 || !IntegerQ[data["precision_actual"]] || data["precision_actual"] Log[10,2]<digits+5 ||
   data["index_names"]=!={"three_pt_parity_index","j_120_index","j_430_index","coordinate_index","deriv_0","deriv_1","monomial_degree"},Throw[$Failed,"blocks"]];
  Do[data[key]=Rationalize[data[key],0],{key,{"j_external","j_internal","delta_minus_x","j_12","j_43","four_pt_struct","pole_list_x","index_values"}}];
  If[!VectorQ[data["j_external"],(IntegerQ[2 #] && #>=0)&] || Length[data["j_external"]]!=4 ||
   !IntegerQ[2 data["j_internal"]] || data["j_internal"]<0 || !MemberQ[{-1,1},data["four_pt_sign"]] ||
   !VectorQ[data["pole_list_x"],NumericQ] || !AllTrue[Lookup[data,{"order","lambda","kept_pole_order"}],(IntegerQ[#] && #>=0)&] ||
   !ListQ[data["index_values"]] || Length[data["index_values"]]!=2 || !ListQ[data["derivs"]] || Length[data["derivs"]]!=2 ||
   !AllTrue[data["index_values"],Function[v,AssociationQ[v] && And@@(ListQ[Lookup[v,#,None]]&/@{"j_120","j_430","coordinates"})]],Throw[$Failed,"blocks"]];
  Do[data[key]=numberString[data[key]],{key,{"delta_12","delta_43","delta_1_plus_2"}}];
  data["derivs"]=Map[N[numberString[#],digits+5]&,data["derivs"],{7}];
  data["Digits"]=digits;data,
 "blocks"];
 If[ans===$Failed,Message[readSpinningBlocks::input,file]];ans];
readSpinningBlocks[___] := (Message[readSpinningBlocks::input,"arguments"];$Failed);
spinBlockPolynomial[data_Association,{a_,b_},{m_Integer,n_Integer},x_Symbol] := Module[{idx,p1,p2,coord,coeff},
 If[m<0 || n<0 || !MemberQ[Range[Abs[data["j_internal"]-data["j_12"]],data["j_internal"]+data["j_12"]],a] ||
  !MemberQ[Range[Abs[data["j_internal"]-data["j_43"]],data["j_internal"]+data["j_43"]],b],Message[spinBlockPolynomial::input];Return[$Failed]];
 idx=Select[Range[Length[data["index_values"]]],
  MemberQ[data["index_values"][[#]]["j_120"],a] && MemberQ[data["index_values"][[#]]["j_430"],b]&];
 If[idx=={},Return[0]];
 If[Length[idx]!=1,Message[spinBlockPolynomial::input];Return[$Failed]];idx=First[idx];
 p1=First@FirstPosition[data["index_values"][[idx]]["j_120"],a];
 p2=First@FirstPosition[data["index_values"][[idx]]["j_430"],b];
 coord=FirstPosition[data["index_values"][[idx]]["coordinates"],"xt"];
 If[MissingQ[coord],Message[spinBlockPolynomial::input];Return[$Failed]];
 coeff=Quiet[Check[data["derivs"][[idx,p1,p2,First[coord],m+1,n+1]],$Failed]];
 If[coeff===$Failed || !VectorQ[coeff,NumericQ],Message[spinBlockPolynomial::input];Return[$Failed]];
 If[IntegerQ[data["j_internal"]],1,I] Sum[coeff[[k]] x^(k-1),{k,Length[coeff]}]];
spinBlockPolynomial[___] := (Message[spinBlockPolynomial::input];$Failed);
spinBlockPrefactor[data_Association,gap_] := {(3-2 Sqrt[2])^gap,data["pole_list_x"]+data["delta_minus_x"]-gap,3-2 Sqrt[2]};
End[];
EndPackage[];
