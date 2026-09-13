SetDirectory[DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]];
Get["spinning-identical.m"];
failed=False;
assert[label_,ok_]:=If[TrueQ[ok],Print[label,"=PASS"],Print[label,"=FAIL"];failed=True];
file=FileNameJoin[{$TemporaryDirectory,"autoboot-blocks-"<>CreateUUID[]<>".json"}];
(* A deliberately simple scalar table: every derivative numerator is 1+2 x.
   Synthetic data tests the interface, not conformal-block physics. *)
raw=<|"j_external"->{0,0,0,0},"j_internal"->0,"delta_minus_x"->0.5,"j_12"->0,"j_43"->0,
 "delta_12"->"0","delta_43"->"0","delta_1_plus_2"->"1.6","four_pt_struct"->{0,0,0,0},"four_pt_sign"->1,
 "order"->8,"lambda"->3,"kept_pole_order"->4,"pole_list_x"->{-1,0},"precision_actual"->300,
 "index_names"->{"three_pt_parity_index","j_120_index","j_430_index","coordinate_index","deriv_0","deriv_1","monomial_degree"},
 "index_values"->{<|"j_120"->{0},"j_430"->{0},"coordinates"->{"xt"}|>,<|"j_120"->{},"j_430"->{},"coordinates"->{"xt"}|>},
 "derivs"->{{{{Table[{"1e0","+2.0"},{m,0,3},{n,0,Floor[(3-m)/2]}]}}},{}}|>;
Export[file,raw,"RawJSON"];data=readSpinningBlocks[file,40];
assert["blocks decimal parser",AssociationQ[data] && data["delta_1_plus_2"]===8/5];
assert["blocks polynomial axes",Chop[spinBlockPolynomial[data,{0,0},{1,1},x]-(1+2 x)]===0];
assert["blocks derivative range",Quiet[spinBlockPolynomial[data,{0,0},{3,1},x]]===$Failed];
assert["shifted prefactor",spinBlockPrefactor[data,1]=={3-2 Sqrt[2],{-3/2,-1/2},3-2 Sqrt[2]}];
minus=data;minus["four_pt_sign"]=-1;
sdp=spinIdenticalSDP[{data,minus},{<|"Spin"->0,"Parity"->1,"Gap"->1|>},x,3];
assert["scalar derivative crossing",AssociationQ[sdp] && Length[sdp["Rows"]]==10 &&
 Chop[sdp["Sectors"]["{0, 1}"]["Matrices"][[3,1,1]]-(4+4 x)]===0];
assert["scalar identity derivative",Abs[sdp["Identity"][[3]]-N[-(32/5) 2^(8/5),40]]<10^-35];
assert["inconsistent table rejected",Quiet[spinIdenticalSDP[{data,Join[minus,<|"delta_1_plus_2"->2|>]},{<|"Spin"->0,"Parity"->1,"Gap"->1|>},x,3]]===$Failed];
assert["conserved external rejected",Quiet[spinIdenticalSDP[{Join[data,<|"j_external"->{1,1,1,1},"delta_1_plus_2"->4|>]},{<|"Spin"->0,"Parity"->1,"Gap"->1|>},x,3]]===$Failed];
assert["unitarity-pole gap rejected",Quiet[spinIdenticalSDP[{data,minus},{<|"Spin"->0,"Parity"->1,"Gap"->1/2|>},x,3]]===$Failed];
bad=raw;bad["delta_12"]="(executed=True;0)";Export[file,bad,"RawJSON"];
assert["untrusted decimal not evaluated",Quiet[readSpinningBlocks[file,40]]===$Failed && !ValueQ[executed]];
bad=raw;bad["precision_actual"]=64;Export[file,bad,"RawJSON"];
assert["insufficient block precision rejected",Quiet[readSpinningBlocks[file,40]]===$Failed];
DeleteFile[file];
Print["SPINNING_BLOCKS=",If[failed,"FAIL","PASS"]];Exit[If[failed,1,0]];
