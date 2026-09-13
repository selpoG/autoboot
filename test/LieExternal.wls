SetDirectory[DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]];
Get["lie-representations.m"];
failed=False;
assert[name_,p_] := If[TrueQ[p],Print[name,"=PASS"],Print[name,"=FAIL"];failed=True];
Do[
 {t,r,a,b,gapCartan,dimension,character,product}=sample;
 name=t<>ToString[r]<>ToString[a];
 assert[name<>" GAP Cartan convention",LieRepresentations`cartan[t,r]===Transpose[gapCartan]];
 assert[name<>" GAP dimension",LieRepresentations`repDimension[t,r,a]===dimension];
 {seconds,ch}=AbsoluteTiming[LieRepresentations`repCharacter[t,r,a]];
 dominant=Sort[KeyValueMap[List,KeySelect[ch,Min[#]>=0&]]];
 assert[name<>" GAP dominant weight multiplicities",dominant===Sort[Transpose[character]]];
 {productSeconds,pr}=AbsoluteTiming[LieRepresentations`repProduct[t,r,a,b]];
 assert[name<>" GAP tensor multiplicities",Sort[Tally[pr]]===Sort[Transpose[product]]];
 If[MemberQ[{{"A",3},{"G",2},{"E",6}},{t,r}],
  assert[name<>" character-convolution cross-check",Sort[pr]===Sort[LieRepresentations`Private`productByCharacters[t,r,a,b]]]];
 Print[name," CHARACTER_SECONDS=",seconds," PRODUCT_SECONDS=",productSeconds];
,{sample,Get["test/fixtures/LieGAP.m"]}];
Exit[If[failed,1,0]];
