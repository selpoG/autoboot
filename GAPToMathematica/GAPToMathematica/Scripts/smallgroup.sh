#!/bin/sh

if [ "$#" != "2" ]; then
   echo "Usage: smallgroup.sh <order> <id>"
   exit 1
fi;

gap -b -q -x 800 << EOI
g := SmallGroup($1, $2);
if not IsPcGroup(g) then g := SmallGroup(1, 1); fi;
ct := Irr(g);
gen := GeneratorsOfGroup(g);
irrs := List(ct, x -> IrreducibleRepresentationsDixon(g, x));
calcG := function()
local x, irr;
Print("(id)=", IdGroup(g));
Print("(cgsize)=", List(ConjugacyClasses(g), x -> Size(x)));
Print("(chartab)=");
for x in ct do Print(List(x)); od;
Print("(gen)=", gen);
Print("(cggen)=", List(ConjugacyClasses(g), x -> Representative(x)));
Print("(irrep)=");
for irr in irrs do Print(List(gen, y -> Image(irr, y))); od;
Print("(elements)=", Elements(g));
return;
end;
output := OutputTextFile("out.txt", false);
SetPrintFormattingStatus(output, false);
OutputLogTo(output);
calcG();
quit;
EOI
