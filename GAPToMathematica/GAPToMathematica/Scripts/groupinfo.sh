#!/bin/sh

if [ "$#" != "1" ]; then
   echo "Usage: groupinfo.sh <group>"
   exit 1
fi;

gap -b -q -x 800 << EOI
id := IdGroup($1);
g := SmallGroup(id[1], id[2]);
if not IsPcGroup(g) then g := SmallGroup(1, 1); fi;
ct := Irr(g);
gen := GeneratorsOfGroup(g);
irrs := List(ct, x -> IrreducibleRepresentationsDixon(g, x:unitary));
calcG := function()
Print("(id)=", IdGroup(g));
Print("(cgsize)=", List(ConjugacyClasses(g), x -> Size(x)));
Print("(chartab)=");
for x in ct do Print(List(x)); od;
Print("(gen)=", gen);
Print("(cggen)=", List(ConjugacyClasses(g), x -> Representative(x)));
Print("(irrep)=");
for irr in irrs do Print(List(gen, y -> Image(irr, y))); od;
return;
end;
output := OutputTextFile("out.txt", false);
SetPrintFormattingStatus(output, false);
OutputLogTo(output);
calcG();
quit;
EOI
