#!/bin/sh

if [ "$#" != "1" ]; then
   echo "Usage: numgroup.sh <order>"
   exit 1
fi;

gap -b -q -x 800 << EOI
output := OutputTextFile("out.txt", false);
SetPrintFormattingStatus(output, false);
OutputLogTo(output);
NumberSmallGroups($1);
quit;
EOI
