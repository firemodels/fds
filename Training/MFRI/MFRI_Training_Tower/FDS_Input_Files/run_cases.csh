#!/bin/csh -f
foreach nn (01)
foreach type (Lavg Savg)
set pbsfile=MFT_$nn\_$type.pbs
cp script.template $pbsfile
perl -p -i.bak -e s/XXXNUM/$nn/g $pbsfile
perl -p -i.bak -e s/YYYTYPE/$type/g $pbsfile
qsub $pbsfile
end
end
