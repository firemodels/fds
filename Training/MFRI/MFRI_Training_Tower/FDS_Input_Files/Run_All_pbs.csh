#!/bin/csh -f
foreach nn (01 02 03 04 05 06 07)
foreach type (min avg max)
set pbsfile=MFT_$nn\_$type.pbs
cp script.template $pbsfile
perl -p -i.bak -e s/XXXNUM/$nn/g $pbsfile
perl -p -i.bak -e s/YYYTYPE/$type/g $pbsfile
echo submitting $pbsfile
qsub $pbsfile
sleep 30
end
end

foreach nn (01)
foreach type (Lavg Savg)
set pbsfile=MFT_$nn\_$type.pbs
cp script.template $pbsfile
perl -p -i.bak -e s/XXXNUM/$nn/g $pbsfile
perl -p -i.bak -e s/YYYTYPE/$type/g $pbsfile
qsub $pbsfile
sleep 30
end
end

foreach nn (01)
foreach type (visavg)
set pbsfile=MFT_$nn\_$type.pbs
cp script.template $pbsfile
perl -p -i.bak -e s/XXXNUM/$nn/g $pbsfile
perl -p -i.bak -e s/YYYTYPE/$type/g $pbsfile
qsub $pbsfile
sleep 30
end
end
