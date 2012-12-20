#!/bin/bash
base=$1
last=$2
underscore="_"
png=".png"

num2seq()
{
  num=$1
  if [[ $num -lt 10 ]] ; then
    seq="000$num"
  elif [[ $num -ge 10 && $num -lt 100 ]] ; then
    seq="00$num"
  elif [[ $num -ge 100 && $num -lt 999 ]] ; then
    seq="0$num"
  else
    seq=$num
  fi
}
echo "RENDERDIR"
echo " frames"
echo "LOADINIFILE"
echo " $base.ini"
index=0
while [[ $index -lt $last ]]; do
  num2seq $index
  file=$base$underscore$seq$png
  if ! [ -e $file ]; then
    echo LOADVOLSMOKEFRAME
    echo " -1 $index"
  fi
  let index=$index+1
done
