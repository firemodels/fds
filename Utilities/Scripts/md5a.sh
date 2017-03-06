#!/bin/bash
FILE=$1
md5 $FILE | awk -F" " '{print $4," ",substr($2,2,length($2)-2)}' 
