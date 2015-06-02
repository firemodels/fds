#!/bin/bash
authors=$(svn log -q $@ | grep -e '^r' | awk 'BEGIN { FS = "|" } ; { print $2 } $1 ~ /r([0-9]+000)/ { print "fetched revision " substr($1, 2) > "/dev/stderr" }' | sort | uniq)
for author in ${authors}; do
echo "${author} = NAME <user@domain>";
done
