#! /bin/sh -e
TMPFILE=/tmp/to_unix$$

tr -d '\015\032' < "$1" > $TMPFILE
mv -f $TMPFILE "$1"
