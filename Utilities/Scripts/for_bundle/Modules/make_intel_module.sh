#!/bin/bash
compilervars=$1
modulefile=$2
echo "#%Module" > $modulefile
perl env2 -from bash -to modulecmd "$compilervars intel64" >> $modulefile
