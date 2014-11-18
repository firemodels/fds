#!/bin/bash

rm -f */*.err
./Run_FDS_Cases.sh -d -m 1

# after cases finish (kill ones that linger after 15 minutes or so)
# look for erros with:

# grep forrtl */*err
# grep ERROR */*err
