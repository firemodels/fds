#!/bin/bash
cd cor
git remote add firemodels git@github.com:firemodels/cor.git
git remote update

cd ../exp
git remote add firemodels git@github.com:firemodels/exp.git
git remote update

cd ../out
git remote add firemodels git@github.com:firemodels/out.git
git remote update

cd ../smv
git remote add firemodels git@github.com:firemodels/smv.git
git remote update

cd ../fds
git remote add firemodels git@github.com:firemodels/fds.git
git remote update

cd ../radcal
git remote add firemodels git@github.com:firemodels/radcal.git
git remote update
