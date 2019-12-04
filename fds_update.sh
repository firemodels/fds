#!/bin/bash
git remote update
git checkout master
git diff firemodels/master
git merge firemodels/master
git push origin master