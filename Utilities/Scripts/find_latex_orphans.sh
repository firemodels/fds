#!/bin/bash

# This script compiles all LaTeX guides, then checks for "orphaned" files.
# These are files that are not referenced in the LaTeX guides.

SVNROOT=`pwd`/../..

cd $SVNROOT/Manuals

#  =============================================================================
#  = Make all guides with the -recorder option to output external dependencies =
#  =============================================================================

cd FDS_Configuration_Management_Plan
echo 'Building FDS Configuration Management Plan'
pdflatex -interaction nonstopmode FDS_Configuration_Management_Plan &> /dev/null
bibtex FDS_Configuration_Management_Plan &> /dev/null
pdflatex -interaction nonstopmode FDS_Configuration_Management_Plan &> /dev/null
pdflatex -interaction nonstopmode -recorder FDS_Configuration_Management_Plan &> /dev/null
if [[ -e FDS_Configuration_Management_Plan.pdf ]];
then
    :
else
    echo 'Warning: Could not create FDS Configuration Management Plan'
fi
cd ..

cd FDS_Technical_Reference_Guide
echo 'Building FDS Technical Reference Guide'
pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide &> /dev/null
bibtex FDS_Technical_Reference_Guide &> /dev/null
pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide &> /dev/null
pdflatex -interaction nonstopmode -recorder FDS_Technical_Reference_Guide &> /dev/null
if [[ -e FDS_Technical_Reference_Guide.pdf ]];
then
    :
else
    echo 'Warning: Could not create FDS Technical Reference Guide'
fi
cd ..

cd FDS_User_Guide
echo 'Building FDS User Guide'
pdflatex -interaction nonstopmode FDS_User_Guide &> /dev/null
bibtex FDS_User_Guide &> /dev/null
pdflatex -interaction nonstopmode FDS_User_Guide &> /dev/null
pdflatex -interaction nonstopmode -recorder FDS_User_Guide &> /dev/null
if [[ -e FDS_User_Guide.pdf ]];
then
    :
else
    echo 'Warning: Could not create FDS User Guide'
fi
cd ..

cd FDS_Validation_Guide
echo 'Building FDS Validation Guide'
pdflatex -interaction nonstopmode FDS_Validation_Guide &> /dev/null
bibtex FDS_Validation_Guide &> /dev/null
pdflatex -interaction nonstopmode FDS_Validation_Guide &> /dev/null
pdflatex -interaction nonstopmode -recorder FDS_Validation_Guide &> /dev/null
if [[ -e FDS_Validation_Guide.pdf ]];
then
    :
else
    echo 'Warning: Could not create FDS Validation Guide'
fi
cd ..

cd FDS_Verification_Guide
echo 'Building FDS Verification Guide'
pdflatex -interaction nonstopmode FDS_Verification_Guide &> /dev/null
bibtex FDS_Verification_Guide &> /dev/null
pdflatex -interaction nonstopmode FDS_Verification_Guide &> /dev/null
pdflatex -interaction nonstopmode -recorder FDS_Verification_Guide &> /dev/null
if [[ -e FDS_Verification_Guide.pdf ]];
then
    :
else
    echo 'Warning: Could not create FDS Verification Guide'
fi
cd ..

cd SMV_Technical_Reference_Guide
echo 'Building SMV Technical Reference Guide'
pdflatex -interaction nonstopmode SMV_Technical_Reference_Guide &> /dev/null
bibtex SMV_Technical_Reference_Guide &> /dev/null
pdflatex -interaction nonstopmode SMV_Technical_Reference_Guide &> /dev/null
pdflatex -interaction nonstopmode -recorder SMV_Technical_Reference_Guide &> /dev/null
if [[ -e SMV_Technical_Reference_Guide.pdf ]];
then
    :
else
    echo 'Warning: Could not create SMV Technical Reference Guide'
fi
cd ..

cd SMV_User_Guide
echo 'Building SMV User Guide'
pdflatex -interaction nonstopmode SMV_User_Guide &> /dev/null
bibtex SMV_User_Guide &> /dev/null
pdflatex -interaction nonstopmode SMV_User_Guide &> /dev/null
pdflatex -interaction nonstopmode -recorder SMV_User_Guide &> /dev/null
if [[ -e SMV_User_Guide.pdf ]];
then
    :
else
    echo 'Warning: Could not create SMV User Guide'
fi
cd ..

cd SMV_Verification_Guide
echo 'Building SMV Verification Guide'
pdflatex -interaction nonstopmode SMV_Verification_Guide &> /dev/null
bibtex SMV_Verification_Guide &> /dev/null
pdflatex -interaction nonstopmode SMV_Verification_Guide &> /dev/null
pdflatex -interaction nonstopmode -recorder SMV_Verification_Guide &> /dev/null
if [[ -e SMV_Verification_Guide.pdf ]];
then
    :
else
    echo 'Warning: Could not create SMV Verification Guide'
fi
cd ..

cd $SVNROOT/Manuals

#  ===========================
#  = Compile lists of images =
#  ===========================

# Compile list of png and pdf images referenced in LaTeX documents
REFERENCED_FILES=`grep -h INPUT */*.fls | grep -E 'pdf|png|eps|jpg' | cut -f2 -d' ' | xargs -n1 basename | sed 's/\..\{3\}$//'`

# Compile list of png and pdf images in Manuals directories
GRAPHICS_FILES=`find . -name *.png -o -name *.pdf -o -name *.eps -o -name *.jpg | xargs -n1 basename | sed 's/\..\{3\}$//'`

#  =================
#  = Compare lists =
#  =================

# See if all graphics are referenced in a LaTeX document
# If not, print "Unused graphics file (with filename)"
for i in $GRAPHICS_FILES
do
    if [[ $REFERENCED_FILES == *"$i"* ]];
        then
            :
    else
        if [[ $i == *Guide* ]];
        then
            :
        else
            echo "Ununsed graphics file: $i"
        fi
    fi
done

