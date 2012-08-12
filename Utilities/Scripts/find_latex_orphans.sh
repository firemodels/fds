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

FIGDIR=/FIGURES
GRAPHICS_FILES_F=`find FDS_Configuration_Management_Plan$FIGDIR FDS_Technical_Reference_Guide$FIGDIR FDS_User_Guide$FIGDIR FDS_Validation_Guide$FIGDIR FDS_Verification_Guide$FIGDIR SMV_Technical_Reference_Guide$FIGDIR SMV_User_Guide$FIGDIR SMV_Verification_Guide$FIGDIR -name *.png -o -name *.pdf -o -name *.eps -o -name *.jpg | xargs -n1 basename | sed 's/\..\{3\}$//'`

FIGDIR=/SCRIPT_FIGURES
GRAPHICS_FILES_SF=`find FDS_Configuration_Management_Plan$FIGDIR FDS_Technical_Reference_Guide$FIGDIR FDS_User_Guide$FIGDIR FDS_Validation_Guide$FIGDIR FDS_Verification_Guide$FIGDIR SMV_Technical_Reference_Guide$FIGDIR SMV_User_Guide$FIGDIR SMV_Verification_Guide$FIGDIR -name *.png -o -name *.pdf -o -name *.eps -o -name *.jpg | xargs -n1 basename | sed 's/\..\{3\}$//'`

# gf: have a loop outputting orphans for each GRPHICS_FILESxxx variable above

#  =================
#  = Compare lists =
#  =================

# See if all graphics are referenced in a LaTeX document
# If not, print "Unused graphics file (with filename)"
for i in $GRAPHICS_FILES_SF
do
    if [[ $REFERENCED_FILES == *"$i"* ]];
        then
            :
    else
        if [[ $i == *Guide* ]] || \
           [[ $i == *Plan* ]]
        then
            :
        else
            echo "Ununsed graphics file in SCRIPT_FIGURES: $i"
        fi
    fi
done

for i in $GRAPHICS_FILES_F
do
    if [[ $REFERENCED_FILES == *"$i"* ]];
        then
            :
    else
        if [[ $i == *Guide* ]] || \
           [[ $i == *Plan* ]]
        then
            :
        else
            echo "Ununsed graphics file in FIGURES: $i"
        fi
    fi
done
