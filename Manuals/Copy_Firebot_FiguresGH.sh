#!/bin/bash
DOWNLOADFIGURES ()
{
  FILE=$1_figures
  echo ***Downloading $FILE.tar.gz
  gh release download FDS_TEST -p $FILE.tar.gz -R github.com/firemodels/test_bundles --clobber
  if [ ! -e $FILE.tar.gz ]; then
    echo ***error: $FILE.tar.gz failed to download
    echo command: gh release download FDS_TEST -p $FILE.tar.gz -R github.com/firemodels/test_bundles --clobber
  else
    echo "   successful"
  fi
}
DOWNLOADFIGURES FDS_TG
DOWNLOADFIGURES FDS_UG
DOWNLOADFIGURES FDS_VALG
DOWNLOADFIGURES FDS_VERG

BASEDIR=`pwd`

rm -rf FIGS
mkdir FIGS
mkdir FIGS/FDS_TG
mkdir FIGS/FDS_UG
mkdir FIGS/FDS_VALG
mkdir FIGS/FDS_VERG

FBTG=$BASEDIR/FIGS/FDS_TG
FBUG=$BASEDIR/FIGS/FDS_UG
FBVG=$BASEDIR/FIGS/FDS_VERG
FBVAL=$BASEDIR/FIGS/FDS_VALG

# Copy Tech Guide Figures
if [ -e $BASEDIR/FDS_TG_figures.tar.gz ]; then
  cd $FBTG
  tar xf $BASEDIR/FDS_TG_figures.tar.gz
  cp * $BASEDIR/FDS_Technical_Reference_Guide/SCRIPT_FIGURES/
  echo FDS Technical Guide figures copied to FDS_Technical_Reference_Guide/SCRIPT_FIGURES
fi

# Copy User's Guide Figures
if [ -e $BASEDIR/FDS_UG_figures.tar.gz ]; then
  cd $FBUG
  tar xf $BASEDIR/FDS_UG_figures.tar.gz
  cp * $BASEDIR/FDS_User_Guide/SCRIPT_FIGURES/
  echo FDS Users Guide figures copied to FDS_User_Reference_Guide/SCRIPT_FIGURES
fi

# Copy Verification Guide Figures
if [ -e $BASEDIR/FDS_VERG_figures.tar.gz ]; then
  cd $FBVG
  tar xf $BASEDIR/FDS_VERG_figures.tar.gz
  cp *.pdf $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
  cp *.png $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
  cp *.tex $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/.
  cp Scatterplots/*.tex $BASEDIR/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/.
  echo FDS Verification Guide figures copied to FDS_Verification_Guide_Reference_Guide/SCRIPT_FIGURES
fi

# Copy Validation Guide Figures
if [ -e $BASEDIR/FDS_VALG_figures.tar.gz ]; then
  cd $FBVAL
  tar xf $BASEDIR/FDS_VALG_figures.tar.gz
  cp -R * $BASEDIR/FDS_Validation_Guide/SCRIPT_FIGURES
  echo FDS Validation Guide Figures copied to FDS_Validation_Guide_Reference_Guide/SCRIPT_FIGURES
fi

