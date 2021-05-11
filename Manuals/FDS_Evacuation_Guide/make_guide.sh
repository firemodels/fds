#!/bin/bash


# Build FDS+Evac Guide

gitrevision=`git describe --long --dirty`
echo "\\newcommand{\\gitrevision}{$gitrevision}" > gitrevision.tex

pdflatex -interaction nonstopmode FDS+Evac_Guide > FDS+Evac_Guide.err
bibtex FDS+Evac_Guide &> FDS+Evac_Guide.err
pdflatex -interaction nonstopmode FDS+Evac_Guide > FDS+Evac_Guide.err
pdflatex -interaction nonstopmode FDS+Evac_Guide > FDS+Evac_Guide.err
pdflatex -interaction nonstopmode FDS+Evac_Guide > FDS+Evac_Guide.err
