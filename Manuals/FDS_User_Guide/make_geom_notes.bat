@echo off
set paper=geom_notes

echo pass 1
pdflatex -interaction nonstopmode %paper% > %paper%.err
bibtex %paper% > %paper%.err
echo pass 2
pdflatex -interaction nonstopmode %paper% > %paper%.err
echo pass 3
pdflatex -interaction nonstopmode %paper% > %paper%.err

find "! LaTeX Error:" %paper%.err
find "Fatal error" %paper%.err
find "Error:" %paper%.err
find "LaTeX Warning:" %paper%.err

echo %paper% build complete
pause

