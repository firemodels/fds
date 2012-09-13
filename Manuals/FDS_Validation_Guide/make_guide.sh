clean_build=1

# Build FDS Validation Guide
pdflatex -interaction nonstopmode FDS_Validation_Guide &> FDS_Validation_Guide.err
bibtex FDS_Validation_Guide &> FDS_Validation_Guide.err
pdflatex -interaction nonstopmode FDS_Validation_Guide &> FDS_Validation_Guide.err
pdflatex -interaction nonstopmode FDS_Validation_Guide &> FDS_Validation_Guide.err

# Scan and report any errors in the LaTeX build process
if [[ `grep -E "Error: pdflatex|Fatal error|! LaTeX Error:" -I FDS_Validation_Guide.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep -E "Error: pdflatex|Fatal error|! LaTeX Error:" -I FDS_Validation_Guide.err
      clean_build=0
fi

# Check for LaTeX warnings (undefined references or duplicate labels)
if [[ `grep -E "undefined|multiply defined|multiply-defined" -I FDS_Validation_Guide.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX warnings detected:"
      grep -E "undefined|multiply defined|multiply-defined" -I FDS_Validation_Guide.err
      clean_build=0
fi

if [[ $clean_build == 0 ]]
   then
      :
   else
      echo "FDS Validation Guide built successfully!"
fi    
