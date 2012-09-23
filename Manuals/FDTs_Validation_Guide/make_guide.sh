clean_build=1

# Build FDTs Validation Guide
pdflatex -interaction nonstopmode FDTs_Validation_Guide &> FDTs_Validation_Guide.err
bibtex FDTs_Validation_Guide &> FDTs_Validation_Guide.err
pdflatex -interaction nonstopmode FDTs_Validation_Guide &> FDTs_Validation_Guide.err
pdflatex -interaction nonstopmode FDTs_Validation_Guide &> FDTs_Validation_Guide.err

# Scan and report any errors in the LaTeX build process
if [[ `grep -E "Error:|Fatal error|! LaTeX Error:" -I FDTs_Validation_Guide.err | grep -v "xpdf supports version 1.5"` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep -E "Error:|Fatal error|! LaTeX Error:" -I FDTs_Validation_Guide.err | grep -v "xpdf supports version 1.5"
      clean_build=0
fi

# Check for LaTeX warnings (undefined references or duplicate labels)
if [[ `grep -E "undefined|multiply defined|multiply-defined" -I FDTs_Validation_Guide.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX warnings detected:"
      grep -E "undefined|multiply defined|multiply-defined" -I FDTs_Validation_Guide.err
      clean_build=0
fi

if [[ $clean_build == 0 ]]
   then
      :
   else
      echo "FDTs Validation Guide built successfully!"
fi    
