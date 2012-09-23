clean_build=1

# Build SMV Verification Guide
pdflatex -interaction nonstopmode SMV_Verification_Guide &> SMV_Verification_Guide.err
bibtex SMV_Verification_Guide &> SMV_Verification_Guide.err
pdflatex -interaction nonstopmode SMV_Verification_Guide &> SMV_Verification_Guide.err
pdflatex -interaction nonstopmode SMV_Verification_Guide &> SMV_Verification_Guide.err

# Scan and report any errors in the LaTeX build process
if [[ `grep -E "Error:|Fatal error|! LaTeX Error:" -I SMV_Verification_Guide.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep -E "Error:|Fatal error|! LaTeX Error:" -I SMV_Verification_Guide.err
      clean_build=0
fi

# Check for LaTeX warnings (undefined references or duplicate labels)
if [[ `grep -E "undefined|multiply defined|multiply-defined" -I SMV_Verification_Guide.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX warnings detected:"
      grep -E "undefined|multiply defined|multiply-defined" -I SMV_Verification_Guide.err
      clean_build=0
fi

if [[ $clean_build == 0 ]]
   then
      :
   else
      echo "SMV Verification Guide built successfully!"
fi    
