clean_build=1

# Build FDS Verification Guide
pdflatex -interaction nonstopmode FDS_Verification_Guide &> FDS_Verification_Guide.err
bibtex FDS_Verification_Guide &> FDS_Verification_Guide.err
pdflatex -interaction nonstopmode FDS_Verification_Guide &> FDS_Verification_Guide.err
pdflatex -interaction nonstopmode FDS_Verification_Guide &> FDS_Verification_Guide.err

# Scan and report any errors in the LaTeX build process
if [[ `grep "! LaTeX Error:" -I FDS_Verification_Guide.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep "! LaTeX Error:" -I FDS_Verification_Guide.err
      clean_build=0
fi

# Check for LaTeX warnings (undefined references or duplicate labels)
if [[ `grep -E "undefined|multiply defined|multiply-defined" -I FDS_Verification_Guide.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX warnings detected:"
      grep -E "undefined|multiply defined|multiply-defined" -I FDS_Verification_Guide.err
      clean_build=0
fi

if [[ $clean_build == 0 ]]
   then
      :
   else
      echo "FDS Verification Guide built successfully!"
fi    
