clean_build=1

# Build FDS Configuration Management Plan
pdflatex -interaction nonstopmode FDS_Configuration_Management_Plan &> FDS_Configuration_Management_Plan.err
bibtex FDS_Configuration_Management_Plan &> FDS_Configuration_Management_Plan.err
pdflatex -interaction nonstopmode FDS_Configuration_Management_Plan &> FDS_Configuration_Management_Plan.err
pdflatex -interaction nonstopmode FDS_Configuration_Management_Plan &> FDS_Configuration_Management_Plan.err

# Scan and report any errors in the LaTeX build process
if [[ `grep -E "Error: pdflatex|Fatal error|! LaTeX Error:" -I FDS_Configuration_Management_Plan.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep -E "Error: pdflatex|Fatal error|! LaTeX Error:" -I FDS_Configuration_Management_Plan.err
      clean_build=0
fi

# Check for LaTeX warnings (undefined references or duplicate labels)
if [[ `grep -E "undefined|multiply defined|multiply-defined" -I FDS_Configuration_Management_Plan.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX warnings detected:"
      grep -E "undefined|multiply defined|multiply-defined" -I FDS_Configuration_Management_Plan.err
      clean_build=0
fi

if [[ $clean_build == 0 ]]
   then
      :
   else
      echo "FDS Configuration Management Plan built successfully!"
fi    
