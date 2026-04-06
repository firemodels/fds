## Model Evaluation Protocol (MEP) for Fire Models Involving Fuels at Liquefied Natural Gas (LNG) Facilities

Under the sponsorship of the Pipeline and Hazardous Materials Safety Administration (PHMSA) of the U.S. Department of Transportation, Sandia National Laboratories has compiled a collection of experimental datasets for use in evaluating models that predict the consequences of accidental fires at LNG facilities. Details can be found in the report written by Anay Luketa, _Model Evaluation Protocol for Fire Models Involving Fuels at Liquefied Natural Gas Facilities (Version 2)_, Sandia Report SAND2026-19309, April 2026. 

The experiments included in this database are also included in the FDS validation suite and the results are reported in the FDS Validation Guide. In addition, Sandia has created an Excel spreadsheet with the basic parameters of each experiment, which is included in this directory under the name `MEPFires_database_spreadsheet_v2.xlsx`. Modelers are expected to report their results in the various sheets that are included in the file.

This directory contains a Python script (`csv2xlsx.py`) that reads the comma-delimited (`.csv`) FDS output files that reside in the `firemodels/out` repository, and transfers snippets of data into the Excel file. A file called `csv2xlsx_parameters.csv` contains a list of the data snippets to be transferred from the FDS output files to the Excel file. Note that the original Excel file is not changed after running the script; but rather a copy is created called `MEPFires_database_spreadsheet_v2_filled.xlsx` which contains the model results.

The columns of `csv2xlsx_parameters.csv` are as follows:

A. The FDS output file

B. The starting row of the FDS data

C. The starting column of the FDS data

D. The number of values, with a negative sign added if the data is to be read across rather than down

E. The name of the tab in the Excel file where the data is to be copied

F. The starting row of the data

G. The starting column of the data

