"""csv2xlsx cut and paste

This script opens a series of comma-delimited .csv files and copies portions of columns into specified
worksheets of a specified Excel (.xlsx) file. It saves the Excel file with a new name. The names of the
new and old Excel files are currently hardwired in the script below.

The names of the csv files are listed in the file called 'csv2xlxs_parameters.csv'. This file contains the full
path to the csv files, followed by the row/column pair where the data starts, the length of the column
snippet, the name of the worksheet in the Excel file, and the row/column pair where the snippet is to be
pasted.

Note that the module used to save the new Excel file does not support wmf (Windows MetaFile) images. Any
wmf images in the original Excel file will be dropped from the new Excel file.

This script requires the Excel parsing module 'openpyxl' and the csv reader/writer 'csv'.
"""

import csv
from openpyxl import load_workbook

wb = load_workbook('LNG_database_validation_comparison_draft.xlsx')

param_file = open('csv2xlsx_parameters.csv','r')
param_reader = csv.reader(param_file,delimiter=',')
params = list(param_reader)

for i in range(0,len(params)):

    orig_row = int(params[i][1]) - 1
    orig_col = int(params[i][2]) - 1
    csv_file = open(params[i][0],'r')
    data = list(csv.reader(csv_file,delimiter=','))
    new_row = int(params[i][5])
    new_col = int(params[i][6])
    length  = int(params[i][3])
    offset  = new_row - orig_row
 
    ws = wb[params[i][4]]

    for row in range(new_row,new_row+length):
        for col in range(new_col,new_col+1):
            _ = ws.cell(column=col, row=row, value=round(float(data[row-offset][orig_col]),2))

    csv_file.close()

wb.save('LNG_database_validation_comparison_modified.xlsx')
