@echo off

:: Name and location of output .tex file with validation SVN statistics
set OUTPUT_TEX_FILE=..\..\Manuals\FDS_Validation_Guide\SCRIPT_FIGURES\ScatterPlots\validation_svn_stats.tex

copy stats_table_template.tex %OUTPUT_TEX_FILE%

