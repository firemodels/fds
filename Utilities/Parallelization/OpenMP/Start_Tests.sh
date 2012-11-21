#! /bin/bash

############### FUNCTIONS START ####################

function select_fds_version {
  FDSVERSION=$(zenity \
        --list \
        --width=1000 \
        --height=500 \
        --checklist \
        --column "" \
        --column "Filename" \
        --column "Description" \
                FALSE "fds5_intel_linux_32_db" "FDS5 non OpenMP version without optimization" \
                FALSE "fds5_openmp_intel_linux_32" "FDS5 OpenMP version without ThreadChecker settings (Release Version)" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck" "FDS5 OpenMP version (complete) with Intel ThreadChecker settings" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_relevant" "FDS5 OpenMP version (only OpenMP relevant files) with Intel ThreadChecker settings" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_pois" "POIS.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_turb" "TURB.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_fire" "FIRE.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_wall" "WALL.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_velo" "VELO.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_divg" "DIVG.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_mass" "MASS.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_radi" "RADI.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_pres" "PRES.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_part" "PART.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" \
                FALSE "fds5_openmp_intel_linux_32_db_TCheck_dump" "DUMP.F90 - FDS5 OpenMP version with Intel ThreadChecker settings for only this file" )
# value after TRUE/FALSE is submitted to FDSVERSION (fds5_openmp_intel_linux_32_db....)
}

function copy_fds_source {
#Copies FDS-Source directory and makefiles in one directory and change to this directory
  cp -r ../../../FDS_Source ./
  cp ../../../FDS_Compilation/makefile FDS_Source/makefile
  cp makefile_openmp_threadcheck FDS_Source/makefile_openmp_threadcheck
  cd FDS_Source
}

function copy_compiled_and_remove_source {
#Copies the compiled file to FDS_OpenMP-Versions directory and deletes the Source directory. Changing to base directory.
  cp "$FDS_COMPILED_FILE" ../FDS_OpenMP-Versions/"$FDS_COMPILED_FILE"
  cd ..
  if [ "$CURRENT_DIRECTORY" = "$PWD" ]; then
    rm -rf FDS_Source
  fi
}


function compile_fds_version {

  for FDS_COMPILED_FILE in $FDSVERSION
    do

      case "$FDS_COMPILED_FILE" in
        fds5_intel_linux_32_db)
          copy_fds_source
          make -f makefile intel_linux_32_db
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32)
          copy_fds_source
          make -f makefile openmp_intel_linux_32
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_with_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_relevant)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_relevant_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_pois)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_pois_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_turb)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_turb_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_fire)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_fire_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_wall)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_wall_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_velo)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_velo_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_divg)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_divg_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_mass)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_mass_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_radi)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_radi_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_pres)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_pres_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_part)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_part_TCheck
          copy_compiled_and_remove_source
          ;;
        fds5_openmp_intel_linux_32_db_TCheck_dump)
          copy_fds_source
          make -f makefile_openmp_threadcheck all_dump_TCheck
          copy_compiled_and_remove_source
          ;;
        *)
          zenity --info --text "No FDS-Version selected"
          ;;
     esac
   done
}


function select_calculation_filename {
### select files which should be calculated by FDS
  FILENAME=$(zenity \
        --list \
        --width=1000 \
        --height=500 \
        --checklist \
        --column "" \
        --column "Filename" \
        --column "Description" \
                FALSE "Test001" "Test1-File" \
                FALSE "Test002" "Test2-File (not existing)" \
                FALSE "Test003" "Test3-File (not existing)" )
# Structure: FALSE "Folder_Name = File_Name without .fds" "Description of the file"
}


function select_delete_mode {
### ask if files (smv, out, sf, bf, s3d...) should be deleted after calculation
DELETE_FILES=$(zenity \
        --list \
        --width=1000 \
        --height=500 \
        --radiolist \
        --column "" \
        --column "YES/NO" \
        --column "Description" \
                FALSE "YES" "Delete files after run" \
                TRUE "NO" "Do not delete files after run" )
}




function run_calculation {
  ### start calculation process for selected files, $NAME is the name of .fds files which should be calculated
  for NAME in $FILENAME
    do

      case "$NAME" in
        Test001)
          run_fds  
          ;;
        *)
          zenity --info --text "Other Value selected" 
          ;;

      esac
    done
}



function do_calculation {

  cp ./FDS_OpenMP-Versions/"$FDS_FILENAME" ./Test_Cases/"$NAME"/"$FDS_FILENAME"
  cd Test_Cases/"$NAME"/
  ./"$FDS_FILENAME" "$NAME".fds >& Screen_Output___"$FDS_FILENAME"___"$NAME".txt
  tcheck_cl -w 200 threadchecker.thr >& Screen_TCheck___"$FDS_FILENAME"___"$NAME".txt
  if [ "$DELETE_FILES" = "YES" ]; then 
    if [ "$CURRENT_DIRECTORY"/Test_Cases/"$NAME" = "$PWD" ]; then
      rm -f "$FDS_FILENAME" *.csv *.out *.smv *.sz *.q *.sf *.bf *.end *.s3d *.thr *.txf fds5*.*
    fi
  fi
  cd ..
  cd ..
}

function run_fds {

  for FDS_FILENAME in $FDSVERSION
  do

    case "$FDS_FILENAME" in
      fds5_intel_linux_32_db)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_relevant)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_pois)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_turb)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_fire)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_wall)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_velo)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_divg)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_mass)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_radi)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_pres)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_part)
        do_calculation
        ;;
      fds5_openmp_intel_linux_32_db_TCheck_dump)
        do_calculation
        ;;
    esac
  done
}




################## END FUNCTIONS #########################

###### START SCRIPT ######

zenity --info --text "First, select the FDS5 versions you want to test. Second, select the .fds-files for testing with the selected FDS5-Version."

DELETE_FILES="NO"
CURRENT_DIRECTORY=$PWD

### safe original file seperator
OIFS=$IFS
### new file seperator is "|"
IFS="|"

select_fds_version
compile_fds_version
select_calculation_filename
select_delete_mode
run_calculation

zenity --info --text "Calculations are finished"

### safe original file separator back
IFS=$OIFS 

