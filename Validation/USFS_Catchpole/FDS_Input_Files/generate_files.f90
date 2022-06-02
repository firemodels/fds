program generate_files

real, dimension(360) :: sigma,beta,delta,U,M,R
character(10) :: test_name(360),endtime,windspeed,vegheight,moistfrac,svratio,masspervolume,vegdensity,heatofreaction,&
                 ignitionend,ignitionp1
character(80) :: command_string
real :: density,heat_of_reaction,ignition_end,ignition_p1
logical :: ex

open(10,file='Test_Matrix.csv',form='formatted')
read(10,*)

do i=1,354
   read(10,*) test_name(i),sigma(i),beta(i),delta(i),U(i),M(i),R(i)
   inquire(file=trim(test_name(i))//'.fds',exist=ex)
   if (ex) then
      write(0,*) i,trim(test_name(i))//'.fds exists'
   endif

   command_string = 'cp template.txt '//trim(test_name(i))//'.fds'
   call execute_command_line(command_string)

   command_string = "perl -pi -e 's/JOBID/"//trim(test_name(i))//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)

   write(endtime,'(f5.0)') min(999.,12./R(i))
   command_string = "perl -pi -e 's/endtime/"//trim(endtime)//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)

   write(windspeed,'(f4.2)') U(i)
   command_string = "perl -pi -e 's/windspeed/"//trim(windspeed)//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)

   write(vegheight,'(f5.3)') delta(i)
   command_string = "perl -pi -e 's/vegheight/"//trim(vegheight)//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)

   write(moistfrac,'(f5.3)') M(i)
   command_string = "perl -pi -e 's/moistfrac/"//trim(moistfrac)//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)

   write(svratio,'(f5.0)') sigma(i)
   command_string = "perl -pi -e 's/svratio/"//trim(svratio)//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)

   if (test_name(i)(1:2)=='MF') then
      density = 442.
      heat_of_reaction = 659.
      ignition_end = 40.
      ignition_p1 = ignition_end + 1
   elseif (test_name(i)(1:4)=='EXSC') then
      density = 398.
      heat_of_reaction = 711.
      ignition_end = 10.
      ignition_p1 = ignition_end + 1
   elseif (test_name(i)(1:4)=='PPMC') then
      density = 510.
      heat_of_reaction = 609.
      ignition_end = 20.
      ignition_p1 = ignition_end + 1
   elseif (test_name(i)(1:4)=='EXMC' .or. test_name(i)(1:2)=='EX') then
      density = 398.
      heat_of_reaction = 711.
      ignition_end = 10.
      ignition_p1 = ignition_end + 1
   endif

   write(masspervolume,'(f5.2)') beta(i)*density
   command_string = "perl -pi -e 's/masspervolume/"//trim(masspervolume)//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)

   write(vegdensity,'(f4.0)') density
   command_string = "perl -pi -e 's/vegdensity/"//trim(vegdensity)//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)

   write(heatofreaction,'(f4.0)') heat_of_reaction
   command_string = "perl -pi -e 's/heatofreaction/"//trim(heatofreaction)//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)

   write(ignitionend,'(f3.0)') ignition_end
   command_string = "perl -pi -e 's/ignitionend/"//trim(ignitionend)//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)

   write(ignitionp1,'(f3.0)') ignition_p1
   command_string = "perl -pi -e 's/ignitionp1/"//trim(ignitionp1)//"/g' "//trim(test_name(i))//".fds"
   call execute_command_line(command_string)
enddo

end program 
