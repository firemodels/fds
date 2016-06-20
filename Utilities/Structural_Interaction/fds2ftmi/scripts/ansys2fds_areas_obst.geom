/prep7																	! Shell file

nsle,S,corner															! Select just the corner nodes (exclude midside nodes)
*get,KCOUNT,kp,,count 													! Get total number of selected kps
*get,ACOUNT,area,,count 												! Get total number of selected areas
*dim,KARRAY,array,KCOUNT,3												! Create NCOUNT x 3 array
*dim,AARRAY,array,ACOUNT,4												! Create ECOUNT x 4 array
*dim,N2V,array,KCOUNT,1													! Create NCOUNT x 1 array
*dim,ELM,array,ACOUNT,1 												! Create ECOUNT x 1 array
*dim,OBST,array,4,3

*vget,N2V(1),kp,,klist													! Fill N2V with KP ID 

*do,I,1,KCOUNT
	*get,KARRAY(I,1),kp,N2V(I),loc,x 									! Fill first column with x-coord.
	*get,KARRAY(I,2),kp,N2V(I),loc,y 									! Fill second column with y-coord.
	*get,KARRAY(I,3),kp,N2V(I),loc,z 									! Fill third column with z-coord.
*enddo

*vget,ELM(1),area,,alist												! Fill ELM column with area ID

*do,I,1,ACOUNT
*get,A,area,ELM(I),loop,1,line,1 										! Get line 1.
*get,AARRAY(I,1),LINE,A,kp,1 											! Fill first column with kp 1 from line 1.
*get,AARRAY(I,2),LINE,A,kp,2 											! Fill second column with kp 2 from line 1.
*get,A,area,ELM(I),loop,1,line,3 										! Get line 2.
*get,AARRAY(I,3),LINE,A,kp,1 											! Fill third column with kp 1 from line 2.
*get,AARRAY(I,4),LINE,A,kp,2 											! Fill forth column with kp 2 from line 2.
*enddo

*do,I,1,ACOUNT															! Translate the kps numbers into the column number.
	*do,J,1,KCOUNT														! This is to avoid search in the N2V vector.
		*if,AARRAY(I,1),EQ,N2V(J),THEN									
		*set,AARRAY(I,1),J 												
		*endif															
		*if,AARRAY(I,2),EQ,N2V(J),THEN 									
		*set,AARRAY(I,2),J
		*endif
		*if,AARRAY(I,3),EQ,N2V(J),THEN
		*set,AARRAY(I,3),J
		*endif
		*if,AARRAY(I,4),EQ,N2V(J),THEN
		*set,AARRAY(I,4),J
		*endif
	*enddo	
*enddo

*cfopen,2fds_input,geo 													! Create file called “2fds_input.geo”

*do,I,1,ACOUNT
	*set,OBST(1,1),KARRAY(AARRAY(I,1),1)								! Saving the coordinates of KPs into matrix.
	*set,OBST(1,2),KARRAY(AARRAY(I,1),2)								! Then choose the max and min to write in OBST line.
	*set,OBST(1,3),KARRAY(AARRAY(I,1),3)
	*set,OBST(2,1),KARRAY(AARRAY(I,2),1)
	*set,OBST(2,2),KARRAY(AARRAY(I,2),2)
	*set,OBST(2,3),KARRAY(AARRAY(I,2),3)
	*set,OBST(3,1),KARRAY(AARRAY(I,3),1)
	*set,OBST(3,2),KARRAY(AARRAY(I,3),2)
	*set,OBST(3,3),KARRAY(AARRAY(I,3),3)	
	*set,OBST(4,1),KARRAY(AARRAY(I,4),1)
	*set,OBST(4,2),KARRAY(AARRAY(I,4),2)
	*set,OBST(4,3),KARRAY(AARRAY(I,4),3)
	*vwrite,min(OBST(1,1),OBST(2,1),OBST(3,1),OBST(4,1)),max(OBST(1,1),OBST(2,1),OBST(3,1),OBST(4,1)),min(OBST(1,2),OBST(2,2),OBST(3,2),OBST(4,2)),max(OBST(1,2),OBST(2,2),OBST(3,2),OBST(4,2)),min(OBST(1,3),OBST(2,3),OBST(3,3),OBST(4,3)),max(OBST(1,3),OBST(2,3),OBST(3,3),OBST(4,3))
&OBST XB=%8.3F,%8.3F,%8.3F,%8.3F,%8.3F,%8.3F, SURF_ID='STEEL', BNDF_OBST=.TRUE.  /
*enddo

finish