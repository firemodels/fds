SAVE

/prep7
*get,HIGHNODE,node,,num,max												! Get highest node number in the model
*get,ETYPE,etyp,,num,max												! Get highest element type number
ETYPE=ETYPE+1
nsla,s,1																! Select nodes on exterior face (surfaces)
ET,5,SURF152															! Define Element Type, 5
TYPE,   5   
MAT,       1
REAL,   
ESYS,       0   
SECNUM, 
!*  
ESURF,0 																! Create dummy surface elements
!*  
esel,S,type,,5 															! Select only surface effect elements

*get,NCOUNT,node,,count 												! Get total number of selected nodes
*get,ECOUNT,elem,,count 												! Get total number of selected elements
*dim,NARRAY,array,NCOUNT,3												! Create NCOUNT x 3 array
*dim,EARRAY,array,ECOUNT,8												! Create ECOUNT x 8 array
*dim,N2V,array,NCOUNT,1													! Create NCOUNT x 1 array
*dim,ELM,array,ECOUNT,1 												! Create ECOUNT x 1 array

*vget,N2V(1),node,,nlist												! Fill N2V with node ID, VERTs ID is the vector line
*do,I,1,NCOUNT
	*get,NARRAY(I,1),node,N2V(I),loc,x 									! Fill first column with x-coord.
	*get,NARRAY(I,2),node,N2V(I),loc,y 									! Fill second column with y-coord.
	*get,NARRAY(I,3),node,N2V(I),loc,z 									! Fill third column with z-coord.
*enddo

*cfopen,nodes,dat 														! Create file called “nodes.dat”
*vwrite,NCOUNT,HIGHNODE															! Writing Nodes
%I %I
*vwrite,N2V(1),NARRAY(1,1),NARRAY(1,2),NARRAY(1,3)						! Write NARRAY to file
%I %8.3F %8.3F %8.3F 
*cfclose																! Close file called “nodes.dat”

*vget,ELM(1),elem,,elist												! Fill ELM column with element ID

*do,I,1,ECOUNT
*get,EARRAY(I,1),elem,ELM(I),node,1 									! Fill first column with node 1.
*get,EARRAY(I,2),elem,ELM(I),node,2 									! Fill second column with node 2.
*get,EARRAY(I,3),elem,ELM(I),node,3 									! Fill third column with node 3.
*get,EARRAY(I,4),elem,ELM(I),node,4 									! Fill fourth column with node 4.
*get,EARRAY(I,5),elem,ELM(I),node,5 									! Fill fifth column with node 5.
*get,EARRAY(I,6),elem,ELM(I),node,6 									! Fill sixth column with node 6.
*get,EARRAY(I,7),elem,ELM(I),node,7 									! Fill seventh column with node 7.
*get,EARRAY(I,8),elem,ELM(I),node,8 									! Fill eighth column with node 8.
*enddo

*cfopen,elements,dat 													! Create file called “elements.dat”
*vwrite,ECOUNT,0,0,ETYPE												! Writing Elements
%8I %8I %8I %8I
*vwrite,ELM(1),EARRAY(1,1),EARRAY(1,2),EARRAY(1,3),EARRAY(1,4),EARRAY(1,5),EARRAY(1,6),EARRAY(1,7),EARRAY(1,8)
%8I %8I %8I %8I %8I %8I %8I %8I %8I
*vwrite
END
*cfclose																! Close file called “elements.dat”

RESUME																	! Delete dummy surface effect elements, Element Type, etc...