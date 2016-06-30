/prep7																	! Tetrahedron file

*get,HIGHNODE,node,,num,max												! Get highest node number in the model
*get,ETYPE,etyp,,num,max												! Get highest element type number
ETYPE=ETYPE+1

nsle,S,corner															! Select just the corner nodes (exclude midside nodes)
*get,NCOUNT,node,,count 												! Get total number of selected nodes
*dim,NARRAY,array,NCOUNT,3												! Create NCOUNT x 3 array
*dim,N2V,array,NCOUNT,1													! Create NCOUNT x 1 array

*vget,N2V(1),node,,nlist												! Fill N2V with node ID, VERTs ID is the vector line
*do,I,1,NCOUNT
	*get,NARRAY(I,1),node,N2V(I),loc,x 									! Fill first column with x-coord.
	*get,NARRAY(I,2),node,N2V(I),loc,y 									! Fill second column with y-coord.
	*get,NARRAY(I,3),node,N2V(I),loc,z 									! Fill third column with z-coord.
*enddo

*cfopen,nodes,dat 														! Create file called “nodes.dat”
*vwrite,NCOUNT,HIGHNODE													! Writing Nodes
%I %I
*vwrite,N2V(1),NARRAY(1,1),NARRAY(1,2),NARRAY(1,3)						! Write NARRAY to file
%I %8.3F %8.3F %8.3F 
*cfclose																! Close file called “nodes.dat”

nsla,S,1																! Select nodes on exterior face (surfaces)
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

*get,ECOUNT,elem,,count 												! Get total number of selected elements
*dim,EARRAY,array,ECOUNT,8												! Create ECOUNT x 8 array
*dim,ELM,array,ECOUNT,1 												! Create ECOUNT x 1 array

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
*vwrite,2*ECOUNT,1,3,ETYPE												! Writing Elements
%8I %8I %8I %8I
*vwrite,ELM(1),EARRAY(1,1),EARRAY(1,2),EARRAY(1,3),EARRAY(1,4),EARRAY(1,5),EARRAY(1,6),EARRAY(1,7),EARRAY(1,8)
%8I %8I %8I %8I %8I %8I %8I %8I %8I

esel,s,type,,1,ETYPE-1,1
ETYPE=ETYPE+1
ensym,,,,all

ET,6,SURF152															! Define Element Type, 6
TYPE,   6   
MAT,       1
REAL,   
ESYS,       0   
SECNUM, 
!*  
ESURF,0 																! Create dummy surface elements
!*  
esel,S,type,,6 															! Select only surface effect elements

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

*vwrite,ELM(1),EARRAY(1,1),EARRAY(1,2),EARRAY(1,3),EARRAY(1,4),EARRAY(1,5),EARRAY(1,6),EARRAY(1,7),EARRAY(1,8)
%8I %8I %8I %8I %8I %8I %8I %8I %8I
*vwrite
END
*cfclose																! Close file called “elements.dat”

esel,S,type,,5,6 

*get,ECOUNT,elem,,count 												! Get total number of selected elements
*dim,EARRAY,array,ECOUNT,8												! Create ECOUNT x 8 array
*dim,ELM,array,ECOUNT,1 												! Create ECOUNT x 1 array

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

*cfopen,2fds_input,geo 													! Create file called “2fds_input.geo”

*vwrite
VERTS=
*vwrite,NARRAY(1,1),NARRAY(1,2),NARRAY(1,3)								! Write NARRAY to file
%8.3F,%8.3F,%8.3F,

*do,I,1,ECOUNT															! Translate the node numbers into VERT numbers
	*do,J,1,NCOUNT														! To write a vector of translation from nodes to VERTs just include 
		*if,EARRAY(I,1),EQ,N2V(J),THEN									! *vwrite,N2V(1)
		*set,EARRAY(I,1),J 												! %I
		*endif															! The value is the node number and the vector line is the VERT number
		*if,EARRAY(I,2),EQ,N2V(J),THEN 									! "N2V(10)=15" means that VERT 10 is node 15
		*set,EARRAY(I,2),J
		*endif
		*if,EARRAY(I,3),EQ,N2V(J),THEN
		*set,EARRAY(I,3),J
		*endif
		*if,EARRAY(I,4),EQ,N2V(J),THEN
		*set,EARRAY(I,4),J
		*endif
	*enddo	
*enddo

*vwrite
FACES=
*do,I,1,ECOUNT/2
*if,EARRAY(I,3),EQ,EARRAY(I,4),then
*vlen,1,1
*vwrite,EARRAY(I,1),EARRAY(I,2),EARRAY(I,3)					 										! Write faces to file (triangles)
%8I,%8I,%8I,
*else
*vlen,1,1
*vwrite,EARRAY(I,1),EARRAY(I,2),EARRAY(I,3),EARRAY(I,1),EARRAY(I,3),EARRAY(I,4)					 	! Write faces to file (quadrilaterals splitted into 2 triangles)
%8I,%8I,%8I, %/%8I,%8I,%8I,
*endif
*enddo	

*do,I,ECOUNT/2+1,ECOUNT
*if,EARRAY(I,3),EQ,EARRAY(I,4),then
*vlen,1,1
*vwrite,EARRAY(I,1),EARRAY(I,2),EARRAY(I,3)					 										! Write faces to file (triangles)
%8I,%8I,%8I,
*else
*vlen,1,1
*vwrite,EARRAY(I,1),EARRAY(I,2),EARRAY(I,4),EARRAY(I,2),EARRAY(I,3),EARRAY(I,4)					 	! Write faces to file (quadrilaterals splitted into 2 triangles)
%8I,%8I,%8I, %/%8I,%8I,%8I,
*endif
*enddo	
RESUME																	! Delete dummy surface effect elements, Element Type, etc...

