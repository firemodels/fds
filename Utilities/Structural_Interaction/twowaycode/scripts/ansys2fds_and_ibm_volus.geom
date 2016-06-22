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
*vwrite,ECOUNT,0,0,ETYPE												! Writing Elements
%8I %8I %8I %8I
*vwrite,ELM(1),EARRAY(1,1),EARRAY(1,2),EARRAY(1,3),EARRAY(1,4),EARRAY(1,5),EARRAY(1,6),EARRAY(1,7),EARRAY(1,8)
%8I %8I %8I %8I %8I %8I %8I %8I %8I
*vwrite
END
*cfclose																! Close file called “elements.dat”

*cfopen,2fds_input,geo 													! Create file called “2fds_input.geo”

*vwrite
&GEOM ID='FEM_MESH',

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
	*enddo	
*enddo

*vwrite
FACES=
*vwrite,EARRAY(1,1),EARRAY(1,2),EARRAY(1,3)					 	! Write faces to file
%8I,%8I,%8I,

RESUME																	! Delete dummy surface effect elements, Element Type, etc...
ALLSEL,ALL 
nsle,S,corner															! Select just the corner nodes (exclude midside nodes)
*get,NCOUNT,node,,count 												! Get total number of selected nodes
*dim,N2V,array,NCOUNT,1													! Create NCOUNT x 1 array
*vget,N2V(1),node,,nlist												! Fill N2V with node ID, VERTs ID is the vector line

*get,VCOUNT,elem,,count 												! Get total number of selected elements
*dim,VARRAY,array,VCOUNT,8												! Create VCOUNT x 8 array
*dim,VELM,array,VCOUNT,1 												! Create VCOUNT x 1 array

*vget,VELM(1),elem,,elist												! Fill ELM column with element ID

*do,I,1,VCOUNT
	*get,VARRAY(I,1),elem,VELM(I),node,1 									! Fill first column with node 1.
	*get,VARRAY(I,2),elem,VELM(I),node,2 									! Fill second column with node 2.
	*get,VARRAY(I,3),elem,VELM(I),node,3 									! Fill third column with node 3.
	*get,VARRAY(I,4),elem,VELM(I),node,4 									! Fill fourth column with node 4.
	*get,VARRAY(I,5),elem,VELM(I),node,5 									! Fill fifth column with node 5.
	*get,VARRAY(I,6),elem,VELM(I),node,6 									! Fill sixth column with node 6.
	*get,VARRAY(I,7),elem,VELM(I),node,7 									! Fill seventh column with node 7.
	*get,VARRAY(I,8),elem,VELM(I),node,8 									! Fill eighth column with node 8.
*enddo

*do,I,1,VCOUNT															! Translate the node numbers into VERT numbers
	*do,J,1,NCOUNT														! To write a vector of translation from nodes to VERTs just include 
		*if,VARRAY(I,1),EQ,N2V(J),THEN									! *vwrite,N2V(1)
		*set,VARRAY(I,1),J 												! %I
		*endif															! The value is the node number and the vector line is the VERT number
		*if,VARRAY(I,2),EQ,N2V(J),THEN 									! "N2V(10)=15" means that VERT 10 is node 15
		*set,VARRAY(I,2),J
		*endif
		*if,VARRAY(I,3),EQ,N2V(J),THEN
		*set,VARRAY(I,3),J
		*endif
		*if,VARRAY(I,4),EQ,N2V(J),THEN
		*set,VARRAY(I,4),J
		*endif
	*enddo	
*enddo

*vwrite
VOLUS=
*vwrite,VARRAY(1,1),VARRAY(1,2),VARRAY(1,3),VARRAY(1,4),			 	! Write VOLUS to file
%8I,%8I,%8I,%8I,	

RESUME																	! Delete dummy surface effect elements, Element Type, etc...

finish