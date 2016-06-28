program barrel
 
WRITE(6,*) ' Enter IJK:'
READ(5,*) IBAR,JBAR,KBAR
WRITE(6,*) ' Enter XB:'
READ(5,*) x1,x2,y1,y2,z1,z2
 
dx  = (x2-x1)/REAL(IBAR)
dy  = (y2-y1)/REAL(JBAR)
dz  = (z2-z1)/REAL(KBAR)
 
zc = -8.6
r0 = 30.9

do k=1,kbar
   z = z1 + (k-.5)*dz
   jf = 0
   do j=1,jbar/2
      y = y1 + (j-.5)*dy
      r  = sqrt(y**2+(z-zc)**2)
      if (r.gt.r0) jf = j
   enddo
   if (jf.gt.0) then
      xs = x1
      xf = x2
      ys = y1
      yf = y1 + dy*jf
      zs = z1 + dz*(k-1)
      zf = z1 + dz*k
      write(8,'(A,6(F5.1,A),A)') "&OBST XB=",xs,",",xf,",",ys,",",yf,",",zs,",",zf,",", " SURF_ID='ROOF' /"
      ys = y2 - dy*jf
      yf = y2
      write(8,'(A,6(F5.1,A),A)') "&OBST XB=",xs,",",xf,",",ys,",",yf,",",zs,",",zf,",", " SURF_ID='ROOF' /"
   endif
enddo
 
STOP
end
