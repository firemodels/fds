      program layer_height
c
      parameter(ntd=20000)
      parameter(ncd=500)
      real d(ntd,ncd),tmp(ntd,ncd),t(ntd),tp(ntd),i1,i2
      real z(0:20),ztc(20)
      integer ind(20),icol(20,20),ntc,ntrees,nn
      real tmpl(ncd),tmph(ncd),zint(ncd),wgt(20)
c
      character(60) infile,outfile
c
      write(6,*) ' Enter number of Trees:'
      read(5,*) ntrees
      write(6,*) ' Enter number of TCs in the tree:'
      read(5,*) ntc
      do n=1,ntc
      write(6,'(a,i2)') ' Enter height (m) and Columns for level ',n
      read(5,*) ztc(n),(icol(n,nn),nn=1,ntrees)
      enddo
      write(6,*) ' Enter weight for each tree:'
      read(5,*) (wgt(nn),nn=1,ntrees)
c
      write(6,*) ' Enter data file name:'
      read(5,'(a)') infile
      write(6,*) ' Enter number of columns in data file:'
      read(5,*) nc
      write(6,*) ' Enter row number where data starts:'
      read(5,*) nr
c
      write(6,*) ' Enter ceiling height (m):'
      read(5,*) z0
c
      write(6,*) ' Enter starting time:'
      read(5,*) t_start
c
      write(6,*) ' Enter name of output file:'
      read(5,'(a)') outfile
c
      open(10,file=infile,status='old',form='formatted')
      do n=1,nr-1
      read(10,*)
      write(6,*) ' Read blank'
      enddo
c
      do 10 i=1,ntd
      read(10,*,end=11) t(i),(d(i,k),k=1,nc-1)
   10 continue
   11 nts = i-1
      close(10)
c
      z(0) = 0.
      do n=1,ntc-1
      z(n) = (ztc(n)+ztc(n+1))/2.
      enddo
      z(ntc) = z0
c
      open(10,file=trim(outfile),form='formatted',status='replace')
c
      write(10,*) '  Time , Height , T_lower , T_upper'
      time_loop: do i=1,nts
      if (t(i).lt.t_start) cycle time_loop
      tmp(i,n) = 0.
      do nn=1,ntrees
      do n=1,ntc                                                                                                
      tmp(i,n) = tmp(i,n) + (273+d(i,icol(n,nn)-1))*wgt(nn)
      enddo                                                                                                     
      enddo                                                                                                     
      i1 = 0.
      i2 = 0.
      do n=1,ntc
      i1 = i1 + tmp(i,n)*(z(n)-z(n-1))
      i2 = i2 + (1./tmp(i,n))*(z(n)-z(n-1))
      enddo
      zint(i)=tmp(i,1)*(i1*i2-z0**2)/(i1+i2*tmp(i,1)**2-2*tmp(i,1)*z0)
      tmpl(i)=tmp(i,1)
      i1 = 0.
      do n=1,ntc
      if (z(n).gt.zint(i)) then
         if (z(n-1).ge.zint(i)) i1 = i1 + tmp(i,n)*(z(n)-z(n-1))
         if (z(n-1).lt.zint(i)) i1 = i1 + tmp(i,n)*(z(n)-zint(i))
         endif
      enddo
      tmph(i) = (1./(z0-zint(i)))*i1
      write(10,'(3(f8.2,a),f8.2)') t(i),',',zint(i),',',
     .     tmpl(i)-273,',',tmph(i)-273
      enddo time_loop
c
      stop
      end
