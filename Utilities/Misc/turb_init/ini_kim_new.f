c
c
c
      program init_rf
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     init_rf "ini_kim_new.f" was modified from "ini_kim.f" to make it
c     work for different size grids without re-compiling. (tms 6/12/00)
c
c     on sgi set: unlimit stack
c     compile on sgi: f77 -r8 -O -o ini_kim_new ini_kim_new.f
c
c     on linux: ulimit -s unlimited
c     compile on intel linux 64: ifort -r8 -mcmodel=medium -O -o turb_init ini_kim_new.f (rjm 11/12/14) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      parameter (nxd = 64, nyd = 64, nzd = 64)
      parameter (nxd = 128, nyd = 128, nzd = 128)

      integer iwk2(0:2*nxd)

      real
     .   wk(2*nxd), wk2(0:2*nxd),
     .   sn(0:nxd), esh(0:nxd),
     .   rindex(nxd), rindex2(nxd)

      real
     .   up(nxd,nyd,nzd), vp(nxd,nyd,nzd), wp(nxd,nyd,nzd),
     .   rad(nxd,nyd,nzd), exx(nxd,nyd,nzd),
     .   ph1(nxd,nyd,nzd), ph2(nxd,nyd,nzd), ph3(nxd,nyd,nzd),
     .   upp(nxd+1,nyd+1,nzd+1),vpp(nxd+1,nyd+1,nzd+1),
     .   wpp(nxd+1,nyd+1,nzd+1),
     .   pp(nxd,nyd,nzd), ekp(nxd,nyd,nzd), src(nxd,nyd,nzd),
     .   cgwk1(nxd,nyd,nzd), cgwk2(nxd,nyd,nzd), cgwk3(nxd,nyd,nzd),
     .   cgwk4(nxd,nyd,nzd), cgwk5(nxd,nyd,nzd), cgwk6(nxd,nyd,nzd),
     .   cgwk7(nxd,nyd,nzd), wtr(nxd,nyd,nzd)

      complex
     .   tmp1(nxd,nyd,nzd), tmp2(nxd,nyd,nzd), tmp3(nxd,nyd,nzd),
     .   tmp4(nxd,nyd,nzd)
	   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(*,*)
c     write(*,*) '   Enter nx, ny, nz: '
      write(*,*) '   Enter nx: '
c     read(*,*) nx, ny, nz
      read(*,*) nx
      ny = nx
      nz = nx
      write(*,*) 'nx,ny,nz:', nx,ny,nz
      write(*,*)

      if ( nx .gt. nxd ) then
       write(*,*) 'FATAL: nx > nxd!'
       stop
      end if

      if ( nx .gt. nxd ) then
       write(*,*) 'FATAL: nx > nxd!'
       stop
      end if

      call ranfld(nx, ny, nz,
     .            iwk2, wk, wk2, sn, esh, rindex, rindex2,
     .            up, vp, wp, rad, exx, ph1, ph2, ph3, 
     .            upp, vpp, wpp, pp, ekp, src,
     .            cgwk1, cgwk2, cgwk3, cgwk4, cgwk5, cgwk6, cgwk7, 
     .            wtr,
     .            tmp1, tmp2, tmp3, tmp4,
     .            mpi,mx,my,mz,ibar,jbar,kbar)

      stop
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
      subroutine ranfld(nx, ny, nz,
     .             iwk2, wk, wk2, sn, esh, rindex, rindex2,
     .             up, vp, wp, rad, exx, ph1, ph2, ph3,
     .             upp, vpp, wpp, pp, ekp, src,
     .             cgwk1, cgwk2, cgwk3, cgwk4, cgwk5, cgwk6, cgwk7,
     .             wtr,
     .             tmp1, tmp2, tmp3, tmp4)

      implicit none

      integer nx, ny, nz, nxh, nyh, nzh

      character*100 fn_iso

      integer
     .   iseed, i, j, k, l, ip1, im1, jp1, jm1, kp1, km1, ifix, 
     .   ip_field

      integer iwk2(0:2*nx)

      real
     .   dx, dy, dz, uref, urms, lnscle, pi, xlen, ylen, zlen, sigma,
     .   umx, umn, vmx, vmn, wmx, wmn,
     .   umean, vmean, wmean, vrms, wrms, ratio, ekmean,
     .   ufac, vfac, wfac, rl11, rl22, rl33, rlkk, 
     .   rt11, rt22, rt33, rtkk,
     .   dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, div,
     .   dudx2, dvdy2, dwdz2, rlamx, rlamy, rlamz, TKE, temp

      real
     .   wk(2*nx),wk2(0:2*nx),
     .   sn(0:nx), esh(0:nx), 
     .   rindex(nx), rindex2(nx)

      real
     .   up(nx,ny,nz), vp(nx,ny,nz), wp(nx,ny,nz), 
     .   rad(nx,ny,nz), exx(nx,ny,nz), 
     .   ph1(nx,ny,nz), ph2(nx,ny,nz), ph3(nx,ny,nz),
     .   upp(nx+1,ny+1,nz+1),vpp(nx+1,ny+1,nz+1),
     .   wpp(nx+1,ny+1,nz+1),
     .   pp(nx,ny,nz), ekp(nx,ny,nz), src(nx,ny,nz),
     .   cgwk1(nx,ny,nz), cgwk2(nx,ny,nz), cgwk3(nx,ny,nz),
     .   cgwk4(nx,ny,nz), cgwk5(nx,ny,nz), cgwk6(nx,ny,nz),
     .   cgwk7(nx,ny,nz), wtr(nx,ny,nz)

      real
     .   rannum

      complex 
     .   tmp1(nx,ny,nz), tmp2(nx,ny,nz), tmp3(nx,ny,nz), tmp4(nx,ny,nz)

      double precision drandm
      double precision dseed

      common/rndo /dseed 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         up(i,j,k) = 0.0
         vp(i,j,k) = 0.0
         wp(i,j,k) = 0.0
         pp(i,j,k) = 0.0
        ekp(i,j,k) = 0.0
        end do
       end do
      end do
c
      do k = 1,nz+1
       do j = 1,ny+1
        do i = 1,nx+1
         upp(i,j,k) = 0.0
         vpp(i,j,k) = 0.0
         wpp(i,j,k) = 0.0
        end do
       end do
      end do
c
      pi = 4.0 * atan(1.0)

c     print*, ' Enter Urms:'
c     read(5,*) urms
c     print*, ' Enter Uref:'
c     read(5,*) uref
c     print*, ' Enter Iseed:'
c     read(5,*) iseed

      iseed = 1000000

      nxh = nx / 2
      nyh = ny / 2
      nzh = nz / 2

      dx = 2.0 * pi / float(nx)
      dy = 2.0 * pi / float(ny)
      dz = 2.0 * pi / float(nz)

      xlen = float(nx) * dx
      ylen = float(ny) * dy
      zlen = float(nz) * dz

      write(*,*) '   Random field Initialization '
      write(*,*) '   =========================== '
      write(*,*) '   generate new field (0) / modify field (1)?'
      read(*,*) ifix
 
      if (ifix .ne. 1 ) then
c
C-------------------------------------------------------------------
C---- warm up the random number generator
C-------------------------------------------------------------------
c      dseed   = 256.d0
c       dseed   = 43756632.d0
       write(*,*)
       write(*,*) ' Enter random seed (dseed): '
       read(*,*) dseed
c       open(unit=10,file='rd_seed.dat')
c       write(10,*) dseed
c       close(10)
       do i = 1,1000
        rannum = drandm(dseed)
       end do
c-------------------------------------------------------------------
c
       call rand3d ( nx, ny, nz, xlen, ylen, zlen, 
     .               uref, urms, ratio, lnscle, iseed, 
     .               up, vp, wp, 
     .               ph1, ph2, ph3, sigma, rad, exx, wk, 
     .               tmp1, tmp2, tmp3, tmp4, 
     .               wtr, sn, esh, rindex, rindex2 )

      elseif ( ifix .eq. 1 ) then

       call fixit( nx, ny, nz, nxh, nyh, nzh, xlen, ylen, zlen, 
     .             up, vp, wp, ekp, 
     .             cgwk1, cgwk2, cgwk3, cgwk4,
     .             tmp1, tmp2, tmp3, tmp4, wk, wk2, iwk2)

      endif
c
c --- Call poisson solver for the pressure field
c
      ip_field = 0
      write(*,*)
      write(*,*) ' Generate pressure field (no=0/yes=1)?'
      read(*,*) ip_field
      if (ip_field .ne. 0) then
       write(*,*) '   Generating Pressure field!'
c      call poisson2(nx,ny,nz,up,vp,wp,pp,cgwk1,src)
c
       call poissonCG(nx,ny,nz,up,vp,wp,pp,src,cgwk1,cgwk2,
     .                cgwk3,cgwk4,cgwk5,cgwk6,cgwk7)
      endif
c
c***********************************
c --- Write out unformatted 
c***********************************
c
c     open(unit=50,file='qphys.dat',form='unformatted')
c     write(50) nx, ny, nz
c     write(50)(((up(i,j,k),i=1,nx),j=1,ny),k=1,nz),
c    .         (((vp(i,j,k),i=1,nx),j=1,ny),k=1,nz),
c    .         (((wp(i,j,k),i=1,nx),j=1,ny),k=1,nz)
c     close(50)
c
c***********************************
c --- Write out for Acusim
c***********************************
c
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         upp(i,j,k) = up(i,j,k)
         vpp(i,j,k) = vp(i,j,k)
         wpp(i,j,k) = wp(i,j,k)
        end do
       end do
      end do
c
      do k = 1,nz
       do j = 1,ny
        upp(nx+1,j,k) = upp(1,j,k)
        vpp(nx+1,j,k) = vpp(1,j,k)
        wpp(nx+1,j,k) = wpp(1,j,k)
       end do
      end do
c
      do k = 1,nz
       do i = 1,nx+1
        upp(i,ny+1,k) = upp(i,1,k)
        vpp(i,ny+1,k) = vpp(i,1,k)
        wpp(i,ny+1,k) = wpp(i,1,k)
       end do
      end do
c
      do j = 1,ny+1
       do i = 1,nx+1
        upp(i,j,nz+1) = upp(i,j,1)
        vpp(i,j,nz+1) = vpp(i,j,1)
        wpp(i,j,nz+1) = wpp(i,j,1)
       end do
      end do
c
c     write(6,*) '   Writing acusim file: kim_ini.dat'
c     open(unit=50,file='kim_ini.dat')
c     write(50,*) nx+1, ny+1, nz+1
c     write(50,*) (((upp(i,j,k),i=1,nx+1),j=1,ny+1),k=1,nz+1),
c    .            (((vpp(i,j,k),i=1,nx+1),j=1,ny+1),k=1,nz+1),
c    .            (((wpp(i,j,k),i=1,nx+1),j=1,ny+1),k=1,nz+1)
c     close(50)
c
c***********************************
c --- Write out for MPSalsa
c***********************************
c
      if ( ifix .ne. 1 ) then
       do k = 1,nz
        do j = 1,ny
         do i = 1,nx
          ekp(i,j,k) = 0.05 * (up(i,j,k)**2+vp(i,j,k)**2+wp(i,j,k)**2)
         end do
        end do
       end do
      endif
      
      write(6,*) '   Writing MPSalsa file: iso_ini.dat'
        
      open(unit=50,file='iso_ini.dat')
      write(50,*) nx, ny, nz
      write(50,52) 0.0, 0.0, 0.0
      write(50,52) nx*dx, ny*dy, nz*dz
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
c           write(50,*) up(i,j,k),vp(i,j,k),wp(i,j,k),pp(i,j,k),ekp(i,j,k)
          write(50,*) up(i,j,k),vp(i,j,k),wp(i,j,k)
        end do
       end do
      end do
      close(50)
52    format(3(f12.7))
c
c --- calculate velocity statistics
c
      urms  = 0.0
      vrms  = 0.0
      wrms  = 0.0
      umean = 0.0
      vmean = 0.0
      wmean = 0.0
      umx   = -10.0e10
      vmx   = -10.0e10
      wmx   = -10.0e10
      umn   =  10.0e10
      vmn   =  10.0e10
      wmn   =  10.0e10

      ekmean = 0.0

      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         umean = umean + up(i,j,k)
         vmean = vmean + vp(i,j,k)
         wmean = wmean + wp(i,j,k)
         urms  = urms + up(i,j,k)**2
         vrms  = vrms + vp(i,j,k)**2
         wrms  = wrms + wp(i,j,k)**2
         umx   = max( umx,up(i,j,k) )
         vmx   = max( vmx,vp(i,j,k) )
         wmx   = max( wmx,wp(i,j,k) )
         umn   = min( umn,up(i,j,k) )
         vmn   = min( vmn,vp(i,j,k) )
         wmn   = min( wmn,wp(i,j,k) )

         ekmean = ekmean + ekp(i,j,k)

        end do
       end do
      end do
      umean = umean / (float(nx)*float(ny)*float(nz))
      vmean = vmean / (float(nx)*float(ny)*float(nz))
      wmean = wmean / (float(nx)*float(ny)*float(nz))
      urms  = sqrt( urms / (float(nx)*float(ny)*float(nz)) )
      vrms  = sqrt( vrms / (float(nx)*float(ny)*float(nz)) )
      wrms  = sqrt( wrms / (float(nx)*float(ny)*float(nz)) )

      TKE   = 0.5 * ( urms**2 + vrms**2 + wrms**2)

      ekmean = ekmean / (float(nx)*float(ny)*float(nz))

      write(6,61) '   umean,vmean,wmean: ',umean,vmean,wmean
      write(6,61) '   urms,vrms,wrms   : ',urms,vrms,wrms
      write(6,61) '   umx,vmx,wmx      : ',umx,vmx,wmx
      write(6,61) '   umn,vmn,wmn      : ',umn,vmn,wmn
61    format(a25,3e12.5)
      write(6,*)  '   TKE              : ', TKE
      write(6,*) '    ekmean           : ', ekmean
c
c --- Estimate the divergence order 2 finite difference
c
      div   = 0.0
      dudx2 = 0.0
      dvdy2 = 0.0
      dwdz2 = 0.0

      do k = 1,nz
       kp1 = k+1
       km1 = k-1
       if (kp1.gt.nz) kp1 = kp1 - nz
       if (km1.lt.1)  km1 = km1 + 1
       do j = 1,nz
        jp1 = j+1
        jm1 = j-1
        if (jp1.gt.ny) jp1 = jp1 - ny
        if (jm1.lt.1)  jm1 = jm1 + 1
        do i = 1,nx
         ip1 = i+1
         im1 = i-1
         if (ip1.gt.nx) ip1 = ip1 - nx
         if (im1.lt.1)  im1 = im1 + 1
c
         dudx = (up(ip1,j,k) - up(im1,j,k))/(2.0*dx)
         dvdx = (vp(ip1,j,k) - vp(im1,j,k))/(2.0*dx)
         dwdx = (wp(ip1,j,k) - wp(im1,j,k))/(2.0*dx)
c
         dudy = (up(i,jp1,k) - up(i,jm1,k))/(2.0*dy)
         dvdy = (vp(i,jp1,k) - vp(i,jm1,k))/(2.0*dy)
         dwdy = (wp(i,jp1,k) - wp(i,jm1,k))/(2.0*dy)
c
         dudz = (up(i,j,kp1) - up(i,j,km1))/(2.0*dz)
         dvdz = (vp(i,j,kp1) - vp(i,j,km1))/(2.0*dz)
         dwdz = (wp(i,j,kp1) - wp(i,j,km1))/(2.0*dz)
c
         div   = div + (dudx + dvdy + dwdz)**2
         dudx2 = dudx2 + dudx*dudx
         dvdy2 = dvdy2 + dvdy*dvdy
         dwdz2 = dwdz2 + dwdz*dwdz

        end do
       end do
      end do
   
      div = sqrt( div / (nx*ny*nz) )
      write(6,*) 'RMS div ', div

      dudx2 = dudx2 / (nx*ny*nz)
      dvdy2 = dvdy2 / (nx*ny*nz)
      dwdz2 = dwdz2 / (nx*ny*nz)

      rlamx = sqrt(urms*urms/dudx2)
      rlamy = sqrt(vrms*vrms/dvdy2)
      rlamz = sqrt(wrms*wrms/dwdz2)

      write(6,*) ' rlamx = ', rlamx
      write(6,*) ' rlamy = ', rlamy
      write(6,*) ' rlamz = ', rlamz
c
c --- Estimate of Integral Length scales
c
      ufac   = 0.0
      vfac   = 0.0
      wfac   = 0.0

      rl11   = 0.0
      rl22   = 0.0
      rl33   = 0.0

      rt11   = 0.0
      rt22   = 0.0
      rt33   = 0.0

      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
c
c --- Longitudinal
c
         do l = 1,nx
          ufac = ( up(i,j,k) - umean ) * ( up(l,j,k) - umean ) * dx
          rl11 = rl11 + ufac
         end do

         do l = 1,ny
          vfac = ( vp(i,j,k) - vmean ) * ( vp(i,l,k) - vmean ) * dy
          rl22 = rl22 + vfac
         end do

         do l = 1,nz
          wfac = ( wp(i,j,k) - wmean ) * ( wp(i,j,l) - wmean ) * dz
          rl33 = rl33 + wfac
         end do
c
c --- Transverse
c
         do l = 1,ny
          ufac = ( up(i,j,k) - umean ) * ( up(i,l,k) - umean ) * dy
          rt11 = rt11 + ufac
         end do

         do l = 1,nz
          vfac = ( vp(i,j,k) - vmean ) * ( vp(i,j,l) - vmean ) * dz
          rt22 = rt22 + vfac
         end do

         do l = 1,nx
          wfac = ( wp(i,j,k) - wmean ) * ( wp(l,j,k) - wmean ) * dx
          rt33 = rt33 + wfac
         end do

        end do
       end do
      end do
c
      rl11 = rl11 / ( urms**2 * nx*ny*nz )
      rl22 = rl22 / ( vrms**2 * nx*ny*nz )
      rl33 = rl33 / ( wrms**2 * nx*ny*nz )
      rlkk = ( rl11 + rl22 + rl33 ) / 3.0
c
      rt11 = rt11 / ( urms**2 * nx*ny*nz )
      rt22 = rt22 / ( vrms**2 * nx*ny*nz )
      rt33 = rt33 / ( wrms**2 * nx*ny*nz )
      rtkk = ( rt11 + rt22 + rt33 ) / 3.0
c
      write(6,*)
      write(6,*) '   Integral scales:'
      write(6,62) '    L11,L22,L33,Lkk: ',rl11,rl22,rl33,rlkk
      write(6,62) '    Lt11,Lt22,Lt33,Ltkk: ',rt11,rt22,rt33,rtkk
      write(6,*)
62    format(a25,4e12.5)
c
c --- Calculate the Energy spectrum of the random velocity field
c
      call fft3d(up, tmp1, nx, ny, nz, wk)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               tmp2(i,j,k) = 0.5 * (tmp1(i,j,k) * conjg(tmp1(i,j,k)))
            end do
         end do
      end do

      call fft3d(vp, tmp1, nx, ny, nz, wk)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               tmp3(i,j,k) = 0.5 * (tmp1(i,j,k) * conjg(tmp1(i,j,k)))
            end do
         end do
      end do

      call fft3d(wp, tmp1, nx, ny, nz, wk)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               tmp4(i,j,k) = 0.5 * (tmp1(i,j,k) * conjg(tmp1(i,j,k)))
            end do
         end do
      end do

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               tmp1(i,j,k) = tmp2(i,j,k) + tmp3(i,j,k) + tmp4(i,j,k)
            end do
         end do
      end do

      call spectrum(tmp1, wk2, nx, ny, nz, iwk2)

      temp = 0.0
      do k = 0,nx
         temp = temp + wk2(k)
      end do
c
c     temp = 0.0
c     do k = 1,nxh
c        temp = temp + 0.5 * ( wk2(k-1) + wk2(k))
c     end do
c
c***********************************
c --- Write out Energy Spectrum
c***********************************
c
      write(*,*)' WARNING:  spectrum not normalized by total energy!'
      write(*,*)' Energy in spectrum, (temp) ', temp
      write(*,*)' SPECTRUM FILE:  ini_kim_spec.dat'

      open(unit=11,file='ini_kim_spec.dat')
      do k = 0,nx
c        write(11,*) k, wk2(k) / temp
         write(11,*) k, wk2(k) 
      end do
      close (11)

      return
      end
c
      subroutine fixit( nx, ny, nz, nxh, nyh, nzh, xlen, ylen, zlen, 
     .                  up, vp, wp, ekp,
     .                  uhat, vhat, what, uvwhat,
     .                  tmp1, tmp2, tmp3, tmp4, 
     .                  wk, wk2, iwk2 )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   
c     Routine fixit reads up,vp,wp and alters the energy spectrum to
c     match a target spectrum, in this case the Comte-Bellot and
c     Corrsin spectrum.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none

      integer kmax, kx, ky, kz, i, j, k, nx, ny, nz, nxh, nyh, nzh, 
     .        numk, ndum, nxdum, nydum, nzdum, knon0

      integer iwk2(0:2*nx)

      real rk, rkx, rky, rkz, ratio, x0, x1, y0, y1, z0, z1,
     .     xlen, ylen, zlen, up, vp, wp, ek, wk2, eratio, dum,
     .     ekp, uhat, vhat, what, uvwhat, eksum, ek_comte, dk,
     .     ek_c_res, ekfac, TKE_mod, TKE_spec

      real wk(2*nx)

      dimension up(nx,ny,nz),vp(nx,ny,nz),wp(nx,ny,nz)
      dimension uhat(nx,ny,nz),vhat(nx,ny,nz),what(nx,ny,nz)
      dimension uvwhat(nx,ny,nz)
      dimension ekp(nx,ny,nz)

      dimension wk2(0:2*nx)
      dimension ek(0:512)

      complex tmp1(nx,ny,nz), tmp2(nx,ny,nz), tmp3(nx,ny,nz), 
     .        tmp4(nx,ny,nz)

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c --- Read up,vp,wp from file
c
      open(unit=50,file='ini_salsa.dat')

      write(*,*) '   Reading: ini_salsa.dat'
      read(50,*) nxdum, nydum, nzdum
      read(50,*) x0, y0, z0
      read(50,*) x1, y1, z1

      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         read(50,*) up(i,j,k), vp(i,j,k), wp(i,j,k)
        end do
       end do
      end do

      close(50)
c
c --- Read target energy spectrum
c
      open(unit=10, file='comte.interp2.dat')
      write(*,*) '   Reading: comte.interp2.dat'
       read(10,*) numk
       do k = 0,numk
        read(10,*) dum, ek(k)

c-rjm   The factor of 100./9. is due to domain scaling
c       needed to match Steve de Bruyn Kops DNS.
c       He uses a domain of 9*2*pi cm.  Since the
c       fft assumes the domain is 2*pi, the energy
c       in the field produced by this program will
c       be multiplied by 9/2pi after the fft to make
c       it dimensionally correct.  The factor of 100.
c       converts from 1/cm to 1/m.

        ek(k) = ek(k)*100./9.

       end do
      close(10)

c
c --- New Energy
c
      TKE_mod = 0.0

      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         TKE_mod = TKE_mod
     .           + 0.5 * (up(i,j,k) * up(i,j,k))
     .           + 0.5 * (vp(i,j,k) * vp(i,j,k))
     .           + 0.5 * (wp(i,j,k) * wp(i,j,k))
        end do
       end do
      end do
      TKE_mod = TKE_mod / (nx*ny*nz)

      write(*,*) 'physical (before fft3d) TKE: ', TKE_mod
c
c --- Transform up,vp,wp
c
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               tmp4(i,j,k) = cmplx(0.0,0.0)
            end do
         end do
      end do

      call fft3d(up, tmp1, nx, ny, nz, wk)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               tmp4(i,j,k) = 0.5 * (tmp1(i,j,k) * conjg(tmp1(i,j,k)))
            end do
         end do
      end do

      call fft3d(vp, tmp2, nx, ny, nz, wk)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               tmp4(i,j,k) = tmp4(i,j,k) 
     .                     + 0.5 * (tmp2(i,j,k) * conjg(tmp2(i,j,k)))
            end do
         end do
      end do

      call fft3d(wp, tmp3, nx, ny, nz, wk)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               tmp4(i,j,k) = tmp4(i,j,k) 
     .                     + 0.5 * (tmp3(i,j,k) * conjg(tmp3(i,j,k)))
            end do
         end do
      end do
c
c --- generate energy spectrum from up,vp,wp
c
      call spectrum(tmp4, wk2, nx, ny, nz, iwk2)
c 
      TKE_spec = 0.0
      knon0 = 0
      do k = 0,nx
       write(6,*) k, wk2(k), ek(k)
       TKE_spec = TKE_spec + wk2(k)
c-rjm  changed to .gt. 1e-10	
       if ( (k.gt.2) .and. ( wk2(k).gt.1e-10 ) ) knon0 = k
      end do
      write(*,*)
      write(*,*) ' k non-zero: ', knon0
      write(*,*) 'Energy in unmodified spectrum: ', TKE_spec
      write(*,*)
c
      ek_c_res = 0.0
      do k = 0,knon0
       ek_c_res = ek_c_res + ek(k)
      end do
      write(*,*)
      write(*,*) '   Comte_Bellot TKE (kmax): ', ek_c_res
      write(*,*)

c
c --- Old Energy
c
      TKE_mod = 0.0

      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         TKE_mod = TKE_mod
     .           + 0.5 * (tmp1(i,j,k) * conjg(tmp1(i,j,k)))
     .           + 0.5 * (tmp2(i,j,k) * conjg(tmp2(i,j,k)))
     .           + 0.5 * (tmp3(i,j,k) * conjg(tmp3(i,j,k)))
        end do
       end do
      end do

      write(*,*)
      write(*,*) 'spectral TKE (just before modifed): ', TKE_mod
      write(*,*)
c
c --- Amplify energy in shells
c
      kmax = nx / 2

      do kx = 1,nx
         rkx = float(kx-1)
         if (rkx .gt. kmax) then
            rkx = kmax+1 - (rkx - kmax+1)
         endif
         do ky = 1,ny
            rky = float(ky-1)
            if (rky .gt. kmax) then
               rky = kmax+1 - (rky - kmax+1)
            endif
            do kz = 1,nz
               rkz = float(kz-1)
               if (rkz .gt. kmax) then
                  rkz = kmax+1 - (rkz - kmax+1)
               endif
               rk = sqrt(rkx*rkx + rkz*rkz + rky*rky)
               k = nint(rk)

               IF ( k .le. knon0 ) THEN
                ratio = sqrt ( ek(k) / wk2(k) )
               ELSE
                ratio = 0.0
               END IF

               tmp1(kx,ky,kz) = ratio * tmp1(kx,ky,kz)
               tmp2(kx,ky,kz) = ratio * tmp2(kx,ky,kz)
               tmp3(kx,ky,kz) = ratio * tmp3(kx,ky,kz)

            end do
         end do
      end do
c
c --- New Energy
c
      TKE_mod = 0.0

      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         TKE_mod = TKE_mod
     .           + 0.5 * real(tmp1(i,j,k) * conjg(tmp1(i,j,k)))
     .           + 0.5 * real(tmp2(i,j,k) * conjg(tmp2(i,j,k)))
     .           + 0.5 * real(tmp3(i,j,k) * conjg(tmp3(i,j,k)))
         tmp4(i,j,k) =  
     .           + 0.5 * (tmp1(i,j,k) * conjg(tmp1(i,j,k)))
     .           + 0.5 * (tmp2(i,j,k) * conjg(tmp2(i,j,k)))
     .           + 0.5 * (tmp3(i,j,k) * conjg(tmp3(i,j,k)))
        end do
       end do
      end do

      write(*,*) 'Modified TKE: ', TKE_mod
c
      call spectrum(tmp4, wk2, nx, ny, nz, iwk2)
c
      TKE_spec = 0.0
      do k = 0,nx
       TKE_spec = TKE_spec + wk2(k)
      end do
      write(*,*)
      write(*,*) 'Energy in modified spectrum: ', TKE_spec
      write(*,*)
c
c --- Transform up,vp,wp back to physical space
c
      call ifft3d(tmp1, up, nx, ny, nz, tmp4, wk)
      call ifft3d(tmp2, vp, nx, ny, nz, tmp4, wk)
      call ifft3d(tmp3, wp, nx, ny, nz, tmp4, wk)
c
c --- New Energy
c
      TKE_mod = 0.0

      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         TKE_mod = TKE_mod
     .           + 0.5 * (up(i,j,k) * up(i,j,k))
     .           + 0.5 * (vp(i,j,k) * vp(i,j,k))
     .           + 0.5 * (wp(i,j,k) * wp(i,j,k))
        end do
       end do
      end do
      TKE_mod = TKE_mod /(nx*ny*nz)

      write(*,*) 'physical TKE (after ifft3d): ', TKE_mod
c
c --- Solve for ksgs
c
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         uhat(i,j,k) = up(i,j,k)*up(i,j,k) 
     .               + vp(i,j,k)*vp(i,j,k)
     .               + wp(i,j,k)*wp(i,j,k)
        end do
       end do
      end do

      call box_fltr(nx, ny, nz, uhat, uvwhat, wk)
      call box_fltr(nx, ny, nz, up, uhat, wk)
      call box_fltr(nx, ny, nz, vp, vhat, wk)
      call box_fltr(nx, ny, nz, wp, what, wk)

      eksum = 0.0
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         ekp(i,j,k) = 0.5 * ( uvwhat(i,j,k) 
     .                      - uhat(i,j,k)*uhat(i,j,k)
     .                      - vhat(i,j,k)*vhat(i,j,k)
     .                      - what(i,j,k)*what(i,j,k) )
         eksum = eksum + ekp(i,j,k)
        end do
       end do 
      end do
      eksum = eksum / (nx*ny*nz)

      ek_comte = 0.0 
      do k = nx+1,numk
       dk = 1.0
       ek_comte = ek_comte + 0.5 * dk * ( ek(k) + ek(k-1) )
      end do

      ratio = ek_comte / eksum

      write(6,*) '   Enter Ksgs factor (0-1):'
      read(*,*)      ekfac

      eksum = 0.0
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         ekp(i,j,k) = ekfac * ratio * ekp(i,j,k)
         eksum      = eksum + ekp(i,j,k)
        end do
       end do
      end do
      eksum = eksum / (nx*ny*nz)

      write(*,*) '   Subgrid ratio: ',ratio
      write(*,*) '   Subgrid ek from comte-bellot spectrum: ',ek_comte
      write(*,*) '   Modified Sugrid energy: ', eksum

      return
      end
c
      subroutine box_fltr(nx, ny, nz, u, uhat, wk)

      dimension u(nx,ny,nz),uhat(nx,ny,nz)
      dimension wk(2*nx)

      rsix = 1.0/6.0

      do i = 1,nx
       do j = 1,ny

        do k = 1,nz
         kp1 = k+1
         if (kp1.gt.nz) kp1 = kp1 - nz
         km1 = k-1
         if (km1.lt.1)  km1 = km1 + nz

         wk(k) = rsix * (u(i,j,km1) + 4.0*u(i,j,k) + u(i,j,kp1))

        end do
        do k = 1,nz
         uhat(i,j,k) = wk(k)
        end do

       end do
      end do
c
      do k = 1,nz
       do i = 1,nx

        do j = 1,ny
         jp1 = j+1
         if (jp1.gt.ny) jp1 = jp1 - ny
         jm1 = j-1
         if (jm1.lt.1)  jm1 = jm1 + ny

         wk(j) = rsix*(uhat(i,jm1,k)+4.0*uhat(i,j,k)+uhat(i,jp1,k))

        end do
        do j = 1,ny
         uhat(i,j,k) = wk(j)
        end do

       end do
      end do
c
      do k = 1,ny
       do j = 1,ny

        do i = 1,nx
         ip1 = i+1
         if (ip1.gt.nx) ip1 = ip1 - nx
         im1 = i-1
         if (im1.lt.1)  im1 = im1 + nx

         wk(i) = rsix*(uhat(im1,j,k)+4.0*uhat(i,j,k)+uhat(ip1,j,k))

        end do
        do i = 1,nx
         uhat(i,j,k) = wk(i)
        end do

       end do
      end do

      return
      end
c
      subroutine poisson2(nx,ny,nz,up,vp,wp,pp,dp,src)

      dimension up(nx,ny,nz),vp(nx,ny,nz),wp(nx,ny,nz),pp(nx,ny,nz)
      dimension src(nx,ny,nz),dp(nx,ny,nz)

      double precision drandm
      double precision dseed

      common/rndo /dseed
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Routine Poisson solves the incompressible, equal spaced Poisson
c     pressure equation to determine pp, the pressure associated with
c     the velocity perturbations.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(6,*)
      write(6,*) '   Poisson Pressure Jacobi Solver'
      write(6,*)

      pi   = 4.0 * atan(1.0)
      dx   = 2.0 * pi / (nx)
      dy   = 2.0 * pi / (ny)
      dz   = 2.0 * pi / (nz)
      r2dx = 1.0 / ( 2.0 * dx )
      r2dy = 1.0 / ( 2.0 * dy )
      r2dz = 1.0 / ( 2.0 * dz )

      sum_div = 0.0
c
c --- Velocity derivative Source term
c
      srcmean = 0.0
      srcrms  = 0.0
      srcmin  = 0.0
      srcmax  = 0.0
c
      do k = 1,nz
       kp1 = k+1
       km1 = k-1
       if (kp1.gt.nz) kp1 = 1
       if (km1.lt.1 ) km1 = nz
       do j = 1,ny
        jp1 = j+1
        jm1 = j-1
        if (jp1.gt.ny) jp1 = jp1 - ny
        if (jm1.lt.1)  jm1 = jm1 + ny
        do i = 1,nx
         ip1 = i+1
         im1 = i-1
         if (ip1.gt.nx) ip1 = ip1 - nx
         if (im1.lt.1)  im1 = im1 + nx

         dudx = ( up(ip1,j,k) - up(im1,j,k) ) * r2dx
         dudy = ( up(i,jp1,k) - up(i,jm1,k) ) * r2dy
         dudz = ( up(i,j,kp1) - up(i,j,km1) ) * r2dz

         dvdx = ( vp(ip1,j,k) - vp(im1,j,k) ) * r2dx
         dvdy = ( vp(i,jp1,k) - vp(i,jm1,k) ) * r2dy
         dvdz = ( vp(i,j,kp1) - vp(i,j,km1) ) * r2dz
 
         dwdx = ( wp(ip1,j,k) - wp(im1,j,k) ) * r2dx
         dwdy = ( wp(i,jp1,k) - wp(i,jm1,k) ) * r2dy
         dwdz = ( wp(i,j,kp1) - wp(i,j,km1) ) * r2dz

         src(i,j,k) = dudx*dudx + dvdy*dvdy + dwdz*dwdz 
     .              + 2.0 * (dudy*dvdx + dwdx*dudz + dvdz*dwdy)

c        r1 = drandm(dseed)
c        src(i,j,k) = 10.0 * (2.0*r1 - 1.0)

c        write(6,*) 'src = ',src(i,j,k)

         sum_div = sum_div + (dudx + dvdy + dwdz)**2

         srcmean = srcmean + src(i,j,k)
         srcmin  = min( srcmin,src(i,j,k) )
         srcmax  = max( srcmax,src(i,j,k) )

        end do
       end do
      end do

      srcmean = srcmean / (nx*ny*nz)

      write(6,*) '   WARNING:  Subtracting mean from src!'
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         srcrms = srcrms + (src(i,j,k) - srcmean)**2
         src(i,j,k) = src(i,j,k) - srcmean
        end do
       end do
      end do

      srcrms = sqrt(srcrms / (nx*ny*nz) )

      write(6,*)
      write(6,*) '   Sum of the div: ', sqrt(sum_div)
      write(6,*) '   Src mean = ',srcmean
      write(6,*) '   Src rms  = ',srcrms
      write(6,*) '   Src min  = ', srcmin
      write(6,*) '   Src max  = ', srcmax
      write(6,*)
c
c --- Iterate Pressure
c
      pref = 0.0
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         pp(i,j,k) = pref
        end do
       end do
      end do
c
      cxyz  = 2.0 * (1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz)) 
      cxyz  = 1.0/cxyz
      cx    = 1.0/(dx*dx)
      cy    = 1.0/(dy*dy)
      cz    = 1.0/(dz*dz)
      coef  = 0.9
      tol   = 1.0e-06
      itmax = 2000
      it    = 0
      dpmax = 1.0e+10

100   continue

      if ( (it.gt.itmax) .or. (dpmax.lt.tol) ) go to 200

      it = it + 1

      do k = 1,nz
       kp1 = k+1
       km1 = k-1
       if (kp1.gt.nz) kp1 = kp1 - nz
       if (km1.lt.1)  km1 = km1 + nz
       do j = 1,ny
        jp1 = j+1
        jm1 = j-1
        if (jp1.gt.ny) jp1 = jp1 - ny
        if (jm1.lt.1)  jm1 = jm1 + ny
        do i = 1,nx
         ip1 = i+1
         im1 = i-1
         if (ip1.gt.nx) ip1 = ip1 - nx
         if (im1.lt.1)  im1 = im1 + nx
 
         dp(i,j,k) = cxyz * ( src(i,j,k)
     .                      + cx * (pp(ip1,j,k) + pp(im1,j,k))
     .                      + cy * (pp(i,jp1,k) + pp(i,jm1,k))
     .                      + cz * (pp(i,j,kp1) + pp(i,j,km1)) )
     .             - pp(i,j,k)
 
        end do
       end do
      end do
c
      dpmax   = 0.0
      rl2norm = 0.0
      pmean   = 0.0
c
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         pp(i,j,k) = pp(i,j,k) + coef * dp(i,j,k)
         pmean     = pmean + pp(i,j,k)
         rl2norm   = rl2norm + dp(i,j,k)*dp(i,j,k)
         dpmax     = max( dpmax,abs(dp(i,j,k)) )
         dp(i,j,k) = 0.0
        end do
       end do
      end do

      pmean = pmean / (nx*ny*nz)

      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         pp(i,j,k) = pp(i,j,k) - pmean
        end do
       end do
      end do

      rl2norm = sqrt ( rl2norm/(nx*ny*nz) )

      if (mod(it,10).eq.0) 
     . write(6,10) 'IT =',it,'  DPMAX = ',dpmax,' L2 = ',rl2norm,
     .             'pmean = ',pmean
10    format(a10,i6,a10,e14.7,a10,e14.7,a10,e14.7)

      go to 100

200   continue

c
c --- Mean Pressure
c
      ppmean = 0.0
      pprms  = 0.0
      ppmn   = 1.0e+10
      ppmx   = 0.0
c    
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         ppmean = ppmean + pp(i,j,k)
         ppmn   = min( ppmn,pp(i,j,k) )
         ppmx   = max( ppmx,pp(i,j,k) )
        end do
       end do
      end do
      ppmean = ppmean / (nx*ny*nz)
c
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         pprms = pprms + (pp(i,j,k) - ppmean)**2
        end do
       end do
      end do
      pprms = sqrt( pprms / (nx*ny*nz) )

      write(6,*)
      write(6,10) 'IT =',it,'  DPMAX = ',dpmax,' L2 = ',rl2norm,
     .            'pmean = ',pmean
      write(6,*) '   Ref. Pressure: ',pref
      write(6,*) '   Mean Pressure: ',ppmean
      write(6,*) '   Max. Pressure: ',ppmx
      write(6,*) '   Min. Pressure: ',ppmn
      write(6,*) '   RMS  Pressure: ',pprms
      write(6,*)
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine poissonCG(nx,ny,nz,up,vp,wp,pp,src,cgwk1,cgwk2,
     .                     cgwk3,cgwk4,cgwk5,cgwk6,cgwk7)

      dimension up(nx,ny,nz),vp(nx,ny,nz),wp(nx,ny,nz),pp(nx,ny,nz)
      dimension src(nx,ny,nz)
      dimension cgwk1(nx,ny,nz),cgwk2(nx,ny,nz),cgwk3(nx,ny,nz)
      dimension cgwk4(nx,ny,nz),cgwk5(nx,ny,nz),cgwk6(nx,ny,nz)
      dimension cgwk7(nx,ny,nz)

      double precision drandm
      double precision dseed

      common/rndo /dseed
      
      common/CG1/ aijk,aim1jk,aip1jk,aijm1k,aijp1k,aijkm1,aijkp1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Routine Poisson solves the incompressible, equal spaced Poisson
c     pressure equation to determine pp, the pressure associated with
c     the velocity perturbations.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(6,*)
      write(6,*) '   Poisson Pressure CG Solver'
      write(6,*)

      pi   = 4.0 * atan(1.0)
      dx   = 2.0 * pi / (nx)
      dy   = 2.0 * pi / (ny)
      dz   = 2.0 * pi / (nz)
      r2dx = 1.0 / ( 2.0 * dx )
      r2dy = 1.0 / ( 2.0 * dy )
      r2dz = 1.0 / ( 2.0 * dz )

      sum_div = 0.0
c
c --- Velocity derivative Source term
c
      srcmean = 0.0
      srcrms  = 0.0
      srcmin  = 0.0
      srcmax  = 0.0
      srcmag  = 1.0e10

      do k = 1,nz
       kp1 = k+1
       km1 = k-1
       if (kp1.gt.nz) kp1 = 1
       if (km1.lt.1 ) km1 = nz
       do j = 1,ny
        jp1 = j+1
        jm1 = j-1
        if (jp1.gt.ny) jp1 = jp1 - ny
        if (jm1.lt.1)  jm1 = jm1 + ny
        do i = 1,nx
         ip1 = i+1
         im1 = i-1
         if (ip1.gt.nx) ip1 = ip1 - nx
         if (im1.lt.1)  im1 = im1 + nx

         dudx = ( up(ip1,j,k) - up(im1,j,k) ) * r2dx
         dudy = ( up(i,jp1,k) - up(i,jm1,k) ) * r2dy
         dudz = ( up(i,j,kp1) - up(i,j,km1) ) * r2dz

         dvdx = ( vp(ip1,j,k) - vp(im1,j,k) ) * r2dx
         dvdy = ( vp(i,jp1,k) - vp(i,jm1,k) ) * r2dy
         dvdz = ( vp(i,j,kp1) - vp(i,j,km1) ) * r2dz

         dwdx = ( wp(ip1,j,k) - wp(im1,j,k) ) * r2dx
         dwdy = ( wp(i,jp1,k) - wp(i,jm1,k) ) * r2dy
         dwdz = ( wp(i,j,kp1) - wp(i,j,km1) ) * r2dz

         src(i,j,k) = dudx*dudx + dvdy*dvdy + dwdz*dwdz
     .              + 2.0 * (dudy*dvdx + dwdx*dudz + dvdz*dwdy)

c        r1 = drandm(dseed)
c        src(i,j,k) = 10.0 * (2.0*r1 - 1.0)

c        write(6,*) 'src = ',src(i,j,k)

         sum_div = sum_div + (dudx + dvdy + dwdz)**2

         srcmean = srcmean + src(i,j,k)
         srcmin  = min( srcmin,src(i,j,k) )
         srcmax  = max( srcmax,src(i,j,k) )

         if ( abs(src(i,j,k)) .lt. srcmag ) then
          srcmag = src(i,j,k)
          imag   = i
          jmag   = j
          kmag   = k
         endif

        end do
       end do
      end do

      srcmean = srcmean / (nx*ny*nz)

      write(6,*) '   WARNING:  Subtracting mean from src!'
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         srcrms = srcrms + (src(i,j,k) - srcmean)**2
         src(i,j,k) = src(i,j,k) - srcmean
        end do
       end do
      end do

      srcrms = sqrt(srcrms / (nx*ny*nz) )

      write(6,*)
      write(6,*) '   RMS of the div: ', sqrt(sum_div/(nx*ny*nz))
      write(6,*) '   Src mean = ',srcmean
      write(6,*) '   Src rms  = ',srcrms
      write(6,*) '   Src min  = ', srcmin
      write(6,*) '   Src max  = ', srcmax
      write(6,*) '   Src mag  = ', srcmag,imag,jmag,kmag
      write(6,*)
c 
c --- Calculate coefficient of A_ij
c
      aijk   =  2.0 * (1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz))
      aim1jk = -1.0/(dx*dx)
      aip1jk = -1.0/(dx*dx)
      aijm1k = -1.0/(dy*dy)
      aijp1k = -1.0/(dy*dy)
      aijkm1 = -1.0/(dz*dz)
      aijkp1 = -1.0/(dz*dz)
c
      tol   = 1.0e-04
      maxit = 2000
c
      call cgsovle(nx,ny,nz,tol,maxit,pp,src,
     .             cgwk1,cgwk2,cgwk3,cgwk4,cgwk5,cgwk6,cgwk7)
c
c --- Mean Pressure
c
      pref   = 0.0
      ppmean = 0.0
      pprms  = 0.0
      ppmn   = 1.0e+10
      ppmx   = 0.0
c
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         ppmean = ppmean + pp(i,j,k)
         ppmn   = min( ppmn,pp(i,j,k) )
         ppmx   = max( ppmx,pp(i,j,k) )
        end do
       end do
      end do
      ppmean = ppmean / (nx*ny*nz)
c
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         pprms = pprms + (pp(i,j,k) - ppmean)**2
        end do
       end do
      end do
      pprms = sqrt( pprms / (nx*ny*nz) )

      write(6,*)
      write(6,*) '   Ref. Pressure: ',pref
      write(6,*) '   Mean Pressure: ',ppmean
      write(6,*) '   Max. Pressure: ',ppmx
      write(6,*) '   Min. Pressure: ',ppmn
      write(6,*) '   RMS  Pressure: ',pprms
      write(6,*)
c
c --- Write out pp function file
c
ccc   open(unit=52,file='pp.data',form='unformatted')
c     write(52) nx,ny,nz
c     write(52) sngl(0.0d0),sngl(0.0d0),sngl(0.0d0),sngl(0.0d0)
c     write(52) ((sngl(up(i,j)), i=1,nx),j=1,ny),
c    .          ((sngl(vp(i,j)), i=1,nx),j=1,ny),
c    .          ((sngl(pp(i,j)), i=1,nx),j=1,ny),
c    .          ((sngl(src(i,j)),i=1,nx),j=1,ny)
ccc   close(52)
c
      return
      end
c
      subroutine cgsovle(nx,ny,nz,tol,maxit,xk,b,
     .                   ax,pk,pkm1,rk,rkm1,zk,zkm1)

      dimension ax(nx,ny,nz), xk(nx,ny,nz), b(nx,ny,nz)

      dimension pk(nx,ny,nz),pkm1(nx,ny,nz)
      dimension rk(nx,ny,nz),rkm1(nx,ny,nz)
      dimension zk(nx,ny,nz),zkm1(nx,ny,nz)

      common/CG1/ aijk,aim1jk,aip1jk,aijm1k,aijp1k,aijkm1,aijkp1

      it = 0
      rl2norm = 1.0e10
c
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         xk(i,j,k) = 0.0
        end do
       end do
      end do
c
      call axmult(nx,ny,nz,xk,ax)
c
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         rk(i,j,k) = b(i,j,k) - ax(i,j,k)
        end do
       end do
      end do

100   continue

      if ((rl2norm .lt. tol) .or. (it.gt.maxit) ) goto 200

      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         zk(i,j,k) = rk(i,j,k) / aijk
        end do
       end do
      end do
c
      it = it + 1
c      
      if ( it .eq. 1 ) then

       do k = 1,nz
        do j = 1,ny
         do i = 1,nx
          pk(i,j,k) = zk(i,j,k)
         end do
        end do
       end do

      else

       beta1 = 0.0
       beta2 = 0.0
       do k = 1,nz
        do j = 1,ny
         do i = 1,nx
          beta1 = beta1 + rk(i,j,k) * zk(i,j,k)
          beta2 = beta2 + rkm1(i,j,k) * zkm1(i,j,k)
         end do
        end do
       end do
       beta = beta1/beta2
c
       do k = 1,nz
        do j = 1,ny
         do i = 1,nx
          pk(i,j,k) = zk(i,j,k) + beta * pk(i,j,k)
         end do
        end do
       end do

      endif
c
      call axmult(nx,ny,nz,pk,ax)
c
      alpha1 = 0.0
      alpha2 = 0.0
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         rkm1(i,j,k) = rk(i,j,k)
         zkm1(i,j,k) = zk(i,j,k)
         alpha1 = alpha1 + rk(i,j,k) * zk(i,j,k)
         alpha2 = alpha2 + pk(i,j,k) * ax(i,j,k)
        end do
       end do
      end do
      alpha = alpha1/alpha2
c
      rl2norm = 0.0
      dxkmax  = 0.0
c
      xkmean = 0.0
      do k = 1,nz
       do j = 1,ny
        do i = 1,nx
         xk(i,j,k) = xk(i,j,k) + alpha * pk(i,j,k)
         rk(i,j,k) = rk(i,j,k) - alpha * ax(i,j,k)

         rl2norm   = rl2norm + ( alpha*pk(i,j,k) )**2
         dxkmax    = max( dxkmax, alpha*pk(i,j,k) )
         xkmean    = xkmean + xk(i,j,k)
        end do
       end do
      end do
c
      rl2norm = sqrt( rl2norm / (nx*ny*nz) )
      xkmean  = xkmean / (nx*ny*nz)
c
      if (mod(it,10).eq.0)
     . write(6,10) 'IT =',it,' DXKMAX = ',dxkmax,' L2 = ',rl2norm,
     . ' xkmean = ',xkmean
10    format(a10,i6,a10,e12.5,a10,e12.5,a10,e12.5)

      goto 100

200   continue

      return
      end
c
      subroutine axmult(nx,ny,nz,x,ax)

      common/CG1/ aijk,aim1jk,aip1jk,aijm1k,aijp1k,aijkm1,aijkp1
      dimension x(nx,ny,nz), ax(nx,ny,nz)

      do k = 1,nz
       kp1 = k+1
       km1 = k-1
       if (kp1.gt.nz) kp1 = kp1 - nz
       if (km1.lt.1)  km1 = km1 + nz
       do j = 1,ny
        jp1 = j+1
        jm1 = j-1
        if (jp1.gt.ny) jp1 = jp1 - ny
        if (jm1.lt.1)  jm1 = jm1 + ny
        do i = 1,nx
         ip1 = i+1
         im1 = i-1
         if (ip1.gt.nx) ip1 = ip1 - nx
         if (im1.lt.1)  im1 = im1 + nx

        ax(i,j,k) = aijk   * x(i,j,k)
     .            + aim1jk * x(im1,j,k)
     .            + aip1jk * x(ip1,j,k)
     .            + aijm1k * x(i,jm1,k)
     .            + aijp1k * x(i,jp1,k)
     .            + aijkm1 * x(i,j,km1)
     .            + aijkp1 * x(i,j,kp1)

        end do
       end do
      end do
       
      return
      end
C
C********************************************************************
      double precision function drandm(dl)
C********************************************************************
C
      double precision dl
      dl = dmod(16807.0d0*dl,2147483647.0d0)
      drandm = dl * 4.6566128752458d-10
C--------------------------------------------------------------------
      return
      end
c
      subroutine rand3d
     .  (nx, ny, nz, xlen, ylen, zlen, 
     .   uref, urms, ratio, lnscle, iseed, 
     .   up, vp, wp, ph1, ph2, ph3, sigma, rad, exx, wk, 
     .   tmp1, tmp2, tmp3, tmp4, 
     .   wtr, sn, esh, rindex, rindex2 )
 
      implicit none

      integer
     .   nx, ny, nz, nxh, nyh, nzh, nxhp, nyhp, nzhp

      integer
     .   i, j, k, m,
     .   irkmx, iseed, ii, jj, kk, im, ip, jm, jp, km, kp,
     .   ishift, jshift, kshift, kseed,
     .   im1, im2, ip1, ip2, jm1, jm2, jp1, jp2,
     .   km1, km2, kp1, kp2, ik, iks

      real
     .   xlen, ylen, zlen, uref, urms, lnscle, pi, dkx, dky, dkz, 
     .   dkx2, dky2, dkz2, a, b, c, d,
     .   rkmax, fac, akx, aky, akz, 
     .   etot, enorm, sigma, eres,
     .   umcon, sp1, rand1, phmax, ur, ui, ratio, 
     .   ran, rk, temp, rk0, x, y, z, vr, vi, dx, dy, dz, sp2, 
     .   rand2, sp3, rand3, wr, wi, r1, rfac, r12dx, r12dy, r12dz

      real wk(2*nx)

      real
     .   up(nx,ny,nz), vp(nx,ny,nz), wp(nx,ny,nz), 
     .   rad(nx,ny,nz), exx(nx,ny,nz), 
     .   ph1(nx,ny,nz), ph2(nx,ny,nz), ph3(nx,ny,nz)
  
      real cpeak

      real sn(0:nx),esh(0:nx)

      real wtr(nx,ny,nz),rindex(nx),rindex2(nx),rkcut

      complex  ai

      complex 
     .   tmp1(nx,ny,nz), tmp2(nx,ny,nz), tmp3(nx,ny,nz), tmp4(nx,ny,nz)

      double precision drandm
      double precision dseed

      common/rndo /dseed
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      nxh  = nx / 2
      nyh  = ny / 2
      nzh  = nz / 2
      nxhp = nx/2+1
      nyhp = ny/2+1
      nzhp = nz/2+1

c   Set up scaling factor (so that up() will have correct urms value)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               exx(i,j,k) = 0.0
            end do
         end do
      end do

c     write random data input info into tape 44

      write(6,60)
c     write(6,61) nx, ny, nz, xlen, ylen, zlen
c     write(6,64) nxh, nyh, nzh
c     write(6,62) uref, urms, lnscle

 60   format(10x, 'The Turbulent Initialization Conditions',/)
 61   format(1x, 'nkx=', i4, 3x, 'nky=', i4, 3x, 'nkz=', i4, 3x, 
     .     'xlen=', f8.5, 3x, 'ylen=', f8.5, 3x, 'zlen=', f8.5)
 62   format(10x,'uref,urms,lnscle',5x,3e13.5)
 64   format(1x, 'nkxh=', i4, 3x, 'nkyh=', i4, 3x, 'nkzh=', i4)

      do i=1,nx
         rindex(i) =float(i)-1.
         rindex2(i)=float(i-1)*float(i-1)
      end do

      pi = 4.0 * atan(1.0)

      dkx = 2.0 * pi / xlen 
      dky = 2.0 * pi / ylen
      dkz = 2.0 * pi / zlen
      dkx2=dkx*dkx
      dky2=dky*dky
      dkz2=dkz*dkz
      rkcut=float(nxh)*sqrt(8./9.)

      do k=1,nzh+1
      do j=1,nyh+1
      do i=1,nxh+1
         wtr(i,j,k)=sqrt(rindex2(i)*dkx2+rindex2(j)*dky2
     1                  +rindex2(k)*dkz2)
      end do
      end do
      end do
      do k = 1,nzh+1
      do j = 1,nyh+1
      do i = 2,nxh
         ii = nx - i + 2
         wtr(ii,j,k)=sqrt(rindex2(i)*dkx2+rindex2(j)*dky2
     1                  +rindex2(k)*dkz2)

      end do
      end do
      end do

      do k = 1,nzh+1
      do j = 2,nyh
      do i = 2,nxh
         ii = nx - i + 2
         jj = ny - j + 2
         wtr(ii,jj,k)=sqrt(rindex2(i)*dkx2+rindex2(j)*dky2
     1                  +rindex2(k)*dkz2)

      end do
      end do
      end do

      do k = 2,nzh
      do j = 1,nyh+1
      do i = 2,nxh
         ii = nx - i + 2
         kk = nz - k + 2
         wtr(ii,j,kk)=sqrt(rindex2(i)*dkx2+rindex2(j)*dky2
     1                  +rindex2(k)*dkz2)

      end do
      end do
      end do

      do k = 2,nzh
      do j = 2,nyh
      do i = 2,nxh
         ii = nx - i + 2
         jj = ny - j + 2
         kk = nz - k + 2
         wtr(ii,jj,kk)=sqrt(rindex2(i)*dkx2+rindex2(j)*dky2
     1                  +rindex2(k)*dkz2)

      end do
      end do
      end do

      do k = 1,nzh+1
      do j = 2,nyh
      do i = 1,nxh+1
         jj = ny - j + 2
         wtr(i,jj,k)=sqrt(rindex2(i)*dkx2+rindex2(j)*dky2
     1                  +rindex2(k)*dkz2)

      end do
      end do
      end do

      do k = 2,nzh
      do j = 1,nyh+1
      do i = 1,nxh+1
         kk = nz - k + 2
         wtr(i,j,kk)=sqrt(rindex2(i)*dkx2+rindex2(j)*dky2
     1                  +rindex2(k)*dkz2)

      end do
      end do
      end do

      do k = 2,nzh
      do j = 2,nyh
      do i = 1,nxh+1
         jj = ny - j + 2
         kk = nz - k + 2
         wtr(i,jj,kk)=sqrt(rindex2(i)*dkx2+rindex2(j)*dky2
     1                  +rindex2(k)*dkz2)

      end do
      end do
      end do

      do k=1,nz
      do j=1,ny
      do i=1,nx
      if (wtr(i,j,k).le.0.0001)wtr(i,j,k)=0.
      if (wtr(i,j,k).gt.0.0001.and.wtr(i,j,k).le.rkcut)wtr(i,j,k)=1.
      if (wtr(i,j,k).gt.rkcut)wtr(i,j,k)=0.
      end do
      end do
      end do

c-rjm      write(6,*) 'cpeak='
c-rjm      read(5,*) cpeak
c-rjm      rk0 = cpeak * max(dkx,dky,dkz)

      do i = 0,nxh
         sn(i) = 0.
      end do

      do k = 1,nz
        kk = k
        if(k.gt.nzhp) kk = nz - k + 2
        km = kk-1
         do j = 1,ny
           jj = j
           if(j.gt.nyhp) jj = ny - j + 2
           jm = jj-1
            do i = 1,nx
              ii = i
              if(i.gt.nxhp) ii = nx - i + 2
              im = ii-1

                iks= im*im+jm*jm+km*km
                rk = sqrt(float(iks))
                do m = 0,nxh
                  if(rk.gt.float(m)-.5.and.rk.le.float(m)+.5 
     1          .and.rk.le.float(nxh)*sqrt(8./9.))
     1            sn(m) = sn(m) + 1.
                end do

            end do
         end do
       end do

c     tsum = 0
c     do i = 0,nxh
c        write(*,*) i,sn(i)
c        tsum = tsum + sn(i)
c     end do
c     write(*,*) tsum
c
c --- esh - E in shell k could be replaced with Comte-Bellot and Corrsin E(k)
c --- 8/9 is the 2/3 rule in 3D
c
c --- Original energy spectrum
c
ccc   do i = 0,nxh
c        rk = float(i)
c        esh(i) = rk**4./(1.+(rk/rk0)**(5./3.+4.))
ccc   end do
c
c --- Comte-Bellot and Corrsin Energy spectrum
c --- Kim's way of doing things.
c     Kim uses t*=0.1, we use seconds, so a factor of 100 must
c     be multiplied to esh(i), this is done at the bottom of the loop.
c     (t.m.smith 07/10/00)
c
       do i = 0,nxh
          rk = float(i)
               if(rk.lt.2.021) then
               a = 1.617
               c = 2.411e-3
               b = 2.021
               d = 4.352e-3
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)

c      added by rjm, 09/04/02
       IF (rk .eq. 0.0) rk = EPSILON(1.0)

               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.2.426) then
               a = 2.021
               c = 4.352e-3
               b = 2.426
               d = 6.093e-3
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.3.234) then
               a = 2.426
               c = 6.093e-3
               b = 3.234
               d = 8.231e-3
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.4.043) then
               a = 3.234
               c = 8.231e-3
               b = 4.043
               d = 8.647e-3
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.5.660) then
               a = 4.043
               c = 8.647e-3
               b = 5.660
               d = 7.190e-3
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.8.085) then
               a = 5.660
               c = 7.190e-3
               b = 8.085
               d = 5.109e-3
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.12.128) then
               a = 8.085
               c = 5.109e-3
               b = 12.128
               d = 3.179e-3
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.16.170) then
               a = 12.128
               c = 3.179e-3
               b = 16.170
               d = 2.271e-3
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.20.213) then
               a = 16.170
               c = 2.271e-3
               b = 20.213
               d = 1.684e-3
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.24.255) then
               a = 20.213
               c = 1.684e-3
               b = 24.255
               d = 1.330e-3
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.32.340) then
               a = 24.255
               c = 1.330e-3
               b = 32.340
               d = 8.893e-4
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.48.510) then
               a = 32.340
               c = 8.893e-4
               b = 48.510
               d = 4.674e-4
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.64.681) then
               a = 48.510
               c = 4.674e-4
               b = 64.681
               d = 2.384e-4
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.80.851) then
               a = 64.681
               c = 2.384e-4
               b = 80.851
               d = 1.404e-4
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.101.063) then
               a = 80.851
               c = 1.404e-4
               b = 101.063
               d = 7.493e-5
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.121.276) then
               a = 101.063
               c = 7.493e-5
               b = 121.276
               d = 4.409e-5
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else if(rk.lt.141.489) then
               a = 121.276
               c = 4.409e-5
               b = 141.489
               d = 2.535e-5
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               go to 999
               else
               a = 141.489
               c = 2.535e-5
               b = 161.701
               d = 1.514e-5
               a = log10(a)
               b = log10(b)
               c = log10(c)
               d = log10(d)
               rk= log10(rk)
               esh(i) =10.**((d-c)/(b-a)*(rk-a)+c)
               end if
  999          continue

               esh(i) = esh(i) * 100.0

       end do
c
      esh(0) = 0.
      do m = 1,nxh
        if(float(m).gt.float(nxh)*sqrt(8./9.)) esh(m) = 0.
      end do

      eres = 0.

      do i = 0,nxh

               rk = sqrt(float(i))

               if(rk.gt.0..and.rk.le.float(nxh)*sqrt(8./9.))
     1            eres=eres+esh(i)

      end do

      write(*,*) 'eres = ',eres

      do k = 1,nzh+1
         akz = dkz * float(k-1)
         km = k-1
         do j = 1,nyh+1
            aky = dky * float(j-1)
            jm = j-1
            do i = 1,nxh+1
               akx = dkx * float(i-1)
               im = i-1
            
               rk = sqrt(aky**2 + akx**2 + akz**2)
               rad(i,j,k) = rk

               do m = 0,nxh
               if( rk.gt.float(m)-.5.and.rk.le.float(m)+.5
     .            .and.rk.le.float(nxh)*sqrt(8./9.) ) then
               ik = m
               exx(i,j,k) = esh(ik)/sn(ik)
               tmp4(i,j,k)= esh(ik)/sn(ik)
               go to 99
               end if
               end do
               exx(i,j,k) = 0.
               tmp4(i,j,k)= 0.
   99          continue

            end do
         end do
      end do

      exx(1,1,1)  = 0.
      tmp4(1,1,1) = 0.

      do k = 1,nzh+1
         do j = 1,nyh+1
            do i = 2,nxh

               ii = nx - i + 2

               exx(ii,j,k)  = exx(i,j,k)
               tmp4(ii,j,k) = tmp4(i,j,k)

            end do
         end do
      end do

      do k = 1,nzh+1
         do j = 2,nyh
            do i = 2,nxh

               ii = nx - i + 2
               jj = ny - j + 2

               exx(ii,jj,k)  = exx(i,j,k)
               tmp4(ii,jj,k) = tmp4(i,j,k)

            end do
         end do
      end do

      do k = 2,nzh
         do j = 1,nyh+1
            do i = 2,nxh

               ii = nx - i + 2
               kk = nz - k + 2

               exx(ii,j,kk)  = exx(i,j,k)
               tmp4(ii,j,kk) = tmp4(i,j,k)

            end do
         end do
      end do

      do k = 2,nzh
         do j = 2,nyh
            do i = 2,nxh

               ii = nx - i + 2
               jj = ny - j + 2
               kk = nz - k + 2

               exx(ii,jj,kk)  = exx(i,j,k)
               tmp4(ii,jj,kk) = tmp4(i,j,k)

            end do
         end do
      end do

      do k = 1,nzh+1
         do j = 2,nyh
            do i = 1,nxh+1

               jj = ny - j + 2

               exx(i,jj,k)  = exx(i,j,k)
               tmp4(i,jj,k) = tmp4(i,j,k)

            end do
         end do
      end do

      do k = 2,nzh
         do j = 1,nyh+1
            do i = 1,nxh+1

               kk = nz - k + 2

               exx(i,j,kk)  = exx(i,j,k)
               tmp4(i,j,kk) = tmp4(i,j,k)

            end do
         end do
      end do

      do k = 2,nzh
         do j = 2,nyh
            do i = 1,nxh+1

               jj = ny - j + 2
               kk = nz - k + 2

               exx(i,jj,kk)  = exx(i,j,k)
               tmp4(i,jj,kk) = tmp4(i,j,k)

            end do
         end do
      end do

      do k=1,nz
      do j=1,ny
      do i=1,nx
         exx(i,j,k)  = wtr(i,j,k)*exx(i,j,k)
         tmp4(i,j,k) = wtr(i,j,k)*tmp4(i,j,k)
      end do
      end do
      end do

c      write(6,*) '   EXX TEST '
c      write(6,121) ((exx(i,j,2),i=1,8),j=1,8)
c 121  format (1x,8(1pe12.5))

c   Normalize target array

      etot = 0.0
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               etot = etot + exx(i,j,k)
            end do
         end do
      end do

      write(*,*) 'etot = ',etot

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
c              exx(i,j,k) = exx(i,j,k) / etot *.5/3.
               exx(i,j,k) = exx(i,j,k) / 3.
            end do
         end do
      end do
 
c   Print out target spectrum

c      write(6,1769)
c1769  format(20x,'exx(i,j,2)')

c      do j = 1,ny
c         write(6,1760) (exx(i,j,2),i=1,nxh)
c      end do

c1760  format(1x,9e13.5)

c   Calculate total energy in target spectrum

      enorm = 0.0
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               enorm = enorm + exx(i,j,k)
            end do
         end do
      end do

c   Convert target spectrum into physical velocities

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               tmp1(i,j,k) = cmplx(0.0,0.0)
               tmp2(i,j,k) = cmplx(0.0,0.0)
               tmp3(i,j,k) = cmplx(0.0,0.0)
            end do
         end do
      end do

      etot = 0.0
      phmax = pi

      do k = 1,nzh+1
         do j = 1,nyh+1
            do i = 1,nxh+1

c   Calculate velocity magnitude

               umcon = sqrt(2.0 * exx(i,j,k))

c   Add in random phase information

c              sp1 = ran(iseed)
c              sp2 = ran(iseed)
c              sp3 = ran(iseed)

               sp1 = drandm(dseed)

               sp2 = drandm(dseed)

               sp3 = drandm(dseed)

c --- range -1 to +1
               rand1 = 2.0 * sp1 - 1.0
               rand2 = 2.0 * sp2 - 1.0
               rand3 = 2.0 * sp3 - 1.0

               ph1(i,j,k) = rand1 * phmax
               ph2(i,j,k) = rand2 * phmax
               ph3(i,j,k) = rand3 * phmax

               if((i.eq.1.and.j.ne.1.and.k.ne.1).or.
     1            (i.ne.1.and.j.eq.1.and.k.ne.1).or.
     1            (i.ne.1.and.j.ne.1.and.k.eq.1).or.
     1            (i.eq.1.and.j.eq.1.and.k.eq.1)) then
                  if(ph1(i,j,k).ge.0.) then
                     ph1(i,j,k) = 0.
                  else
                     ph1(i,j,k) = pi
                  end if
                  if(ph2(i,j,k).ge.0.) then
                     ph2(i,j,k) = 0.
                  else
                     ph2(i,j,k) = pi
                  end if
                  if(ph3(i,j,k).ge.0.) then
                     ph3(i,j,k) = 0.
                  else
                     ph3(i,j,k) = pi
                  end if
               end if

               ur = umcon * cos(ph1(i,j,k))
               ui =-umcon * sin(ph1(i,j,k))
               tmp1(i,j,k) = cmplx(ur,ui)

               vr = umcon * cos(ph2(i,j,k))
               vi =-umcon * sin(ph2(i,j,k))
               tmp2(i,j,k) = cmplx(vr,vi)

               wr = umcon * cos(ph3(i,j,k))
               wi =-umcon * sin(ph3(i,j,k))
               tmp3(i,j,k) = cmplx(wr,wi)

               etot = etot + 0.5 * umcon**2

            end do
         end do
      end do

      tmp1(1,1,1) = 0.
      tmp2(1,1,1) = 0.
      tmp3(1,1,1) = 0.

      do k = 1,nzh+1
         do j = 1,nyh+1
            do i = 2,nxh

               ii = nx - i + 2

               umcon = sqrt(2.0 * exx(i,j,k))
               ur = umcon * cos(ph1(i,j,k))
               ui =-umcon * sin(ph1(i,j,k))
               vr = umcon * cos(ph2(i,j,k))
               vi =-umcon * sin(ph2(i,j,k))
               wr = umcon * cos(ph3(i,j,k))
               wi =-umcon * sin(ph3(i,j,k))

               tmp1(ii,j,k) = cmplx(ur,-ui)
               tmp2(ii,j,k) = cmplx(vr,-vi)
               tmp3(ii,j,k) = cmplx(wr,-wi)
               etot = etot + 0.5 * umcon**2

            end do
         end do
      end do

      do k = 1,nzh+1
         do j = 2,nyh
            do i = 2,nxh

               ii = nx - i + 2
               jj = ny - j + 2

               umcon = sqrt(2.0 * exx(i,j,k))
               ur = umcon * cos(ph1(i,j,k))
               ui =-umcon * sin(ph1(i,j,k))
               vr = umcon * cos(ph2(i,j,k))
               vi =-umcon * sin(ph2(i,j,k))
               wr = umcon * cos(ph3(i,j,k))
               wi =-umcon * sin(ph3(i,j,k))

               tmp1(ii,jj,k) = cmplx(ur,ui)
               tmp2(ii,jj,k) = cmplx(vr,vi)
               tmp3(ii,jj,k) = cmplx(wr,wi)
               etot = etot + 0.5 * umcon**2

            end do
         end do
      end do

      do k = 2,nzh
         do j = 1,nyh+1
            do i = 2,nxh

               ii = nx - i + 2
               kk = nz - k + 2

               umcon = sqrt(2.0 * exx(i,j,k))
               ur = umcon * cos(ph1(i,j,k))
               ui =-umcon * sin(ph1(i,j,k))
               vr = umcon * cos(ph2(i,j,k))
               vi =-umcon * sin(ph2(i,j,k))
               wr = umcon * cos(ph3(i,j,k))
               wi =-umcon * sin(ph3(i,j,k))

               tmp1(ii,j,kk) = cmplx(ur,ui)
               tmp2(ii,j,kk) = cmplx(vr,vi)
               tmp3(ii,j,kk) = cmplx(wr,wi)
               etot = etot + 0.5 * umcon**2

            end do
         end do
      end do

      do k = 2,nzh
         do j = 2,nyh
            do i = 2,nxh

               ii = nx - i + 2
               jj = ny - j + 2
               kk = nz - k + 2

               umcon = sqrt(2.0 * exx(i,j,k))
               ur = umcon * cos(ph1(i,j,k))
               ui =-umcon * sin(ph1(i,j,k))
               vr = umcon * cos(ph2(i,j,k))
               vi =-umcon * sin(ph2(i,j,k))
               wr = umcon * cos(ph3(i,j,k))
               wi =-umcon * sin(ph3(i,j,k))

               tmp1(ii,jj,kk) = cmplx(ur,-ui)
               tmp2(ii,jj,kk) = cmplx(vr,-vi)
               tmp3(ii,jj,kk) = cmplx(wr,-wi)
               etot = etot + 0.5 * umcon**2

            end do
         end do
      end do

      do k = 1,nzh+1
         do j = 2,nyh
            do i = 1,nxh+1

               jj = ny - j + 2

               umcon = sqrt(2.0 * exx(i,j,k))
               ur = umcon * cos(ph1(i,j,k))
               ui =-umcon * sin(ph1(i,j,k))
               vr = umcon * cos(ph2(i,j,k))
               vi =-umcon * sin(ph2(i,j,k))
               wr = umcon * cos(ph3(i,j,k))
               wi =-umcon * sin(ph3(i,j,k))

               tmp1(i,jj,k) = cmplx(ur,-ui)
               tmp2(i,jj,k) = cmplx(vr,-vi)
               tmp3(i,jj,k) = cmplx(wr,-wi)
               etot = etot + 0.5 * umcon**2

            end do
         end do
      end do

      do k = 2,nzh
         do j = 1,nyh+1
            do i = 1,nxh+1

               kk = nz - k + 2

               umcon = sqrt(2.0 * exx(i,j,k))
               ur = umcon * cos(ph1(i,j,k))
               ui =-umcon * sin(ph1(i,j,k))
               vr = umcon * cos(ph2(i,j,k))
               vi =-umcon * sin(ph2(i,j,k))
               wr = umcon * cos(ph3(i,j,k))
               wi =-umcon * sin(ph3(i,j,k))

               tmp1(i,j,kk) = cmplx(ur,-ui)
               tmp2(i,j,kk) = cmplx(vr,-vi)
               tmp3(i,j,kk) = cmplx(wr,-wi)
               etot = etot + 0.5 * umcon**2

            end do
         end do
      end do

      do k = 2,nzh
         do j = 2,nyh
            do i = 1,nxh+1

               jj = ny - j + 2
               kk = nz - k + 2

               umcon = sqrt(2.0 * exx(i,j,k))
               ur = umcon * cos(ph1(i,j,k))
               ui =-umcon * sin(ph1(i,j,k))
               vr = umcon * cos(ph2(i,j,k))
               vi =-umcon * sin(ph2(i,j,k))
               wr = umcon * cos(ph3(i,j,k))
               wi =-umcon * sin(ph3(i,j,k))

               tmp1(i,jj,kk) = cmplx(ur,ui)
               tmp2(i,jj,kk) = cmplx(vr,vi)
               tmp3(i,jj,kk) = cmplx(wr,wi)
               etot = etot + 0.5 * umcon**2

            end do
         end do
      end do

      do k=1,nz
      do j=1,ny
      do i=1,nx
         tmp1(i,j,k) = wtr(i,j,k)*tmp1(i,j,k)
         tmp2(i,j,k) = wtr(i,j,k)*tmp2(i,j,k)
         tmp3(i,j,k) = wtr(i,j,k)*tmp3(i,j,k)
      end do
      end do
      end do

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               exx(i,j,k) = tmp4(i,j,k)
            end do
         end do
      end do

c   Print out total energies

      write(6,71) enorm, etot
 71   format(' Exx total = ',e15.7,' Energy in U (Spectral)= ',e15.7)

c   Calculate velocity in physical space

      call ifft3d(tmp1, up, nx, ny, nz, tmp4, wk)
      call ifft3d(tmp2, vp, nx, ny, nz, tmp4, wk)
      call ifft3d(tmp3, wp, nx, ny, nz, tmp4, wk)

c     call fft3d(up, tmp1, nx, ny, nz, wk)
c     call ifft3d(tmp1, up, nx, ny, nz, tmp4, wk)
c   compute energy in physical space

      etot = 0.0

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
            etot = etot+.5*(up(i,j,k)**2+vp(i,j,k)**2+wp(i,j,k)**2)
            end do
         end do
      end do

      etot = etot / (float(nx)*float(ny)*float(nz))
      print*, ' Energy (physical)=', etot

      ratio = 3.*enorm / etot
      write(6,*) '   ratio = ',ratio

c   Read scaling factor
      write(6,*)
      write(6,*) ' Enter scale factor (sigma):'
      read(*,*)   sigma

c   Scale by SIGMA

      ratio = sqrt(ratio)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               up(i,j,k) = sigma*ratio * up(i,j,k)
               vp(i,j,k) = sigma*ratio * vp(i,j,k)
               wp(i,j,k) = sigma*ratio * wp(i,j,k)
            end do
         end do
      end do
 
      return
      end
c
      subroutine fft3d(v, vht, il, jl, kl, z2)

c  This routine performs an FFT on the real array v and places the
c  result in the complex array vht

c     implicit none

      integer id, il, jl, kl, i, j, k, ilm1, jlm1, klm1, nn(3)

      real tke
      real v(il,jl,kl), z2(2*il)

      complex vht(il,jl,kl)

c     Convert v to spectral space

      tke = 0.0

      do k = 1,kl
         do j = 1,jl
            do i = 1,il
               vht(i,j,k) = cmplx(v(i,j,k),0.0)
               tke = tke + 0.5 * v(i,j,k)**2
            end do
         end do
      end do

      tke = tke / float(il * jl * kl)

      print*, ' FFT3D internal check-'
      print*, ' TKEave (physical) =', tke

      nn(1) = il
      nn(2) = jl
      nn(3) = kl

      call fourt(vht,nn,3,-1,1,z2)

      do k = 1,kl
         do j = 1,jl
            do i = 1,il
               vht(i,j,k) = vht(i,j,k) / nn(1) / nn(2) / nn(3)
            end do
         end do
      end do

      tke = 0.0

      do k = 1,kl
         do j = 1,jl
            do i = 1,il
               tke = tke + 0.5 * real(vht(i,j,k) * conjg(vht(i,j,k)))
            end do
         end do
      end do

      print*, ' TKEave (spectral) =', tke

      return
      end
c
      subroutine ifft3d(vht, v, il, jl, kl, vhtmp, z2)

c  This routine performs an inverse 3-D FFT

      integer il, jl, kl, i, j, k, nn(3)

      real tke,tkep,ratio
      real v(il,jl,kl), z2(2*il)

      complex vht(il,jl,kl), vhtmp(il,jl,kl)

       tke = 0.0

       do k = 1,kl
          do j = 1,jl
             do i = 1,il
                tke = tke + 0.5 * real(vht(i,j,k) * conjg(vht(i,j,k)))
             end do
          end do
       end do

       print*, ' IFFT3D internal check-'
       print*, ' TKEave (spectral) =', tke

      nn(1) = il
      nn(2) = jl
      nn(3) = kl

      do k = 1,kl
         do j = 1,jl
            do i = 1,il
               vhtmp(i,j,k) = vht(i,j,k) * nn(1) * nn(2) * nn(3) 
            end do
         end do
      end do

      call fourt(vhtmp,nn,3,1,1,z2)
      
      tkep = 0.0

      do k = 1,kl
         do j = 1,jl
            do i = 1,il
               v(i,j,k) = real(vhtmp(i,j,k)) / nn(1) / nn(2) / nn(3)
               tkep = tkep + 0.5 * v(i,j,k)**2
            end do
         end do
      end do

       tkep = tkep / float(il * jl * kl)

       print*, ' TKEave (physical) =', tkep

      return
      end
c
      subroutine spectrum(vht, vt, il, jl, kl, num)

      implicit none

      integer
     .   kmax, il, jl, kl, kx, ky, kz, k, k1, ksum, ilm1, jlm1, klm1

       integer num(0:2*il)

      real
     .   rk, temp, vsum, rkx, rky, rkz, anum, etot

      real
     .   vt(0:2*il)

      complex
     .   vht(il,jl,kl)

      kmax = il / 2

      do k = 0,2*il
       vt(k) = 0.0
       num(k) = 0
      end do

      anum = 0.0
      etot = 0.0

      do kx = 1,il
         rkx = float(kx-1)
         if (rkx .gt. kmax) then
            rkx = kmax+1 - (rkx - kmax+1)
         endif
         do ky = 1,jl
            rky = float(ky-1)
            if (rky .gt. kmax) then
               rky = kmax+1 - (rky - kmax+1)
            endif
            do kz = 1,kl
               rkz = float(kz-1)
               if (rkz .gt. kmax) then
                  rkz = kmax+1 - (rkz - kmax+1)
               endif
               rk     = sqrt(rkx*rkx + rkz*rkz + rky*rky)
               k      = nint(rk)
               num(k) = num(k) + 1
               temp   = real(vht(kx,ky,kz) * conjg(vht(kx,ky,kz)))
               etot   = etot + sqrt(temp)
               vt(k)  = vt(k) + sqrt(temp)
c              anum   = anum + 0.5 * temp
c              vt(k)  = vt(k) + 0.5 * temp
            end do
         end do
      end do

      print*, ' '
      print*, ' Spectrum Internal Check-'
      print*, ' Total Energy in 3D field = ', etot
      print*, ' k, Num(k), vt(k)'
      ksum = 0
      vsum = 0.0
      do k = 0, il
       print*, k, num(k), vt(k)
       ksum = ksum + num(k)
       vsum = vsum + vt(k)
      end do

      print*, ' ksum: ', ksum
      print*, ' Total Energy in spectrum: ', vsum
      print*, ' '

      return
      end
c
      subroutine fourt(data,nn,ndim,isign,iform,work)
c
c     the cooley-tukey fast fourier transform in usasi basic fortran
c     transform(j1,j2,,,,) = sum(data(i1,i2,,,,)*w1**((i2-1)*(j2-1))
c                                 *w2**((i2-1)*(j2-1))*,,,),
c     where i1 and j1 run from 1 to nn(1) and w1=exp(isign*2*pi=
c     sqrt(-1)/nn(1)), etc.  there is no limit on the dimensionality
c     (number of subscripts) of the data array.  if an inverse
c     transform (isign=+1) is performed upon an array of transformed
c     (isign=-1) data, the original data will reappear.
c     multiplied by nn(1)*nn(2)*,,,  the array of input data must be
c     in complex format.  however, if all imaginary parts are zero (i.e.
c     the data are disguised real) running time is cut up to forty per-
c     cent.  (for fastest transform of real data, nn(1) should be even.)
c     the transform values are always complex and are returned in the
c     original array of data, replacing the input data.  the length
c     of each dimension of the data array may be any integer.  the
c     program runs faster on composite integers than on primes, and is
c     particularly fast on numbers rich in factors of two.
c
c     timing is in fact given by the following formula.  let ntot be the
c     total number of points (real or complex) in the data array, that
c     is, ntot=nn(1)*nn(2)*...  decompose ntot into its prime factors,
c     such as 2**k2 * 3**k3 * 5**k5 * ...  let sum2 be the sum of all
c     the factors of two in ntot, that is, sum2 = 2*k2.  let sumf be
c     the sum of all other factors of ntot, that is, sumf = 3*k3*5*k5*..
c     the time taken by a multidimensional transform on these ntot data
c     is t = t0 + ntot*(t1+t2*sum2+t3*sumf).  on the cdc 3300 (floating
c     point add time = six microseconds), t = 3000 + ntot*(600+40*sum2+
c     175*sumf) microseconds on complex data.
c
c     implementation of the definition by summation will run in a time
c     proportional to ntot*(nn(1)+nn(2)+...).  for highly composite ntot
c     the savings offered by this program can be dramatic.  a one-dimen-
c     sional array 4000 in length will be transformed in 4000*(600+
c     40*(2+2+2+2+2)+175*(5+5+5)) = 14.5 seconds versus about 4000*
c     4000*175 = 2800 seconds for the straightforward technique.
c
c     the fast fourier transform places three restrictions upon the
c     data.
c     1.  the number of input data and the number of transform values
c     must be the same.
c     2.  both the input data and the transform values must represent
c     equispaced points in their respective domains of time and
c     frequency.  calling these spacings deltat and deltaf, it must be
c     true that deltaf=2*pi/(nn(i)*deltat).  of course, deltat need not
c     be the same for every dimension.
c     3.  conceptually at least, the input data and the transform output
c     represent single cycles of periodic functions.
c
c     the calling sequence is--
c     call fourt(data,nn,ndim,isign,iform,work)
c
c     data is the array used to hold the real and imaginary parts
c     of the data on input and the transform values on output.  it
c     is a multidimensional floating point array, with the real and
c     imaginary parts of a datum stored immediately adjacent in storage
c     (such as fortran iv places them).  normal fortran ordering is
c     expected, the first subscript changing fastest.  the dimensions
c     are given in the integer array nn, of length ndim.  isign is -1
c     to indicate a forward transform (exponential sign is -) and +1
c     for an inverse transform (sign is +).  iform is +1 if the data are
c     complex, 0 if the data are real.  if it is 0, the imaginary
c     parts of the data must be set to zero.  as explained above, the
c     transform values are always complex and are stored in array data.
c     work is an array used for working storage.  it is floating point
c     real, one dimensional of length equal to twice the largest array
c     dimension nn(i) that is not a power of two.  if all nn(i) are
c     powers of two, it is not needed and may be replaced by zero in the
c     calling sequence.  thus, for a one-dimensional array, nn(1) odd,
c     work occupies as many storage locations as data.  if supplied,
c     work must not be the same array as data.  all subscripts of all
c     arrays begin at one.
c
c     example 1.  three-dimensional forward fourier transform of a
c     complex array dimensioned 32 by 25 by 13 in fortran iv.
c     dimension data(32,25,13),work(50),nn(3)
c     complex data
c     data nn/32,25,13/
c     do 1 i=1,32
c     do 1 j=1,25
c     do 1 k=1,13
c  1  data(i,j,k)=complex value
c     call fourt(data,nn,3,-1,1,work)
c
c     example 2.  one-dimensional forward transform of a real array of
c     length 64 in fortran ii,
c     dimension data(2,64)
c     do 2 i=1,64
c     data(1,i)=real part
c  2  data(2,i)=0.
c     call fourt(data,64,1,-1,0,0)
c
c     there are no error messages or error halts in this program.  the
c     program returns immediately if ndim or any nn(i) is less than one.
c
c     program by norman brenner from the basic program by charles
c     rader,  june 1967.  the idea for the digit reversal was
c     suggested by ralph alter.
c
c     this is the fastest and most versatile version of the fft known
c     to the author.  a program called four2 is available that also
c     performs the fast fourier transform and is written in usasi basic
c     fortran.  it is about one third as long and restricts the
c     dimensions of the input array (which must be complex) to be powers
c     of two.  another program, called four1, is one tenth as long and
c     runs two thirds as fast on a one-dimensional complex array whose
c     length is a power of two.
c
c     reference--
c     ieee audio transactions (june 1967), special issue on the fft.

      dimension data(*),nn(1),ifact(32),work(1)

      data twopi/6.2831853071796/,rthlf/0.70710678118655/
      data nprev/0/,np0/0/
c the following call is for gathering statistics on library use at ncar
c     call q8qst4( 4hxlib      , 5hfourt     ,5hfourt  ,10hversion  9)
      if(ndim-1)920,1,1
1     ntot=2
      do 2 idim=1,ndim
      if(nn(idim))920,920,2
2     ntot=ntot*nn(idim)
c
c     main loop for each dimension
c
      np1=2
      do 910 idim=1,ndim
      n=nn(idim)
      np2=np1*n
      if(n-1)920,900,5
c
c     is n a power of two and if not, what are its factors
c
5     m=n
      ntwo=np1
      if=1
      idiv=2
10    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv)50,11,11
11    if(irem)20,12,20
12    ntwo=ntwo+ntwo
      ifact(if)=idiv
      if=if+1
      m=iquot
      go to 10
20    idiv=3
      inon2=if
30    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv)60,31,31
31    if(irem)40,32,40
32    ifact(if)=idiv
      if=if+1
      m=iquot
      go to 30
40    idiv=idiv+2
      go to 30
50    inon2=if
      if(irem)60,51,60
51    ntwo=ntwo+ntwo
      go to 70
60    ifact(if)=m
c
c     separate four cases--
c        1. complex transform or real transform for the 4th, 9th,etc.
c           dimensions.
c        2. real transform for the 2nd or 3rd dimension.  method--
c           transform half the data, supplying the other half by con-
c           jugate symmetry.
c        3. real transform for the 1st dimension, n odd.  method--
c           set the imaginary parts to zero.
c        4. real transform for the 1st dimension, n even.  method--
c           transform a complex array of length n/2 whose real parts
c           are the even numbered real values and whose imaginary parts
c           are the odd numbered real values.  separate and supply
c           the second half by conjugate symmetry.
c
70    icase=1
      ifmin=1
      i1rng=np1
      if(idim-4)71,100,100
71    if(iform)72,72,100
72    icase=2
      i1rng=np0*(1+nprev/2)
      if(idim-1)73,73,100
73    icase=3
      i1rng=np1
      if(ntwo-np1)100,100,74
74    icase=4
      ifmin=2
      ntwo=ntwo/2
      n=n/2
      np2=np2/2
      ntot=ntot/2
      i=1
      do 80 j=1,ntot
      data(j)=data(i)
80    i=i+2
c
c     shuffle data by bit reversal, since n=2**k.  as the shuffling
c     can be done by simple interchange, no working array is needed
c
100   if(ntwo-np2)200,110,110
110   np2hf=np2/2
      j=1
      do 150 i2=1,np2,np1
      if(j-i2)120,130,130
120   i1max=i2+np1-2
      do 125 i1=i2,i1max,2
      do 125 i3=i1,ntot,np2
      j3=j+i3-i2
      tempr=data(i3)
      tempi=data(i3+1)
      data(i3)=data(j3)
      data(i3+1)=data(j3+1)
      data(j3)=tempr
125   data(j3+1)=tempi
130   m=np2hf
140   if(j-m)150,150,145
145   j=j-m
      m=m/2
      if(m-np1)150,140,140
150   j=j+m
      go to 300
c
c     shuffle data by digit reversal for general n
c
200   nwork=2*n
      do 270 i1=1,np1,2
      do 270 i3=i1,ntot,np2
      j=i3
      do 260 i=1,nwork,2
      if(icase-3)210,220,210
210   work(i)=data(j)
      work(i+1)=data(j+1)
      go to 230
220   work(i)=data(j)
      work(i+1)=0.
230   ifp2=np2
      if=ifmin
240   ifp1=ifp2/ifact(if)
      j=j+ifp1
      if(j-i3-ifp2)260,250,250
250   j=j-ifp2
      ifp2=ifp1
      if=if+1
      if(ifp2-np1)260,260,240
260   continue
      i2max=i3+np2-np1
      i=1
      do 270 i2=i3,i2max,np1
      data(i2)=work(i)
      data(i2+1)=work(i+1)
270   i=i+2
c
c     main loop for factors of two.  perform fourier transforms of
c     length four, with one of length two if needed.  the twiddle factor
c     w=exp(isign*2*pi*sqrt(-1)*m/(4*mmax)).  check for w=isign*sqrt(-1)
c     and repeat for w=w*(1+isign*sqrt(-1))/sqrt(2).
c
300   if(ntwo-np1)600,600,305
305   np1tw=np1+np1
      ipar=ntwo/np1
310   if(ipar-2)350,330,320
320   ipar=ipar/4
      go to 310
330   do 340 i1=1,i1rng,2
      do 340 k1=i1,ntot,np1tw
      k2=k1+np1
      tempr=data(k2)
      tempi=data(k2+1)
      data(k2)=data(k1)-tempr
      data(k2+1)=data(k1+1)-tempi
      data(k1)=data(k1)+tempr
340   data(k1+1)=data(k1+1)+tempi
350   mmax=np1
360   if(mmax-ntwo/2)370,600,600
370   lmax=max0(np1tw,mmax/2)
      do 570 l=np1,lmax,np1tw
      m=l
      if(mmax-np1)420,420,380
380   theta=-twopi*float(l)/float(4*mmax)
      if(isign)400,390,390
390   theta=-theta
400   wr=cos(theta)
      wi=sin(theta)
410   w2r=wr*wr-wi*wi
      w2i=2.*wr*wi
      w3r=w2r*wr-w2i*wi
      w3i=w2r*wi+w2i*wr
420   do 530 i1=1,i1rng,2
      kmin=i1+ipar*m
      if(mmax-np1)430,430,440
430   kmin=i1
440   kdif=ipar*mmax
450   kstep=4*kdif
      if(kstep-ntwo)460,460,530
460   do 520 k1=kmin,ntot,kstep
      k2=k1+kdif
      k3=k2+kdif
      k4=k3+kdif
      if(mmax-np1)470,470,480
470   u1r=data(k1)+data(k2)
      u1i=data(k1+1)+data(k2+1)
      u2r=data(k3)+data(k4)
      u2i=data(k3+1)+data(k4+1)
      u3r=data(k1)-data(k2)
      u3i=data(k1+1)-data(k2+1)
      if(isign)471,472,472
471   u4r=data(k3+1)-data(k4+1)
      u4i=data(k4)-data(k3)
      go to 510
472   u4r=data(k4+1)-data(k3+1)
      u4i=data(k3)-data(k4)
      go to 510
480   t2r=w2r*data(k2)-w2i*data(k2+1)
      t2i=w2r*data(k2+1)+w2i*data(k2)
      t3r=wr*data(k3)-wi*data(k3+1)
      t3i=wr*data(k3+1)+wi*data(k3)
      t4r=w3r*data(k4)-w3i*data(k4+1)
      t4i=w3r*data(k4+1)+w3i*data(k4)
      u1r=data(k1)+t2r
      u1i=data(k1+1)+t2i
      u2r=t3r+t4r
      u2i=t3i+t4i
      u3r=data(k1)-t2r
      u3i=data(k1+1)-t2i
      if(isign)490,500,500
490   u4r=t3i-t4i
      u4i=t4r-t3r
      go to 510
500   u4r=t4i-t3i
      u4i=t3r-t4r
510   data(k1)=u1r+u2r
      data(k1+1)=u1i+u2i
      data(k2)=u3r+u4r
      data(k2+1)=u3i+u4i
      data(k3)=u1r-u2r
      data(k3+1)=u1i-u2i
      data(k4)=u3r-u4r
520   data(k4+1)=u3i-u4i
      kdif=kstep
      kmin=4*(kmin-i1)+i1
      go to 450
530   continue
      m=m+lmax
      if(m-mmax)540,540,570
540   if(isign)550,560,560
550   tempr=wr
      wr=(wr+wi)*rthlf
      wi=(wi-tempr)*rthlf
      go to 410
560   tempr=wr
      wr=(wr-wi)*rthlf
      wi=(tempr+wi)*rthlf
      go to 410
570   continue
      ipar=3-ipar
      mmax=mmax+mmax
      go to 360
c
c     main loop for factors not equal to two.  apply the twiddle factor
c     w=exp(isign*2*pi*sqrt(-1)*(j1-1)*(j2-j1)/(ifp1+ifp2)), then
c     perform a fourier transform of length ifact(if), making use of
c     conjugate symmetries.
c
600   if(ntwo-np2)605,700,700
605   ifp1=ntwo
      if=inon2
      np1hf=np1/2
610   ifp2=ifact(if)*ifp1
      j1min=np1+1
      if(j1min-ifp1)615,615,640
615   do 635 j1=j1min,ifp1,np1
      theta=-twopi*float(j1-1)/float(ifp2)
      if(isign)625,620,620
620   theta=-theta
625   wstpr=cos(theta)
      wstpi=sin(theta)
      wr=wstpr
      wi=wstpi
      j2min=j1+ifp1
      j2max=j1+ifp2-ifp1
      do 635 j2=j2min,j2max,ifp1
      i1max=j2+i1rng-2
      do 630 i1=j2,i1max,2
      do 630 j3=i1,ntot,ifp2
      tempr=data(j3)
      data(j3)=data(j3)*wr-data(j3+1)*wi
630   data(j3+1)=tempr*wi+data(j3+1)*wr
      tempr=wr
      wr=wr*wstpr-wi*wstpi
635   wi=tempr*wstpi+wi*wstpr
640   theta=-twopi/float(ifact(if))
      if(isign)650,645,645
645   theta=-theta
650   wstpr=cos(theta)
      wstpi=sin(theta)
      j2rng=ifp1*(1+ifact(if)/2)
      do 695 i1=1,i1rng,2
      do 695 i3=i1,ntot,np2
      j2max=i3+j2rng-ifp1
      do 690 j2=i3,j2max,ifp1
      j1max=j2+ifp1-np1
      do 680 j1=j2,j1max,np1
      j3max=j1+np2-ifp2
      do 680 j3=j1,j3max,ifp2
      jmin=j3-j2+i3
      jmax=jmin+ifp2-ifp1
      i=1+(j3-i3)/np1hf
      if(j2-i3)655,655,665
655   sumr=0.
      sumi=0.
      do 660 j=jmin,jmax,ifp1
 659  sumr=sumr+data(j)
660   sumi=sumi+data(j+1)
      work(i)=sumr
      work(i+1)=sumi
      go to 680
665   iconj=1+(ifp2-2*j2+i3+j3)/np1hf
      j=jmax
      sumr=data(j)
      sumi=data(j+1)
      oldsr=0.
      oldsi=0.
      j=j-ifp1
670   tempr=sumr
      tempi=sumi
      sumr=twowr*sumr-oldsr+data(j)
      sumi=twowr*sumi-oldsi+data(j+1)
      oldsr=tempr
      oldsi=tempi
      j=j-ifp1
      if(j-jmin)675,675,670
675   tempr=wr*sumr-oldsr+data(j)
      tempi=wi*sumi
      work(i)=tempr-tempi
      work(iconj)=tempr+tempi
      tempr=wr*sumi-oldsi+data(j+1)
      tempi=wi*sumr
      work(i+1)=tempr+tempi
      work(iconj+1)=tempr-tempi
680   continue
      if(j2-i3)685,685,686
685   wr=wstpr
      wi=wstpi
      go to 690
686   tempr=wr
      wr=wr*wstpr-wi*wstpi
      wi=tempr*wstpi+wi*wstpr
690   twowr=wr+wr
      i=1
      i2max=i3+np2-np1
      do 695 i2=i3,i2max,np1
      data(i2)=work(i)
      data(i2+1)=work(i+1)
695   i=i+2
      if=if+1
      ifp1=ifp2
      if(ifp1-np2)610,700,700
c
c     complete a real transform in the 1st dimension, n even, by con-
c     jugate symmetries.
c
700   go to (900,800,900,701),icase
701   nhalf=n
      n=n+n
      theta=-twopi/float(n)
      if(isign)703,702,702
702   theta=-theta
703   wstpr=cos(theta)
      wstpi=sin(theta)
      wr=wstpr
      wi=wstpi
      imin=3
      jmin=2*nhalf-1
      go to 725
710   j=jmin
      do 720 i=imin,ntot,np2
      sumr=(data(i)+data(j))/2.
      sumi=(data(i+1)+data(j+1))/2.
      difr=(data(i)-data(j))/2.
      difi=(data(i+1)-data(j+1))/2.
      tempr=wr*sumi+wi*difr
      tempi=wi*sumi-wr*difr
      data(i)=sumr+tempr
      data(i+1)=difi+tempi
      data(j)=sumr-tempr
      data(j+1)=-difi+tempi
720   j=j+np2
      imin=imin+2
      jmin=jmin-2
      tempr=wr
      wr=wr*wstpr-wi*wstpi
      wi=tempr*wstpi+wi*wstpr
725   if(imin-jmin)710,730,740
730   if(isign)731,740,740
731   do 735 i=imin,ntot,np2
735   data(i+1)=-data(i+1)
740   np2=np2+np2
      ntot=ntot+ntot
      j=ntot+1
      imax=ntot/2+1
745   imin=imax-2*nhalf
      i=imin
      go to 755
750   data(j)=data(i)
      data(j+1)=-data(i+1)
755   i=i+2
      j=j-2
      if(i-imax)750,760,760
760   data(j)=data(imin)-data(imin+1)
      data(j+1)=0.
      if(i-j)770,780,780
765   data(j)=data(i)
      data(j+1)=data(i+1)
770   i=i-2
      j=j-2
      if(i-imin)775,775,765
775   data(j)=data(imin)+data(imin+1)
      data(j+1)=0.
      imax=imin
      go to 745
780   data(1)=data(1)+data(2)
      data(2)=0.
      go to 900
c
c     complete a real transform for the 2nd or 3rd dimension by
c     conjugate symmetries.
c
800   if(i1rng-np1)805,900,900
805   do 860 i3=1,ntot,np2
      i2max=i3+np2-np1
      do 860 i2=i3,i2max,np1
      imin=i2+i1rng
      imax=i2+np1-2
      jmax=2*i3+np1-imin
      if(i2-i3)820,820,810
810   jmax=jmax+np2
820   if(idim-2)850,850,830
830   j=jmax+np0
      do 840 i=imin,imax,2
      data(i)=data(j)
      data(i+1)=-data(j+1)
840   j=j-2
850   j=jmax
      do 860 i=imin,imax,np0
      data(i)=data(j)
      data(i+1)=-data(j+1)
860   j=j-np0
c
c     end of loop on each dimension
c
900   np0=np1
      np1=np2
910   nprev=n
920   return
      end

