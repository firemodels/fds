      MODULE POIS
C
C Poisson Solver Routines
C
      USE PREC
      IMPLICIT REAL(EB) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PRIVATE
      REAL(EB) scale
      INTEGER kappa,nmax,ikpwr
      LOGICAL tpose,nocopy,outary
C
      PUBLIC H3CZIS,H3CZSS,H2CZSS,H2CYSS,H3CSSS,H2CZIS,H3CSIS,H2CYIS
C
C
      CONTAINS
C
C
      SUBROUTINE h3czis(xs,xf,l,lbdcnd,ys,yf,m,mbdcnd,zs,zf,n,nbdcnd,
     *                  h,elmbda,ldimf,mdimf,ierror,save)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB) save(-3:*),h(0:*)

C                               CHECK FOR INVALID INPUT

      ierror = 0

      IF (xs.GT.xf) THEN
          ierror = ierror + 1
          save(ierror) = 1.
      END IF

      IF (l.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 2.
      END IF

      IF (lbdcnd.LT.0 .OR. lbdcnd.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 3.
      END IF

      IF (ys.GT.yf) THEN
          ierror = ierror + 1
          save(ierror) = 4.
      END IF

      IF (m.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 5.
      END IF

      IF (mbdcnd.LT.0 .OR. mbdcnd.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 6.
      END IF

      IF (zs.GT.zf) THEN
          ierror = ierror + 1
          save(ierror) = 7.
      END IF

      IF (n.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 8.
      END IF

      IF (nbdcnd.LT.0 .OR. nbdcnd.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 9.
      END IF

      IF (ldimf.LT.l) THEN
          ierror = ierror + 1
          save(ierror) = 10.
      END IF

      IF (mdimf.LT.m) THEN
          ierror = ierror + 1
          save(ierror) = 11.
      END IF

      IF (ierror.NE.0) THEN
          RETURN
      ELSE
          save(1) = 0.
      END IF

C                               DEFINE GRID PARAMETERS

      dx = (xf-xs)/l
      dy = (yf-ys)/m
      dz = (zf-zs)/n
      dlxsqr = 1./dx**2
      dlysqr = 1./dy**2
      dlzsqr = 1./dz**2
      lp = lbdcnd + 1
      mp = mbdcnd + 1
      np = nbdcnd + 1

C                               ALLOCATE SAVE ARRAY

      ia = 12
      ib = ia + l
      ic = ib + l
      id = ic + l
      is = id + l

C                               DEFINE THE A,B,C COEFFICIENTS
C                               IN SAVE ARRAY

      DO 100 i = 0,l - 1
          hm = .5*(h(i)+h(i+1))
          hp = .5*(h(i+1)+h(i+2))
          save(ia+i) = dlxsqr/(h(i+1)*hm)
          save(ic+i) = dlxsqr/(h(i+1)*hp)
          save(ib+i) = - (save(ia+i)+save(ic+i)) + elmbda
          save(id+i) = dlysqr
  100 CONTINUE

      select case(lp)
      case(2)
      save(ib) = save(ib) - save(ia)
      save(ic-1) = save(ic-1) - save(id-1)
      case(3)
      save(ib) = save(ib) - save(ia)
      save(ic-1) = save(ic-1) + save(id-1)
      case(4)
      save(ib) = save(ib) + save(ia)
      save(ic-1) = save(ic-1) + save(id-1)
      case(5)
      save(ib) = save(ib) + save(ia)
      save(ic-1) = save(ic-1) - save(id-1)
      end select


C                               DETERMINE WHETHER OR NOT BOUNDARY
C                               CONDITION IS PERIODIC IN X

      IF (lbdcnd.EQ.0) THEN
          lperod = 0
      ELSE
          lperod = 1
      END IF

C                               INITIALIZE SOLVER ROUTINE S3CFIS.
C
      CALL s3cfis(lperod,l,mbdcnd,m,nbdcnd,n,dlzsqr,save(ia),save(ib),
     *            save(ic),save(id),ldimf,mdimf,ir,save(is))

C                               TEST ERROR FLAG FROM S3CFIS FOR
C                               INTERNAL ERROR

      IF (ir.NE.0) THEN
          save(1) = 99.
          ierror = 1

          RETURN
      END IF

C                               SAVE PARAMETERS FOR H3CZSS IN SAVE ARRAY

      save(2) = dx
      save(3) = l
      save(4) = lp
      save(5) = dy
      save(6) = m
      save(7) = mp
      save(8) = dz
      save(9) = n
      save(10) = np
      save(11) = elmbda

      RETURN
      END SUBROUTINE h3czis

      SUBROUTINE h3czss(bdxs,bdxf,bdys,bdyf,bdzs,bdzf,ldimf,mdimf,f,
     *                  pertrb,save,w,h)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB) bdxs(mdimf,*),bdxf(mdimf,*),bdys(ldimf,*),bdyf(ldimf,*),
     *         bdzs(ldimf,*),bdzf(ldimf,*),f(ldimf,mdimf,*),save(-3:*),
     *         w(*),h(0:*)
C
C                               CHECK VALUE OF IERROR (=SAVE(1)).
C                               IF NON-ZERO, RETURN.

      IF (save(1).NE.0.) RETURN

C                               GET PARAMETERS FOR H3CZSS FROM SAVE
C                               ARRAY WHERE THEY WERE STORED IN
C                               INITIALIZATION SUBROUTINE H3CZIS.

      dx = save(2)
      l = save(3)
      lp = save(4)
      dy = save(5)
      m = save(6)
      mp = save(7)
      dz = save(8)
      n = save(9)
      np = save(10)
      elmbda = save(11)

      dlxrcp = 1./dx
      twdxsq = 2./dx**2
      dlyrcp = 1./dy
      twdysq = 2./dy**2
      dlzrcp = 1./dz
      twdzsq = 2./dz**2

C                               ALLOCATE SAVE ARRAY

      ia = 12
      ib = ia + l
      ic = ib + l
      id = ic + l
      is = 12 + 4*l


C                               ENTER BOUNDARY DATA FOR X-BOUNDARIES

      if (lp.eq.2 .or. lp.eq.3) then
      do k = 1,n
      do j = 1,m
      f(1,j,k) = f(1,j,k) - 2.*bdxs(j,k)*save(ia)
      enddo
      enddo
      endif

      if (lp.eq.4 .or. lp.eq.5) then
      do k = 1,n
      do j = 1,m
      f(1,j,k) = f(1,j,k) + save(ia)*dx*bdxs(j,k)
      enddo
      enddo
      endif

      if (lp.eq.2 .or. lp.eq.5) then
      do k = 1,n
      do j = 1,m
      f(l,j,k) = f(l,j,k) - 2.*bdxf(j,k)*save(id-1)
      enddo
      enddo
      endif

      if (lp.eq.3 .or. lp.eq.4) then
      do k = 1,n
      do j = 1,m
      f(l,j,k) = f(l,j,k) - save(id-1)*dx*bdxf(j,k)
      enddo
      enddo
      endif
 
C                               ENTER BOUNDARY DATA FOR Y-BOUNDARIES

      if (mp.eq.2 .or. mp.eq.3) then
      DO k = 1,n
      DO i = 1,l
      f(i,1,k) = f(i,1,k) - bdys(i,k)*twdysq
      enddo
      enddo
      endif

      if (mp.eq.4 .or. mp.eq.5) then
      DO k = 1,n
      DO i = 1,l
      f(i,1,k) = f(i,1,k) + bdys(i,k)*dlyrcp
      enddo
      enddo
      endif

      if (mp.eq.2 .or. mp.eq.5) then
      DO k = 1,n
      DO i = 1,l
      f(i,m,k) = f(i,m,k) - bdyf(i,k)*twdysq
      enddo
      enddo
      endif

      if (mp.eq.3 .or. mp.eq.4) then
      DO k = 1,n
      DO i = 1,l
      f(i,m,k) = f(i,m,k) - bdyf(i,k)*dlyrcp
      enddo
      enddo
      endif

C                               ENTER BOUNDARY DATA FOR Z-BOUNDARIES

      if (np.eq.2 .or. np.eq.3) then
      DO j = 1,m
      DO i = 1,l
      f(i,j,1) = f(i,j,1) - bdzs(i,j)*twdzsq
      enddo
      enddo
      endif

      if (np.eq.4 .or. np.eq.5) then
      DO j = 1,m
      DO i = 1,l
      f(i,j,1) = f(i,j,1) + bdzs(i,j)*dlzrcp
      enddo
      enddo
      endif

      if (np.eq.2 .or. np.eq.5) then
      DO j = 1,m
      DO i = 1,l
      f(i,j,n) = f(i,j,n) - bdzf(i,j)*twdzsq
      enddo
      enddo
      endif

      if (np.eq.3 .or. np.eq.4) then
      DO j = 1,m
      DO i = 1,l
      f(i,j,n) = f(i,j,n) - bdzf(i,j)*dlzrcp
      enddo
      enddo
      endif

      pertrb = 0.
      pert   = 0.
      ising = 0

C                               FOR SINGULAR PROBLEMS ADJUST DATA TO
C                               INSURE A SOLUTION WILL EXIST.  GO THRU
C                               THIS CODE TWICE: ISING=1 FOR CALCULATING
C                               PERTRB; ISING=2 FOR NORMALIZING SOLUTION
C                               AFTER IT IS COMPUTED.

      SELECT CASE(lp)
      CASE(1)   ; GOTO 630
      CASE(2:3) ; GOTO 750
      CASE(4)   ; GOTO 630
      CASE(5)   ; GOTO 750
      END SELECT

  630 CONTINUE
      SELECT CASE(mp)
      CASE(1)   ; GOTO 640
      CASE(2:3) ; GOTO 750
      CASE(4)   ; GOTO 640
      CASE(5)   ; GOTO 750
      END SELECT

  640 CONTINUE
      SELECT CASE(np)
      CASE(1)   ; GOTO 650
      CASE(2:3) ; GOTO 750
      CASE(4)   ; GOTO 650
      CASE(5)   ; GOTO 750
      END SELECT

  650 CONTINUE
      IF (elmbda.NE.0.) GO TO 750
      ising = 1
  660 CONTINUE
      pert = 0.
      DO 670 i = 1,l
          w(i) = 0.
  670 CONTINUE
      DO 700 k = 1,n
          DO 690 j = 1,m
              DO 680 i = 1,l
                  w(i) = w(i) + f(i,j,k)
  680         CONTINUE
  690     CONTINUE
  700 CONTINUE
      s1 = 0.
      s3 = 0.
      DO 710 i = 1,l
          s3 = s3 + h(i)
          s1 = s1 + h(i)*w(i)
  710 CONTINUE

      s3   = s3*m*n
      pert = s1/s3


C                               ADJUST F ARRAY BY PERT

      DO 740 k = 1,n
          DO 730 j = 1,m
              DO 720 i = 1,l
                  f(i,j,k) = f(i,j,k) - pert
  720         CONTINUE
  730     CONTINUE
  740 CONTINUE
  750 CONTINUE

C                               IF NORMALIZING SOLUTION, RESTORE PERTRB
C                               AND JUMP TO END

      IF (ising.EQ.2) THEN
          pertrb = prtsav

          GO TO 800

      END IF

      prtsav = pert

C                               SOLVE THE EQUATION

      CALL s3cfss(ldimf,mdimf,f,save(is),w)

C                               IF A SINGULAR PROBLEM,
C                               RE-NORMALIZE SOLUTION (ISING=2)

      IF (ising.EQ.1) THEN
          ising = 2

          GO TO 660

      END IF

  800 CONTINUE
      RETURN
      END SUBROUTINE h3czss
C
C
      SUBROUTINE s3cfis(lperod,l,mperod,m,nperod,n,scal,a,b,c,d,ldimf,
     *                  mdimf,ierror,save)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   a(l),b(l),c(l),d(l),save(-3:*)

C                               CHECK FOR INVALID INPUT

      ierror = 0

      IF (lperod.NE.0 .AND. lperod.NE.1) THEN
          ierror = 1
          save(1) = 1.
      END IF

      IF (l.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 2.
      END IF

      IF (mperod.LT.0 .AND. mperod.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 3.
      END IF

      IF (m.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 4.
      END IF

      IF (nperod.LT.0 .AND. nperod.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 5.
      END IF

      IF (n.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 6.
      END IF

      IF (ldimf.LT.l) THEN
          ierror = ierror + 1
          save(ierror) = 7.
      END IF

      IF (mdimf.LT.m) THEN
          ierror = ierror + 1
          save(ierror) = 8.
      END IF

      IF (lperod.EQ.0) THEN
          DO 100 i = 1,l

              IF (a(i).NE.a(1)) GO TO 110
              IF (b(i).NE.b(1)) GO TO 110
              IF (c(i).NE.a(1)) GO TO 110
              IF (d(i).NE.d(1)) GO TO 110

  100     CONTINUE

          GO TO 120

  110     CONTINUE
          ierror = ierror + 1
          save(ierror) = 9.
      END IF

  120 CONTINUE

      IF (ierror.NE.0) THEN
          RETURN
      ELSE
          save(1) = ierror
      END IF

      ldimfc=l
      if (ldimf.GT.l .AND. MOD(l,2).EQ.0) ldimfc=l+1

      igrid = 2

      CALL fsh00s(igrid,lperod,l,mperod,m,nperod,n,ldimfc,
     *            scal,a,b,c,d,save)

      RETURN
      END SUBROUTINE s3cfis
C
C
      SUBROUTINE s3cfss(ldimf,mdimf,f,save,w)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   f(ldimf,mdimf,*),save(-3:*),w(*)

C                               CHECK VALUE OF IERROR (=SAVE(1)).
C                               IF NON-ZERO, RETURN.

      IF (save(1).NE.0.) RETURN

      CALL fsh02s(ldimf,mdimf,f,save,w)

      RETURN
      END SUBROUTINE s3cfss
C
C
      SUBROUTINE fsh00s(igrid,lperod,l,mperod,m,nperod,n,ldimfc,c2,
     *                  a,b,c,d,save)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   a(l),b(l),c(l),d(l),save(-3:*)

C                               THIS SUBROUTINE INITIALIZES FFT SOLVER

C                               ALLOCATE SAVE ARRAY

      ia = 12
      ib = ia + l
      ic = ib + l
      icfy = ic + l

      IF (igrid.EQ.1) THEN
          icfz = icfy + 2*m
          ifctrd = icfz + 2*n
      ELSE
          icfz = icfy + 4*m
          ifctrd = icfz + 4*n
      END IF

      iwsy = ifctrd + ldimfc*m*n
      iwsz = iwsy + m + 16
      izrt = iwsz + n + 16

C                               COPY COEFFICIENT ARRAYS A,B, AND C INTO
C                               SAVE ARRAY.  A COPY OF B IS MADE BECAUSE
C                               IN THE NEXT LEVEL ROUTINE, BOUNDARY
C                               ELEMENTS OF B MAY BE CHANGED.

      DO 100 i = 0,l - 1
          save(ia+i) = a(i+1)
          save(ib+i) = b(i+1)
          save(ic+i) = c(i+1)
  100 CONTINUE
      lp = lperod + 1
      mp = mperod + 1
      np = nperod + 1

C                               CALL LOWER LEVEL INITIALIZATION ROUTINE

      CALL fsh01s(igrid,l,lp,m,mp,d,n,np,ldimfc,c2,save(ia),save(ib),
     *            save(ic),save(icfy),save(icfz),save(ifctrd),
     *            save(iwsy),save(iwsz),save(izrt))

C                               SAVE PARAMETERS FOR SUBROUTINE SOLVER

      save(2) = l
      save(3) = lp
      save(4) = m
      save(5) = mp
      save(6) = n
      save(7) = np
      save(8) = icfz
      save(9) = iwsz
      save(10) = izrt
      save(11) = igrid

      RETURN
      END SUBROUTINE fsh00s
C
C
      SUBROUTINE fsh01s(igrid,l,lp,m,mp,d,n,np,ldimfc,c2,a,b,c,cfy,cfz,
     *                  fctrd,wsavey,wsavez,zrt)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB) a(l),b(l),c(l),d(l),cfy(4*m),cfz(4*n),fctrd(ldimfc,n,m),
     *         wsavey(m+16),wsavez(n+16),zrt(n)

C                               INITIALIZATION ROUTINE FOR FFT SOLVERS

      pi = fsh19s()
      eps = fsh20s()

      IF (l.GT.1 .AND. lp.EQ.1) THEN
          lh = (l+1)/2
          lodd = 1

          IF (2*lh.EQ.l) lodd = 2
          c(lh-1) = 0.
          a(lh) = 0.
          c(lh) = 2.*c(lh)

       IF (lodd.eq.1) then
          b(lh-1) = b(lh-1) - a(lh-1)
          b(l) = b(l) + a(l)
          endif

          if (lodd.eq.2) a(l) = c(lh)
      END IF

C                               COMPUTE TRANSFORM ROOTS FOR J-DIRECTION

      IF (m.EQ.1) THEN
          cfy(1) = 0.
      ELSE

          IF (igrid.EQ.1) THEN
              mrdel = ((mp-1)* (mp-3)* (mp-5))/3
              del = pi/ (2.* (m+mrdel))
          ELSE
              del = pi/ (2*m)
          END IF

          if (mp.eq.1) then
          cfy(1) = 0.
          cfy(m) = -4.
          DO j = 2,m - 1,2
          cfy(j) = -4.*sin(j*del)**2
          cfy(j+1) = cfy(j)
          enddo
          endif

          if (mp.eq.2) then
          DO j = 1,m
          cfy(j) = -4.*sin(j*del)**2
          enddo
          endif

          if (mp.eq.3 .or. mp.eq.5) then
          DO j = 1,m
          cfy(j) = -4.*sin((j-.5)*del)**2
          enddo
          endif

          if (mp.eq.4) then
          DO j = 1,m
          cfy(j) = -4.*sin((j-1)*del)**2
          enddo
          endif

      END IF

C                               COMPUTE TRANSFORM ROOTS IN K-DIRECTION

      IF (n.EQ.1) THEN
          zrt(1) = 0.
      ELSE

          IF (igrid.EQ.1) THEN
              nrdel = ((np-1)* (np-3)* (np-5))/3
              del = pi/ (2.* (n+nrdel))
          ELSE
              del = pi/ (2*n)
          END IF

          if (np.eq.1) then
          zrt(1) = 0.
          zrt(n) = -4.*c2
          DO k = 2,n - 1,2
          zrt(k) = -4.*c2*sin(k*del)**2
          zrt(k+1) = zrt(k)
          enddo
          endif

          if (np.eq.2) then
          DO k = 1,n
          zrt(k) = -4.*c2*sin(k*del)**2
          enddo
          endif

          if (np.eq.3 .or. np.eq.5) then
          DO k = 1,n
          zrt(k) = -4.*c2*sin((k-.5)*del)**2
          enddo
          endif

          if (np.eq.4) then
          DO k = 1,n
          zrt(k) = -4.*c2*sin((k-1)*del)**2
          enddo
          endif

      END IF

      IF (l.GT.1) THEN

C                               FACTOR M*N TRIDIAGONAL SYSTEMS.
C                               FIRST, DO THE POSSIBLY SINGULAR
C                               CASE CORRESPONDING TO J = K = 1.

          fctrd(1,1,1) = 1./ (b(1)+d(1)*cfy(1)+zrt(1))
          DO 310 i = 2,l - 1
              fctrd(i,1,1) = 1./ (b(i)+d(i)*cfy(1)+zrt(1)-
     *                       a(i)*c(i-1)*fctrd(i-1,1,1))
  310     CONTINUE

C                               IF TRIDIAGONAL SYSTEM
C                               (...,A(I),B(I),C(I),...) IS SINGULAR
C                               THEN FCTRD(1,1,l) IS 1./0.  IF
C                               DENOMINATOR IS WITHIN ROUND-OFF OF 0,
C                               SET FCTRD(1,1,l) ARBITRARILY

          den = b(l) + d(l)*cfy(1) + zrt(1) - a(l)*c(l-1)*fctrd(l-1,1,1)
          bmax = abs(b(1))
          DO 320 i = 2,l
              bmax = max(bmax,abs(b(i)))
  320     CONTINUE

          IF (abs(den/bmax).LE.10.*eps) den = bmax
          fctrd(l,1,1) = 1./den

C                               FACTOR CASES J=1, K=2,...,N.

          DO 330 k = 2,n
              fctrd(1,k,1) = 1./ (b(1)+d(1)*cfy(1)+zrt(k))
  330     CONTINUE
          DO 350 i = 2,l
              DO 340 k = 2,n
                  fctrd(i,k,1) = 1./ (b(i)+d(i)*cfy(1)+zrt(k)-
     *                           a(i)*c(i-1)*fctrd(i-1,k,1))
  340         CONTINUE
  350     CONTINUE

C                               FACTOR CASES K=1, J=2,...,M.

          DO 360 j = 2,m
              fctrd(1,1,j) = 1./ (b(1)+d(1)*cfy(j)+zrt(1))
  360     CONTINUE
          DO 380 i = 2,l
              DO 370 j = 2,m
                  fctrd(i,1,j) = 1./ (b(i)+d(i)*cfy(j)+zrt(1)-
     *                           a(i)*c(i-1)*fctrd(i-1,1,j))
  370         CONTINUE
  380     CONTINUE

C                               FACTOR REMAINING CASES.

          DO 400 k = 2,n
              DO 390 j = 2,m
                  fctrd(1,k,j) = 1./ (b(1)+d(1)*cfy(j)+zrt(k))
  390         CONTINUE
  400     CONTINUE
          DO 430 i = 2,l
              DO 420 k = 2,n
                  DO 410 j = 2,m
                      fctrd(i,k,j) = 1./ (b(i)+d(i)*cfy(j)+zrt(k)-
     *                               a(i)*c(i-1)*fctrd(i-1,k,j))
  410             CONTINUE
  420         CONTINUE
  430     CONTINUE
      END IF

C                               INITIALIZE FFT TRANSFORMS AND
C                               PRE-PROCESSING COEFFICIENTS IN J

      IF (m.NE.1) THEN

          SELECT CASE(igrid)
          CASE(1)   ; GOTO 440
          CASE(2)   ; GOTO 450
          END SELECT

  440     CONTINUE
          SELECT CASE(mp)
          CASE(1)   ; GOTO 460
          CASE(2)   ; GOTO 470
          CASE(3)   ; GOTO 480
          CASE(4)   ; GOTO 500
          CASE(5)   ; GOTO 510
          END SELECT

  450     CONTINUE
          SELECT CASE(mp)
          CASE(1)   ; GOTO 460
          CASE(2)   ; GOTO 480
          CASE(3)   ; GOTO 490
          CASE(4)   ; GOTO 510
          CASE(5)   ; GOTO 520
          END SELECT

  460     CONTINUE
          CALL vsrfti(m,wsavey)

          GO TO 530

  470     CONTINUE
          CALL vsinti(m,cfy,wsavey)

          GO TO 530

  480     CONTINUE
ceb          CALL vssini(m,cfy(1),cfy(m+1),wsavey)
          CALL VSCOSI(m,cfy(1),cfy(m+1),wsavey)

          GO TO 530

  490     CONTINUE
ceb          CALL vssnqi(m,cfy(1),cfy(m+1),cfy(2*m+1),cfy(3*m+1),wsavey)
          CALL VSCSQI(m,cfy(1),cfy(m+1),cfy(2*m+1),cfy(3*m+1),wsavey)

          GO TO 530

  500     CONTINUE
          CALL vcosti(m,cfy,wsavey)

          GO TO 530

  510     CONTINUE
          CALL vscosi(m,cfy(1),cfy(m+1),wsavey)

          GO TO 530

  520     CONTINUE
          CALL vscsqi(m,cfy(1),cfy(m+1),cfy(2*m+1),cfy(3*m+1),wsavey)
  530     CONTINUE
      END IF

C                               INITIALIZE FFT TRANSFORMS AND
C                               PRE-PROCESSING COEFFICIENTS IN K

      IF (n.NE.1) THEN

          SELECT CASE(igrid)
          CASE(1)   ; GOTO 540
          CASE(2)   ; GOTO 550
          END SELECT

  540     CONTINUE
          SELECT CASE(np)
          CASE(1)   ; GOTO 560
          CASE(2)   ; GOTO 570
          CASE(3)   ; GOTO 580
          CASE(4)   ; GOTO 600
          CASE(5)   ; GOTO 610
          END SELECT

  550     CONTINUE
          SELECT CASE(np)
          CASE(1)   ; GOTO 560
          CASE(2)   ; GOTO 580
          CASE(3)   ; GOTO 590
          CASE(4)   ; GOTO 610
          CASE(5)   ; GOTO 620
          END SELECT

  560     CONTINUE
          CALL vsrfti(n,wsavez)

          GO TO 630

  570     CONTINUE
          CALL vsinti(n,cfz,wsavez)

          GO TO 630

  580     CONTINUE
ceb          CALL vssini(n,cfz(1),cfz(n+1),wsavez)
          CALL VSCOSI(n,cfz(1),cfz(n+1),wsavez)

          GO TO 630

  590     CONTINUE
ceb          CALL vssnqi(n,cfz(1),cfz(n+1),cfz(2*n+1),cfz(3*n+1),wsavez)
          CALL VSCSQI(n,cfz(1),cfz(n+1),cfz(2*n+1),cfz(3*n+1),wsavez)

          GO TO 630

  600     CONTINUE
          CALL vcosti(n,cfz,wsavez)

          GO TO 630

  610     CONTINUE
          CALL vscosi(n,cfz(1),cfz(n+1),wsavez)

          GO TO 630

  620     CONTINUE
          CALL vscsqi(n,cfz(1),cfz(n+1),cfz(2*n+1),cfz(3*n+1),wsavez)
  630     CONTINUE
      END IF

      RETURN
      END SUBROUTINE fsh01s
C
C
      SUBROUTINE fsh02s(ldimf,mdimf,f,save,w)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   f(ldimf,mdimf,*),save(-3:*),w(*)

C                               RETRIEVE CONSTANTS FROM SAVE ARRAY

      l = save(2)
      lp = save(3)
      m = save(4)
      mp = save(5)
      n = save(6)
      np = save(7)
      igrid = save(11)

C                               ALLOCATION OF SAVE ARRAY

      ia = 12
      ic = ia + 2*l
      icfy = ic + l

      IF (igrid.EQ.1) THEN
          icfz = icfy + 2*m
          ifctrd = icfz + 2*n
      ELSE
          icfz = icfy + 4*m
          ifctrd = icfz + 4*n
      END IF

      ldimfc=l
      if (ldimf.gt.l .and. mod(l,2).eq.0) ldimfc=l+1
      iwsy = ifctrd + ldimfc*m*n
      iwsz = iwsy + m + 16
      ldimft=l

      IF (ldimf.EQ.l .AND. mdimf.EQ.m) THEN

C                               NO HOLES IN DATA ARRAY, SO CALL SOLVER

          CALL fsh03s(igrid,l,lp,m,mp,n,np,ldimfc,ldimft,f,save(icfy),
     *                save(icfz),
     *                w,save(ia),save(ic),save(ifctrd),save(iwsy),
     *                save(iwsz))
      ELSE
      if (ldimf.gt.l .and. mod(l,2).eq.0) ldimft=l+1

C                               PACK DATA ARRAY, CALL SOLVER,
C                               AND THEN UNPACK SOLUTION ARRAY

          CALL fsh04s(l,m,n,ldimf,mdimf,ldimft,f,w)
          CALL fsh03s(igrid,l,lp,m,mp,n,np,ldimfc,ldimft,w,save(icfy),
     *                save(icfz),f,
     *                save(ia),save(ic),save(ifctrd),save(iwsy),
     *                save(iwsz))
          CALL fsh05s(l,m,n,ldimf,mdimf,ldimft,f,w)
      END IF

      RETURN
      END SUBROUTINE fsh02s
C
C
      SUBROUTINE fsh03s(igrid,l,lp,m,mp,n,np,ldimfc,ldimft,f,cfy,cfz,
     *                  ft,a,c,fctrd,wsavey,wsavez)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   a(l),c(l),f(ldimft,m,n),cfy(4*m),cfz(4*n),
     *          ft(ldimft,m,n),
     *          fctrd(ldimfc,m,n),wsavey(m+16),wsavez(n+16)
      LOGICAL datary,datasw

C                               ZERO OUT BOTTOM PLANE OF ARRAY FT

      DO 99 k=1,n
      DO 98 j=1,m
      ft(ldimft,j,k)=0.
   98 end do
   99 end do

      nocopy=.TRUE.
      datary=.TRUE.
      scale=1.
      ifwrd = 1
  100 CONTINUE

      IF (n.NE.1) THEN
      tpose=.false.
      if (ifwrd.eq.2) tpose=.true.

C                               TRANSFORM IN Z

        if (datary) then
            call fsh26s(igrid,ifwrd,np,l,n,m,ldimft,f,ft,cfz,wsavez)
            datary=outary
        else
            call fsh26s(igrid,ifwrd,np,l,n,m,ldimft,ft,f,cfz,wsavez)
            datary=.not.outary
        endif

      END IF
      SELECT CASE(ifwrd)
      CASE(1)   ; GOTO 490
      CASE(2)   ; GOTO 510
      END SELECT

  490 continue
      IF (m.NE.1) THEN
      tpose=.true.
      if (ifwrd.eq.2) tpose=.false.

C                               TRANSFORM Y

C
        if (datary) then
            call fsh26s(igrid,ifwrd,mp,l,m,n,ldimft,f,ft,cfy,wsavey)
            datary=outary
        else
            call fsh26s(igrid,ifwrd,mp,l,m,n,ldimft,ft,f,cfy,wsavey)
            datary=.not.outary
        endif 

      END IF

      SELECT CASE(ifwrd)
      CASE(1)   ; GOTO 285
      CASE(2)   ; GOTO 100
      END SELECT
  285 continue

          IF (l.GT.1) THEN

C                               SOLVE TRIDIAGONAL SYSTEMS IN X THAT WERE
C                               PREVIOUSLY FACTORED IN FSH01S

C                               CALL VECTORIZED TRIDIAGONAL SOLVER

              datasw=.false.
              if (np.eq.1) datasw=.not.datasw
              if (mp.eq.1) datasw=.not.datasw
              if (datary) then
                 if (datasw) then
              CALL fsh06s(l,lp,m*n,ldimfc,ldimft,scale,a,c,f,ft,fctrd)
                 datary=.false.
                 else
              CALL fsh06s(l,lp,m*n,ldimfc,ldimft,scale,a,c,f,f,fctrd)
                 endif
              else
                 if (datasw) then
              CALL fsh06s(l,lp,m*n,ldimfc,ldimft,scale,a,c,ft,f,fctrd)
                 datary=.true.
                 else
              CALL fsh06s(l,lp,m*n,ldimfc,ldimft,scale,a,c,ft,ft,fctrd)
                 endif
              endif

          END IF

          ifwrd = 2

          GO TO 490

  510 continue

      if (.not.datary) then
         do 513 k=1,n
         do 512 j=1,m
         do 511 i=1,l
            f(i,j,k)=ft(i,j,k)
  511       end do
  512       end do
  513       end do
      endif
      RETURN
      END SUBROUTINE fsh03s
C
C
      SUBROUTINE fsh04s(l,m,n,ldimf,mdimf,ldimg,f,g)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   f(ldimf,mdimf,n),g(ldimg,m,n)

C                              THIS SUBROUTINE PACKS THE SUB-ARRAY
C                              F(I,J,K), I=1,...,L, J=1,...,M, K=1,...,N
C                              INTO THE ARRAY G.

      DO 120 k = 1,n
          DO 110 j = 1,m
              DO 100 i = 1,l
                  g(i,j,k) = f(i,j,k)
  100         CONTINUE
  110     CONTINUE
  120 CONTINUE
      if (ldimg.gt.l) then
         do 201 k=1,n
         do 200 j=1,m
            g(ldimg,j,k)=0.
  200    continue
  201    continue
      endif
      RETURN
      END SUBROUTINE fsh04s
C
C
      SUBROUTINE fsh05s(l,m,n,ldimf,mdimf,ldimg,f,g)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   f(ldimf,mdimf,n),g(ldimg,m,n)

C                               THIS SUBROUTINE EXPANDS THE ARRAY G OF
C                               DIMENSION L X M X N INTO THE ARRAY F OF
C                               DIMENSION LDIMF X MDIMF X N.

      DO 120 k = 1,n
          DO 110 j = 1,m
              DO 100 i = 1,l
                  f(i,j,k) = g(i,j,k)
  100         CONTINUE
  110     CONTINUE
  120 CONTINUE
      RETURN
      END SUBROUTINE fsh05s
C
C
      SUBROUTINE fsh06s(l,lp,m,ldimfc,ldimft,scale,a,c,f,ft,fctrd)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   a(l),c(l),f(ldimft,m),ft(ldimft,m),fctrd(ldimfc,m)

C                               THIS SUBROUTINE SOLVES M TRIDIAGONAL
C                               SYSTEMS OF ORDER L THAT HAVE
C                               FACTORIZATION STORED IN ARRAY FCTRD AND
C                               RIGHT SIDES IN FT

      scale=scale*scale

      IF (lp.eq.1) THEN 
         lh=(l+1)/2
         lq=(lh-1)/2

         DO 610 i=1,lq
            do 600 j=1,m
               f(i,j)=f(i,j)-f(2*lh-i,j)
               f(lh-i,j)=f(lh-i,j)-f(lh+i,j)
               f(lh+i,j)=f(lh-i,j)+2.*f(lh+i,j)
               f(2*lh-i,j)=f(i,j)+2.*f(2*lh-i,j)
  600       continue
            call sswap(m,f(i,1),ldimft,f(lh-i,1),ldimft)
  610    continue
         DO 620 j=1,m
            f(lh,j)=2.*f(lh,j)
  620       end do
         IF (mod(l,2).eq.0)  THEN
            DO 630 j=1,m
               f(l,j)=2.*f(l,j)
  630          end do
         ENDIF
         if (mod(lh,2).eq.0) then
            do 640 j=1,m
               f(  lh/2,j)=f(lh/2,j)-f(3*lh/2,j)
               f(3*lh/2,j)=f(lh/2,j)+2.*f(3*lh/2,j)
  640       continue
         endif
         scale=0.5*scale
      ENDIF
C                               FORWARD SUBSTITUTION

         DO 120 j = 1,m
             ft(1,j) = scale*f(1,j)*fctrd(1,j)
         DO 110 i = 2,l
                 ft(i,j) = (scale*f(i,j)-a(i)*ft(i-1,j))*fctrd(i,j)
  110        CONTINUE
  120    CONTINUE

C                               BACKWARD SUBSTITUTION

         DO 140 j = 1,m
         DO 130 i = l - 1,1,-1
                 ft(i,j) = ft(i,j) - c(i)*fctrd(i,j)*ft(i+1,j)
  130        CONTINUE
  140    CONTINUE

      IF (lp.eq.1) THEN

         do 710 i=1,lq
            do 700 j=1,m
               ft(i,j)=ft(lh+i,j)+ft(i,j)
               ft(lh-i,j)=ft(2*lh-i,j)+ft(lh-i,j)
               ft(2*lh-i,j)=2.*ft(2*lh-i,j)-ft(lh-i,j)
               ft(  lh+i,j)=2.*ft(lh+i,j)-ft(i,j)
  700       CONTINUE
            call sswap(m,ft(i,1),ldimft,ft(lh-i,1),ldimft)
  710    CONTINUE
         if (mod(lh,2).eq.0) then
             do 720 j=1,m
                ft(  lh/2,j)=ft(3*lh/2,j)+ft(lh/2,j)
                ft(3*lh/2,j)=2.*ft(3*lh/2,j)-ft(lh/2,j)
  720        continue
         endif

      ENDIF 
         
      RETURN
      END SUBROUTINE fsh06s
C
C
      FUNCTION FSH19S()

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+


C                               THIS FUNCTION FURNISHES PI TO MACHINE
C                               PREC

C                               PI=3.14159265358979323846264338327950288

      fsh19s = 3.141592653589793_EB
      RETURN
      END FUNCTION FSH19S
C
C
      FUNCTION fsh20s()

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C                               THIS FUNCTION FURNISHES THE MACHINE
C                               PREC EPS--THE SMALLEST POSITIVE
C                               MACHINE NUMBER SATISFYING

C                               FL(1.+EPS) > 1.

c     fsh20s = 2.4E-7_EB
      fsh20s = 4.E-15_EB

      RETURN
      END FUNCTION fsh20s
C
C
      SUBROUTINE fsh26s(igrid,ifwrd,mp,l,m,n,ldimft,f,ft,cfy,wsavey)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   f(ldimft,m,n), ft(ldimft,m,n), cfy(4*m), wsavey(m+16)

          SELECT CASE(igrid)
          CASE(1)   ; GOTO 110
          CASE(2)   ; GOTO 120
          END SELECT

  110     CONTINUE
          SELECT CASE(mp)
          CASE(1)   ; GOTO 130
          CASE(2)   ; GOTO 160
          CASE(3)   ; GOTO 170
          CASE(4)   ; GOTO 220
          CASE(5)   ; GOTO 230
          END SELECT
  120     CONTINUE
          SELECT CASE(mp)
          CASE(1)   ; GOTO 130
          CASE(2)   ; GOTO 180
          CASE(3)   ; GOTO 210
          CASE(4)   ; GOTO 240
          CASE(5)   ; GOTO 270
          END SELECT
  130     CONTINUE
          SELECT CASE(ifwrd)
          CASE(1)   ; GOTO 140
          CASE(2)   ; GOTO 150
          END SELECT
  140     CONTINUE
          CALL vsrftf(f,l,m,n,ldimft,ft,wsavey)
          GO TO 280
  150     CONTINUE
          CALL vsrftb(f,l,m,n,ldimft,ft,wsavey)
          GO TO 280
  160     CONTINUE
          CALL vsint(f,l,m,n,ldimft,ft,cfy,wsavey)
          GO TO 280
  170     CONTINUE
          SELECT CASE(ifwrd)
          CASE(1)   ; GOTO 200
          CASE(2)   ; GOTO 190
          END SELECT
  180     CONTINUE
          SELECT CASE(ifwrd)
          CASE(1)   ; GOTO 190
          CASE(2)   ; GOTO 200
          END SELECT
  190     CONTINUE
          CALL vssinf(f,l,m,n,ldimft,ft,cfy(1),cfy(m+1),wsavey)
          GO TO 280
  200     CONTINUE
          CALL vssinb(f,l,m,n,ldimft,ft,cfy(1),cfy(m+1),wsavey)
          GO TO 280
  210     CONTINUE
          CALL vssinq(f,l,m,n,ldimft,ft,cfy(1),cfy(m+1),
     *                cfy(2*m+1),cfy(3*m+1),wsavey)
          GO TO 280
  220     CONTINUE
          CALL vcost(f,l,m,n,ldimft,ft,cfy,wsavey)
          GO TO 280
  230     CONTINUE
          SELECT CASE(ifwrd)
          CASE(1)   ; GOTO 260
          CASE(2)   ; GOTO 250
          END SELECT
  240     CONTINUE
          SELECT CASE(ifwrd)
          CASE(1)   ; GOTO 250
          CASE(2)   ; GOTO 260
          END SELECT
  250     CONTINUE
          CALL vscosf(f,l,m,n,ldimft,ft,cfy(1),cfy(m+1),wsavey)
          GO TO 280
  260     CONTINUE
          CALL vscosb(f,l,m,n,ldimft,ft,cfy(1),cfy(m+1),wsavey)
          GO TO 280
  270     CONTINUE
          CALL vscosq(f,l,m,n,ldimft,ft,cfy(1),cfy(m+1),
     *                cfy(2*m+1),cfy(3*m+1),wsavey)
  280     CONTINUE
      
      return
      END SUBROUTINE fsh26s
C
C
      SUBROUTINE VCOST(x,l,m,n,LDIMX,xt,c,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)         X(LDIMX*N,M),xt(lDIMX*n,m),c(m),WSAVE(M+15)
C
C     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
C           = .FALSE. IF TRANSFORMING THIRD INDEX
C
      mM1 = m-1
      mS2 = m/2
      scale=scale*sqrt(0.5_EB)
      IF (TPOSE) THEN
          call vcost1(l,m,n,LDIMX,ms2,c,xt,x)
      ELSE
         do 401 i=1,ldimx*n
         xt(i,m)=  x(i,1)-x(i,m)
         Xt(I,1) = X(I,1)+X(I,m)
  401    end do
         DO 104 j=2,mS2
         jC = m+1-j
         do 402 i=1,ldimx*n
         xt(i,m) = xt(i,m)+c(jC)*(x(i,j)-x(i,jc))
         Xt(I,j) =  (X(I,j)+X(I,jC))-(c(j)*(x(i,j)-x(i,jc)))
         Xt(I,jC) = (X(I,j)+X(I,jC))+(c(j)*(x(i,j)-x(i,jc)))
  402    continue
  104    CONTINUE
         IF (MOD(M,2) .eq. 0) go to 404
         do 403 i=1,ldimx*n
         Xt(I,mS2+1) = 2.*X(I,ms2+1)
  403    end do
  404    continue
      end if
      IF (M.GT.3) THEN
  450    CALL vRFFTF (LDIMX*n,mM1,Xt,LDIMX*n,x,WSAVE)
         if (outary) then
            call vcosta(l,m,n*LDIMX,mm1,XT(1,M),x,xt)
         else
            call vcosta(l,m,n*LDIMX,mm1,XT(1,M),xt,x)
         endif
      ELSEIF (M.EQ.2) THEN
         OUTARY=.FALSE.
      ELSE
         DO 302 I=1,Ldimx*n
            X(I,2)=XT(I,3)
            X(I,1)=XT(I,1)+XT(I,2)
            X(I,3)=XT(I,1)-XT(I,2)
  302       end do
         OUTARY=.TRUE.
         SCALE=SCALE*SQRT(0.5_EB)
      ENDIF
      return
      end SUBROUTINE VCOST
C
C
      SUBROUTINE VCOST1(L,M,N,LDIMX,MS2,c,XT,X)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   X(LDIMX,M,N),XT(LDIMX,N,M),c(*)
      do 405 k=1,n
      do 401 i=1,l
      xt(i,k,m)=x(i,1,K)-x(i,M,K)
      Xt(I,k,1) = X(I,1,K)+X(I,M,K)
  401 end do
  405 end do
      DO 104 j=2,mS2
      jC = m+1-j
      do 412 k=1,n
      do 402 i=1,l
      xt(i,k,m) = xt(i,k,m)+c(jC)*(x(i,j,k)-x(i,jc,k))
      Xt(I,k,j) =  (X(I,j,k)+X(I,jC,k))-(c(j)*(x(i,j,k)-x(i,jc,k)))
      Xt(I,k,jC) = (X(I,j,k)+X(I,jC,k))+(c(j)*(x(i,j,k)-x(i,jc,k)))
  402 continue
  412 continue
  104 CONTINUE
      IF (MOD(M,2) .eq. 0) go to 404
      do 407 k=1,n
      do 403 i=1,l
      Xt(I,k,mS2+1) = 2.*X(I,MS2+1,k)
  403 end do
  407 end do
  404 continue
      RETURN
      END SUBROUTINE VCOST1
C
C
      SUBROUTINE VCOSTA(l,m,LDIMX,mm1,PL,x,xt)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   x(LDIMX,m),xt(LDIMX,m),PL(LDIMX)
      do 405 i=1,ldimx
      x(i,1) = xt(i,1)
      X(I,2) = PL(i)
  405 end do
      DO 105 j=4,m,2
      do 406 i=1,ldimx
      X(I,j) = X(I,j-2)-Xt(I,j-1)
      X(I,j-1) = xt(i,j-2)
  406 end do
  105 CONTINUE
      IF (MOD(M,2) .eq. 0) go to 409
      do 407 i=1,ldimx
      X(I,m) = xt(i,mm1)
  407 end do
  409 continue
      RETURN
      END SUBROUTINE VCOSTA
C
C
      SUBROUTINE vCOSTI(N,c,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

C
C     *******   NOTE   ******
C
C     SPECIAL DOCUMENTATION FOR USE OF INTERNAL PARAMETER TPOSE.  IF
C     TPOSE = .TRUE., THEN ALL PREPROCESSING ROUTINES TRANSPOSE SECOND
C                     AND THIRD INDICES IN PREPARATION FOR TRANSFORM OF
C                     ORIGINAL SECOND INDEX.  IF
C     TPOSE = .FALSE., THEN ALL PREPROCESSING ROUTINES PREPARE ONLY FOR
C                      TRANSFORM OF THIRD INDEX.  TO USE IN THIS CASE,
C                      CALL ROUTINE WITH M AND N INTERCHANGED.
C
      REAL(EB)    c(n), WSAVE(n+15)

C                               INITIALIZE NOCOPY AND TPOSE TO DEFAULT
C                               VALUES

      nocopy = .FALSE.
      tpose = .TRUE.

      PI = fsh19s()
      IF (N .LE. 3) RETURN
      NP1 = N+1
      NS2 = N/2
      DT = PI/(n-1)
C
      DO 101 K=2,NS2
      KC = NP1-K
      c(K) = 2.*SIN((k-1)*dt)
      c(KC) = 2.*COS((K-1)*DT)
  101 CONTINUE
      CALL vRFFTI (N-1,WSAVE)
      RETURN
      END SUBROUTINE vCOSTI
C
C
      SUBROUTINE VRFFTF (M,N,R,MDIMR,rt,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)         R(MDIMR,N)  ,RT(M,N)    ,WSAVE(N+15)
      IF (N .EQ. 1) RETURN
      CALL VRFTF1 (M,N,R,MDIMR,rt,WSAVE(1),WSAVE(N+1))
      RETURN
      END SUBROUTINE VRFFTF
C
C
      SUBROUTINE VRFFTI (N,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)         WSAVE(N+15)

C                               INITIALIZE NOCOPY AND TPOSE TO DEFAULT
C                               VALUES

      nocopy = .FALSE.
      tpose = .TRUE.

      IF (N .LE. 1) RETURN
      CALL VRFTI1 (N,WSAVE(1),WSAVE(N+1))
      RETURN
      END SUBROUTINE VRFFTI
C
C
      SUBROUTINE VRFTF1 (M,N,C,MDIMC,ch,WA,FAC)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)         C(MDIMC,N) ,Ch(M,N)  ,WA(N)   ,FAC(15)
C
      NF = FAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = FAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL VRADF4 (M,IDO,L1,C,MDIMC,ch,m,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL VRADF4 (M,IDO,L1,CH,m,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL VRADF2 (M,IDO,L1,C,MDIMC,ch,m,WA(IW))
         GO TO 110
  103    CALL VRADF2 (M,IDO,L1,CH,m,C,MDIMC,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL VRADF3 (M,IDO,L1,C,MDIMC,ch,m,WA(IW),WA(IX2))
         GO TO 110
  105    CALL VRADF3 (M,IDO,L1,CH,M,c,mDIMC,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
      CALL VRADF5(M,IDO,L1,C,MDIMC,ch,m,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107 CALL VRADF5 (M,IDO,L1,CH,m,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL VRADFG (M,IDO,IP,L1,IDL1,C,C,C,mdimc,CH,CH,M,WA(IW))
         NA = 1
         GO TO 110
  109    CALL VRADFG (M,IDO,IP,L1,IDL1,CH,CH,CH,m,C,C,MDIMC,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      outary=.TRUE.
      if (nocopy) then
         scale=scale*sqrt(1.0_EB/real(n,EB))
         if (na.eq.0) then
            outary=.FALSE.
         endif
      else
         SCALE=SQRT(1.0_EB/REAL(N,EB))
         IF (NA .EQ. 1) GO TO 113
         DO 115 J=1,N
         DO 112 I=1,M
            C(I,J) = SCALE*CH(I,J)
  112    CONTINUE
  115    CONTINUE
         RETURN
  113    DO 116 J=1,N
         DO 114 I=1,M
            C(I,J)=SCALE*C(I,J)
  114    CONTINUE
  116    CONTINUE
      endif
      RETURN
      END SUBROUTINE VRFTF1
C
C
      SUBROUTINE VRFTI1 (N,WA,FAC)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)         WA(N)      ,FAC(15)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4.LT.0) GOTO 102
      IF (J-4.EQ.0) GOTO 102
      IF (J-4.GT.0) GOTO 103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR.LT.0) GOTO 101
      IF (NR.EQ.0) GOTO 105
      IF (NR.GT.0) GOTO 101
  105 NF = NF+1
      FAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         FAC(IB+2) = FAC(IB+1)
  106 CONTINUE
      FAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      FAC(1) = N
      FAC(2) = NF
      TPI = 2.*fsh19s()
      ARGH = TPI/REAL(N,EB)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = FAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = REAL(LD,EB)*ARGH
            FI = 0.
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END SUBROUTINE VRFTI1
C
C
      SUBROUTINE VSCOSB(F,L,M,N,LDIMF,FT,C1,C2,WORK)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF*N,M),FT(LDIMF*N,M),C1(M),C2(M),WORK(M+15)
C
C     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
C           = .FALSE. IF TRANSFORMING THIRD INDEX
C
C
C     PREPROCESSING
C
         IF (TPOSE) THEN
            CALL VSCSB1(L,M,N,LDIMF,F,FT,C1,C2)
         ELSE
            DO 100 I=1,Ldimf*n
               FT(I,1)=.5*F(I,1)
  100          end do
            DO 201 J=2,M
            DO 200 I=1,Ldimf*n
               FT(I,J) =(C1(J)*F(I,J)+C2(J)*F(I,M-J+2))
  200          end do
  201          end do
         ENDIF
C
C     REAL,PERIODIC ANALYSIS
C
      CALL VRFFTF(LDIMF*N,M,FT,LDIMF*N,F,WORK)
C
C     POSTPROCESSING
C
      SCALE=SQRT(2.0_EB)*SCALE
      IF (OUTARY) THEN
         CALL VSCSBA(L,M,N*LDIMF,FT,F)
      ELSE
         CALL VSCSBA(L,M,N*LDIMF,F,FT)
      ENDIF
      RETURN
      END SUBROUTINE VSCOSB
C
C
      SUBROUTINE VSCOSF(F,L,M,N,LDIMF,FT,C1,C2,WORK)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF*N,M),FT(LDIMF*N,M),C1(M),C2(M),WORK(M+15)
C
C     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
C           = .FALSE. IF TRANSFORMING THIRD INDEX
C
C
C     PREPROCESSING
C
      SCALE=SQRT(2.0_EB)*SCALE
      IF (TPOSE) THEN
         CALL VSCSF1(L,M,N,LDIMF,F,FT)
      ELSE
         DO 100 I=1,Ldimf*n
            FT(I,1)=F(I,1)
  100       end do
         CONTINUE
         IF (MOD(M,2).EQ.0) THEN
            DO 130 I=1,Ldimf*n
               FT(I,M)=F(I,M)
  130          end do
         ENDIF
         DO 201 J=2,M-1,2
         DO 200 I=1,Ldimf*n
            FT(I,J)   =  0.5*(F(I,J)+F(I,J+1))
            FT(I,J+1) = -0.5*(F(I,J)-F(I,J+1))
  200       end do
  201       end do
      ENDIF
C
C     REAL,PERIODIC SYNTHESIS
C
      CALL VRFFTB(LDIMF*N,M,FT,LDIMF*N,F,WORK)
C
C     POSTPROCESSING
C
      IF (OUTARY) THEN
         CALL VSCSFA(L,M,N*LDIMF,F,FT,C1,C2)
      ELSE
         CALL VSCSFA(L,M,N*LDIMF,FT,F,C1,C2)
      ENDIF
      RETURN
      END SUBROUTINE VSCOSF
C
C
      SUBROUTINE VSCOSI(N,C1,C2,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

c      ENTRY VSSINI(N,C1,C2,WSAVE)
      REAL(EB)   C1(N),C2(N),WSAVE(N+15)

      PI=fsh19s()
      DX=PI/(2*N)
C
C     GENERATE A(I)+-B(I)
C
      DO 100 I=1,N
         C=COS((I-1)*DX)
         S=SIN((I-1)*DX)
         C1(I)=.5*(S+C)
         C2(I)=.5*(S-C)
  100    end do
C
C     INITIALIZE VRFFTPK ROUTINE
C
      CALL VRFFTI(N,WSAVE)
      RETURN
      END SUBROUTINE VSCOSI
C
C
      SUBROUTINE VSCOSQ(F,L,M,N,LDIMF,FT,C1,C2,C3,C4,WORK)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF*N,M),FT(LDIMF*N,M),
     *          C1(M),C2(M),C3(M),C4(M),WORK(M+15)
C
C     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
C           = .FALSE. IF TRANSFORMING THIRD INDEX
C
C
C     PREPROCESSING
C
      IF (TPOSE) THEN
         CALL VSCSQ1(L,M,N,LDIMF,F,FT,C1,C2)
      ELSE
         DO 100 I=1,Ldimf*n
            FT(I,1)=F(I,1)
  100       end do
         IF (MOD(M,2).EQ.0) THEN
            DO 110 I=1,Ldimf*n
               FT(I,M)=-F(I,M)
  110          end do
         ENDIF
         DO 201 J=2,M-1,2
            JBY2=J/2
            DO 200 I=1,Ldimf*n
               FT(I,J)   =  F(I,J+1)*C1(JBY2)+F(I,J)*C2(JBY2)
               FT(I,J+1) = -F(I,J+1)*C2(JBY2)+F(I,J)*C1(JBY2)
  200          end do
  201          end do
      ENDIF
C
C     REAL,PERIODIC SYNTHESIS
C
      CALL VRFFTB(LDIMF*N,M,FT,LDIMF*N,F,WORK)
C
C     POSTPROCESSING
C
      IF (OUTARY) THEN
         DO 301 J=1,M
         DO 300 I=1,Ldimf*n
            F(I,J)=C4(J)*FT(I,J)+C3(J)*FT(I,M+1-J)
  300       end do
  301       end do
      ELSE
         DO 311 J=1,M
         DO 310 I=1,Ldimf*n
            FT(I,J)=C4(J)*F(I,J)+C3(J)*F(I,M+1-J)
  310       end do
  311       end do
      ENDIF
      RETURN
      END SUBROUTINE VSCOSQ
C
C
      SUBROUTINE VSCSB1(L,M,N,LDIMF,F,FT,C1,C2)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M,N), FT(LDIMF,N,M), C1(M), C2(M)
      DO 101 K=1,N
      DO 100 I=1,L
         FT(I,K,1)=.5*F(I,1,K)
  100    end do
  101    end do
      DO 202 K=1,N
      DO 201 J=2,M
      DO 200 I=1,L
         FT(I,K,J) =(C1(J)*F(I,J,K)+C2(J)*F(I,M-J+2,K))
 200     end do
 201     end do
 202     end do
      RETURN
      END SUBROUTINE VSCSB1
C
C
      SUBROUTINE VSCSBA(L,M,LDIMF,FT,F)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M), FT(LDIMF,M)
         DO 300 I=1,Ldimf
            F(I,1) = FT(I,1)
  300       end do
      IF (MOD(M,2).EQ.0) THEN
         DO 310 I=1,Ldimf
            F(I,M) = FT(I,M)
 310        end do
      ENDIF
      DO 401 J=2,M-1,2
         DO 400 I=1,Ldimf
            F(I,J)=FT(I,J)-FT(I,J+1)
            F(I,J+1)=FT(I,J)+FT(I,J+1)
  400       end do
  401       end do
      RETURN
      END SUBROUTINE VSCSBA
C
C
      SUBROUTINE VSCSF1(L,M,N,LDIMF,F,FT)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M,N), FT(LDIMF,N,M)
      DO 120 K=1,N
      DO 100 I=1,L
         FT(I,K,1)=F(I,1,K)
  100  end do
  120 CONTINUE
      IF (MOD(M,2).EQ.0) THEN
         DO 131 K=1,N
         DO 130 I=1,L
            FT(I,K,M)=F(I,M,K)
  130       end do
  131       end do
      ENDIF
      DO 202 K=1,N
      DO 201 J=2,M-1,2
      DO 200 I=1,L
         FT(I,K,J)   =  0.5*(F(I,J,K)+F(I,J+1,K))
         FT(I,K,J+1) = -0.5*(F(I,J,K)-F(I,J+1,K))
  200    end do
  201    end do
  202    end do
      RETURN
      END SUBROUTINE VSCSF1
C
C
      SUBROUTINE VSCSFA(L,M,LDIMF,F,FT,C1,C2)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M),FT(LDIMF,M),C1(M),C2(M)
      DO 301 J=2,M
      DO 300 I=1,Ldimf
         F(I,J)=C1(J)*FT(I,J)-C2(J)*FT(I,M+2-J)
  300    end do
  301    end do
      DO 400 I=1,Ldimf
         F(I,1)=FT(I,1)
  400    end do
      RETURN
      END SUBROUTINE VSCSFA
C
C
      SUBROUTINE VSCSQ1(L,M,N,LDIMF,F,FT,C1,C2)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M,N), FT(LDIMF,N,M), C1(M), C2(M)
      DO 101 K=1,N
      DO 100 I=1,L
         FT(I,K,1)=F(I,1,K)
  100    end do
  101    end do
      IF (MOD(M,2).EQ.0) THEN
         DO 111 K=1,N
         DO 110 I=1,L
            FT(I,K,M)=-F(I,M,K)
  110       end do
  111       end do
      ENDIF
      DO 202 J=2,M-1,2
         JBY2=J/2
         DO 201 K=1,N
         DO 200 I=1,L
            FT(I,K,J)   =  F(I,J+1,K)*C1(JBY2)+F(I,J,K)*C2(JBY2)
            FT(I,K,J+1) = -F(I,J+1,K)*C2(JBY2)+F(I,J,K)*C1(JBY2)
  200       end do
  201       end do
  202       end do
      RETURN
      END SUBROUTINE VSCSQ1
C
C
      SUBROUTINE VSCSQI(N,C1,C2,C3,C4,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

c      ENTRY VSSNQI(N,C1,C2,C3,C4,WSAVE)
      REAL(EB)   C1(N),C2(N),C3(N),C4(N),WSAVE(N+15)

      PI=fsh19s()
      DX=PI/N
      SCALE=SQRT(.5_EB)
C
C     GENERATE A(I)+-B(I)
C
      DO 100 I=1,(N-1)/2
         C=COS(I*DX)
         S=SIN(I*DX)
         C1(I)=.5*(S+C)
         C2(I)=.5*(C-S)
  100    end do
C
      DX=PI/(2*N)
      DO 200 I=1,N
         C=COS((I-.5)*DX)
         S=SIN((I-.5)*DX)
         C3(I)=SCALE*(C+S)
         C4(I)=SCALE*(C-S)
  200    end do
C
C     INITIALIZE VRFFTPK ROUTINE
C
      CALL VRFFTI(N,WSAVE)
      RETURN
      END SUBROUTINE VSCSQI
C
C
      SUBROUTINE vSINT(x,l,m,N,ldimx,xt,c,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)    X(ldimx*n,m),xt(ldimx*n,m+1),c(m/2),WSAVE(M+16)
C
C     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
C           = .FALSE. IF TRANSFORMING THIRD INDEX
C
      sqrt2i=sqrt(.5_EB)
      MODm = MOD(m,2)
  103 mP1 = m+1
      mS2 = m/2
      do 409 i=1,ldimx*n
      xt(i,1) = 0.
      xt(i,m+1)=0.
  409 end do
C                         Zero out last plane because for sine 
C                         it doesn't get done in fsh02
      if (tpose) then
         call vsint1(l,m,n,ldimx,modm,ms2,c,xt,x)
      else
         DO 105 J=1,MS2
         JC=M+1-J
         DO 104 I=1,Ldimx*n
           xt(i,j+1) = (X(I,J)-X(I,JC))+c(j)*(X(I,J)+X(I,JC))
           xt(i,jC+1) = c(j)*(X(I,J)+X(I,JC))-(X(I,J)-X(I,JC))
  104    CONTINUE
  105    CONTINUE
         IF (MODm .eq. 0) go to 405
         do 504 i=1,ldimx*n
           xt(i,mS2+2) = 4.*X(I,MS2+1)
  504      end do
  405    continue
      endif
      CALL vRFFTF (ldimx*n,mP1,xt,ldimx*n,x,WSAVE)
      scale=scale*sqrt2i
      if (outary) then
         call vsinta(m,ldimx*n,modm,x,xt)
      else
         call vsinta(m,ldimx*n,modm,xt,x)
      endif
      return
      END SUBROUTINE vSINT
C
C
      SUBROUTINE VSINT1(L,M,N,ldimx,MODM,MS2,C,XT,X)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   XT(Ldimx,N,M+1), X(Ldimx,M,N), C(*)
         DO 403 j=1,mS2
         jC = m+1-j
         do 402 k=1,n
         do 401 i=1,l
           xt(i,k,j+1) = (X(I,j,K)-X(I,jC,k))+c(j)*(X(I,j,K)+X(I,jC,k))
           xt(i,k,jC+1) = c(j)*(X(I,j,K)+X(I,jC,k))-(X(I,j,K)-X(I,jC,k))
  401      end do
  402      end do
  403      end do
      IF (MODm .eq. 0) go to 405
         do 406 k=1,n
         do 404 i=1,l
            xt(i,k,mS2+2) = 4.*X(I,MS2+1,k)
  404       end do
  406       end do
  405 continue
      return
      END SUBROUTINE VSINT1
C
C
      SUBROUTINE VSINTA(m,ldimx,modm,x,xt)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   x(ldimx,m), xt(ldimx,m+1)
      do 406 i=1,ldimx
      X(I,1) = 0.5*xt(I,1)
  406 end do
      DO 105 j=2,m-1,2
      do 407 i=1,ldimx
         X(I,j) = -xt(i,j+1)
         X(I,j+1) = X(I,j-1)+xt(i,j)
  407    end do
  105 CONTINUE
      IF (MODm .NE. 0) go to 900
      do 408 i=1,ldimx
      X(I,M) = -xt(i,m+1)
  408 continue
  900 continue
      RETURN
      END SUBROUTINE VSINTA
C
C
      SUBROUTINE vSINTI(N,c,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

C     *******   NOTE   ******
C
C     SPECIAL DOCUMENTATION FOR USE OF INTERNAL PARAMETER TPOSE.  IF
C     TPOSE = .TRUE., THEN ALL PREPROCESSING ROUTINES TRANSPOSE SECOND
C                     AND THIRD INDICES IN PREPARATION FOR TRANSFORM OF
C                     ORIGINAL SECOND INDEX.  IF
C     TPOSE = .FALSE., THEN ALL PREPROCESSING ROUTINES PREPARE ONLY FOR
C                      TRANSFORM OF THIRD INDEX.  TO USE IN THIS CASE,
C                      CALL ROUTINE WITH M AND N INTERCHANGED.
C
      REAL(EB)   c(n/2),wsave (N+16)

C                               INITIALIZE NOCOPY AND TPOSE TO DEFAULT
C                               VALUES

      nocopy = .FALSE.
      tpose = .TRUE.

      IF (N .LE. 1) RETURN
      pi=fsh19s()
      NP1 = N+1
      NS2 = N/2
      DT = PI/REAL(NP1,EB)
      DO 101 K=1,ns2
      c(K) = 2.*SIN(K*DT)
  101 CONTINUE
      CALL vRFFTI (NP1,WSAVE)
      RETURN
      END SUBROUTINE vSINTI
C
C
      SUBROUTINE VSRFTB(F,L,M,N,ldimf,FT,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

C
C     *******   NOTE   ******
C
C     SPECIAL DOCUMENTATION FOR USE OF INTERNAL PARAMETER TPOSE.  IF
C     TPOSE = .TRUE., THEN ALL PREPROCESSING ROUTINES TRANSPOSE SECOND
C                     AND THIRD INDICES IN PREPARATION FOR TRANSFORM OF
C                     ORIGINAL SECOND INDEX.  IF
C     TPOSE = .FALSE., THEN ALL PREPROCESSING ROUTINES PREPARE ONLY FOR
C                      TRANSFORM OF THIRD INDEX.  TO USE IN THIS CASE,
C                      CALL ROUTINE WITH M AND N INTERCHANGED.
C
C
C     VSFFTPK, VERSION 2, JUNE 1988
C
      REAL(EB)   F(Ldimf,M,N),FT(Ldimf,N,M),WSAVE(M+15)
C
C     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
C           = .FALSE. IF TRANSFORMING THIRD INDEX
C
      IF (TPOSE) THEN
C
C        RE-ORDER INPUT
C
         DO 202 K=1,N
         DO 201 J=1,M
         DO 200 I=1,L
            FT(I,K,J)=F(I,J,K)
  200    CONTINUE
  201    CONTINUE
  202    CONTINUE
C
C        REAL, PERIODIC TRANSFORM
C
         CALL VRFFTB(Ldimf*N,M,FT,Ldimf*n,F,WSAVE)
         outary=.not.outary
         if (.not.nocopy) then
            CALL VSRTB1(L,M,N*LDIMF,F,FT)
         endif
      ELSE
         CALL VRFFTB(Ldimf*N,M,F,Ldimf*n,FT,WSAVE)
      ENDIF
      RETURN
      END SUBROUTINE VSRFTB
C
C
      SUBROUTINE VSRFTF(F,L,M,N,ldimf,FT,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(Ldimf,N,M),FT(Ldimf,N,M),WSAVE(M+15)
C
C     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
C           = .FALSE. IF TRANSFORMING THIRD INDEX
C
      IF (TPOSE) THEN
C
C        RE-ORDER INPUT
C
         CALL VSRTF1(L,M,N,LDIMF,F,FT)
C
C        REAL, PERIODIC TRANSFORM
C
         CALL VRFFTF(Ldimf*N,M,FT,Ldimf*n,F,WSAVE)
         outary=.not.outary
         if (.not.nocopy) then
            CALL VSRTB1(L,M,N*LDIMF,F,FT)
         endif
      ELSE
         CALL VRFFTF(Ldimf*N,M,F,Ldimf*n,FT,WSAVE)
      ENDIF
      RETURN
      END SUBROUTINE VSRFTF
C
C
      SUBROUTINE VSRFTI(N,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

C
C     *******   NOTE   ******
C
C     SPECIAL DOCUMENTATION FOR USE OF INTERNAL PARAMETER TPOSE.  IF
C     TPOSE = .TRUE., THEN ALL PREPROCESSING ROUTINES TRANSPOSE SECOND
C                     AND THIRD INDICES IN PREPARATION FOR TRANSFORM OF
C                     ORIGINAL SECOND INDEX.  IF
C     TPOSE = .FALSE., THEN ALL PREPROCESSING ROUTINES PREPARE ONLY FOR
C                      TRANSFORM OF THIRD INDEX.  TO USE IN THIS CASE,
C                      CALL ROUTINE WITH M AND N INTERCHANGED.
C
      REAL(EB)   WSAVE(N+15)
C
C     INITIALIZE VRFFTPK ROUTINE
C
      CALL VRFFTI(N,WSAVE)
      RETURN
      END SUBROUTINE VSRFTI
C
C
      SUBROUTINE VSRTB1(L,M,LDIMF,F,FT)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M),FT(LDIMF,M)
      DO 301 J=1,M
      DO 300 I=1,Ldimf
      F(I,J)=FT(I,J)
  300 CONTINUE
  301 CONTINUE
      RETURN
      END SUBROUTINE VSRTB1
C
C
      SUBROUTINE VSRTF1(L,M,N,LDIMF,F,FT)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M,N),FT(LDIMF,N,M)
      DO 102 K=1,N
      DO 101 J=1,M
      DO 100 I=1,L
      FT(I,K,J)=F(I,J,K)
  100 end do
  101 end do
  102 end do
      RETURN
      END SUBROUTINE VSRTF1
C
C
      SUBROUTINE VSSINB(F,L,M,N,LDIMF,FT,C1,C2,WORK)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF*N,M),FT(LDIMF*N,M),C1(M),C2(M),WORK(M+15)
C
C     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
C           = .FALSE. IF TRANSFORMING THIRD INDEX
C
C
C     PREPROCESSING
C
      IF (TPOSE) THEN
         CALL VSSNB1(L,M,N,LDIMF,F,FT,C1,C2)
      ELSE
         DO 101 J=2,M
         DO 100 I=1,Ldimf*n
            FT(I,J)=C1(J)*F(I,J-1)-C2(J)*F(I,M-J+1)
  100       end do
  101       end do
         DO 200 I=1,Ldimf*n
            FT(I,1) = .5*F(I,M)
  200       end do
      ENDIF
C
C     REAL,PERIODIC ANALYSIS
C
      CALL VRFFTF(LDIMF*N,M,FT,LDIMF*N,F,WORK)
C
C     POSTPROCESSING
C
      SCALE=SQRT(2.0_EB)*SCALE
      IF (OUTARY) THEN
         CALL VSSNBA(L,M,N*LDIMF,F,FT)
      ELSE
         CALL VSSNBA(L,M,N*LDIMF,FT,F)
      ENDIF
      RETURN
      END SUBROUTINE VSSINB
C
C
      SUBROUTINE VSSINF(F,L,M,N,LDIMF,FT,C1,C2,WORK)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF*N,M),FT(LDIMF*N,M),C1(M),C2(M),WORK(M+15)
C
C     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
C           = .FALSE. IF TRANSFORMING THIRD INDEX
C
C
C     PREPROCESSING
C
      SCALE=SQRT(2.0_EB)*SCALE
      IF (TPOSE) THEN
         CALL VSSNF1(L,M,N,LDIMF,F,FT)
      ELSE
         DO 100 I=1,Ldimf*n
            FT(I,1)=F(I,1)
  100       end do
         IF (MOD(M,2).EQ.0) THEN
            DO 110 I=1,Ldimf*n
               FT(I,M)=-F(I,M)
  110          end do
         ENDIF
         DO 201 J=2,M-1,2
         DO 200 I=1,Ldimf*n
            FT(I,J)   = 0.5*(F(I,J+1)-F(I,J))
            FT(I,J+1) =-0.5*(F(I,J+1)+F(I,J))
  200       end do
  201       end do
      ENDIF
C
C     REAL,PERIODIC SYNTHESIS
C
      CALL VRFFTB(LDIMF*N,M,FT,LDIMF*N,F,WORK)
C
C     POSTPROCESSING
C
      IF (OUTARY) THEN
         CALL VSSNFA(L,M,N*LDIMF,F,FT,C1,C2)
      ELSE
         CALL VSSNFA(L,M,N*LDIMF,FT,F,C1,C2)
      ENDIF
      RETURN
      END SUBROUTINE VSSINF
C
C
      SUBROUTINE VSSINQ(F,L,M,N,LDIMF,FT,C1,C2,C3,C4,WORK)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF*N,M),FT(LDIMF*N,M),
     *          C1(M),C2(M),C3(M),C4(M),WORK(M+15)
C
C     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
C           = .FALSE. IF TRANSFORMING THIRD INDEX
C
C
C     PREPROCESSING
C
      IF (TPOSE) THEN
         CALL VSSNQ1(L,M,N,LDIMF,F,FT,C1,C2)
      ELSE
         DO 100 I=1,Ldimf*n
            FT(I,1)=F(I,1)
  100       end do
         IF (MOD(M,2).EQ.0) THEN
            DO 110 I=1,Ldimf*n
               FT(I,M)=F(I,M)
  110          end do
         ENDIF
         DO 201 J=2,M-1,2
            JBY2=J/2
            DO 200 I=1,Ldimf*n
               FT(I,J)   =   F(I,J+1)*C1(JBY2)-F(I,J)*C2(JBY2)
               FT(I,J+1) = F(I,J+1)*(-C2(JBY2))-F(I,J)*C1(JBY2)
  200          end do
  201          end do
      ENDIF
C
C     REAL,PERIODIC SYNTHESIS
C
      CALL VRFFTB(LDIMF*N,M,FT,LDIMF*N,F,WORK)
C
C     POSTPROCESSING
C
      IF (OUTARY) THEN
         DO 301 J=1,M
         DO 300 I=1,Ldimf*n
            F(I,J)=C3(J)*FT(I,J)-C4(J)*FT(I,M+1-J)
  300       end do
  301       end do
      ELSE
         DO 311 J=1,M
         DO 310 I=1,Ldimf*n
            FT(I,J)=C3(J)*F(I,J)-C4(J)*F(I,M+1-J)
  310       end do
  311       end do
      ENDIF
      RETURN
      END SUBROUTINE VSSINQ
C
C
      SUBROUTINE VSSNB1(L,M,N,LDIMF,F,FT,C1,C2)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M,N), FT(LDIMF,N,M), C1(M), C2(M)
         DO 102 J=2,M
         DO 101 K=1,N
         DO 100 I=1,L
            FT(I,K,J)=C1(J)*F(I,J-1,K)-C2(J)*F(I,M-J+1,K)
  100       end do
  101       end do
  102       end do
         DO 201 K=1,N
         DO 200 I=1,L
            FT(I,K,1) = .5*F(I,M,K)
  200       end do
  201       end do
      RETURN
      END SUBROUTINE VSSNB1
C
C
      SUBROUTINE VSSNBA(L,M,LDIMF,F,FT)
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M), FT(LDIMF,M)
      DO 300 I=1,Ldimf
         F(I,1) = FT(I,1)
  300    end do
      IF (MOD(M,2).EQ.0) THEN
         DO 310 I=1,Ldimf
            F(I,M) = -FT(I,M)
  310       end do
      ENDIF
      DO 401 J=2,M-1,2
      DO 400 I=1,Ldimf
         F(I,J)    = -(FT(I,J+1)+FT(I,J))
         F(I,J+1)  =   FT(I,J)-FT(I,J+1)
  400    end do
  401    end do
      RETURN
      END SUBROUTINE VSSNBA
C
C
      SUBROUTINE VSSNF1(L,M,N,LDIMF,F,FT)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M,N), FT(LDIMF,N,M)
      DO 101 K=1,N
      DO 100 I=1,L
         FT(I,K,1)=F(I,1,K)
  100    end do
  101    end do
      IF (MOD(M,2).EQ.0) THEN
         DO 111 K=1,N
         DO 110 I=1,L
            FT(I,K,M)=-F(I,M,K)
  110       end do
  111       end do
      ENDIF
      DO 202 J=2,M-1,2
      DO 201 K=1,N
      DO 200 I=1,L
         FT(I,K,J)   = 0.5*(F(I,J+1,K)-F(I,J,K))
         FT(I,K,J+1) =-0.5*(F(I,J+1,K)+F(I,J,K))
  200    end do
  201    end do
  202    end do
      RETURN
      END SUBROUTINE VSSNF1
C
C
      SUBROUTINE VSSNFA(L,M,LDIMF,F,FT,C1,C2)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M), FT(LDIMF,M), C1(M), C2(M)
      DO 301 J=2,M
      DO 300 I=1,Ldimf
         F(I,J-1)=C1(J)*FT(I,J)+C2(J)*FT(I,M+2-J)
  300    end do
  301    end do
      DO 400 I=1,Ldimf
         F(I,M)=FT(I,1)
  400    end do
      RETURN
      END SUBROUTINE VSSNFA
C
C
      SUBROUTINE VSSNQ1(L,M,N,LDIMF,F,FT,C1,C2)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)   F(LDIMF,M,N), FT(LDIMF,N,M), C1(M), C2(M)
      DO 101 K=1,N
      DO 100 I=1,L
         FT(I,K,1)=F(I,1,K)
  100    end do
  101    end do
      IF (MOD(M,2).EQ.0) THEN
         DO 111 K=1,N
         DO 110 I=1,L
            FT(I,K,M)=F(I,M,K)
  110       end do
  111       end do
      ENDIF
      DO 202 J=2,M-1,2
         JBY2=J/2
         DO 201 K=1,N
         DO 200 I=1,L
            FT(I,K,J)   =   F(I,J+1,K)*C1(JBY2)-F(I,J,K)*C2(JBY2)
            FT(I,K,J+1) = -(F(I,J+1,K)*C2(JBY2)+F(I,J,K)*C1(JBY2))
  200       end do
  201       end do
  202       end do
      RETURN
      END SUBROUTINE VSSNQ1
C
C
      SUBROUTINE VRADF2 (MP,IDO,L1,CC,mdimc,CH,mdimch,WA1)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)     CH(Mdimch,IDO,2,L1)  ,CC(MDIMC,IDO,L1,2)     ,
     1                WA1(IDO)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+CC(M,1,K,2)
         CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2.LT.0) GOTO 107
      IF (IDO-2.EQ.0) GOTO 105
      IF (IDO-2.GT.0) GOTO 102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            DO 1003 M=1,MP
            CH(M,I,1,K) = CC(M,I,K,1)+(WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))
            CH(M,IC,2,K) = (WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-CC(M,I,K,1)
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+(WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)-(WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         DO 1006 M=1,MP
         CH(M,1,2,K) = -CC(M,IDO,K,2)
         CH(M,IDO,1,K) = CC(M,IDO,K,1)
 1006    CONTINUE
  106 CONTINUE
  107 RETURN
      END SUBROUTINE VRADF2
C
C
      SUBROUTINE VRADF3 (MP,IDO,L1,CC,mdimc,CH,MDIMCh,WA1,WA2)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)     CH(MDIMCh,IDO,3,L1)  ,CC(MDIMC,IDO,L1,3)     ,
     1                WA1(IDO)     ,WA2(IDO)
      ARG=2.*fsh19s()/3.
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,2)+CC(M,1,K,3))
         CH(M,1,3,K) = TAUI*(CC(M,1,K,3)-CC(M,1,K,2))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TAUR*
     1      (CC(M,1,K,2)+CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))+(TAUI*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*
     1       CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))-(TAUI*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*
     1       CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))))+(TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))))
            CH(M,IC,2,K) = (TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))))-(CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END SUBROUTINE VRADF3
C
C
      SUBROUTINE VRADF4 (MP,IDO,L1,CC,mdimc,CH,Mdimch,WA1,WA2,WA3)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)      CC(MDIMC,IDO,L1,4)   ,CH(Mdimch,IDO,4,L1)     ,
     1                WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)
      HSQT2=SQRT(2.0_EB)/2.
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = (CC(M,1,K,2)+CC(M,1,K,4))
     1      +(CC(M,1,K,1)+CC(M,1,K,3))
         CH(M,IDO,4,K) = (CC(M,1,K,1)+CC(M,1,K,3))
     1      -(CC(M,1,K,2)+CC(M,1,K,4))
         CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,3)
         CH(M,1,3,K) = CC(M,1,K,4)-CC(M,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2.LT.0) GOTO 107
      IF (IDO-2.EQ.0) GOTO 105
      IF (IDO-2.GT.0) GOTO 102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            DO 1003 M=1,MP
            CH(M,I-1,1,K) = ((WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4)))+(CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+
     1       WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))+(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,4,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))-(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))+(CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I,3,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,2,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3)))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         DO 1006 M=1,MP
            CH(M,IDO,1,K) = (HSQT2*(CC(M,IDO,K,2)-CC(M,IDO,K,4)))+
     1       CC(M,IDO,K,1)
            CH(M,IDO,3,K) = CC(M,IDO,K,1)-(HSQT2*(CC(M,IDO,K,2)-
     1       CC(M,IDO,K,4)))
            CH(M,1,2,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))-
     1       CC(M,IDO,K,3)
            CH(M,1,4,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))+
     1       CC(M,IDO,K,3)
 1006    CONTINUE
  106 CONTINUE
  107 RETURN
      END SUBROUTINE VRADF4 
C
C
      SUBROUTINE VRADF5 (MP,IDO,L1,CC,mdimc,CH,MDIMCh,WA1,WA2,WA3,WA4)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)    CC(MDIMC,IDO,L1,5)    ,CH(MDIMCh,IDO,5,L1)     ,
     1           WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
      ARG=2.*fsh19s()/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,5)+CC(M,1,K,2))+
     1    (CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TR11*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR12*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,3,K) = TI11*(CC(M,1,K,5)-CC(M,1,K,2))+TI12*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
         CH(M,IDO,4,K) = CC(M,1,K,1)+TR12*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR11*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,5,K) = TI12*(CC(M,1,K,5)-CC(M,1,K,2))-TI11*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5)))+((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I-1,3,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1      ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1       +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))+TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4)))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1     ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1      +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))-(TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI11*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,2,K) = (TI11*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I-1,5,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))+(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))-(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I,5,K) = (CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI12*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,4,K) = (TI12*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END SUBROUTINE VRADF5 
C
C
      SUBROUTINE VRADFG(MP,IDO,IP,L1,IDL1,CC,C1,C2,mdimc,
     *                  CH,CH2,MDIMCh,WA)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)       CH(MDIMCh,IDO,L1,IP)   ,CC(MDIMC,IDO,IP,L1)  ,
     1            C1(MDIMC,IDO,L1,IP)    ,C2(MDIMC,IDL1,IP),
     2                CH2(MDIMCh,IDL1,IP)           ,WA(IDO)
      TPI=2.*fsh19s()
      ARG = TPI/REAL(IP,EB)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         DO 1001 M=1,MP
         CH2(M,IK,1) = C2(M,IK,1)
 1001    CONTINUE
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            DO 1002 M=1,MP
            CH(M,1,K,J) = C1(M,1,K,J)
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               DO 1004 M=1,MP
               CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ)
     1           *C1(M,I,K,J)
               CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ)
     1           *C1(M,I-1,K,J)
 1004          CONTINUE
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               DO 1008 M=1,MP
               CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ)
     1           *C1(M,I,K,J)
               CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ)
     1           *C1(M,I-1,K,J)
 1008          CONTINUE
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               DO 1012 M=1,MP
               C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
               C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
               C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
 1012          CONTINUE
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               DO 1016 M=1,MP
               C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
               C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
               C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
 1016          CONTINUE
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         DO 1020 M=1,MP
         C2(M,IK,1) = CH2(M,IK,1)
 1020    CONTINUE
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            DO 1022 M=1,MP
            C1(M,1,K,J) = CH(M,1,K,J)+CH(M,1,K,JC)
            C1(M,1,K,JC) = CH(M,1,K,JC)-CH(M,1,K,J)
 1022       CONTINUE
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.
      AI1 = 0.
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            DO 1024 M=1,MP
            CH2(M,IK,L) = C2(M,IK,1)+AR1*C2(M,IK,2)
            CH2(M,IK,LC) = AI1*C2(M,IK,IP)
 1024       CONTINUE
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               DO 1025 M=1,MP
               CH2(M,IK,L) = CH2(M,IK,L)+AR2*C2(M,IK,J)
               CH2(M,IK,LC) = CH2(M,IK,LC)+AI2*C2(M,IK,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            DO 1028 M=1,MP
            CH2(M,IK,1) = CH2(M,IK,1)+C2(M,IK,J)
 1028       CONTINUE
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            DO 1030 M=1,MP
            CC(M,I,1,K) = CH(M,I,K,1)
 1030       CONTINUE
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            DO 1033 M=1,MP
            CC(M,I,1,K) = CH(M,I,K,1)
 1033       CONTINUE
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            DO 1036 M=1,MP
            CC(M,IDO,J2-2,K) = CH(M,1,K,J)
            CC(M,1,J2-1,K) = CH(M,1,K,JC)
 1036       CONTINUE
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               DO 1038 M=1,MP
               CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
               CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
               CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
 1038          CONTINUE
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               DO 1042 M=1,MP
               CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
               CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
               CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
 1042          CONTINUE
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END SUBROUTINE VRADFG
C
C
      SUBROUTINE VRFFTB(M,N,R,MDIMR,rt,WSAVE)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)       R(MDIMR,N),RT(M,N),WSAVE(N+15)
      IF (N .EQ. 1) RETURN
      CALL VRFTB1 (M,N,R,MDIMR,rt,WSAVE(1),WSAVE(N+1))
      RETURN
      END SUBROUTINE VRFFTB
C
C
      SUBROUTINE VRFTB1 (M,N,C,MDIMC,ch,WA,FAC)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)         CH(M,N), C(MDIMC,N), WA(N) ,FAC(15)
      NF = FAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = FAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL VRADB4 (M,IDO,L1,C,mdimc,CH,M,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL VRADB4 (M,IDO,L1,CH,m,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL VRADB2 (M,IDO,L1,C,mdimc,CH,M,WA(IW))
         GO TO 105
  104    CALL VRADB2 (M,IDO,L1,CH,m,C,MDIMC,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL VRADB3 (M,IDO,L1,C,mdimc,CH,M,WA(IW),WA(IX2))
         GO TO 108
  107    CALL VRADB3 (M,IDO,L1,CH,m,C,MDIMC,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
      CALL VRADB5 (M,IDO,L1,C,mdimc,CH,M,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110 CALL VRADB5 (M,IDO,L1,CH,m,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL VRADBG (M,IDO,IP,L1,IDL1,C,C,C,mdimc,CH,CH,M,WA(IW))
         GO TO 114
  113    CALL VRADBG (M,IDO,IP,L1,IDL1,CH,CH,CH,m,C,C,MDIMC,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      outary=.TRUE.
      if (nocopy) then
         scale=scale*sqrt(1.0_EB/real(n,EB))
         if (na.eq.1) then
            outary=.FALSE.
         endif
      else
         SCALE=SQRT(1.0_EB/REAL(N,EB))
         IF (NA .EQ. 0) GO TO 118
         DO 127 J=1,N
         DO 117 I=1,M
            C(I,J) = SCALE*CH(I,J)
  117    CONTINUE
  127    CONTINUE
         RETURN
  118    DO 129 J=1,N
         DO 119 I=1,M
            C(I,J)=SCALE*C(I,J)
  119    CONTINUE
  129    CONTINUE
      endif
      RETURN
      END SUBROUTINE VRFTB1
C
C
      SUBROUTINE VRADB2 (MP,IDO,L1,CC,mdimc,CH,MDIMCH,WA1)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)    CC(MDIMC,IDO,2,L1)    ,CH(MDIMCH,IDO,L1,2),
     1                WA1(IDO)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+CC(M,IDO,2,K)
         CH(M,1,K,2) = CC(M,1,1,K)-CC(M,IDO,2,K)
 1001     CONTINUE
  101 CONTINUE
      IF (IDO-2.LT.0) GOTO 107
      IF (IDO-2.EQ.0) GOTO 105
      IF (IDO-2.GT.0) GOTO 102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+CC(M,IC-1,2,K)
            CH(M,I,K,1) = CC(M,I,1,K)-CC(M,IC,2,K)
            CH(M,I-1,K,2) = WA1(I-2)*(CC(M,I-1,1,K)-CC(M,IC-1,2,K))
     1      -WA1(I-1)*(CC(M,I,1,K)+CC(M,IC,2,K))
            CH(M,I,K,2) = WA1(I-2)*(CC(M,I,1,K)+CC(M,IC,2,K))+WA1(I-1)
     1      *(CC(M,I-1,1,K)-CC(M,IC-1,2,K))
 1002          CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
          DO 1003 M=1,MP
         CH(M,IDO,K,1) = CC(M,IDO,1,K)+CC(M,IDO,1,K)
         CH(M,IDO,K,2) = -(CC(M,1,2,K)+CC(M,1,2,K))
 1003     CONTINUE
  106 CONTINUE
  107 RETURN
      END SUBROUTINE VRADB2
C
C
      SUBROUTINE VRADB3 (MP,IDO,L1,CC,mdimc,CH,MDIMCh,WA1,WA2)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)    CC(MDIMC,IDO,3,L1)    ,CH(MDIMCh,IDO,L1,3),
     1                WA1(IDO)   ,WA2(IDO)
      ARG=2.*fsh19s()/3.
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+2.*CC(M,IDO,2,K)
         CH(M,1,K,2) = CC(M,1,1,K)+(2.*TAUR)*CC(M,IDO,2,K)
     1   -(2.*TAUI)*CC(M,1,3,K)
         CH(M,1,K,3) = CC(M,1,1,K)+(2.*TAUR)*CC(M,IDO,2,K)
     1   +2.*TAUI*CC(M,1,3,K)
 1001     CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
            CH(M,I,K,1) = CC(M,I,1,K)+(CC(M,I,3,K)-CC(M,IC,2,K))
            CH(M,I-1,K,2) = WA1(I-2)*
     1 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-
     * (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
     2                   -WA1(I-1)*
     3 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))+
     * (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
            CH(M,I,K,2) = WA1(I-2)*
     4 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))+
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
     5                  +WA1(I-1)*
     6 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
              CH(M,I-1,K,3) = WA2(I-2)*
     7 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))+
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
     8                      -WA2(I-1)*
     9 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))-
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
            CH(M,I,K,3) = WA2(I-2)*
     1 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))-
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
     2                 +WA2(I-1)*
     3 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))+
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
 1002          CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END SUBROUTINE VRADB3
C
C
      SUBROUTINE VRADB4 (MP,IDO,L1,CC,mdimc,CH,Mdimch,WA1,WA2,WA3)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)    CC(MDIMC,IDO,4,L1)  ,CH(MDIMCh,IDO,L1,4)    ,
     1                WA1(IDO)  ,WA2(IDO)  ,WA3(IDO)
      SQRT2=SQRT(2.0_EB)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,3) = (CC(M,1,1,K)+CC(M,IDO,4,K))
     1   -(CC(M,IDO,2,K)+CC(M,IDO,2,K))
         CH(M,1,K,1) = (CC(M,1,1,K)+CC(M,IDO,4,K))
     1   +(CC(M,IDO,2,K)+CC(M,IDO,2,K))
         CH(M,1,K,4) = (CC(M,1,1,K)-CC(M,IDO,4,K))
     1   +(CC(M,1,3,K)+CC(M,1,3,K))
         CH(M,1,K,2) = (CC(M,1,1,K)-CC(M,IDO,4,K))
     1   -(CC(M,1,3,K)+CC(M,1,3,K))
 1001     CONTINUE
  101 CONTINUE
      IF (IDO-2.LT.0) GOTO 107
      IF (IDO-2.EQ.0) GOTO 105
      IF (IDO-2.GT.0) GOTO 102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = (CC(M,I-1,1,K)+CC(M,IC-1,4,K))
     1      +(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
            CH(M,I,K,1) = (CC(M,I,1,K)-CC(M,IC,4,K))
     1      +(CC(M,I,3,K)-CC(M,IC,2,K))
            CH(M,I-1,K,2)=WA1(I-2)*((CC(M,I-1,1,K)-CC(M,IC-1,4,K))
     1      -(CC(M,I,3,K)+CC(M,IC,2,K)))-WA1(I-1)
     1      *((CC(M,I,1,K)+CC(M,IC,4,K))+(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))
            CH(M,I,K,2)=WA1(I-2)*((CC(M,I,1,K)+CC(M,IC,4,K))
     1      +(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))+WA1(I-1)
     1      *((CC(M,I-1,1,K)-CC(M,IC-1,4,K))-(CC(M,I,3,K)+CC(M,IC,2,K)))
            CH(M,I-1,K,3)=WA2(I-2)*((CC(M,I-1,1,K)+CC(M,IC-1,4,K))
     1      -(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-WA2(I-1)
     1      *((CC(M,I,1,K)-CC(M,IC,4,K))-(CC(M,I,3,K)-CC(M,IC,2,K)))
            CH(M,I,K,3)=WA2(I-2)*((CC(M,I,1,K)-CC(M,IC,4,K))
     1      -(CC(M,I,3,K)-CC(M,IC,2,K)))+WA2(I-1)
     1      *((CC(M,I-1,1,K)+CC(M,IC-1,4,K))-(CC(M,I-1,3,K)
     1      +CC(M,IC-1,2,K)))
            CH(M,I-1,K,4)=WA3(I-2)*((CC(M,I-1,1,K)-CC(M,IC-1,4,K))
     1      +(CC(M,I,3,K)+CC(M,IC,2,K)))-WA3(I-1)
     1     *((CC(M,I,1,K)+CC(M,IC,4,K))-(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))
            CH(M,I,K,4)=WA3(I-2)*((CC(M,I,1,K)+CC(M,IC,4,K))
     1      -(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))+WA3(I-1)
     1      *((CC(M,I-1,1,K)-CC(M,IC-1,4,K))+(CC(M,I,3,K)+CC(M,IC,2,K)))
 1002          CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
               DO 1003 M=1,MP
         CH(M,IDO,K,1) = (CC(M,IDO,1,K)+CC(M,IDO,3,K))
     1   +(CC(M,IDO,1,K)+CC(M,IDO,3,K))
         CH(M,IDO,K,2) = SQRT2*((CC(M,IDO,1,K)-CC(M,IDO,3,K))
     1   -(CC(M,1,2,K)+CC(M,1,4,K)))
         CH(M,IDO,K,3) = (CC(M,1,4,K)-CC(M,1,2,K))
     1   +(CC(M,1,4,K)-CC(M,1,2,K))
         CH(M,IDO,K,4) = -SQRT2*((CC(M,IDO,1,K)-CC(M,IDO,3,K))
     1   +(CC(M,1,2,K)+CC(M,1,4,K)))
 1003          CONTINUE
  106 CONTINUE
  107 RETURN
      END SUBROUTINE VRADB4
C
C
      SUBROUTINE VRADB5 (MP,IDO,L1,CC,mdimc,CH,MDIMCh,WA1,WA2,WA3,WA4)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)    CC(MDIMC,IDO,5,L1)    ,CH(MDIMCh,IDO,L1,5),
     1             WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
      ARG=2.*fsh19s()/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
      DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+2.*CC(M,IDO,2,K)+2.*CC(M,IDO,4,K)
         CH(M,1,K,2) = (CC(M,1,1,K)+TR11*2.*CC(M,IDO,2,K)
     1   +TR12*2.*CC(M,IDO,4,K))-(TI11*2.*CC(M,1,3,K)
     1   +TI12*2.*CC(M,1,5,K))
         CH(M,1,K,3) = (CC(M,1,1,K)+TR12*2.*CC(M,IDO,2,K)
     1   +TR11*2.*CC(M,IDO,4,K))-(TI12*2.*CC(M,1,3,K)
     1   -TI11*2.*CC(M,1,5,K))
         CH(M,1,K,4) = (CC(M,1,1,K)+TR12*2.*CC(M,IDO,2,K)
     1   +TR11*2.*CC(M,IDO,4,K))+(TI12*2.*CC(M,1,3,K)
     1   -TI11*2.*CC(M,1,5,K))
         CH(M,1,K,5) = (CC(M,1,1,K)+TR11*2.*CC(M,IDO,2,K)
     1   +TR12*2.*CC(M,IDO,4,K))+(TI11*2.*CC(M,1,3,K)
     1   +TI12*2.*CC(M,1,5,K))
 1001          CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
      DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +(CC(M,I-1,5,K)+CC(M,IC-1,4,K))
            CH(M,I,K,1) = CC(M,I,1,K)+(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +(CC(M,I,5,K)-CC(M,IC,4,K))
            CH(M,I-1,K,2) = WA1(I-2)*((CC(M,I-1,1,K)+TR11*
     1      (CC(M,I-1,3,K)+CC(M,IC-1,2,K))+TR12
     1      *(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA1(I-1)*((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))+(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,2) = WA1(I-2)*((CC(M,I,1,K)+TR11*(CC(M,I,3,K)
     1      -CC(M,IC,2,K))+TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI11*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))+TI12
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))+WA1(I-1)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)
     1      +CC(M,IC-1,2,K))+TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))
     1      -(TI11*(CC(M,I,3,K)+CC(M,IC,2,K))+TI12
     1      *(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,3) = WA2(I-2)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1     -WA2(I-1)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,3) = WA2(I-2)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA2(I-1)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,4) = WA3(I-2)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA3(I-1)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      -(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,4) = WA3(I-2)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      -(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA3(I-1)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,5) = WA4(I-2)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA4(I-1)
     1      *((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))-(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,5) = WA4(I-2)
     1      *((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))-(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA4(I-1)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
 1002          CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END SUBROUTINE VRADB5 
C
C
      SUBROUTINE VRADBG (MP,IDO,IP,L1,IDL1,CC,C1,C2,mdimc,
     *                   CH,CH2,MDIMCh,WA)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

      REAL(EB)      CH(MDIMCh,IDO,L1,IP)    ,CC(MDIMC,IDO,IP,L1) ,
     1           C1(MDIMC,IDO,L1,IP)     ,C2(MDIMC,IDL1,IP),
     2                CH2(MDIMCh,IDL1,IP)       ,WA(IDO)
      TPI=2.*fsh19s()
      ARG = TPI/REAL(IP,EB)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            DO 1001 M=1,MP
            CH(M,I,K,1) = CC(M,I,1,K)
 1001       CONTINUE
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            DO 1004 M=1,MP
            CH(M,I,K,1) = CC(M,I,1,K)
 1004       CONTINUE
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            DO 1007 M=1,MP
            CH(M,1,K,J) = CC(M,IDO,J2-2,K)+CC(M,IDO,J2-2,K)
            CH(M,1,K,JC) = CC(M,1,J2-1,K)+CC(M,1,J2-1,K)
 1007       CONTINUE
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               DO 1009 M=1,MP
               CH(M,I-1,K,J) = CC(M,I-1,2*J-1,K)+CC(M,IC-1,2*J-2,K)
               CH(M,I-1,K,JC) = CC(M,I-1,2*J-1,K)-CC(M,IC-1,2*J-2,K)
               CH(M,I,K,J) = CC(M,I,2*J-1,K)-CC(M,IC,2*J-2,K)
               CH(M,I,K,JC) = CC(M,I,2*J-1,K)+CC(M,IC,2*J-2,K)
 1009          CONTINUE
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               DO 1013 M=1,MP
               CH(M,I-1,K,J) = CC(M,I-1,2*J-1,K)+CC(M,IC-1,2*J-2,K)
               CH(M,I-1,K,JC) = CC(M,I-1,2*J-1,K)-CC(M,IC-1,2*J-2,K)
               CH(M,I,K,J) = CC(M,I,2*J-1,K)-CC(M,IC,2*J-2,K)
               CH(M,I,K,JC) = CC(M,I,2*J-1,K)+CC(M,IC,2*J-2,K)
 1013          CONTINUE
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            DO 1017 M=1,MP
            C2(M,IK,L) = CH2(M,IK,1)+AR1*CH2(M,IK,2)
            C2(M,IK,LC) = AI1*CH2(M,IK,IP)
 1017       CONTINUE
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               DO 1018 M=1,MP
               C2(M,IK,L) = C2(M,IK,L)+AR2*CH2(M,IK,J)
               C2(M,IK,LC) = C2(M,IK,LC)+AI2*CH2(M,IK,JC)
 1018          CONTINUE
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            DO 1021 M=1,MP
            CH2(M,IK,1) = CH2(M,IK,1)+CH2(M,IK,J)
 1021       CONTINUE
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            DO 1023 M=1,MP
            CH(M,1,K,J) = C1(M,1,K,J)-C1(M,1,K,JC)
            CH(M,1,K,JC) = C1(M,1,K,J)+C1(M,1,K,JC)
 1023       CONTINUE
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               DO 1025 M=1,MP
               CH(M,I-1,K,J) = C1(M,I-1,K,J)-C1(M,I,K,JC)
               CH(M,I-1,K,JC) = C1(M,I-1,K,J)+C1(M,I,K,JC)
               CH(M,I,K,J) = C1(M,I,K,J)+C1(M,I-1,K,JC)
               CH(M,I,K,JC) = C1(M,I,K,J)-C1(M,I-1,K,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               DO 1029 M=1,MP
               CH(M,I-1,K,J) = C1(M,I-1,K,J)-C1(M,I,K,JC)
               CH(M,I-1,K,JC) = C1(M,I-1,K,J)+C1(M,I,K,JC)
               CH(M,I,K,J) = C1(M,I,K,J)+C1(M,I-1,K,JC)
               CH(M,I,K,JC) = C1(M,I,K,J)-C1(M,I-1,K,JC)
 1029          CONTINUE
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         DO 1033 M=1,MP
         C2(M,IK,1) = CH2(M,IK,1)
 1033    CONTINUE
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            DO 1034 M=1,MP
            C1(M,1,K,J) = CH(M,1,K,J)
 1034       CONTINUE
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               DO 1036 M=1,MP
               C1(M,I-1,K,J) = WA(IDIJ-1)*CH(M,I-1,K,J)-WA(IDIJ)*
     1          CH(M,I,K,J)
               C1(M,I,K,J) = WA(IDIJ-1)*CH(M,I,K,J)+WA(IDIJ)*
     1          CH(M,I-1,K,J)
 1036          CONTINUE
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               DO 1040 M=1,MP
               C1(M,I-1,K,J) = WA(IDIJ-1)*CH(M,I-1,K,J)-WA(IDIJ)*
     1          CH(M,I,K,J)
               C1(M,I,K,J) = WA(IDIJ-1)*CH(M,I,K,J)+WA(IDIJ)*
     1          CH(M,I-1,K,J)
 1040          CONTINUE
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END SUBROUTINE VRADBG
C
C
      SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
C***BEGIN PROLOGUE  SSWAP
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A5
C***KEYWORDS  BLAS,INTERCHANGE,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Interchange s.p vectors
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C       SX  input vector SY (unchanged if N .LE. 0)
C       SY  input vector SX (unchanged if N .LE. 0)
C
C     Interchange single precision SX and single precision SY.
C     For I = 0 to N-1, interchange  SX(LX+I*INCX) and SY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSWAP
C
      REAL(EB) SX(1),SY(1),STEMP1,STEMP2,STEMP3
C***FIRST EXECUTABLE STATEMENT  SSWAP
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) THEN
        IF(INCX-1.LT.0) GOTO 5
        IF(INCX-1.EQ.0) GOTO 20
        IF(INCX-1.GT.0) GOTO 60
        ENDIF
    5 CONTINUE
C
C       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        STEMP1 = SX(IX)
        SX(IX) = SY(IY)
        SY(IY) = STEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        STEMP1 = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        STEMP1 = SX(I)
        STEMP2 = SX(I+1)
        STEMP3 = SX(I+2)
        SX(I) = SY(I)
        SX(I+1) = SY(I+1)
        SX(I+2) = SY(I+2)
        SY(I) = STEMP1
        SY(I+1) = STEMP2
        SY(I+2) = STEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
C
C     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
      NS = N*INCX
        DO 70 I=1,NS,INCX
        STEMP1 = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP1
   70   CONTINUE
      RETURN
      END SUBROUTINE SSWAP
C
C
      SUBROUTINE h3csis(rs,rf,l,lbdcnd,ts,tf,m,mbdcnd,ps,pf,n,nbdcnd,
     *                  elmbda,ldimf,mdimf,ierror,save,w,hx,hy)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   save(-3:*),w(*),hx(0:*),hy(0:*)


      pi = fsh19s()

C                               CHECK FOR INVALID INPUT

      ierror = 0

cmcg  IF (rs.LT.0.0) THEN
cmcg      ierror = ierror + 1
cmcg      save(ierror) = 1.
cmcg  END IF

      IF (rf.LE.rs) THEN
          ierror = ierror + 1
          save(ierror) = 2.
      END IF

      IF (l.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 3.
      END IF

      IF (lbdcnd.LT.1 .OR. lbdcnd.GT.6) THEN
          ierror = ierror + 1
          save(ierror) = 4.
      END IF

cmcg  IF (rs.EQ.0.0 .AND. lbdcnd.LE.4) THEN
cmcg      ierror = ierror + 1
cmcg      save(ierror) = 5.
cmcg  END IF

      IF (rs.NE.0.0 .AND. lbdcnd.GE.5) THEN
          ierror = ierror + 1
          save(ierror) = 6.
      END IF

      IF (tf.LE.ts) THEN
          ierror = ierror + 1
          save(ierror) = 7.
      END IF

      IF (m.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 8.
      END IF

      IF (mbdcnd.LT.1 .OR. mbdcnd.GT.9) THEN
          ierror = ierror + 1
          save(ierror) = 9.
      END IF

cmcg  IF (ts.EQ.0. .AND. (mbdcnd.EQ.3.OR.mbdcnd.EQ.4.OR.
cmcg *    mbdcnd.EQ.8)) THEN
cmcg      ierror = ierror + 1
cmcg      save(ierror) = 10.
cmcg  END IF

      IF (tf.EQ.pi .AND. (mbdcnd.EQ.2.OR.mbdcnd.EQ.3.OR.
     *    mbdcnd.EQ.6)) THEN
          ierror = ierror + 1
          save(ierror) = 11.
      END IF

      IF (pf.LE.ps) THEN
          ierror = ierror + 1
          save(ierror) = 12.
      END IF

      IF (n.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 13.
      END IF

      IF (nbdcnd.LT.0 .OR. nbdcnd.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 14.
      END IF

      IF (ldimf.LT.l) THEN
          ierror = ierror + 1
          save(ierror) = 15.
      END IF

      IF (mdimf.LT.m) THEN
          ierror = ierror + 1
          save(ierror) = 16.
      END IF

      IF (mbdcnd.GT.4 .AND. (nbdcnd.EQ.1.OR.nbdcnd.EQ.2.OR.
     *    nbdcnd.EQ.4)) THEN
          ierror = ierror + 1
          save(ierror) = 17.
      END IF

      IF (lbdcnd.GT.4 .AND. .NOT. (mbdcnd.EQ.3.OR.mbdcnd.EQ.6.OR.
     *    mbdcnd.EQ.8.OR.mbdcnd.EQ.9)) THEN
          ierror = ierror + 1
          save(ierror) = 18.
      END IF

      IF (lbdcnd.GT.4 .AND. .NOT. (nbdcnd.EQ.0..OR.nbdcnd.EQ.3)) THEN
          ierror = ierror + 1
          save(ierror) = 19.
      END IF

      IF (ierror.NE.0) THEN
          RETURN
      ELSE
          save(1) = ierror
      END IF

C                               DEFINE GRID PARAMETERS

      dr = (rf-rs)/float(l)
      drby2 = dr/2.
      drsqr = 1./ (dr**2)
      dt = (tf-ts)/float(m)
      dtby2 = dt/2.
      dtsqr = 1./ (dt**2)
      dp = (pf-ps)/float(n)
      dpsqr = 1./ (dp**2)

C                               ALLOCATE SAVE ARRAY

      ial = 13
      ibl = ial + l
      icl = ibl + l
      idl = icl + l
      isl = idl + l
      iam = isl + l
      ibm = iam + m
      icm = ibm + m
      idm = icm + m
      ism = idm + m
      isvps = ism + m

C                               DEFINE AL,BL,CL,DL,SL ARRAYS IN
C                               ARRAY SAVE.  SL IS SYMMETRIZER FOR
C                               R-OPERATOR.
      DO 100 i = 1,l
          hxm = .5*(hx(i-1)+hx(i))
          hxp = .5*(hx(i)+hx(i+1))
cmcg      ri = rs + (i-.5)*dr
cmcg      save(idl+i) = 1./ri**2
          save(idl+i) = 1.
cmcg      save(ial+i) = drsqr* (ri-drby2)**2
cmcg      save(icl+i) = drsqr* (ri+drby2)**2
          save(ial+i) = drsqr/(hx(i)*hxm)
          save(icl+i) = drsqr/(hx(i)*hxp)
          save(ibl+i) = - (save(ial+i)+save(icl+i)) + elmbda/save(idl+i)
cmcg      save(isl+i) = ri**2
          save(isl+i) = 1.
  100 CONTINUE

C                               DEFINE BOUNDARY COEFFICIENTS

      SELECT CASE(lbdcnd)
      CASE(1:2) ; GOTO 110
      CASE(3:6) ; GOTO 120
      END SELECT

  110 CONTINUE
      save(ibl+1) = save(ibl+1) - save(ial+1)

      GO TO 130

  120 CONTINUE
      save(ibl+1) = save(ibl+1) + save(ial+1)
  130 CONTINUE

      SELECT CASE(lbdcnd)
      CASE(1)   ; GOTO 140
      CASE(2:3) ; GOTO 150
      CASE(4:5) ; GOTO 140
      CASE(6)   ; GOTO 150
      END SELECT

  140 CONTINUE
      save(ibl+l) = save(ibl+l) - save(icl+l)

      GO TO 160

  150 CONTINUE
      save(ibl+l) = save(ibl+l) + save(icl+l)
  160 CONTINUE

C                               DEFINE ARRAYS AM,BM,CM,DM,SM.  SM IS
C                               SYMMETRIZER FOR THETA-OPERATOR.
      DO 170 j = 1,m
          hym = .5*(hy(j-1)+hy(j))
          hyp = .5*(hy(j)+hy(j+1))
cmcg      tj = ts + (j-.5)*dt
cmcg      save(ism+j) = sin(tj)
          save(ism+j) = 1.     
cmcg      save(idm+j) = dpsqr/sin(tj)**2
          save(idm+j) = dpsqr
cmcg      save(iam+j) = dtsqr*sin(tj-dtby2)/sin(tj)
cmcg      save(icm+j) = dtsqr*sin(tj+dtby2)/sin(tj)
          save(iam+j) = dtsqr/(hy(j)*hym)
          save(icm+j) = dtsqr/(hy(j)*hyp)
          save(ibm+j) = - (save(iam+j)+save(icm+j))
  170 CONTINUE

C                               DEFINE BOUNDARY COEFFICIENTS

      SELECT CASE(mbdcnd)
      CASE(1:2) ; GOTO 180
      CASE(3:6) ; GOTO 190
      CASE(7)   ; GOTO 180
      CASE(8:9) ; GOTO 190
      END SELECT

  180 CONTINUE
      save(ibm+1) = save(ibm+1) - save(iam+1)

      GO TO 200

  190 CONTINUE
      save(ibm+1) = save(ibm+1) + save(iam+1)
  200 CONTINUE

      SELECT CASE(mbdcnd)
      CASE(1)   ; GOTO 210
      CASE(2:3) ; GOTO 220
      CASE(4:5) ; GOTO 210
      CASE(6:9) ; GOTO 220
      END SELECT

  210 CONTINUE
      save(ibm+m) = save(ibm+m) - save(icm+m)

      GO TO 230

  220 CONTINUE
      save(ibm+m) = save(ibm+m) + save(icm+m)
  230 CONTINUE

C                               INITIALIZE SOLVER ROUTINE S3CCIS

      CALL s3ccis(l,save(ial+1),save(ibl+1),save(icl+1),
     *            m,save(iam+1:iam+m),
     *            save(ibm+1),save(icm+1),save(idm+1),nbdcnd,n,ldimf,
     *            mdimf,ierr1,save(isvps+1),w)

C                               TEST ERROR FLAG FROM S3CCIS FOR
C                               INTERNAL ERROR

      IF (ierr1.NE.0) THEN
          save(1) = 99.
          ierror = 1
          RETURN
      END IF

C                               SCALE RADIAL COEFFICIENTS
      DO 240 i = 1,l
          save(ial+i) = save(ial+i)*save(idl+i)
          save(icl+i) = save(icl+i)*save(idl+i)
          save(ibl+i) = save(ibl+i)*save(idl+i)
  240 CONTINUE

C                               COMPUTE SCALING FOR SINGULAR PROBLEMS

      sum = 0.
      DO 250 j = 1,m
cmcg      sum = sum + save(ism+j)
          sum = sum + save(ism+j)*hy(j)
  250 CONTINUE

      s3 = 0.
      DO 260 i = 1,l
cmcg      s3 = s3 + save(isl+i)
          s3 = s3 + save(isl+i)*hx(i)
  260 CONTINUE
      save(11) = sum*s3*n

C                               RESTORE ARRAY DM

      DO 270 j = 1,m
          save(idm+j) = save(idm+j)/dpsqr
  270 CONTINUE

C                               SAVE PARAMETERS FOR HS3SPH IN SAVE ARRAY

      save(2) = dr
      save(3) = l
      save(4) = lbdcnd
      save(5) = dt
      save(6) = m
      save(7) = mbdcnd
      save(8) = dp
      save(9) = n
      save(10) = nbdcnd
      save(12) = elmbda

      save(-1) = kappa
      save(-2) = nmax 
      save(-3) = ikpwr

      RETURN
      END SUBROUTINE h3csis
C
C
      SUBROUTINE h3csss(bdrs,bdrf,bdts,bdtf,bdps,bdpf,ldimf,mdimf,f,
     *                  pertrb,save,w,hx,hy)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB) bdrs(mdimf,*),bdrf(mdimf,*),bdts(ldimf,*),bdtf(ldimf,*),
     *         bdps(ldimf,*),bdpf(ldimf,*),f(ldimf,mdimf,*),save(-3:*),
     *         w(*),hx(0:*),hy(0:*)


C                               CHECK VALUE OF IERROR (=SAVE(1)).
C                               IF NON-ZERO, RETURN.

      IF (save(1).NE.0.) RETURN

C                               GET PARAMETERS FOR H3CSSS FROM SAVE
C                               ARRAY WHERE THEY WERE STORED IN
C                               INITIALIZATION SUBROUTINE H3CSIS.

      kappa = save(-1)   ! Extra variables added to save array
      nmax  = save(-2)  
      ikpwr = save(-3) 
c
      dr = save(2)
      l = save(3)
      lbdcnd = save(4)

      dt = save(5)
      m = save(6)

      mbdcnd = save(7)

      dp = save(8)
      dpr = 1./dp
      dpsqr = 2./ (dp**2)
      n = save(9)
      nbdcnd = save(10)

      elmbda = save(12)

C                               ALLOCATE SAVE ARRAY

      ial = 13
      ibl = ial + l
      icl = ibl + l
      idl = icl + l
      isl = idl + l
      iam = isl + l
      ibm = iam + m
      icm = ibm + m
      idm = icm + m
      ism = idm + m
      isvps = ism + m

C                               ENTER BOUNDARY DATA FOR R-BOUNDARIES

      SELECT CASE(lbdcnd)
      CASE(1:2) ; GOTO 180
      CASE(3:4) ; GOTO 210
      CASE(5:6) ; GOTO 240
      END SELECT

  180 CONTINUE
      DO 200 k = 1,n
          DO 190 j = 1,m
              f(1,j,k) = f(1,j,k) - 2.*save(ial+1)*bdrs(j,k)
  190     CONTINUE
  200 CONTINUE

      GO TO 240

  210 CONTINUE
      DO 230 k = 1,n
          DO 220 j = 1,m
              f(1,j,k) = f(1,j,k) + dr*save(ial+1)*bdrs(j,k)
  220     CONTINUE
  230 CONTINUE

  240 CONTINUE
      SELECT CASE(lbdcnd)
      CASE(1)   ; GOTO 250
      CASE(2:3) ; GOTO 280
      CASE(4:5) ; GOTO 250
      CASE(6)   ; GOTO 280
      END SELECT

  250 CONTINUE
      DO 270 k = 1,n
          DO 260 j = 1,m
              f(l,j,k) = f(l,j,k) - 2.*save(icl+l)*bdrf(j,k)
  260     CONTINUE
  270 CONTINUE

      GO TO 310

  280 CONTINUE
      DO 300 k = 1,n
          DO 290 j = 1,m
              f(l,j,k) = f(l,j,k) - dr*save(icl+l)*bdrf(j,k)
  290     CONTINUE
  300 CONTINUE
  310 CONTINUE

C                               ENTER BOUNDARY DATA FOR THETA-BOUNDARIES

      SELECT CASE(mbdcnd)
      CASE(1:2) ; GOTO 320
      CASE(3:4) ; GOTO 350
      CASE(5:6) ; GOTO 380
      CASE(7)   ; GOTO 320
      CASE(8)   ; GOTO 350
      CASE(9)   ; GOTO 380
      END SELECT

  320 CONTINUE
      DO 340 k = 1,n
          DO 330 i = 1,l
              f(i,1,k) = f(i,1,k) - 2.*save(iam+1)*save(idl+i)*bdts(i,k)
  330     CONTINUE
  340 CONTINUE

      GO TO 380

  350 CONTINUE
      scal = dt*save(iam+1)
      DO 370 k = 1,n
          DO 360 i = 1,l
              f(i,1,k) = f(i,1,k) + scal*save(idl+i)*bdts(i,k)
  360     CONTINUE
  370 CONTINUE

  380 CONTINUE
      SELECT CASE(mbdcnd)
      CASE(1)   ; GOTO 390
      CASE(2:3) ; GOTO 420
      CASE(4:5) ; GOTO 390
      CASE(6)   ; GOTO 420
      CASE(7:9) ; GOTO 450
      END SELECT

  390 CONTINUE
      DO 410 k = 1,n
          DO 400 i = 1,l
              f(i,m,k) = f(i,m,k) - 2.*save(icm+m)*save(idl+i)*bdtf(i,k)
  400     CONTINUE
  410 CONTINUE

      GO TO 450

  420 CONTINUE
      scal = dt*save(icm+m)
      DO 440 k = 1,n
          DO 430 i = 1,l
              f(i,m,k) = f(i,m,k) - scal*save(idl+i)*bdtf(i,k)
  430     CONTINUE
  440 CONTINUE
  450 CONTINUE

C                               ENTER BOUNDARY DATA FOR PHI-BOUNDARIES

      SELECT CASE(nbdcnd+1)
      CASE(1)   ; GOTO 590
      CASE(2:3) ; GOTO 460
      CASE(4:5) ; GOTO 490
      END SELECT

  460 CONTINUE
      DO 480 j = 1,m
          scal = dpsqr*save(idm+j)
          DO 470 i = 1,l
              f(i,j,1) = f(i,j,1) - scal*save(idl+i)*bdps(i,j)
  470     CONTINUE
  480 CONTINUE

      GO TO 520

  490 CONTINUE
      DO 510 j = 1,m
          scal = dpr*save(idm+j)
          DO 500 i = 1,l
              f(i,j,1) = f(i,j,1) + scal*save(idl+i)*bdps(i,j)
  500     CONTINUE
  510 CONTINUE

  520 CONTINUE
      SELECT CASE(nbdcnd+1)
      CASE(1)   ; GOTO 590
      CASE(2)   ; GOTO 530
      CASE(3:4) ; GOTO 560
      CASE(5)   ; GOTO 530
      END SELECT

  530 CONTINUE
      DO 550 j = 1,m
          scal = dpsqr*save(idm+j)
          DO 540 i = 1,l
              f(i,j,n) = f(i,j,n) - scal*save(idl+i)*bdpf(i,j)
  540     CONTINUE
  550 CONTINUE

      GO TO 590

  560 CONTINUE
      DO 580 j = 1,m
          scal = dpr*save(idm+j)
          DO 570 i = 1,l
              f(i,j,n) = f(i,j,n) - scal*save(idl+i)*bdpf(i,j)
  570     CONTINUE
  580 CONTINUE
  590 CONTINUE

      pertrb = 0.
      ising = 0

C                               FOR SINGULAR PROBLEMS ADJUST DATA TO
C                               INSURE A SOLUTION WILL EXIST.  GO THRU
C                               THIS CODE TWICE: ISING=1 FOR CALCULATING
C                               PERTRB; ISING=2 FOR NORMALIZING SOLUTION
C                               AFTER IT IS COMPUTED.

      SELECT CASE(lbdcnd)
      CASE(1:2) ; GOTO 770
      CASE(3)   ; GOTO 620
      CASE(4:5) ; GOTO 770
      CASE(6)   ; GOTO 620
      END SELECT

  620 CONTINUE
      SELECT CASE(mbdcnd)
      CASE(1:2) ; GOTO 770
      CASE(3)   ; GOTO 630
      CASE(4:5) ; GOTO 770
      CASE(6)   ; GOTO 630
      CASE(7)   ; GOTO 770
      CASE(8:9) ; GOTO 630
      END SELECT

  630 CONTINUE
      SELECT CASE(nbdcnd+1)
      CASE(1)   ; GOTO 640
      CASE(2:3) ; GOTO 770
      CASE(4)   ; GOTO 640
      CASE(5)   ; GOTO 770
      END SELECT

  640 CONTINUE
      IF (elmbda.eq.0.) then
         goto 650
         ELSE
         goto 770
         ENDIF
  650 CONTINUE
      ising = 1
  660 CONTINUE
      DO 670 i = 1,l
          w(i) = 0.
  670 CONTINUE
      DO 700 k = 1,n
          DO 690 j = 1,m
              DO 680 i = 1,l
cmcg              w(i) = w(i) + save(ism+j)*f(i,j,k)
                  w(i) = w(i) + hy(j)*save(ism+j)*f(i,j,k)
  680         CONTINUE
  690     CONTINUE
  700 CONTINUE

      sum = 0.
      DO 710 i = 1,l
cmcg      sum = sum + w(i)*save(isl+i)
          sum = sum + w(i)*save(isl+i)*hx(i)
  710 CONTINUE

      pert = sum/save(11)
C                               ADJUST F ARRAY BY PERT

      DO 740 k = 1,n
          DO 730 j = 1,m
              DO 720 i = 1,l
                  f(i,j,k) = f(i,j,k) - pert
  720         CONTINUE
  730     CONTINUE
  740 CONTINUE

C                               IF NORMALIZING SOLUTION, RESTORE PERTRB
C                               AND JUMP TO END

      IF (ising.EQ.2) THEN
          pertrb = prtsav

          GO TO 850

      END IF

      prtsav = pert

  770 CONTINUE

C                               SCALE RIGHT SIDE OF EQUATION BEFORE
C                               CALL TO S3CCSS

      DO 800 k = 1,n
          DO 790 j = 1,m
              DO 780 i = 1,l
                  f(i,j,k) = f(i,j,k)*save(isl+i)
  780         CONTINUE
  790     CONTINUE
  800 CONTINUE

C                               SOLVE SYSTEM USING S3CCSS

      CALL s3ccss(ldimf,mdimf,f,save(isvps+1),w)

C                               IF A SINGULAR PROBLEM,
C                               RE-NORMALIZE SOLUTION (ISING=2)

      IF (ising.EQ.1) THEN
          ising = 2

          GO TO 660

      END IF

  850 CONTINUE

      RETURN
      END SUBROUTINE h3csss
C
C
      SUBROUTINE s3ccis(l,al,bl,cl,m,am,bm,cm,dm,nperod,n,ldimf,mdimf,
     *                  ierror,save,w)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB) al(l),bl(l),cl(l),am(m),bm(m),cm(m),dm(m),save(-3:*),w(*)



C                               CHECK FOR INVALID INPUT

      ierror = 0

      IF (l.LT.2) THEN
          ierror = ierror + 1
          save(ierror) = 1.
      END IF

      IF (m.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 2.
      END IF

      IF (ldimf.LT.l) THEN
          ierror = ierror + 1
          save(ierror) = 3.
      END IF

      DO 100 j = 2,m

          IF (am(j)*cm(j-1).LT.0.0) THEN
              ierror = ierror + 1
              save(ierror) = 5.

              GO TO 110

          END IF

  100 CONTINUE
  110 CONTINUE

      IF (mdimf.LT.m) THEN
          ierror = ierror + 1
          save(ierror) = 6.
      END IF

      IF (nperod.LT.0 .AND. nperod.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 7.
      END IF

      IF (n.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 8.
      END IF

      IF (ierror.NE.0) THEN
          RETURN
      ELSE
          save(1) = ierror
      END IF

      igrid = 2

      CALL fsh21s(igrid,l,al,bl,cl,m,am,bm,cm,dm,nperod,n,ierror,save,w)

      RETURN
      END SUBROUTINE s3ccis
C
C
      SUBROUTINE s3ccss(ldimf,mdimf,f,save,w)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   f(ldimf,mdimf,*),save(-3:*),w(*)


C                               CHECK VALUE OF IERROR (=SAVE(1)).
C                               IF NON-ZERO, RETURN.

      IF (save(1).NE.0.) RETURN

      CALL fsh15s(ldimf,mdimf,f,save,w)

      RETURN
      END SUBROUTINE s3ccss
C
C
      SUBROUTINE fsh15s(ldimy,mdimy,y,save,w)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   y(*),save(-3:*),w(*)

C                               SOLVER BASED ON CYCLIC REDUCTION AND
C                               FAST FOURIER TRANSFORMS

      l = save(2)
C     ml = save(3)              this variable is never used
      m = save(4)
      n = save(5)
      np = save(6)
      igrid = save(7)

C                               ALLOCATE SAVE ARRAY

      ial = 8
      ibl = ial + l
      icl = ibl + l
      iam = icl + l
      icm = iam + m
      icfz = icm + m

      IF (igrid.EQ.1) THEN
          iwsz = icfz + 2*n
      ELSE
          iwsz = icfz + 4*n
      END IF

      ib = iwsz + n + 16
      icf = ib + ((kappa-2)*ikpwr+kappa+6)*n

      IF (ldimy.EQ.l .AND. mdimy.EQ.m) THEN

C                               DATA ARRAY HAS NO HOLES, SO CALL SOLVER

          CALL fsh16s(igrid,l,m,n,np,save(ial),save(ibl),save(icl),
     *                save(iam),save(icm),save(icfz),save(iwsz),
     *                save(ib),save(icf),y,w,w(1+m),w(1+m+l*m))

      ELSE

C                               PACK DATA ARRAY, CALL SOLVER, AND UNPACK

          CALL fsh04s(l,m,n,ldimy,mdimy,l,y,w)

          IF (n.GT.1) THEN
              leny = l*m* (n+1)
          ELSE
              leny = l*m
          END IF

          CALL fsh16s(igrid,l,m,n,np,save(ial),save(ibl),save(icl),
     *                save(iam),save(icm),save(icfz),save(iwsz),
     *                save(ib),save(icf),w,w(1+leny),w(1+leny+m),y)

          CALL fsh05s(l,m,n,ldimy,mdimy,l,y,w)

      END IF

      RETURN
      END SUBROUTINE fsh15s
C
C
      SUBROUTINE fsh16s(igrid,l,m,n,np,al,bl,cl,am,cm,cfz,wsavez,b,
     *                  coef,f,w1,w2,ft)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+


      REAL(EB)   al(l),bl(l),cl(l),am(m),cm(m),cfz(4*n),wsavez(n+16),
     *          b(*),coef(*),f(l,m,n),w1(m),w2(l*m),ft(l,m,n)
      LOGICAL datary

C                               BEGIN SOLUTION

      nocopy=.TRUE.
      datary=.TRUE.
      scale=1.
      ifwrd = 1
      tpose=.FALSE.
      ldimft=l
  100 CONTINUE

      IF (n.NE.1) THEN

C                               TRANSFORM IN Z

        iptr=7+2*(ifwrd-1)
        IF (datary) THEN
            call fsh26s(igrid,ifwrd,np,l,n,m,ldimft,f,ft,cfz,wsavez)
            datary=outary
        ELSE
            call fsh26s(igrid,ifwrd,np,l,n,m,ldimft,ft,f,cfz,wsavez)
            datary=.not.outary
        ENDIF

      END IF

      IF (ifwrd.NE.1) GO TO 900
      ib = 1
      icf = 1
      DO 890 k = 1,n
         IF (datary) THEN
            CALL fsh17s(am,cm,l,al,bl,cl,l,f(1,1,k),
     *                  b(ib),coef(icf),w1,w2,ft)
         ELSE
            CALL fsh17s(am,cm,l,al,bl,cl,l,ft(1,1,k),
     *                  b(ib),coef(icf),w1,w2,f)
         ENDIF
         ib = ib + (kappa-2)*ikpwr + kappa + 6
         icf = icf + (2*kappa-4.5)*ikpwr + kappa + 8
  890 CONTINUE
      ifwrd = 2

      GO TO 100

  900 CONTINUE

      IF (datary) THEN
         DO k=1,n
         DO j=1,m
         DO i=1,l
            f(i,j,k)=scale*f(i,j,k)
         ENDDO
         ENDDO
         ENDDO
      ELSE
         DO k=1,n
         DO j=1,m
         DO i=1,l
            f(i,j,k)=scale*ft(i,j,k)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE fsh16s
C
C
      SUBROUTINE fsh17s(an,cn,m,am,bm,cm,idimy,y,b,coef,rt,ws,d)

C                               SUBROUTINE FSH17S SOLVES THE LINEAR
C                               SYSTEM BY VECTORIZED CYCLIC REDUCTION

      REAL(EB)   an(nmax),cn(nmax),am(m),bm(m),cm(m),b(*),
     *          y(idimy,nmax),ws(nmax,m),rt(nmax),coef(*),d(nmax,m)

C                               LET KAPPA = LOG2(NMAX) + 1, THEN
C LENGTH OF B ARRAY = (KAPPA-2)*2**(KAPPA+1) + KAPPA + 5
C LENGTH OF COEF ARRAY = (2*KAPPA-4.5)*2**(KAPPA+1) + KAPPA + 8


C                               BEGIN REDUCTION PHASE

      kdo = kappa - 1
      if = 2**kappa
      ic = 0
      DO 210 ir = 0,kappa - 2
          irm1 = ir - 1
          i2 = 2**ir
          i3 = i2 + i2/2
          i4 = i2 + i2
          is = 0
          DO 130 i = i2,nmax,i4
              CALL fsh10s(i,ir,im2,nm2)
              DO 100 ip = 1,nm2
                  rt(is+ip) = b(im2+ip-1)
  100         CONTINUE
              DO 120 ip = 1,nm2
                  DO 110 j = 1,m
                      ws(is+ip,j) = y(j,i)
  110             CONTINUE
  120         CONTINUE
              is = is + nm2
  130     CONTINUE
          CALL fsh18s(is,m,am,bm,cm,rt,ws,d)
          is = 0
          DO 200 i = i4,nmax,i4
              ipi2 = i + i2
              imi2 = i - i2
              CALL fsh10s(imi2,ir,im2,nm2)
              DO 160 ip = 1,nm2
                  DO 150 j = 1,m
                      y(j,i) = y(j,i) + coef(ic+ip)*ws(is+ip,j)
  150             CONTINUE
  160         CONTINUE
              ic = ic + nm2

              IF (ipi2.GT.nmax) GO TO 200
              is = is + nm2
              CALL fsh10s(ipi2,ir,ip2,np2)
              DO 190 ip = 1,np2
                  DO 180 j = 1,m
                      y(j,i) = y(j,i) + coef(ic+ip)*ws(is+ip,j)
  180             CONTINUE
  190         CONTINUE
              ic = ic + np2
  200     CONTINUE
  210 CONTINUE

C                               BEGIN BACK SUBSTITUTION PHASE

      DO 530 ir = kdo,0,-1
          irm1 = ir - 1
          i2 = 2**ir
          i1 = i2/2
          i4 = i2 + i2

          IF (ir.EQ.kdo) GO TO 430

          IF (ir.NE.0) GO TO 280
          DO 270 i = 1,nmax,2

              IF (i.NE.1) GO TO 230
              DO 220 j = 1,m
                  y(j,i) = y(j,i) - cn(i)*y(j,i+1)
  220         CONTINUE

              GO TO 270

  230         CONTINUE
              IF (i.NE.nmax) GO TO 250
              DO 240 j = 1,m
                  y(j,i) = y(j,i) - an(i)*y(j,i-1)
  240         CONTINUE

              GO TO 270

  250         CONTINUE
              DO 260 j = 1,m
                  y(j,i) = y(j,i) - an(i)*y(j,i-1) - cn(i)*y(j,i+1)
  260         CONTINUE
  270     CONTINUE

          GO TO 430

  280     CONTINUE
          is = 0
          DO 360 i = i2,nmax,i4
              imi1 = i - i1
              ipi1 = i + i1
              imi2 = i - i2
              ipi2 = i + i2

              IF (i.EQ.i2) GO TO 320
              CALL fsh10s(imi1,irm1,im1,nm1)
              DO 290 ip = 1,nm1
                  rt(is+ip) = b(im1+ip-1)
  290         CONTINUE
              DO 310 ip = 1,nm1
                  DO 300 j = 1,m
                      ws(is+ip,j) = y(j,imi2)
  300             CONTINUE
  310         CONTINUE
              is = is + nm1

  320         CONTINUE
              IF (ipi2.GT.nmax) GO TO 360
              CALL fsh10s(ipi1,irm1,ip1,np1)
              DO 330 ip = 1,np1
                  rt(is+ip) = b(ip1+ip-1)
  330         CONTINUE
              DO 350 ip = 1,np1
                  DO 340 j = 1,m
                      ws(is+ip,j) = y(j,ipi2)
  340             CONTINUE
  350         CONTINUE
              is = is + np1
  360     CONTINUE
          CALL fsh18s(is,m,am,bm,cm,rt,ws,d)
          is = 0
          DO 420 i = i2,nmax,i4
              imi1 = i - i1
              ipi1 = i + i1

              IF (i.EQ.i2) GO TO 390
              CALL fsh10s(imi1,irm1,im1,nm1)
              DO 380 ip = 1,nm1
                  DO 370 j = 1,m
                      y(j,i) = y(j,i) + coef(ic+ip)*ws(is+ip,j)
  370             CONTINUE
  380         CONTINUE
              ic = ic + nm1
              is = is + nm1

  390         CONTINUE
              IF (i+i2.GT.nmax) GO TO 420
              CALL fsh10s(ipi1,irm1,ip1,np1)
              DO 410 ip = 1,np1
                  DO 400 j = 1,m
                      y(j,i) = y(j,i) + coef(ic+ip)*ws(is+ip,j)
  400             CONTINUE
  410         CONTINUE
              ic = ic + np1
              is = is + np1
  420     CONTINUE
  430     CONTINUE
          is = 0
          DO 470 i = i2,nmax,i4
              CALL fsh10s(i,ir,iz,nz)
              DO 440 ip = 1,nz
                  rt(is+ip) = b(iz+ip-1)
  440         CONTINUE
              DO 460 ip = 1,nz
                  DO 450 j = 1,m
                      ws(is+ip,j) = y(j,i)
  450             CONTINUE
  460         CONTINUE
              is = is + nz
  470     CONTINUE
          CALL fsh18s(is,m,am,bm,cm,rt,ws,d)
          is = 0
          DO 520 i = i2,nmax,i4
              CALL fsh10s(i,ir,iz,nz)
              DO 490 j = 1,m
                  y(j,i) = coef(ic+1)*ws(is+1,j)
  490         CONTINUE
              DO 510 ip = 2,nz
                  DO 500 j = 1,m
                      y(j,i) = y(j,i) + coef(ic+ip)*ws(is+ip,j)
  500             CONTINUE
  510         CONTINUE
              ic = ic + nz
              is = is + nz
  520     CONTINUE
  530 CONTINUE
      RETURN
      END SUBROUTINE fsh17s
C
C
      SUBROUTINE fsh18s(ideg,m,a,b,c,xl,y,d)

      REAL(EB)   a(m),b(m),c(m),y(nmax,m),d(nmax,m),xl(ideg)

C                              FSH18S IS A VECTORIZED TRIDIAGONAL SOLVER
C                              THAT ALSO DOES FACTORIZATION

      mm1 = m - 1
      DO 100 k = 1,ideg
          d(k,1) = c(1)/ (b(1)-xl(k))
          y(k,1) = y(k,1)/ (b(1)-xl(k))
  100 CONTINUE
      DO 120 i = 2,mm1
          DO 110 k = 1,ideg
              d(k,i) = c(i)/ (b(i)-xl(k)-a(i)*d(k,i-1))
              y(k,i) = (y(k,i)-a(i)*y(k,i-1))/
     *                 (b(i)-xl(k)-a(i)*d(k,i-1))
  110     CONTINUE
  120 CONTINUE
      DO 130 k = 1,ideg
          d(k,m) = y(k,m) - a(m)*y(k,mm1)
          y(k,m) = b(m) - xl(k) - a(m)*d(k,mm1)
  130 CONTINUE
      DO 140 k = 1,ideg

C                               Y(K,M) = CVMGZ(0.,D(K,M)/Y(K,M),Y(K,M))
C                               ON A CRAY-1

          IF (y(k,m).NE.0.) y(k,m) = d(k,m)/y(k,m)
  140 CONTINUE
      DO 160 i = m - 1,1,-1
          DO 150 k = 1,ideg
              y(k,i) = y(k,i) - d(k,i)*y(k,i+1)
  150     CONTINUE
  160 CONTINUE
      RETURN
      END SUBROUTINE fsh18s
C
C
      SUBROUTINE fsh21s(igrid,l,al,bl,cl,m,am,bm,cm,dm,nperod,n,ierror,
     *                  save,w)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+


      REAL(EB) al(l),bl(l),cl(l),am(m),bm(m),cm(m),dm(m),save(-3:*),w(*)

C                               THIS ROUTINE INITIALIZES SOLVER BASED ON
C                               CYCLIC REDUCTION AND FAST FOURIER TRANS.

C                               COMPUTE SMALLEST INTEGER KAPPA SUCH THAT
C                                         2**KAPPA .GE. M+1

      kappa = 2
      ml = 4
  110 CONTINUE

      IF (ml.LT.m+1) THEN
          kappa = kappa + 1
          ml = 2*ml

          GO TO 110

      END IF

      ikpwr = 2*ml
      ml = ml - 1

      nmax = m
      np = nperod + 1
      pi = fsh19s()

C                               ALLOCATE SAVE ARRAY

      ial = 8
      ibl = ial + l
      icl = ibl + l
      iam = icl + l
      icm = iam + m
      icfz = icm + m

      IF (igrid.EQ.1) THEN
          iwsz = icfz + 2*n
      ELSE
          iwsz = icfz + 4*n
      END IF

      ib = iwsz + n + 16
      icf = ib + ((kappa-2)*ikpwr+kappa+6)*n

C                               COMPUTE TRANSFORM ROOTS IN K-DIRECTION
C                               AND STORE IN W(1),...,W(N)

      IF (n.EQ.1) THEN

          w(1) = 0.

      ELSE

          IF (igrid.EQ.1) THEN
              nrdel = ((np-1)* (np-3)* (np-5))/3
              del = pi/ (2.* (n+nrdel))
          ELSE
              del = pi/ (2*n)
          END IF

          if (np.eq.1) then
          w(1) = 0.
          w(n) = -4.
          DO k = 2,n - 1,2
              w(k) = -4.*sin(k*del)**2
              w(k+1) = w(k)
              enddo
          endif

          if (np.eq.2) then
          DO k = 1,n
              w(k) = -4.*sin(k*del)**2
              enddo
          endif

          if (np.eq.3 .or. np.eq.5) then
          DO k = 1,n
              w(k) = -4.*sin((k-.5)*del)**2
              enddo
          endif

          if (np.eq.4) then
          DO k = 1,n
              w(k) = -4.*sin((k-1)*del)**2
              enddo
          endif

      END IF

C                               INITIALIZE FFT TRANSFORMS AND
C                               PRE-PROCESSING COEFFICIENTS IN K

      IF (n.NE.1) THEN

          SELECT CASE(igrid)
          CASE(1)   ; GOTO 210
          CASE(2)   ; GOTO 220
          END SELECT

  210     CONTINUE
          SELECT CASE(np)
          CASE(1)   ; GOTO 230
          CASE(2)   ; GOTO 240
          CASE(3)   ; GOTO 250
          CASE(4)   ; GOTO 270
          CASE(5)   ; GOTO 280
          END SELECT

  220     CONTINUE
          SELECT CASE(np)
          CASE(1)   ; GOTO 230
          CASE(2)   ; GOTO 250
          CASE(3)   ; GOTO 260
          CASE(4)   ; GOTO 280
          CASE(5)   ; GOTO 290
          END SELECT

  230     CONTINUE
          CALL vsrfti(n,save(iwsz))

          GO TO 300

  240     CONTINUE
          CALL vsinti(n,save(icfz),save(iwsz))

          GO TO 300

  250     CONTINUE
          CALL VSCOSI(n,save(icfz),save(icfz+n),save(iwsz))

          GO TO 300

  260     CONTINUE
ceb          CALL vssnqi(n,save(icfz),save(icfz+n),save(icfz+2*n),
ceb     *                save(icfz+3*n),save(iwsz))
          CALL VSCSQI(n,save(icfz),save(icfz+n),save(icfz+2*n),
     *                save(icfz+3*n),save(iwsz))

          GO TO 300

  270     CONTINUE
          CALL vcosti(n,save(icfz),save(iwsz))

          GO TO 300

  280     CONTINUE
          CALL vscosi(n,save(icfz),save(icfz+n),save(iwsz))

          GO TO 300

  290     CONTINUE
          CALL vscsqi(n,save(icfz),save(icfz+n),save(icfz+2*n),
     *                save(icfz+3*n),save(iwsz))
  300     CONTINUE
      END IF

      DO 320 k = 1,n

C                               COMPUTE NEW BM ARRAY AND STORE IN
C                               W(N+1),...,W(N+M)
          DO 310 j = 1,m
              w(n+j) = bm(j) + w(k)*dm(j)
  310     CONTINUE

C                               SUBROUTINE FSH07S COMPUTES THE ROOTS OF
C                               THE B POLYNOMIALS

          CALL fsh07s(ierr1,am,w(n+1),cm,save(ib),w(m+n+1),
     *                w(3*m+n+1))


          IF (ierr1.NE.0) THEN
              ierror = 1
              save(1) = ierr1
              RETURN
          END IF

C                               FSH08S COMPUTES COEFFICIENTS OF PARTIAL
C                               FRACTION EXPANSIONS

          CALL fsh08s(am,cm,save(ib),save(icf),w(m+n+1))

          ib = ib + (kappa-2)*ikpwr + kappa + 6
          icf = icf + (2*kappa-4.5)*ikpwr + kappa + 8
  320 CONTINUE

C                               SAVE QUANTITIES FOR USE IN SOLVERS

      save(2) = l
C     save(3) = ml              this variable never used in solver
      save(4) = m
      save(5) = n
      save(6) = np
      save(7) = igrid
C
      DO 330 i = 1,l
          save(ial+i-1) = al(i)
          save(ibl+i-1) = bl(i)
          save(icl+i-1) = cl(i)
  330 CONTINUE
C
      DO 340 j = 1,m
          save(iam+j-1) = am(j)
          save(icm+j-1) = cm(j)
  340 CONTINUE

      RETURN
      END SUBROUTINE fsh21s
C
C
      SUBROUTINE fsh07s(ierror,an,bn,cn,b,ah,bh)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C                               FSH07S COMPUTES THE ROOTS OF THE B
C                               POLYNOMIALS USING SUBROUTINE FSH14S
C                               WHICH IS A MODIFICATION THE EISPACK
C                               SUBROUTINE TQLRAT.  IERROR IS SET TO 4
C                               IF EITHER FSH14S FAILS OR IF
C                               A(J+1)*C(J) .LT. 0 FOR SOME J.
C                               AH,BH ARE TEMPORARY WORK ARRAYS.

      REAL(EB)   an(*),bn(*),cn(*),b(*),ah(*),bh(*)

      ierror = 0
      DO 100 j = 2,nmax
          arg = an(j)*cn(j-1)
          b(j) = sign(sqrt(arg),an(j))
  100 CONTINUE
      if = 2**kappa
      kdo = kappa - 1
      DO 160 l = 1,kdo
          ir = l - 1
          i2 = 2**ir
          i4 = i2 + i2
          ipl = i4 - 1
          ifd = if - i4
          DO 150 i = i4,ifd,i4
              CALL fsh10s(i,l,ib,nb)
              IF (nb.le.0) goto 160
              IF (nb.gt.0) goto 110
  110         CONTINUE
              js = i - ipl
              jf = js + nb - 1
              ls = 0
              DO 120 j = js,jf
                  ls = ls + 1
                  bh(ls) = bn(j)
                  ah(ls) = b(j)
  120         CONTINUE
              CALL fsh14s(nb,bh,ah,ierror)
              IF (ierror.eq.0) then
                 goto 130
                 else
                 goto 190
                 endif
  130         CONTINUE
              lh = ib - 1
              DO 140 j = 1,nb
                  lh = lh + 1
                  b(lh) = -bh(j)
  140         CONTINUE
  150     CONTINUE
  160 CONTINUE
      DO 170 j = 1,nmax
          b(j) = -bn(j)
  170 CONTINUE
      RETURN
  190 CONTINUE
      ierror = 4
      RETURN
      END SUBROUTINE fsh07s
C
C
      SUBROUTINE fsh08s(an,cn,b,coef,t)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C                               SUBROUTINE FSH08S GENERATES COEFFICIENTS
C                               FOR PARTIAL FRACTION EXPANSION OF MATRIX
C                               POLYNOMIAL PRODUCTS
      REAL(EB)   an(nmax),cn(nmax),b(*),coef(*),t(nmax,4),dum(0:0)

      if = 2**kappa
      kdo = kappa - 1
      is = 1
      DO 130 ir = 0,kappa - 2
          irm1 = ir - 1
          i2 = 2**ir
          i3 = i2 + i2/2
          i4 = i2 + i2
          DO 120 i = i4,nmax,i4
              imi2 = i - i2
              imi3 = i - i3
              CALL fsh09s(i,ir,idxa,na)
              CALL fsh10s(imi2,ir,im2,nm2)
              CALL fsh10s(imi3,irm1,im3,nm3)
              if (im3.lt.1) im3=1
              CALL fsh12s(ir,na,an(idxa),nm2,b(im2),nm3,b(im3),0,dum,
     *                    coef(is),t,t(1,2),t(1,3),t(1,4))
              is = is + nm2
              ipi2 = i + i2

              IF (ipi2.GT.nmax) GO TO 120
              ipi3 = i + i3
              CALL fsh11s(i,ir,idxc,nc)
              CALL fsh10s(ipi2,ir,ip2,np2)
              CALL fsh10s(ipi3,irm1,ip3,np3)
              if (ip3.le.0) ip3 = 1   ! Prevent out of bounds
              CALL fsh12s(ir,nc,cn(idxc),np2,b(ip2),np3,b(ip3),0,dum,
     *                    coef(is),t,t(1,2),t(1,3),t(1,4))
              is = is + np2
  120     CONTINUE
  130 CONTINUE

C                               BEGIN BACK SUBSTITUTION PHASE

      DO 190 ir = kdo,0,-1
          irm1 = ir - 1
          i2 = 2**ir
          i1 = i2/2
          i4 = i2 + i2

          IF (ir.EQ.kdo) GO TO 160

          IF (ir.EQ.0) GO TO 160
          DO 150 i = i2,nmax,i4

              IF (i.EQ.i2) GO TO 140
              imi1 = i - i1
              CALL fsh09s(i,ir,idxa,na)
              CALL fsh10s(imi1,irm1,im1,nm1)
              CALL fsh12s(irm1,na,an(idxa),nm1,b(im1),0,dum,0,dum,
     *                    coef(is),t,t(1,2),t(1,3),t(1,4))
              is = is + nm1

  140         CONTINUE
              IF (i+i2.GT.nmax) GO TO 150
              ipi1 = i + i1
              CALL fsh11s(i,ir,idxc,nc)
              CALL fsh10s(ipi1,irm1,ip1,np1)
              CALL fsh12s(irm1,nc,cn(idxc),np1,b(ip1),0,dum,0,dum,
     *                    coef(is),t,t(1,2),t(1,3),t(1,4))
              is = is + np1
  150     CONTINUE
  160     CONTINUE
          DO 180 i = i2,nmax,i4
              imi1 = i - i1
              ipi1 = i + i1
              CALL fsh10s(i,ir,iz,nz)
              CALL fsh10s(imi1,irm1,im1,nm1)
              CALL fsh10s(ipi1,irm1,ip1,np1)
              if (im1.lt.1) im1=1
              if (ip1.lt.1) ip1=1
              CALL fsh12s(ir,0,dum,nz,b(iz),nm1,b(im1),np1,b(ip1),
     *                    coef(is),t,t(1,2),t(1,3),t(1,4))
              is = is + nz
  180     CONTINUE
  190 CONTINUE
      RETURN
      END SUBROUTINE fsh08s
C
C
      SUBROUTINE fsh09s(i,ir,idxa,na)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+


      na = 2**ir
      idxa = i - na + 1

      IF (i-nmax.le.0) goto 110
      IF (i-nmax.gt.0) goto 100
  100 CONTINUE
      na = 0
  110 CONTINUE
      RETURN
      END SUBROUTINE fsh09s
C
C
      SUBROUTINE fsh10s(i,ir,idx,idp)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+


C                               B(IDX) IS THE LOCATION OF THE FIRST ROOT
C                               OF THE B(I,IR) POLYNOMIAL


      idp = 0
      idx = 0

      IF (ir.lt.0) goto 160
      IF (ir.eq.0) goto 100
      IF (ir.gt.0) goto 120

  100 CONTINUE
      IF (i-nmax.le.0) goto 110
      IF (i-nmax.gt.0) goto 160
  110 CONTINUE
      idx = i
      idp = 1
      RETURN
  120 CONTINUE
      izh = 2**ir
      id = i - izh - izh
      idx = id + id + (ir-1)*ikpwr + ir + (ikpwr-i)/izh + 4
      ipl = izh - 1
      idp = izh + izh - 1

      IF (i-ipl-nmax.le.0) goto 140
      IF (i-ipl-nmax.gt.0) goto 130
  130 CONTINUE
      idp = 0
      RETURN

  140 CONTINUE
      IF (i+ipl-nmax.le.0) goto 160
      IF (i+ipl-nmax.gt.0) goto 150
  150 CONTINUE
      idp = nmax + ipl - i + 1
  160 CONTINUE
      RETURN
      END SUBROUTINE fsh10s
C
C
      SUBROUTINE fsh11s(i,ir,idxc,nc)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+


      nc = 2**ir
      idxc = i

      IF (idxc+nc-1-nmax.le.0) goto 110
      IF (idxc+nc-1-nmax.gt.0) goto 100
  100 CONTINUE
      nc = 0
  110 CONTINUE
      RETURN
      END SUBROUTINE fsh11s
C
C
      SUBROUTINE fsh12s(ir,na,a,nc,c,nb1,b1,nb2,b2,coef,w1,w2,w3,w4)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

C                               SUBROUTINE FSH12S ACTUALLY GENERATES THE
C                               COEFFICIENTS

      REAL(EB) a(na),c(nc),b1(nb1),b2(nb2),coef(nc),w1(*),w2(*),w3(*),
     *         w4(*)

C                               MAXIMUM LENGTH OF W1, W2, W3, W4 IS NM.


      nb = nb1 + nb2
      DO 110 ip = 1,nb1
          w1(ip) = b1(ip)
  110 CONTINUE
      DO 120 ip = 1,nb2
          w1(nb1+ip) = b2(ip)
  120 CONTINUE

C                               MERGE TWO STRINGS

      CALL fsh13s(w1,0,nb1,nb1,nb2,nb)
      DO 130 ip = 1,nb
          w4(ip) = w1(nb+ip)
  130 CONTINUE
      m = max0(na,nc,nb)
      DO 140 i = 1,na
          w1(i) = a(i)
  140 CONTINUE
      DO 150 i = na + 1,m
          w1(i) = 1.
  150 CONTINUE
      DO 160 i = nb + 1,m
          w2(i) = 1.
  160 CONTINUE
      DO 170 i = nc + 1,m
          w3(i) = 1.
  170 CONTINUE
      sign = 1.

      IF (ir.EQ.0 .AND. na.EQ.1) sign = -1.
      DO 210 k = 1,nc
          DO 180 i = 1,nb
              w2(i) = c(k) - w4(i)
  180     CONTINUE
          DO 190 i = 1,nc
              w3(i) = c(k) - c(i)
  190     CONTINUE
          w3(k) = 1.
          coef(k) = sign
          DO 200 i = 1,m
              coef(k) = coef(k)*w1(i)*w2(i)/w3(i)
  200     CONTINUE
  210 CONTINUE
      RETURN
      END SUBROUTINE fsh12s
C
C
      SUBROUTINE fsh13s(tcos,i1,m1,i2,m2,i3)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB)   tcos(*)

C                               THIS SUBROUTINE MERGES TWO ASCENDING
C                               STRINGS OF NUMBERS IN THE ARRAY TCOS.
C                               THE FIRST STRING IS OF LENGTH M1 AND
C                               STARTS AT TCOS(I1+1).  THE SECOND
C                               STRING IS OF LENGTH M2 AND STARTS AT
C                               TCOS(I2+1).  THE MERGED STRING GOES INTO
C                               TCOS(I3+1).  IT IS ALMOST A STRAIGHT
C                               COPY FROM THE ORIGINAL FISHPAK VERSION.

      j1 = 1
      j2 = 1
      j = i3

      IF (m1+m2.EQ.0) GO TO 180

      IF (m1.EQ.0) GO TO 160

      IF (m2.EQ.0) GO TO 130
  100 CONTINUE
      j = j + 1
      x = tcos(i1+j1)
      y = tcos(i2+j2)

      IF (x.LE.y) THEN
          tcos(j) = x
          j1 = j1 + 1

          IF (j1.GT.m1) GO TO 150

          GO TO 100

      ELSE
          tcos(j) = y
          j2 = j2 + 1
      END IF

      IF (j2.LE.m2) GO TO 100

      IF (j1.GT.m1) GO TO 180
  130 CONTINUE
      k = j - j1 + 1

      DO 140 j = j1,m1
          tcos(k+j) = tcos(j+i1)
  140 CONTINUE

      GO TO 180

  150 CONTINUE

      IF (j2.GT.m2) GO TO 180
  160 CONTINUE
      k = j - j2 + 1

      DO 170 j = j2,m2
          tcos(k+j) = tcos(j+i2)
  170 CONTINUE
  180 CONTINUE
      RETURN
      END SUBROUTINE fsh13s
C
C
      SUBROUTINE fsh14s(n,d,e2,ierr)
C
      REAL(EB)   d(n),e2(n)
C
C                               THIS SUBROUTINE IS A MODIFICATION OF THE
C                               EISPACK SUBROUTINE TQLRAT ALGORITHM 464,
C                               COMM. ACM 16, 689(1973) BY REINSCH.

C THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.

C                               ON INPUT-
C
C                               N IS THE ORDER OF THE MATRIX,
C
C D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C E2 CONTAINS THE                SUBDIAGONAL ELEMENTS OF THE
C INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C                               ON OUTPUT-
C
C D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C                               THE SMALLEST EIGENVALUES,
C
C                               E2 HAS BEEN DESTROYED,
C
C                               IERR IS SET TO
C                               ZERO       FOR NORMAL RETURN,
C J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                               DETERMINED AFTER 30 ITERATIONS.
C
C ********** EPS IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C THE RELATIVE PREC OF FLOATING POINT ARITHMETIC.

      ierr = 0

      IF (n.EQ.1) GO TO 240

      eps = fsh20s()
      eps = 3.55e-15

      DO 100 i = 2,n
          e2(i-1) = e2(i)*e2(i)
  100 CONTINUE

      f = 0.0
      b = 0.0
      e2(n) = 0.0

      DO 210 l = 1,n
          j = 0
          h = eps* (abs(d(l))+sqrt(e2(l)))

          IF (b.GT.h) GO TO 110
          b = h
          c = b*b

C            LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT

  110     CONTINUE
          DO 120 m = l,n

              IF (e2(m).LE.c) GO TO 130

C            E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C            THROUGH THE BOTTOM OF THE LOOP

  120     CONTINUE

  130     CONTINUE
          IF (m.EQ.l) GO TO 170

  140     CONTINUE
          IF (j.EQ.30) GO TO 230
          j = j + 1

C                               FORM SHIFT

          l1 = l + 1
          s = sqrt(e2(l))
          g = d(l)
          p = (d(l1)-g)/ (2.0*s)
          r = sqrt(p*p+1.0_EB)
          d(l) = s/ (p+sign(r,p))
          h = g - d(l)

          DO 150 i = l1,n
              d(i) = d(i) - h
  150     CONTINUE

          f = f + h

C                               RATIONAL QL TRANSFORMATION

          g = d(m)

          IF (g.EQ.0.0) g = b
          h = g
          s = 0.0
          mml = m - l

C                               FOR I=M-1 STEP -1 UNTIL L DO --

          DO 160 ii = 1,mml
              i = m - ii
              p = g*h
              r = p + e2(i)
              e2(i+1) = s*r
              s = e2(i)/r
              d(i+1) = h + s* (h+d(i))
              g = d(i) - e2(i)/g

              IF (g.EQ.0.0) g = b
              h = g*p/r
  160     CONTINUE
C
          e2(l) = s*g
          d(l) = h

C                               GUARD AGAINST UNDERFLOWED H

          IF (h.EQ.0.0) GO TO 170

          IF (abs(e2(l)).LE.abs(c/h)) GO TO 170
          e2(l) = h*e2(l)

          IF (e2(l).NE.0.0) GO TO 140
  170     CONTINUE
          p = d(l) + f

C                               ORDER EIGENVALUES

          IF (l.EQ.1) GO TO 190

C                               FOR I=L STEP -1 UNTIL 2 DO --

          DO 180 ii = 2,l
              i = l + 2 - ii

              IF (p.GE.d(i-1)) GO TO 200
              d(i) = d(i-1)
  180     CONTINUE

  190     CONTINUE
          i = 1
  200     CONTINUE
          d(i) = p
  210 CONTINUE

      IF (abs(d(n)).GE.abs(d(1))) GO TO 240
      nhalf = n/2
      DO 220 i = 1,nhalf
          ntop = n - i
          dhold = d(i)
          d(i) = d(ntop+1)
          d(ntop+1) = dhold
  220 CONTINUE

      GO TO 240

C                               SET ERROR -- NO CONVERGENCE TO AN
C                               EIGENVALUE AFTER 30 ITERATIONS

  230 CONTINUE
      ierr = l
  240 CONTINUE
      RETURN
      END SUBROUTINE fsh14s
C
C
      SUBROUTINE h2czis(xs,xf,l,lbdcnd,ys,yf,m,mbdcnd,h,
     *                  elmbda,ldimf,ierror,save)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB) h(0:*)
      REAL(EB) save(-3:*)



C                               CHECK FOR INVALID INPUT

      ierror = 0

      IF (xs.GT.xf) THEN
          ierror = ierror + 1
          save(ierror) = 1.
      END IF

      IF (l.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 2.
      END IF

      IF (lbdcnd.LT.0 .OR. lbdcnd.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 3.
      END IF

      IF (ys.GT.yf) THEN
          ierror = ierror + 1
          save(ierror) = 4.
      END IF

      IF (m.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 5.
      END IF

      IF (mbdcnd.LT.0 .OR. mbdcnd.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 6.
      END IF

      IF (ldimf.LT.l) THEN
          ierror = ierror + 1
          save(ierror) = 7.
      END IF

      IF (ierror.NE.0) THEN
          RETURN
      ELSE
          save(1) = 0.
      END IF

C                               DEFINE GRID PARAMETERS

      dx = (xf-xs)/l
      dy = (yf-ys)/m
      dlxsqr = 1./dx**2
      dlysqr = 1./dy**2
      lp = lbdcnd + 1
      mp = mbdcnd + 1

C                               ALLOCATE SAVE ARRAY

      ia = 9
      ib = ia + l
      ic = ib + l
      id = ic + l
      is = id + l

C                               DEFINE THE A,B,C,D COEFFICIENTS
C                               IN SAVE ARRAY
      DO 100 i = 0,l - 1
          hm = .5*(h(i)+h(i+1))
          hp = .5*(h(i+1)+h(i+2))
          save(ia+i) = dlxsqr/(h(i+1)*hm)
          save(ic+i) = dlxsqr/(h(i+1)*hp)
          save(ib+i) = - (save(ia+i)+save(ic+i)) + elmbda
          save(id+i) = dlysqr
  100 CONTINUE

      select case(lp)
      case(2)
      save(ib) = save(ib) - save(ia)
      save(ic-1) = save(ic-1) - save(id-1)
      case(3)
      save(ib) = save(ib) - save(ia)
      save(ic-1) = save(ic-1) + save(id-1)
      case(4)
      save(ib) = save(ib) + save(ia)
      save(ic-1) = save(ic-1) + save(id-1)
      case(5)
      save(ib) = save(ib) + save(ia)
      save(ic-1) = save(ic-1) - save(id-1)
      end select

C                               DETERMINE WHETHER OR NOT BOUNDARY
C                               CONDITION IS PERIODIC IN X

      IF (lbdcnd.EQ.0) THEN
          lperod = 0
      ELSE
          lperod = 1
      END IF

C                               INITIALIZE SOLVER ROUTINE S2CFIS.
C
      CALL s2cfis(lperod,l,mbdcnd,m,save(ia),save(ib),
     *            save(ic),save(id),ldimf,ir,save(is))

C                               TEST ERROR FLAG FROM S2CFIS FOR
C                               INTERNAL ERROR

      IF (ir.NE.0) THEN
          save(1) = 99.
          ierror = 1

          RETURN
      END IF

C                               SAVE PARAMETERS FOR H2CCSS IN SAVE ARRAY

      save(2) = dx
      save(3) = l
      save(4) = lp
      save(5) = dy
      save(6) = m
      save(7) = mp
      save(8) = elmbda

      RETURN
      END SUBROUTINE h2czis
C
C
      SUBROUTINE h2czss(bdxs,bdxf,bdys,bdyf,ldimf,f,pertrb,save,w,h)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB) bdxs(*),bdxf(*),bdys(*),bdyf(*),f(ldimf,*),save(-3:*),
     .         w(*),h(*)


C                               CHECK VALUE OF IERROR (=SAVE(1)).
C                               IF NON-ZERO, RETURN.

      IF (save(1).NE.0.) RETURN

C                               GET PARAMETERS FOR H2CCSS FROM SAVE
C                               ARRAY WHERE THEY WERE STORED IN
C                               INITIALIZATION SUBROUTINE H2CCIS.

      dx = save(2)
      l = save(3)
      lp = save(4)
      dy = save(5)
      m = save(6)
      mp = save(7)
      elmbda = save(8)

      dlxrcp = 1./dx
      twdxsq = 2./dx**2
      dlyrcp = 1./dy
      twdysq = 2./dy**2

C                               ALLOCATE SAVE ARRAY

      ia = 9
      ib = ia + l
      ic = ib + l
      id = ic + l
      is = 9 + 4*l


C                               ENTER BOUNDARY DATA FOR X-BOUNDARIES

      if (lp.eq.2 .or. lp.eq.3) then
      DO j = 1,m
      f(1,j) = f(1,j) - 2.*bdxs(j)*save(ia)
      enddo
      endif

      if (lp.eq.4 .or. lp.eq.5) then
      DO j = 1,m
      f(1,j) = f(1,j) + save(ia)*dx*bdxs(j)
      enddo
      endif

      if (lp.eq.2 .or. lp.eq.5) then
      DO j = 1,m
      f(l,j) = f(l,j) - 2.*bdxf(j)*save(id-1)
      enddo
      endif

      if (lp.eq.3 .or. lp.eq.4) then
      DO j = 1,m
      f(l,j) = f(l,j) - save(id-1)*dx*bdxf(j)
      enddo
      endif

C                               ENTER BOUNDARY DATA FOR Y-BOUNDARIES

      if (mp.eq.2 .or. mp.eq.3) then
      DO i = 1,l
      f(i,1) = f(i,1) - bdys(i)*twdysq
      enddo
      endif

      if (mp.eq.4 .or. mp.eq.5) then
      DO i = 1,l
      f(i,1) = f(i,1) + bdys(i)*dlyrcp
      enddo
      endif

      if (mp.eq.2 .or. mp.eq.5) then
      DO i = 1,l
      f(i,m) = f(i,m) - bdyf(i)*twdysq
      enddo
      endif

      if (mp.eq.3 .or. mp.eq.4) then
      DO i = 1,l
      f(i,m) = f(i,m) - bdyf(i)*dlyrcp
      enddo
      endif

      pert=0.
      pertrb = 0.
      ising = 0

C                               FOR SINGULAR PROBLEMS ADJUST DATA TO
C                               INSURE A SOLUTION WILL EXIST.  GO THRU
C                               THIS CODE TWICE: ISING=1 FOR CALCULATING
C                               PERTRB; ISING=2 FOR NORMALIZING SOLUTION
C                               AFTER IT IS COMPUTED.

      SELECT CASE(lp)
      CASE(1)   ; GOTO 630
      CASE(2:3) ; GOTO 750
      CASE(4)   ; GOTO 630
      CASE(5)   ; GOTO 750
      END SELECT

  630 CONTINUE
      SELECT CASE(mp)
      CASE(1)   ; GOTO 640
      CASE(2:3) ; GOTO 750
      CASE(4)   ; GOTO 640
      CASE(5)   ; GOTO 750
      END SELECT

  640 CONTINUE

      IF (elmbda.NE.0.) GO TO 750
      ising = 1
  660 CONTINUE
      pert = 0.
      DO 670 i = 1,l
          w(i) = 0.
  670 CONTINUE
          DO 690 j = 1,m
              DO 680 i = 1,l
                  w(i) = w(i) + f(i,j)
  680         CONTINUE
  690     CONTINUE
      s1 = 0.
      s3 = 0.
      DO 710 i = 1,l
          s3 = s3 + h(i+1)
          s1 = s1 + h(i+1)*w(i)
  710 CONTINUE

      s3   = s3*m
      pert = s1/s3

C                               ADJUST F ARRAY BY PERT

          DO 730 j = 1,m
              DO 720 i = 1,l
                  f(i,j) = f(i,j) - pert
  720         CONTINUE
  730     CONTINUE
  750 CONTINUE

C                               IF NORMALIZING SOLUTION, RESTORE PERTRB
C                               AND JUMP TO END

      IF (ising.EQ.2) THEN
          pertrb = prtsav

          GO TO 800

      END IF

      prtsav = pert

C                               SOLVE THE EQUATION

      CALL s2cfss(ldimf,f,save(is),w)

C                               IF A SINGULAR PROBLEM,
C                               RE-NORMALIZE SOLUTION (ISING=2)

      IF (ising.EQ.1) THEN
          ising = 2

          GO TO 660

      END IF

  800 CONTINUE
      RETURN
      END SUBROUTINE h2czss
C
C
      SUBROUTINE s2cfis(lperod,l,mperod,m,a,b,c,d,ldimf,ierror,save)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB) a(l),b(l),c(l),d(l),save(-3:*)


C                               CHECK FOR INVALID INPUT

      ierror = 0

      IF (lperod.NE.0 .AND. lperod.NE.1) THEN
          ierror = 1
          save(1) = 1.
      END IF

      IF (l.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 2.
      END IF

      IF (mperod.LT.0 .AND. mperod.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 3.
      END IF

      IF (m.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 4.
      END IF

      IF (ldimf.LT.l) THEN
          ierror = ierror + 1
          save(ierror) = 5.
      END IF

      IF (lperod.EQ.0) THEN
          DO 100 i = 1,l

              IF (a(i).NE.a(1)) GO TO 110
              IF (b(i).NE.b(1)) GO TO 110
              IF (c(i).NE.a(1)) GO TO 110
              IF (d(i).NE.d(1)) GO TO 110

  100     CONTINUE

          GO TO 120

  110     CONTINUE
          ierror = ierror + 1
          save(ierror) = 6.
      END IF

  120 CONTINUE

      IF (ierror.NE.0) THEN
          RETURN
      ELSE
          save(1) = ierror
      END IF

      ldimfc=l
      IF (ldimf.GT.l .AND. mod(l,2).EQ.0) ldimfc=l+1
      igrid = 2
      n = 1
      nperod = 0
      scal = 0.

      CALL fsh00s(igrid,lperod,l,mperod,m,nperod,n,ldimfc,
     *            scal,a,b,c,d,save)

      RETURN
      END SUBROUTINE s2cfis
C
C
      SUBROUTINE s2cfss(ldimf,f,save,w)

C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      REAL(EB) f(ldimf,*),save(-3:*),w(*)


C                               CHECK VALUE OF IERROR (=SAVE(1)).
C                               IF NON-ZERO, RETURN.

      IF (save(1).NE.0.) RETURN

C                               DEBUG PRINTOUTS UNDER CONTROL OF
C                               PARAMETER LVLPRN

      m = save(4)

      CALL fsh02s(ldimf,m,f,save,w)

      RETURN
      END SUBROUTINE s2cfss
C
C
      SUBROUTINE h2cyis(rs,rf,l,lbdcnd,zs,zf,n,nbdcnd,
     *                  elmbda,xmu,ldimf,ierror,save)
C
C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      DIMENSION save(-3:*)


C                               CHECK FOR INVALID INPUT

      ierror = 0

      IF (rs.LT.0.0) THEN
          ierror = ierror + 1
          save(ierror) = 1.
      END IF

      IF (rf.LE.rs) THEN
          ierror = ierror + 1
          save(ierror) = 2.
      END IF

      IF (l.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 3.
      END IF

      IF (lbdcnd.LT.1 .OR. lbdcnd.GT.6) THEN
          ierror = ierror + 1
          save(ierror) = 4.
      END IF

      IF (rs.EQ.0.0 .AND. lbdcnd.LE.4) THEN
          ierror = ierror + 1
          save(ierror) = 5.
      END IF

      IF (rs.EQ.0.0 .AND. xmu.NE.0.) THEN
          ierror = ierror + 1
          save(ierror) = 6.
      END IF

      IF (zf.LE.zs) THEN
          ierror = ierror + 1
          save(ierror) = 7.
      END IF

      IF (n.LT.3) THEN
          ierror = ierror + 1
          save(ierror) = 8.
      END IF

      IF (nbdcnd.LT.0 .OR. nbdcnd.GT.4) THEN
          ierror = ierror + 1
          save(ierror) = 9.
      END IF

      IF (ldimf.LT.l) THEN
          ierror = ierror + 1
          save(ierror) = 10.
      END IF

      IF (ierror.NE.0) THEN
          RETURN
      ELSE
          save(1) = ierror
      END IF

C                               DEFINE GRID PARAMETERS

      dz = (zf-zs)/float(n)
      dzsqr = 1./ (dz**2)
      np = nbdcnd + 1

      dr = (rf-rs)/float(l)
      drby2 = dr/2.
      drsqr = 1./ (dr**2)

C                               ALLOCATE SAVE ARRAY

      ia = 13
      ib = ia + l
      ic = ib + l
      id = ic + l
      is = id + l
      isvps = is + l
C
C                               DEFINE A,B,C,D,S ARRAYS IN ARRAY SAVE.
C                               S ARRAY IS SYMMETRIZER FOR SING. CASES

CDIR$ IVDEP
      DO 100 i = 1,l
          ri = rs + (i-.5)*dr
          save(id+i) = dzsqr
          save(ia+i) = drsqr* (ri-drby2)/ri
          save(ic+i) = drsqr* (ri+drby2)/ri
          save(ib+i) = - (save(ia+i)+save(ic+i)) + elmbda + xmu/ri**2
          save(is+i) = ri
  100 CONTINUE

C                               DEFINE BOUNDARY COEFFICIENTS

      SELECT CASE(lbdcnd)
      CASE(1:2) ; GOTO 110
      CASE(3:6) ; GOTO 120
      END SELECT

  110 CONTINUE
      save(ib+1) = save(ib+1) - save(ia+1)

      GO TO 130

  120 CONTINUE
      save(ib+1) = save(ib+1) + save(ia+1)
  130 CONTINUE

      SELECT CASE(lbdcnd)
      CASE(1)   ; GOTO 140
      CASE(2:3) ; GOTO 150
      CASE(4:5) ; GOTO 140
      CASE(6)   ; GOTO 150
      END SELECT

  140 CONTINUE
      save(ib+l) = save(ib+l) - save(ic+l)

      GO TO 160

  150 CONTINUE
      save(ib+l) = save(ib+l) + save(ic+l)
  160 CONTINUE

C                               INITIALIZE SOLVER ROUTINE S2CFIS

      CALL s2cfis(1,l,nbdcnd,n,save(ia+1),save(ib+1),
     *            save(ic+1),save(id+1),ldimf,ir,save(isvps+1))

C                               TEST ERROR FLAG FROM S2CFIS FOR
C                               INTERNAL ERROR

      IF (ir.NE.0) THEN
          save(1) = 99.
          ierror = 1

          RETURN
      END IF

C                               SAVE PARAMETERS FOR H3CYSS IN SAVE ARRAY

      save(2) = dr
      save(3) = l
      save(4) = lbdcnd
C     save(5) = drsqr           this parameter is not used
      save(6) = dz
      save(7) = n
      save(8) = np
      save(9) = dzsqr
      save(10) = 1./dz
      save(11) = elmbda
      save(12) = xmu

      RETURN
      END SUBROUTINE h2cyis
C
C
      SUBROUTINE h2cyss(bdrs,bdrf,bdzs,bdzf,ldimf,f,pertrb,save,w)
C
C +--------------------------------------------------------------------+
C |                                                                    |
C |                       COPYRIGHT (C) 1989 BY                        |
C |                          ROLAND A. SWEET                           |
C |                        ALL RIGHTS RESERVED                         |
C |                                                                    |
C +--------------------------------------------------------------------+

      DIMENSION bdrs(*),bdrf(*),bdzs(*),bdzf(*),f(ldimf,*),save(-3:*),
     .          w(*)


C                               CHECK VALUE OF IERROR (=SAVE(1)).
C                               IF NON-ZERO, RETURN.

      IF (save(1).NE.0.) RETURN

C                               GET PARAMETERS FOR H2CYSS FROM SAVE
C                               ARRAYWHERE THEY WERE STORED IN
C                               INITIALIZATION SUBROUTINE H2CYIS.

      dr = save(2)
      l = save(3)
      lbdcnd = save(4)
C     drsqr = save(5)      these two parameters are not used
C     dz = save(6)

      n = save(7)
      np = save(8)
      dzsqr = save(9)
      dzr = save(10)

      elmbda = save(11)
      xmu = save(12)

C                               ALLOCATE SAVE ARRAY

      ia = 13
      ib = ia + l
      ic = ib + l
      id = ic + l
      is = id + l
      isvps = is + l

C                               ENTER BOUNDARY DATA FOR R-BOUNDARIES.

      if (lbdcnd.eq.1 .or. lbdcnd.eq.2) then
      DO k = 1,n
      f(1,k) = f(1,k) - 2.*save(ia+1)*bdrs(k)
      enddo
      endif

      if (lbdcnd.eq.3 .or. lbdcnd.eq.4) then
      DO k = 1,n
      f(1,k) = f(1,k) + dr*save(ia+1)*bdrs(k)
      enddo
      endif

      if (lbdcnd.eq.4 .or. lbdcnd.eq.5) then
      DO k = 1,n
      f(l,k) = f(l,k) - 2.*save(ic+l)*bdrf(k)
      enddo
      endif

      if (lbdcnd.eq.2 .or. lbdcnd.eq.3 .or. lbdcnd.eq.6) then
      DO k = 1,n
      f(l,k) = f(l,k) - dr*save(ic+l)*bdrf(k)
      enddo
      endif

C                               ENTER BOUNDARY DATA FOR Z-BOUNDARIES.

      if (np.eq.2 .or. np.eq.3) then
      DO i = 1,l
      f(i,1) = f(i,1) - 2.*dzsqr*bdzs(i)
      enddo
      endif

      if (np.eq.4 .or. np.eq.5) then
      DO i = 1,l
      f(i,1) = f(i,1) + dzr*bdzs(i)
      enddo
      endif

      if (np.eq.2 .or. np.eq.5) then
      DO i = 1,l
      f(i,n) = f(i,n) - 2.*dzsqr*bdzf(i)
      enddo
      endif

      if (np.eq.3 .or. np.eq.4) then
      DO i = 1,l
      f(i,n) = f(i,n) - dzr*bdzf(i)
      enddo
      endif

      pertrb = 0.
      ising = 0

C                               FOR SINGULAR PROBLEMS ADJUST DATA TO
C                               INSURE A SOLUTION WILL EXIST.  GO THRU
C                               THIS CODE TWICE: ISING=1 FOR CALCULATING
C                               PERTRB; ISING=2 FOR NORMALIZING SOLUTION
C                               AFTER IT IS COMPUTED.

      SELECT CASE(lbdcnd)
      CASE(1:2) ; GOTO 770
      CASE(3)   ; GOTO 620
      CASE(4:5) ; GOTO 770
      CASE(6)   ; GOTO 620
      END SELECT

  620 CONTINUE
      SELECT CASE(np)
      CASE(1)   ; GOTO 640
      CASE(2:3) ; GOTO 770
      CASE(4)   ; GOTO 640
      CASE(5)   ; GOTO 770
      END SELECT

  640 CONTINUE
      IF (elmbda.eq.0.) THEN
         GOTO 650
         ELSE
         GOTO 770
         endif

  650 CONTINUE
      IF (xmu.eq.0.) THEN
         GOTO 665
         ELSE
         GOTO 770
         endif
  665 CONTINUE
      ising = 1
  660 CONTINUE
      w(1:l) = 0.
      DO k = 1,n
      DO i = 1,l
          w(i) = w(i) + f(i,k)
          enddo
          enddo

      s3 = 0.
      s1 = 0.
      DO i = 1,l
          s3 = s3 + save(is+i)
          s1 = s1 + w(i)*save(is+i)
          enddo
      s3 = s3*n
      pert = s1/s3
C
      DO k = 1,n
      DO i = 1,l
         f(i,k) = f(i,k) - pert
         enddo
         enddo
C
C                               IF NORMALIZING SOLUTION, RESTORE PERTRB
C                               AND JUMP TO END

      IF (ising.EQ.2) THEN
          pertrb = prtsav
          GO TO 800
          ENDIF

      prtsav = pert

  770 CONTINUE

C                               SOLVE SYSTEM USING S2CFSS

      CALL s2cfss(ldimf,f,save(isvps+1),w)

C                               IF A SINGULAR PROBLEM,
C                               RE-NORMALIZE SOLUTION (ISING=2)

      IF (ising.EQ.1) THEN
          ising = 2

          GO TO 660

      END IF

  800 CONTINUE
      RETURN
      END SUBROUTINE h2cyss
C
      END MODULE POIS
