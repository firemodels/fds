
!*********************************************************
!*****   FDS5-Evac: Modules and misc routines        *****
!*********************************************************
!
! VTT Technical Research Centre of Finland            2008
!
! Author:  Timo Korhonen
! Date:    13.8.2008
! FDS Version:  5.2.0
! Evac Version: 2.0.0
!
! This file (ieva.f90) contains:
! * dcdflib.f90 (Netlib cumulative density function library)
! * stat.f90 (random number generation etc, by T.K./Pfs3.0)
! * miscellaneous Evac routines, like initializations
!
!*********************************************************
Module DCDFLIB
  Implicit None

  Private
  Public cdfbet,cdfgam,cdfnor,gamma_ieva

  !
  !*********************************************************
  !*****             DCDFLIB                           *****
  !*********************************************************
  !
  ! The Fortran90 version form:
  ! http://people.scs.fsu.edu/~burkardt/f_src/dcdflib/dcdflib.html
  !
  ! Unused functions commented out by Timo Korhonen, 2017 (not all dcdflib needed in stat.f90 routines)
  ! Implicit integer = real conversions made explicit, integer = INT(real) by Timo Korhonen, 2017
  ! Unused 'label continue' statements removed by Timo Korhonen, 2017
  ! Assigned goto statements removed by Timo Korhonen, 2008.
  ! Statement functions replaced with internal function by Timo Korhonen, 2008.
  !
  ! Original F77 version from Netlib: http://www.netlib.org/random:
  !
  !                                    DCDFLIB
  !
  !             Library of Fortran Routines for Cumulative Distribution
  !                  Functions, Inverses, and Other Parameters
  !
  !                                 (February, 1994)
  !                     Summary Documentation of Each Routine
  !
  !                             Compiled and Written by:
  !
  !                                  Barry W. Brown
  !                                   James Lovato
  !                                   Kathy Russell
  !
  !     Department of Biomathematics, Box 237
  !     The University of Texas, M.D. Anderson Cancer Center
  !     1515 Holcombe Boulevard
  !     Houston, TX      77030
  !
  !  This work was supported by grant CA-16672 from the National Cancer Institute.
  !
  !
  !                           SUMMARY OF DCDFLIB
  !
  ! This  library  contains routines  to compute  cumulative  distribution
  ! functions, inverses, and    parameters  of the  distribution  for  the
  ! following set of statistical distributions:
  !
  !     (1) Beta
  !     (2) Binomial
  !     (3) Chi-square
  !     (4) Noncentral Chi-square
  !     (5) F
  !     (6) Noncentral F
  !     (7) Gamma
  !     (8) Negative Binomial
  !     (9) Normal
  !     (10) Poisson
  !     (11) Student's t
  !
  ! Given values of all but one parameter of a distribution, the other is
  ! computed.  These calculations are done with  FORTRAN Double Precision
  ! variables.
  !*********************************************************
  !
  !
  !*********************************************************

!!$  Real ( kind = 8 ) :: algdiv, alnrel, apser, bcorr, beta, beta_asym, &
!!$       beta_frac, beta_log, beta_pser, beta_rcomp, beta_rcomp1, &
!!$       beta_up, dbetrm, dexpm1, dinvnr, dlanor, dstrem, dt1, &
!!$       error_f, error_fc, esum, eval_pol, exparg, fpser, &
!!$       gam1, gamma, gamma_ln1, gamma_log, gsumln, ipmpar, &
!!$       psi, rcomp, rexp, rlog, rlog1, stvaln

Contains

  Function algdiv ( a, b )

    !*****************************************************************************80
    !
    !! ALGDIV computes ln ( Gamma ( B ) / Gamma ( A + B ) ) when 8 <= B.
    !
    !  Discussion:
    !
    !    In this algorithm, DEL(X) is the function defined by
    !
    !      ln ( Gamma(X) ) = ( X - 0.5 ) * ln ( X ) - X + 0.5 * ln ( 2 * PI )
    !                      + DEL(X).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, define the arguments.
    !
    !    Output, real ( kind = 8 ) ALGDIV, the value of ln(Gamma(B)/Gamma(A+B)).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) algdiv
    !    Real ( kind = 8 ) alnrel
    Real ( kind = 8 ) b
    Real ( kind = 8 ) c
    Real ( kind = 8 ), Parameter :: c0 =  0.833333333333333D-01
    Real ( kind = 8 ), Parameter :: c1 = -0.277777777760991D-02
    Real ( kind = 8 ), Parameter :: c2 =  0.793650666825390D-03
    Real ( kind = 8 ), Parameter :: c3 = -0.595202931351870D-03
    Real ( kind = 8 ), Parameter :: c4 =  0.837308034031215D-03
    Real ( kind = 8 ), Parameter :: c5 = -0.165322962780713D-02
    Real ( kind = 8 ) d
    Real ( kind = 8 ) h
    Real ( kind = 8 ) s11
    Real ( kind = 8 ) s3
    Real ( kind = 8 ) s5
    Real ( kind = 8 ) s7
    Real ( kind = 8 ) s9
    Real ( kind = 8 ) t
    Real ( kind = 8 ) u
    Real ( kind = 8 ) v
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    Real ( kind = 8 ) x2

    If ( b < a ) Then
       h = b / a
       c = 1.0D+00 / ( 1.0D+00 + h )
       x = h / ( 1.0D+00 + h )
       d = a + ( b - 0.5D+00 )
    Else
       h = a / b
       c = h / ( 1.0D+00 + h )
       x = 1.0D+00 / ( 1.0D+00 + h )
       d = b + ( a - 0.5D+00 )
    End If
    !
    !  Set SN = (1 - X**N)/(1 - X).
    !
    x2 = x * x
    s3 = 1.0D+00 + ( x + x2 )
    s5 = 1.0D+00 + ( x + x2 * s3 )
    s7 = 1.0D+00 + ( x + x2 * s5 )
    s9 = 1.0D+00 + ( x + x2 * s7 )
    s11 = 1.0D+00 + ( x + x2 * s9 )
    !
    !  Set W = DEL(B) - DEL(A + B).
    !
    t = ( 1.0D+00 / b )**2
    w = (((( &
         c5 * s11  * t &
         + c4 * s9 ) * t &
         + c3 * s7 ) * t &
         + c2 * s5 ) * t &
         + c1 * s3 ) * t &
         + c0

    w = w * ( c / b )
    !
    !  Combine the results.
    !
    u = d * alnrel ( a / b )
    v = a * ( Log ( b ) - 1.0D+00 )

    If ( v < u ) Then
       algdiv = ( w - v ) - u
    Else
       algdiv = ( w - u ) - v
    End If

    Return
  End Function algdiv

  Function alnrel ( a )

    !*****************************************************************************80
    !
    !! ALNREL evaluates the function ln ( 1 + A ).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, the argument.
    !
    !    Output, real ( kind = 8 ) ALNREL, the value of ln ( 1 + A ).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) alnrel
    Real ( kind = 8 ), Parameter :: p1 = -0.129418923021993D+01
    Real ( kind = 8 ), Parameter :: p2 =  0.405303492862024D+00
    Real ( kind = 8 ), Parameter :: p3 = -0.178874546012214D-01
    Real ( kind = 8 ), Parameter :: q1 = -0.162752256355323D+01
    Real ( kind = 8 ), Parameter :: q2 =  0.747811014037616D+00
    Real ( kind = 8 ), Parameter :: q3 = -0.845104217945565D-01
    Real ( kind = 8 ) t
    Real ( kind = 8 ) t2
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x

    If ( Abs ( a ) <= 0.375D+00 ) Then

       t = a / ( a +  2.0D+00  )
       t2 = t * t

       w = ((( p3 * t2 + p2 ) * t2 + p1 ) * t2 + 1.0D+00 ) &
            / ((( q3 * t2 + q2 ) * t2 + q1 ) * t2 + 1.0D+00 )

       alnrel =  2.0D+00  * t * w

    Else

       x = 1.0D+00 + Real ( a, kind = 8 )
       alnrel = Log ( x )

    End If

    Return
  End Function alnrel

  Function apser ( a, b, x, eps )

    !*****************************************************************************80
    !
    !! APSER computes the incomplete beta ratio I(SUB(1-X))(B,A).
    !
    !  Discussion:
    !
    !    APSER is used only for cases where
    !
    !      A <= min ( EPS, EPS * B ),
    !      B * X <= 1, and
    !      X <= 0.5.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, X, the parameters of the
    !    incomplete beta ratio.
    !
    !    Input, real ( kind = 8 ) EPS, a tolerance.
    !
    !    Output, real ( kind = 8 ) APSER, the computed value of the
    !    incomplete beta ratio.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) aj
    Real ( kind = 8 ) apser
    Real ( kind = 8 ) b
    Real ( kind = 8 ) bx
    Real ( kind = 8 ) c
    Real ( kind = 8 ) eps
    Real ( kind = 8 ), Parameter :: g = 0.577215664901533D+00
    Real ( kind = 8 ) j
    !    Real ( kind = 8 ) psi
    Real ( kind = 8 ) s
    Real ( kind = 8 ) t
    Real ( kind = 8 ) tol
    Real ( kind = 8 ) x

    bx = b * x
    t = x - bx

    If ( b * eps <= 0.02D+00 ) Then
       c = Log ( x ) + psi ( b ) + g + t
    Else
       c = Log ( bx ) + g + t
    End If

    tol = 5.0D+00 * eps * Abs ( c )
    j = 1.0D+00
    s = 0.0D+00

    Do

       j = j + 1.0D+00
       t = t * ( x - bx / j )
       aj = t / j
       s = s + aj

       If ( Abs ( aj ) <= tol ) Then
          Exit
       End If

    End Do

    apser = -a * ( c + s )

    Return
  End Function apser

  Function bcorr ( a0, b0 )

    !*****************************************************************************80
    !
    !! BCORR evaluates DEL(A0) + DEL(B0) - DEL(A0 + B0).
    !
    !  Discussion:
    !
    !    The function DEL(A) is a remainder term that is used in the expression:
    !
    !      ln ( Gamma ( A ) ) = ( A - 0.5 ) * ln ( A )
    !        - A + 0.5 * ln ( 2 * PI ) + DEL ( A ),
    !
    !    or, in other words, DEL ( A ) is defined as:
    !
    !      DEL ( A ) = ln ( Gamma ( A ) ) - ( A - 0.5 ) * ln ( A )
    !        + A + 0.5 * ln ( 2 * PI ).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A0, B0, the arguments.
    !    It is assumed that 8 <= A0 and 8 <= B0.
    !
    !    Output, real ( kind = 8 ) BCORR, the value of the function.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) a0
    Real ( kind = 8 ) b
    Real ( kind = 8 ) b0
    Real ( kind = 8 ) bcorr
    Real ( kind = 8 ) c
    Real ( kind = 8 ), Parameter :: c0 =  0.833333333333333D-01
    Real ( kind = 8 ), Parameter :: c1 = -0.277777777760991D-02
    Real ( kind = 8 ), Parameter :: c2 =  0.793650666825390D-03
    Real ( kind = 8 ), Parameter :: c3 = -0.595202931351870D-03
    Real ( kind = 8 ), Parameter :: c4 =  0.837308034031215D-03
    Real ( kind = 8 ), Parameter :: c5 = -0.165322962780713D-02
    Real ( kind = 8 ) h
    Real ( kind = 8 ) s11
    Real ( kind = 8 ) s3
    Real ( kind = 8 ) s5
    Real ( kind = 8 ) s7
    Real ( kind = 8 ) s9
    Real ( kind = 8 ) t
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    Real ( kind = 8 ) x2

    a = Min ( a0, b0 )
    b = Max ( a0, b0 )

    h = a / b
    c = h / ( 1.0D+00 + h )
    x = 1.0D+00 / ( 1.0D+00 + h )
    x2 = x * x
    !
    !  Set SN = (1 - X**N)/(1 - X)
    !
    s3 = 1.0D+00 + ( x + x2 )
    s5 = 1.0D+00 + ( x + x2 * s3 )
    s7 = 1.0D+00 + ( x + x2 * s5 )
    s9 = 1.0D+00 + ( x + x2 * s7 )
    s11 = 1.0D+00 + ( x + x2 * s9 )
    !
    !  Set W = DEL(B) - DEL(A + B)
    !
    t = ( 1.0D+00 / b )**2

    w = (((( &
         c5 * s11  * t &
         + c4 * s9 ) * t &
         + c3 * s7 ) * t &
         + c2 * s5 ) * t &
         + c1 * s3 ) * t &
         + c0

    w = w * ( c / b )
    !
    !  Compute  DEL(A) + W.
    !
    t = ( 1.0D+00 / a )**2

    bcorr = ((((( &
         c5   * t &
         + c4 ) * t &
         + c3 ) * t &
         + c2 ) * t &
         + c1 ) * t &
         + c0 ) / a + w

    Return
  End Function bcorr

  Function beta ( a, b )

    !*****************************************************************************80
    !
    !! BETA evaluates the beta function.
    !
    !  Modified:
    !
    !    03 December 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, the arguments of the beta function.
    !
    !    Output, real ( kind = 8 ) BETA, the value of the beta function.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) b
    Real ( kind = 8 ) beta
    !    Real ( kind = 8 ) beta_log

    beta = Exp ( beta_log ( a, b ) )

    Return
  End Function beta

  Function beta_asym ( a, b, lambda, eps )

    !*****************************************************************************80
    !
    !! BETA_ASYM computes an asymptotic expansion for IX(A,B), for large A and B.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, the parameters of the function.
    !    A and B should be nonnegative.  It is assumed that both A and B
    !    are greater than or equal to 15.
    !
    !    Input, real ( kind = 8 ) LAMBDA, the value of ( A + B ) * Y - B.
    !    It is assumed that 0 <= LAMBDA.
    !
    !    Input, real ( kind = 8 ) EPS, the tolerance.
    !
    Implicit None

    Integer, Parameter :: num = 20

    Real ( kind = 8 ) a
    Real ( kind = 8 ) a0(num+1)
    Real ( kind = 8 ) b
    Real ( kind = 8 ) b0(num+1)
    !    Real ( kind = 8 ) bcorr
    Real ( kind = 8 ) beta_asym
    Real ( kind = 8 ) bsum
    Real ( kind = 8 ) c(num+1)
    Real ( kind = 8 ) d(num+1)
    Real ( kind = 8 ) dsum
    Real ( kind = 8 ), Parameter :: e0 = 1.12837916709551D+00
    Real ( kind = 8 ), Parameter :: e1 = 0.353553390593274D+00
    Real ( kind = 8 ) eps
    !    Real ( kind = 8 ) error_fc
    Real ( kind = 8 ) f
    Real ( kind = 8 ) h
    Real ( kind = 8 ) h2
    Real ( kind = 8 ) hn
    Integer i
    Integer j
    Real ( kind = 8 ) j0
    Real ( kind = 8 ) j1
    Real ( kind = 8 ) lambda
    Integer m
    Integer mm1
    Integer mmj
    Integer n
    Integer np1
    Real ( kind = 8 ) r
    Real ( kind = 8 ) r0
    Real ( kind = 8 ) r1
    !    Real ( kind = 8 ) rlog1
    Real ( kind = 8 ) s
    Real ( kind = 8 ) sum1
    Real ( kind = 8 ) t
    Real ( kind = 8 ) t0
    Real ( kind = 8 ) t1
    Real ( kind = 8 ) u
    Real ( kind = 8 ) w
    Real ( kind = 8 ) w0
    Real ( kind = 8 ) z
    Real ( kind = 8 ) z0
    Real ( kind = 8 ) z2
    Real ( kind = 8 ) zn
    Real ( kind = 8 ) znm1

    beta_asym = 0.0D+00

    If ( a < b ) Then
       h = a / b
       r0 = 1.0D+00 / ( 1.0D+00 + h )
       r1 = ( b - a ) / b
       w0 = 1.0D+00 / Sqrt ( a * ( 1.0D+00 + h ))
    Else
       h = b / a
       r0 = 1.0D+00 / ( 1.0D+00 + h )
       r1 = ( b - a ) / a
       w0 = 1.0D+00 / Sqrt ( b * ( 1.0D+00 + h ))
    End If

    f = a * rlog1 ( - lambda / a ) + b * rlog1 ( lambda / b )
    t = Exp ( - f )
    If ( t == 0.0D+00 ) Then
       Return
    End If

    z0 = Sqrt ( f )
    z = 0.5D+00 * ( z0 / e1 )
    z2 = f + f

    a0(1) = ( 2.0D+00 / 3.0D+00 ) * r1
    c(1) = -0.5D+00 * a0(1)
    d(1) = -c(1)
    j0 = ( 0.5D+00 / e0 ) * error_fc ( 1, z0 )
    j1 = e1
    sum1 = j0 + d(1) * w0 * j1

    s = 1.0D+00
    h2 = h * h
    hn = 1.0D+00
    w = w0
    znm1 = z
    zn = z2

    Do n = 2, num, 2

       hn = h2 * hn
       a0(n) = 2.0D+00 * r0 * ( 1.0D+00 + h * hn ) &
            / ( n + 2.0D+00 )
       np1 = n + 1
       s = s + hn
       a0(np1) = 2.0D+00 * r1 * s / ( n + 3.0D+00 )

       Do i = n, np1

          r = -0.5D+00 * ( i + 1.0D+00 )
          b0(1) = r * a0(1)
          Do m = 2, i
             bsum = 0.0D+00
             mm1 = m - 1
             Do j = 1, mm1
                mmj = m - j
                bsum = bsum + ( j * r - mmj ) * a0(j) * b0(mmj)
             End Do
             b0(m) = r * a0(m) + bsum / m
          End Do

          c(i) = b0(i) / ( i + 1.0D+00 )

          dsum = 0.0
          Do j = 1, i-1
             dsum = dsum + d(i-j) * c(j)
          End Do
          d(i) = - ( dsum + c(i) )

       End Do

       j0 = e1 * znm1 + ( n - 1.0D+00 ) * j0
       j1 = e1 * zn + n * j1
       znm1 = z2 * znm1
       zn = z2 * zn
       w = w0 * w
       t0 = d(n) * w * j0
       w = w0 * w
       t1 = d(np1) * w * j1
       sum1 = sum1 + ( t0 + t1 )

       If ( ( Abs ( t0 ) + Abs ( t1 )) <= eps * sum1 ) Then
          u = Exp ( - bcorr ( a, b ) )
          beta_asym = e0 * t * u * sum1
          Return
       End If

    End Do

    u = Exp ( - bcorr ( a, b ) )
    beta_asym = e0 * t * u * sum1

    Return
  End Function beta_asym

  Function beta_frac ( a, b, x, y, lambda, eps )

    !*****************************************************************************80
    !
    !! BETA_FRAC evaluates a continued fraction expansion for IX(A,B).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, the parameters of the function.
    !    A and B should be nonnegative.  It is assumed that both A and
    !    B are greater than 1.
    !
    !    Input, real ( kind = 8 ) X, Y.  X is the argument of the
    !    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
    !
    !    Input, real ( kind = 8 ) LAMBDA, the value of ( A + B ) * Y - B.
    !
    !    Input, real ( kind = 8 ) EPS, a tolerance.
    !
    !    Output, real ( kind = 8 ) BETA_FRAC, the value of the continued
    !    fraction approximation for IX(A,B).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) alpha
    Real ( kind = 8 ) an
    Real ( kind = 8 ) anp1
    Real ( kind = 8 ) b
    Real ( kind = 8 ) beta
    Real ( kind = 8 ) beta_frac
    !    Real ( kind = 8 ) beta_rcomp
    Real ( kind = 8 ) bn
    Real ( kind = 8 ) bnp1
    Real ( kind = 8 ) c
    Real ( kind = 8 ) c0
    Real ( kind = 8 ) c1
    Real ( kind = 8 ) e
    Real ( kind = 8 ) eps
    Real ( kind = 8 ) lambda
    Real ( kind = 8 ) n
    Real ( kind = 8 ) p
    Real ( kind = 8 ) r
    Real ( kind = 8 ) r0
    Real ( kind = 8 ) s
    Real ( kind = 8 ) t
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    Real ( kind = 8 ) y
    Real ( kind = 8 ) yp1

    beta_frac = beta_rcomp ( a, b, x, y )

    If ( beta_frac == 0.0D+00 ) Then
       Return
    End If

    c = 1.0D+00 + lambda
    c0 = b / a
    c1 = 1.0D+00 + 1.0D+00 / a
    yp1 = y + 1.0D+00

    n = 0.0D+00
    p = 1.0D+00
    s = a + 1.0D+00
    an = 0.0D+00
    bn = 1.0D+00
    anp1 = 1.0D+00
    bnp1 = c / c1
    r = c1 / c
    !
    !  Continued fraction calculation.
    !
    Do

       n = n + 1.0D+00
       t = n / a
       w = n * ( b - n ) * x
       e = a / s
       alpha = ( p * ( p + c0 ) * e * e ) * ( w * x )
       e = ( 1.0D+00 + t ) / ( c1 + t + t )
       beta = n + w / s + e * ( c + n * yp1 )
       p = 1.0D+00 + t
       s = s +  2.0D+00
       !
       !  Update AN, BN, ANP1, and BNP1.
       !
       t = alpha * an + beta * anp1
       an = anp1
       anp1 = t
       t = alpha * bn + beta * bnp1
       bn = bnp1
       bnp1 = t

       r0 = r
       r = anp1 / bnp1

       If ( Abs ( r - r0 ) <= eps * r ) Then
          beta_frac = beta_frac * r
          Exit
       End If
       !
       !  Rescale AN, BN, ANP1, and BNP1.
       !
       an = an / bnp1
       bn = bn / bnp1
       anp1 = r
       bnp1 = 1.0D+00

    End Do

    Return
  End Function beta_frac

  Subroutine beta_grat ( a, b, x, y, w, eps, ierr )

    !*****************************************************************************80
    !
    !! BETA_GRAT evaluates an asymptotic expansion for IX(A,B).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, the parameters of the function.
    !    A and B should be nonnegative.  It is assumed that 15 <= A
    !    and B <= 1, and that B is less than A.
    !
    !    Input, real ( kind = 8 ) X, Y.  X is the argument of the
    !    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
    !
    !    Input/output, real ( kind = 8 ) W, a quantity to which the
    !    result of the computation is to be added on output.
    !
    !    Input, real ( kind = 8 ) EPS, a tolerance.
    !
    !    Output, integer IERR, an error flag, which is 0 if no error
    !    was detected.
    !
    Implicit None

    Real ( kind = 8 ) a
    !    Real ( kind = 8 ) algdiv
    !    Real ( kind = 8 ) alnrel
    Real ( kind = 8 ) b
    Real ( kind = 8 ) bm1
    Real ( kind = 8 ) bp2n
    Real ( kind = 8 ) c(30)
    Real ( kind = 8 ) cn
    Real ( kind = 8 ) coef
    Real ( kind = 8 ) d(30)
    Real ( kind = 8 ) dj
    Real ( kind = 8 ) eps
    !    Real ( kind = 8 ) gam1
    Integer i
    Integer ierr
    Real ( kind = 8 ) j
    Real ( kind = 8 ) l
    Real ( kind = 8 ) lnx
    Integer n
    Real ( kind = 8 ) n2
    Real ( kind = 8 ) nu
    Real ( kind = 8 ) p
    Real ( kind = 8 ) q
    Real ( kind = 8 ) r
    Real ( kind = 8 ) s
    Real ( kind = 8 ) sum1
    Real ( kind = 8 ) t
    Real ( kind = 8 ) t2
    Real ( kind = 8 ) u
    Real ( kind = 8 ) v
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    Real ( kind = 8 ) y
    Real ( kind = 8 ) z

    bm1 = ( b - 0.5D+00 ) - 0.5D+00
    nu = a + 0.5D+00 * bm1

    If ( y <= 0.375D+00 ) Then
       lnx = alnrel ( - y )
    Else
       lnx = Log ( x )
    End If

    z = -nu * lnx

    If ( b * z == 0.0D+00 ) Then
       ierr = 1
       Return
    End If
    !
    !  Computation of the expansion.
    !
    !  Set R = EXP(-Z)*Z**B/GAMMA(B)
    !
    r = b * ( 1.0D+00 + gam1 ( b ) ) * Exp ( b * Log ( z ))
    r = r * Exp ( a * lnx ) * Exp ( 0.5D+00 * bm1 * lnx )
    u = algdiv ( b, a ) + b * Log ( nu )
    u = r * Exp ( - u )

    If ( u == 0.0D+00 ) Then
       ierr = 1
       Return
    End If

    Call gamma_rat1 ( b, z, r, p, q, eps )

    v = 0.25D+00 * ( 1.0D+00 / nu )**2
    t2 = 0.25D+00 * lnx * lnx
    l = w / u
    j = q / r
    sum1 = j
    t = 1.0D+00
    cn = 1.0D+00
    n2 = 0.0D+00

    Do n = 1, 30

       bp2n = b + n2
       j = ( bp2n * ( bp2n + 1.0D+00 ) * j &
            + ( z + bp2n + 1.0D+00 ) * t ) * v
       n2 = n2 +  2.0D+00
       t = t * t2
       cn = cn / ( n2 * ( n2 + 1.0D+00 ))
       c(n) = cn
       s = 0.0D+00

       coef = b - n
       Do i = 1, n-1
          s = s + coef * c(i) * d(n-i)
          coef = coef + b
       End Do

       d(n) = bm1 * cn + s / n
       dj = d(n) * j
       sum1 = sum1 + dj

       If ( sum1 <= 0.0D+00 ) Then
          ierr = 1
          Return
       End If

       If ( Abs ( dj ) <= eps * ( sum1 + l ) ) Then
          ierr = 0
          w = w + u * sum1
          Return
       End If

    End Do

    ierr = 0
    w = w + u * sum1

    Return
  End Subroutine beta_grat

  Subroutine beta_inc ( a, b, x, y, w, w1, ierr )

    !*****************************************************************************80
    !
    !! BETA_INC evaluates the incomplete beta function IX(A,B).
    !
    !  Author:
    !
    !    Alfred Morris,
    !    Naval Surface Weapons Center,
    !    Dahlgren, Virginia.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, the parameters of the function.
    !    A and B should be nonnegative.
    !
    !    Input, real ( kind = 8 ) X, Y.  X is the argument of the
    !    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
    !
    !    Output, real ( kind = 8 ) W, W1, the values of IX(A,B) and
    !    1-IX(A,B).
    !
    !    Output, integer IERR, the error flag.
    !    0, no error was detected.
    !    1, A or B is negative;
    !    2, A = B = 0;
    !    3, X < 0 or 1 < X;
    !    4, Y < 0 or 1 < Y;
    !    5, X + Y /= 1;
    !    6, X = A = 0;
    !    7, Y = B = 0.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) a0
    !    Real ( kind = 8 ) apser
    Real ( kind = 8 ) b
    Real ( kind = 8 ) b0
    !    Real ( kind = 8 ) beta_asym
    !    Real ( kind = 8 ) beta_frac
    !    Real ( kind = 8 ) beta_pser
    !    Real ( kind = 8 ) beta_up
    Real ( kind = 8 ) eps
    !    Real ( kind = 8 ) fpser
    Integer ierr
    Integer ierr1
    Integer ind
    Real ( kind = 8 ) lambda
    Integer n
    Real ( kind = 8 ) t
    Real ( kind = 8 ) w
    Real ( kind = 8 ) w1
    Real ( kind = 8 ) x
    Real ( kind = 8 ) x0
    Real ( kind = 8 ) y
    Real ( kind = 8 ) y0
    Real ( kind = 8 ) z

    eps = Epsilon ( eps )
    w = 0.0D+00
    w1 = 0.0D+00

    If ( a < 0.0D+00 .Or. b < 0.0D+00 ) Then
       ierr = 1
       Return
    End If

    If ( a == 0.0D+00 .And. b == 0.0D+00 ) Then
       ierr = 2
       Return
    End If

    If ( x < 0.0D+00 .Or. 1.0D+00 < x ) Then
       ierr = 3
       Return
    End If

    If ( y < 0.0D+00 .Or. 1.0D+00 < y ) Then
       ierr = 4
       Return
    End If

    z = ( ( x + y ) - 0.5D+00 ) - 0.5D+00

    If ( 3.0D+00 * eps < Abs ( z ) ) Then
       ierr = 5
       Return
    End If

    ierr = 0

    If ( x == 0.0D+00 ) Then
       w = 0.0D+00
       w1 = 1.0D+00
       If ( a == 0.0D+00 ) Then
          ierr = 6
       End If
       Return
    End If

    If ( y == 0.0D+00 ) Then
       If ( b == 0.0D+00 ) Then
          ierr = 7
          Return
       End If
       w = 1.0D+00
       w1 = 0.0D+00
       Return
    End If

    If ( a == 0.0D+00 ) Then
       w = 1.0D+00
       w1 = 0.0D+00
       Return
    End If

    If ( b == 0.0D+00 ) Then
       w = 0.0D+00
       w1 = 1.0D+00
       Return
    End If

    eps = Max ( eps, 1.0D-15 )

    If ( Max ( a, b ) < 0.001D+00 * eps ) Then
       go to 260
    End If

    ind = 0
    a0 = a
    b0 = b
    x0 = x
    y0 = y

    If ( 1.0D+00 < Min ( a0, b0 ) ) Then
       go to 40
    End If
    !
    !  Procedure for A0 <= 1 or B0 <= 1
    !
    If ( 0.5D+00 < x ) Then
       ind = 1
       a0 = b
       b0 = a
       x0 = y
       y0 = x
    End If

    If ( b0 < Min ( eps, eps * a0 ) ) Then
       go to 90
    End If

    If ( a0 < Min ( eps, eps * b0 ) .And. b0 * x0 <= 1.0D+00 ) Then
       go to 100
    End If

    If ( 1.0D+00 < Max ( a0, b0 ) ) Then
       go to 20
    End If

    If ( Min ( 0.2D+00, b0 ) <= a0 ) Then
       go to 110
    End If

    If ( x0**a0 <= 0.9D+00 ) Then
       go to 110
    End If

    If ( 0.3D+00 <= x0 ) Then
       go to 120
    End If

    n = 20
    go to 140

20  Continue

    If ( b0 <= 1.0D+00 ) Then
       go to 110
    End If

    If ( 0.3D+00 <= x0 ) Then
       go to 120
    End If

    If ( 0.1D+00 <= x0 ) Then
       go to 30
    End If

    If ( (x0 * b0 )**a0 <= 0.7D+00 ) Then
       go to 110
    End If

30  Continue

    If ( 15.0D+00 < b0 ) Then
       go to 150
    End If

    n = 20
    go to 140
    !
    !  PROCEDURE for 1 < A0 and 1 < B0.
    !
40  Continue

    If ( a <= b ) Then
       lambda = a - ( a + b ) * x
    Else
       lambda = ( a + b ) * y - b
    End If

    If ( lambda < 0.0D+00 ) Then
       ind = 1
       a0 = b
       b0 = a
       x0 = y
       y0 = x
       lambda = Abs ( lambda )
    End If

! 70  Continue

    If ( b0 < 40.0D+00 .And. b0 * x0 <= 0.7D+00 ) Then
       go to 110
    End If

    If ( b0 < 40.0D+00 ) Then
       go to 160
    End If

    If ( b0 < a0 ) Then
       go to 80
    End If

    If ( a0 <= 100.0D+00 ) Then
       go to 130
    End If

    If ( 0.03D+00 * a0 < lambda ) Then
       go to 130
    End If

    go to 200

80  Continue

    If ( b0 <= 100.0D+00 ) Then
       go to 130
    End If

    If ( 0.03D+00 * b0 < lambda ) Then
       go to 130
    End If

    go to 200
    !
    !  Evaluation of the appropriate algorithm.
    !
90  Continue

    w = fpser ( a0, b0, x0, eps )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250

100 Continue

    w1 = apser ( a0, b0, x0, eps )
    w = 0.5D+00 + ( 0.5D+00 - w1 )
    go to 250

110 Continue

    w = beta_pser ( a0, b0, x0, eps )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250

120 Continue

    w1 = beta_pser ( b0, a0, y0, eps )
    w = 0.5D+00 + ( 0.5D+00 - w1 )
    go to 250

130 Continue

    w = beta_frac ( a0, b0, x0, y0, lambda, 15.0D+00 * eps )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250

140 Continue

    w1 = beta_up ( b0, a0, y0, x0, n, eps )
    b0 = b0 + n

150 Continue

    Call beta_grat ( b0, a0, y0, x0, w1, 15.0D+00 * eps, ierr1 )
    w = 0.5D+00 + ( 0.5D+00 - w1 )
    go to 250

160 Continue

    n = INT(b0)
    b0 = b0 - n

    If ( b0 == 0.0D+00 ) Then
       n = n - 1
       b0 = 1.0D+00
    End If

! 170 Continue

    w = beta_up ( b0, a0, y0, x0, n, eps )

    If ( x0 <= 0.7D+00 ) Then
       w = w + beta_pser ( a0, b0, x0, eps )
       w1 = 0.5D+00 + ( 0.5D+00 - w )
       go to 250
    End If

    If ( a0 <= 15.0D+00 ) Then
       n = 20
       w = w + beta_up ( a0, b0, x0, y0, n, eps )
       a0 = a0 + n
    End If

! 190 Continue

    Call beta_grat ( a0, b0, x0, y0, w, 15.0D+00 * eps, ierr1 )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250

200 Continue

    w = beta_asym ( a0, b0, lambda, 100.0D+00 * eps )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250
    !
    !  Termination of the procedure.
    !
250 Continue

    If ( ind /= 0 ) Then
       t = w
       w = w1
       w1 = t
    End If

    Return
    !
    !  Procedure for A and B < 0.001 * EPS
    !
260 Continue

    w = b / ( a + b )
    w1 = a / ( a + b )

    Return
  End Subroutine beta_inc

!!$  Subroutine beta_inc_values ( n_data, a, b, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! BETA_INC_VALUES returns some values of the incomplete Beta function.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    The incomplete Beta function may be written
!!$    !
!!$    !      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
!!$    !                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
!!$    !
!!$    !    Thus,
!!$    !
!!$    !      BETA_INC(A,B,0.0) = 0.0
!!$    !      BETA_INC(A,B,1.0) = 1.0
!!$    !
!!$    !    Note that in Mathematica, the expressions:
!!$    !
!!$    !      BETA[A,B]   = Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
!!$    !      BETA[X,A,B] = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
!!$    !
!!$    !    and thus, to evaluate the incomplete Beta function requires:
!!$    !
!!$    !      BETA_INC(A,B,X) = BETA[X,A,B] / BETA[A,B]
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    17 February 2004
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !    Karl Pearson,
!!$    !    Tables of the Incomplete Beta Function,
!!$    !    Cambridge University Press, 1968.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, real ( kind = 8 ) A, B, X, the arguments of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 30
!!$
!!$    Real ( kind = 8 ) a
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: a_vec = (/ &
!!$         0.5D+00,  0.5D+00,  0.5D+00,  1.0D+00, &
!!$         1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
!!$         2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
!!$         2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
!!$         2.0D+00,  5.5D+00, 10.0D+00, 10.0D+00, &
!!$         10.0D+00, 10.0D+00, 20.0D+00, 20.0D+00, &
!!$         20.0D+00, 20.0D+00, 20.0D+00, 30.0D+00, &
!!$         30.0D+00, 40.0D+00 /)
!!$    Real ( kind = 8 ) b
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: b_vec = (/ &
!!$         0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
!!$         0.5D+00,  0.5D+00,  0.5D+00,  1.0D+00, &
!!$         2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
!!$         2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
!!$         2.0D+00,  5.0D+00,  0.5D+00,  5.0D+00, &
!!$         5.0D+00, 10.0D+00,  5.0D+00, 10.0D+00, &
!!$         10.0D+00, 20.0D+00, 20.0D+00, 10.0D+00, &
!!$         10.0D+00, 20.0D+00 /)
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.0637686D+00, 0.2048328D+00, 1.0000000D+00, 0.0D+00,       &
!!$         0.0050126D+00, 0.0513167D+00, 0.2928932D+00, 0.5000000D+00, &
!!$         0.028D+00,     0.104D+00,     0.216D+00,     0.352D+00,     &
!!$         0.500D+00,     0.648D+00,     0.784D+00,     0.896D+00,     &
!!$         0.972D+00,     0.4361909D+00, 0.1516409D+00, 0.0897827D+00, &
!!$         1.0000000D+00, 0.5000000D+00, 0.4598773D+00, 0.2146816D+00, &
!!$         0.9507365D+00, 0.5000000D+00, 0.8979414D+00, 0.2241297D+00, &
!!$         0.7586405D+00, 0.7001783D+00 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         0.01D+00, 0.10D+00, 1.00D+00, 0.0D+00,  &
!!$         0.01D+00, 0.10D+00, 0.50D+00, 0.50D+00, &
!!$         0.1D+00,  0.2D+00,  0.3D+00,  0.4D+00,  &
!!$         0.5D+00,  0.6D+00,  0.7D+00,  0.8D+00,  &
!!$         0.9D+00,  0.50D+00, 0.90D+00, 0.50D+00, &
!!$         1.00D+00, 0.50D+00, 0.80D+00, 0.60D+00, &
!!$         0.80D+00, 0.50D+00, 0.60D+00, 0.70D+00, &
!!$         0.80D+00, 0.70D+00 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       a = 0.0D+00
!!$       b = 0.0D+00
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       a = a_vec(n_data)
!!$       b = b_vec(n_data)
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine beta_inc_values

  Function beta_log ( a0, b0 )

    !*****************************************************************************80
    !
    !! BETA_LOG evaluates the logarithm of the beta function.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A0, B0, the parameters of the function.
    !    A0 and B0 should be nonnegative.
    !
    !    Output, real ( kind = 8 ) BETA_LOG, the value of the logarithm
    !    of the Beta function.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) a0
    !    Real ( kind = 8 ) algdiv
    !    Real ( kind = 8 ) alnrel
    Real ( kind = 8 ) b
    Real ( kind = 8 ) b0
    !    Real ( kind = 8 ) bcorr
    Real ( kind = 8 ) beta_log
    Real ( kind = 8 ) c
    Real ( kind = 8 ), Parameter :: e = 0.918938533204673D+00
    !    Real ( kind = 8 ) gamma_log
    !    Real ( kind = 8 ) gsumln
    Real ( kind = 8 ) h
    Integer i
    Integer n
    Real ( kind = 8 ) u
    Real ( kind = 8 ) v
    Real ( kind = 8 ) w
    Real ( kind = 8 ) z

    a = Min ( a0, b0 )
    b = Max ( a0, b0 )

    If ( 8.0D+00 <= a ) Then
       go to 100
    End If

    If ( 1.0D+00 <= a ) Then
       go to 20
    End If
    !
    !  Procedure when A < 1
    !
    If ( b < 8.0D+00 ) Then
       beta_log = gamma_log ( a ) + ( gamma_log ( b ) - gamma_log ( a + b ) )
    Else
       beta_log = gamma_log ( a ) + algdiv ( a, b )
    End If

    Return
    !
    !  Procedure when 1 <= A < 8
    !
20  Continue

    If ( 2.0D+00 < a ) Then
       go to 40
    End If

    If ( b <= 2.0D+00 ) Then
       beta_log = gamma_log ( a ) + gamma_log ( b ) - gsumln ( a, b )
       Return
    End If

    w = 0.0D+00
    If ( b < 8.0D+00 ) Then
       go to 60
    End If

    beta_log = gamma_log ( a ) + algdiv ( a, b )
    Return
    !
    !  Reduction of A when B <= 1000.
    !
40  Continue

    If ( 1000.0D+00 < b ) Then
       go to 80
    End If

    n = INT(a - 1.0D+00)
    w = 1.0D+00
    Do i = 1, n
       a = a - 1.0D+00
       h = a / b
       w = w * ( h / ( 1.0D+00 + h ) )
    End Do
    w = Log ( w )

    If ( 8.0D+00 <= b ) Then
       beta_log = w + gamma_log ( a ) + algdiv ( a, b )
       Return
    End If
    !
    !  Reduction of B when B < 8.
    !
60  Continue

    n = INT(b - 1.0D+00)
    z = 1.0D+00
    Do i = 1, n
       b = b - 1.0D+00
       z = z * ( b / ( a + b ))
    End Do

    beta_log = w + Log ( z ) + ( gamma_log ( a ) + ( gamma_log ( b ) &
         - gsumln ( a, b ) ) )
    Return
    !
    !  Reduction of A when 1000 < B.
    !
80  Continue

    n = INT(a - 1.0D+00)
    w = 1.0D+00
    Do i = 1, n
       a = a - 1.0D+00
       w = w * ( a / ( 1.0D+00 + a / b ))
    End Do

    beta_log = ( Log ( w ) - n * Log ( b ) ) + ( gamma_log ( a ) + algdiv ( a, b ) )

    Return
    !
    !  Procedure when 8 <= A.
    !
100 Continue

    w = bcorr ( a, b )
    h = a / b
    c = h / ( 1.0D+00 + h )
    u = - ( a - 0.5D+00 ) * Log ( c )
    v = b * alnrel ( h )

    If ( v < u ) Then
       beta_log = ((( -0.5D+00 * Log ( b ) + e ) + w ) - v ) - u
    Else
       beta_log = ((( -0.5D+00 * Log ( b ) + e ) + w ) - u ) - v
    End If

    Return
  End Function beta_log

  Function beta_pser ( a, b, x, eps )

    !*****************************************************************************80
    !
    !! BETA_PSER uses a power series expansion to evaluate IX(A,B)(X).
    !
    !  Discussion:
    !
    !    BETA_PSER is used when B <= 1 or B*X <= 0.7.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, the parameters.
    !
    !    Input, real ( kind = 8 ) X, the point where the function
    !    is to be evaluated.
    !
    !    Input, real ( kind = 8 ) EPS, the tolerance.
    !
    !    Output, real ( kind = 8 ) BETA_PSER, the approximate value of IX(A,B)(X).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) a0
    !    Real ( kind = 8 ) algdiv
    Real ( kind = 8 ) apb
    Real ( kind = 8 ) b
    Real ( kind = 8 ) b0
    !    Real ( kind = 8 ) beta_log
    Real ( kind = 8 ) beta_pser
    Real ( kind = 8 ) c
    Real ( kind = 8 ) eps
    !    Real ( kind = 8 ) gam1
    !    Real ( kind = 8 ) gamma_ln1
    Integer i
    Integer m
    Real ( kind = 8 ) n
    Real ( kind = 8 ) sum1
    Real ( kind = 8 ) t
    Real ( kind = 8 ) tol
    Real ( kind = 8 ) u
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    Real ( kind = 8 ) z

    beta_pser = 0.0D+00

    If ( x == 0.0D+00 ) Then
       Return
    End If
    !
    !  Compute the factor X**A/(A*BETA(A,B))
    !
    a0 = Min ( a, b )

    If ( 1.0D+00 <= a0 ) Then

       z = a * Log ( x ) - beta_log ( a, b )
       beta_pser = Exp ( z ) / a

    Else

       b0 = Max ( a, b )

       If ( b0 <= 1.0D+00 ) Then

          beta_pser = x**a
          If ( beta_pser == 0.0D+00 ) Then
             Return
          End If

          apb = a + b

          If ( apb <= 1.0D+00 ) Then
             z = 1.0D+00 + gam1 ( apb )
          Else
             u = a + b - 1.0D+00
             z = ( 1.0D+00 + gam1 ( u ) ) / apb
          End If

          c = ( 1.0D+00 + gam1 ( a ) ) &
               * ( 1.0D+00 + gam1 ( b ) ) / z
          beta_pser = beta_pser * c * ( b / apb )

       Else If ( b0 < 8.0D+00 ) Then

          u = gamma_ln1 ( a0 )
          m = INT(b0 - 1.0D+00)

          c = 1.0D+00
          Do i = 1, m
             b0 = b0 - 1.0D+00
             c = c * ( b0 / ( a0 + b0 ))
          End Do

          u = Log ( c ) + u
          z = a * Log ( x ) - u
          b0 = b0 - 1.0D+00
          apb = a0 + b0

          If ( apb <= 1.0D+00 ) Then
             t = 1.0D+00 + gam1 ( apb )
          Else
             u = a0 + b0 - 1.0D+00
             t = ( 1.0D+00 + gam1 ( u ) ) / apb
          End If

          beta_pser = Exp ( z ) * ( a0 / a ) &
               * ( 1.0D+00 + gam1 ( b0 )) / t

       Else If ( 8.0D+00 <= b0 ) Then

          u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
          z = a * Log ( x ) - u
          beta_pser = ( a0 / a ) * Exp ( z )

       End If

    End If

    If ( beta_pser == 0.0D+00 .Or. a <= 0.1D+00 * eps ) Then
       Return
    End If
    !
    !  Compute the series.
    !
    sum1 = 0.0D+00
    n = 0.0D+00
    c = 1.0D+00
    tol = eps / a

    Do

       n = n + 1.0D+00
       c = c * ( 0.5D+00 + ( 0.5D+00 - b / n ) ) * x
       w = c / ( a + n )
       sum1 = sum1 + w

       If ( Abs ( w ) <= tol ) Then
          Exit
       End If

    End Do

    beta_pser = beta_pser * ( 1.0D+00 + a * sum1 )

    Return
  End Function beta_pser

  Function beta_rcomp ( a, b, x, y )

    !*****************************************************************************80
    !
    !! BETA_RCOMP evaluates X**A * Y**B / Beta(A,B).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, the parameters of the Beta function.
    !    A and B should be nonnegative.
    !
    !    Input, real ( kind = 8 ) X, Y, define the numerator of the fraction.
    !
    !    Output, real ( kind = 8 ) BETA_RCOMP, the value of X**A * Y**B / Beta(A,B).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) a0
    !    Real ( kind = 8 ) algdiv
    !    Real ( kind = 8 ) alnrel
    Real ( kind = 8 ) apb
    Real ( kind = 8 ) b
    Real ( kind = 8 ) b0
    !    Real ( kind = 8 ) bcorr
    !    Real ( kind = 8 ) beta_log
    Real ( kind = 8 ) beta_rcomp
    Real ( kind = 8 ) c
    Real ( kind = 8 ), Parameter :: const = 0.398942280401433D+00
    Real ( kind = 8 ) e
    !    Real ( kind = 8 ) gam1
    !    Real ( kind = 8 ) gamma_ln1
    Real ( kind = 8 ) h
    Integer i
    Real ( kind = 8 ) lambda
    Real ( kind = 8 ) lnx
    Real ( kind = 8 ) lny
    Integer n
    !    Real ( kind = 8 ) rlog1
    Real ( kind = 8 ) t
    Real ( kind = 8 ) u
    Real ( kind = 8 ) v
    Real ( kind = 8 ) x
    Real ( kind = 8 ) x0
    Real ( kind = 8 ) y
    Real ( kind = 8 ) y0
    Real ( kind = 8 ) z

    beta_rcomp = 0.0D+00
    If ( x == 0.0D+00 .Or. y == 0.0D+00 ) Then
       Return
    End If

    a0 = Min ( a, b )

    If ( a0 < 8.0D+00 ) Then

       If ( x <= 0.375D+00 ) Then
          lnx = Log ( x )
          lny = alnrel ( - x )
       Else If ( y <= 0.375D+00 ) Then
          lnx = alnrel ( - y )
          lny = Log ( y )
       Else
          lnx = Log ( x )
          lny = Log ( y )
       End If

       z = a * lnx + b * lny

       If ( 1.0D+00 <= a0 ) Then
          z = z - beta_log ( a, b )
          beta_rcomp = Exp ( z )
          Return
       End If
       !
       !  Procedure for A < 1 or B < 1
       !
       b0 = Max ( a, b )

       If ( b0 <= 1.0D+00 ) Then

          beta_rcomp = Exp ( z )
          If ( beta_rcomp == 0.0D+00 ) Then
             Return
          End If

          apb = a + b

          If ( apb <= 1.0D+00 ) Then
             z = 1.0D+00 + gam1 ( apb )
          Else
             u = a + b - 1.0D+00
             z = ( 1.0D+00 + gam1 ( u ) ) / apb
          End If

          c = ( 1.0D+00 + gam1 ( a ) ) &
               * ( 1.0D+00 + gam1 ( b ) ) / z
          beta_rcomp = beta_rcomp * ( a0 * c ) &
               / ( 1.0D+00 + a0 / b0 )

       Else If ( b0 < 8.0D+00 ) Then

          u = gamma_ln1 ( a0 )
          n = INT(b0 - 1.0D+00)

          c = 1.0D+00
          Do i = 1, n
             b0 = b0 - 1.0D+00
             c = c * ( b0 / ( a0 + b0 ))
          End Do
          u = Log ( c ) + u

          z = z - u
          b0 = b0 - 1.0D+00
          apb = a0 + b0

          If ( apb <= 1.0D+00 ) Then
             t = 1.0D+00 + gam1 ( apb )
          Else
             u = a0 + b0 - 1.0D+00
             t = ( 1.0D+00 + gam1 ( u ) ) / apb
          End If

          beta_rcomp = a0 * Exp ( z ) * ( 1.0D+00 + gam1 ( b0 ) ) / t

       Else If ( 8.0D+00 <= b0 ) Then

          u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
          beta_rcomp = a0 * Exp ( z - u )

       End If

    Else

       If ( a <= b ) Then
          h = a / b
          x0 = h / ( 1.0D+00 + h )
          y0 = 1.0D+00 / (  1.0D+00 + h )
          lambda = a - ( a + b ) * x
       Else
          h = b / a
          x0 = 1.0D+00 / ( 1.0D+00 + h )
          y0 = h / ( 1.0D+00 + h )
          lambda = ( a + b ) * y - b
       End If

       e = -lambda / a

       If ( Abs ( e ) <= 0.6D+00 ) Then
          u = rlog1 ( e )
       Else
          u = e - Log ( x / x0 )
       End If

       e = lambda / b

       If ( Abs ( e ) <= 0.6D+00 ) Then
          v = rlog1 ( e )
       Else
          v = e - Log ( y / y0 )
       End If

       z = Exp ( - ( a * u + b * v ) )
       beta_rcomp = const * Sqrt ( b * x0 ) * z * Exp ( - bcorr ( a, b ))

    End If

    Return
  End Function beta_rcomp

  Function beta_rcomp1 ( mu, a, b, x, y )

    !*****************************************************************************80
    !
    !! BETA_RCOMP1 evaluates exp(MU) * X**A * Y**B / Beta(A,B).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, integer MU, ?
    !
    !    Input, real ( kind = 8 ) A, B, the parameters of the Beta function.
    !    A and B should be nonnegative.
    !
    !    Input, real ( kind = 8 ) X, Y, ?
    !
    !    Output, real ( kind = 8 ) BETA_RCOMP1, the value of
    !    exp(MU) * X**A * Y**B / Beta(A,B).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) a0
    !    Real ( kind = 8 ) algdiv
    !    Real ( kind = 8 ) alnrel
    Real ( kind = 8 ) apb
    Real ( kind = 8 ) b
    Real ( kind = 8 ) b0
    !    Real ( kind = 8 ) bcorr
    !    Real ( kind = 8 ) beta_log
    Real ( kind = 8 ) beta_rcomp1
    Real ( kind = 8 ) c
    Real ( kind = 8 ), Parameter :: const = 0.398942280401433D+00
    Real ( kind = 8 ) e
    !    Real ( kind = 8 ) esum
    !    Real ( kind = 8 ) gam1
    !    Real ( kind = 8 ) gamma_ln1
    Real ( kind = 8 ) h
    Integer i
    Real ( kind = 8 ) lambda
    Real ( kind = 8 ) lnx
    Real ( kind = 8 ) lny
    Integer mu
    Integer n
    !    Real ( kind = 8 ) rlog1
    Real ( kind = 8 ) t
    Real ( kind = 8 ) u
    Real ( kind = 8 ) v
    Real ( kind = 8 ) x
    Real ( kind = 8 ) x0
    Real ( kind = 8 ) y
    Real ( kind = 8 ) y0
    Real ( kind = 8 ) z

    a0 = Min ( a, b )

    If ( 8.0D+00 <= a0 ) Then
       go to 130
    End If

    If ( x <= 0.375D+00 ) Then
       lnx = Log ( x )
       lny = alnrel ( - x )
    Else If ( y <= 0.375D+00 ) Then
       lnx = alnrel ( - y )
       lny = Log ( y )
    Else
       lnx = Log ( x )
       lny = Log ( y )
    End If

    z = a * lnx + b * lny

    If ( 1.0D+00 <= a0 ) Then
       z = z - beta_log ( a, b )
       beta_rcomp1 = esum ( mu, z )
       Return
    End If
    !
    !  Procedure for A < 1 OR B < 1
    !
! 40  Continue

    b0 = Max ( a, b )

    If ( 8.0D+00 <= b0 ) Then
       go to 120
    End If

    If ( 1.0D+00 < b0 ) Then
       go to 70
    End If
    !
    !  Algorithm for B0 <= 1
    !
    beta_rcomp1 = esum ( mu, z )
    If ( beta_rcomp1 == 0.0D+00 ) Then
       Return
    End If

    apb = a + b

    If ( apb <= 1.0D+00 ) Then
       z = 1.0D+00 + gam1 ( apb )
    Else
       u = Real ( a, kind = 8 ) + Real ( b, kind = 8 ) - 1.0D+00
       z = ( 1.0D+00 + gam1 ( u )) / apb
    End If

    c = ( 1.0D+00 + gam1 ( a ) ) &
         * ( 1.0D+00 + gam1 ( b ) ) / z
    beta_rcomp1 = beta_rcomp1 * ( a0 * c ) / ( 1.0D+00 + a0 / b0 )
    Return
    !
    !  Algorithm for 1 < B0 < 8
    !
70  Continue

    u = gamma_ln1 ( a0 )
    n = INT(b0 - 1.0D+00)

    c = 1.0D+00
    Do i = 1, n
       b0 = b0 - 1.0D+00
       c = c * ( b0 / ( a0 + b0 ) )
    End Do
    u = Log ( c ) + u

    z = z - u
    b0 = b0 - 1.0D+00
    apb = a0 + b0

    If ( apb <= 1.0D+00 ) Then
       t = 1.0D+00 + gam1 ( apb )
    Else
       u = a0 + b0 - 1.0D+00
       t = ( 1.0D+00 + gam1 ( u ) ) / apb
    End If

    beta_rcomp1 = a0 * esum ( mu, z ) &
         * ( 1.0D+00 + gam1 ( b0 ) ) / t
    Return
    !
    !  Algorithm for 8 <= B0.
    !
120 Continue

    u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
    beta_rcomp1 = a0 * esum ( mu, z-u )
    Return
    !
    !  Procedure for 8 <= A and 8 <= B.
    !
130 Continue

    If ( a <= b ) Then
       h = a / b
       x0 = h / ( 1.0D+00 + h )
       y0 = 1.0D+00 / ( 1.0D+00 + h )
       lambda = a - ( a + b ) * x
    Else
       h = b / a
       x0 = 1.0D+00 / ( 1.0D+00 + h )
       y0 = h / ( 1.0D+00 + h )
       lambda = ( a + b ) * y - b
    End If

    e = -lambda / a

    If ( Abs ( e ) <= 0.6D+00 ) Then
       u = rlog1 ( e )
    Else
       u = e - Log ( x / x0 )
    End If

    e = lambda / b

    If ( Abs ( e ) <= 0.6D+00 ) Then
       v = rlog1 ( e )
    Else
       v = e - Log ( y / y0 )
    End If

    z = esum ( mu, - ( a * u + b * v ))
    beta_rcomp1 = const * Sqrt ( b * x0 ) * z * Exp ( - bcorr ( a, b ) )

    Return
  End Function beta_rcomp1

  Function beta_up ( a, b, x, y, n, eps )

    !*****************************************************************************80
    !
    !! BETA_UP evaluates IX(A,B) - IX(A+N,B) where N is a positive integer.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, the parameters of the function.
    !    A and B should be nonnegative.
    !
    !    Input, real ( kind = 8 ) X, Y, ?
    !
    !    Input, integer N, ?
    !
    !    Input, real ( kind = 8 ) EPS, the tolerance.
    !
    !    Output, real ( kind = 8 ) BETA_UP, the value of IX(A,B) - IX(A+N,B).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) ap1
    Real ( kind = 8 ) apb
    Real ( kind = 8 ) b
    !    Real ( kind = 8 ) beta_rcomp1
    Real ( kind = 8 ) beta_up
    Real ( kind = 8 ) d
    Real ( kind = 8 ) eps
    !    Real ( kind = 8 ) exparg
    Integer i
    Integer k
    Real ( kind = 8 ) l
    Integer mu
    Integer n
    Real ( kind = 8 ) r
    Real ( kind = 8 ) t
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    Real ( kind = 8 ) y
    !
    !  Obtain the scaling factor EXP(-MU) AND
    !  EXP(MU)*(X**A*Y**B/BETA(A,B))/A
    !
    apb = a + b
    ap1 = a + 1.0D+00
    mu = 0
    d = 1.0D+00

    If ( n /= 1 ) Then

       If ( 1.0D+00 <= a ) Then

          If ( 1.1D+00 * ap1 <= apb ) Then
             mu = INT(Abs ( exparg ( 1 ) ))
             k = INT(exparg ( 0 ))
             If ( k < mu ) Then
                mu = k
             End If
             t = mu
             d = Exp ( - t )
          End If

       End If

    End If

    beta_up = beta_rcomp1 ( mu, a, b, x, y ) / a

    If ( n == 1 .Or. beta_up == 0.0D+00 ) Then
       Return
    End If

    w = d
    !
    !  Let K be the index of the maximum term.
    !
    k = 0

    If ( 1.0D+00 < b ) Then

       If ( y <= 0.0001D+00 ) Then

          k = n - 1

       Else

          r = ( b - 1.0D+00 ) * x / y - a

          If ( 1.0D+00 <= r ) Then
             k = n - 1
             t = n - 1
             If ( r < t ) Then
                k = INT(r)
             End If
          End If

       End If
       !
       !  Add the increasing terms of the series.
       !
       Do i = 1, k
          l = i - 1
          d = ( ( apb + l ) / ( ap1 + l ) ) * x * d
          w = w + d
       End Do

    End If
    !
    !  Add the remaining terms of the series.
    !
    Do i = k+1, n-1
       l = i - 1
       d = ( ( apb + l ) / ( ap1 + l ) ) * x * d
       w = w + d
       If ( d <= eps * w ) Then
          beta_up = beta_up * w
          Return
       End If
    End Do

    beta_up = beta_up * w

    Return
  End Function beta_up

!!$  Subroutine binomial_cdf_values ( n_data, a, b, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    CDF(X)(A,B) is the probability of at most X successes in A trials,
!!$    !    given that the probability of success on a single trial is B.
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    27 May 2001
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !    Daniel Zwillinger,
!!$    !    CRC Standard Mathematical Tables and Formulae,
!!$    !    30th Edition, CRC Press, 1996, pages 651-652.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, integer A, real ( kind = 8 ) B, integer X, the arguments
!!$    !    of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 17
!!$
!!$    Integer a
!!$    Integer, Save, Dimension ( n_max ) :: a_vec = (/ &
!!$         2,  2,  2,  2, &
!!$         2,  4,  4,  4, &
!!$         4, 10, 10, 10, &
!!$         10, 10, 10, 10, &
!!$         10 /)
!!$    Real ( kind = 8 ) b
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: b_vec = (/ &
!!$         0.05D+00, 0.05D+00, 0.05D+00, 0.50D+00, &
!!$         0.50D+00, 0.25D+00, 0.25D+00, 0.25D+00, &
!!$         0.25D+00, 0.05D+00, 0.10D+00, 0.15D+00, &
!!$         0.20D+00, 0.25D+00, 0.30D+00, 0.40D+00, &
!!$         0.50D+00 /)
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.9025D+00, 0.9975D+00, 1.0000D+00, 0.2500D+00, &
!!$         0.7500D+00, 0.3164D+00, 0.7383D+00, 0.9492D+00, &
!!$         0.9961D+00, 0.9999D+00, 0.9984D+00, 0.9901D+00, &
!!$         0.9672D+00, 0.9219D+00, 0.8497D+00, 0.6331D+00, &
!!$         0.3770D+00 /)
!!$    Integer n_data
!!$    Integer x
!!$    Integer, Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         0, 1, 2, 0, &
!!$         1, 0, 1, 2, &
!!$         3, 4, 4, 4, &
!!$         4, 4, 4, 4, &
!!$         4 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       a = 0
!!$       b = 0.0D+00
!!$       x = 0
!!$       fx = 0.0D+00
!!$    Else
!!$       a = a_vec(n_data)
!!$       b = b_vec(n_data)
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine binomial_cdf_values

  Subroutine cdfbet ( which, p, q, x, y, a, b, status, bound )

    !*****************************************************************************80
    !
    !! CDFBET evaluates the CDF of the Beta Distribution.
    !
    !  Discussion:
    !
    !    This routine calculates any one parameter of the beta distribution
    !    given the others.
    !
    !    The value P of the cumulative distribution function is calculated
    !    directly by code associated with the reference.
    !
    !    Computation of the other parameters involves a seach for a value that
    !    produces the desired value of P.  The search relies on the
    !    monotonicity of P with respect to the other parameters.
    !
    !    The beta density is proportional to t^(A-1) * (1-t)^(B-1).
    !
    !  Modified:
    !
    !    09 June 2004
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, integer WHICH, indicates which of the next four argument
    !    values is to be calculated from the others.
    !    1: Calculate P and Q from X, Y, A and B;
    !    2: Calculate X and Y from P, Q, A and B;
    !    3: Calculate A from P, Q, X, Y and B;
    !    4: Calculate B from P, Q, X, Y and A.
    !
    !    Input/output, real ( kind = 8 ) P, the integral from 0 to X of the
    !    chi-square distribution.  Input range: [0, 1].
    !
    !    Input/output, real ( kind = 8 ) Q, equals 1-P.  Input range: [0, 1].
    !
    !    Input/output, real ( kind = 8 ) X, the upper limit of integration
    !    of the beta density.  If it is an input value, it should lie in
    !    the range [0,1].  If it is an output value, it will be searched for
    !    in the range [0,1].
    !
    !    Input/output, real ( kind = 8 ) Y, equal to 1-X.  If it is an input
    !    value, it should lie in the range [0,1].  If it is an output value,
    !    it will be searched for in the range [0,1].
    !
    !    Input/output, real ( kind = 8 ) A, the first parameter of the beta
    !    density.  If it is an input value, it should lie in the range
    !    (0, +infinity).  If it is an output value, it will be searched
    !    for in the range [1D-300,1D300].
    !
    !    Input/output, real ( kind = 8 ) B, the second parameter of the beta
    !    density.  If it is an input value, it should lie in the range
    !    (0, +infinity).  If it is an output value, it will be searched
    !    for in the range [1D-300,1D300].
    !
    !    Output, integer STATUS, reports the status of the computation.
    !     0, if the calculation completed correctly;
    !    -I, if the input parameter number I is out of range;
    !    +1, if the answer appears to be lower than lowest search bound;
    !    +2, if the answer appears to be higher than greatest search bound;
    !    +3, if P + Q /= 1;
    !    +4, if X + Y /= 1.
    !
    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
    !    If STATUS is negative, then this is the value exceeded by parameter I.
    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ), Parameter :: atol = 1.0D-50
    Real ( kind = 8 ) b
    Real ( kind = 8 ) bound
    Real ( kind = 8 ) ccum
    Real ( kind = 8 ) cum
    Real ( kind = 8 ) fx
    Real ( kind = 8 ), Parameter :: inf = 1.0D+300
    Real ( kind = 8 ) p
    Real ( kind = 8 ) q
    Logical qhi
    Logical qleft
    Integer status
    Real ( kind = 8 ), Parameter :: tol = 1.0D-08
    Integer which
    Real ( kind = 8 ) x
    Real ( kind = 8 ) xhi
    Real ( kind = 8 ) xlo
    Real ( kind = 8 ) y

    If ( which < 1 ) Then
       bound = 1.0D+00
       status = -1
       Write ( *, '(a)' ) ' '
       Write ( *, '(a)' ) 'CDFBET - Fatal error!'
       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
       Return
    End If

    If ( 4 < which ) Then
       bound = 4.0D+00
       status = -1
       Write ( *, '(a)' ) ' '
       Write ( *, '(a)' ) 'CDFBET - Fatal error!'
       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
       Return
    End If
    !
    !  Unless P is to be computed, make sure it is legal.
    !
    If ( which /= 1 ) Then
       If ( p < 0.0D+00 ) Then
          bound = 0.0D+00
          status = -2
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter P is out of range.'
          Return
       Else If ( 1.0D+00 < p ) Then
          bound = 1.0D+00
          status = -2
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter P is out of range.'
          Return
       End If
    End If
    !
    !  Unless Q is to be computed, make sure it is legal.
    !
    If ( which /= 1 ) Then
       If ( q < 0.0D+00 ) Then
          bound = 0.0D+00
          status = -3
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
          Return
       Else If ( 1.0D+00 < q ) Then
          bound = 1.0D+00
          status = -3
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
          Return
       End If
    End If
    !
    !  Unless X is to be computed, make sure it is legal.
    !
    If ( which /= 2 ) Then
       If ( x < 0.0D+00 ) Then
          bound = 0.0D+00
          status = -4
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter X is out of range.'
          Return
       Else If ( 1.0D+00 < x ) Then
          bound = 1.0D+00
          status = -4
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter X is out of range.'
          Return
       End If
    End If
    !
    !  Unless Y is to be computed, make sure it is legal.
    !
    If ( which /= 2 ) Then
       If ( y < 0.0D+00 ) Then
          bound = 0.0D+00
          status = -5
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter Y is out of range.'
          Return
       Else If ( 1.0D+00 < y ) Then
          bound = 1.0D+00
          status = -5
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter Y is out of range.'
          Return
       End If
    End If
    !
    !  Unless A is to be computed, make sure it is legal.
    !
    If ( which /= 3 ) Then
       If ( a <= 0.0D+00 ) Then
          bound = 0.0D+00
          status = -6
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter A is out of range.'
          Return
       End If
    End If
    !
    !  Unless B is to be computed, make sure it is legal.
    !
    If ( which /= 4 ) Then
       If ( b <= 0.0D+00 ) Then
          bound = 0.0D+00
          status = -7
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter B is out of range.'
          Return
       End If
    End If
    !
    !  Check that P + Q = 1.
    !
    If ( which /= 1 ) Then
       If ( 3.0D+00 * Epsilon ( p ) < Abs ( ( p + q ) - 1.0D+00 ) ) Then
          status = 3
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  P + Q /= 1.'
          Return
       End If
    End If
    !
    !  Check that X + Y = 1.
    !
    If ( which /= 2 ) Then
       If ( 3.0D+00 * Epsilon ( x ) < Abs ( ( x + y ) - 1.0D+00 ) ) Then
          status = 4
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
          Write ( *, '(a)' ) '  X + Y /= 1.'
          Return
       End If
    End If
    !
    !  Compute P and Q.
    !
    If ( which == 1 ) Then

       Call cumbet ( x, y, a, b, p, q )
       status = 0
       !
       !  Compute X and Y.
       !
    Else If ( which == 2 ) Then

       ! Call dstzr ( 0.0D+00, 1.0D+00, atol, tol )
       Call dzror ( status, atol, tol, 0.0D+00, 1.0D+00, qleft, qhi, .TRUE. )

       If ( p <= q ) Then

          status = 0
          fx = 0.0D+00
          Call dzror ( status, x, fx, xlo, xhi, qleft, qhi, .FALSE. )
          y = 1.0D+00 - x

          Do

             If ( status /= 1 ) Then
                Exit
             End If

             Call cumbet ( x, y, a, b, cum, ccum )
             fx = cum - p
             Call dzror ( status, x, fx, xlo, xhi, qleft, qhi, .FALSE. )
             y = 1.0D+00 - x

          End Do

       Else

          status = 0
          fx = 0.0D+00
          Call dzror ( status, y, fx, xlo, xhi, qleft, qhi, .FALSE. )
          x = 1.0D+00 - y

          Do

             If ( status /= 1 ) Then
                Exit
             End If

             Call cumbet ( x, y, a, b, cum, ccum )
             fx = ccum - q
             Call dzror ( status, y, fx, xlo, xhi, qleft, qhi, .FALSE. )
             x = 1.0D+00 - y

          End Do

       End If

       If ( status == -1 ) Then
          If ( qleft ) Then
             status = 1
             bound = 0.0D+00
             Write ( *, '(a)' ) ' '
             Write ( *, '(a)' ) 'CDFBET - Warning!'
             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
             Write ( *, '(a)' ) '  the search bound.'
          Else
             status = 2
             bound = 1.0D+00
             Write ( *, '(a)' ) ' '
             Write ( *, '(a)' ) 'CDFBET - Warning!'
             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
             Write ( *, '(a)' ) '  the search bound.'
          End If
       End If
       !
       !  Compute A.
       !
    Else If ( which == 3 ) Then

       a = 5.0D+00
       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
       Call dinvr ( status, a, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
       status = 0
       Call dinvr ( status, a, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )

       Do

          If ( status /= 1 ) Then
             Exit
          End If

          Call cumbet ( x, y, a, b, cum, ccum )

          If ( p <= q ) Then
             fx = cum - p
          Else
             fx = ccum - q
          End If

          Call dinvr ( status, a, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )

       End Do

       If ( status == -1 ) Then

          If ( qleft ) Then
             status = 1
             bound = 0.0D+00
             Write ( *, '(a)' ) ' '
             Write ( *, '(a)' ) 'CDFBET - Warning!'
             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
             Write ( *, '(a)' ) '  the search bound.'
          Else
             status = 2
             bound = inf
             Write ( *, '(a)' ) ' '
             Write ( *, '(a)' ) 'CDFBET - Warning!'
             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
             Write ( *, '(a)' ) '  the search bound.'
          End If

       End If
       !
       !  Compute B.
       !
    Else If ( which == 4 ) Then

       b = 5.0D+00
       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
       Call dinvr ( status, b, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
       status = 0
       Call dinvr ( status, b, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )

       Do

          If ( status /= 1 ) Then
             Exit
          End If

          Call cumbet ( x, y, a, b, cum, ccum )

          If ( p <= q ) Then
             fx = cum - p
          Else
             fx = ccum - q
          End If

          Call dinvr ( status, b, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )

       End Do

       If ( status == -1 ) Then
          If ( qleft ) Then
             status = 1
             bound = 0.0D+00
             Write ( *, '(a)' ) ' '
             Write ( *, '(a)' ) 'CDFBET - Warning!'
             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
             Write ( *, '(a)' ) '  the search bound.'
          Else
             status = 2
             bound = inf
             Write ( *, '(a)' ) ' '
             Write ( *, '(a)' ) 'CDFBET - Warning!'
             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
             Write ( *, '(a)' ) '  the search bound.'
          End If
       End If

    End If

    Return
  End Subroutine cdfbet

!!$  Subroutine cdfbin ( which, p, q, s, xn, pr, ompr, status, bound )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CDFBIN evaluates the CDF of the Binomial distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine calculates any one parameter of the binomial distribution
!!$    !    given the others.
!!$    !
!!$    !    The value P of the cumulative distribution function is calculated
!!$    !    directly.
!!$    !
!!$    !    Computation of the other parameters involves a seach for a value that
!!$    !    produces the desired value of P.  The search relies on the
!!$    !    monotonicity of P with respect to the other parameters.
!!$    !
!!$    !    P is the probablility of S or fewer successes in XN binomial trials,
!!$    !    each trial having an individual probability of success of PR.
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    09 June 2004
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.5.24.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, integer WHICH, indicates which of argument values is to
!!$    !    be calculated from the others.
!!$    !    1: Calculate P and Q from S, XN, PR and OMPR;
!!$    !    2: Calculate S from P, Q, XN, PR and OMPR;
!!$    !    3: Calculate XN from P, Q, S, PR and OMPR;
!!$    !    4: Calculate PR and OMPR from P, Q, S and XN.
!!$    !
!!$    !    Input/output, real ( kind = 8 ) P, the cumulation, from 0 to S,
!!$    !    of the binomial distribution.  If P is an input value, it should
!!$    !    lie in the range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) Q, equal to 1-P.  If Q is an input
!!$    !    value, it should lie in the range [0,1].  If Q is an output value,
!!$    !    it will lie in the range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) S, the number of successes observed.
!!$    !    Whether this is an input or output value, it should lie in the
!!$    !    range [0,XN].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) XN, the number of binomial trials.
!!$    !    If this is an input value it should lie in the range: (0, +infinity).
!!$    !    If it is an output value it will be searched for in the
!!$    !    range [1.0D-300, 1.0D+300].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) PR, the probability of success in each
!!$    !    binomial trial.  Whether this is an input or output value, it should
!!$    !    lie in the range: [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) OMPR, equal to 1-PR.  Whether this is an
!!$    !    input or output value, it should lie in the range [0,1].  Also, it should
!!$    !    be the case that PR + OMPR = 1.
!!$    !
!!$    !    Output, integer STATUS, reports the status of the computation.
!!$    !     0, if the calculation completed correctly;
!!$    !    -I, if the input parameter number I is out of range;
!!$    !    +1, if the answer appears to be lower than lowest search bound;
!!$    !    +2, if the answer appears to be higher than greatest search bound;
!!$    !    +3, if P + Q /= 1;
!!$    !    +4, if PR + OMPR /= 1.
!!$    !
!!$    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
!!$    !    If STATUS is negative, then this is the value exceeded by parameter I.
!!$    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ), Parameter :: atol = 1.0D-50
!!$    Real ( kind = 8 ) bound
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Parameter :: inf = 1.0D+300
!!$    Real ( kind = 8 ) ompr
!!$    Real ( kind = 8 ) p
!!$    Real ( kind = 8 ) pr
!!$    Real ( kind = 8 ) q
!!$    Logical qhi
!!$    Logical qleft
!!$    Real ( kind = 8 ) s
!!$    Integer status
!!$    Real ( kind = 8 ), Parameter :: tol = 1.0D-08
!!$    Integer which
!!$    Real ( kind = 8 ) xhi
!!$    Real ( kind = 8 ) xlo
!!$    Real ( kind = 8 ) xn
!!$    !
!!$    !  Check the arguments.
!!$    !
!!$    If ( which < 1 ) Then
!!$       bound = 1.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$
!!$    If ( 4 < which ) Then
!!$       bound = 4.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$    !
!!$    !  Unless P is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( p < 0.0D+00 ) Then
!!$          status = -2
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < p ) Then
!!$          status = -2
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless Q is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( q < 0.0D+00 ) Then
!!$          status = -3
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < q ) Then
!!$          status = -3
!!$          bound = 1.0
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless XN is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 3 ) Then
!!$       If ( xn <= 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -5
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter XN is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless S is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 2 ) Then
!!$       If ( s < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -4
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter S is out of range.'
!!$          Return
!!$       Else If ( which /= 3 .And. xn < s ) Then
!!$          bound = xn
!!$          status = -4
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter S is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless PR is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 4 ) Then
!!$       If ( pr < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -6
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter PR is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < pr ) Then
!!$          bound = 1.0D+00
!!$          status = -6
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter PR is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless OMPR is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 4 ) Then
!!$       If ( ompr < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -7
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter OMPR is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < ompr ) Then
!!$          bound = 1.0D+00
!!$          status = -7
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter OMPR is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Check that P + Q = 1.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( 3.0D+00 * Epsilon ( p ) < Abs ( ( p + q ) - 1.0D+00 ) ) Then
!!$          status = 3
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  P + Q /= 1.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Check that PR + OMPR = 1.
!!$    !
!!$    If ( which /= 4 ) Then
!!$       If ( 3.0D+00 * Epsilon ( 1.0D+00 ) &
!!$            < Abs ( ( pr + ompr ) - 1.0D+00 ) ) Then
!!$          status = 4
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFBET - Fatal error!'
!!$          Write ( *, '(a)' ) '  PR + OMPR /= 1.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Calculate P and Q.
!!$    !
!!$    If ( which == 1 ) Then
!!$
!!$       Call cumbin ( s, xn, pr, ompr, p, q )
!!$       status = 0
!!$       !
!!$       !  Calculate S.
!!$       !
!!$    Else If ( which == 2 ) Then
!!$
!!$       s = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, xn, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, s, fx, qleft, qhi, 0.0D+00, xn, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       fx = 0.0D+00
!!$
!!$       Call dinvr ( status, s, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumbin ( s, xn, pr, ompr, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, s, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFBIN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = xn
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFBIN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate XN.
!!$       !
!!$    Else If ( which == 3 ) Then
!!$
!!$       xn = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, xn, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, xn, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumbin ( s, xn, pr, ompr, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, xn, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFBIN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$             Return
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFBIN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$             Return
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate PR and OMPR.
!!$       !
!!$    Else If ( which == 4 ) Then
!!$
!!$       ! Call dstzr ( 0.0D+00, 1.0D+00, atol, tol )
!!$       Call dzror ( status, atol, tol, 0.0D+00, 1.0D+00, qleft, qhi, .TRUE. )
!!$
!!$       If ( p <= q ) Then
!!$
!!$          status = 0
!!$          Call dzror ( status, pr, fx, xlo, xhi, qleft, qhi, .FALSE. )
!!$          ompr = 1.0D+00 - pr
!!$
!!$          Do
!!$
!!$             If ( status /= 1 ) Then
!!$                Exit
!!$             End If
!!$
!!$             Call cumbin ( s, xn, pr, ompr, cum, ccum )
!!$             fx = cum - p
!!$             Call dzror ( status, pr, fx, xlo, xhi, qleft, qhi, .FALSE. )
!!$             ompr = 1.0D+00 - pr
!!$
!!$          End Do
!!$
!!$       Else
!!$
!!$          status = 0
!!$          Call dzror ( status, ompr, fx, xlo, xhi, qleft, qhi, .FALSE. )
!!$          pr = 1.0D+00 - ompr
!!$
!!$          Do
!!$
!!$             If ( status /= 1 ) Then
!!$                Exit
!!$             End If
!!$
!!$             Call cumbin ( s, xn, pr, ompr, cum, ccum )
!!$             fx = ccum - q
!!$             Call dzror ( status, ompr, fx, xlo, xhi, qleft, qhi, .FALSE. )
!!$             pr = 1.0D+00 - ompr
!!$
!!$          End Do
!!$
!!$       End If
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFBIN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = 1.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFBIN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$
!!$    End If
!!$
!!$    Return
!!$  End Subroutine cdfbin

!!$  Subroutine cdfchi ( which, p, q, x, df, status, bound )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CDFCHI evaluates the CDF of the chi square distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine calculates any one parameter of the chi square distribution
!!$    !    given the others.
!!$    !
!!$    !    The value P of the cumulative distribution function is calculated
!!$    !    directly.
!!$    !
!!$    !    Computation of the other parameters involves a seach for a value that
!!$    !    produces the desired value of P.  The search relies on the
!!$    !    monotonicity of P with respect to the other parameters.
!!$    !
!!$    !    The CDF of the chi square distribution can be evaluated
!!$    !    within Mathematica by commands such as:
!!$    !
!!$    !      Needs["Statistics`ContinuousDistributions`"]
!!$    !      CDF [ ChiSquareDistribution [ DF ], X ]
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.4.19.
!!$    !
!!$    !    Stephen Wolfram,
!!$    !    The Mathematica Book,
!!$    !    Fourth Edition,
!!$    !    Wolfram Media / Cambridge University Press, 1999.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, integer WHICH, indicates which argument is to be calculated
!!$    !    from the others.
!!$    !    1: Calculate P and Q from X and DF;
!!$    !    2: Calculate X from P, Q and DF;
!!$    !    3: Calculate DF from P, Q and X.
!!$    !
!!$    !    Input/output, real ( kind = 8 ) P, the integral from 0 to X of
!!$    !    the chi-square distribution.  If this is an input value, it should
!!$    !    lie in the range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) Q, equal to 1-P.  If Q is an input
!!$    !    value, it should lie in the range [0,1].  If Q is an output value,
!!$    !    it will lie in the range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) X, the upper limit of integration
!!$    !    of the chi-square distribution.  If this is an input
!!$    !    value, it should lie in the range: [0, +infinity).  If it is an output
!!$    !    value, it will be searched for in the range: [0,1.0D+300].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) DF, the degrees of freedom of the
!!$    !    chi-square distribution.  If this is an input value, it should lie
!!$    !    in the range: (0, +infinity).  If it is an output value, it will be
!!$    !    searched for in the range: [ 1.0D-300, 1.0D+300].
!!$    !
!!$    !    Output, integer STATUS, reports the status of the computation.
!!$    !     0, if the calculation completed correctly;
!!$    !    -I, if the input parameter number I is out of range;
!!$    !    +1, if the answer appears to be lower than lowest search bound;
!!$    !    +2, if the answer appears to be higher than greatest search bound;
!!$    !    +3, if P + Q /= 1;
!!$    !    +10, an error was returned from CUMGAM.
!!$    !
!!$    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
!!$    !    If STATUS is negative, then this is the value exceeded by parameter I.
!!$    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ), Parameter :: atol = 1.0D-50
!!$    Real ( kind = 8 ) bound
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) df
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Parameter :: inf = 1.0D+300
!!$    Real ( kind = 8 ) p
!!$    Real ( kind = 8 ) porq
!!$    Real ( kind = 8 ) q
!!$    Logical qhi
!!$    Logical qleft
!!$    Integer status
!!$    Real ( kind = 8 ), Parameter :: tol = 1.0D-08
!!$    Integer which
!!$    Real ( kind = 8 ) x
!!$    !
!!$    !  Check the arguments.
!!$    !
!!$    If ( which < 1 ) Then
!!$       bound = 1.0
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFCHI - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$
!!$    If ( 3 < which ) Then
!!$       bound = 3.0
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFCHI - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$    !
!!$    !  Unless P is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( p < 0.0D+00 ) Then
!!$          status = -2
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < p ) Then
!!$          status = -2
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless Q is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( q < 0.0D+00 ) Then
!!$          status = -3
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < q ) Then
!!$          status = -3
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless X is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 2 ) Then
!!$       If ( x < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -4
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter X is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless DF is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 3 ) Then
!!$       If ( df <= 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -5
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter DF is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Check that P + Q = 1.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( 3.0D+00 * Epsilon ( p ) < Abs ( ( p + q ) - 1.0D+00 ) ) Then
!!$          status = 3
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHI - Fatal error!'
!!$          Write ( *, '(a)' ) '  P + Q /= 1.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Select the minimum of P and Q.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       porq = Min ( p, q )
!!$    End If
!!$    !
!!$    !  Calculate P and Q.
!!$    !
!!$    If ( which == 1 ) Then
!!$
!!$       status = 0
!!$       Call cumchi ( x, df, p, q )
!!$
!!$       If ( 1.5D+00 < porq ) Then
!!$          status = 10
!!$          Return
!!$       End If
!!$       !
!!$       !  Calculate X.
!!$       !
!!$    Else If ( which == 2 ) Then
!!$
!!$       x = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, x, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       fx = 0.0D+00
!!$       Call dinvr ( status, x, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumchi ( x, df, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          If ( 1.5D+00 < fx + porq ) Then
!!$             status = 10
!!$             Return
!!$          End If
!!$
!!$          Call dinvr ( status, x, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFCHI - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFCHI - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$
!!$       End If
!!$       !
!!$       !  Calculate DF.
!!$       !
!!$    Else If ( which == 3 ) Then
!!$
!!$       df = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, df, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, df, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumchi ( x, df, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          If ( 1.5D+00 < fx + porq ) Then
!!$             status = 10
!!$             Return
!!$          End If
!!$
!!$          Call dinvr ( status, df, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFCHI - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFCHI - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$
!!$    End If
!!$
!!$    Return
!!$  End Subroutine cdfchi

!!$  Subroutine cdfchn ( which, p, q, x, df, pnonc, status, bound )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CDFCHN evaluates the CDF of the Noncentral Chi-Square.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine calculates any one parameter of the noncentral chi-square
!!$    !    distribution given values for the others.
!!$    !
!!$    !    The value P of the cumulative distribution function is calculated
!!$    !    directly.
!!$    !
!!$    !    Computation of the other parameters involves a seach for a value that
!!$    !    produces the desired value of P.  The search relies on the
!!$    !    monotonicity of P with respect to the other parameters.
!!$    !
!!$    !    The computation time required for this routine is proportional
!!$    !    to the noncentrality parameter (PNONC).  Very large values of
!!$    !    this parameter can consume immense computer resources.  This is
!!$    !    why the search range is bounded by 10,000.
!!$    !
!!$    !    The CDF of the noncentral chi square distribution can be evaluated
!!$    !    within Mathematica by commands such as:
!!$    !
!!$    !      Needs["Statistics`ContinuousDistributions`"]
!!$    !      CDF[ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.5.25.
!!$    !
!!$    !    Stephen Wolfram,
!!$    !    The Mathematica Book,
!!$    !    Fourth Edition,
!!$    !    Wolfram Media / Cambridge University Press, 1999.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, integer WHICH, indicates which argument is to be calculated
!!$    !    from the others.
!!$    !    1: Calculate P and Q from X, DF and PNONC;
!!$    !    2: Calculate X from P, DF and PNONC;
!!$    !    3: Calculate DF from P, X and PNONC;
!!$    !    4: Calculate PNONC from P, X and DF.
!!$    !
!!$    !    Input/output, real ( kind = 8 ) P, the integral from 0 to X of
!!$    !    the noncentral chi-square distribution.  If this is an input
!!$    !    value, it should lie in the range: [0, 1.0-1.0D-16).
!!$    !
!!$    !    Input/output, real ( kind = 8 ) Q, is generally not used by this
!!$    !    subroutine and is only included for similarity with other routines.
!!$    !    However, if P is to be computed, then a value will also be computed
!!$    !    for Q.
!!$    !
!!$    !    Input, real ( kind = 8 ) X, the upper limit of integration of the
!!$    !    noncentral chi-square distribution.  If this is an input value, it
!!$    !    should lie in the range: [0, +infinity).  If it is an output value,
!!$    !    it will be sought in the range: [0,1.0D+300].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) DF, the number of degrees of freedom
!!$    !    of the noncentral chi-square distribution.  If this is an input value,
!!$    !    it should lie in the range: (0, +infinity).  If it is an output value,
!!$    !    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) PNONC, the noncentrality parameter of
!!$    !    the noncentral chi-square distribution.  If this is an input value, it
!!$    !    should lie in the range: [0, +infinity).  If it is an output value,
!!$    !    it will be searched for in the range: [0,1.0D+4]
!!$    !
!!$    !    Output, integer STATUS, reports on the calculation.
!!$    !    0, if calculation completed correctly;
!!$    !    -I, if input parameter number I is out of range;
!!$    !    1, if the answer appears to be lower than the lowest search bound;
!!$    !    2, if the answer appears to be higher than the greatest search bound.
!!$    !
!!$    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
!!$    !    If STATUS is negative, then this is the value exceeded by parameter I.
!!$    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ), Parameter :: atol=1.0D-50
!!$    Real ( kind = 8 ) bound
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) df
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Parameter :: inf=1.0D+300
!!$    Real ( kind = 8 ) p
!!$    Real ( kind = 8 ) pnonc
!!$    Real ( kind = 8 ) q
!!$    Logical qhi
!!$    Logical qleft
!!$    Integer status
!!$    Real ( kind = 8 ), Parameter :: tent4=1.0D+04
!!$    Real ( kind = 8 ), Parameter :: tol=1.0D-08
!!$    Integer which
!!$    Real ( kind = 8 ) x
!!$    !
!!$    !  Check the arguments.
!!$    !
!!$    If ( which < 1 ) Then
!!$       bound = 1.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFCHN - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$
!!$    If ( 4 < which ) Then
!!$       bound = 4.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFCHN - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$    !
!!$    !  Unless P is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( p < 0.0D+00 ) Then
!!$          status = -2
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < p ) Then
!!$          status = -2
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless X is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 2 ) Then
!!$       If ( x < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -4
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter X is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless DF is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 3 ) Then
!!$       If ( df <= 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -5
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter DF is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless PNONC is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 4 ) Then
!!$       If ( pnonc < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -6
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFCHN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter PNONC is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Calculate P and Q.
!!$    !
!!$    If ( which == 1 ) Then
!!$
!!$       Call cumchn ( x, df, pnonc, p, q )
!!$       status = 0
!!$       !
!!$       !  Calculate X.
!!$       !
!!$    Else If ( which == 2 ) Then
!!$
!!$       x = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, x, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       fx = 0.0D+00
!!$       Call dinvr ( status, x, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$          Call cumchn ( x, df, pnonc, cum, ccum )
!!$          fx = cum - p
!!$          Call dinvr ( status, x, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFCHN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFCHN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate DF.
!!$       !
!!$    Else If ( which == 3 ) Then
!!$
!!$       df = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, df, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, df, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$          Call cumchn ( x, df, pnonc, cum, ccum )
!!$          fx = cum - p
!!$          Call dinvr ( status, df, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$       End Do
!!$
!!$       If ( status == -1 )Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFCHN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFCHN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate PNONC.
!!$       !
!!$    Else If ( which == 4 ) Then
!!$
!!$       pnonc = 5.0D+00
!!$       !Call dstinv ( 0.0D+00, tent4, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, pnonc, fx, qleft, qhi, 0.0D+00, tent4, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, pnonc, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumchn ( x, df, pnonc, cum, ccum )
!!$          fx = cum - p
!!$          Call dinvr ( status, pnonc, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFCHN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = tent4
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFCHN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$
!!$    End If
!!$
!!$    Return
!!$  End Subroutine cdfchn

!!$  Subroutine cdff ( which, p, q, f, dfn, dfd, status, bound )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CDFF evaluates the CDF of the F distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine calculates any one parameter of the F distribution
!!$    !    given the others.
!!$    !
!!$    !    The value P of the cumulative distribution function is calculated
!!$    !    directly.
!!$    !
!!$    !    Computation of the other parameters involves a seach for a value that
!!$    !    produces the desired value of P.  The search relies on the
!!$    !    monotonicity of P with respect to the other parameters.
!!$    !
!!$    !    The value of the cumulative F distribution is not necessarily
!!$    !    monotone in either degrees of freedom.  There thus may be two
!!$    !    values that provide a given CDF value.  This routine assumes
!!$    !    monotonicity and will find an arbitrary one of the two values.
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.6.2.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, integer WHICH, indicates which argument is to be calculated
!!$    !    from the others.
!!$    !    1: Calculate P and Q from F, DFN and DFD;
!!$    !    2: Calculate F from P, Q, DFN and DFD;
!!$    !    3: Calculate DFN from P, Q, F and DFD;
!!$    !    4: Calculate DFD from P, Q, F and DFN.
!!$    !
!!$    !    Input/output, real ( kind = 8 ) P, the integral from 0 to F of
!!$    !    the F-density.  If it is an input value, it should lie in the
!!$    !    range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) Q, equal to 1-P.  If Q is an input
!!$    !    value, it should lie in the range [0,1].  If Q is an output value,
!!$    !    it will lie in the range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) F, the upper limit of integration
!!$    !    of the F-density.  If this is an input value, it should lie in the
!!$    !    range [0, +infinity).  If it is an output value, it will be searched
!!$    !    for in the range [0,1.0D+300].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) DFN, the number of degrees of
!!$    !    freedom of the numerator sum of squares.  If this is an input value,
!!$    !    it should lie in the range: (0, +infinity).  If it is an output value,
!!$    !    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) DFD, the number of degrees of freedom
!!$    !    of the denominator sum of squares.  If this is an input value, it should
!!$    !    lie in the range: (0, +infinity).  If it is an output value, it will
!!$    !    be searched for in the  range: [ 1.0D-300, 1.0D+300].
!!$    !
!!$    !    Output, integer STATUS, reports the status of the computation.
!!$    !     0, if the calculation completed correctly;
!!$    !    -I, if the input parameter number I is out of range;
!!$    !    +1, if the answer appears to be lower than lowest search bound;
!!$    !    +2, if the answer appears to be higher than greatest search bound;
!!$    !    +3, if P + Q /= 1.
!!$    !
!!$    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
!!$    !    If STATUS is negative, then this is the value exceeded by parameter I.
!!$    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ), Parameter :: atol = 1.0D-50
!!$    Real ( kind = 8 ) bound
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) dfd
!!$    Real ( kind = 8 ) dfn
!!$    Real ( kind = 8 ) f
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Parameter :: inf = 1.0D+300
!!$    Real ( kind = 8 ) p
!!$    Real ( kind = 8 ) q
!!$    Logical qhi
!!$    Logical qleft
!!$    Integer status
!!$    Real ( kind = 8 ), Parameter :: tol = 1.0D-08
!!$    Integer which
!!$    !
!!$    !  Check the arguments.
!!$    !
!!$    If ( which < 1 ) Then
!!$       bound = 1.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFF - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$
!!$    If ( 4 < which ) Then
!!$       bound = 4.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFF - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$    !
!!$    !  Unless P is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( p < 0.0D+00 ) Then
!!$          status = -2
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFF - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < p ) Then
!!$          status = -2
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFF - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless Q is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( q < 0.0D+00 ) Then
!!$          status = -3
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFF - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < q ) Then
!!$          status = -3
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFF - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless F is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 2 ) Then
!!$       If ( f < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -4
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFF - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter F is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless DFN is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 3 ) Then
!!$       If ( dfn <= 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -5
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFF - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter DFN is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless DFD is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 4 ) Then
!!$       If ( dfd <= 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -6
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFF - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter DFD is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Check that P + Q = 1.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( 3.0D+00 * Epsilon ( p ) < Abs ( ( p + q ) - 1.0D+00 ) ) Then
!!$          status = 3
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFF - Fatal error!'
!!$          Write ( *, '(a)' ) '  P + Q /= 1.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Calculate P and Q.
!!$    !
!!$    If ( which == 1 ) Then
!!$
!!$       Call cumf ( f, dfn, dfd, p, q )
!!$       status = 0
!!$       !
!!$       !  Calculate F.
!!$       !
!!$    Else If ( which == 2 ) Then
!!$
!!$       f = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, f, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       fx = 0.0D+00
!!$       Call dinvr ( status, f, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumf ( f, dfn, dfd, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, f, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFF - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFF - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate DFN.
!!$       !
!!$    Else If ( which == 3 ) Then
!!$
!!$       dfn = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, dfn, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, dfn, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumf ( f, dfn, dfd, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, dfn, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFF - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$             Return
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFF - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$             Return
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate DFD.
!!$       !
!!$    Else If ( which == 4 ) Then
!!$
!!$       dfd = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, dfd, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, dfd, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumf ( f, dfn, dfd, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, dfd, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFF - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFF - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$
!!$    End If
!!$
!!$    Return
!!$  End Subroutine cdff

!!$  Subroutine cdffnc ( which, p, q, f, dfn, dfd, pnonc, status, bound )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CDFFNC evaluates the CDF of the Noncentral F distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine originally used 1.0D+300 as the upper bound for the
!!$    !    interval in which many of the missing parameters are to be sought.
!!$    !    Since the underlying rootfinder routine needs to evaluate the
!!$    !    function at this point, it is no surprise that the program was
!!$    !    experiencing overflows.  A less extravagant upper bound
!!$    !    is being tried for now!
!!$    !
!!$    !    This routine calculates any one parameter of the Noncentral F distribution
!!$    !    given the others.
!!$    !
!!$    !    The value P of the cumulative distribution function is calculated
!!$    !    directly.
!!$    !
!!$    !    Computation of the other parameters involves a seach for a value that
!!$    !    produces the desired value of P.  The search relies on the
!!$    !    monotonicity of P with respect to the other parameters.
!!$    !
!!$    !    The computation time required for this routine is proportional
!!$    !    to the noncentrality parameter PNONC.  Very large values of
!!$    !    this parameter can consume immense computer resources.  This is
!!$    !    why the search range is bounded by 10,000.
!!$    !
!!$    !    The value of the cumulative noncentral F distribution is not
!!$    !    necessarily monotone in either degree of freedom.  There thus
!!$    !    may be two values that provide a given CDF value.  This routine
!!$    !    assumes monotonicity and will find an arbitrary one of the two
!!$    !    values.
!!$    !
!!$    !    The CDF of the noncentral F distribution can be evaluated
!!$    !    within Mathematica by commands such as:
!!$    !
!!$    !      Needs["Statistics`ContinuousDistributions`"]
!!$    !      CDF [ NoncentralFRatioDistribution [ DFN, DFD, PNONC ], X ]
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    15 June 2004
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.6.20.
!!$    !
!!$    !    Stephen Wolfram,
!!$    !    The Mathematica Book,
!!$    !    Fourth Edition,
!!$    !    Wolfram Media / Cambridge University Press, 1999.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, integer WHICH, indicates which argument is to be calculated
!!$    !    from the others.
!!$    !    1: Calculate P and Q from F, DFN, DFD and PNONC;
!!$    !    2: Calculate F from P, Q, DFN, DFD and PNONC;
!!$    !    3: Calculate DFN from P, Q, F, DFD and PNONC;
!!$    !    4: Calculate DFD from P, Q, F, DFN and PNONC;
!!$    !    5: Calculate PNONC from P, Q, F, DFN and DFD.
!!$    !
!!$    !    Input/output, real ( kind = 8 ) P, the integral from 0 to F of
!!$    !    the noncentral F-density.  If P is an input value it should
!!$    !    lie in the range [0,1) (Not including 1!).
!!$    !
!!$    !    Dummy, real ( kind = 8 ) Q, is not used by this subroutine,
!!$    !    and is only included for similarity with the other routines.
!!$    !    Its input value is not checked.  If P is to be computed, the
!!$    !    Q is set to 1 - P.
!!$    !
!!$    !    Input/output, real ( kind = 8 ) F, the upper limit of integration
!!$    !    of the noncentral F-density.  If this is an input value, it should
!!$    !    lie in the range: [0, +infinity).  If it is an output value, it
!!$    !    will be searched for in the range: [0,1.0D+30].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) DFN, the number of degrees of freedom
!!$    !    of the numerator sum of squares.  If this is an input value, it should
!!$    !    lie in the range: (0, +infinity).  If it is an output value, it will
!!$    !    be searched for in the range: [ 1.0, 1.0D+30].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) DFD, the number of degrees of freedom
!!$    !    of the denominator sum of squares.  If this is an input value, it should
!!$    !    be in range: (0, +infinity).  If it is an output value, it will be
!!$    !    searched for in the range [1.0, 1.0D+30].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) PNONC, the noncentrality parameter
!!$    !    If this is an input value, it should be nonnegative.
!!$    !    If it is an output value, it will be searched for in the range: [0,1.0D+4].
!!$    !
!!$    !    Output, integer STATUS, reports the status of the computation.
!!$    !     0, if the calculation completed correctly;
!!$    !    -I, if the input parameter number I is out of range;
!!$    !    +1, if the answer appears to be lower than lowest search bound;
!!$    !    +2, if the answer appears to be higher than greatest search bound;
!!$    !    +3, if P + Q /= 1.
!!$    !
!!$    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
!!$    !    If STATUS is negative, then this is the value exceeded by parameter I.
!!$    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ), Parameter :: atol = 1.0D-50
!!$    Real ( kind = 8 ) bound
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) dfd
!!$    Real ( kind = 8 ) dfn
!!$    Real ( kind = 8 ) f
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Parameter :: inf = 1.0D+30
!!$    Real ( kind = 8 ) p
!!$    Real ( kind = 8 ) pnonc
!!$    Real ( kind = 8 ) q
!!$    Logical qhi
!!$    Logical qleft
!!$    Integer status
!!$    Real ( kind = 8 ), Parameter :: tent4 = 1.0D+04
!!$    Real ( kind = 8 ), Parameter :: tol = 1.0D-08
!!$    Integer which
!!$    !
!!$    !  Check the arguments.
!!$    !
!!$    If ( which < 1 ) Then
!!$       bound = 1.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFFNC - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$
!!$    If ( 5 < which ) Then
!!$       bound = 5.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFFNC - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$    !
!!$    !  Unless P is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( p < 0.0D+00 ) Then
!!$          status = -2
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFFNC - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < p ) Then
!!$          status = -2
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFFNC - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless F is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 2 ) Then
!!$       If ( f < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -4
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFFNC - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter F is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless DFN is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 3 ) Then
!!$       If ( dfn <= 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -5
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFFNC - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter DFN is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless DFD is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 4 ) Then
!!$       If ( dfd <= 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -6
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFFNC - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter DFD is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless PNONC is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 5 ) Then
!!$       If ( pnonc < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -7
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFFNC - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter PNONC is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    ! Calculate P and Q.
!!$    !
!!$    If ( which == 1 ) Then
!!$
!!$       Call cumfnc ( f, dfn, dfd, pnonc, p, q )
!!$       status = 0
!!$       !
!!$       !  Calculate F.
!!$       !
!!$    Else If ( which == 2 ) Then
!!$
!!$       f = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, f, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       fx = 0.0D+00
!!$       Call dinvr ( status, f, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumfnc ( f, dfn, dfd, pnonc, cum, ccum )
!!$          fx = cum - p
!!$          Call dinvr ( status, f, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFFNC - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$             Return
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFFNC - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$             Return
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate DFN.
!!$       !
!!$    Else If ( which == 3 ) Then
!!$
!!$       dfn = 5.0D+00
!!$       ! Call dstinv ( 1.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, dfn, fx, qleft, qhi, 1.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, dfn, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumfnc ( f, dfn, dfd, pnonc, cum, ccum )
!!$          fx = cum - p
!!$
!!$          Call dinvr ( status, dfn, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFFNC - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFFNC - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate DFD.
!!$       !
!!$    Else If ( which == 4 ) Then
!!$
!!$       dfd = 5.0D+00
!!$       ! Call dstinv ( 1.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, dfd, fx, qleft, qhi, 1.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, dfd, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumfnc ( f, dfn, dfd, pnonc, cum, ccum )
!!$          fx = cum - p
!!$          Call dinvr ( status, dfd, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFFNC - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFFNC - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate PNONC.
!!$       !
!!$    Else If ( which == 5 ) Then
!!$
!!$       pnonc = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, tent4, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, pnonc, fx, qleft, qhi, 0.0D+00, tent4, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, pnonc, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumfnc ( f, dfn, dfd, pnonc, cum, ccum )
!!$          fx = cum - p
!!$
!!$          Call dinvr ( status, pnonc, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFFNC - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = tent4
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFFNC - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$
!!$    End If
!!$
!!$    Return
!!$  End Subroutine cdffnc

  Subroutine cdfgam ( which, p, q, x, shape, scale, status, bound )

    !*****************************************************************************80
    !
    !! CDFGAM evaluates the CDF of the Gamma Distribution.
    !
    !  Discussion:
    !
    !    This routine calculates any one parameter of the Gamma distribution
    !    given the others.
    !
    !    The cumulative distribution function P is calculated directly.
    !
    !    Computation of the other parameters involves a seach for a value that
    !    produces the desired value of P.  The search relies on the
    !    monotonicity of P with respect to the other parameters.
    !
    !    The gamma density is proportional to T**(SHAPE - 1) * EXP(- SCALE * T)
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 654:
    !    Computation of the incomplete gamma function ratios and their inverse,
    !    ACM Transactions on Mathematical Software,
    !    Volume 12, 1986, pages 377-393.
    !
    !  Parameters:
    !
    !    Input, integer WHICH, indicates which argument is to be calculated
    !    from the others.
    !    1: Calculate P and Q from X, SHAPE and SCALE;
    !    2: Calculate X from P, Q, SHAPE and SCALE;
    !    3: Calculate SHAPE from P, Q, X and SCALE;
    !    4: Calculate SCALE from P, Q, X and SHAPE.
    !
    !    Input/output, real ( kind = 8 ) P, the integral from 0 to X of the
    !    Gamma density.  If this is an input value, it should lie in the
    !    range: [0,1].
    !
    !    Input/output, real ( kind = 8 ) Q, equal to 1-P.  If Q is an input
    !    value, it should lie in the range [0,1].  If Q is an output value,
    !    it will lie in the range [0,1].
    !
    !    Input/output, real ( kind = 8 ) X, the upper limit of integration of
    !    the Gamma density.  If this is an input value, it should lie in the
    !    range: [0, +infinity).  If it is an output value, it will lie in
    !    the range: [0,1E300].
    !
    !    Input/output, real ( kind = 8 ) SHAPE, the shape parameter of the
    !    Gamma density.  If this is an input value, it should lie in the range:
    !    (0, +infinity).  If it is an output value, it will be searched for
    !    in the range: [1.0D-300,1.0D+300].
    !
    !    Input/output, real ( kind = 8 ) SCALE, the scale parameter of the
    !    Gamma density.  If this is an input value, it should lie in the range
    !    (0, +infinity).  If it is an output value, it will be searched for
    !    in the range: (1.0D-300,1.0D+300].
    !
    !    Output, integer STATUS, reports the status of the computation.
    !     0, if the calculation completed correctly;
    !    -I, if the input parameter number I is out of range;
    !    +1, if the answer appears to be lower than lowest search bound;
    !    +2, if the answer appears to be higher than greatest search bound;
    !    +3, if P + Q /= 1;
    !    +10, if the Gamma or inverse Gamma routine cannot compute the answer.
    !    This usually happens only for X and SHAPE very large (more than 1.0D+10.
    !
    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
    !    If STATUS is negative, then this is the value exceeded by parameter I.
    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
    !
    Implicit None

    Real ( kind = 8 ), Parameter :: atol = 1.0D-50
    Real ( kind = 8 ) bound
    Real ( kind = 8 ) ccum
    Real ( kind = 8 ) cum
    Real ( kind = 8 ) fx
    Integer ierr
    Real ( kind = 8 ), Parameter :: inf=1.0D+300
    Real ( kind = 8 ) p
    Real ( kind = 8 ) porq
    Real ( kind = 8 ) q
    Logical qhi
    Logical qleft
    Real ( kind = 8 ) scale
    Real ( kind = 8 ) shape
    Real ( kind = 8 ), Parameter :: tol = 1.0D-08
    Integer status,which
    Real ( kind = 8 ) x
    Real ( kind = 8 ) xscale
    Real ( kind = 8 ) xx
    !
    !  Check the arguments.
    !
    If ( which < 1 ) Then
       bound = 1.0D+00
       status = -1
       Write ( *, '(a)' ) ' '
       Write ( *, '(a)' ) 'CDFGAM - Fatal error!'
       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
       Return
    End If

    If ( 4 < which ) Then
       bound = 4.0D+00
       status = -1
       Write ( *, '(a)' ) ' '
       Write ( *, '(a)' ) 'CDFGAM - Fatal error!'
       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
       Return
    End If
    !
    !  Unless P is to be computed, make sure it is legal.
    !
    If ( which /= 1 ) Then
       If ( p < 0.0D+00 ) Then
          status = -2
          bound = 0.0D+00
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFGAM - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter P is out of range.'
          Return
       Else If ( 1.0D+00 < p ) Then
          status = -2
          bound = 1.0D+00
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFGAM - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter P is out of range.'
          Return
       End If
    End If
    !
    !  Unless Q is to be computed, make sure it is legal.
    !
    If ( which /= 1 ) Then
       If ( q < 0.0D+00 ) Then
          status = -3
          bound = 0.0D+00
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFGAM - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
          Return
       Else If ( 1.0D+00 < q ) Then
          status = -3
          bound = 1.0
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFGAM - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
          Return
       End If
    End If
    !
    !  Unless X is to be computed, make sure it is legal.
    !
    If ( which /= 2 ) Then
       If ( x < 0.0D+00 ) Then
          bound = 0.0D+00
          status = -4
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFGAM - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter X is out of range.'
          Return
       End If
    End If
    !
    !  Unless SHAPE is to be computed, make sure it is legal.
    !
    If ( which /= 3 ) Then
       If ( shape <= 0.0D+00 ) Then
          bound = 0.0D+00
          status = -5
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFGAM - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter SHAPE is out of range.'
          Return
       End If
    End If
    !
    !  Unless SCALE is to be computed, make sure it is legal.
    !
    If ( which /= 4 ) Then
       If ( scale <= 0.0D+00 ) Then
          bound = 0.0D+00
          status = -6
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFGAM - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter SCALE is out of range.'
          Return
       End If
    End If
    !
    !  Check that P + Q = 1.
    !
    If ( which /= 1 ) Then
       If ( 3.0D+00 * Epsilon ( p ) < Abs ( ( p + q ) - 1.0D+0 ) ) Then
          status = 3
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFGAM - Fatal error!'
          Write ( *, '(a)' ) '  P + Q /= 1.'
          Return
       End If
    End If
    !
    !  Select the minimum of P or Q.
    !
    If ( which /= 1 ) Then
       porq = Min ( p, q )
    End If
    !
    !  Calculate P and Q.
    !
    If ( which == 1 ) Then

       status = 0
       xscale = x * scale
       Call cumgam ( xscale, shape, p, q )

       If ( 1.5D+00 < porq ) Then
          status = 10
       End If
       !
       !  Calculate X.
       !
    Else If ( which == 2 ) Then

       Call gamma_inc_inv ( shape, xx, -1.0D+00, p, q, ierr )

       If ( ierr < 0.0D+00 ) Then
          status = 10
          Return
       End If

       x = xx / scale
       status = 0
       !
       !  Calculate SHAPE.
       !
    Else If ( which == 3 ) Then

       shape = 5.0D+00
       xscale = x * scale
       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
       Call dinvr ( status, shape, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
       status = 0
       fx = 0.0D+00
       Call dinvr ( status, shape, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )

       Do

          If ( status /= 1 ) Then
             Exit
          End If

          Call cumgam ( xscale, shape, cum, ccum )

          If ( p <= q ) Then
             fx = cum - p
          Else
             fx = ccum - q
          End If

          If ( p <= q .And. 1.5D+00 < cum ) Then
             status = 10
             Return
          Else If ( q < p .And. 1.5D+00 < ccum ) Then
             status = 10
             Return
          End If

          Call dinvr ( status, shape, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )

       End Do

       If ( status == -1 ) Then
          If ( qleft ) Then
             status = 1
             bound = 0.0D+00
             Write ( *, '(a)' ) ' '
             Write ( *, '(a)' ) 'CDFGAM - Warning!'
             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
             Write ( *, '(a)' ) '  the search bound.'
          Else
             status = 2
             bound = inf
             Write ( *, '(a)' ) ' '
             Write ( *, '(a)' ) 'CDFGAM - Warning!'
             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
             Write ( *, '(a)' ) '  the search bound.'
          End If
       End If
       !
       !  Calculate SCALE.
       !
    Else If ( which == 4 ) Then

       Call gamma_inc_inv ( shape, xx, -1.0D+00, p, q, ierr )

       If ( ierr < 0.0D+00 ) Then
          status = 10
       Else
          scale = xx / x
          status = 0
       End If

    End If

    Return
  End Subroutine cdfgam

!!$  Subroutine cdfnbn ( which, p, q, f, s, pr, ompr, status, bound )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CDFNBN evaluates the CDF of the Negative Binomial distribution
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine calculates any one parameter of the negative binomial
!!$    !    distribution given values for the others.
!!$    !
!!$    !    The cumulative negative binomial distribution returns the
!!$    !    probability that there will be F or fewer failures before the
!!$    !    S-th success in binomial trials each of which has probability of
!!$    !    success PR.
!!$    !
!!$    !    The individual term of the negative binomial is the probability of
!!$    !    F failures before S successes and is
!!$    !    Choose( F, S+F-1 ) * PR^(S) * (1-PR)^F
!!$    !
!!$    !    Computation of other parameters involve a seach for a value that
!!$    !    produces the desired value of P.  The search relies on the
!!$    !    monotonicity of P with respect to the other parameters.
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.5.26.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, integer WHICH, indicates which argument is to be calculated
!!$    !    from the others.
!!$    !    1: Calculate P and Q from F, S, PR and OMPR;
!!$    !    2: Calculate F from P, Q, S, PR and OMPR;
!!$    !    3: Calculate S from P, Q, F, PR and OMPR;
!!$    !    4: Calculate PR and OMPR from P, Q, F and S.
!!$    !
!!$    !    Input/output, real ( kind = 8 ) P, the cumulation from 0 to F of
!!$    !    the negative binomial distribution.  If P is an input value, it
!!$    !    should lie in the range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) Q, equal to 1-P.  If Q is an input
!!$    !    value, it should lie in the range [0,1].  If Q is an output value,
!!$    !    it will lie in the range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) F, the upper limit of cumulation of
!!$    !    the binomial distribution.  There are F or fewer failures before
!!$    !    the S-th success.  If this is an input value, it may lie in the
!!$    !    range [0,+infinity), and if it is an output value, it will be searched
!!$    !    for in the range [0,1.0D+300].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) S, the number of successes.
!!$    !    If this is an input value, it should lie in the range: [0, +infinity).
!!$    !    If it is an output value, it will be searched for in the range:
!!$    !    [0, 1.0D+300].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) PR, the probability of success in each
!!$    !    binomial trial.  Whether an input or output value, it should lie in the
!!$    !    range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) OMPR, the value of (1-PR).  Whether an
!!$    !    input or output value, it should lie in the range [0,1].
!!$    !
!!$    !    Output, integer STATUS, reports the status of the computation.
!!$    !     0, if the calculation completed correctly;
!!$    !    -I, if the input parameter number I is out of range;
!!$    !    +1, if the answer appears to be lower than lowest search bound;
!!$    !    +2, if the answer appears to be higher than greatest search bound;
!!$    !    +3, if P + Q /= 1;
!!$    !    +4, if PR + OMPR /= 1.
!!$    !
!!$    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
!!$    !    If STATUS is negative, then this is the value exceeded by parameter I.
!!$    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ), Parameter :: atol = 1.0D-50
!!$    Real ( kind = 8 ) bound
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) f
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Parameter :: inf = 1.0D+300
!!$    Real ( kind = 8 ) ompr
!!$    Real ( kind = 8 ) p
!!$    Real ( kind = 8 ) pr
!!$    Real ( kind = 8 ) q
!!$    Logical qhi
!!$    Logical qleft
!!$    Real ( kind = 8 ) s
!!$    Integer status
!!$    Real ( kind = 8 ), Parameter :: tol = 1.0D-08
!!$    Integer which
!!$    Real ( kind = 8 ) xhi
!!$    Real ( kind = 8 ) xlo
!!$    !
!!$    !  Check the arguments.
!!$    !
!!$    If ( which < 1 ) Then
!!$       bound = 1.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$
!!$    If ( 4 < which ) Then
!!$       bound = 4.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$    !
!!$    !  Unless P is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( p < 0.0D+00 ) Then
!!$          status = -2
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < p ) Then
!!$          status = -2
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless Q is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( q < 0.0D+00 ) Then
!!$          status = -3
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < q ) Then
!!$          status = -3
!!$          bound = 1.0
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless F is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 2 ) Then
!!$       If ( f < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -4
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter F is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless S is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 3 ) Then
!!$       If ( s < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -5
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter S is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless PR is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 4 ) Then
!!$       If ( pr < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -6
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter PR is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < pr ) Then
!!$          bound = 1.0D+00
!!$          status = -6
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter PR is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless OMPR is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 4 ) Then
!!$       If ( ompr < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -7
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter OMPR is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < ompr ) Then
!!$          bound = 1.0D+00
!!$          status = -7
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter OMPR is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Check that P + Q = 1.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( 3.0D+00 * Epsilon ( p ) < Abs ( ( p + q ) - 1.0D+00 ) ) Then
!!$          status = 3
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  P + Q /= 1.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Check that PR + OMPR = 1.
!!$    !
!!$    If ( which /= 4 ) Then
!!$       If ( 3.0D+00 * Epsilon ( pr ) < Abs ( ( pr + ompr ) - 1.0D+00 ) ) Then
!!$          status = 4
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFNBN - Fatal error!'
!!$          Write ( *, '(a)' ) '  PR + OMPR /= 1.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Calculate P and Q.
!!$    !
!!$    If ( which == 1 ) Then
!!$
!!$       Call cumnbn ( f, s, pr, ompr, p, q )
!!$       status = 0
!!$       !
!!$       !  Calculate F.
!!$       !
!!$    Else If ( which == 2 ) Then
!!$
!!$       f = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, f, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       fx = 0.0D+00
!!$       Call dinvr ( status, f, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumnbn ( f, s, pr, ompr, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, f, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFNBN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFNBN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate S.
!!$       !
!!$    Else If ( which == 3 ) Then
!!$
!!$       s = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, s, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, s, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumnbn ( f, s, pr, ompr, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, s, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFNBN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFNBn - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate PR and OMPR.
!!$       !
!!$    Else If ( which == 4 ) Then
!!$
!!$       ! Call dstzr ( 0.0D+00, 1.0D+00, atol, tol )
!!$       Call dzror ( status, atol, tol, 0.0D+00, 1.0D+00, qleft, qhi, .TRUE. )
!!$
!!$       If ( p <= q ) Then
!!$
!!$          status = 0
!!$          Call dzror ( status, pr, fx, xlo, xhi, qleft, qhi, .FALSE. )
!!$          ompr = 1.0D+00 - pr
!!$
!!$          Do
!!$             If ( status /= 1 ) Then
!!$                Exit
!!$             End If
!!$             Call cumnbn ( f, s, pr, ompr, cum, ccum )
!!$             fx = cum - p
!!$             Call dzror ( status, pr, fx, xlo, xhi, qleft, qhi, .FALSE. )
!!$             ompr = 1.0D+00 - pr
!!$          End Do
!!$
!!$       Else
!!$
!!$          status = 0
!!$          Call dzror ( status, ompr, fx, xlo, xhi, qleft, qhi, .FALSE. )
!!$          pr = 1.0D+00 - ompr
!!$
!!$          Do
!!$
!!$             If ( status /= 1 ) Then
!!$                Exit
!!$             End If
!!$
!!$             Call cumnbn ( f, s, pr, ompr, cum, ccum )
!!$             fx = ccum - q
!!$             Call dzror ( status, ompr, fx, xlo, xhi, qleft, qhi, .FALSE. )
!!$             pr = 1.0D+00 - ompr
!!$
!!$          End Do
!!$
!!$       End If
!!$
!!$       If ( status == -1 ) Then
!!$
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFNBN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = 1.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFNBN - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$
!!$       End If
!!$
!!$    End If
!!$
!!$    Return
!!$  End Subroutine cdfnbn

  Subroutine cdfnor ( which, p, q, x, mean, sd, status, bound )

    !*****************************************************************************80
    !
    !! CDFNOR evaluates the CDF of the Normal distribution.
    !
    !  Discussion:
    !
    !    A slightly modified version of ANORM from SPECFUN
    !    is used to calculate the cumulative standard normal distribution.
    !
    !    The rational functions from pages 90-95 of Kennedy and Gentle
    !    are used as starting values to a Newton iteration which
    !    compute the inverse standard normal.  Therefore no searches are
    !    necessary for any parameter.
    !
    !    For X < -15, the asymptotic expansion for the normal is used  as
    !    the starting value in finding the inverse standard normal.
    !
    !    The normal density is proportional to
    !    exp ( - 0.5D+00 * (( X - MEAN)/SD)**2)
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions
    !    1966, Formula 26.2.12.
    !
    !    William Cody,
    !    Algorithm 715:
    !    SPECFUN - A Portable FORTRAN Package of
    !    Special Function Routines and Test Drivers,
    !    ACM Transactions on Mathematical Software,
    !    Volume 19, Number 1, pages 22-32, 1993.
    !
    !    William Kennedy, James Gentle,
    !    Statistical Computing,
    !    Marcel Dekker, NY, 1980,
    !    QA276.4 K46
    !
    !  Parameters:
    !
    !    Input, integer WHICH, indicates which argument is to be calculated
    !    from the others.
    !    1: Calculate P and Q from X, MEAN and SD;
    !    2: Calculate X from P, Q, MEAN and SD;
    !    3: Calculate MEAN from P, Q, X and SD;
    !    4: Calculate SD from P, Q, X and MEAN.
    !
    !    Input/output, real ( kind = 8 ) P, the integral from -infinity to X
    !    of the Normal density.  If this is an input or output value, it will
    !    lie in the range [0,1].
    !
    !    Input/output, real ( kind = 8 ) Q, equal to 1-P.  If Q is an input
    !    value, it should lie in the range [0,1].  If Q is an output value,
    !    it will lie in the range [0,1].
    !
    !    Input/output, real ( kind = 8 ) X, the upper limit of integration of
    !    the Normal density.
    !
    !    Input/output, real ( kind = 8 ) MEAN, the mean of the Normal density.
    !
    !    Input/output, real ( kind = 8 ) SD, the standard deviation of the
    !    Normal density.  If this is an input value, it should lie in the
    !    range (0,+infinity).
    !
    !    Output, integer STATUS, the status of the calculation.
    !    0, if calculation completed correctly;
    !    -I, if input parameter number I is out of range;
    !    1, if answer appears to be lower than lowest search bound;
    !    2, if answer appears to be higher than greatest search bound;
    !    3, if P + Q /= 1.
    !
    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
    !    If STATUS is negative, then this is the value exceeded by parameter I.
    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
    !
    Implicit None

    Real ( kind = 8 ) bound
!    Real ( kind = 8 ) dinvnr
    Real ( kind = 8 ) mean
    Real ( kind = 8 ) p
    Real ( kind = 8 ) q
    Real ( kind = 8 ) sd
    Integer status
    Integer which
    Real ( kind = 8 ) x
    Real ( kind = 8 ) z
    !
    !  Check the arguments.
    !
    status = 0

    If ( which < 1 ) Then
       status = -1
       bound = 1.0D+00
       Write ( *, '(a)' ) ' '
       Write ( *, '(a)' ) 'CDFNOR - Fatal error!'
       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
       Return
    Else If ( 4 < which ) Then
       status = -1
       bound = 4.0D+00
       Write ( *, '(a)' ) ' '
       Write ( *, '(a)' ) 'CDFNOR - Fatal error!'
       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
       Return
    End If
    !
    !  Unless P is to be computed, make sure it is legal.
    !
    If ( which /= 1 ) Then
       If ( p < 0.0D+00 ) Then
          status = -2
          bound = 0.0D+00
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFNOR - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter P is out of range.'
          Return
       Else If ( 1.0D+00 < p ) Then
          status = -2
          bound = 1.0D+00
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFNOR - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter P is out of range.'
          Return
       End If
    End If
    !
    !  Unless Q is to be computed, make sure it is legal.
    !
    If ( which /= 1 ) Then
       If ( q < 0.0D+00 ) Then
          status = -3
          bound = 0.0D+00
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFNOR - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
          Return
       Else If ( 1.0D+00 < q ) Then
          status = -3
          bound = 1.0D+00
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFNOR - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
          Return
       End If
    End If
    !
    !  Check that P + Q = 1.
    !
    If ( which /= 1 ) Then
       If ( 3.0D+00 * Epsilon ( p ) < Abs ( ( p + q ) - 1.0D+00 ) ) Then
          status = 3
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFNOR - Fatal error!'
          Write ( *, '(a)' ) '  P + Q /= 1.'
          Return
       End If
    End If

    If ( which /= 4 ) Then
       If ( sd <= 0.0D+00 ) Then
          bound = 0.0D+00
          status = -6
          Write ( *, '(a)' ) ' '
          Write ( *, '(a)' ) 'CDFNOR - Fatal error!'
          Write ( *, '(a)' ) '  Input parameter SD is out of range.'
          Return
       End If
    End If
    !
    !  Calculate P and Q.
    !
    If ( which == 1 ) Then

       z = ( x - mean ) / sd
       Call cumnor ( z, p, q )
       !
       !  Calculate X.
       !
    Else If ( which == 2 ) Then

       z = dinvnr ( p, q )
       x = sd * z + mean
       !
       !  Calculate MEAN.
       !
    Else If ( which == 3 ) Then

       z = dinvnr ( p, q )
       mean = x - sd * z
       !
       !  Calculate SD.
       !
    Else If ( which == 4 ) Then

       z = dinvnr ( p, q )
       sd = ( x - mean ) / z

    End If

    Return
  End Subroutine cdfnor

!!$  Subroutine cdfpoi ( which, p, q, s, xlam, status, bound )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CDFPOI evaluates the CDF of the Poisson distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine calculates any one parameter of the Poisson distribution
!!$    !    given the others.
!!$    !
!!$    !    The value P of the cumulative distribution function is calculated
!!$    !    directly.
!!$    !
!!$    !    Computation of other parameters involve a seach for a value that
!!$    !    produces the desired value of P.  The search relies on the
!!$    !    monotonicity of P with respect to the other parameters.
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.4.21.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, integer WHICH, indicates which argument is to be calculated
!!$    !    from the others.
!!$    !    1: Calculate P and Q from S and XLAM;
!!$    !    2: Calculate A from P, Q and XLAM;
!!$    !    3: Calculate XLAM from P, Q and S.
!!$    !
!!$    !    Input/output, real ( kind = 8 ) P, the cumulation from 0 to S of the
!!$    !    Poisson density.  Whether this is an input or output value, it will
!!$    !    lie in the range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) Q, equal to 1-P.  If Q is an input
!!$    !    value, it should lie in the range [0,1].  If Q is an output value,
!!$    !    it will lie in the range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) S, the upper limit of cumulation of
!!$    !    the Poisson CDF.  If this is an input value, it should lie in
!!$    !    the range: [0, +infinity).  If it is an output value, it will be
!!$    !    searched for in the range: [0,1.0D+300].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) XLAM, the mean of the Poisson
!!$    !    distribution.  If this is an input value, it should lie in the range
!!$    !    [0, +infinity).  If it is an output value, it will be searched for
!!$    !    in the range: [0,1E300].
!!$    !
!!$    !    Output, integer STATUS, reports the status of the computation.
!!$    !     0, if the calculation completed correctly;
!!$    !    -I, if the input parameter number I is out of range;
!!$    !    +1, if the answer appears to be lower than lowest search bound;
!!$    !    +2, if the answer appears to be higher than greatest search bound;
!!$    !    +3, if P + Q /= 1.
!!$    !
!!$    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
!!$    !    If STATUS is negative, then this is the value exceeded by parameter I.
!!$    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ), Parameter :: atol = 1.0D-50
!!$    Real ( kind = 8 ) bound
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Parameter :: inf = 1.0D+300
!!$    Real ( kind = 8 ) p
!!$    Real ( kind = 8 ) q
!!$    Logical qhi
!!$    Logical qleft
!!$    Real ( kind = 8 ) s
!!$    Integer status
!!$    Real ( kind = 8 ), Parameter :: tol = 1.0D-08
!!$    Integer which
!!$    Real ( kind = 8 ) xlam
!!$    !
!!$    !  Check the arguments.
!!$    !
!!$    If ( which < 1 ) Then
!!$       bound = 1.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFPOI - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$
!!$    If ( 3 < which ) Then
!!$       bound = 3.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFPOI - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$    !
!!$    !  Unless P is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( p < 0.0D+00 ) Then
!!$          status = -2
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFPOI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < p ) Then
!!$          status = -2
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFPOI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless Q is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( q < 0.0D+00 ) Then
!!$          status = -3
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFPOI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < q ) Then
!!$          status = -3
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFPOI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless S is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 2 ) Then
!!$       If ( s < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -4
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFPOI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter S is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless XLAM is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 3 ) Then
!!$       If ( xlam < 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -5
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFPOI - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter XLAM is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Check that P + Q = 1.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( 3.0D+00 * Epsilon ( p ) < Abs ( ( p + q ) - 1.0D+00 ) ) Then
!!$          status = 3
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFPOI - Fatal error!'
!!$          Write ( *, '(a)' ) '  P + Q /= 1.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Calculate P and Q.
!!$    !
!!$    If ( which == 1 ) Then
!!$
!!$       Call cumpoi ( s, xlam, p, q )
!!$       status = 0
!!$       !
!!$       !  Calculate S.
!!$       !
!!$    Else If ( which == 2 ) Then
!!$
!!$       s = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, s, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       fx = 0.0D+00
!!$       Call dinvr ( status, s, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumpoi ( s, xlam, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, s, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFPOI - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFPOI - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate XLAM.
!!$       !
!!$    Else If ( which == 3 ) Then
!!$
!!$       xlam = 5.0D+00
!!$       ! Call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, xlam, fx, qleft, qhi, 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, xlam, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumpoi ( s, xlam, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, xlam, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFPOI - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFPOI - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$
!!$    End If
!!$
!!$    Return
!!$  End Subroutine cdfpoi

!!$  Subroutine cdft ( which, p, q, t, df, status, bound )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CDFT evaluates the CDF of the T distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine calculates any one parameter of the T distribution
!!$    !    given the others.
!!$    !
!!$    !    The value P of the cumulative distribution function is calculated
!!$    !    directly.
!!$    !
!!$    !    Computation of other parameters involve a seach for a value that
!!$    !    produces the desired value of P.   The search relies on the
!!$    !    monotonicity of P with respect to the other parameters.
!!$    !
!!$    !    The original version of this routine allowed the search interval
!!$    !    to extend from -1.0D+300 to +1.0D+300, which is fine until you
!!$    !    try to evaluate a function at such a point!
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.5.27.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, integer WHICH, indicates which argument is to be calculated
!!$    !    from the others.
!!$    !    1 : Calculate P and Q from T and DF;
!!$    !    2 : Calculate T from P, Q and DF;
!!$    !    3 : Calculate DF from P, Q and T.
!!$    !
!!$    !    Input/output, real ( kind = 8 ) P, the integral from -infinity to T of
!!$    !    the T-density.  Whether an input or output value, this will lie in the
!!$    !    range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) Q, equal to 1-P.  If Q is an input
!!$    !    value, it should lie in the range [0,1].  If Q is an output value,
!!$    !    it will lie in the range [0,1].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) T, the upper limit of integration of
!!$    !    the T-density.  If this is an input value, it may have any value.
!!$    !    It it is an output value, it will be searched for in the range
!!$    !    [ -1.0D+30, 1.0D+30 ].
!!$    !
!!$    !    Input/output, real ( kind = 8 ) DF, the number of degrees of freedom
!!$    !    of the T distribution.  If this is an input value, it should lie
!!$    !    in the range: (0 , +infinity).  If it is an output value, it will be
!!$    !    searched for in the range: [1, 1.0D+10].
!!$    !
!!$    !    Output, integer STATUS, reports the status of the computation.
!!$    !     0, if the calculation completed correctly;
!!$    !    -I, if the input parameter number I is out of range;
!!$    !    +1, if the answer appears to be lower than lowest search bound;
!!$    !    +2, if the answer appears to be higher than greatest search bound;
!!$    !    +3, if P + Q /= 1.
!!$    !
!!$    !    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
!!$    !    If STATUS is negative, then this is the value exceeded by parameter I.
!!$    !    if STATUS is 1 or 2, this is the search bound that was exceeded.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ), Parameter :: atol = 1.0D-50
!!$    Real ( kind = 8 ) bound
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) df
!!$    !    Real ( kind = 8 ) dt1
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Parameter :: inf = 1.0D+30
!!$    Real ( kind = 8 ), Parameter :: maxdf = 1.0D+10
!!$    Real ( kind = 8 ) p
!!$    Real ( kind = 8 ) q
!!$    Logical qhi
!!$    Logical qleft
!!$    Integer status
!!$    Real ( kind = 8 ) t
!!$    Real ( kind = 8 ), Parameter :: tol = 1.0D-08
!!$    Integer which
!!$    !
!!$    !  Check the arguments.
!!$    !
!!$    If ( which < 1 ) Then
!!$       bound = 1.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFT - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$
!!$    If ( 3 < which ) Then
!!$       bound = 3.0D+00
!!$       status = -1
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CDFT - Fatal error!'
!!$       Write ( *, '(a)' ) '  Input parameter WHICH is out of range.'
!!$       Return
!!$    End If
!!$    !
!!$    !  Unless P is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( p < 0.0D+00 ) Then
!!$          status = -2
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFT - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < p ) Then
!!$          status = -2
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFT - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter P is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless Q is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( q < 0.0D+00 ) Then
!!$          status = -3
!!$          bound = 0.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFT - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       Else If ( 1.0D+00 < q ) Then
!!$          status = -3
!!$          bound = 1.0D+00
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFT - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter Q is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Unless DF is to be computed, make sure it is legal.
!!$    !
!!$    If ( which /= 3 ) Then
!!$       If ( df <= 0.0D+00 ) Then
!!$          bound = 0.0D+00
!!$          status = -5
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFT - Fatal error!'
!!$          Write ( *, '(a)' ) '  Input parameter DF is out of range.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Check that P + Q = 1.
!!$    !
!!$    If ( which /= 1 ) Then
!!$       If ( 3.0D+00 * Epsilon ( 1.0D+00 ) &
!!$            < Abs ( ( p + q ) - 1.0D+00 ) ) Then
!!$          status = 3
!!$          Write ( *, '(a)' ) ' '
!!$          Write ( *, '(a)' ) 'CDFT - Fatal error!'
!!$          Write ( *, '(a)' ) '  P + Q /= 1.'
!!$          Return
!!$       End If
!!$    End If
!!$    !
!!$    !  Calculate P and Q.
!!$    !
!!$    If ( which == 1 ) Then
!!$
!!$       Call cumt ( t, df, p, q )
!!$       status = 0
!!$       !
!!$       !  Calculate T.
!!$       !
!!$    Else If ( which == 2 ) Then
!!$
!!$       t = dt1 ( p, q, df )
!!$       ! Call dstinv ( -inf, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, t, fx, qleft, qhi, -inf, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       fx = 0.0D+00
!!$       Call dinvr ( status, t, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumt ( t, df, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, t, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft )Then
!!$             status = 1
!!$             bound = -inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFT - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = inf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFT - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$       !
!!$       !  Calculate DF.
!!$       !
!!$    Else If ( which == 3 ) Then
!!$
!!$       df = 5.0D+00
!!$       ! Call dstinv ( 1.0D+00, maxdf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
!!$       Call dinvr ( status, df, fx, qleft, qhi, 1.0D+00, maxdf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol, .TRUE. )
!!$       status = 0
!!$       Call dinvr ( status, df, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       Do
!!$
!!$          If ( status /= 1 ) Then
!!$             Exit
!!$          End If
!!$
!!$          Call cumt ( t, df, cum, ccum )
!!$
!!$          If ( p <= q ) Then
!!$             fx = cum - p
!!$          Else
!!$             fx = ccum - q
!!$          End If
!!$
!!$          Call dinvr ( status, df, fx, qleft, qhi, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, .FALSE. )
!!$
!!$       End Do
!!$
!!$       If ( status == -1 ) Then
!!$          If ( qleft ) Then
!!$             status = 1
!!$             bound = 0.0D+00
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFT - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be lower than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          Else
!!$             status = 2
!!$             bound = maxdf
!!$             Write ( *, '(a)' ) ' '
!!$             Write ( *, '(a)' ) 'CDFT - Warning!'
!!$             Write ( *, '(a)' ) '  The desired answer appears to be higher than'
!!$             Write ( *, '(a)' ) '  the search bound.'
!!$          End If
!!$       End If
!!$
!!$    End If
!!$
!!$    Return
!!$  End Subroutine cdft

!!$  Subroutine chi_noncentral_cdf_values ( n_data, x, lambda, df, cdf )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CHI_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    The CDF of the noncentral chi square distribution can be evaluated
!!$    !    within Mathematica by commands such as:
!!$    !
!!$    !      Needs["Statistics`ContinuousDistributions`"]
!!$    !      CDF [ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    12 June 2004
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Stephen Wolfram,
!!$    !    The Mathematica Book,
!!$    !    Fourth Edition,
!!$    !    Wolfram Media / Cambridge University Press, 1999.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, real ( kind = 8 ) X, the argument of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) LAMBDA, the noncentrality parameter.
!!$    !
!!$    !    Output, integer DF, the number of degrees of freedom.
!!$    !
!!$    !    Output, real ( kind = 8 ) CDF, the noncentral chi CDF.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 27
!!$
!!$    Real ( kind = 8 ) cdf
!!$    Real, Save, Dimension ( n_max ) :: cdf_vec = (/ &
!!$         0.839944E+00, 0.695906E+00, 0.535088E+00, &
!!$         0.764784E+00, 0.620644E+00, 0.469167E+00, &
!!$         0.307088E+00, 0.220382E+00, 0.150025E+00, &
!!$         0.307116E-02, 0.176398E-02, 0.981679E-03, &
!!$         0.165175E-01, 0.202342E-03, 0.498448E-06, &
!!$         0.151325E-01, 0.209041E-02, 0.246502E-03, &
!!$         0.263684E-01, 0.185798E-01, 0.130574E-01, &
!!$         0.583804E-01, 0.424978E-01, 0.308214E-01, &
!!$         0.105788E+00, 0.794084E-01, 0.593201E-01 /)
!!$    Integer df
!!$    Integer, Save, Dimension ( n_max ) :: df_vec = (/ &
!!$         1,   2,   3, &
!!$         1,   2,   3, &
!!$         1,   2,   3, &
!!$         1,   2,   3, &
!!$         60,  80, 100, &
!!$         1,   2,   3, &
!!$         10,  10,  10, &
!!$         10,  10,  10, &
!!$         10,  10,  10 /)
!!$    Real ( kind = 8 ) lambda
!!$    Real, Save, Dimension ( n_max ) :: lambda_vec = (/ &
!!$         0.5E+00,  0.5E+00,  0.5E+00, &
!!$         1.0E+00,  1.0E+00,  1.0E+00, &
!!$         5.0E+00,  5.0E+00,  5.0E+00, &
!!$         20.0E+00, 20.0E+00, 20.0E+00, &
!!$         30.0E+00, 30.0E+00, 30.0E+00, &
!!$         5.0E+00,  5.0E+00,  5.0E+00, &
!!$         2.0E+00,  3.0E+00,  4.0E+00, &
!!$         2.0E+00,  3.0E+00,  4.0E+00, &
!!$         2.0E+00,  3.0E+00,  4.0E+00 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real, Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         3.000E+00,  3.000E+00,  3.000E+00, &
!!$         3.000E+00,  3.000E+00,  3.000E+00, &
!!$         3.000E+00,  3.000E+00,  3.000E+00, &
!!$         3.000E+00,  3.000E+00,  3.000E+00, &
!!$         60.000E+00, 60.000E+00, 60.000E+00, &
!!$         0.050E+00,  0.050E+00,  0.050E+00, &
!!$         4.000E+00,  4.000E+00,  4.000E+00, &
!!$         5.000E+00,  5.000E+00,  5.000E+00, &
!!$         6.000E+00,  6.000E+00,  6.000E+00 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       x = 0.0D+00
!!$       lambda = 0.0D+00
!!$       df = 0
!!$       cdf = 0.0D+00
!!$    Else
!!$       x = x_vec(n_data)
!!$       lambda = lambda_vec(n_data)
!!$       df = df_vec(n_data)
!!$       cdf = cdf_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine chi_noncentral_cdf_values

!!$  Subroutine chi_square_cdf_values ( n_data, a, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    The value of CHI_CDF ( DF, X ) can be evaluated in Mathematica by
!!$    !    commands like:
!!$    !
!!$    !      Needs["Statistics`ContinuousDistributions`"]
!!$    !      CDF[ChiSquareDistribution[DF], X ]
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    11 June 2004
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !    Stephen Wolfram,
!!$    !    The Mathematica Book,
!!$    !    Fourth Edition,
!!$    !    Wolfram Media / Cambridge University Press, 1999.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, integer A, real ( kind = 8 ) X, the arguments of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 21
!!$
!!$    Integer a
!!$    Integer, Save, Dimension ( n_max ) :: a_vec = (/ &
!!$         1,  2,  1,  2, &
!!$         1,  2,  3,  4, &
!!$         1,  2,  3,  4, &
!!$         5,  3,  3,  3, &
!!$         3,  3, 10, 10, &
!!$         10 /)
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.0796557D+00, 0.00498752D+00, 0.112463D+00,    0.00995017D+00, &
!!$         0.472911D+00,  0.181269D+00,   0.0597575D+00,   0.0175231D+00, &
!!$         0.682689D+00,  0.393469D+00,   0.198748D+00,    0.090204D+00, &
!!$         0.0374342D+00, 0.427593D+00,   0.608375D+00,    0.738536D+00, &
!!$         0.828203D+00,  0.88839D+00,    0.000172116D+00, 0.00365985D+00, &
!!$         0.0185759D+00 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         0.01D+00, 0.01D+00, 0.02D+00, 0.02D+00, &
!!$         0.40D+00, 0.40D+00, 0.40D+00, 0.40D+00, &
!!$         1.00D+00, 1.00D+00, 1.00D+00, 1.00D+00, &
!!$         1.00D+00, 2.00D+00, 3.00D+00, 4.00D+00, &
!!$         5.00D+00, 6.00D+00, 1.00D+00, 2.00D+00, &
!!$         3.00D+00 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       a = 0
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       a = a_vec(n_data)
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine chi_square_cdf_values

  Subroutine cumbet ( x, y, a, b, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMBET evaluates the cumulative incomplete beta distribution.
    !
    !  Discussion:
    !
    !    This routine calculates the CDF to X of the incomplete beta distribution
    !    with parameters A and B.  This is the integral from 0 to x
    !    of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
    !
    !  Reference:
    !
    !    Armido Didonato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios.
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, Number 3, September 1992, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the upper limit of integration.
    !
    !    Input, real ( kind = 8 ) Y, the value of 1-X.
    !
    !    Input, real ( kind = 8 ) A, B, the parameters of the distribution.
    !
    !    Output, real ( kind = 8 ) CUM, CCUM, the values of the cumulative
    !    density function and complementary cumulative density function.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) b
    Real ( kind = 8 ) ccum
    Real ( kind = 8 ) cum
    Integer ierr
    Real ( kind = 8 ) x
    Real ( kind = 8 ) y

    If ( x <= 0.0D+00 ) Then

       cum = 0.0
       ccum = 1.0D+00

    Else If ( y <= 0.0D+00 ) Then

       cum = 1.0D+00
       ccum = 0.0

    Else

       Call beta_inc ( a, b, x, y, cum, ccum, ierr )

    End If

    Return
  End Subroutine cumbet

!!$  Subroutine cumbin ( s, xn, pr, ompr, cum, ccum )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CUMBIN evaluates the cumulative binomial distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine returns the probability of 0 to S successes in XN binomial
!!$    !    trials, each of which has a probability of success, PR.
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.5.24.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, real ( kind = 8 ) S, the upper limit of summation.
!!$    !
!!$    !    Input, real ( kind = 8 ) XN, the number of trials.
!!$    !
!!$    !    Input, real ( kind = 8 ) PR, the probability of success in one trial.
!!$    !
!!$    !    Input, real ( kind = 8 ) OMPR, equals ( 1 - PR ).
!!$    !
!!$    !    Output, real ( kind = 8 ) CUM, the cumulative binomial distribution.
!!$    !
!!$    !    Output, real ( kind = 8 ) CCUM, the complement of the cumulative
!!$    !    binomial distribution.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) ompr
!!$    Real ( kind = 8 ) pr
!!$    Real ( kind = 8 ) s
!!$    Real ( kind = 8 ) xn
!!$
!!$    If ( s < xn ) Then
!!$
!!$       Call cumbet ( pr, ompr, s + 1.0D+00, xn - s, ccum, cum )
!!$
!!$    Else
!!$
!!$       cum = 1.0D+00
!!$       ccum = 0.0D+00
!!$
!!$    End If
!!$
!!$    Return
!!$  End Subroutine cumbin

! Subroutine cumchi ( x, df, cum, ccum )

!   !*****************************************************************************80
!   !
!   !! CUMCHI evaluates the cumulative chi-square distribution.
!   !
!   !  Parameters:
!   !
!   !    Input, real ( kind = 8 ) X, the upper limit of integration.
!   !
!   !    Input, real ( kind = 8 ) DF, the egrees of freedom of the
!   !    chi-square distribution.
!   !
!   !    Output, real ( kind = 8 ) CUM, the cumulative chi-square distribution.
!   !
!   !    Output, real ( kind = 8 ) CCUM, the complement of the cumulative
!   !    chi-square distribution.
!   !
!   Implicit None

!   Real ( kind = 8 ) a
!   Real ( kind = 8 ) ccum
!   Real ( kind = 8 ) cum
!   Real ( kind = 8 ) df
!   Real ( kind = 8 ) x
!   Real ( kind = 8 ) xx

!   a = df * 0.5D+00
!   xx = x * 0.5D+00

!   Call cumgam ( xx, a, cum, ccum )

!   Return
! End Subroutine cumchi

!!$  Subroutine cumchn ( x, df, pnonc, cum, ccum )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CUMCHN evaluates the cumulative noncentral chi-square distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine calculates the cumulative noncentral chi-square
!!$    !    distribution, i.e., the probability that a random variable
!!$    !    which follows the noncentral chi-square distribution, with
!!$    !    noncentrality parameter PNONC and continuous degrees of
!!$    !    freedom DF, is less than or equal to X.
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.4.25.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, real ( kind = 8 ) X, the upper limit of integration.
!!$    !
!!$    !    Input, real ( kind = 8 ) DF, the number of degrees of freedom.
!!$    !
!!$    !    Input, real ( kind = 8 ) PNONC, the noncentrality parameter of
!!$    !    the noncentral chi-square distribution.
!!$    !
!!$    !    Output, real ( kind = 8 ) CUM, CCUM, the CDF and complementary
!!$    !    CDF of the noncentral chi-square distribution.
!!$    !
!!$    !  Local Parameters:
!!$    !
!!$    !    Local, real ( kind = 8 ) EPS, the convergence criterion.  The sum
!!$    !    stops when a term is less than EPS*SUM.
!!$    !
!!$    !    Local, integer NTIRED, the maximum number of terms to be evaluated
!!$    !    in each sum.
!!$    !
!!$    !    Local, logical QCONV, is TRUE if convergence was achieved, that is,
!!$    !    the program did not stop on NTIRED criterion.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ) adj
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) centaj
!!$    Real ( kind = 8 ) centwt
!!$    Real ( kind = 8 ) chid2
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) df
!!$    Real ( kind = 8 ) dfd2
!!$    !Timo: Real ( kind = 8 ) dg
!!$    Real ( kind = 8 ), Parameter :: eps = 0.00001D+00
!!$    !    Real ( kind = 8 ) gamma_log
!!$    Integer i
!!$    Integer icent
!!$    Integer iterb
!!$    Integer iterf
!!$    Real ( kind = 8 ) lcntaj
!!$    Real ( kind = 8 ) lcntwt
!!$    Real ( kind = 8 ) lfact
!!$    Integer, Parameter :: ntired = 1000
!!$    Real ( kind = 8 ) pcent
!!$    Real ( kind = 8 ) pnonc
!!$    Real ( kind = 8 ) pterm
!!$    !Timo: Logical qsmall
!!$    Real ( kind = 8 ) sum1
!!$    Real ( kind = 8 ) sumadj
!!$    Real ( kind = 8 ) term
!!$    Real ( kind = 8 ) wt
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ) xnonc
!!$    !Timo: Statement function qsmall converted to an internal function
!!$    !Timo: qsmall(xx) = sum1 < 1.0D-20 .Or. xx < eps * sum1
!!$    !Timo: dg(ii) = df +  2.0D+00  * Real ( ii, kind = 8 )
!!$    !Timo: Statement function dg converted to an internal function
!!$    !Timo: dg(xi) = df +  2.0D+00  * xi
!!$
!!$    If ( x <= 0.0D+00 ) Then
!!$       cum = 0.0D+00
!!$       ccum = 1.0D+00
!!$       Return
!!$    End If
!!$    !
!!$    !  When the noncentrality parameter is (essentially) zero,
!!$    !  use cumulative chi-square distribution
!!$    !
!!$    If ( pnonc <= 1.0D-10 ) Then
!!$       Call cumchi ( x, df, cum, ccum )
!!$       Return
!!$    End If
!!$
!!$    xnonc = pnonc /  2.0D+00
!!$    !
!!$    !  The following code calculates the weight, chi-square, and
!!$    !  adjustment term for the central term in the infinite series.
!!$    !  The central term is the one in which the poisson weight is
!!$    !  greatest.  The adjustment term is the amount that must
!!$    !  be subtracted from the chi-square to move up two degrees
!!$    !  of freedom.
!!$    !
!!$    icent = Int ( xnonc )
!!$    If ( icent == 0 ) Then
!!$       icent = 1
!!$    End If
!!$
!!$    chid2 = x /  2.0D+00
!!$    !
!!$    !  Calculate central weight term.
!!$    !
!!$    lfact = gamma_log ( Real ( icent + 1, kind = 8 ) )
!!$    lcntwt = - xnonc + icent * Log ( xnonc ) - lfact
!!$    centwt = Exp ( lcntwt )
!!$    !
!!$    !  Calculate central chi-square.
!!$    !
!!$    Call cumchi ( x, dg(Real(icent,kind=8)), pcent, ccum )
!!$    !
!!$    !  Calculate central adjustment term.
!!$    !
!!$    dfd2 = dg(Real(icent,kind=8)) /  2.0D+00
!!$    lfact = gamma_log ( 1.0D+00 + dfd2 )
!!$    lcntaj = dfd2 * Log ( chid2 ) - chid2 - lfact
!!$    centaj = Exp ( lcntaj )
!!$    sum1 = centwt * pcent
!!$    !
!!$    !  Sum backwards from the central term towards zero.
!!$    !  Quit whenever either
!!$    !  (1) the zero term is reached, or
!!$    !  (2) the term gets small relative to the sum, or
!!$    !  (3) More than NTIRED terms are totaled.
!!$    !
!!$    iterb = 0
!!$    sumadj = 0.0D+00
!!$    adj = centaj
!!$    wt = centwt
!!$    i = icent
!!$    term = 0.0D+00
!!$
!!$    Do
!!$
!!$       dfd2 = dg(Real(i,kind=8)) /  2.0D+00
!!$       !
!!$       !  Adjust chi-square for two fewer degrees of freedom.
!!$       !  The adjusted value ends up in PTERM.
!!$       !
!!$       adj = adj * dfd2 / chid2
!!$       sumadj = sumadj + adj
!!$       pterm = pcent + sumadj
!!$       !
!!$       !  Adjust Poisson weight for J decreased by one.
!!$       !
!!$       wt = wt * ( i / xnonc )
!!$       term = wt * pterm
!!$       sum1 = sum1 + term
!!$       i = i - 1
!!$       iterb = iterb + 1
!!$
!!$       If ( ntired < iterb .Or. qsmall ( term ) .Or. i == 0 ) Then
!!$          Exit
!!$       End If
!!$
!!$    End Do
!!$
!!$    iterf = 0
!!$    !
!!$    !  Now sum forward from the central term towards infinity.
!!$    !  Quit when either
!!$    !    (1) the term gets small relative to the sum, or
!!$    !    (2) More than NTIRED terms are totaled.
!!$    !
!!$    sumadj = centaj
!!$    adj = centaj
!!$    wt = centwt
!!$    i = icent
!!$    !
!!$    !  Update weights for next higher J.
!!$    !
!!$    Do
!!$
!!$       wt = wt * ( xnonc / ( i + 1 ) )
!!$       !
!!$       !  Calculate PTERM and add term to sum.
!!$       !
!!$       pterm = pcent - sumadj
!!$       term = wt * pterm
!!$       sum1 = sum1 + term
!!$       !
!!$       !  Update adjustment term for DF for next iteration.
!!$       !
!!$       i = i + 1
!!$       dfd2 = dg(Real(i,kind=8)) /  2.0D+00
!!$       adj = adj * chid2 / dfd2
!!$       sumadj = sumadj + adj
!!$       iterf = iterf + 1
!!$
!!$       If ( ntired < iterf .Or. qsmall ( term ) ) Then
!!$          Exit
!!$       End If
!!$
!!$    End Do
!!$
!!$    cum = sum1
!!$    ccum = 0.5D+00 + ( 0.5D+00 - cum )
!!$
!!$    Return
!!$  Contains
!!$
!!$    Logical Function qsmall(xx)
!!$      Real ( kind = 8 ), Intent(In) :: xx
!!$      qsmall = sum1 < 1.0D-20 .Or. xx < eps * sum1
!!$      Return
!!$    End Function qsmall
!!$
!!$    Real ( kind = 8 ) Function dg(xi)
!!$      Real ( kind = 8 ), Intent(In) :: xi
!!$      dg = df +  2.0D+00  * xi
!!$      Return
!!$    End Function dg
!!$
!!$  End Subroutine cumchn

! Subroutine cumf ( f, dfn, dfd, cum, ccum )

!   !*****************************************************************************80
!   !
!   !! CUMF evaluates the cumulative F distribution.
!   !
!   !  Discussion:
!   !
!   !    This routine computes the integral from 0 to F of the F density with DFN
!   !    numerator and DFD denominator degrees of freedom.
!   !
!   !  Reference:
!   !
!   !    Milton Abramowitz, Irene Stegun,
!   !    Handbook of Mathematical Functions
!   !    1966, Formula 26.5.28.
!   !
!   !  Parameters:
!   !
!   !    Input, real ( kind = 8 ) F, the upper limit of integration.
!   !
!   !    Input, real ( kind = 8 ) DFN, DFD, the number of degrees of
!   !    freedom for the numerator and denominator.
!   !
!   !    Output, real ( kind = 8 ) CUM, CCUM, the value of the F CDF and
!   !    the complementary F CDF.
!   !
!   Implicit None

!   Real ( kind = 8 ) ccum
!   Real ( kind = 8 ) cum
!   Real ( kind = 8 ) dfd
!   Real ( kind = 8 ) dfn
!   Real ( kind = 8 ) dsum
!   Real ( kind = 8 ) f
!   Integer ierr
!   Real ( kind = 8 ) prod
!   Real ( kind = 8 ) xx
!   Real ( kind = 8 ) yy

!   If ( f <= 0.0D+00 ) Then
!      cum = 0.0D+00
!      ccum = 1.0D+00
!      Return
!   End If

!   prod = dfn * f
!   !
!   !  XX is such that the incomplete beta with parameters
!   !  DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
!   !
!   !  YY is 1 - XX
!   !
!   !  Calculate the smaller of XX and YY accurately.
!   !
!   dsum = dfd + prod
!   xx = dfd / dsum

!   If ( 0.5D+00 < xx ) Then
!      yy = prod / dsum
!      xx = 1.0D+00 - yy
!   Else
!      yy = 1.0D+00 - xx
!   End If

!   Call beta_inc ( 0.5D+00*dfd, 0.5D+00*dfn, xx, yy, ccum, cum, ierr )

!   Return
! End Subroutine cumf

!!$  Subroutine cumfnc ( f, dfn, dfd, pnonc, cum, ccum )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CUMFNC evaluates the cumulative noncentral F distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine computes the noncentral F distribution with DFN and DFD
!!$    !    degrees of freedom and noncentrality parameter PNONC.
!!$    !
!!$    !    The series is calculated backward and forward from J = LAMBDA/2
!!$    !    (this is the term with the largest Poisson weight) until
!!$    !    the convergence criterion is met.
!!$    !
!!$    !    The sum continues until a succeeding term is less than EPS
!!$    !    times the sum (or the sum is less than 1.0D-20).  EPS is
!!$    !    set to 1.0D-4 in a data statement which can be changed.
!!$    !
!!$    !
!!$    !    The original version of this routine allowed the input values
!!$    !    of DFN and DFD to be negative (nonsensical) or zero (which
!!$    !    caused numerical overflow.)  I have forced both these values
!!$    !    to be at least 1.
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    15 June 2004
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.5.16, 26.6.17, 26.6.18, 26.6.20.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, real ( kind = 8 ) F, the upper limit of integration.
!!$    !
!!$    !    Input, real ( kind = 8 ) DFN, DFD, the number of degrees of freedom
!!$    !    in the numerator and denominator.  Both DFN and DFD must be positive,
!!$    !    and normally would be integers.  This routine requires that they
!!$    !    be no less than 1.
!!$    !
!!$    !    Input, real ( kind = 8 ) PNONC, the noncentrality parameter.
!!$    !
!!$    !    Output, real ( kind = 8 ) CUM, CCUM, the noncentral F CDF and
!!$    !    complementary CDF.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ) adn
!!$    Real ( kind = 8 ) arg1
!!$    Real ( kind = 8 ) aup
!!$    Real ( kind = 8 ) b
!!$    Real ( kind = 8 ) betdn
!!$    Real ( kind = 8 ) betup
!!$    Real ( kind = 8 ) centwt
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) dfd
!!$    Real ( kind = 8 ) dfn
!!$    Real ( kind = 8 ) dnterm
!!$    Real ( kind = 8 ) dsum
!!$    Real ( kind = 8 ) dummy
!!$    Real ( kind = 8 ), Parameter :: eps = 0.0001D+00
!!$    Real ( kind = 8 ) f
!!$    !    Real ( kind = 8 ) gamma_log
!!$    Integer i
!!$    Integer icent
!!$    Integer ierr
!!$    Real ( kind = 8 ) pnonc
!!$    Real ( kind = 8 ) prod
!!$    !Timo: Logical qsmall
!!$    Real ( kind = 8 ) sum1
!!$    Real ( kind = 8 ) upterm
!!$    Real ( kind = 8 ) xmult
!!$    Real ( kind = 8 ) xnonc
!!$    Real ( kind = 8 ) xx
!!$    Real ( kind = 8 ) yy
!!$
!!$    !Timo: Statement function qsmall converted to an internal function
!!$    !Timo: qsmall(x) = sum1 < 1.0D-20 .Or. x < eps * sum1
!!$
!!$    If ( f <= 0.0D+00 ) Then
!!$       cum = 0.0D+00
!!$       ccum = 1.0D+00
!!$       Return
!!$    End If
!!$
!!$    If ( dfn < 1.0D+00 ) Then
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CUMFNC - Fatal error!'
!!$       Write ( *, '(a)' ) '  DFN < 1.'
!!$       Stop
!!$    End If
!!$
!!$    If ( dfd < 1.0D+00 ) Then
!!$       Write ( *, '(a)' ) ' '
!!$       Write ( *, '(a)' ) 'CUMFNC - Fatal error!'
!!$       Write ( *, '(a)' ) '  DFD < 1.'
!!$       Stop
!!$    End If
!!$    !
!!$    !  Handle case in which the noncentrality parameter is essentially zero.
!!$    !
!!$    If ( pnonc < 1.0D-10 ) Then
!!$       Call cumf ( f, dfn, dfd, cum, ccum )
!!$       Return
!!$    End If
!!$
!!$    xnonc = pnonc /  2.0D+00
!!$    !
!!$    !  Calculate the central term of the Poisson weighting factor.
!!$    !
!!$    icent = Int ( xnonc )
!!$
!!$    If ( icent == 0 ) Then
!!$       icent = 1
!!$    End If
!!$    !
!!$    !  Compute central weight term.
!!$    !
!!$    centwt = Exp ( -xnonc + icent * Log ( xnonc ) &
!!$         - gamma_log ( Real ( icent + 1, kind = 8  ) ) )
!!$    !
!!$    !  Compute central incomplete beta term.
!!$    !  Ensure that minimum of arg to beta and 1 - arg is computed accurately.
!!$    !
!!$    prod = dfn * f
!!$    dsum = dfd + prod
!!$    yy = dfd / dsum
!!$
!!$    If ( 0.5D+00 < yy ) Then
!!$       xx = prod / dsum
!!$       yy = 1.0D+00 - xx
!!$    Else
!!$       xx = 1.0D+00 - yy
!!$    End If
!!$
!!$    arg1 = 0.5D+00 * dfn + Real ( icent, kind = 8 )
!!$    Call beta_inc ( arg1, 0.5D+00*dfd, xx, yy, betdn, &
!!$         dummy, ierr )
!!$
!!$    adn = dfn / 2.0D+00 + Real ( icent, kind = 8 )
!!$    aup = adn
!!$    b = dfd / 2.0D+00
!!$    betup = betdn
!!$    sum1 = centwt * betdn
!!$    !
!!$    !  Now sum terms backward from ICENT until convergence or all done.
!!$    !
!!$    xmult = centwt
!!$    i = icent
!!$    dnterm = Exp ( gamma_log ( adn + b ) &
!!$         - gamma_log ( adn + 1.0D+00 ) &
!!$         - gamma_log ( b ) + adn * Log ( xx ) + b * Log ( yy ) )
!!$
!!$    Do
!!$
!!$       If ( qsmall ( xmult * betdn ) .Or. i <= 0 ) Then
!!$          Exit
!!$       End If
!!$
!!$       xmult = xmult * ( Real ( i, kind = 8 ) / xnonc )
!!$       i = i - 1
!!$       adn = adn - 1.0D+00
!!$       dnterm = ( adn + 1.0D+00 ) / ( ( adn + b ) * xx ) * dnterm
!!$       betdn = betdn + dnterm
!!$       sum1 = sum1 + xmult * betdn
!!$
!!$    End Do
!!$
!!$    i = icent + 1
!!$    !
!!$    !  Now sum forward until convergence.
!!$    !
!!$    xmult = centwt
!!$
!!$    If ( ( aup - 1.0D+00 + b ) == 0 ) Then
!!$
!!$       upterm = Exp ( - gamma_log ( aup ) - gamma_log ( b ) &
!!$            + ( aup - 1.0D+00 ) * Log ( xx ) + b * Log ( yy ) )
!!$
!!$    Else
!!$
!!$       upterm = Exp ( gamma_log ( aup - 1.0D+00 + b ) - gamma_log ( aup ) &
!!$            - gamma_log ( b ) + ( aup - 1.0D+00 ) * Log ( xx ) + b * Log ( yy ) )
!!$
!!$    End If
!!$
!!$    Do
!!$
!!$       xmult = xmult * ( xnonc / i )
!!$       i = i + 1
!!$       aup = aup + 1.0D+00
!!$       upterm = ( aup + b -  2.0D+00  ) * xx / ( aup - 1.0D+00 ) * upterm
!!$       betup = betup - upterm
!!$       sum1 = sum1 + xmult * betup
!!$
!!$       If ( qsmall ( xmult * betup ) ) Then
!!$          Exit
!!$       End If
!!$
!!$    End Do
!!$
!!$    cum = sum1
!!$    ccum = 0.5D+00 + ( 0.5D+00 - cum )
!!$
!!$    Return
!!$  Contains
!!$
!!$    Logical Function qsmall(x)
!!$      Real ( kind = 8 ), Intent(In) :: x
!!$      qsmall = sum1 < 1.0D-20 .Or. x < eps * sum1
!!$      Return
!!$    End Function qsmall
!!$
!!$  End Subroutine cumfnc

  Subroutine cumgam ( x, a, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMGAM evaluates the cumulative incomplete gamma distribution.
    !
    !  Discussion:
    !
    !    This routine computes the cumulative distribution function of the
    !    incomplete gamma distribution, i.e., the integral from 0 to X of
    !
    !      (1/GAM(A))*EXP(-T)*T**(A-1) DT
    !
    !    where GAM(A) is the complete gamma function of A, i.e.,
    !
    !      GAM(A) = integral from 0 to infinity of EXP(-T)*T**(A-1) DT
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the upper limit of integration.
    !
    !    Input, real ( kind = 8 ) A, the shape parameter of the incomplete
    !    Gamma distribution.
    !
    !    Output, real ( kind = 8 ) CUM, CCUM, the incomplete Gamma CDF and
    !    complementary CDF.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) ccum
    Real ( kind = 8 ) cum
    Real ( kind = 8 ) x

    If ( x <= 0.0D+00 ) Then

       cum = 0.0D+00
       ccum = 1.0D+00

    Else

       Call gamma_inc ( a, x, cum, ccum, 0 )

    End If

    Return
  End Subroutine cumgam

!!$  Subroutine cumnbn ( f, s, pr, ompr, cum, ccum )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CUMNBN evaluates the cumulative negative binomial distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine returns the probability that there will be F or
!!$    !    fewer failures before there are S successes, with each binomial
!!$    !    trial having a probability of success PR.
!!$    !
!!$    !    Prob(# failures = F | S successes, PR)  =
!!$    !                        ( S + F - 1 )
!!$    !                        (            ) * PR^S * (1-PR)^F
!!$    !                        (      F     )
!!$    !
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions
!!$    !    1966, Formula 26.5.26.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, real ( kind = 8 ) F, the number of failures.
!!$    !
!!$    !    Input, real ( kind = 8 ) S, the number of successes.
!!$    !
!!$    !    Input, real ( kind = 8 ) PR, OMPR, the probability of success on
!!$    !    each binomial trial, and the value of (1-PR).
!!$    !
!!$    !    Output, real ( kind = 8 ) CUM, CCUM, the negative binomial CDF,
!!$    !    and the complementary CDF.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) f
!!$    Real ( kind = 8 ) ompr
!!$    Real ( kind = 8 ) pr
!!$    Real ( kind = 8 ) s
!!$
!!$    Call cumbet ( pr, ompr, s, f+1.D+00, cum, ccum )
!!$
!!$    Return
!!$  End Subroutine cumnbn

  Subroutine cumnor ( arg, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMNOR computes the cumulative normal distribution.
    !
    !  Discussion:
    !
    !    This function evaluates the normal distribution function:
    !
    !                              / x
    !                     1       |       -t*t/2
    !          P(x) = ----------- |      e       dt
    !                 sqrt(2 pi)  |
    !                             /-oo
    !
    !    This transportable program uses rational functions that
    !    theoretically approximate the normal distribution function to
    !    at least 18 significant decimal digits.  The accuracy achieved
    !    depends on the arithmetic system, the compiler, the intrinsic
    !    functions, and proper selection of the machine dependent
    !    constants.
    !
    !  Author:
    !
    !    William Cody
    !    Mathematics and Computer Science Division
    !    Argonne National Laboratory
    !    Argonne, IL 60439
    !
    !  Reference:
    !
    !    William Cody,
    !    Rational Chebyshev approximations for the error function,
    !    Mathematics of Computation,
    !    1969, pages 631-637.
    !
    !    William Cody,
    !    Algorithm 715:
    !    SPECFUN - A Portable FORTRAN Package of Special Function Routines
    !    and Test Drivers,
    !    ACM Transactions on Mathematical Software,
    !    Volume 19, Number 1, 1993, pages 22-32.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) ARG, the upper limit of integration.
    !
    !    Output, real ( kind = 8 ) CUM, CCUM, the Normal density CDF and
    !    complementary CDF.
    !
    !  Local Parameters:
    !
    !    Local, real ( kind = 8 ) EPS, the argument below which anorm(x)
    !    may be represented by 0.5 and above which  x*x  will not underflow.
    !    A conservative value is the largest machine number X
    !    such that   1.0D+00 + X = 1.0D+00   to machine precision.
    !
    Implicit None

    Real ( kind = 8 ), Parameter, Dimension ( 5 ) :: a = (/ &
         2.2352520354606839287D+00, &
         1.6102823106855587881D+02, &
         1.0676894854603709582D+03, &
         1.8154981253343561249D+04, &
         6.5682337918207449113D-02 /)
    Real ( kind = 8 ) arg
    Real ( kind = 8 ), Parameter, Dimension ( 4 ) :: b = (/ &
         4.7202581904688241870D+01, &
         9.7609855173777669322D+02, &
         1.0260932208618978205D+04, &
         4.5507789335026729956D+04 /)
    Real ( kind = 8 ), Parameter, Dimension ( 9 ) :: c = (/ &
         3.9894151208813466764D-01, &
         8.8831497943883759412D+00, &
         9.3506656132177855979D+01, &
         5.9727027639480026226D+02, &
         2.4945375852903726711D+03, &
         6.8481904505362823326D+03, &
         1.1602651437647350124D+04, &
         9.8427148383839780218D+03, &
         1.0765576773720192317D-08 /)
    Real ( kind = 8 ) ccum
    Real ( kind = 8 ) cum
    Real ( kind = 8 ), Parameter, Dimension ( 8 ) :: d = (/ &
         2.2266688044328115691D+01, &
         2.3538790178262499861D+02, &
         1.5193775994075548050D+03, &
         6.4855582982667607550D+03, &
         1.8615571640885098091D+04, &
         3.4900952721145977266D+04, &
         3.8912003286093271411D+04, &
         1.9685429676859990727D+04 /)
    Real ( kind = 8 ) del
    Real ( kind = 8 ) eps
    Integer i
    Real ( kind = 8 ), Parameter, Dimension ( 6 ) :: p = (/ &
         2.1589853405795699D-01, &
         1.274011611602473639D-01, &
         2.2235277870649807D-02, &
         1.421619193227893466D-03, &
         2.9112874951168792D-05, &
         2.307344176494017303D-02 /)
    Real ( kind = 8 ), Parameter, Dimension ( 5 ) :: q = (/ &
         1.28426009614491121D+00, &
         4.68238212480865118D-01, &
         6.59881378689285515D-02, &
         3.78239633202758244D-03, &
         7.29751555083966205D-05 /)
    Real ( kind = 8 ), Parameter :: root32 = 5.656854248D+00
    Real ( kind = 8 ), Parameter :: sixten = 16.0D+00
    Real ( kind = 8 ) temp
    Real ( kind = 8 ), Parameter :: sqrpi = 3.9894228040143267794D-01
    Real ( kind = 8 ), Parameter :: thrsh = 0.66291D+00
    Real ( kind = 8 ) x
    Real ( kind = 8 ) xden
    Real ( kind = 8 ) xnum
    Real ( kind = 8 ) y
    Real ( kind = 8 ) xsq
    !
    !  Machine dependent constants
    !
    eps = Epsilon ( 1.0D+00 ) * 0.5D+00

    x = arg
    y = Abs ( x )

    If ( y <= thrsh ) Then
       !
       !  Evaluate  anorm  for  |X| <= 0.66291
       !
       If ( eps < y ) Then
          xsq = x * x
       Else
          xsq = 0.0D+00
       End If

       xnum = a(5) * xsq
       xden = xsq
       Do i = 1, 3
          xnum = ( xnum + a(i) ) * xsq
          xden = ( xden + b(i) ) * xsq
       End Do
       cum = x * ( xnum + a(4) ) / ( xden + b(4) )
       temp = cum
       cum = 0.5D+00 + temp
       ccum = 0.5D+00 - temp
       !
       !  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
       !
    Else If ( y <= root32 ) Then

       xnum = c(9) * y
       xden = y
       Do i = 1, 7
          xnum = ( xnum + c(i) ) * y
          xden = ( xden + d(i) ) * y
       End Do
       cum = ( xnum + c(8) ) / ( xden + d(8) )
       xsq = Aint ( y * sixten ) / sixten
       del = ( y - xsq ) * ( y + xsq )
       cum = Exp ( - xsq * xsq * 0.5D+00 ) * Exp ( -del * 0.5D+00 ) * cum
       ccum = 1.0D+00 - cum

       If ( 0.0D+00 < x ) Then
          Call r8_swap ( cum, ccum )
       End If
       !
       !  Evaluate ANORM for sqrt(32) < |X|.
       !
    Else

       cum = 0.0D+00
       xsq = 1.0D+00 / ( x * x )
       xnum = p(6) * xsq
       xden = xsq
       Do i = 1, 4
          xnum = ( xnum + p(i) ) * xsq
          xden = ( xden + q(i) ) * xsq
       End Do

       cum = xsq * ( xnum + p(5) ) / ( xden + q(5) )
       cum = ( sqrpi - cum ) / y
       xsq = Aint ( x * sixten ) / sixten
       del = ( x - xsq ) * ( x + xsq )
       cum = Exp ( - xsq * xsq * 0.5D+00 ) &
            * Exp ( - del * 0.5D+00 ) * cum
       ccum = 1.0D+00 - cum

       If ( 0.0D+00 < x ) Then
          Call r8_swap ( cum, ccum )
       End If

    End If

    If ( cum < Tiny ( cum ) ) Then
       cum = 0.0D+00
    End If

    If ( ccum < Tiny ( ccum ) ) Then
       ccum = 0.0D+00
    End If

    Return
  End Subroutine cumnor

!!$  Subroutine cumpoi ( s, xlam, cum, ccum )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CUMPOI evaluates the cumulative Poisson distribution.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    This routine returns the probability of S or fewer events in a Poisson
!!$    !    distribution with mean XLAM.
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    Formula 26.4.21.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, real ( kind = 8 ) S, the upper limit of cumulation of the
!!$    !    Poisson density function.
!!$    !
!!$    !    Input, real ( kind = 8 ) XLAM, the mean of the Poisson distribution.
!!$    !
!!$    !    Output, real ( kind = 8 ) CUM, CCUM, the Poisson density CDF and
!!$    !    complementary CDF.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) chi
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) df
!!$    Real ( kind = 8 ) s
!!$    Real ( kind = 8 ) xlam
!!$
!!$    df =  2.0D+00  * ( s + 1.0D+00 )
!!$    chi =  2.0D+00  * xlam
!!$
!!$    Call cumchi ( chi, df, ccum, cum )
!!$
!!$    Return
!!$  End Subroutine cumpoi

!!$  Subroutine cumt ( t, df, cum, ccum )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! CUMT evaluates the cumulative T distribution.
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    Formula 26.5.27.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input, real ( kind = 8 ) T, the upper limit of integration.
!!$    !
!!$    !    Input, real ( kind = 8 ) DF, the number of degrees of freedom of
!!$    !    the T distribution.
!!$    !
!!$    !    Output, real ( kind = 8 ) CUM, CCUM, the T distribution CDF and
!!$    !    complementary CDF.
!!$    !
!!$    Implicit None
!!$
!!$    Real ( kind = 8 ) a
!!$    Real ( kind = 8 ) ccum
!!$    Real ( kind = 8 ) cum
!!$    Real ( kind = 8 ) df
!!$    Real ( kind = 8 ) oma
!!$    Real ( kind = 8 ) t
!!$    Real ( kind = 8 ) xx
!!$    Real ( kind = 8 ) yy
!!$
!!$    xx = df / ( df + t**2 )
!!$    yy = t**2 / ( df + t**2 )
!!$
!!$    Call cumbet ( xx, yy, 0.5D+00*df, 0.5D+00, a, oma )
!!$
!!$    If ( t <= 0.0D+00 ) Then
!!$       cum = 0.5D+00 * a
!!$       ccum = oma + cum
!!$    Else
!!$       ccum = 0.5D+00 * a
!!$       cum = oma + ccum
!!$    End If
!!$
!!$    Return
!!$  End Subroutine cumt

  Function dbetrm ( a, b )

    !*****************************************************************************80
    !
    !! DBETRM computes the Sterling remainder for the complete beta function.
    !
    !  Discussion:
    !
    !    Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
    !    where Lgamma is the log of the (complete) gamma function
    !
    !    Let ZZ be approximation obtained if each log gamma is approximated
    !    by Sterling's formula, i.e.,
    !
    !      Sterling(Z) = log ( sqrt ( 2 * PI ) ) + ( Z - 0.5 ) * log ( Z ) - Z
    !
    !    The Sterling remainder is Log(Beta(A,B)) - ZZ.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, the parameters of the Beta function.
    !
    !    Output, real ( kind = 8 ) DBETRM, the Sterling remainder.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) b
    Real ( kind = 8 ) dbetrm
    !    Real ( kind = 8 ) dstrem
    !
    !  Try to sum from smallest to largest
    !
    dbetrm = -dstrem ( a + b )
    dbetrm = dbetrm + dstrem ( Max ( a, b ) )
    dbetrm = dbetrm + dstrem ( Min ( a, b ) )

    Return
  End Function dbetrm

  Function dexpm1 ( x )

    !*****************************************************************************80
    !
    !! DEXPM1 evaluates the function EXP(X) - 1.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the value at which exp(X)-1 is desired.
    !
    !    Output, real ( kind = 8 ) DEXPM1, the value of exp(X)-1.
    !
    Implicit None

    Real ( kind = 8 ) bot
    Real ( kind = 8 ) dexpm1
    Real ( kind = 8 ), Parameter :: p1 =  0.914041914819518D-09
    Real ( kind = 8 ), Parameter :: p2 =  0.238082361044469D-01
    Real ( kind = 8 ), Parameter :: q1 = -0.499999999085958D+00
    Real ( kind = 8 ), Parameter :: q2 =  0.107141568980644D+00
    Real ( kind = 8 ), Parameter :: q3 = -0.119041179760821D-01
    Real ( kind = 8 ), Parameter :: q4 =  0.595130811860248D-03
    Real ( kind = 8 ) top
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x

    If ( Abs ( x ) <= 0.15D+00 ) Then

       top = ( p2 * x + p1 ) * x + 1.0D+00
       bot = ((( q4 * x + q3 ) * x + q2 ) * x + q1 ) * x + 1.0D+00
       dexpm1 = x * ( top / bot )

    Else

       w = Exp ( x )

       If ( x <= 0.0D+00 ) Then
          dexpm1 = ( w - 0.5D+00 ) - 0.5D+00
       Else
          dexpm1 = w * ( 0.5D+00 &
               + ( 0.5D+00 - 1.0D+00 / w ))
       End If

    End If

    Return
  End Function dexpm1

  Function dinvnr ( p, q )

    !*****************************************************************************80
    !
    !! DINVNR computes the inverse of the normal distribution.
    !
    !  Discussion:
    !
    !    This routine returns X such that
    !
    !      CUMNOR(X) = P,
    !
    !    that is, so that
    !
    !      P = integral ( -Infinity <= T <= X ) exp(-U*U/2)/sqrt(2*PI) dU
    !
    !    The rational function on page 95 of Kennedy and Gentle is used as a
    !    starting value for the Newton method of finding roots.
    !
    !  Reference:
    !
    !    William Kennedy, James Gentle,
    !    Statistical Computing,
    !    Marcel Dekker, NY, 1980,
    !    QA276.4 K46
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) P, Q, the probability, and the complementary
    !    probability.
    !
    !    Output, real ( kind = 8 ) DINVNR, the argument X for which the
    !    Normal CDF has the value P.
    !
    Implicit None

    Real ( kind = 8 ) ccum
    Real ( kind = 8 ) cum
    Real ( kind = 8 ) dinvnr
    Real ( kind = 8 ) dx
    Real ( kind = 8 ), Parameter :: eps = 1.0D-13
    Integer i
    Integer, Parameter :: maxit = 100
    Real ( kind = 8 ) p
    Real ( kind = 8 ) pp
    Real ( kind = 8 ) q
    Real ( kind = 8 ), Parameter :: r2pi = 0.3989422804014326D+00
    Real ( kind = 8 ) strtx
    !    Real ( kind = 8 ) stvaln
    Real ( kind = 8 ) xcur

    pp = Min ( p, q )
    strtx = stvaln ( pp )
    xcur = strtx
    !
    !  Newton iterations.
    !
    Do i = 1, maxit

       Call cumnor ( xcur, cum, ccum )
       dx = ( cum - pp ) / ( r2pi * Exp ( -0.5D+00 * xcur * xcur ) )
       xcur = xcur - dx

       If ( Abs ( dx / xcur ) < eps ) Then
          If ( p <= q ) Then
             dinvnr = xcur
          Else
             dinvnr = -xcur
          End If
          Return
       End If

    End Do

    If ( p <= q ) Then
       dinvnr = strtx
    Else
       dinvnr = -strtx
    End If

    Return
  End Function dinvnr

  ! Subroutine dinvr ( status, x, fx, qleft, qhi )
  Subroutine dinvr ( status, x, fx, qleft, qhi, zsmall, zbig, zabsst, zrelst, zstpmu, zabsto, zrelto, is_entry )

    !*****************************************************************************80
    !
    !! DINVR bounds the zero of the function and invokes DZROR.
    !
    !  Discussion:
    !
    !    This routine seeks to find bounds on a root of the function and
    !    invokes ZROR to perform the zero finding.  STINVR must have been
    !    called before this routine in order to set its parameters.
    !
    !  Reference:
    !
    !    JCP Bus, TJ Dekker,
    !    Two Efficient Algorithms with Guaranteed Convergence for
    !    Finding a Zero of a Function,
    !    ACM Transactions on Mathematical Software,
    !    Volume 1, Number 4, pages 330-345, 1975.
    !
    !  Parameters:
    !
    !    Input/output, integer STATUS.  At the beginning of a zero finding
    !    problem, STATUS should be set to 0 and INVR invoked.  The value
    !    of parameters other than X will be ignored on this call.
    !    If INVR needs the function to be evaluated, it will set STATUS to 1
    !    and return.  The value of the function should be set in FX and INVR
    !    again called without changing any of its other parameters.
    !    If INVR finishes without error, it returns with STATUS 0, and X an
    !    approximate root of F(X).
    !    If INVR cannot bound the function, it returns a negative STATUS and
    !    sets QLEFT and QHI.
    !
    !    Output, real ( kind = 8 ) X, the value at which F(X) is to be evaluated.
    !
    !    Input, real ( kind = 8 ) FX, the value of F(X) calculated by the user
    !    on the previous call, when INVR returned with STATUS = 1.
    !
    !    Output, logical QLEFT, is defined only if QMFINV returns FALSE.  In that
    !    case, QLEFT is TRUE if the stepping search terminated unsucessfully
    !    at SMALL, and FALSE if the search terminated unsucessfully at BIG.
    !
    !    Output, logical QHI, is defined only if QMFINV returns FALSE.  In that
    !    case, it is TRUE if Y < F(X) at the termination of the search and FALSE
    !    if F(X) < Y.
    !
    Implicit None

    ! Timo Korhonen 2014: Remove the Entry dstzr
    LOGICAL, INTENT(IN) :: is_entry

    Real ( kind = 8 ), Save :: absstp
    Real ( kind = 8 ), Save :: abstol
    Real ( kind = 8 ), Save :: big
    Real ( kind = 8 ), Save :: fbig
    Real ( kind = 8 ), Save :: fsmall
    Real ( kind = 8 ) fx
    Integer, Save :: i99999
    Logical, Save :: qbdd
    Logical, Save :: qcond
    Logical, Save :: qdum1
    Logical, Save :: qdum2
    Logical qhi
    Logical, Save :: qincr
    Logical qleft
    Logical, Save :: qlim
    Logical, Save :: qup
    Real ( kind = 8 ), Save :: relstp
    Real ( kind = 8 ), Save :: reltol
    Real ( kind = 8 ), Save :: small
    Integer status
    Real ( kind = 8 ), Save :: step
    Real ( kind = 8 ), Save :: stpmul
    Real ( kind = 8 ) x
    Real ( kind = 8 ), Save :: xhi
    Real ( kind = 8 ), Save :: xlb
    Real ( kind = 8 ), Save :: xlo
    Real ( kind = 8 ), Save :: xsave
    Real ( kind = 8 ), Save :: xub
    Real ( kind = 8 ), Save :: yy
    Real ( kind = 8 ) zabsst
    Real ( kind = 8 ) zabsto
    Real ( kind = 8 ) zbig
    Real ( kind = 8 ) zrelst
    Real ( kind = 8 ) zrelto
    Real ( kind = 8 ) zsmall
    Real ( kind = 8 ) zstpmu
!Timo:    Save

    ! Timo Korhonen 2014: Remove the Entry dstzr
    If ( is_entry ) Then
       Go to 280 ! Entry dstinv
    End If

    If ( 0 < status ) Then
       ! go to i99999
       Select Case(i99999)
          Case(10)
             Go To 10
          Case(20)
             Go To 20
          Case(90)
             Go To 90
          Case(130)
             Go To 130
          Case(200)
             Go To 200
          Case(270)
             Go To 270
       End Select
    End If

    qcond = .Not. ( small <= x .And. x <= big )

    If ( .Not. ( small <= x .And. x <= big ) ) Then
       Write ( *, '(a)' ) ' '
       Write ( *, '(a)' ) 'DINVR - Fatal error!'
       Write ( *, '(a)' ) '  The values SMALL, X, BIG are not monotone.'
       Stop
    End If

    xsave = x
    !
    !  See that SMALL and BIG bound the zero and set QINCR.
    !
    x = small
    !
    !  GET-function-VALUE
    !
    ! Assign 10 to i99999
    i99999 = 10
    status = 1
    Return

10  Continue

    fsmall = fx
    x = big
    !
    !  GET-function-VALUE
    !
    ! Assign 20 to i99999
    i99999 = 20
    status = 1
    Return

20  Continue

    fbig = fx

    qincr = ( fsmall < fbig )

    If ( fsmall <= fbig ) Then

       If ( 0.0D+00 < fsmall ) Then
          status = -1
          qleft = .True.
          qhi = .True.
          Return
       End If

       If ( fbig < 0.0D+00 ) Then
          status = -1
          qleft = .False.
          qhi = .False.
          Return
       End If

    Else If ( fbig < fsmall ) Then

       If ( fsmall < 0.0D+00 ) Then
          status = -1
          qleft = .True.
          qhi = .False.
          Return
       End If

       If ( 0.0D+00 < fbig ) Then
          status = -1
          qleft = .False.
          qhi = .True.
          Return
       End If

    End If

    x = xsave
    step = Max ( absstp, relstp * Abs ( x ) )
    !
    !  YY = F(X) - Y
    !  GET-function-VALUE
    !
    ! Assign 90 to i99999
    i99999 = 90
    status = 1
    Return

90  Continue

    yy = fx

    If ( yy == 0.0D+00 ) Then
       status = 0
       Return
    End If

! 100 Continue

    qup = ( qincr .And. ( yy < 0.0D+00 ) ) .Or. &
         ( .Not. qincr .And. ( 0.0D+00 < yy ) )
    !
    !  Handle case in which we must step higher.
    !
    If (.Not. qup ) Then
       go to 170
    End If

    xlb = xsave
    xub = Min ( xlb + step, big )
    go to 120

110 Continue

    If ( qcond ) Then
       go to 150
    End If
    !
    !  YY = F(XUB) - Y
    !
120 Continue

    x = xub
    !
    !  GET-function-VALUE
    !
    ! Assign 130 to i99999
    i99999 = 130
    status = 1
    Return

130 Continue

    yy = fx
    qbdd = ( qincr .And. ( 0.0D+00 <= yy ) ) .Or. &
         ( .Not. qincr .And. ( yy <= 0.0D+00 ) )
    qlim = ( big <= xub )
    qcond = qbdd .Or. qlim

    If ( .Not. qcond ) Then
       step = stpmul * step
       xlb = xub
       xub = Min ( xlb + step, big )
    End If

    go to 110

150 Continue

    If ( qlim .And. .Not. qbdd ) Then
       status = -1
       qleft = .False.
       qhi = .Not. qincr
       x = big
       Return
    End If

! 160 Continue

    go to 240
    !
    !  Handle the case in which we must step lower.
    !
170 Continue

    xub = xsave
    xlb = Max ( xub - step, small )
    go to 190

180 Continue

    If ( qcond ) Then
       go to 220
    End If
    !
    !  YY = F(XLB) - Y
    !
190 Continue

    x = xlb
    !
    !  GET-function-VALUE
    !
    ! Assign 200 to i99999
    i99999 = 200
    status = 1
    Return

200 Continue

    yy = fx
    qbdd = ( qincr .And. ( yy <= 0.0D+00 ) ) .Or. &
         ( .Not. qincr .And. ( 0.0D+00 <= yy ) )
    qlim = xlb <= small
    qcond = qbdd .Or. qlim

    If ( .Not. qcond ) Then
       step = stpmul * step
       xub = xlb
       xlb = Max ( xub - step, small )
    End If

    go to 180

220 Continue

    If ( qlim .And. ( .Not. qbdd ) ) Then
       status = -1
       qleft = .True.
       qhi = qincr
       x = small
       Return
    End If

! 230 Continue
240 Continue

    ! Call dstzr ( xlb, xub, abstol, reltol )
    Call dzror ( status, abstol, reltol, xlb, xub, qleft, qhi, .TRUE. )
    !
    !  If we reach here, XLB and XUB bound the zero of F.
    !
    status = 0
    go to 260

250 Continue

    If ( status /= 1 ) Then
       x = xlo
       status = 0
       Return
    End If

260 Continue

    Call dzror ( status, x, fx, xlo, xhi, qdum1, qdum2, .FALSE. )

    If ( status /= 1 ) Then
       go to 250
    End If
    !
    !  GET-function-VALUE
    !
    ! Assign 270 to i99999
    i99999 = 270
    status = 1
    Return

270 Continue
    go to 250

    ! Entry dstinv ( zsmall, zbig, zabsst, zrelst, zstpmu, zabsto, zrelto )
280 Continue
    !*****************************************************************************80
    !
    !! DSTINV SeT INverse finder - Reverse Communication
    !
    !  Discussion:
    !
    !    This routine is given a monotone function F, and a value Y,
    !    and seeks an argument value X such that F(X) = Y.
    !
    !    This routine uses reverse communication -- see invr.
    !    This routine sets quantities needed by INVR.
    !
    !    F must be a monotone function, the results of QMFINV are
    !    otherwise undefined.  QINCR must be TRUE if F is nondecreasing
    !    and FALSE if F is nonincreasing.
    !
    !    QMFINV will return TRUE if and only if F(SMALL) and
    !    F(BIG) bracket Y, i. e.,
    !      QINCR is TRUE and F(SMALL) <= Y <= F(BIG) or
    !      QINCR is FALSE and F(BIG) <= Y <= F(SMALL)
    !
    !    if QMFINV returns TRUE, then the X returned satisfies
    !    the following condition.  let
    !      TOL(X) = MAX ( ABSTOL, RELTOL * ABS ( X ) )
    !    then if QINCR is TRUE,
    !      F(X-TOL(X)) <= Y <= F(X+TOL(X))
    !    and if QINCR is FALSE
    !      F(X-TOL(X)) >= Y >= F(X+TOL(X))
    !
    !    Compares F(X) with Y for the input value of X then uses QINCR
    !    to determine whether to step left or right to bound the
    !    desired x.  the initial step size is
    !
    !      max ( ABSSTP, RELSTP * ABS ( S ) )
    !
    !    for the input value of X.
    !
    !    Iteratively steps right or left until it bounds X.
    !    At each step which doesn't bound X, the step size is doubled.
    !    The routine is careful never to step beyond SMALL or BIG.  If
    !    it hasn't bounded X at SMALL or BIG, QMFINV returns FALSE
    !    after setting QLEFT and QHI.
    !
    !    If X is successfully bounded then Algorithm R of the paper
    !    Bus and Dekker is employed to find the zero of the function F(X)-Y.
    !    This is routine QRZERO.
    !
    !  Reference:
    !
    !    JCP Bus, TJ Dekker,
    !    Two Efficient Algorithms with Guaranteed Convergence for
    !    Finding a Zero of a Function,
    !    ACM Transactions on Mathematical Software,
    !    Volume 1, Number 4, pages 330-345, 1975.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) ZSMALL, ZBIG, the left and right endpoints
    !    of the interval to be searched for a solution.
    !
    !    Input, real ( kind = 8 ) ZABSST, ZRELSTP, the initial step size in
    !    the search is max ( ZABSST, ZRELST * abs ( X ) ).
    !
    !    Input, real ( kind = 8 ) STPMUL.  When a step doesn't bound the zero,
    !    the stepsize is multiplied by STPMUL and another step taken.  A
    !    popular value is 2.0.
    !
    !    Input, real ( kind = 8 ) ABSTOL, RELTOL, two numbers that determine
    !    the accuracy of the solution
    !
    small = zsmall
    big = zbig
    absstp = zabsst
    relstp = zrelst
    stpmul = zstpmu
    abstol = zabsto
    reltol = zrelto

    Return
  End Subroutine dinvr

  Function dlanor ( x )

    !*****************************************************************************80
    !
    !! DLANOR evaluates the logarithm of the asymptotic Normal CDF.
    !
    !  Discussion:
    !
    !    This routine computes the logarithm of the cumulative normal distribution
    !    from abs ( x ) to infinity for  5 <= abs ( X ).
    !
    !    The relative error at X = 5 is about 0.5D-5.
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions
    !    1966, Formula 26.2.12.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the value at which the Normal CDF is to be
    !    evaluated.  It is assumed that 5 <= abs ( X ).
    !
    !    Output, real ( kind = 8 ) DLANOR, the logarithm of the asymptotic
    !    Normal CDF.
    !
    Implicit None

    !    Real ( kind = 8 ) alnrel
    Real ( kind = 8 ) approx
    Real ( kind = 8 ), Save, Dimension ( 0:11 ) :: coef = (/ &
         -1.0D+00,  3.0D+00,  -15.0D+00,  105.0D+00,  -945.0D+00,  &
         10395.0D+00, -135135.0D+00,  2027025.0D+00,  -34459425.0D+00, &
         654729075.0D+00, -13749310575D+00,  316234143225.0D+00 /)
    Real ( kind = 8 ) correc
    Real ( kind = 8 ), Parameter :: dlsqpi = 0.91893853320467274177D+00
    !    Real ( kind = 8 ) eval_pol
    Real ( kind = 8 ) dlanor
    Real ( kind = 8 ) x
    Real ( kind = 8 ) xx
    Real ( kind = 8 ) xx2

    xx = Abs ( x )

    If ( Abs ( x ) < 5.0D+00 ) Then
       Write ( *, '(a)' ) ' '
       Write ( *, '(a)' ) 'DLANOR - Fatal error!'
       Write ( *, '(a)' ) '  The argument X is too small.'
    End If

    approx = - dlsqpi - 0.5D+00 * x**2 - Log ( Abs ( x ) )

    xx2 = xx * xx
    correc = eval_pol ( coef, 11, 1.0D+00 / xx2 ) / xx2
    correc = alnrel ( correc )

    dlanor = approx + correc

    Return
  End Function dlanor

  Function dstrem ( z )

    !*****************************************************************************80
    !
    !! DSTREM computes the Sterling remainder ln ( Gamma ( Z ) ) - Sterling ( Z ).
    !
    !  Discussion:
    !
    !    This routine returns
    !
    !      ln ( Gamma ( Z ) ) - Sterling ( Z )
    !
    !    where Sterling(Z) is Sterling's approximation to ln ( Gamma ( Z ) ).
    !
    !    Sterling(Z) = ln ( sqrt ( 2 * PI ) ) + ( Z - 0.5 ) * ln ( Z ) - Z
    !
    !    If 6 <= Z, the routine uses 9 terms of a series in Bernoulli numbers,
    !    with values calculated using Maple.
    !
    !    Otherwise, the difference is computed explicitly.
    !
    !  Modified:
    !
    !    14 June 2004
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) Z, the value at which the Sterling
    !    remainder is to be calculated.  Z must be positive.
    !
    !    Output, real ( kind = 8 ) DSTREM, the Sterling remainder.
    !
    Implicit None

    Integer, Parameter :: ncoef = 9

    Real ( kind = 8 ), Parameter, Dimension ( 0:ncoef ) :: coef = (/ &
         0.0D+00, &
         0.0833333333333333333333333333333D+00, &
         -0.00277777777777777777777777777778D+00, &
         0.000793650793650793650793650793651D+00, &
         -0.000595238095238095238095238095238D+00, &
         0.000841750841750841750841750841751D+00, &
         -0.00191752691752691752691752691753D+00, &
         0.00641025641025641025641025641026D+00, &
         -0.0295506535947712418300653594771D+00, &
         0.179644372368830573164938490016D+00 /)
    Real ( kind = 8 ) dstrem
    !    Real ( kind = 8 ) eval_pol
    !    Real ( kind = 8 ) gamma_log
    Real ( kind = 8 ), Parameter :: hln2pi = 0.91893853320467274178D+00
    Real ( kind = 8 ) sterl
    Real ( kind = 8 ) z

    If ( z <= 0.0D+00 ) Then
       Write ( *, '(a)' ) ' '
       Write ( *, '(a)' ) 'DSTREM - Fatal error!'
       Write ( *, '(a)' ) '  Zero or negative argument Z.'
       Stop
    End If

    If ( 6.0D+00 < z ) Then
       dstrem = eval_pol ( coef, ncoef, 1.0D+00 / z**2 ) * z
    Else
       sterl = hln2pi + ( z - 0.5D+00 ) * Log ( z ) - z
       dstrem = gamma_log ( z ) - sterl
    End If

    Return
  End Function dstrem

  Function dt1 ( p, q, df )

    !*****************************************************************************80
    !
    !! DT1 computes an approximate inverse of the cumulative T distribution.
    !
    !  Discussion:
    !
    !    This routine returns the inverse of the T distribution function, that is,
    !    the integral from 0 to INVT of the T density is P.  This is an
    !    initial approximation.
    !
    !    Thanks to Charles Katholi for pointing out that the RESHAPE
    !    function should not use a range in the "SHAPE" field (0:4,4),
    !    but simply the number of rows and columns (5,4), 04 May 2006.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) P, Q, the value whose inverse from the
    !    T distribution CDF is desired, and the value (1-P).
    !
    !    Input, real ( kind = 8 ) DF, the number of degrees of freedom of the
    !    T distribution.
    !
    !    Output, real ( kind = 8 ) DT1, the approximate value of X for which
    !    the T density CDF with DF degrees of freedom has value P.
    !
    Implicit None

    Real ( kind = 8 ), Dimension(0:4,4) :: coef = Reshape ( (/ &
         1.0D+00,     1.0D+00,    0.0D+00,   0.0D+00,  0.0D+00, &
         3.0D+00,    16.0D+00,    5.0D+00,   0.0D+00,  0.0D+00, &
         -15.0D+00,    17.0D+00,   19.0D+00,   3.0D+00,  0.0D+00, &
         -945.0D+00, -1920.0D+00, 1482.0D+00, 776.0D+00, 79.0D+00/), (/ 5, 4 /) )
    Real ( kind = 8 ), Parameter, Dimension ( 4 ) :: denom = (/ &
         4.0D+00, 96.0D+00, 384.0D+00, 92160.0D+00 /)
    Real ( kind = 8 ) denpow
    !    Real ( kind = 8 ) eval_pol
    Real ( kind = 8 ) df
    !    Real ( kind = 8 ) dinvnr
    Real ( kind = 8 ) dt1
    Integer i
    Integer, Parameter, Dimension ( 4 ) :: ideg = (/ 1, 2, 3, 4 /)
    Real ( kind = 8 ) p
    Real ( kind = 8 ) q
    Real ( kind = 8 ) sum1
    Real ( kind = 8 ) term
    Real ( kind = 8 ) x
    Real ( kind = 8 ) xp
    Real ( kind = 8 ) xx

    x = Abs ( dinvnr ( p, q ) )
    xx = x * x

    sum1 = x
    denpow = 1.0D+00
    Do i = 1, 4
       term = eval_pol ( coef(0,i), ideg(i), xx ) * x
       denpow = denpow * df
       sum1 = sum1 + term / ( denpow * denom(i) )
    End Do

    If ( 0.5D+00 <= p ) Then
       xp = sum1
    Else
       xp = -sum1
    End If

    dt1 = xp

    Return
  End Function dt1

  Subroutine dzror ( status, x, fx, xlo, xhi, qleft, qhi, is_entry )

    !*****************************************************************************80
    !
    !! DZROR seeks a zero of a function, using reverse communication.
    !
    !  Discussion:
    !
    !    This routine performs the zero finding.  STZROR must have been called
    !    before this routine in order to set its parameters.
    !
    !  Modified:
    !
    !    09 June 2004
    !
    !  Reference:
    !
    !    JCP Bus, TJ Dekker,
    !    Two Efficient Algorithms with Guaranteed Convergence for
    !    Finding a Zero of a Function,
    !    ACM Transactions on Mathematical Software,
    !    Volume 1, Number 4, pages 330-345, 1975.
    !
    !  Parameters:
    !
    !     STATUS <--> At the beginning of a zero finding problem, STATUS
    !                 should be set to 0 and ZROR invoked.  (The value
    !                 of other parameters will be ignored on this call.)
    !
    !                 When ZROR needs the function evaluated, it will set
    !                 STATUS to 1 and return.  The value of the function
    !                 should be set in FX and ZROR again called without
    !                 changing any of its other parameters.
    !
    !                 When ZROR has finished without error, it will return
    !                 with STATUS 0.  In that case (XLO,XHI) bound the answe
    !
    !                 If ZROR finds an error (which implies that F(XLO)-Y an
    !                 F(XHI)-Y have the same sign, it returns STATUS -1.  In
    !                 this case, XLO and XHI are undefined.
    !                         integer STATUS
    !
    !    Output, real ( kind = 8 ) X, the value of X at which F(X) is to
    !    be evaluated.
    !
    !    Input, real ( kind = 8 ) FX, the value of F(X), which must be calculated
    !    by the user when ZROR has returned on the previous call with STATUS = 1.
    !
    !    Output, real ( kind = 8 ) XLO, XHI, are lower and upper bounds for the
    !    solution when ZROR returns with STATUS = 0.
    !
    !    Output, logical QLEFT,is TRUE if the stepping search terminated
    !    unsucessfully at XLO.  If it is FALSE, the search terminated
    !    unsucessfully at XHI.
    !
    !    Output, logical QHI, is TRUE if Y < F(X) at the termination of the
    !    search and FALSE if F(X) < Y at the termination of the search.
    !
    Implicit None

    ! Timo Korhonen 2014: Remove the Entry dstzr
    LOGICAL, INTENT(IN) :: is_entry

    Real ( kind = 8 ), Save :: a
    Real ( kind = 8 ), Save :: abstol
    Real ( kind = 8 ), Save :: b ! save?
    Real ( kind = 8 ), Save :: c ! save?
    Real ( kind = 8 ), Save :: d ! save?
    Integer, Save :: ext
    Real ( kind = 8 ), Save :: fa
    Real ( kind = 8 ), Save :: fb
    Real ( kind = 8 ), Save :: fc ! save?
    Real ( kind = 8 ), Save :: fd ! save?
    Real ( kind = 8 ) fda
    Real ( kind = 8 ) fdb
    Logical, Save :: first
    !Timo: Real ( kind = 8 ) ftol
    Real ( kind = 8 ) fx
    Integer, Save :: i99999
    Real ( kind = 8 ), Save :: m ! save?
    Real ( kind = 8 ), Save :: mb ! save?
    Real ( kind = 8 ), Save :: p ! save?
    Real ( kind = 8 ), Save :: q ! save?
    Logical qhi
    Logical qleft
    Logical qrzero
    Real ( kind = 8 ), Save :: reltol
    Integer status
    Real ( kind = 8 ), Save :: tol ! save?
    Real ( kind = 8 ), Save :: w ! save?
    Real ( kind = 8 ) x
    Real ( kind = 8 ) xhi
    Real ( kind = 8 ) xlo
    Real ( kind = 8 ), Save :: xxhi = 0.0D+00
    Real ( kind = 8 ), Save :: xxlo = 0.0D+00
    ! Timo Korhonen 2014: Remove the Entry dstzr
    ! Real ( kind = 8 ) zabstl
    ! Real ( kind = 8 ) zreltl
    ! Real ( kind = 8 ) zxhi
    ! Real ( kind = 8 ) zxlo

    ! Save ! Gnu compiler says that this Save is redundant
    !Timo: Statement function ftol converted to an internal function
    !Timo: ftol(zx) = 0.5D+00 * Max ( abstol, reltol * Abs ( zx ) )

    ! Timo Korhonen 2014: Remove the Entry dstzr
    If ( is_entry ) Then
       Go to 250 ! Entry dstzr
    End If

    If ( 0 < status ) Then
       go to 280
    End If

    xlo = xxlo
    xhi = xxhi
    b = xlo
    x = xlo
    !
    !     GET-function-VALUE
    !
    ! Assign 10 to i99999
    i99999 = 10
    go to 270

10  Continue

    fb = fx
    xlo = xhi
    a = xlo
    x = xlo
    !
    !     GET-function-VALUE
    !
    ! Assign 20 to i99999
    i99999 = 20
    go to 270
    !
    !  Check that F(ZXLO) < 0 < F(ZXHI)  or F(ZXLO) > 0 > F(ZXHI)
    !
20  Continue

    If (.Not. ( fb < 0.0D+00 ) ) Then
       go to 40
    End If

    If (.Not. ( fx < 0.0D+00 ) ) Then
       go to 30
    End If

    status = -1
    qleft = fx < fb
    qhi = .False.
    Return

30  Continue
40  Continue

    If (.Not. ( 0.0D+00 < fb )) Then
       go to 60
    End If

    If (.Not. ( 0.0D+00 < fx )) Then
       go to 50
    End If

    status = -1
    qleft = ( fb < fx )
    qhi = .True.
    Return

50  Continue
60  Continue

    fa = fx
    first = .True.

70  Continue

    c = a
    fc = fa
    ext = 0

80  Continue

    If ( Abs ( fc ) < Abs ( fb ) ) Then

       If ( c == a ) Then
          d = a
          fd = fa
       End If

       a = b
       fa = fb
       xlo = c
       b = xlo
       fb = fc
       c = a
       fc = fa

    End If

    tol = ftol ( xlo )
    m = ( c + b ) * 0.5D+00
    mb = m - b

    If (.Not. ( tol < Abs ( mb ) ) ) Then
       go to 240
    End If

    If ( 3 < ext ) Then
       w = mb
       go to 190
    End If

! 110 Continue

    tol = Sign ( tol, mb )
    p = ( b - a ) * fb
    !
    !  I had to insert a rudimentary check on the divisions here
    !  to avoid ninny errors, JVB, 09 June 2004.
    !
    If ( first ) Then

       q = fa - fb
       first = .False.

    Else

       If ( d == b ) Then
          fdb = 1.0
       Else
          fdb = ( fd - fb ) / ( d - b )
       End If

       If ( d == a ) Then
          fda = 1.0
       Else
          fda = ( fd - fa ) / ( d - a )
       End If

       p = fda * p
       q = fdb * fa - fda * fb

    End If

! 130 Continue

    If ( p < 0.0D+00 ) Then
       p = -p
       q = -q
    End If

! 140 Continue

    If ( ext == 3 ) Then
       p = p *  2.0D+00
    End If

    If (.Not. ( ( p * 1.0D+00 ) == 0.0D+00 .Or. p <= ( q * tol ) ) ) Then
       go to 150
    End If

    w = tol
    go to 180

150 Continue

    If ( p < mb * q ) Then
       w = p / q
    Else
       w = mb
    End If

180 Continue
190 Continue

    d = a
    fd = fa
    a = b
    fa = fb
    b = b + w
    xlo = b
    x = xlo
    !
    !     GET-function-VALUE
    !
    ! Assign 200 to i99999
    i99999 = 200
    go to 270

200 Continue

    fb = fx

    If ( 0.0D+00 <= fc * fb ) Then

       go to 70

    Else

       If ( w == mb ) Then
          ext = 0
       Else
          ext = ext + 1
       End If

       go to 80

    End If

240 Continue

    xhi = c
    qrzero = ( 0.0D+00 <= fc .And. fb <= 0.0D+00 ) .Or. &
         ( fc < 0.0D+00 .And. fb >= 0.0D+00 )

    If ( qrzero ) Then
       status = 0
    Else
       status = -1
    End If

    Return

    ! Timo Korhonen 2014: Remove the Entry dstzr
    ! Entry dstzr ( zxlo, zxhi, zabstl, zreltl )
    ! Subroutine dzror ( status, x, fx, xlo, xhi, qleft, qhi, is_entry )
    ! Use (xlo=zxlo, xhi=zxhi, x=zabstl, fx=zreltl)
250 Continue ! Entry dstzr

    !*****************************************************************************80
    !
    !! DSTZR - SeT ZeRo finder - Reverse communication version
    !
    !  Discussion:
    !
    !    This routine sets quantities needed by ZROR.  The function of ZROR
    !    and the quantities set is given here.
    !
    !    Given a function F, find XLO such that F(XLO) = 0.
    !
    !     Input condition. F is a real ( kind = 8 ) function of a single
    !     real ( kind = 8 ) argument and XLO and XHI are such that
    !          F(XLO)*F(XHI)  <=  0.0
    !
    !     If the input condition is met, QRZERO returns .TRUE.
    !     and output values of XLO and XHI satisfy the following
    !          F(XLO)*F(XHI)  <= 0.
    !          ABS ( F(XLO) ) <= ABS ( F(XHI) )
    !          ABS ( XLO - XHI ) <= TOL(X)
    !     where
    !          TOL(X) = MAX ( ABSTOL, RELTOL * ABS ( X ) )
    !
    !     If this algorithm does not find XLO and XHI satisfying
    !     these conditions then QRZERO returns .FALSE.  This
    !     implies that the input condition was not met.
    !
    !  Reference:
    !
    !    JCP Bus, TJ Dekker,
    !    Two Efficient Algorithms with Guaranteed Convergence for
    !    Finding a Zero of a Function,
    !    ACM Transactions on Mathematical Software,
    !    Volume 1, Number 4, pages 330-345, 1975.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) XLO, XHI, the left and right endpoints of the
    !    interval to be searched for a solution.
    !
    !    Input, real ( kind = 8 ) ABSTOL, RELTOL, two numbers that determine
    !    the accuracy of the solution.
    !
    xxlo = xlo
    xxhi = xhi
    abstol = x
    reltol = fx
    !    xxlo = zxlo
    !    xxhi = zxhi
    !    abstol = zabstl
    !    reltol = zreltl
    Return
    !
    !     TO GET-function-VALUE
    !
270 status = 1
    Return

280 Continue
    ! go to i99999
    Select Case(i99999)
       Case(10)
          Go To 10
       Case(20)
          Go To 20
       Case(200)
          Go To 200
    End Select

  Contains
    Real ( kind = 8 ) Function ftol(zx)
      Real ( kind = 8 ), Intent(In) :: zx
      ftol = 0.5D+00 * Max ( abstol, reltol * Abs ( zx ) )
    End Function ftol

  End Subroutine dzror

!!$  Subroutine erf_values ( n_data, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! ERF_VALUES returns some values of the ERF or "error" function.
!!$    !
!!$    !  Definition:
!!$    !
!!$    !    ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    17 April 2001
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, real ( kind = 8 ) X, the argument of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 21
!!$
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.0000000000D+00, 0.1124629160D+00, 0.2227025892D+00, 0.3286267595D+00, &
!!$         0.4283923550D+00, 0.5204998778D+00, 0.6038560908D+00, 0.6778011938D+00, &
!!$         0.7421009647D+00, 0.7969082124D+00, 0.8427007929D+00, 0.8802050696D+00, &
!!$         0.9103139782D+00, 0.9340079449D+00, 0.9522851198D+00, 0.9661051465D+00, &
!!$         0.9763483833D+00, 0.9837904586D+00, 0.9890905016D+00, 0.9927904292D+00, &
!!$         0.9953222650D+00 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         0.0D+00, 0.1D+00, 0.2D+00, 0.3D+00, &
!!$         0.4D+00, 0.5D+00, 0.6D+00, 0.7D+00, &
!!$         0.8D+00, 0.9D+00, 1.0D+00, 1.1D+00, &
!!$         1.2D+00, 1.3D+00, 1.4D+00, 1.5D+00, &
!!$         1.6D+00, 1.7D+00, 1.8D+00, 1.9D+00, &
!!$         2.0D+00 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine erf_values

  Function error_f ( x )

    !*****************************************************************************80
    !
    !! ERROR_F evaluates the error function.
    !
    !  Discussion:
    !
    !    Since some compilers already supply a routine named ERF which evaluates
    !    the error function, this routine has been given a distinct, if
    !    somewhat unnatural, name.
    !
    !    The function is defined by:
    !
    !      ERF(X) = ( 2 / sqrt ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( - T**2 ) dT.
    !
    !    Properties of the function include:
    !
    !      Limit ( X -> -Infinity ) ERF(X) =          -1.0;
    !                               ERF(0) =           0.0;
    !                               ERF(0.476936...) = 0.5;
    !      Limit ( X -> +Infinity ) ERF(X) =          +1.0.
    !
    !      0.5D+00 * ( ERF(X/sqrt(2)) + 1 ) = Normal_01_CDF(X)
    !
    !  Modified:
    !
    !    17 November 2006
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument.
    !
    !    Output, real ( kind = 8 ) ERF, the value of the error function at X.
    !
    Implicit None

    Real ( kind = 8 ), Parameter, Dimension ( 5 ) :: a = (/ &
         0.771058495001320D-04, &
         -0.133733772997339D-02, &
         0.323076579225834D-01, &
         0.479137145607681D-01, &
         0.128379167095513D+00 /)
    Real ( kind = 8 ) ax
    Real ( kind = 8 ), Parameter, Dimension ( 3 ) :: b = (/ &
         0.301048631703895D-02, &
         0.538971687740286D-01, &
         0.375795757275549D+00 /)
    Real ( kind = 8 ) bot
    Real ( kind = 8 ), Parameter :: c = 0.564189583547756D+00
    Real ( kind = 8 ) error_f
    Real ( kind = 8 ), Dimension ( 8 ) :: p = (/   &
         -1.36864857382717D-07, 5.64195517478974D-01, &
         7.21175825088309D+00, 4.31622272220567D+01, &
         1.52989285046940D+02, 3.39320816734344D+02, &
         4.51918953711873D+02, 3.00459261020162D+02 /)
    Real ( kind = 8 ), Dimension ( 8 ) :: q = (/ &
         1.00000000000000D+00, 1.27827273196294D+01, &
         7.70001529352295D+01, 2.77585444743988D+02, &
         6.38980264465631D+02, 9.31354094850610D+02, &
         7.90950925327898D+02, 3.00459260956983D+02 /)
    Real ( kind = 8 ), Dimension ( 5 ) :: r = (/ &
         2.10144126479064D+00, 2.62370141675169D+01, &
         2.13688200555087D+01, 4.65807828718470D+00, &
         2.82094791773523D-01 /)
    Real ( kind = 8 ), Parameter, Dimension ( 4 ) :: s = (/ &
         9.41537750555460D+01, 1.87114811799590D+02, &
         9.90191814623914D+01, 1.80124575948747D+02 /)
    Real ( kind = 8 ) t
    Real ( kind = 8 ) top
    Real ( kind = 8 ) x
    Real ( kind = 8 ) x2

    ax = Abs ( x )

    If ( ax <= 0.5D+00 ) Then

       t = x * x

       top = (((( a(1)   * t &
            + a(2) ) * t &
            + a(3) ) * t &
            + a(4) ) * t &
            + a(5) ) + 1.0D+00

       bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0D+00
       error_f = ax * ( top / bot )

    Else If ( ax <= 4.0D+00 ) Then

       top = (((((( p(1)   * ax &
            + p(2) ) * ax &
            + p(3) ) * ax &
            + p(4) ) * ax &
            + p(5) ) * ax &
            + p(6) ) * ax &
            + p(7) ) * ax &
            + p(8)

       bot = (((((( q(1) * ax + q(2) ) * ax + q(3) ) * ax + q(4) ) * ax &
            + q(5) ) * ax + q(6) ) * ax + q(7) ) * ax + q(8)

       error_f = 0.5D+00 &
            + ( 0.5D+00 - Exp ( - x * x ) * top / bot )

    Else If ( ax < 5.8D+00 ) Then

       x2 = x * x
       t = 1.0D+00 / x2

       top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)

       bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t &
            + 1.0D+00

       error_f = ( c - top / ( x2 * bot )) / ax
       error_f = 0.5D+00 &
            + ( 0.5D+00 - Exp ( - x2 ) * error_f )

    Else

       error_f = 1.0D+00

    End If

    If ( x < 0.0D+00 ) Then
       error_f = -error_f
    End If

    Return
  End Function error_f

  Function error_fc ( ind, x )

    !*****************************************************************************80
    !
    !! ERROR_FC evaluates the complementary error function.
    !
    !  Modified:
    !
    !    09 December 1999
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, integer IND, chooses the scaling.
    !    If IND is nonzero, then the value returned has been multiplied by
    !    EXP(X*X).
    !
    !    Input, real ( kind = 8 ) X, the argument of the function.
    !
    !    Output, real ( kind = 8 ) ERROR_FC, the value of the complementary
    !    error function.
    !
    Implicit None

    Real ( kind = 8 ), Dimension ( 5 ) :: a = (/ &
         0.771058495001320D-04,  -0.133733772997339D-02, &
         0.323076579225834D-01,   0.479137145607681D-01, &
         0.128379167095513D+00 /)
    Real ( kind = 8 ) ax
    Real ( kind = 8 ), Dimension(3) :: b = (/ &
         0.301048631703895D-02, &
         0.538971687740286D-01, &
         0.375795757275549D+00 /)
    Real ( kind = 8 ) bot
    Real ( kind = 8 ), Parameter :: c = 0.564189583547756D+00
    Real ( kind = 8 ) e
    Real ( kind = 8 ) error_fc
    !    Real ( kind = 8 ) exparg
    Integer ind
    Real ( kind = 8 ), Dimension ( 8 ) :: p = (/ &
         -1.36864857382717D-07, 5.64195517478974D-01, &
         7.21175825088309D+00, 4.31622272220567D+01, &
         1.52989285046940D+02, 3.39320816734344D+02, &
         4.51918953711873D+02, 3.00459261020162D+02 /)
    Real ( kind = 8 ), Dimension ( 8 ) :: q = (/  &
         1.00000000000000D+00, 1.27827273196294D+01, &
         7.70001529352295D+01, 2.77585444743988D+02, &
         6.38980264465631D+02, 9.31354094850610D+02, &
         7.90950925327898D+02, 3.00459260956983D+02 /)
    Real ( kind = 8 ), Dimension ( 5 ) :: r = (/ &
         2.10144126479064D+00, 2.62370141675169D+01, &
         2.13688200555087D+01, 4.65807828718470D+00, &
         2.82094791773523D-01 /)
    Real ( kind = 8 ), Dimension ( 4 ) :: s = (/ &
         9.41537750555460D+01, 1.87114811799590D+02, &
         9.90191814623914D+01, 1.80124575948747D+02 /)
    Real ( kind = 8 ) t
    Real ( kind = 8 ) top
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    !
    !  ABS ( X ) <= 0.5
    !
    ax = Abs ( x )

    If ( ax <= 0.5D+00 ) Then

       t = x * x

       top = (((( a(1) * t + a(2) ) * t + a(3) ) * t + a(4) ) * t + a(5) ) &
            + 1.0D+00

       bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0D+00

       error_fc = 0.5D+00 + ( 0.5D+00 &
            - x * ( top / bot ) )

       If ( ind /= 0 ) Then
          error_fc = Exp ( t ) * error_fc
       End If

       Return

    End If
    !
    !  0.5 < abs ( X ) <= 4
    !
    If ( ax <= 4.0D+00 ) Then

       top = (((((( p(1) * ax + p(2)) * ax + p(3)) * ax + p(4)) * ax &
            + p(5)) * ax + p(6)) * ax + p(7)) * ax + p(8)

       bot = (((((( q(1) * ax + q(2)) * ax + q(3)) * ax + q(4)) * ax &
            + q(5)) * ax + q(6)) * ax + q(7)) * ax + q(8)

       error_fc = top / bot
       !
       !  4 < ABS ( X )
       !
    Else

       If ( x <= -5.6D+00 ) Then

          If ( ind == 0 ) Then
             error_fc =  2.0D+00
          Else
             error_fc =  2.0D+00  * Exp ( x * x )
          End If

          Return

       End If

       If ( ind == 0 ) Then

          If ( 100.0D+00 < x ) Then
             error_fc = 0.0D+00
             Return
          End If

          If ( -exparg ( 1 ) < x * x ) Then
             error_fc = 0.0D+00
             Return
          End If

       End If

       t = ( 1.0D+00 / x )**2

       top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)

       bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t &
            + 1.0D+00

       error_fc = ( c - t * top / bot ) / ax

    End If
    !
    !  Final assembly.
    !
    If ( ind /= 0 ) Then

       If ( x < 0.0D+00 ) Then
          error_fc =  2.0D+00  * Exp ( x * x ) - error_fc
       End If

    Else

       w = x * x
       t = w
       e = w - t
       error_fc = (( 0.5D+00 &
            + ( 0.5D+00 - e ) ) * Exp ( - t ) ) * error_fc

       If ( x < 0.0D+00 ) Then
          error_fc =  2.0D+00  - error_fc
       End If

    End If

    Return
  End Function error_fc

  Function esum ( mu, x )

    !*****************************************************************************80
    !
    !! ESUM evaluates exp ( MU + X ).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, integer MU, part of the argument.
    !
    !    Input, real ( kind = 8 ) X, part of the argument.
    !
    !    Output, real ( kind = 8 ) ESUM, the value of exp ( MU + X ).
    !
    Implicit None

    Real ( kind = 8 ) esum
    Integer mu
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x

    If ( x <= 0.0D+00 ) Then
       If ( 0 <= mu ) Then
          w = mu + x
          If ( w <= 0.0D+00 ) Then
             esum = Exp ( w )
             Return
          End If
       End If
    Else If ( 0.0D+00 < x ) Then
       If ( mu <= 0 ) Then
          w = mu + x
          If ( 0.0D+00 <= w ) Then
             esum = Exp ( w )
             Return
          End If
       End If
    End If

    w = mu
    esum = Exp ( w ) * Exp ( x )

    Return
  End Function esum

  Function eval_pol ( a, n, x )

    !*****************************************************************************80
    !
    !! EVAL_POL evaluates a polynomial at X.
    !
    !  Discussion:
    !
    !    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
    !
    !  Modified:
    !
    !    15 December 1999
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A(0:N), coefficients of the polynomial.
    !
    !    Input, integer N, length of A.
    !
    !    Input, real ( kind = 8 ) X, the point at which the polynomial
    !    is to be evaluated.
    !
    !    Output, real ( kind = 8 ) EVAL_POL, the value of the polynomial at X.
    !
    Implicit None

    Integer n

    Real ( kind = 8 ) a(0:n)
    Real ( kind = 8 ) eval_pol
    Integer i
    Real ( kind = 8 ) term
    Real ( kind = 8 ) x

    term = a(n)
    Do i = n - 1, 0, -1
       term = term * x + a(i)
    End Do

    eval_pol = term

    Return
  End Function eval_pol

  Function exparg ( l )

    !*****************************************************************************80
    !
    !! EXPARG returns the largest or smallest legal argument for EXP.
    !
    !  Discussion:
    !
    !    Only an approximate limit for the argument of EXP is desired.
    !
    !  Modified:
    !
    !    09 December 1999
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, integer L, indicates which limit is desired.
    !    If L = 0, then the largest positive argument for EXP is desired.
    !    Otherwise, the largest negative argument for EXP for which the
    !    result is nonzero is desired.
    !
    !    Output, real ( kind = 8 ) EXPARG, the desired value.
    !
    Implicit None

    Integer b
    Real ( kind = 8 ) exparg
    !    Integer ipmpar
    Integer l
    Real ( kind = 8 ) lnb
    Integer m
    !
    !  Get the arithmetic base.
    !
    b = ipmpar(4)
    !
    !  Compute the logarithm of the arithmetic base.
    !
    If ( b == 2 ) Then
       lnb = 0.69314718055995D+00
    Else If ( b == 8 ) Then
       lnb = 2.0794415416798D+00
    Else If ( b == 16 ) Then
       lnb = 2.7725887222398D+00
    Else
       lnb = Log ( Real ( b, kind = 8 ) )
    End If

    If ( l /= 0 ) Then
       m = ipmpar(9) - 1
       exparg = 0.99999D+00 * ( m * lnb )
    Else
       m = ipmpar(10)
       exparg = 0.99999D+00 * ( m * lnb )
    End If

    Return
  End Function exparg

!!$  Subroutine f_cdf_values ( n_data, a, b, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! F_CDF_VALUES returns some values of the F CDF test function.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    The value of F_CDF ( DFN, DFD, X ) can be evaluated in Mathematica by
!!$    !    commands like:
!!$    !
!!$    !      Needs["Statistics`ContinuousDistributions`"]
!!$    !      CDF[FRatioDistribution[ DFN, DFD ], X ]
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    11 June 2004
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !    Stephen Wolfram,
!!$    !    The Mathematica Book,
!!$    !    Fourth Edition,
!!$    !    Wolfram Media / Cambridge University Press, 1999.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, integer A, integer B, real ( kind = 8 ) X, the arguments
!!$    !    of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 20
!!$
!!$    Integer a
!!$    Integer, Save, Dimension ( n_max ) :: a_vec = (/ &
!!$         1, 1, 5, 1, &
!!$         2, 4, 1, 6, &
!!$         8, 1, 3, 6, &
!!$         1, 1, 1, 1, &
!!$         2, 3, 4, 5 /)
!!$    Integer b
!!$    Integer, Save, Dimension ( n_max ) :: b_vec = (/ &
!!$         1,  5,  1,  5, &
!!$         10, 20,  5,  6, &
!!$         16,  5, 10, 12, &
!!$         5,  5,  5,  5, &
!!$         5,  5,  5,  5 /)
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.500000D+00, 0.499971D+00, 0.499603D+00, 0.749699D+00, &
!!$         0.750466D+00, 0.751416D+00, 0.899987D+00, 0.899713D+00, &
!!$         0.900285D+00, 0.950025D+00, 0.950057D+00, 0.950193D+00, &
!!$         0.975013D+00, 0.990002D+00, 0.994998D+00, 0.999000D+00, &
!!$         0.568799D+00, 0.535145D+00, 0.514343D+00, 0.500000D+00 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         1.00D+00,  0.528D+00, 1.89D+00,  1.69D+00, &
!!$         1.60D+00,  1.47D+00,  4.06D+00,  3.05D+00, &
!!$         2.09D+00,  6.61D+00,  3.71D+00,  3.00D+00, &
!!$         10.01D+00, 16.26D+00, 22.78D+00, 47.18D+00, &
!!$         1.00D+00,  1.00D+00,  1.00D+00,  1.00D+00 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       a = 0
!!$       b = 0
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       a = a_vec(n_data)
!!$       b = b_vec(n_data)
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine f_cdf_values

!!$  Subroutine f_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    The value of NONCENTRAL_F_CDF ( DFN, DFD, LAMDA, X ) can be evaluated
!!$    !    in Mathematica by commands like:
!!$    !
!!$    !      Needs["Statistics`ContinuousDistributions`"]
!!$    !      CDF[NoncentralFRatioDistribution[ DFN, DFD, LAMBDA ], X ]
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    12 June 2004
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !    Stephen Wolfram,
!!$    !    The Mathematica Book,
!!$    !    Fourth Edition,
!!$    !    Wolfram Media / Cambridge University Press, 1999.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, integer A, integer B, real ( kind = 8 ) LAMBDA, the
!!$    !    parameters of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) X, the argument of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 22
!!$
!!$    Integer a
!!$    Integer, Save, Dimension ( n_max ) :: a_vec = (/ &
!!$         1,  1,  1,  1, &
!!$         1,  1,  1,  1, &
!!$         1,  1,  2,  2, &
!!$         3,  3,  4,  4, &
!!$         5,  5,  6,  6, &
!!$         8, 16 /)
!!$    Integer b
!!$    Integer, Save, Dimension ( n_max ) :: b_vec = (/ &
!!$         1,  5,  5,  5, &
!!$         5,  5,  5,  5, &
!!$         5,  5,  5, 10, &
!!$         5,  5,  5,  5, &
!!$         1,  5,  6, 12, &
!!$         16,  8 /)
!!$    Real ( kind = 8 ) fx
!!$    Real, Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.500000E+00, 0.636783E+00, 0.584092E+00, 0.323443E+00, &
!!$         0.450119E+00, 0.607888E+00, 0.705928E+00, 0.772178E+00, &
!!$         0.819105E+00, 0.317035E+00, 0.432722E+00, 0.450270E+00, &
!!$         0.426188E+00, 0.337744E+00, 0.422911E+00, 0.692767E+00, &
!!$         0.363217E+00, 0.421005E+00, 0.426667E+00, 0.446402E+00, &
!!$         0.844589E+00, 0.816368E+00 /)
!!$    Real ( kind = 8 ) lambda
!!$    Real, Save, Dimension ( n_max ) :: lambda_vec = (/ &
!!$         0.00E+00,  0.000E+00, 0.25E+00,  1.00E+00, &
!!$         1.00E+00,  1.00E+00,  1.00E+00,  1.00E+00, &
!!$         1.00E+00,  2.00E+00,  1.00E+00,  1.00E+00, &
!!$         1.00E+00,  2.00E+00,  1.00E+00,  1.00E+00, &
!!$         0.00E+00,  1.00E+00,  1.00E+00,  1.00E+00, &
!!$         1.00E+00,  1.00E+00 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real, Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         1.00E+00,  1.00E+00, 1.00E+00,  0.50E+00, &
!!$         1.00E+00,  2.00E+00, 3.00E+00,  4.00E+00, &
!!$         5.00E+00,  1.00E+00, 1.00E+00,  1.00E+00, &
!!$         1.00E+00,  1.00E+00, 1.00E+00,  2.00E+00, &
!!$         1.00E+00,  1.00E+00, 1.00E+00,  1.00E+00, &
!!$         2.00E+00,  2.00E+00 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       a = 0
!!$       b = 0
!!$       lambda = 0.0D+00
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       a = a_vec(n_data)
!!$       b = b_vec(n_data)
!!$       lambda = lambda_vec(n_data)
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine f_noncentral_cdf_values

  Function fpser ( a, b, x, eps )

    !*****************************************************************************80
    !
    !! FPSER evaluates IX(A,B)(X) for very small B.
    !
    !  Discussion:
    !
    !    This routine is appropriate for use when
    !
    !      B < min ( EPS, EPS * A )
    !
    !    and
    !
    !      X <= 0.5.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, parameters of the function.
    !
    !    Input, real ( kind = 8 ) X, the point at which the function is to
    !    be evaluated.
    !
    !    Input, real ( kind = 8 ) EPS, a tolerance.
    !
    !    Output, real ( kind = 8 ) FPSER, the value of IX(A,B)(X).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) an
    Real ( kind = 8 ) b
    Real ( kind = 8 ) c
    Real ( kind = 8 ) eps
    !    Real ( kind = 8 ) exparg
    Real ( kind = 8 ) fpser
    Real ( kind = 8 ) s
    Real ( kind = 8 ) t
    Real ( kind = 8 ) tol
    Real ( kind = 8 ) x

    fpser = 1.0D+00

    If ( 1.0D-03 * eps < a ) Then
       fpser = 0.0D+00
       t = a * Log ( x )
       If ( t < exparg ( 1 ) ) Then
          Return
       End If
       fpser = Exp ( t )
    End If
    !
    !  1/B(A,B) = B
    !
    fpser = ( b / a ) * fpser
    tol = eps / a
    an = a + 1.0D+00
    t = x
    s = t / an

    Do

       an = an + 1.0D+00
       t = x * t
       c = t / an
       s = s + c

       If ( Abs ( c ) <= tol ) Then
          Exit
       End If

    End Do

    fpser = fpser * ( 1.0D+00 + a * s )

    Return
  End Function fpser

  Function gam1 ( a )

    !*****************************************************************************80
    !
    !! GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5 <= A <= 1.5
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, forms the argument of the Gamma function.
    !
    !    Output, real ( kind = 8 ) GAM1, the value of 1 / GAMMA ( A + 1 ) - 1.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) bot
    Real ( kind = 8 ) d
    Real ( kind = 8 ) gam1
    Real ( kind = 8 ), Parameter, Dimension ( 7 ) :: p = (/ &
         0.577215664901533D+00, -0.409078193005776D+00, &
         -0.230975380857675D+00,  0.597275330452234D-01, &
         0.766968181649490D-02, -0.514889771323592D-02, &
         0.589597428611429D-03 /)
    Real ( kind = 8 ), Dimension ( 5 ) :: q = (/ &
         0.100000000000000D+01, 0.427569613095214D+00, &
         0.158451672430138D+00, 0.261132021441447D-01, &
         0.423244297896961D-02 /)
    Real ( kind = 8 ), Dimension ( 9 ) :: r = (/ &
         -0.422784335098468D+00, -0.771330383816272D+00, &
         -0.244757765222226D+00,  0.118378989872749D+00, &
         0.930357293360349D-03, -0.118290993445146D-01, &
         0.223047661158249D-02,  0.266505979058923D-03, &
         -0.132674909766242D-03 /)
    Real ( kind = 8 ), Parameter :: s1 = 0.273076135303957D+00
    Real ( kind = 8 ), Parameter :: s2 = 0.559398236957378D-01
    Real ( kind = 8 ) t
    Real ( kind = 8 ) top
    Real ( kind = 8 ) w

    d = a - 0.5D+00

    If ( 0.0D+00 < d ) Then
       t = d - 0.5D+00
    Else
       t = a
    End If

    If ( t == 0.0D+00 ) Then

       gam1 = 0.0D+00

    Else If ( 0.0D+00 < t ) Then

       top = (((((    &
            p(7)   &
            * t + p(6) ) &
            * t + p(5) ) &
            * t + p(4) ) &
            * t + p(3) ) &
            * t + p(2) ) &
            * t + p(1)

       bot = ((( q(5) * t + q(4) ) * t + q(3) ) * t + q(2) ) * t &
            + 1.0D+00

       w = top / bot

       If ( d <= 0.0D+00 ) Then
          gam1 = a * w
       Else
          gam1 = ( t / a ) * ( ( w - 0.5D+00 ) &
               - 0.5D+00 )
       End If

    Else If ( t < 0.0D+00 ) Then

       top = (((((((  &
            r(9)   &
            * t + r(8) ) &
            * t + r(7) ) &
            * t + r(6) ) &
            * t + r(5) ) &
            * t + r(4) ) &
            * t + r(3) ) &
            * t + r(2) ) &
            * t + r(1)

       bot = ( s2 * t + s1 ) * t + 1.0D+00
       w = top / bot

       If ( d <= 0.0D+00 ) Then
          gam1 = a * ( ( w + 0.5D+00 ) + 0.5D+00 )
       Else
          gam1 = t * w / a
       End If

    End If

    Return
  End Function gam1

  Function gamma_ieva ( a )

    !*****************************************************************************80
    !
    !! GAMMA evaluates the gamma function.
    !
    !  Author:
    !
    !    Alfred Morris,
    !    Naval Surface Weapons Center,
    !    Dahlgren, Virginia.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, the argument of the Gamma function.
    !
    !    Output, real ( kind = 8 ) GAMMA, the value of the Gamma function.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) bot
    Real ( kind = 8 ), Parameter :: d = 0.41893853320467274178D+00
    !    Real ( kind = 8 ) exparg
    Real ( kind = 8 ) g
    Real ( kind = 8 ) :: gamma_ieva
    Integer i
    Integer j
    Real ( kind = 8 ) lnx
    Integer m
    Integer n
    Real ( kind = 8 ), Dimension ( 7 ) :: p = (/ &
         0.539637273585445D-03, 0.261939260042690D-02, &
         0.204493667594920D-01, 0.730981088720487D-01, &
         0.279648642639792D+00, 0.553413866010467D+00, &
         1.0D+00 /)
    Real ( kind = 8 ), Parameter :: pi = 3.1415926535898D+00
    Real ( kind = 8 ), Dimension ( 7 ) :: q = (/ &
         -0.832979206704073D-03,  0.470059485860584D-02, &
         0.225211131035340D-01, -0.170458969313360D+00, &
         -0.567902761974940D-01,  0.113062953091122D+01, &
         1.0D+00 /)
    Real ( kind = 8 ), Parameter :: r1 =  0.820756370353826D-03
    Real ( kind = 8 ), Parameter :: r2 = -0.595156336428591D-03
    Real ( kind = 8 ), Parameter :: r3 =  0.793650663183693D-03
    Real ( kind = 8 ), Parameter :: r4 = -0.277777777770481D-02
    Real ( kind = 8 ), Parameter :: r5 =  0.833333333333333D-01
    Real ( kind = 8 ) s
    Real ( kind = 8 ) t
    Real ( kind = 8 ) top
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    Real ( kind = 8 ) z

    gamma_ieva = 0.0D+00
    x = a

    If ( Abs ( a ) < 15.0D+00 ) Then
       !
       !  Evaluation of GAMMA(A) for |A| < 15
       !
       t = 1.0D+00
       m = Int ( a ) - 1
       !
       !  Let T be the product of A-J when 2 <= A.
       !
       If ( 0 <= m ) Then

          Do j = 1, m
             x = x - 1.0D+00
             t = x * t
          End Do

          x = x - 1.0D+00
          !
          !  Let T be the product of A+J WHEN A < 1
          !
       Else

          t = a

          If ( a <= 0.0D+00 ) Then

             m = - m - 1

             Do j = 1, m
                x = x + 1.0D+00
                t = x * t
             End Do

             x = ( x + 0.5D+00 ) + 0.5D+00
             t = x * t
             If ( t == 0.0D+00 ) Then
                Return
             End If

          End If
          !
          !  Check if 1/T can overflow.
          !
          If ( Abs ( t ) < 1.0D-30 ) Then
             If ( 1.0001D+00 < Abs ( t ) * Huge ( t ) ) Then
                gamma_ieva = 1.0D+00 / t
             End If
             Return
          End If

       End If
       !
       !  Compute Gamma(1 + X) for 0 <= X < 1.
       !
       top = p(1)
       bot = q(1)
       Do i = 2, 7
          top = top * x + p(i)
          bot = bot * x + q(i)
       End Do

       gamma_ieva = top / bot
       !
       !  Termination.
       !
       If ( 1.0D+00 <= a ) Then
          gamma_ieva = gamma_ieva * t
       Else
          gamma_ieva = gamma_ieva / t
       End If
       !
       !  Evaluation of Gamma(A) FOR 15 <= ABS ( A ).
       !
    Else

       If ( 1000.0D+00 <= Abs ( a ) ) Then
          Return
       End If

       If ( a <= 0.0D+00 ) Then

          x = -a
          n = INT(x)
          t = x - n

          If ( 0.9D+00 < t ) Then
             t = 1.0D+00 - t
          End If

          s = Sin ( pi * t ) / pi

          If ( Mod ( n, 2 ) == 0 ) Then
             s = -s
          End If

          If ( s == 0.0D+00 ) Then
             Return
          End If

       End If
       !
       !  Compute the modified asymptotic sum.
       !
       t = 1.0D+00 / ( x * x )

       g = (((( r1 * t + r2 ) * t + r3 ) * t + r4 ) * t + r5 ) / x

       lnx = Log ( x )
       !
       !  Final assembly.
       !
       z = x
       g = ( d + g ) + ( z - 0.5D+00 ) &
            * ( lnx - 1.0D+00 )
       w = g
       t = g - Real ( w, kind = 8 )

       If ( 0.99999D+00 * exparg ( 0 ) < w ) Then
          Return
       End If

       gamma_ieva = Exp ( w )* ( 1.0D+00 + t )

       If ( a < 0.0D+00 ) Then
          gamma_ieva = ( 1.0D+00 / ( gamma_ieva * s ) ) / x
       End If

    End If

    Return
  End Function gamma_ieva

  Subroutine gamma_inc ( a, x, ans, qans, ind )

    !*****************************************************************************80
    !
    !! GAMMA_INC evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
    !
    !  Modified:
    !
    !    16 April 2005
    !
    !  Author:
    !
    !    Alfred Morris,
    !    Naval Surface Weapons Center,
    !    Dahlgren, Virginia.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, X, the arguments of the incomplete
    !    gamma ratio.  A and X must be nonnegative.  A and X cannot
    !    both be zero.
    !
    !    Output, real ( kind = 8 ) ANS, QANS.  On normal output,
    !    ANS = P(A,X) and QANS = Q(A,X).  However, ANS is set to 2 if
    !    A or X is negative, or both are 0, or when the answer is
    !    computationally indeterminate because A is extremely large
    !    and X is very close to A.
    !
    !    Input, integer IND, indicates the accuracy request:
    !    0, as much accuracy as possible.
    !    1, to within 1 unit of the 6-th significant digit,
    !    otherwise, to within 1 unit of the 3rd significant digit.
    !
    !  Local Parameters:
    !
    !     ALOG10 = LN(10)
    !     RT2PIN = 1/SQRT(2*PI)
    !     RTPI   = SQRT(PI)
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) a2n
    Real ( kind = 8 ) a2nm1
    Real ( kind = 8 ) acc
    Real ( kind = 8 ), Dimension ( 3 ) :: acc0 = (/ &
         5.0D-15, 5.0D-07, 5.0D-04 /)
    Real ( kind = 8 ), Parameter :: alog10 = 2.30258509299405D+00
    Real ( kind = 8 ) am0
    Real ( kind = 8 ) amn
    Real ( kind = 8 ) an
    Real ( kind = 8 ) an0
    Real ( kind = 8 ) ans
    Real ( kind = 8 ) apn
    Real ( kind = 8 ) b2n
    Real ( kind = 8 ) b2nm1
    Real ( kind = 8 ) big(3)
    Real ( kind = 8 ) c
    Real ( kind = 8 ) c0
    Real ( kind = 8 ) c1
    Real ( kind = 8 ) c2
    Real ( kind = 8 ) c3
    Real ( kind = 8 ) c4
    Real ( kind = 8 ) c5
    Real ( kind = 8 ) c6
    Real ( kind = 8 ) cma
    Real ( kind = 8 ) d0(13)
    Real ( kind = 8 ) d1(12)
    Real ( kind = 8 ) d2(10)
    Real ( kind = 8 ) d3(8)
    Real ( kind = 8 ) d4(6)
    Real ( kind = 8 ) d5(4)
    Real ( kind = 8 ) d6(2)
    Real ( kind = 8 ) d10
    Real ( kind = 8 ) d20
    Real ( kind = 8 ) d30
    Real ( kind = 8 ) d40
    Real ( kind = 8 ) d50
    Real ( kind = 8 ) d60
    Real ( kind = 8 ) d70
    Real ( kind = 8 ) e
    Real ( kind = 8 ) e0
    Real ( kind = 8 ) e00(3)
    !    Real ( kind = 8 ) error_f
    !    Real ( kind = 8 ) error_fc
    Real ( kind = 8 ) g
    !    Real ( kind = 8 ) gam1
    !    Real ( kind = 8 ) gamma
    Real ( kind = 8 ) h
    Integer i
    Integer ind
    Integer iop
    Real ( kind = 8 ) j
    Real ( kind = 8 ) l
    Integer m
    Integer n
    Integer n_max
    Real ( kind = 8 ) qans
    Real ( kind = 8 ) r
    !    Real ( kind = 8 ) rexp
    !    Real ( kind = 8 ) rlog
    Real ( kind = 8 ), Parameter :: rt2pin = 0.398942280401433D+00
    Real ( kind = 8 ) rta
    Real ( kind = 8 ), Parameter :: rtpi = 1.77245385090552D+00
    Real ( kind = 8 ) rtx
    Real ( kind = 8 ) s
    Real ( kind = 8 ) sum1
    Real ( kind = 8 ) t
    Real ( kind = 8 ) t1
    Real ( kind = 8 ) tol
    Real ( kind = 8 ) twoa
    Real ( kind = 8 ) u
    Real ( kind = 8 ) w
    Real ( kind = 8 ) wk(20)
    Real ( kind = 8 ) x
    Real ( kind = 8 ) x0
    Real ( kind = 8 ) x00(3)
    Real ( kind = 8 ) y
    Real ( kind = 8 ) z

    Data big(1)/20.0D+00/,big(2)/14.0D+00/,big(3)/10.0D+00/
    Data e00(1)/0.25D-03/,e00(2)/0.25D-01/,e00(3)/0.14D+00/
    Data x00(1)/31.0D+00/,x00(2)/17.0D+00/,x00(3)/9.7D+00/
    Data d0(1)/0.833333333333333D-01/
    Data d0(2)/-0.148148148148148D-01/
    Data d0(3)/0.115740740740741D-02/,d0(4)/0.352733686067019D-03/
    Data d0(5)/-0.178755144032922D-03/,d0(6)/0.391926317852244D-04/
    Data d0(7)/-0.218544851067999D-05/,d0(8)/-0.185406221071516D-05/
    Data d0(9)/0.829671134095309D-06/,d0(10)/-0.176659527368261D-06/
    Data d0(11)/0.670785354340150D-08/,d0(12)/0.102618097842403D-07/
    Data d0(13)/-0.438203601845335D-08/
    Data d10/-0.185185185185185D-02/,d1(1)/-0.347222222222222D-02/
    Data d1(2)/0.264550264550265D-02/,d1(3)/-0.990226337448560D-03/
    Data d1(4)/0.205761316872428D-03/,d1(5)/-0.401877572016461D-06/
    Data d1(6)/-0.180985503344900D-04/,d1(7)/0.764916091608111D-05/
    Data d1(8)/-0.161209008945634D-05/,d1(9)/0.464712780280743D-08/
    Data d1(10)/0.137863344691572D-06/,d1(11)/-0.575254560351770D-07/
    Data d1(12)/0.119516285997781D-07/
    Data d20/0.413359788359788D-02/,d2(1)/-0.268132716049383D-02/
    Data d2(2)/0.771604938271605D-03/,d2(3)/0.200938786008230D-05/
    Data d2(4)/-0.107366532263652D-03/,d2(5)/0.529234488291201D-04/
    Data d2(6)/-0.127606351886187D-04/,d2(7)/0.342357873409614D-07/
    Data d2(8)/0.137219573090629D-05/,d2(9)/-0.629899213838006D-06/
    Data d2(10)/0.142806142060642D-06/
    Data d30/0.649434156378601D-03/,d3(1)/0.229472093621399D-03/
    Data d3(2)/-0.469189494395256D-03/,d3(3)/0.267720632062839D-03/
    Data d3(4)/-0.756180167188398D-04/,d3(5)/-0.239650511386730D-06/
    Data d3(6)/0.110826541153473D-04/,d3(7)/-0.567495282699160D-05/
    Data d3(8)/0.142309007324359D-05/
    Data d40/-0.861888290916712D-03/,d4(1)/0.784039221720067D-03/
    Data d4(2)/-0.299072480303190D-03/,d4(3)/-0.146384525788434D-05/
    Data d4(4)/0.664149821546512D-04/,d4(5)/-0.396836504717943D-04/
    Data d4(6)/0.113757269706784D-04/
    Data d50/-0.336798553366358D-03/,d5(1)/-0.697281375836586D-04/
    Data d5(2)/0.277275324495939D-03/,d5(3)/-0.199325705161888D-03/
    Data d5(4)/0.679778047793721D-04/
    Data d60/0.531307936463992D-03/,d6(1)/-0.592166437353694D-03/
    Data d6(2)/0.270878209671804D-03/
    Data d70 / 0.344367606892378D-03/

    e = Epsilon ( 1.0D+00 )

    If ( a < 0.0D+00 .Or. x < 0.0D+00 ) Then
       ans = 2.0D+00
       Return
    End If

    If ( a == 0.0D+00 .And. x == 0.0D+00 ) Then
       ans = 2.0D+00
       Return
    End If

    If ( a * x == 0.0D+00 ) Then
       If ( x <= a ) Then
          ans = 0.0D+00
          qans = 1.0D+00
       Else
          ans = 1.0D+00
          qans = 0.0D+00
       End If
       Return
    End If

    iop = ind + 1
    If ( iop /= 1 .And. iop /= 2 ) iop = 3
    acc = Max ( acc0(iop), e )
    e0 = e00(iop)
    x0 = x00(iop)
    !
    !  Select the appropriate algorithm.
    !
    If ( 1.0D+00 <= a ) Then
       go to 10
    End If

    If ( a == 0.5D+00 ) Then
       go to 390
    End If

    If ( x < 1.1D+00 ) Then
       go to 160
    End If

    t1 = a * Log ( x ) - x
    u = a * Exp ( t1 )

    If ( u == 0.0D+00 ) Then
       ans = 1.0D+00
       qans = 0.0D+00
       Return
    End If

    r = u * ( 1.0D+00 + gam1 ( a ) )
    go to 250

10  Continue

    If ( big(iop) <= a ) Then
       go to 30
    End If

    If ( x < a .Or. x0 <= x ) Then
       go to 20
    End If

    twoa = a + a
    m = Int ( twoa )

    If ( twoa == Real ( m, kind = 8 ) ) Then
       i = m / 2
       If ( a == Real ( i, kind = 8 ) ) Then
          go to 210
       End If
       go to 220
    End If

20  Continue

    t1 = a * Log ( x ) - x
    r = Exp ( t1 ) / gamma ( a )
    go to 40

30  Continue

    l = x / a

    If ( l == 0.0D+00 ) Then
       ans = 0.0D+00
       qans = 1.0D+00
       Return
    End If

    s = 0.5D+00 + ( 0.5D+00 - l )
    z = rlog ( l )
    If ( 700.0D+00 / a <= z ) Then
       go to 410
    End If

    y = a * z
    rta = Sqrt ( a )

    If ( Abs ( s ) <= e0 / rta ) Then
       go to 330
    End If

    If ( Abs ( s ) <= 0.4D+00 ) Then
       go to 270
    End If

    t = ( 1.0D+00 / a )**2
    t1 = ((( 0.75D+00 * t - 1.0D+00 ) * t + 3.5D+00 ) &
         * t - 105.0D+00 ) / ( a * 1260.0D+00 )
    t1 = t1 - y
    r = rt2pin * rta * Exp ( t1 )

40  Continue

    If ( r == 0.0D+00 ) Then
       If ( x <= a ) Then
          ans = 0.0D+00
          qans = 1.0D+00
       Else
          ans = 1.0D+00
          qans = 0.0D+00
       End If
       Return
    End If

    If ( x <= Max ( a, alog10 ) ) Then
       go to 50
    End If

    If ( x < x0 ) Then
       go to 250
    End If

    go to 100
    !
    !  Taylor series for P/R.
    !
50  Continue

    apn = a + 1.0D+00
    t = x / apn
    wk(1) = t

    n = 20

    Do i = 2, 20
       apn = apn + 1.0D+00
       t = t * ( x / apn )
       If ( t <= 1.0D-03 ) Then
          n = i
          Exit
       End If
       wk(i) = t
    End Do

    sum1 = t

    tol = 0.5D+00 * acc

    Do

       apn = apn + 1.0D+00
       t = t * ( x / apn )
       sum1 = sum1 + t

       If ( t <= tol ) Then
          Exit
       End If

    End Do

    n_max = n - 1
    Do m = 1, n_max
       n = n - 1
       sum1 = sum1 + wk(n)
    End Do

    ans = ( r / a ) * ( 1.0D+00 + sum1 )
    qans = 0.5D+00 + ( 0.5D+00 - ans )
    Return
    !
    !  Asymptotic expansion.
    !
100 Continue

    amn = a - 1.0D+00
    t = amn / x
    wk(1) = t

    n = 20

    Do i = 2, 20
       amn = amn - 1.0D+00
       t = t * ( amn / x )
       If ( Abs ( t ) <= 1.0D-03 ) Then
          n = i
          Exit
       End If
       wk(i) = t
    End Do

    sum1 = t

    Do

       If ( Abs ( t ) <= acc ) Then
          Exit
       End If

       amn = amn - 1.0D+00
       t = t * ( amn / x )
       sum1 = sum1 + t

    End Do

    n_max = n - 1
    Do m = 1, n_max
       n = n - 1
       sum1 = sum1 + wk(n)
    End Do
    qans = ( r / x ) * ( 1.0D+00 + sum1 )
    ans = 0.5D+00 + ( 0.5D+00 - qans )
    Return
    !
    !  Taylor series for P(A,X)/X**A
    !
160 Continue

    an = 3.0D+00
    c = x
    sum1 = x / ( a + 3.0D+00 )
    tol = 3.0D+00 * acc / ( a + 1.0D+00 )

    Do

       an = an + 1.0D+00
       c = -c * ( x / an )
       t = c / ( a + an )
       sum1 = sum1 + t

       If ( Abs ( t ) <= tol ) Then
          Exit
       End If

    End Do

    j = a * x * ( ( sum1 / 6.0D+00 - 0.5D+00 / &
         ( a +  2.0D+00  ) ) * x + 1.0D+00 &
         / ( a + 1.0D+00 ) )

    z = a * Log ( x )
    h = gam1 ( a )
    g = 1.0D+00 + h

    If ( x < 0.25D+00 ) Then
       go to 180
    End If

    If ( a < x / 2.59D+00 ) Then
       go to 200
    End If

    go to 190

180 Continue

    If ( -0.13394D+00 < z ) Then
       go to 200
    End If

190 Continue

    w = Exp ( z )
    ans = w * g * ( 0.5D+00 + ( 0.5D+00 - j ))
    qans = 0.5D+00 + ( 0.5D+00 - ans )
    Return

200 Continue

    l = rexp ( z )
    w = 0.5D+00 + ( 0.5D+00 + l )
    qans = ( w * j - l ) * g - h

    If ( qans < 0.0D+00 ) Then
       ans = 1.0D+00
       qans = 0.0D+00
       Return
    End If

    ans = 0.5D+00 + ( 0.5D+00 - qans )
    Return
    !
    !  Finite sums for Q when 1 <= A and 2*A is an integer.
    !
210 Continue

    sum1 = Exp ( - x )
    t = sum1
    n = 1
    c = 0.0D+00
    go to 230

220 Continue

    rtx = Sqrt ( x )
    sum1 = error_fc ( 0, rtx )
    t = Exp ( -x ) / ( rtpi * rtx )
    n = 0
    c = -0.5D+00

230 Continue

    Do While ( n /= i )
       n = n + 1
       c = c + 1.0D+00
       t = ( x * t ) / c
       sum1 = sum1 + t
    End Do

! 240 Continue

    qans = sum1
    ans = 0.5D+00 + ( 0.5D+00 - qans )
    Return
    !
    !  Continued fraction expansion.
    !
250 Continue

    tol = Max ( 5.0D+00 * e, acc )
    a2nm1 = 1.0D+00
    a2n = 1.0D+00
    b2nm1 = x
    b2n = x + ( 1.0D+00 - a )
    c = 1.0D+00

    Do

       a2nm1 = x * a2n + c * a2nm1
       b2nm1 = x * b2n + c * b2nm1
       am0 = a2nm1 / b2nm1
       c = c + 1.0D+00
       cma = c - a
       a2n = a2nm1 + cma * a2n
       b2n = b2nm1 + cma * b2n
       an0 = a2n / b2n

       If ( Abs ( an0 - am0 ) < tol * an0 ) Then
          Exit
       End If

    End Do

    qans = r * an0
    ans = 0.5D+00 + ( 0.5D+00 - qans )
    Return
    !
    !  General Temme expansion.
    !
270 Continue

    If ( Abs ( s ) <= 2.0D+00 * e .And. 3.28D-03 < a * e * e ) Then
       ans =  2.0D+00
       Return
    End If

    c = Exp ( - y )
    w = 0.5D+00 * error_fc ( 1, Sqrt ( y ) )
    u = 1.0D+00 / a
    z = Sqrt ( z + z )

    If ( l < 1.0D+00 ) Then
       z = -z
    End If

    If ( iop < 2 ) Then

       If ( Abs ( s ) <= 1.0D-03 ) Then

          c0 = ((((((     &
               d0(7)   &
               * z + d0(6) ) &
               * z + d0(5) ) &
               * z + d0(4) ) &
               * z + d0(3) ) &
               * z + d0(2) ) &
               * z + d0(1) ) &
               * z - 1.0D+00 / 3.0D+00

          c1 = (((((      &
               d1(6)   &
               * z + d1(5) ) &
               * z + d1(4) ) &
               * z + d1(3) ) &
               * z + d1(2) ) &
               * z + d1(1) ) &
               * z + d10

          c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20

          c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30

          c4 = ( d4(2) * z + d4(1) ) * z + d40
          c5 = ( d5(2) * z + d5(1) ) * z + d50
          c6 = d6(1) * z + d60

          t = (((((( d70 &
               * u + c6 ) &
               * u + c5 ) &
               * u + c4 ) &
               * u + c3 ) &
               * u + c2 ) &
               * u + c1 ) &
               * u + c0

       Else

          c0 = (((((((((((( &
               d0(13)   &
               * z + d0(12) ) &
               * z + d0(11) ) &
               * z + d0(10) ) &
               * z + d0(9)  ) &
               * z + d0(8)  ) &
               * z + d0(7)  ) &
               * z + d0(6)  ) &
               * z + d0(5)  ) &
               * z + d0(4)  ) &
               * z + d0(3)  ) &
               * z + d0(2)  ) &
               * z + d0(1)  ) &
               * z - 1.0D+00 / 3.0D+00

          c1 = ((((((((((( &
               d1(12) &
               * z + d1(11) &
               ) * z + d1(10) &
               ) * z + d1(9)  &
               ) * z + d1(8)  &
               ) * z + d1(7)  &
               ) * z + d1(6)  &
               ) * z + d1(5)  &
               ) * z + d1(4)  &
               ) * z + d1(3)  &
               ) * z + d1(2)  &
               ) * z + d1(1)  &
               ) * z + d10

          c2 = ((((((((( &
               d2(10) &
               * z + d2(9) &
               ) * z + d2(8) &
               ) * z + d2(7) &
               ) * z + d2(6) &
               ) * z + d2(5) &
               ) * z + d2(4) &
               ) * z + d2(3) &
               ) * z + d2(2) &
               ) * z + d2(1) &
               ) * z + d20

          c3 = ((((((( &
               d3(8) &
               * z + d3(7) &
               ) * z + d3(6) &
               ) * z + d3(5) &
               ) * z + d3(4) &
               ) * z + d3(3) &
               ) * z + d3(2) &
               ) * z + d3(1) &
               ) * z + d30

          c4 = ((((( d4(6)*z+d4(5))*z+d4(4))*z+d4(3))*z+d4(2))*z+d4(1))*z + d40

          c5 = (((d5(4)*z+d5(3))*z+d5(2))*z+d5(1))*z + d50

          c6 = ( d6(2) * z + d6(1) ) * z + d60

          t = ((((((   &
               d70    &
               * u + c6 ) &
               * u + c5 ) &
               * u + c4 ) &
               * u + c3 ) &
               * u + c2 ) &
               * u + c1 ) &
               * u + c0

       End If

    Else If ( iop == 2 ) Then

       c0 = (((((      &
            d0(6)   &
            * z + d0(5) ) &
            * z + d0(4) ) &
            * z + d0(3) ) &
            * z + d0(2) ) &
            * z + d0(1) ) &
            * z - 1.0D+00 / 3.0D+00

       c1 = ((( d1(4) * z + d1(3) ) * z + d1(2) ) * z + d1(1) ) * z + d10
       c2 = d2(1) * z + d20
       t = ( c2 * u + c1 ) * u + c0

    Else If ( 2 < iop ) Then

       t = (( d0(3) * z + d0(2) ) * z + d0(1) ) * z - 1.0D+00 / 3.0D+00

    End If

310 Continue

    If ( 1.0D+00 <= l ) Then
       qans = c * ( w + rt2pin * t / rta )
       ans = 0.5D+00 + ( 0.5D+00 - qans )
    Else
       ans = c * ( w - rt2pin * t / rta )
       qans = 0.5D+00 + ( 0.5D+00 - ans )
    End If

    Return
    !
    !  Temme expansion for L = 1
    !
330 Continue

    If ( 3.28D-03 < a * e * e ) Then
       ans =  2.0D+00
       Return
    End If

    c = 0.5D+00 + ( 0.5D+00 - y )
    w = ( 0.5D+00 - Sqrt ( y ) &
         * ( 0.5D+00 &
         + ( 0.5D+00 - y / 3.0D+00 ) ) / rtpi ) / c
    u = 1.0D+00 / a
    z = Sqrt ( z + z )

    If ( l < 1.0D+00 ) Then
       z = -z
    End If

    If ( iop < 2 ) Then

       c0 = ((((((     &
            d0(7)   &
            * z + d0(6) ) &
            * z + d0(5) ) &
            * z + d0(4) ) &
            * z + d0(3) ) &
            * z + d0(2) ) &
            * z + d0(1) ) &
            * z - 1.0D+00 / 3.0D+00

       c1 = (((((      &
            d1(6)   &
            * z + d1(5) ) &
            * z + d1(4) ) &
            * z + d1(3) ) &
            * z + d1(2) ) &
            * z + d1(1) ) &
            * z + d10

       c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20

       c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30

       c4 = ( d4(2) * z + d4(1) ) * z + d40
       c5 = ( d5(2) * z + d5(1) ) * z + d50
       c6 = d6(1) * z + d60

       t = (((((( d70 &
            * u + c6 ) &
            * u + c5 ) &
            * u + c4 ) &
            * u + c3 ) &
            * u + c2 ) &
            * u + c1 ) &
            * u + c0

    Else If ( iop == 2 ) Then

       c0 = ( d0(2) * z + d0(1) ) * z - 1.0D+00 / 3.0D+00
       c1 = d1(1) * z + d10
       t = ( d20 * u + c1 ) * u + c0

    Else If ( 2 < iop ) Then

       t = d0(1) * z - 1.0D+00 / 3.0D+00

    End If

    go to 310
    !
    !  Special cases
    !
390 Continue

    If ( x < 0.25D+00 ) Then
       ans = error_f ( Sqrt ( x ) )
       qans = 0.5D+00 + ( 0.5D+00 - ans )
    Else
       qans = error_fc ( 0, Sqrt ( x ) )
       ans = 0.5D+00 + ( 0.5D+00 - qans )
    End If

    Return

410 Continue

    If ( Abs ( s ) <= 2.0D+00 * e ) Then
       ans =  2.0D+00
       Return
    End If

    If ( x <= a ) Then
       ans = 0.0D+00
       qans = 1.0D+00
    Else
       ans = 1.0D+00
       qans = 0.0D+00
    End If

    Return
  End Subroutine gamma_inc

  Subroutine gamma_inc_inv ( a, x, x0, p, q, ierr )

    !*****************************************************************************80
    !
    !! GAMMA_INC_INV computes the inverse incomplete gamma ratio function.
    !
    !  Discussion:
    !
    !    The routine is given positive A, and nonnegative P and Q where P + Q = 1.
    !    The value X is computed with the property that P(A,X) = P and Q(A,X) = Q.
    !    Schroder iteration is employed.  The routine attempts to compute X
    !    to 10 significant digits if this is possible for the particular computer
    !    arithmetic being used.
    !
    !  Author:
    !
    !    Alfred Morris,
    !    Naval Surface Weapons Center,
    !    Dahlgren, Virginia.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, the parameter in the incomplete gamma
    !    ratio.  A must be positive.
    !
    !    Output, real ( kind = 8 ) X, the computed point for which the
    !    incomplete gamma functions have the values P and Q.
    !
    !    Input, real ( kind = 8 ) X0, an optional initial approximation
    !    for the solution X.  If the user does not want to supply an
    !    initial approximation, then X0 should be set to 0, or a negative
    !    value.
    !
    !    Input, real ( kind = 8 ) P, Q, the values of the incomplete gamma
    !    functions, for which the corresponding argument is desired.
    !
    !    Output, integer IERR, error flag.
    !    0, the solution was obtained. Iteration was not used.
    !    0 < K, The solution was obtained. IERR iterations were performed.
    !    -2, A <= 0
    !    -3, No solution was obtained. The ratio Q/A is too large.
    !    -4, P + Q /= 1
    !    -6, 20 iterations were performed. The most recent value obtained
    !        for X is given.  This cannot occur if X0 <= 0.
    !    -7, Iteration failed. No value is given for X.
    !        This may occur when X is approximately 0.
    !    -8, A value for X has been obtained, but the routine is not certain
    !        of its accuracy.  Iteration cannot be performed in this
    !        case. If X0 <= 0, this can occur only when P or Q is
    !        approximately 0. If X0 is positive then this can occur when A is
    !        exceedingly close to X and A is extremely large (say A >= 1.E20).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ), Parameter :: a0 = 3.31125922108741D+00
    Real ( kind = 8 ), Parameter :: a1 = 11.6616720288968D+00
    Real ( kind = 8 ), Parameter :: a2 = 4.28342155967104D+00
    Real ( kind = 8 ), Parameter :: a3 = 0.213623493715853D+00
    !    Real ( kind = 8 ) alnrel
    Real ( kind = 8 ) am1
    Real ( kind = 8 ) amax
    Real ( kind = 8 ), Dimension(2) :: amin = (/ &
         500.0D+00, 100.0D+00 /)
    Real ( kind = 8 ) ap1
    Real ( kind = 8 ) ap2
    Real ( kind = 8 ) ap3
    Real ( kind = 8 ) apn
    Real ( kind = 8 ) b
    Real ( kind = 8 ), Parameter :: b1 = 6.61053765625462D+00
    Real ( kind = 8 ), Parameter :: b2 = 6.40691597760039D+00
    Real ( kind = 8 ), Parameter :: b3 = 1.27364489782223D+00
    Real ( kind = 8 ), Parameter :: b4 = .036117081018842D+00
    Real ( kind = 8 ), Dimension ( 2 ) :: bmin = (/ &
         1.0D-28, 1.0D-13 /)
    Real ( kind = 8 ), Parameter :: c = 0.577215664901533D+00
    Real ( kind = 8 ) c1
    Real ( kind = 8 ) c2
    Real ( kind = 8 ) c3
    Real ( kind = 8 ) c4
    Real ( kind = 8 ) c5
    Real ( kind = 8 ) d
    Real ( kind = 8 ), Dimension ( 2 ) :: dmin = (/ &
         1.0D-06, 1.0D-04 /)
    Real ( kind = 8 ) e
    Real ( kind = 8 ) e2
    Real ( kind = 8 ), Dimension ( 2 ) :: emin = (/ &
         2.0D-03, 6.0D-03 /)
    Real ( kind = 8 ) eps
    Real ( kind = 8 ), Dimension ( 2 ) :: eps0 = (/ &
         1.0D-10, 1.0D-08 /)
    Real ( kind = 8 ) g
    !    Real ( kind = 8 ) gamma_log
    !    Real ( kind = 8 ) gamma_ln1
    !    Real ( kind = 8 ) gamma
    Real ( kind = 8 ) h
    Real ( kind = 8 ), Parameter :: half = 0.5D+00
    Integer ierr
    Integer iop
    Real ( kind = 8 ), Parameter :: ln10 = 2.302585D+00
    Real ( kind = 8 ) p
    Real ( kind = 8 ) pn
    Real ( kind = 8 ) q
    Real ( kind = 8 ) qg
    Real ( kind = 8 ) qn
    Real ( kind = 8 ) r
    !    Real ( kind = 8 ) rcomp
    Real ( kind = 8 ) rta
    Real ( kind = 8 ) s
    Real ( kind = 8 ) s2
    Real ( kind = 8 ) sum1
    Real ( kind = 8 ) t
    Real ( kind = 8 ), Parameter :: tol = 1.0D-05
    Real ( kind = 8 ), Parameter :: two =  2.0D+00
    Real ( kind = 8 ) u
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    Real ( kind = 8 ) x0
    Real ( kind = 8 ) xn
    Real ( kind = 8 ) y
    Real ( kind = 8 ) z

    e = Epsilon ( e )

    x = 0.0D+00

    If ( a <= 0.0D+00 ) Then
       ierr = -2
       Return
    End If

    t = p + q - 1.0D+00

    If ( e < Abs ( t ) ) Then
       ierr = -4
       Return
    End If

    ierr = 0

    If ( p == 0.0D+00 ) Then
       Return
    End If

    If ( q == 0.0D+00 ) Then
       x = Huge ( x )
       Return
    End If

    If ( a == 1.0D+00 ) Then
       If ( 0.9D+00 <= q ) Then
          x = -alnrel ( - p )
       Else
          x = -Log ( q )
       End If
       Return
    End If

    e2 = two * e
    amax = 0.4D-10 / ( e * e )

    If ( 1.0D-10 < e ) Then
       iop = 2
    Else
       iop = 1
    End If

    eps = eps0(iop)
    xn = x0

    If ( 0.0D+00 < x0 ) Then
       go to 160
    End If
    !
    !  Selection of the initial approximation XN of X when A < 1.
    !
    If ( 1.0D+00 < a ) Then
       go to 80
    End If

    g = gamma ( a + 1.0D+00 )
    qg = q * g

    If ( qg == 0.0D+00 ) Then
       x = Huge ( x )
       ierr = -8
       Return
    End If

    b = qg / a

    If ( 0.6D+00 * a < qg ) Then
       go to 40
    End If

    If ( a < 0.30D+00 .And. 0.35D+00 <= b ) Then
       t = Exp ( - ( b + c ) )
       u = t * Exp ( t )
       xn = t * Exp ( u )
       go to 160
    End If

    If ( 0.45D+00 <= b ) Then
       go to 40
    End If

    If ( b == 0.0D+00 ) Then
       x = Huge ( x )
       ierr = -8
       Return
    End If

    y = -Log ( b )
    s = half + ( half - a )
    z = Log ( y )
    t = y - s * z

    If ( 0.15D+00 <= b ) Then
       xn = y - s * Log ( t ) - Log ( 1.0D+00 + s / ( t + 1.0D+00 ) )
       go to 220
    End If

    If ( 0.01D+00 < b ) Then
       u = ( ( t + two * ( 3.0D+00 - a ) ) * t &
            + ( two - a ) * ( 3.0D+00 - a )) / &
            ( ( t + ( 5.0D+00 - a ) ) * t + two )
       xn = y - s * Log ( t ) - Log ( u )
       go to 220
    End If

30  Continue

    c1 = -s * z
    c2 = -s * ( 1.0D+00 + c1 )

    c3 = s * (( half * c1 &
         + ( two - a ) ) * c1 + ( 2.5D+00 - 1.5D+00 * a ) )

    c4 = -s * ((( c1 / 3.0D+00 + ( 2.5D+00 - 1.5D+00 * a ) ) * c1 &
         + ( ( a - 6.0D+00 ) * a + 7.0D+00 ) ) &
         * c1 + ( ( 11.0D+00 * a - 46.0D+00 ) * a + 47.0D+00 ) / 6.0D+00 )

    c5 = -s * (((( - c1 / 4.0D+00 + ( 11.0D+00 * a - 17.0D+00 ) / 6.0D+00 ) * c1 &
         + ( ( -3.0D+00 * a + 13.0D+00 ) * a - 13.0D+00 ) ) * c1 &
         + half &
         * ( ( ( two * a - 25.0D+00 ) * a + 72.0D+00 ) &
         * a - 61.0D+00 ) ) * c1 &
         + ( ( ( 25.0D+00 * a - 195.0D+00 ) * a &
         + 477.0D+00 ) * a - 379.0D+00 ) / 12.0D+00 )

    xn = (((( c5 / y + c4 ) / y + c3 ) / y + c2 ) / y + c1 ) + y

    If ( 1.0D+00 < a ) Then
       go to 220
    End If

    If ( bmin(iop) < b ) Then
       go to 220
    End If

    x = xn
    Return

40  Continue

    If ( b * q <= 1.0D-08 ) Then
       xn = Exp ( - ( q / a + c ))
    Else If ( 0.9D+00 < p ) Then
       xn = Exp ( ( alnrel ( - q ) + gamma_ln1 ( a )) / a )
    Else
       xn = Exp ( Log ( p * g ) / a )
    End If

    If ( xn == 0.0D+00 ) Then
       ierr = -3
       Return
    End If

    t = half + ( half - xn / ( a + 1.0D+00 ))
    xn = xn / t
    go to 160
    !
    !  Selection of the initial approximation XN of X when 1 < A.
    !
80  Continue

    If ( 0.5D+00 < q ) Then
       w = Log ( p )
    Else
       w = Log ( q )
    End If

    t = Sqrt ( - two * w )

    s = t - ((( a3 * t + a2 ) * t + a1 ) * t + a0 ) / (((( &
         b4 * t + b3 ) * t + b2 ) * t + b1 ) * t + 1.0D+00 )

    If ( 0.5D+00 < q ) Then
       s = -s
    End If

    rta = Sqrt ( a )
    s2 = s * s

    xn = a + s * rta + ( s2 - 1.0D+00 ) / 3.0D+00 + s * ( s2 - 7.0D+00 ) &
         / ( 36.0D+00 * rta ) - ( ( 3.0D+00 * s2 + 7.0D+00 ) * s2 - 16.0D+00 ) &
         / ( 810.0D+00 * a ) + s * (( 9.0D+00 * s2 + 256.0D+00 ) * s2 - 433.0D+00 ) &
         / ( 38880.0D+00 * a * rta )

    xn = Max ( xn, 0.0D+00 )

    If ( amin(iop) <= a ) Then

       x = xn
       d = half + ( half - x / a )

       If ( Abs ( d ) <= dmin(iop) ) Then
          Return
       End If

    End If

! 110 Continue

    If ( p <= 0.5D+00 ) Then
       go to 130
    End If

    If ( xn < 3.0D+00 * a ) Then
       go to 220
    End If

    y = - ( w + gamma_log ( a ) )
    d = Max ( two, a * ( a - 1.0D+00 ) )

    If ( ln10 * d <= y ) Then
       s = 1.0D+00 - a
       z = Log ( y )
       go to 30
    End If

! 120 Continue

    t = a - 1.0D+00
    xn = y + t * Log ( xn ) - alnrel ( -t / ( xn + 1.0D+00 ) )
    xn = y + t * Log ( xn ) - alnrel ( -t / ( xn + 1.0D+00 ) )
    go to 220

130 Continue

    ap1 = a + 1.0D+00

    If ( 0.70D+00 * ap1 < xn ) Then
       go to 170
    End If

    w = w + gamma_log ( ap1 )

    If ( xn <= 0.15 * ap1 ) Then
       ap2 = a + two
       ap3 = a + 3.0D+00
       x = Exp ( ( w + x ) / a )
       x = Exp ( ( w + x - Log ( 1.0D+00 + ( x / ap1 ) &
            * ( 1.0D+00 + x / ap2 ) ) ) / a )
       x = Exp ( ( w + x - Log ( 1.0D+00 + ( x / ap1 ) &
            * ( 1.0D+00 + x / ap2 ) ) ) / a )
       x = Exp ( ( w + x - Log ( 1.0D+00 + ( x / ap1 ) &
            * ( 1.0D+00 + ( x / ap2 ) &
            * ( 1.0D+00 + x / ap3 ) ) ) ) / a )
       xn = x

       If ( xn <= 1.0D-02 * ap1 ) Then
          If ( xn <= emin(iop) * ap1 ) Then
             Return
          End If
          go to 170
       End If

    End If

    apn = ap1
    t = xn / apn
    sum1 = 1.0D+00 + t

    Do

       apn = apn + 1.0D+00
       t = t * ( xn / apn )
       sum1 = sum1 + t

       If ( t <= 1.0D-04 ) Then
          Exit
       End If

    End Do

    t = w - Log ( sum1 )
    xn = Exp ( ( xn + t ) / a )
    xn = xn * ( 1.0D+00 - ( a * Log ( xn ) - xn - t ) / ( a - xn ) )
    go to 170
    !
    !  Schroder iteration using P.
    !
160 Continue

    If ( 0.5D+00 < p ) Then
       go to 220
    End If

170 Continue

    If ( p <= 1.0D+10 * Tiny ( p ) ) Then
       x = xn
       ierr = -8
       Return
    End If

    am1 = ( a - half ) - half

180 Continue

    If ( amax < a ) Then
       d = half + ( half - xn / a )
       If ( Abs ( d ) <= e2 ) Then
          x = xn
          ierr = -8
          Return
       End If
    End If

! 190 Continue

    If ( 20 <= ierr ) Then
       ierr = -6
       Return
    End If

    ierr = ierr + 1
    Call gamma_inc ( a, xn, pn, qn, 0 )

    If ( pn == 0.0D+00 .Or. qn == 0.0D+00 ) Then
       x = xn
       ierr = -8
       Return
    End If

    r = rcomp ( a, xn )

    If ( r == 0.0D+00 ) Then
       x = xn
       ierr = -8
       Return
    End If

    t = ( pn - p ) / r
    w = half * ( am1 - xn )

    If ( Abs ( t ) <= 0.1D+00 .And. Abs ( w * t ) <= 0.1D+00 ) Then
       go to 200
    End If

    x = xn * ( 1.0D+00 - t )

    If ( x <= 0.0D+00 ) Then
       ierr = -7
       Return
    End If

    d = Abs ( t )
    go to 210

200 Continue

    h = t * ( 1.0D+00 + w * t )
    x = xn * ( 1.0D+00 - h )

    If ( x <= 0.0D+00 ) Then
       ierr = -7
       Return
    End If

    If ( 1.0D+00 <= Abs ( w ) .And. Abs ( w ) * t * t <= eps ) Then
       Return
    End If

    d = Abs ( h )

210 Continue

    xn = x

    If ( d <= tol ) Then

       If ( d <= eps ) Then
          Return
       End If

       If ( Abs ( p - pn ) <= tol * p ) Then
          Return
       End If

    End If

    go to 180
    !
    !  Schroder iteration using Q.
    !
220 Continue

    If ( q <= 1.0D+10 * Tiny ( q ) ) Then
       x = xn
       ierr = -8
       Return
    End If

    am1 = ( a - half ) - half

230 Continue

    If ( amax < a ) Then
       d = half + ( half - xn / a )
       If ( Abs ( d ) <= e2 ) Then
          x = xn
          ierr = -8
          Return
       End If
    End If

    If ( 20 <= ierr ) Then
       ierr = -6
       Return
    End If

    ierr = ierr + 1
    Call gamma_inc ( a, xn, pn, qn, 0 )

    If ( pn == 0.0D+00 .Or. qn == 0.0D+00 ) Then
       x = xn
       ierr = -8
       Return
    End If

    r = rcomp ( a, xn )

    If ( r == 0.0D+00 ) Then
       x = xn
       ierr = -8
       Return
    End If

    t = ( q - qn ) / r
    w = half * ( am1 - xn )

    If ( Abs ( t ) <= 0.1 .And. Abs ( w * t ) <= 0.1 ) Then
       go to 250
    End If

    x = xn * ( 1.0D+00 - t )

    If ( x <= 0.0D+00 ) Then
       ierr = -7
       Return
    End If

    d = Abs ( t )
    go to 260

250 Continue

    h = t * ( 1.0D+00 + w * t )
    x = xn * ( 1.0D+00 - h )

    If ( x <= 0.0D+00 ) Then
       ierr = -7
       Return
    End If

    If ( 1.0D+00 <= Abs ( w ) .And. Abs ( w ) * t * t <= eps ) Then
       Return
    End If

    d = Abs ( h )

260 Continue

    xn = x

    If ( tol < d ) Then
       go to 230
    End If

    If ( d <= eps ) Then
       Return
    End If

    If ( Abs ( q - qn ) <= tol * q ) Then
       Return
    End If

    go to 230
  End Subroutine gamma_inc_inv

!!$  Subroutine gamma_inc_values ( n_data, a, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    The (normalized) incomplete Gamma function P(A,X) is defined as:
!!$    !
!!$    !      PN(A,X) = 1/GAMMA(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
!!$    !
!!$    !    With this definition, for all A and X,
!!$    !
!!$    !      0 <= PN(A,X) <= 1
!!$    !
!!$    !    and
!!$    !
!!$    !      PN(A,INFINITY) = 1.0
!!$    !
!!$    !    Mathematica can compute this value as
!!$    !
!!$    !      1 - GammaRegularized[A,X]
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    08 May 2001
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, real ( kind = 8 ) A, X, the arguments of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 20
!!$
!!$    Real ( kind = 8 ) a
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: a_vec = (/ &
!!$         0.1D+00,  0.1D+00,  0.1D+00,  0.5D+00, &
!!$         0.5D+00,  0.5D+00,  1.0D+00,  1.0D+00, &
!!$         1.0D+00,  1.1D+00,  1.1D+00,  1.1D+00, &
!!$         2.0D+00,  2.0D+00,  2.0D+00,  6.0D+00, &
!!$         6.0D+00, 11.0D+00, 26.0D+00, 41.0D+00 /)
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.7420263D+00, 0.9119753D+00, 0.9898955D+00, 0.2931279D+00, &
!!$         0.7656418D+00, 0.9921661D+00, 0.0951626D+00, 0.6321206D+00, &
!!$         0.9932621D+00, 0.0757471D+00, 0.6076457D+00, 0.9933425D+00, &
!!$         0.0091054D+00, 0.4130643D+00, 0.9931450D+00, 0.0387318D+00, &
!!$         0.9825937D+00, 0.9404267D+00, 0.4863866D+00, 0.7359709D+00 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         3.1622777D-02, 3.1622777D-01, 1.5811388D+00, 7.0710678D-02, &
!!$         7.0710678D-01, 3.5355339D+00, 0.1000000D+00, 1.0000000D+00, &
!!$         5.0000000D+00, 1.0488088D-01, 1.0488088D+00, 5.2440442D+00, &
!!$         1.4142136D-01, 1.4142136D+00, 7.0710678D+00, 2.4494897D+00, &
!!$         1.2247449D+01, 1.6583124D+01, 2.5495098D+01, 4.4821870D+01 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       a = 0.0D+00
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       a = a_vec(n_data)
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine gamma_inc_values

  Function gamma_ln1 ( a )

    !*****************************************************************************80
    !
    !! GAMMA_LN1 evaluates ln ( Gamma ( 1 + A ) ), for -0.2 <= A <= 1.25.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, defines the argument of the function.
    !
    !    Output, real ( kind = 8 ) GAMMA_LN1, the value of ln ( Gamma ( 1 + A ) ).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) bot
    Real ( kind = 8 ) gamma_ln1
    Real ( kind = 8 ), Parameter :: p0 =  0.577215664901533D+00
    Real ( kind = 8 ), Parameter :: p1 =  0.844203922187225D+00
    Real ( kind = 8 ), Parameter :: p2 = -0.168860593646662D+00
    Real ( kind = 8 ), Parameter :: p3 = -0.780427615533591D+00
    Real ( kind = 8 ), Parameter :: p4 = -0.402055799310489D+00
    Real ( kind = 8 ), Parameter :: p5 = -0.673562214325671D-01
    Real ( kind = 8 ), Parameter :: p6 = -0.271935708322958D-02
    Real ( kind = 8 ), Parameter :: q1 =  0.288743195473681D+01
    Real ( kind = 8 ), Parameter :: q2 =  0.312755088914843D+01
    Real ( kind = 8 ), Parameter :: q3 =  0.156875193295039D+01
    Real ( kind = 8 ), Parameter :: q4 =  0.361951990101499D+00
    Real ( kind = 8 ), Parameter :: q5 =  0.325038868253937D-01
    Real ( kind = 8 ), Parameter :: q6 =  0.667465618796164D-03
    Real ( kind = 8 ), Parameter :: r0 = 0.422784335098467D+00
    Real ( kind = 8 ), Parameter :: r1 = 0.848044614534529D+00
    Real ( kind = 8 ), Parameter :: r2 = 0.565221050691933D+00
    Real ( kind = 8 ), Parameter :: r3 = 0.156513060486551D+00
    Real ( kind = 8 ), Parameter :: r4 = 0.170502484022650D-01
    Real ( kind = 8 ), Parameter :: r5 = 0.497958207639485D-03
    Real ( kind = 8 ), Parameter :: s1 = 0.124313399877507D+01
    Real ( kind = 8 ), Parameter :: s2 = 0.548042109832463D+00
    Real ( kind = 8 ), Parameter :: s3 = 0.101552187439830D+00
    Real ( kind = 8 ), Parameter :: s4 = 0.713309612391000D-02
    Real ( kind = 8 ), Parameter :: s5 = 0.116165475989616D-03
    Real ( kind = 8 ) top
    Real ( kind = 8 ) x

    If ( a < 0.6D+00 ) Then

       top = (((((  &
            p6   &
            * a + p5 ) &
            * a + p4 ) &
            * a + p3 ) &
            * a + p2 ) &
            * a + p1 ) &
            * a + p0

       bot = (((((  &
            q6   &
            * a + q5 ) &
            * a + q4 ) &
            * a + q3 ) &
            * a + q2 ) &
            * a + q1 ) &
            * a + 1.0D+00

       gamma_ln1 = -a * ( top / bot )

    Else

       x = ( a - 0.5D+00 ) - 0.5D+00

       top = ((((( r5 * x + r4 ) * x + r3 ) * x + r2 ) * x + r1 ) * x + r0 )

       bot = ((((( s5 * x + s4 ) * x + s3 ) * x + s2 ) * x + s1 ) * x + 1.0D+00 )

       gamma_ln1 = x * ( top / bot )

    End If

    Return
  End Function gamma_ln1

  Function gamma_log ( a )

    !*****************************************************************************80
    !
    !! GAMMA_LOG evaluates ln ( Gamma ( A ) ) for positive A.
    !
    !  Author:
    !
    !    Alfred Morris,
    !    Naval Surface Weapons Center,
    !    Dahlgren, Virginia.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, the argument of the function.
    !    A should be positive.
    !
    !    Output, real ( kind = 8 ), GAMMA_LOG, the value of ln ( Gamma ( A ) ).
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ), Parameter :: c0 =  0.833333333333333D-01
    Real ( kind = 8 ), Parameter :: c1 = -0.277777777760991D-02
    Real ( kind = 8 ), Parameter :: c2 =  0.793650666825390D-03
    Real ( kind = 8 ), Parameter :: c3 = -0.595202931351870D-03
    Real ( kind = 8 ), Parameter :: c4 =  0.837308034031215D-03
    Real ( kind = 8 ), Parameter :: c5 = -0.165322962780713D-02
    Real ( kind = 8 ), Parameter :: d  =  0.418938533204673D+00
    Real ( kind = 8 ) gamma_log
    !    Real ( kind = 8 ) gamma_ln1
    Integer i
    Integer n
    Real ( kind = 8 ) t
    Real ( kind = 8 ) w

    If ( a <= 0.8D+00 ) Then

       gamma_log = gamma_ln1 ( a ) - Log ( a )

    Else If ( a <= 2.25D+00 ) Then

       t = ( a - 0.5D+00 ) - 0.5D+00
       gamma_log = gamma_ln1 ( t )

    Else If ( a < 10.0D+00 ) Then

       n = INT(a - 1.25D+00)
       t = a
       w = 1.0D+00
       Do i = 1, n
          t = t - 1.0D+00
          w = t * w
       End Do

       gamma_log = gamma_ln1 ( t - 1.0D+00 ) + Log ( w )

    Else

       t = ( 1.0D+00 / a )**2

       w = ((((( c5 * t + c4 ) * t + c3 ) * t + c2 ) * t + c1 ) * t + c0 ) / a

       gamma_log = ( d + w ) + ( a - 0.5D+00 ) &
            * ( Log ( a ) - 1.0D+00 )

    End If

    Return
  End Function gamma_log

  Subroutine gamma_rat1 ( a, x, r, p, q, eps )

    !*****************************************************************************80
    !
    !! GAMMA_RAT1 evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, X, the parameters of the functions.
    !    It is assumed that A <= 1.
    !
    !    Input, real ( kind = 8 ) R, the value exp(-X) * X**A / Gamma(A).
    !
    !    Output, real ( kind = 8 ) P, Q, the values of P(A,X) and Q(A,X).
    !
    !    Input, real ( kind = 8 ) EPS, the tolerance.
    !
    Implicit None

    Real ( kind = 8 ) a
    Real ( kind = 8 ) a2n
    Real ( kind = 8 ) a2nm1
    Real ( kind = 8 ) am0
    Real ( kind = 8 ) an
    Real ( kind = 8 ) an0
    Real ( kind = 8 ) b2n
    Real ( kind = 8 ) b2nm1
    Real ( kind = 8 ) c
    Real ( kind = 8 ) cma
    Real ( kind = 8 ) eps
    !    Real ( kind = 8 ) error_f
    !    Real ( kind = 8 ) error_fc
    Real ( kind = 8 ) g
    !    Real ( kind = 8 ) gam1
    Real ( kind = 8 ) h
    Real ( kind = 8 ) j
    Real ( kind = 8 ) l
    Real ( kind = 8 ) p
    Real ( kind = 8 ) q
    Real ( kind = 8 ) r
    !    Real ( kind = 8 ) rexp
    Real ( kind = 8 ) sum1
    Real ( kind = 8 ) t
    Real ( kind = 8 ) tol
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    Real ( kind = 8 ) z

    If ( a * x == 0.0D+00 ) Then

       If ( x <= a ) Then
          p = 0.0D+00
          q = 1.0D+00
       Else
          p = 1.0D+00
          q = 0.0D+00
       End If

       Return
    End If

    If ( a == 0.5D+00 ) Then

       If ( x < 0.25D+00 ) Then
          p = error_f ( Sqrt ( x ) )
          q = 0.5D+00 + ( 0.5D+00 - p )
       Else
          q = error_fc ( 0, Sqrt ( x ) )
          p = 0.5D+00 + ( 0.5D+00 - q )
       End If

       Return

    End If
    !
    !  Taylor series for P(A,X)/X**A
    !
    If ( x < 1.1D+00 ) Then

       an = 3.0
       c = x
       sum1 = x / ( a + 3.0D+00 )
       tol = 0.1D+00 * eps / ( a + 1.0D+00 )

       Do

          an = an + 1.0D+00
          c = -c * ( x / an )
          t = c / ( a + an )
          sum1 = sum1 + t

          If ( Abs ( t ) <= tol ) Then
             Exit
          End If

       End Do

       j = a * x * ( ( sum1 / 6.0D+00 - 0.5D+00 &
            / ( a +  2.0D+00  ) ) &
            * x + 1.0D+00 / ( a + 1.0D+00 ) )

       z = a * Log ( x )
       h = gam1 ( a )
       g = 1.0D+00 + h

       If ( x < 0.25D+00 ) Then
          go to 30
       End If

       If ( a < x / 2.59D+00 ) Then
          go to 50
       Else
          go to 40
       End If

30     Continue

       If ( -0.13394D+00 < z ) Then
          go to 50
       End If

40     Continue

       w = Exp ( z )
       p = w * g * ( 0.5D+00 + ( 0.5D+00 - j ))
       q = 0.5D+00 + ( 0.5D+00 - p )
       Return

50     Continue

       l = rexp ( z )
       w = 0.5D+00 + ( 0.5D+00 + l )
       q = ( w * j - l ) * g - h

       If  ( q < 0.0D+00 ) Then
          p = 1.0D+00
          q = 0.0D+00
       Else
          p = 0.5D+00 + ( 0.5D+00 - q )
       End If
       !
       !  Continued fraction expansion.
       !
    Else

       a2nm1 = 1.0D+00
       a2n = 1.0D+00
       b2nm1 = x
       b2n = x + ( 1.0D+00 - a )
       c = 1.0D+00

       Do

          a2nm1 = x * a2n + c * a2nm1
          b2nm1 = x * b2n + c * b2nm1
          am0 = a2nm1 / b2nm1
          c = c + 1.0D+00
          cma = c - a
          a2n = a2nm1 + cma * a2n
          b2n = b2nm1 + cma * b2n
          an0 = a2n / b2n

          If ( Abs ( an0 - am0 ) < eps * an0 ) Then
             Exit
          End If

       End Do

       q = r * an0
       p = 0.5D+00 + ( 0.5D+00 - q )

    End If

    Return
  End Subroutine gamma_rat1

!!$  Subroutine gamma_values ( n_data, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! GAMMA_VALUES returns some values of the Gamma function.
!!$    !
!!$    !  Definition:
!!$    !
!!$    !    Gamma(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) exp(-T) dT
!!$    !
!!$    !  Recursion:
!!$    !
!!$    !    Gamma(X+1) = X * Gamma(X)
!!$    !
!!$    !  Restrictions:
!!$    !
!!$    !    0 < X ( a software restriction).
!!$    !
!!$    !  Special values:
!!$    !
!!$    !    GAMMA(0.5) = sqrt(PI)
!!$    !
!!$    !    For N a positive integer, GAMMA(N+1) = N!, the standard factorial.
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    17 April 2001
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, real ( kind = 8 ) X, the argument of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 18
!!$
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         4.590845D+00,     2.218160D+00,     1.489192D+00,     1.164230D+00, &
!!$         1.0000000000D+00, 0.9513507699D+00, 0.9181687424D+00, 0.8974706963D+00, &
!!$         0.8872638175D+00, 0.8862269255D+00, 0.8935153493D+00, 0.9086387329D+00, &
!!$         0.9313837710D+00, 0.9617658319D+00, 1.0000000000D+00, 3.6288000D+05, &
!!$         1.2164510D+17,    8.8417620D+30 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         0.2D+00,  0.4D+00,  0.6D+00,  0.8D+00, &
!!$         1.0D+00,  1.1D+00,  1.2D+00,  1.3D+00, &
!!$         1.4D+00,  1.5D+00,  1.6D+00,  1.7D+00, &
!!$         1.8D+00,  1.9D+00,  2.0D+00, 10.0D+00, &
!!$         20.0D+00, 30.0D+00 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine gamma_values

  Function gsumln ( a, b )

    !*****************************************************************************80
    !
    !! GSUMLN evaluates the function ln(Gamma(A + B)).
    !
    !  Discussion:
    !
    !    GSUMLN is used for 1 <= A <= 2 and 1 <= B <= 2
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, values whose sum is the argument of
    !    the Gamma function.
    !
    !    Output, real ( kind = 8 ) GSUMLN, the value of ln(Gamma(A+B)).
    !
    Implicit None

    Real ( kind = 8 ) a
    !    Real ( kind = 8 ) alnrel
    Real ( kind = 8 ) b
    !    Real ( kind = 8 ) gamma_ln1
    Real ( kind = 8 ) gsumln
    Real ( kind = 8 ) x

    x = a + b - 2.0D+00

    If ( x <= 0.25D+00 ) Then
       gsumln = gamma_ln1 ( 1.0D+00 + x )
    Else If ( x <= 1.25D+00 ) Then
       gsumln = gamma_ln1 ( x ) + alnrel ( x )
    Else
       gsumln = gamma_ln1 ( x - 1.0D+00 ) + Log ( x * ( 1.0D+00 + x ) )
    End If

    Return
  End Function gsumln

  Function ipmpar ( i )

    !*****************************************************************************80
    !
    !! IPMPAR returns integer machine constants.
    !
    !  Discussion:
    !
    !    Input arguments 1 through 3 are queries about integer arithmetic.
    !    We assume integers are represented in the N-digit, base A form
    !
    !      sign * ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
    !
    !    where 0 <= X(0:N-1) < A.
    !
    !    Then:
    !
    !      IPMPAR(1) = A, the base of integer arithmetic;
    !      IPMPAR(2) = N, the number of base A digits;
    !      IPMPAR(3) = A**N - 1, the largest magnitude.
    !
    !    It is assumed that the single and real ( kind = 8 ) floating
    !    point arithmetics have the same base, say B, and that the
    !    nonzero numbers are represented in the form
    !
    !      sign * (B**E) * (X(1)/B + ... + X(M)/B**M)
    !
    !    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
    !    EMIN <= E <= EMAX.
    !
    !    Input argument 4 is a query about the base of real arithmetic:
    !
    !      IPMPAR(4) = B, the base of single and real ( kind = 8 ) arithmetic.
    !
    !    Input arguments 5 through 7 are queries about single precision
    !    floating point arithmetic:
    !
    !     IPMPAR(5) = M, the number of base B digits for single precision.
    !     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
    !     IPMPAR(7) = EMAX, the largest exponent E for single precision.
    !
    !    Input arguments 8 through 10 are queries about real ( kind = 8 )
    !    floating point arithmetic:
    !
    !     IPMPAR(8) = M, the number of base B digits for real ( kind = 8 ).
    !     IPMPAR(9) = EMIN, the smallest exponent E for real ( kind = 8 ).
    !     IPMPAR(10) = EMAX, the largest exponent E for real ( kind = 8 ).
    !
    !  Reference:
    !
    !    Phyllis Fox, Andrew Hall, Norman Schryer,
    !    Algorithm 528,
    !    Framework for a Portable FORTRAN Subroutine Library,
    !    ACM Transactions on Mathematical Software,
    !    Volume 4, 1978, pages 176-188.
    !
    !  Parameters:
    !
    !    Input, integer I, the index of the desired constant.
    !
    !    Output, integer IPMPAR, the value of the desired constant.
    !
    Implicit None

    Integer i
    Integer imach(10)
    Integer ipmpar
    !
    !     MACHINE CONSTANTS FOR AMDAHL MACHINES.
    !
    !     data imach( 1) /   2 /
    !     data imach( 2) /  31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /  16 /
    !     data imach( 5) /   6 /
    !     data imach( 6) / -64 /
    !     data imach( 7) /  63 /
    !     data imach( 8) /  14 /
    !     data imach( 9) / -64 /
    !     data imach(10) /  63 /
    !
    !     Machine constants for the AT&T 3B SERIES, AT&T
    !     PC 7300, AND AT&T 6300.
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    24 /
    !     data imach( 6) /  -125 /
    !     data imach( 7) /   128 /
    !     data imach( 8) /    53 /
    !     data imach( 9) / -1021 /
    !     data imach(10) /  1024 /
    !
    !     Machine constants for the BURROUGHS 1700 SYSTEM.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   33 /
    !     data imach( 3) / 8589934591 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   24 /
    !     data imach( 6) / -256 /
    !     data imach( 7) /  255 /
    !     data imach( 8) /   60 /
    !     data imach( 9) / -256 /
    !     data imach(10) /  255 /
    !
    !     Machine constants for the BURROUGHS 5700 SYSTEM.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   39 /
    !     data imach( 3) / 549755813887 /
    !     data imach( 4) /    8 /
    !     data imach( 5) /   13 /
    !     data imach( 6) /  -50 /
    !     data imach( 7) /   76 /
    !     data imach( 8) /   26 /
    !     data imach( 9) /  -50 /
    !     data imach(10) /   76 /
    !
    !     Machine constants for the BURROUGHS 6700/7700 SYSTEMS.
    !
    !     data imach( 1) /      2 /
    !     data imach( 2) /     39 /
    !     data imach( 3) / 549755813887 /
    !     data imach( 4) /      8 /
    !     data imach( 5) /     13 /
    !     data imach( 6) /    -50 /
    !     data imach( 7) /     76 /
    !     data imach( 8) /     26 /
    !     data imach( 9) / -32754 /
    !     data imach(10) /  32780 /
    !
    !     Machine constants for the CDC 6000/7000 SERIES
    !     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
    !     ARITHMETIC (NOS OPERATING SYSTEM).
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   48 /
    !     data imach( 3) / 281474976710655 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   48 /
    !     data imach( 6) / -974 /
    !     data imach( 7) / 1070 /
    !     data imach( 8) /   95 /
    !     data imach( 9) / -926 /
    !     data imach(10) / 1070 /
    !
    !     Machine constants for the CDC CYBER 995 64 BIT
    !     ARITHMETIC (NOS/VE OPERATING SYSTEM).
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    63 /
    !     data imach( 3) / 9223372036854775807 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    48 /
    !     data imach( 6) / -4096 /
    !     data imach( 7) /  4095 /
    !     data imach( 8) /    96 /
    !     data imach( 9) / -4096 /
    !     data imach(10) /  4095 /
    !
    !     Machine constants for the CRAY 1, XMP, 2, AND 3.
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    63 /
    !     data imach( 3) / 9223372036854775807 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    47 /
    !     data imach( 6) / -8189 /
    !     data imach( 7) /  8190 /
    !     data imach( 8) /    94 /
    !     data imach( 9) / -8099 /
    !     data imach(10) /  8190 /
    !
    !     Machine constants for the data GENERAL ECLIPSE S/200.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   15 /
    !     data imach( 3) / 32767 /
    !     data imach( 4) /   16 /
    !     data imach( 5) /    6 /
    !     data imach( 6) /  -64 /
    !     data imach( 7) /   63 /
    !     data imach( 8) /   14 /
    !     data imach( 9) /  -64 /
    !     data imach(10) /   63 /
    !
    !     Machine constants for the HARRIS 220.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   23 /
    !     data imach( 3) / 8388607 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   23 /
    !     data imach( 6) / -127 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   38 /
    !     data imach( 9) / -127 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the HONEYWELL 600/6000
    !     AND DPS 8/70 SERIES.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   35 /
    !     data imach( 3) / 34359738367 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   27 /
    !     data imach( 6) / -127 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   63 /
    !     data imach( 9) / -127 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the HP 2100
    !     3 WORD real ( kind = 8 ) OPTION WITH FTN4
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   15 /
    !     data imach( 3) / 32767 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   23 /
    !     data imach( 6) / -128 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   39 /
    !     data imach( 9) / -128 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the HP 2100
    !     4 WORD real ( kind = 8 ) OPTION WITH FTN4
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   15 /
    !     data imach( 3) / 32767 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   23 /
    !     data imach( 6) / -128 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   55 /
    !     data imach( 9) / -128 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the HP 9000.
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    24 /
    !     data imach( 6) /  -126 /
    !     data imach( 7) /   128 /
    !     data imach( 8) /    53 /
    !     data imach( 9) / -1021 /
    !     data imach(10) /  1024 /
    !
    !     Machine constants for the IBM 360/370 SERIES,
    !     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
    !     5/7/9 AND THE SEL SYSTEMS 85/86.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /   16 /
    !     data imach( 5) /    6 /
    !     data imach( 6) /  -64 /
    !     data imach( 7) /   63 /
    !     data imach( 8) /   14 /
    !     data imach( 9) /  -64 /
    !     data imach(10) /   63 /
    !
    !     Machine constants for the IBM PC.
    !
    !      data imach(1)/2/
    !      data imach(2)/31/
    !      data imach(3)/2147483647/
    !      data imach(4)/2/
    !      data imach(5)/24/
    !      data imach(6)/-125/
    !      data imach(7)/128/
    !      data imach(8)/53/
    !      data imach(9)/-1021/
    !      data imach(10)/1024/
    !
    !     Machine constants for the MACINTOSH II - ABSOFT
    !     MACFORTRAN II.
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    24 /
    !     data imach( 6) /  -125 /
    !     data imach( 7) /   128 /
    !     data imach( 8) /    53 /
    !     data imach( 9) / -1021 /
    !     data imach(10) /  1024 /
    !
    !     Machine constants for the MICROVAX - VMS FORTRAN.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   24 /
    !     data imach( 6) / -127 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   56 /
    !     data imach( 9) / -127 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the PDP-11 FORTRAN SUPPORTING
    !     32-BIT integer ARITHMETIC.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   24 /
    !     data imach( 6) / -127 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   56 /
    !     data imach( 9) / -127 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the SEQUENT BALANCE 8000.
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    24 /
    !     data imach( 6) /  -125 /
    !     data imach( 7) /   128 /
    !     data imach( 8) /    53 /
    !     data imach( 9) / -1021 /
    !     data imach(10) /  1024 /
    !
    !     Machine constants for the SILICON GRAPHICS IRIS-4D
    !     SERIES (MIPS R3000 PROCESSOR).
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    24 /
    !     data imach( 6) /  -125 /
    !     data imach( 7) /   128 /
    !     data imach( 8) /    53 /
    !     data imach( 9) / -1021 /
    !     data imach(10) /  1024 /
    !
    !     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
    !     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
    !     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
    !
    Data imach( 1) /     2 /
    Data imach( 2) /    31 /
    Data imach( 3) / 2147483647 /
    Data imach( 4) /     2 /
    Data imach( 5) /    24 /
    Data imach( 6) /  -125 /
    Data imach( 7) /   128 /
    Data imach( 8) /    53 /
    Data imach( 9) / -1021 /
    Data imach(10) /  1024 /
    !
    !     Machine constants for the UNIVAC 1100 SERIES.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   35 /
    !     data imach( 3) / 34359738367 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   27 /
    !     data imach( 6) / -128 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   60 /
    !     data imach( 9) /-1024 /
    !     data imach(10) / 1023 /
    !
    !     Machine constants for the VAX 11/780.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   24 /
    !     data imach( 6) / -127 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   56 /
    !     data imach( 9) / -127 /
    !     data imach(10) /  127 /
    !
    ipmpar = imach(i)

    Return
  End Function ipmpar

!!$  Subroutine negative_binomial_cdf_values ( n_data, f, s, p, cdf )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    Assume that a coin has a probability P of coming up heads on
!!$    !    any one trial.  Suppose that we plan to flip the coin until we
!!$    !    achieve a total of S heads.  If we let F represent the number of
!!$    !    tails that occur in this process, then the value of F satisfies
!!$    !    a negative binomial PDF:
!!$    !
!!$    !      PDF(F,S,P) = Choose ( F from F+S-1 ) * P**S * (1-P)**F
!!$    !
!!$    !    The negative binomial CDF is the probability that there are F or
!!$    !    fewer failures upon the attainment of the S-th success.  Thus,
!!$    !
!!$    !      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    07 June 2004
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    FC Powell,
!!$    !    Statistical Tables for Sociology, Biology and Physical Sciences,
!!$    !    Cambridge University Press, 1982.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, integer F, the maximum number of failures.
!!$    !
!!$    !    Output, integer S, the number of successes.
!!$    !
!!$    !    Output, real ( kind = 8 ) P, the probability of a success on one trial.
!!$    !
!!$    !    Output, real ( kind = 8 ) CDF, the probability of at most F failures
!!$    !    before the S-th success.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 27
!!$
!!$    Real ( kind = 8 ) cdf
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: cdf_vec = (/ &
!!$         0.6367D+00, 0.3633D+00, 0.1445D+00, &
!!$         0.5000D+00, 0.2266D+00, 0.0625D+00, &
!!$         0.3438D+00, 0.1094D+00, 0.0156D+00, &
!!$         0.1792D+00, 0.0410D+00, 0.0041D+00, &
!!$         0.0705D+00, 0.0109D+00, 0.0007D+00, &
!!$         0.9862D+00, 0.9150D+00, 0.7472D+00, &
!!$         0.8499D+00, 0.5497D+00, 0.2662D+00, &
!!$         0.6513D+00, 0.2639D+00, 0.0702D+00, &
!!$         1.0000D+00, 0.0199D+00, 0.0001D+00 /)
!!$    Integer f
!!$    Integer, Save, Dimension ( n_max ) :: f_vec = (/ &
!!$         4,  3,  2, &
!!$         3,  2,  1, &
!!$         2,  1,  0, &
!!$         2,  1,  0, &
!!$         2,  1,  0, &
!!$         11, 10,  9, &
!!$         17, 16, 15, &
!!$         9,  8,  7, &
!!$         2,  1,  0 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) p
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: p_vec = (/ &
!!$         0.50D+00, 0.50D+00, 0.50D+00, &
!!$         0.50D+00, 0.50D+00, 0.50D+00, &
!!$         0.50D+00, 0.50D+00, 0.50D+00, &
!!$         0.40D+00, 0.40D+00, 0.40D+00, &
!!$         0.30D+00, 0.30D+00, 0.30D+00, &
!!$         0.30D+00, 0.30D+00, 0.30D+00, &
!!$         0.10D+00, 0.10D+00, 0.10D+00, &
!!$         0.10D+00, 0.10D+00, 0.10D+00, &
!!$         0.01D+00, 0.01D+00, 0.01D+00 /)
!!$    Integer s
!!$    Integer, Save, Dimension ( n_max ) :: s_vec = (/ &
!!$         4, 5, 6, &
!!$         4, 5, 6, &
!!$         4, 5, 6, &
!!$         4, 5, 6, &
!!$         4, 5, 6, &
!!$         1, 2, 3, &
!!$         1, 2, 3, &
!!$         1, 2, 3, &
!!$         0, 1, 2 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       f = 0
!!$       s = 0
!!$       p = 0.0D+00
!!$       cdf = 0.0D+00
!!$    Else
!!$       f = f_vec(n_data)
!!$       s = s_vec(n_data)
!!$       p = p_vec(n_data)
!!$       cdf = cdf_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine negative_binomial_cdf_values

!!$  Subroutine normal_01_cdf_values ( n_data, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    In Mathematica, the function can be evaluated by:
!!$    !
!!$    !      Needs["Statistics`ContinuousDistributions`"]
!!$    !      dist = NormalDistribution [ 0, 1 ]
!!$    !      CDF [ dist, x ]
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    28 August 2004
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !    Stephen Wolfram,
!!$    !    The Mathematica Book,
!!$    !    Fourth Edition,
!!$    !    Wolfram Media / Cambridge University Press, 1999.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, real ( kind = 8 ) X, the argument of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 17
!!$
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.5000000000000000D+00, &
!!$         0.5398278372770290D+00, &
!!$         0.5792597094391030D+00, &
!!$         0.6179114221889526D+00, &
!!$         0.6554217416103242D+00, &
!!$         0.6914624612740131D+00, &
!!$         0.7257468822499270D+00, &
!!$         0.7580363477769270D+00, &
!!$         0.7881446014166033D+00, &
!!$         0.8159398746532405D+00, &
!!$         0.8413447460685429D+00, &
!!$         0.9331927987311419D+00, &
!!$         0.9772498680518208D+00, &
!!$         0.9937903346742239D+00, &
!!$         0.9986501019683699D+00, &
!!$         0.9997673709209645D+00, &
!!$         0.9999683287581669D+00 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         0.0000000000000000D+00, &
!!$         0.1000000000000000D+00, &
!!$         0.2000000000000000D+00, &
!!$         0.3000000000000000D+00, &
!!$         0.4000000000000000D+00, &
!!$         0.5000000000000000D+00, &
!!$         0.6000000000000000D+00, &
!!$         0.7000000000000000D+00, &
!!$         0.8000000000000000D+00, &
!!$         0.9000000000000000D+00, &
!!$         0.1000000000000000D+01, &
!!$         0.1500000000000000D+01, &
!!$         0.2000000000000000D+01, &
!!$         0.2500000000000000D+01, &
!!$         0.3000000000000000D+01, &
!!$         0.3500000000000000D+01, &
!!$         0.4000000000000000D+01 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine normal_01_cdf_values

!!$  Subroutine normal_cdf_values ( n_data, mu, sigma, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! NORMAL_CDF_VALUES returns some values of the Normal CDF.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    In Mathematica, the function can be evaluated by:
!!$    !
!!$    !      Needs["Statistics`ContinuousDistributions`"]
!!$    !      dist = NormalDistribution [ mu, sigma ]
!!$    !      CDF [ dist, x ]
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    05 August 2004
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !    Stephen Wolfram,
!!$    !    The Mathematica Book,
!!$    !    Fourth Edition,
!!$    !    Wolfram Media / Cambridge University Press, 1999.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, real ( kind = 8 ) MU, the mean of the distribution.
!!$    !
!!$    !    Output, real ( kind = 8 ) SIGMA, the variance of the distribution.
!!$    !
!!$    !    Output, real ( kind = 8 ) X, the argument of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 12
!!$
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.5000000000000000D+00, &
!!$         0.9772498680518208D+00, &
!!$         0.9999683287581669D+00, &
!!$         0.9999999990134124D+00, &
!!$         0.6914624612740131D+00, &
!!$         0.6305586598182364D+00, &
!!$         0.5987063256829237D+00, &
!!$         0.5792597094391030D+00, &
!!$         0.6914624612740131D+00, &
!!$         0.5000000000000000D+00, &
!!$         0.3085375387259869D+00, &
!!$         0.1586552539314571D+00 /)
!!$    Real ( kind = 8 ) mu
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: mu_vec = (/ &
!!$         0.1000000000000000D+01, &
!!$         0.1000000000000000D+01, &
!!$         0.1000000000000000D+01, &
!!$         0.1000000000000000D+01, &
!!$         0.1000000000000000D+01, &
!!$         0.1000000000000000D+01, &
!!$         0.1000000000000000D+01, &
!!$         0.1000000000000000D+01, &
!!$         0.2000000000000000D+01, &
!!$         0.3000000000000000D+01, &
!!$         0.4000000000000000D+01, &
!!$         0.5000000000000000D+01 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) sigma
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: sigma_vec = (/ &
!!$         0.5000000000000000D+00, &
!!$         0.5000000000000000D+00, &
!!$         0.5000000000000000D+00, &
!!$         0.5000000000000000D+00, &
!!$         0.2000000000000000D+01, &
!!$         0.3000000000000000D+01, &
!!$         0.4000000000000000D+01, &
!!$         0.5000000000000000D+01, &
!!$         0.2000000000000000D+01, &
!!$         0.2000000000000000D+01, &
!!$         0.2000000000000000D+01, &
!!$         0.2000000000000000D+01 /)
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         0.1000000000000000D+01, &
!!$         0.2000000000000000D+01, &
!!$         0.3000000000000000D+01, &
!!$         0.4000000000000000D+01, &
!!$         0.2000000000000000D+01, &
!!$         0.2000000000000000D+01, &
!!$         0.2000000000000000D+01, &
!!$         0.2000000000000000D+01, &
!!$         0.3000000000000000D+01, &
!!$         0.3000000000000000D+01, &
!!$         0.3000000000000000D+01, &
!!$         0.3000000000000000D+01 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       mu = 0.0D+00
!!$       sigma = 0.0D+00
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       mu = mu_vec(n_data)
!!$       sigma = sigma_vec(n_data)
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine normal_cdf_values

!!$  Subroutine poisson_cdf_values ( n_data, a, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! POISSON_CDF_VALUES returns some values of the Poisson CDF.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    CDF(X)(A) is the probability of at most X successes in unit time,
!!$    !    given that the expected mean number of successes is A.
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    28 May 2001
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !    Daniel Zwillinger,
!!$    !    CRC Standard Mathematical Tables and Formulae,
!!$    !    30th Edition, CRC Press, 1996, pages 653-658.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, real ( kind = 8 ) A, integer X, the arguments of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 21
!!$
!!$    Real ( kind = 8 ) a
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: a_vec = (/ &
!!$         0.02D+00, 0.10D+00, 0.10D+00, 0.50D+00, &
!!$         0.50D+00, 0.50D+00, 1.00D+00, 1.00D+00, &
!!$         1.00D+00, 1.00D+00, 2.00D+00, 2.00D+00, &
!!$         2.00D+00, 2.00D+00, 5.00D+00, 5.00D+00, &
!!$         5.00D+00, 5.00D+00, 5.00D+00, 5.00D+00, &
!!$         5.00D+00 /)
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.980D+00, 0.905D+00, 0.995D+00, 0.607D+00, &
!!$         0.910D+00, 0.986D+00, 0.368D+00, 0.736D+00, &
!!$         0.920D+00, 0.981D+00, 0.135D+00, 0.406D+00, &
!!$         0.677D+00, 0.857D+00, 0.007D+00, 0.040D+00, &
!!$         0.125D+00, 0.265D+00, 0.441D+00, 0.616D+00, &
!!$         0.762D+00 /)
!!$    Integer n_data
!!$    Integer x
!!$    Integer, Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         0, 0, 1, 0, &
!!$         1, 2, 0, 1, &
!!$         2, 3, 0, 1, &
!!$         2, 3, 0, 1, &
!!$         2, 3, 4, 5, &
!!$         6 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       a = 0.0D+00
!!$       x = 0
!!$       fx = 0.0D+00
!!$    Else
!!$       a = a_vec(n_data)
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine poisson_cdf_values

  Function psi ( xx )

    !*****************************************************************************80
    !
    !! PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
    !
    !  Discussion:
    !
    !    The main computation involves evaluation of rational Chebyshev
    !    approximations.  PSI was written at Argonne National Laboratory
    !    for FUNPACK, and subsequently modified by A. H. Morris of NSWC.
    !
    !  Reference:
    !
    !    William Cody, Anthony Strecok, Henry Thacher,
    !    Chebyshev Approximations for the Psi Function,
    !    Mathematics of Computation,
    !    Volume 27, 1973, pages 123-127.
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) XX, the argument of the psi function.
    !
    !    Output, real ( kind = 8 ) PSI, the value of the psi function.  PSI
    !    is assigned the value 0 when the psi function is undefined.
    !
    Implicit None

    Real ( kind = 8 ) aug
    Real ( kind = 8 ) den
    Real ( kind = 8 ), Parameter :: dx0 = &
         1.461632144968362341262659542325721325D+00
    Integer i
    !    Integer ipmpar
    Integer m
    Integer n
    Integer nq
    Real ( kind = 8 ), Parameter, Dimension ( 7 ) :: p1 = (/ &
         0.895385022981970D-02, &
         0.477762828042627D+01, &
         0.142441585084029D+03, &
         0.118645200713425D+04, &
         0.363351846806499D+04, &
         0.413810161269013D+04, &
         0.130560269827897D+04/)
    Real ( kind = 8 ), Dimension ( 4 ) :: p2 = (/ &
         -0.212940445131011D+01, &
         -0.701677227766759D+01, &
         -0.448616543918019D+01, &
         -0.648157123766197D+00 /)
    Real ( kind = 8 ), Parameter :: piov4 = 0.785398163397448D+00
    Real ( kind = 8 ) psi
    !
    !  Coefficients for rational approximation of
    !  PSI(X) / (X - X0),  0.5D+00 <= X <= 3.0D+00
    !
    Real ( kind = 8 ), Dimension ( 6 ) :: q1 = (/ &
         0.448452573429826D+02, &
         0.520752771467162D+03, &
         0.221000799247830D+04, &
         0.364127349079381D+04, &
         0.190831076596300D+04, &
         0.691091682714533D-05 /)
    Real ( kind = 8 ), Dimension ( 4 ) :: q2 = (/ &
         0.322703493791143D+02, &
         0.892920700481861D+02, &
         0.546117738103215D+02, &
         0.777788548522962D+01 /)
    Real ( kind = 8 ) sgn
    Real ( kind = 8 ) upper
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x
    Real ( kind = 8 ) xmax1
    Real ( kind = 8 ) xmx0
    Real ( kind = 8 ) xsmall
    Real ( kind = 8 ) xx
    Real ( kind = 8 ) z
    !
    !  XMAX1 is the largest positive floating point constant with entirely
    !  integer representation.  It is also used as negative of lower bound
    !  on acceptable negative arguments and as the positive argument beyond which
    !  psi may be represented as LOG(X).
    !
    xmax1 = Real ( ipmpar(3), kind = 8 )
    xmax1 = Min ( xmax1, 1.0D+00 / Epsilon ( xmax1 ) )
    !
    !  XSMALL is the absolute argument below which PI*COTAN(PI*X)
    !  may be represented by 1/X.
    !
    xsmall = 1.0D-09

    x = xx
    aug = 0.0D+00

    If ( x == 0.0D+00 ) Then
       psi = 0.0D+00
       Return
    End If
    !
    !  X < 0.5,  Use reflection formula PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
    !
    If ( x < 0.5D+00 ) Then
       !
       !  0 < ABS ( X ) <= XSMALL.  Use 1/X as a substitute for PI*COTAN(PI*X)
       !
       If ( Abs ( x ) <= xsmall ) Then
          aug = -1.0D+00 / x
          go to 40
       End If
       !
       !  Reduction of argument for cotangent.
       !
       w = -x
       sgn = piov4

       If ( w <= 0.0D+00 ) Then
          w = -w
          sgn = -sgn
       End If
       !
       !  Make an error exit if X <= -XMAX1
       !
       If ( xmax1 <= w ) Then
          psi = 0.0D+00
          Return
       End If

       nq = Int ( w )
       w = w - Real ( nq, kind = 8 )
       nq = Int ( w * 4.0D+00 )
       w = 4.0D+00 * ( w - Real ( nq, kind = 8 ) * 0.25D+00 )
       !
       !  W is now related to the fractional part of 4.0D+00 * X.
       !  Adjust argument to correspond to values in first
       !  quadrant and determine sign.
       !
       n = nq / 2
       If ( n + n /= nq ) Then
          w = 1.0D+00 - w
       End If

       z = piov4 * w
       m = n / 2

       If ( m + m /= n ) Then
          sgn = -sgn
       End If
       !
       !  Determine final value for -PI * COTAN(PI*X).
       !
       n = ( nq + 1 ) / 2
       m = n / 2
       m = m + m

       If ( m == n ) Then

          If ( z == 0.0D+00 ) Then
             psi = 0.0D+00
             Return
          End If

          aug = 4.0D+00 * sgn * ( Cos(z) / Sin(z) )

       Else

          aug = 4.0D+00 * sgn * ( Sin(z) / Cos(z) )

       End If

40     Continue

       x = 1.0D+00 - x

    End If
    !
    !  0.5 <= X <= 3
    !
    If ( x <= 3.0D+00 ) Then

       den = x
       upper = p1(1) * x

       Do i = 1, 5
          den = ( den + q1(i) ) * x
          upper = ( upper + p1(i+1) ) * x
       End Do

       den = ( upper + p1(7) ) / ( den + q1(6) )
       xmx0 = Real ( x, kind = 8 ) - dx0
       psi = den * xmx0 + aug
       !
       !  3 < X < XMAX1
       !
    Else If ( x < xmax1 ) Then

       w = 1.0D+00 / x**2
       den = w
       upper = p2(1) * w

       Do i = 1, 3
          den = ( den + q2(i) ) * w
          upper = ( upper + p2(i+1) ) * w
       End Do

       aug = upper / ( den + q2(4) ) - 0.5D+00 / x + aug
       psi = aug + Log ( x )
       !
       !  XMAX1 <= X
       !
    Else

       psi = aug + Log ( x )

    End If

    Return
  End Function psi

!!$  Subroutine psi_values ( n_data, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! PSI_VALUES returns some values of the Psi or Digamma function.
!!$    !
!!$    !  Discussion:
!!$    !
!!$    !    PSI(X) = d LN ( GAMMA ( X ) ) / d X = GAMMA'(X) / GAMMA(X)
!!$    !
!!$    !    PSI(1) = - Euler's constant.
!!$    !
!!$    !    PSI(X+1) = PSI(X) + 1 / X.
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    17 May 2001
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, real ( kind = 8 ) X, the argument of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 11
!!$
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         -0.5772156649D+00, -0.4237549404D+00, -0.2890398966D+00, &
!!$         -0.1691908889D+00, -0.0613845446D+00, -0.0364899740D+00, &
!!$         0.1260474528D+00,  0.2085478749D+00,  0.2849914333D+00, &
!!$         0.3561841612D+00,  0.4227843351D+00 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         1.0D+00,  1.1D+00,  1.2D+00,  &
!!$         1.3D+00,  1.4D+00,  1.5D+00,  &
!!$         1.6D+00,  1.7D+00,  1.8D+00,  &
!!$         1.9D+00,  2.0D+00 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine psi_values

  Subroutine r8_swap ( x, y )

    !*****************************************************************************80
    !
    !! R8_SWAP swaps two R8 values.
    !
    !  Modified:
    !
    !    01 May 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
    !    Y have been interchanged.
    !
    Implicit None

    Real ( kind = 8 ) x
    Real ( kind = 8 ) y
    Real ( kind = 8 ) z

    z = x
    x = y
    y = z

    Return
  End Subroutine r8_swap

  Function rcomp ( a, x )

    !*****************************************************************************80
    !
    !! RCOMP evaluates exp(-X) * X**A / Gamma(A).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, X, arguments of the quantity to be computed.
    !
    !    Output, real ( kind = 8 ) RCOMP, the value of exp(-X) * X**A / Gamma(A).
    !
    !  Local parameters:
    !
    !    RT2PIN = 1/SQRT(2*PI)
    !
    Implicit None

    Real ( kind = 8 ) a
    !    Real ( kind = 8 ) gam1
    !    Real ( kind = 8 ) gamma
    Real ( kind = 8 ) rcomp
    !    Real ( kind = 8 ) rlog
    Real ( kind = 8 ), Parameter :: rt2pin = 0.398942280401433D+00
    Real ( kind = 8 ) t
    Real ( kind = 8 ) t1
    Real ( kind = 8 ) u
    Real ( kind = 8 ) x

    If ( a < 20.0D+00 ) Then

       t = a * Log ( x ) - x

       If ( a < 1.0D+00 ) Then
          rcomp = ( a * Exp ( t ) ) * ( 1.0D+00 + gam1 ( a ) )
       Else
          rcomp = Exp ( t ) / gamma ( a )
       End If

    Else

       u = x / a

       If ( u == 0.0D+00 ) Then
          rcomp = 0.0D+00
       Else
          t = ( 1.0D+00 / a )**2
          t1 = ((( 0.75D+00 * t - 1.0D+00 ) * t + 3.5D+00 ) * t - 105.0D+00 ) &
               / ( a * 1260.0D+00 )
          t1 = t1 - a * rlog ( u )
          rcomp = rt2pin * Sqrt ( a ) * Exp ( t1 )
       End If

    End If

    Return
  End Function rcomp

  Function rexp ( x )

    !*****************************************************************************80
    !
    !! REXP evaluates the function EXP(X) - 1.
    !
    !  Modified:
    !
    !    09 December 1999
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the function.
    !
    !    Output, real ( kind = 8 ) REXP, the value of EXP(X)-1.
    !
    Implicit None

    Real ( kind = 8 ), Parameter :: p1 =  0.914041914819518D-09
    Real ( kind = 8 ), Parameter :: p2 =  0.238082361044469D-01
    Real ( kind = 8 ), Parameter :: q1 = -0.499999999085958D+00
    Real ( kind = 8 ), Parameter :: q2 =  0.107141568980644D+00
    Real ( kind = 8 ), Parameter :: q3 = -0.119041179760821D-01
    Real ( kind = 8 ), Parameter :: q4 =  0.595130811860248D-03
    Real ( kind = 8 ) rexp
    Real ( kind = 8 ) w
    Real ( kind = 8 ) x

    If ( Abs ( x ) <= 0.15D+00 ) Then

       rexp = x * ( ( ( p2 * x + p1 ) * x + 1.0D+00 ) &
            / ( ( ( ( q4 * x + q3 ) * x + q2 ) * x + q1 ) * x + 1.0D+00 ) )

    Else

       w = Exp ( x )

       If ( x <= 0.0D+00 ) Then
          rexp = ( w - 0.5D+00 ) - 0.5D+00
       Else
          rexp = w * ( 0.5D+00 + ( 0.5D+00 - 1.0D+00 / w ) )
       End If

    End If

    Return
  End Function rexp

  Function rlog ( x )

    !*****************************************************************************80
    !
    !! RLOG computes X - 1 - LN(X).
    !
    !  Modified:
    !
    !    06 August 2004
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the function.
    !
    !    Output, real ( kind = 8 ) RLOG, the value of the function.
    !
    Implicit None

    Real ( kind = 8 ), Parameter :: a  =  0.566749439387324D-01
    Real ( kind = 8 ), Parameter :: b  =  0.456512608815524D-01
    Real ( kind = 8 ), Parameter :: half = 0.5D+00
    Real ( kind = 8 ), Parameter :: p0 =  0.333333333333333D+00
    Real ( kind = 8 ), Parameter :: p1 = -0.224696413112536D+00
    Real ( kind = 8 ), Parameter :: p2 =  0.620886815375787D-02
    Real ( kind = 8 ), Parameter :: q1 = -0.127408923933623D+01
    Real ( kind = 8 ), Parameter :: q2 =  0.354508718369557D+00
    Real ( kind = 8 ) r
    Real ( kind = 8 ) rlog
    Real ( kind = 8 ) t
    Real ( kind = 8 ), Parameter :: two =  2.0D+00
    Real ( kind = 8 ) u
    Real ( kind = 8 ) w
    Real ( kind = 8 ) w1
    Real ( kind = 8 ) x

    If ( x < 0.61D+00 ) Then

       r = ( x - 0.5D+00 ) - 0.5D+00
       rlog = r - Log ( x )

    Else If ( x < 1.57D+00 ) Then

       If ( x < 0.82D+00 ) Then

          u = x - 0.7D+00
          u = u / 0.7D+00
          w1 = a - u * 0.3D+00

       Else If ( x < 1.18D+00 ) Then

          u = ( x - half ) - half
          w1 = 0.0D+00

       Else If ( x < 1.57D+00 ) Then

          u = 0.75D+00 * x - 1.0D+00
          w1 = b + u / 3.0D+00

       End If

       r = u / ( u + two )
       t = r * r
       w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
       rlog = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

    Else If ( 1.57D+00 <= x ) Then

       r = ( x - half ) - half
       rlog = r - Log ( x )

    End If

    Return
  End Function rlog

  Function rlog1 ( x )

    !*****************************************************************************80
    !
    !! RLOG1 evaluates the function X - ln ( 1 + X ).
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument.
    !
    !    Output, real ( kind = 8 ) RLOG1, the value of X - ln ( 1 + X ).
    !
    Implicit None

    Real ( kind = 8 ), Parameter :: a = 0.566749439387324D-01
    Real ( kind = 8 ), Parameter :: b = 0.456512608815524D-01
    Real ( kind = 8 ) h
    Real ( kind = 8 ), Parameter :: half = 0.5D+00
    Real ( kind = 8 ), Parameter :: p0 = 0.333333333333333D+00
    Real ( kind = 8 ), Parameter :: p1 = -0.224696413112536D+00
    Real ( kind = 8 ), Parameter :: p2 = 0.620886815375787D-02
    Real ( kind = 8 ), Parameter :: q1 = -0.127408923933623D+01
    Real ( kind = 8 ), Parameter :: q2 = 0.354508718369557D+00
    Real ( kind = 8 ) r
    Real ( kind = 8 ) rlog1
    Real ( kind = 8 ) t
    Real ( kind = 8 ), Parameter :: two =  2.0D+00
    Real ( kind = 8 ) w
    Real ( kind = 8 ) w1
    Real ( kind = 8 ) x

    If ( x < -0.39D+00 ) Then

       w = ( x + half ) + half
       rlog1 = x - Log ( w )

    Else If ( x < -0.18D+00 ) Then

       h = x + 0.3D+00
       h = h / 0.7D+00
       w1 = a - h * 0.3D+00

       r = h / ( h + 2.0D+00 )
       t = r * r
       w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
       rlog1 = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

    Else If ( x <= 0.18D+00 ) Then

       h = x
       w1 = 0.0D+00

       r = h / ( h + two )
       t = r * r
       w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
       rlog1 = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

    Else If ( x <= 0.57D+00 ) Then

       h = 0.75D+00 * x - 0.25D+00
       w1 = b + h / 3.0D+00

       r = h / ( h + 2.0D+00 )
       t = r * r
       w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
       rlog1 = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

    Else

       w = ( x + half ) + half
       rlog1 = x - Log ( w )

    End If

    Return
  End Function rlog1

!!$  Subroutine student_cdf_values ( n_data, a, x, fx )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! STUDENT_CDF_VALUES returns some values of the Student CDF.
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    02 June 2001
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Reference:
!!$    !
!!$    !    Milton Abramowitz, Irene Stegun,
!!$    !    Handbook of Mathematical Functions,
!!$    !    US Department of Commerce, 1964.
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!!$    !    first call.  On each call, the routine increments N_DATA by 1, and
!!$    !    returns the corresponding data; when there is no more data, the
!!$    !    output value of N_DATA will be 0 again.
!!$    !
!!$    !    Output, integer A, real ( kind = 8 ) X, the arguments of the function.
!!$    !
!!$    !    Output, real ( kind = 8 ) FX, the value of the function.
!!$    !
!!$    Implicit None
!!$
!!$    Integer, Parameter :: n_max = 13
!!$
!!$    Integer a
!!$    Integer, Save, Dimension ( n_max ) :: a_vec = (/ &
!!$         1, 2, 3, 4, &
!!$         5, 2, 5, 2, &
!!$         5, 2, 3, 4, &
!!$         5 /)
!!$    Real ( kind = 8 ) fx
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: fx_vec = (/ &
!!$         0.60D+00, 0.60D+00, 0.60D+00, 0.60D+00, &
!!$         0.60D+00, 0.75D+00, 0.75D+00, 0.95D+00, &
!!$         0.95D+00, 0.99D+00, 0.99D+00, 0.99D+00, &
!!$         0.99D+00 /)
!!$    Integer n_data
!!$    Real ( kind = 8 ) x
!!$    Real ( kind = 8 ), Save, Dimension ( n_max ) :: x_vec = (/ &
!!$         0.325D+00, 0.289D+00, 0.277D+00, 0.271D+00, &
!!$         0.267D+00, 0.816D+00, 0.727D+00, 2.920D+00, &
!!$         2.015D+00, 6.965D+00, 4.541D+00, 3.747D+00, &
!!$         3.365D+00 /)
!!$
!!$    If ( n_data < 0 ) Then
!!$       n_data = 0
!!$    End If
!!$
!!$    n_data = n_data + 1
!!$
!!$    If ( n_max < n_data ) Then
!!$       n_data = 0
!!$       a = 0
!!$       x = 0.0D+00
!!$       fx = 0.0D+00
!!$    Else
!!$       a = a_vec(n_data)
!!$       x = x_vec(n_data)
!!$       fx = fx_vec(n_data)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine student_cdf_values

  Function stvaln ( p )

    !*****************************************************************************80
    !
    !! STVALN provides starting values for the inverse of the normal distribution.
    !
    !  Discussion:
    !
    !    The routine returns an X for which it is approximately true that
    !      P = CUMNOR(X),
    !    that is,
    !      P = Integral ( -infinity < U <= X ) exp(-U*U/2)/sqrt(2*PI) dU.
    !
    !  Reference:
    !
    !    William Kennedy, James Gentle,
    !    Statistical Computing,
    !    Marcel Dekker, NY, 1980, page 95,
    !    QA276.4 K46
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) P, the probability whose normal deviate
    !    is sought.
    !
    !    Output, real ( kind = 8 ) STVALN, the normal deviate whose probability
    !    is approximately P.
    !
    Implicit None

    !    Real ( kind = 8 ) eval_pol
    Real ( kind = 8 ) p
    Real ( kind = 8 ) sgn
    Real ( kind = 8 ) stvaln
    Real ( kind = 8 ), Parameter, Dimension(0:4) :: xden = (/ &
         0.993484626060D-01, &
         0.588581570495D+00, &
         0.531103462366D+00, &
         0.103537752850D+00, &
         0.38560700634D-02 /)
    Real ( kind = 8 ), Parameter, Dimension(0:4) :: xnum = (/ &
         -0.322232431088D+00, &
         -1.000000000000D+00, &
         -0.342242088547D+00, &
         -0.204231210245D-01, &
         -0.453642210148D-04 /)
    Real ( kind = 8 ) y
    Real ( kind = 8 ) z

    If ( p <= 0.5D+00 ) Then

       sgn = -1.0D+00
       z = p

    Else

       sgn = 1.0D+00
       z = 1.0D+00 - p

    End If

    y = Sqrt ( -2.0D+00 * Log ( z ) )
    stvaln = y + eval_pol ( xnum, 4, y ) / eval_pol ( xden, 4, y )
    stvaln = sgn * stvaln

    Return
  End Function stvaln

!!$  Subroutine timestamp ( )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
!!$    !
!!$    !  Example:
!!$    !
!!$    !    May 31 2001   9:45:54.872 AM
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    15 March 2003
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    None
!!$    !
!!$    Implicit None
!!$
!!$    Character ( len = 40 ) string
!!$
!!$    Call timestring ( string )
!!$
!!$    Write ( *, '(a)' ) Trim ( string )
!!$
!!$    Return
!!$  End Subroutine timestamp
!!$
!!$  Subroutine timestring ( string )
!!$
!!$    !*****************************************************************************80
!!$    !
!!$    !! TIMESTRING writes the current YMDHMS date into a string.
!!$    !
!!$    !  Example:
!!$    !
!!$    !    STRING = 'May 31 2001   9:45:54.872 AM'
!!$    !
!!$    !  Modified:
!!$    !
!!$    !    15 March 2003
!!$    !
!!$    !  Author:
!!$    !
!!$    !    John Burkardt
!!$    !
!!$    !  Parameters:
!!$    !
!!$    !    Output, character ( len = * ) STRING, contains the date information.
!!$    !    A character length of 40 should always be sufficient.
!!$    !
!!$    Implicit None
!!$
!!$    Character ( len = 8 ) ampm
!!$    Integer d
!!$    Character ( len = 8 ) date
!!$    Integer h
!!$    Integer m
!!$    Integer mm
!!$    Character ( len = 9 ), Parameter, Dimension(12) :: month = (/ &
!!$         'January  ', 'February ', 'March    ', 'April    ', &
!!$         'May      ', 'June     ', 'July     ', 'August   ', &
!!$         'September', 'October  ', 'November ', 'December ' /)
!!$    Integer n
!!$    Integer s
!!$    Character ( len = * ) string
!!$    Character ( len = 10 ) time
!!$    Integer values(8)
!!$    Integer y
!!$    Character ( len = 5 ) zone
!!$
!!$    Call Date_and_time ( date, time, zone, values )
!!$
!!$    y = values(1)
!!$    m = values(2)
!!$    d = values(3)
!!$    h = values(5)
!!$    n = values(6)
!!$    s = values(7)
!!$    mm = values(8)
!!$
!!$    If ( h < 12 ) Then
!!$       ampm = 'AM'
!!$    Else If ( h == 12 ) Then
!!$       If ( n == 0 .And. s == 0 ) Then
!!$          ampm = 'Noon'
!!$       Else
!!$          ampm = 'PM'
!!$       End If
!!$    Else
!!$       h = h - 12
!!$       If ( h < 12 ) Then
!!$          ampm = 'PM'
!!$       Else If ( h == 12 ) Then
!!$          If ( n == 0 .And. s == 0 ) Then
!!$             ampm = 'Midnight'
!!$          Else
!!$             ampm = 'AM'
!!$          End If
!!$       End If
!!$    End If
!!$
!!$    Write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
!!$         Trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, Trim ( ampm )
!!$
!!$    Return
!!$  End Subroutine timestring

End Module DCDFLIB

!*********************************************************
!*****   FDS5-Evac Statistical subroutine library    *****
!*********************************************************
!
! VTT Technical Research Centre of Finland            2006
!
! Author:  Timo Korhonen
! Date:    22.5.2006
! Version: 4.99
!
! Uses: dcdflib.f90, which is included in this file.
!
!*********************************************************
! File contains following subroutines:
!
! * RandomNumbers
!   - generates random numbers from the specified distributions
!   - Many subprograms are included: different distributions for
!     MC-sampling (xxx_rnd).
!   - User given density: now a discrete density can also be given
!
!
!*********************************************************

Module STAT

  Use DCDFLIB ! ieva.f90
  Use PRECISION_PARAMETERS ! prec.f90g
  Use GLOBAL_CONSTANTS, GC_Gamma => Gamma    ! cons.f90
  Use COMP_FUNCTIONS, Only: SHUTDOWN         ! func.f90
  Use MATH_FUNCTIONS, Only: AFILL ! func.f90
  Use MEMORY_FUNCTIONS, Only: ChkMemErr ! func.f90

  Implicit None
  !
  Private
  ! Public subprograms (called from the main program)
  Public RandomNumbers, UserRandomNumbers
  !
  Integer :: IZERO

Contains
  !
  !
  !*********************************************************
  !***** Subroutine RandomNumbers Starts               *****
  !*********************************************************
  !
  ! INPUTS:
  ! n_rnd: Number of random numbers returned
  ! n_par: how many parameters for the distribution
  ! RandomPara(n_par): Parameters of the distribution
  ! RandomType:   1) Uniform, 2) ...
  !
  ! OUTPUTS:
  ! RandomVector(nmc): The random sample
  !
  Subroutine RandomNumbers( n_rnd, n_par, RandomType, RandomPara, rnd_vec)
    Implicit None

    !----------------------------------------------
    ! Argument definitions
    !----------------------------------------------
    ! Input variables (by value)
    Integer, Intent(In) :: n_rnd, n_par, RandomType
    Real(EB), Intent(In) :: RandomPara(n_par)
    ! Output variables (by reference)
    Real(EB), Intent(Out) :: rnd_vec(n_rnd)

    !----------------------------------------------
    ! Local variables
    !----------------------------------------------
    !----------------------------------------------
    !
    ! 1: uniform (TESTED: OK)
    ! 2==>6: beta (TESTED: OK)
    ! 3: gamma  (TESTED: OK)
    ! 4: normal (TESTED: OK)
    ! 5: lognormal (TESTED: OK)
    !    17th April 2009: high end truncated and shifted lognormal (TESTED: OK)
    ! 6==>2: Truncated normal (TESTED: OK)
    ! 7: Triangular (TESTED: OK)
    ! 8: Weibull (TESTED: OK) (alpha=1: Exponential)
    ! 9: Gumbel (TESTED: OK)

    Select Case (RandomType)
    Case (1)    ! Uniform distribution
       ! Parameters: (ave,min,max) ave not used
       Call Uniform_rnd(n_rnd, RandomPara(2), RandomPara(3), rnd_vec)
    Case (2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       Call TNormal_rnd(n_rnd, RandomPara(1), RandomPara(2), RandomPara(3), RandomPara(4), rnd_vec)
    Case (3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       Call Gamma_rnd(n_rnd, RandomPara(2), RandomPara(3), rnd_vec)
    Case (4)   ! Normal
       ! Parameters: (ave,sigma)
       Call Normal_rnd(n_rnd, RandomPara(1), RandomPara(2), rnd_vec)
    Case (5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       ! 17th April 2009: para3=cutoff high side, para4=shift
       Call LogNormal_rnd(n_rnd, RandomPara(1), RandomPara(2), RandomPara(3), RandomPara(4), rnd_vec)
    Case (6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       Call Beta_rnd(n_rnd, RandomPara(2), RandomPara(3), rnd_vec)
    Case (7)   ! Triangular
       ! Parameters: (peak,min,max)
       Call Triang_rnd(n_rnd, RandomPara(1), RandomPara(2), RandomPara(3), rnd_vec)
    Case(8)    ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda) ave not used
       Call Weibull_rnd(n_rnd, RandomPara(2), RandomPara(3), rnd_vec)
    Case(9)    ! Gumbel
       ! Parameters: (ave,alpha) ave not used
       Call Gumbel_rnd(n_rnd, RandomPara(2), rnd_vec)
    Case Default
       ! Character(80) MESSAGE, Call SHUTDOWN(MESSAGE)
       Call SHUTDOWN('ieva.f90: Subroutine RandomNumbers')
       ! Stop
    End Select
  End Subroutine RandomNumbers
  ! ============================================================
  ! RandomNumbers ENDS HERE.
  ! ============================================================

  Subroutine UserRandomNumbers(n_rnd, mode, npts, x_user, f, rnd_vec)
    ! INPUTS:
    !   n_rnd: Number of random numbers returned
    !   mode:  1) density given
    !          2) distribution given
    !          3) discrete distribution (as histogram)
    !   npts:  number of points in the input function
    !   x_user:    the points on the x-axis
    !   f: the density or distribution values
    ! OUTPUTS:
    !   RandomVector(n_rnd): The random sample
    Implicit None

    ! Passed Variables
    Integer, Intent(IN) :: n_rnd, mode, npts
    Real(EB), Intent(IN) :: x_user(npts), f(npts)
    Real(EB), Intent(OUT) :: rnd_vec(n_rnd)
    ! rnd_vec:  32 bit random numbers as a vector
    ! n_rnd:   vector lenght

    ! Local variables
    Integer i
    Real(EB) p, x
    Real(EB), Allocatable :: yi(:), invdist(:)
    !    External Interpolate1d, Dens2InvDistr

    Call Random_number(rnd_vec)
    !    Call ranmar(rnd_vec,n_rnd)        ! Uniform (0,1) on the y axis

    ! Calculate next the inverse cdf
    Allocate( yi(npts),STAT=IZERO )
    Call ChkMemErr('STAT','UserRandomNumbers yi',IZERO)
    Allocate( invdist(npts),STAT=IZERO )
    Call ChkMemErr('STAT','UserRandomNumbers invdist',IZERO)

    Call Dens2InvDistr(mode, npts, x_user, f, yi, invdist)

    Do i = 1, n_rnd
       p = rnd_vec(i)
       If ( mode == 3 ) Then   ! Discrete density
          ! Next if statement is needed if the first point of the given
          ! density is not equal to zero (interpolation routine can not
          ! interpolate below the first value...)
          If ( p < yi(1) ) Then
             rnd_vec(i) = x_user(1)
          Else
             Call Interpolate1d(npts, yi, invdist, p, x)
             ! Now a discrete distribution
             rnd_vec(i) = x_user( Max( 1,Min( npts,Int(1.0_EB + x) ) ) )
          End If
       Else
          Call Interpolate1d(npts, yi, invdist, p, x)
          rnd_vec(i) = x
       End If
    End Do

    Deallocate( invdist)
    Deallocate( yi)

  End Subroutine UserRandomNumbers

  !
  Subroutine Beta_rnd(n_rnd,a,b,rnd_vec)
    Implicit None
    ! Passed Variables
    Integer n_rnd
    Real(EB) a, b
    Real(EB) rnd_vec(n_rnd)
    ! a: first parameter, b: second parameter, n_rnd: vector lenght
    Integer i, imode, status
    Real(EB) p,q,x,y, bound
    !
    Call Random_number(rnd_vec)
    imode = 2   ! calculate inverse cdf
    Do i = 1, n_rnd
       p = rnd_vec(i)
       q = 1.0_EB-p
       Call cdfbet(imode,p,q,x,y,a,b,status,bound)
       rnd_vec(i) = x
    End Do
  End Subroutine Beta_rnd

  Subroutine Gamma_rnd(n_rnd,a,b,rnd_vec)
    Implicit None
    ! Passed Variables
    Integer n_rnd
    Real(EB) a, b
    Real(EB) rnd_vec(n_rnd)
    ! a: first parameter, b: second parameter, n_rnd: vector lenght
    Integer i, imode, status
    Real(EB) p,q,x, bound, binv
    !
    Call Random_number(rnd_vec)
    imode = 2   ! calculate inverse cdf
    Do i = 1, n_rnd
       p = rnd_vec(i)
       q = 1.0_EB-p
       binv = 1.0_EB/b
       Call cdfgam(imode,p,q,x,a,binv,status,bound)
       rnd_vec(i) = x
    End Do
  End Subroutine Gamma_rnd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Lognormal distribution
  ! If X has a lognormal distribution, then log(X) is normally distributed.
  ! Here the logarithm is the natural logarithm, that is to base e, sometimes
  ! denoted as ln.  To generate random variates from this distribution, generate
  ! a random deviate from the normal distribution with mean and variance equal
  ! to the mean and variance of the logarithms of X, then take its exponential.

  ! Relationship between the mean & variance of log(X) and the mean & variance
  ! of X, when X has a lognormal distribution.
  ! Let m = mean of log(X), and s^2 = variance of log(X)
  ! Then
  ! mean of X     = exp(m + 0.5s^2)
  ! variance of X = (mean(X))^2.[exp(s^2) - 1]

  ! In the reverse direction (rarely used)
  ! variance of log(X) = log[1 + var(X)/(mean(X))^2]
  ! mean of log(X)     = log(mean(X) - 0.5var(log(X))

  ! N.B. The above formulae relate to population parameters; they will only be
  !      approximate if applied to sample values.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine LogNormal_rnd(n_rnd,a,b,cutoff,shift,rnd_vec)
    Implicit None
    ! Passed Variables
    Integer n_rnd
    Real(EB) a, b, cutoff, shift
    Real(EB) rnd_vec(n_rnd)
    ! rnd_vec:  32 bit random numbers as a vector
    ! a: first parameter, b: second parameter, n_rnd: vector lenght
    Integer i, imode, status
    Real(EB) p,q,x,bound, pmax, pmin, xmax
    !
    Call Random_number(rnd_vec)
    imode = 1   ! calculate cdf
    xmax = Log(cutoff-shift) ! Natural logarithm, i.e., base is 'e'
    Call cdfnor(imode,p,q,xmax,a,b,status,bound)
    pmax = p
    pmin = 0.0_EB
    imode = 2   ! calculate inverse cdf
    Do i = 1, n_rnd
       p = rnd_vec(i)*(pmax-pmin) + pmin
       q = 1.0_EB-p
       Call cdfnor(imode,p,q,x,a,b,status,bound)
       x = Min(x,xmax)
       rnd_vec(i) = Exp(x)+shift    ! lognormal
    End Do
  End Subroutine LogNormal_rnd

  Subroutine TNormal_rnd(n_rnd,a,b,xmin,xmax,rnd_vec)
    ! Truncated normal distribution
    Implicit None
    ! Passed Variables
    Integer n_rnd
    Real(EB) a, b, xmin, xmax
    Real(EB) rnd_vec(n_rnd)
    ! rnd_vec:  32 bit random numbers as a vector
    ! a: first parameter, b: second parameter, n_rnd: vector lenght
    Integer i, imode, status
    Real(EB) p,q,x,bound, pmax, pmin
    !
    Call Random_number(rnd_vec)
    imode = 1   ! calculate cdf
    x = xmin
    Call cdfnor(imode,p,q,x,a,b,status,bound)
    pmin = p
    x = xmax
    Call cdfnor(imode,p,q,x,a,b,status,bound)
    pmax = p
    imode = 2   ! calculate inverse cdf
    Do i = 1, n_rnd
       p = rnd_vec(i)*(pmax-pmin) + pmin
       q = 1.0_EB-p
       Call cdfnor(imode,p,q,x,a,b,status,bound)
       rnd_vec(i) = Min(Max(x,xmin),xmax)
       ! rnd_vec(i) = Min(x,xmax)
    End Do
  End Subroutine TNormal_rnd

  Subroutine Normal_rnd(n_rnd,a,b,rnd_vec)
    Implicit None
    ! Passed Variables
    Integer n_rnd
    Real(EB) a, b
    Real(EB) rnd_vec(n_rnd)
    ! rnd_vec:  32 bit random numbers as a vector
    ! a: first parameter, b: second parameter, n_rnd: vector lenght
    Integer i, imode, status
    Real(EB) p,q,x,bound
    !
    Call Random_number(rnd_vec)
    imode = 2   ! calculate inverse cdf
    Do i = 1, n_rnd
       p = rnd_vec(i)
       q = 1.0_EB-p
       Call cdfnor(imode,p,q,x,a,b,status,bound)
       rnd_vec(i) = x
    End Do
  End Subroutine Normal_rnd

  Subroutine Triang_rnd(n_rnd,peak,a,b,rnd_vec)
    Implicit None
    ! Passed Variables
    Integer n_rnd
    Real(EB) a, b, peak
    Real(EB) rnd_vec(n_rnd)
    ! rnd_vec:  32 bit random numbers as a vector
    ! a: first parameter, b: second parameter, n_rnd: vector lenght
    Integer i, imode
    Real(EB) p,x
    !
    Call Random_number(rnd_vec)
    imode = 2   ! calculate inverse cdf
    Do i = 1, n_rnd
       p = rnd_vec(i)
       Call cdftri(imode,p,x,peak,a,b)
       rnd_vec(i) = x
    End Do
  End Subroutine Triang_rnd

  Subroutine cdftri(imode,p,x,peak,a,b)
    ! Cumulative (and inverse) triangular distribution
    ! Written by T.K. (25/2/2003), works correctly
    Implicit None
    Integer imode
    Real(EB) p, x, peak, a, b
    ! Local variables
    Real(EB) y_peak, p_peak, xkk1, xkk2, y2
    y_peak   = 2.0_EB / (b-a)  ! maximum of the density
    p_peak = (peak-a)/(b-a)  ! cdf at the 'peak'
    xkk1 = y_peak/(peak-a)   ! derivate 1
    xkk2 = y_peak/(b-peak)  ! derivate 2
    If ( imode == 1 ) Then   ! cdf
       If ( x <= a ) Then
          p = 0.0_EB
       Else If ( x <= peak ) Then
          p = 0.5_EB*(x-a)*(x-a)*xkk1
       Else If ( x < b ) Then
          y2 = y_peak - (x-peak)*xkk2
          p = p_peak + (x-peak)*0.5_EB*(y_peak+y2)
       Else
          p = 1.0_EB
       End If
    Else If (imode == 2 ) Then   ! inverse cdf
       If ( p <= 0.0_EB ) Then
          x = a
       Else If ( p <= p_peak ) Then
          x = a + Sqrt( 2.0_EB*p / xkk1 )
       Else If ( p < 1.0_EB ) Then
          x = b - Sqrt( 2.0_EB*(1.0_EB-p)/xkk2)
       Else
          x = b
       End If
    Else
       Call SHUTDOWN('ieva.f90: Subroutine cdftri')
       ! Stop
    End If
  End Subroutine cdftri

  Subroutine Weibull_rnd(n_rnd,a,b,rnd_vec)
    Implicit None
    ! Passed Variables
    Integer n_rnd
    Real(EB) a, b
    Real(EB) rnd_vec(n_rnd)
    ! rnd_vec:  32 bit random numbers as a vector
    ! a: first parameter, b: second parameter, n_rnd: vector lenght
    Integer i, imode
    Real(EB) p,x
    !
    Call Random_number(rnd_vec)
    imode = 2   ! calculate inverse cdf
    Do i = 1, n_rnd
       p = rnd_vec(i)
       Call cdfwei(imode,p,x,a,b)
       rnd_vec(i) = x
    End Do
  End Subroutine Weibull_rnd

  Subroutine cdfwei(imode,p,x,alpha,lambda)
    ! Cumulative (and inverse) Weibull distribution
    ! Written by T.K. (25/2/2003)
    ! f(x) = alpha*lambda*(lamda*x)^(alpha-1)*exp(-(lambda*x)^alpha)
    Implicit None
    Integer imode
    Real(EB) p, x, alpha, lambda

    If ( (alpha <= 0.0_EB) .Or. (lambda <= 0.0_EB ) ) Then
       Call SHUTDOWN('ieva.f90: Subroutine cdfwei 1')
       ! Stop
    End If
    If ( imode == 1 ) Then   ! cdf
       If ( x <= 0.0_EB ) Then
          p = 0.0_EB  !  Weibull defined: x > 0
       Else
          p = 1.0_EB - Exp(-((lambda*x)**alpha))
       End If
    Else If (imode == 2 ) Then   ! inverse cdf
       If ( p <= 0.0_EB ) Then
          x = 0
       Else If ( p < 1.0_EB ) Then
          x = (1.0_EB/lambda)*((-Log(1.0_EB-p))**(1.0_EB/alpha))
       Else
          x = Huge(p)  ! largest positive number
       End If
    Else
       Call SHUTDOWN('ieva.f90: Subroutine cdfwei 2')
       ! Stop
    End If
  End Subroutine cdfwei

  Subroutine Gumbel_rnd(n_rnd,a,rnd_vec)
    Implicit None
    ! Passed Variables
    Integer n_rnd
    Real(EB) a
    Real(EB) rnd_vec(n_rnd)
    ! rnd_vec:  32 bit random numbers as a vector
    ! a: first parameter, b: second parameter, n_rnd: vector lenght
    Integer i, imode
    Real(EB) p,x
    !
    Call Random_number(rnd_vec)
    imode = 2   ! calculate inverse cdf
    Do i = 1, n_rnd
       p = rnd_vec(i)
       Call cdfgum(imode,p,x,a)
       rnd_vec(i) = x
    End Do
  End Subroutine Gumbel_rnd

  Subroutine cdfgum(imode,p,x,alpha)
    ! Cumulative (and inverse) Gumbel distribution
    ! Written by T.K. (25/2/2003)
    ! f(x) = alpha*exp(-alpha*x)*exp( -exp(-alpha*x) )
    Implicit None
    Integer imode
    Real(EB) p, x, alpha
    ! Local variables

    If ( alpha <= 0.0_EB ) Then
       Call SHUTDOWN('ieva.f90: Subroutine cdfgum 1')
       ! Stop
    End If
    If ( imode == 1 ) Then   ! cdf
       If ( x < 0.0_EB ) Then
          p = Exp( -1.0_EB/Exp(alpha*x) ) ! numerically better
          !          p = Exp( -1.0_EB*Exp(-alpha*x) )
       Else
          p = Exp( -1.0_EB*Exp(-alpha*x) )
       End If
    Else If (imode == 2 ) Then   ! inverse cdf
       If ( p <= 0.0_EB ) Then
          x = -Huge(p)  ! smallest possible value
       Else If ( p < 1.0_EB ) Then
          x = (-1.0_EB/alpha)*Log(-1.0_EB*Log(p))
       Else
          x = Huge(p)  ! largest positive number
       End If
    Else
       Call SHUTDOWN('ieva.f90: Subroutine cdfgum 2')
       ! Stop
    End If
  End Subroutine cdfgum

  !-------------------------------------------------------------
  Subroutine Uniform_rnd(n_rnd,a,b,rnd_vec)
    Implicit None
    ! Passed Variables
    Integer n_rnd
    Real(EB) a, b
    Real(EB) rnd_vec(n_rnd)
    ! rnd_vec:  32 bit random numbers as a vector
    ! a: lower bound, b: upper bound, n_rnd: vector lenght
    !
    Call Random_number(rnd_vec)
    rnd_vec(:) = rnd_vec(:)*(b-a) + a  ! Inverse CDF to get value on the x axis
  End Subroutine Uniform_rnd

  !-------------------------------------------------------------

  Subroutine Dens2InvDistr(mode, npts, x_user, f, yi, invdist)
    Implicit None
    !
    ! Calculates the cumulative distribution function from
    ! the probability density f(x). The inverse is also
    ! calculated.
    !
    ! Assumptions: f(x) >= 0 for all x
    !              f(x)  > 0 for some x
    ! These assumptions are not checked, because they should
    ! hold, if f(x) is a probability density function.
    !
    !----------------------------------------------
    ! Argument definitions
    !----------------------------------------------
    Integer mode     ! 1: density given, 2: distr. given
    Integer npts     ! number of points on the x-axis
    Real(EB) x_user(npts), f(npts), yi(npts), invdist(npts)
    !----------------------------------------------
    ! Local variables
    !----------------------------------------------
    Integer i
    Real(EB) sum, dx

    invdist = 0.0_EB

    If ( mode == 1 ) Then
       ! Calculate the distribution (trapezoidal rule)
       sum  = 0.0_EB
       Do i = 2, npts
          dx      = x_user(i) - x_user(i-1)
          sum     = sum + 0.5_EB*dx*( f(i) + f(i-1) )
          invdist(i) = sum         ! the distribution function
          ! invdist used as tmp array to save memory
       End Do
       invdist = (1.0_EB/sum)*invdist   ! now it is normalized
       invdist = Min(invdist,1.0_EB)    ! one or less (do not trust numerical accuracy
       ! of the normalization)
       ! Make the inverse of the distribution
       yi      = invdist
       invdist = x_user
    Else If ( mode == 3 ) Then
       ! Now a discrere distribution density is given:
       ! It is assumed that before the first point the density is zero and
       ! after the last given point the density is also zero
       sum  = 0.0_EB
       Do i = 1, npts
          sum     = sum +  f(i)
          invdist(i) = sum         ! the distribution function
          ! invdist used as tmp array to save memory
       End Do
       invdist = (1.0_EB/sum)*invdist   ! now it is normalized
       invdist = Min(invdist,1.0_EB)    ! one or less (do not trust numerical accuracy
       ! of the normalization)
       ! Make the inverse of the distribution
       yi      = invdist
       Do i = 1, npts
          invdist(i) = Dble(i)
       End Do
    Else
       ! Make the inverse of the distribution
       yi      = f
       invdist = x_user
    End If

  End Subroutine Dens2InvDistr

  Subroutine Interpolate1d(nx,x,y,x_user,ans)
    Integer nx
    Real(EB)    x(nx), y(nx), x_user, ans
    Integer jl, jm, ju
    !
    ! Find the index of the value just below X_USER
    !
    jl = 0
    ju = nx+1
    Do While ( (ju-jl) > 1 )
       jm = ( ju + jl )/2
       If ( (x(nx) > x(1) ) .Eqv. ( x_user > x(jm) ) ) Then
          jl = jm
       Else
          ju = jm
       End If
    End Do
    !
    ! Interpolate between JL and JL+1
    !
    If ( jl >= nx ) Then
       ans = y(nx)
    Else If ( jl < 1 ) Then
       ans = y(1)
    Else
       ans = y(jl) + (x_user-x(jl)) * (y(jl+1)-y(jl)) &
            / (x(jl+1)-x(jl)+1.0e-30_EB*(x(2)-x(1)))
    End If
  End Subroutine Interpolate1d

End Module STAT
