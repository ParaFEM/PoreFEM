! *******************************************************************
!  *                                                                 *
!  *                           function gamma                        *
!  *                                                                 *
!  *******************************************************************
!  Single Precision Version 1.0
!  Written by Gordon A. Fenton, TUNS, 1991
!
!  PURPOSE  to return the Gamma (factorial) function
!
!  This routine calculates the definite integral
!
!                  inf
!    gamma(z) = int  t**(z-1)*exp(-t) dt
!                  0
!
!  which is the Gamma function, for argument z. The method used
!  is a series expansion followed by a recursion relationship. The
!  coefficients used in the series expansion are obtained from Abramowitz
!  and Stegun ("Handbook of Mathematical Functions", US Dept. of Commerce,
!  10th printing, 1972), pg 256. The result should have a relative accuracy
!  of better than 1.e-14. If gamma(zz) would result in a number bigger
!  than the maximum single precision number, "big", then "big" is
!  returned. This occurs also if zz is a negative integer.
!  Arguments are as follows;
!
!     zz    real value for which the gamma function is desired. (input)
!
!  NOTES: 1) "eps" is the machine epsilon (smallest number which can be added
!            to 1.0 without being lost)
!         2) "big" is the biggest floating point number
!         3) set IDEBUG = 1 if error messages should be reported
!--------------------------------------------------------------------------
      real function gammaf( zz )
      parameter (IDEBUG = 1)
      dimension c(26)
      data c/1.0000000000000000,  0.5772156649015329,     &
            -0.6558780715202538, -0.0420026350340952,     &
             0.1665386113822915, -0.0421977345555443,     &
            -0.0096219715278770,  0.0072189432466630,     &
            -0.0011651675918591, -0.0002152416741149,     &
             0.0001280502823882, -0.0000201348547807,     &
            -0.0000012504934821,  0.0000011330272320,     &
            -0.0000002056338417,  0.0000000061160950,     &
             0.0000000050020075, -0.0000000011812746,     &
             0.0000000001043427,  0.0000000000077823,     &
            -0.0000000000036968,  0.0000000000005100,     &
            -0.0000000000000206, -0.0000000000000054,     &
             0.0000000000000014,  0.0000000000000001/
      data zero/0./, eps/0.119209e-06/, one/1./, big/1.701411733e+38/
      data pi/3.141592653589793238462643/

   1  format(a,e13.6,a)
!----------------------------- start executable statements
      z = abs(zz)

      q = z + eps
      n = int( q )
      f = z - float(n)
!					is z an integer?
      if( f .lt. q*eps ) then
         if( zz .le. zero ) go to 20
         n = n - 1
         f = one
         a = one
      else
!					series expansion for f < 1.0
         a = f*(c(22) + f*(c(23) + f*(c(24) + f*(c(25) + f*c(26)))))
         a = f*(c(17) + f*(c(18) + f*(c(19) + f*(c(20) + f*(c(21)+a)))))
         a = f*(c(12) + f*(c(13) + f*(c(14) + f*(c(15) + f*(c(16)+a)))))
         a = f*(c( 7) + f*(c( 8) + f*(c( 9) + f*(c(10) + f*(c(11)+a)))))
         a = f*(c( 2) + f*(c( 3) + f*(c( 4) + f*(c( 5) + f*(c( 6)+a)))))
         a = one/(f*(c( 1) + a))
      endif
!					apply recursion and watch for overflow
      if( n .gt. 0 ) then
         a = a*f
         do 10 i = 1, n-1
            q = f + float(i)
            if( a .gt. (big/q) ) go to 20
            a = a*q
  10     continue
      endif

      if( zz .lt. zero ) then
         gammaf = -pi/(z*a*sin(pi*z))
      else
         gammaf = a
      endif

      return
!				overflow error
  20  if( IDEBUG .eq. 1 )                                             &
         write(0,1)'Overflow calculating Gamma(',zz,') in GAMMA.'
      gammaf = big
      return
      end
