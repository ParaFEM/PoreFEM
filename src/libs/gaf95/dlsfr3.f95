!  *******************************************************************
!  *                                                                 *
!  *                     Function dlsfr3                             *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 1.31
!  Written by Gordon A. Fenton, TUNS, Mar. 15, 1995
!  Latest Update: May 9, 2001
!
!  PURPOSE  returns the covariance between two points in a separable 3-D
!           fractional Gaussian noise field
!
!  Returns the covariance between two points in the field separated by the
!  lag vector {X,Y,Z}. The covariance function is a separable fractional
!  Gaussian noise, given by
!
!         B(X,Y,Z) = var*r(X)*r(Y)*r(Z)
!
!  for two points separated by distance {X,Y,Z}, where r(X), r(Y), and r(Z)
!  are the correlation functions given by
!
!               1               2H       2H           2H
!     r(X) = -------- [ |X + pb|  -  2|X|  +  |X - pb|  ]		(1)
!                 2H
!             2*pb
!
!
!               1               2G       2G           2G
!     r(Y) = -------- [ |Y + pb|  -  2|Y|  +  |Y - pb|  ]		(2)
!                 2G
!             2*pb
!
!
!               1               2F       2F           2F
!     r(Z) = -------- [ |Z + pb|  -  2|Z|  +  |Z - pb|  ]		(2)
!                 2F
!             2*pb
!
!  and where `pb' is the length over which the fractional Brownian motion is
!  averaged in order to make this, the derivative process, exist. Normally
!  `pb' is selected to be quite small (of the order of the size of the
!  discretization interval). The parameters `F', `G', and `H' are the
!  self-similar parameters, or Hurst exponents, (0.5 < F,G,H < 1) with
!  F = G = H = 1 giving a total correlated field and F = G = H = 0.5 giving
!  a white noise field (for sufficiently small pb). Finally var is the
!  point variance of the process.
!
!  The parameters var, pb, H, G, and F are brought in via the common
!  block /dparam/.
!
!  If var < 0, then this function computes the variance of a local average
!  of the random field, |var|*V(X,Y,Z), averaged over a domain having side
!  dimensions {X,Y,Z}. For a separable covariance function, the variance
!  function is also separable,
!
!     V(X,Y,Z) = V(X)*V(Y)*V(Z)
!
!  where (X,Y,Z) gives the dimensions of the averaging region. The 1-D
!  variance function corresponding to Eq. (1) is given by
!
!                        r        r      r            r
!                |X + pb| - 2( |X| + |pb| ) + |X - pb|
!       V(X) = -----------------------------------------		(3)
!                      2                  2H
!                     X * (2H+1)*(2H+2)*pb
!
!  where r = (2*H+2). Notice that as X goes to zero the above becomes 0/0, but
!  its limiting value is 1.0. A similar relationship holds for the variance
!  function corresponding to (2).
!
!  Parameters to this process are brought in through the common block
!  /dparam/ and are described as follows;
!
!     var	the point variance of the fGn process.
!
!      pb	the averaging length discussed above. This is assumed common
!		for all three directions.
!
!       H	the Hurst exponent governing the process in the X direction.
!		In this implementation, 0.5 < H < 1; values of H near 1 yield
!		a correlation structure which remains very high (and thus a
!		covariance matrix which may be nearly singular). Values of H
!		near 0.5 yield a band- limited white noise process.
!
!       G	the Hurst exponent governing the process in the Y direction.
!		In this implementation, 0.5 < G < 1; values of G near 1 yield
!		a correlation structure which remains very high (and thus a
!		covariance matrix which may be nearly singular). Values of G
!		near 0.5 yield a band- limited white noise process.
!
!       F	the Hurst exponent governing the process in the Z direction.
!		In this implementation, 0.5 < F < 1; values of F near 1 yield
!		a correlation structure which remains very high (and thus a
!		covariance matrix which may be nearly singular). Values of F
!		near 0.5 yield a band- limited white noise process.
!
!  Arguments to this function are just the X, Y, and Z elements of the
!  lag vector between points in the field for which the covariance is
!  desired (or the dimensions of the physical averaging region, if var < 0).
!
!  REVISION HISTORY:
!  1.1	added lvarfn flag to end of /dparam/ common block. Set it to
!	true if the variance function is to be returned. Otherwise this
!	function returns the covariance. (May 1/00)
!  1.2	eliminated lvarfn -- now return covariances only if var < 0 (Apr 9/01)
!  1.3	reversed default - now return covariances if var > 0. (Apr 11/01)
!  1.31	revised above documentation to reflect revision 1.3 (May 9/01)
!---------------------------------------------------------------------- 
      real*8 function dlsfr3( X, Y, Z )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, H, G, F
      data zero/0.d0/, one/1.d0/, two/2.d0/
      abs(q)  = dabs(q)

      if( var .lt. zero ) then			! return variance function
!							X-direction var func
         if( X .eq. zero ) then
            v1 = one
         else
            e0 = two*H
            e1 = e0 + one
            e2 = e0 + two
            t1 = abs(X+pb)**e2-two*(abs(X)**e2+pb**e2)+abs(X-pb)**e2
            v1 = t1/(X*X*e1*e2*pb**e0)
         endif
!							Y-direction var func
         if( Y .eq. zero ) then
            v2 = one
         else
            e0 = two*G
            e1 = e0 + one
            e2 = e0 + two
            t1 = abs(Y+pb)**e2-two*(abs(Y)**e2+pb**e2)+abs(Y-pb)**e2
            v2 = t1/(Y*Y*e1*e2*pb**e0)
         endif
!						Z-direction var func
         if( Z .eq. zero ) then
            v3 = one
         else
            e0 = two*F
            e1 = e0 + one
            e2 = e0 + two
            t1 = abs(Z+pb)**e2-two*(abs(Z)**e2+pb**e2)+abs(Z-pb)**e2
            v3 = t1/(Z*Z*e1*e2*pb**e0)
         endif
!							final variance function
         dlsfr3 = -var*v1*v2*v3
      else					! return covariance
         eH = two*H
         eG = two*G
         eF = two*F
         pH = two*(pb**eH)
         pG = two*(pb**eG)
         pF = two*(pb**eF)
         tH = abs(X+pb)**eH-two*(abs(X)**eH)+abs(X-pb)**eH
         tG = abs(Y+pb)**eG-two*(abs(Y)**eG)+abs(Y-pb)**eG
         tF = abs(Y+pb)**eF-two*(abs(Y)**eF)+abs(Y-pb)**eF
         dlsfr2 = var*tH*tG*tF/(pH*pG*pF)
      endif

      return
      end
