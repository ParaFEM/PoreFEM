!  *******************************************************************
!  *                                                                 *
!  *                     Function dlafs3                             *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 1.4
!  Written by Gordon A. Fenton, TUNS, Mar. 15, 1995
!  Latest Update: Jan 11, 2005
!
!  PURPOSE  returns the covariance between two points in a partially
!           separable 3-D fractional Gaussian noise field.
!
!  Returns the covariance between points in a 3-D fractional Gaussian noise
!  process which is isotropic in the X-Y plane. The covariance function is
!  partially separable, and is given by
!
!         B(X,Y,Z) = var*r(X,Y)*r(Z)
!
!  where r(Z) is the correlation function in the z-direction,
!
!               1               2G       2G           2G
!     r(Z) = -------- [ |Z + pb|  -  2|Z|  +  |Z - pb|  ]		(1)
!                 2G
!             2*pb
!
!  while, in the X-Y plane, the isotropic correlation function r(X,Y) is
!  given by
!
!                     1               2H       2H           2H
!  r(X,Y) = r(s) = -------- [ |s + pb|  -  2|s|  +  |s - pb|  ]		(2)
!                       2H
!                   2*pb
!
!  where s = sqrt(X*X + Y*Y), and `pb' is the length over which the fractional
!  Brownian motion is averaged in order to make this, the derivative process,
!  exist. Normally `pb' is selected to be quite small (of the order of the
!  size of the discretization interval). The parameters `G' and `H' are the
!  self-similar parameters, or Hurst exponents, (0.5 < G,H < 1) with
!  G = H = 1 giving a total correlated field and G = H = 0.5 giving a white
!  noise field (for arbitrarily small pb).
!
!  The parameters var, pb, H, and G are brought in through the common block
!  /dparam/. X, Y, and Z are the function arguments and are the elements
!  of the lag vector separating the two points in the field for which the
!  covariance is desired.
!
!  If var < 0, then this function computes the variance of a local average
!  of the random field, averaged over a domain having side dimensions {X,Y,Z}.
!  In this case, the variance function (which is also separable)
!
!     V(X,Y,Z) = V(X,Y)*V(Z)
!
!  is used, and the function returns |var|*V(X,Y,Z). The 1-D
!  variance function corresponding to Eq. (1) is given by
!
!                        r        r      r            r
!                |Z + pb| - 2( |Z| + |pb| ) + |Z - pb|
!       V(Z) = -----------------------------------------		(3)
!                      2                  2H
!                     Z * (2H+1)*(2H+2)*pb
!
!  where r = (2*G+2). Notice that as Z goes to zero the above becomes 0/0, but
!  its limiting value is 1.0.
!
!  The variance function which yields a close approximation to Eq. (2)c  is given by
!
!                      r          r        r              r
!              |D + pb|  -  2( |D|  +  |pb|  )  + |D - pb|
!   V(X,Y) = --------------------------------------------------		(4)
!                          2                  2G
!                         D * (2G+1)*(2G+2)*pb
!
!  where r = (2H+2) and D = |X| + |Y| (the 1-norm of {X,Y}). The actual
!  covariance function corresponding to (4) is equal to (2) along the
!  X,Y coordinate axes, and remains close to (2) for G > 0.8. For G < 0.7,
!  the covariance in the diagonal direction is too high, but since the
!  covariance in such fields dies off very rapidly, the error is
!  negligible.
!  
!  Parameters to this process are brought in through the common block
!  /dparam/ and are described as follows;
!
!     var	the point variance of the fGn process.
!
!      pb	the averaging length discussed above.
!
!       H	the Hurst exponent governing the process in the X,Y plane.
!		In this implementation, 0.5 < H < 1; values of H near 1 yield
!		a correlation structure which remains very high (and thus a
!		covariance matrix which may be nearly singular). Values of H
!		near 0.5 yield a band- limited white noise process.
!
!       G	the Hurst exponent governing the process in the Z direction.
!		In this implementation, 0.5 < G < 1; values of G near 1 yield
!		a correlation structure which remains very high (and thus a
!		covariance matrix which may be nearly singular). Values of G
!		near 0.5 yield a band- limited white noise process.
!
!      da	dummy placeholder provided so that the common block /dparam/
!		has the same form for a variety of variance functions.
!
!  Arguments to this function are just the X, Y, and Z dimensions of the
!  lag distance (or the physical averaging region).
!
!  REVISION HISTORY:
!  1.1	added lvarfn flag to end of /dparam/ common block. Set it to
!	true if the variance function is to be returned. Otherwise this
!	function returns the covariance. (Apr 29/00)
!  1.11	minor changes to documentation above (Jul 14/00)
!  1.2	eliminated lvarfn -- now return covariances only if var < 0 (Apr 9/01)
!  1.3	reversed default - now return covariances if var > 0 (Apr 11/01)
!  1.31	revised above documentation to reflect revision 1.3 (May 9/01)
!  1.4	made X-Y the isotropic plane (Z is usually vertical) (Jan 11/05)
!---------------------------------------------------------------------- 
      real*8 function dlafs3( X, Y, Z )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, H, G, da
      data zero/0.d0/, one/1.d0/, two/2.d0/
      abs(q)  = dabs(q)
      sqrt(q) = dsqrt(q)

      if( var .lt. zero ) then			! return variance function
!						Z-direction variance function
         if( Z .eq. zero ) then
            v1 = one
         else
            e0 = two*G
            e1 = e0 + one
            e2 = e0 + two
            t1 = abs(Z+pb)**e2 - two*(abs(Z)**e2+pb**e2) + abs(Z-pb)**e2
            v1 = t1/(Z*Z*e1*e2*(pb**e0))
         endif
!						X-Y plane variance function
         if( (X .eq. zero) .and. (Y .eq. zero) ) then
            v2 = one
         else
            D  = abs(Y) + abs(X)
            e0 = two*H
            e1 = e0 + one
            e2 = e0 + two
            t1 = (D+pb)**e2 - two*(D**e2 + pb**e2) + abs(D-pb)**e2
            v2 = t1/(D*D*e1*e2*pb**e0)
         endif
         dlafs3 = -var*v1*v2
      else					! return covariance
!						Z-direction covariance
         if( Z .eq. zero ) then
            v1 = one
         else
            e0 = two*G
            t1 = abs(Z+pb)**e0 - two*(abs(Z)**e0) + abs(Z-pb)**e0
            v1 = t1/(two*(pb**e0))
         endif
!						Y-direction covariance
         if( (Y .eq. zero) .and. (X .eq. zero) ) then
            v2 = one
         else
            e0 = two*H
            D  = sqrt(Y*Y + X*X)
            t1 = (D+pb)**e0 - two*(D**e0) + abs(D-pb)**e0
            v2 = t1/(two*(pb**e0))
         endif
         dlafs3 = var*v1*v2
      endif

      return
      end
