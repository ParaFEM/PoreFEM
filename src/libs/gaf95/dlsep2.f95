!  *******************************************************************
!  *                                                                 *
!  *                     Function dlsep2                             *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 1.51
!  Written by Gordon A. Fenton, Sept. 1989
!  Latest Update: May 9, 2001
!
!  PURPOSE  returns the covariance between two points in a 2-D process having
!           a separable Markovian covariance function
!
!  This function returns the covariance between two points, separated by lag
!  vector {X,Y}, in a 2-D random field having separable Markovian covariance
!  function. The covariance function for this field is given by
!
!      B(X,Y) = var * exp{ -(2|X|/dthx) } * exp{ - (2|Y|/dthy) }
!
!  where var is the point variance, and dthx and dthy are the scales of
!  fluctuation in the x- and y-directions, respectively.
!  NOTE: this function is NOT isotropic even if dthx = dthy.
!  The parameters var, dthx, and dthy are brought in via the
!  common block /dparam/.
!
!  If var < 0, then this function returns the variance of a local average
!  of the process, |var|*V(X,Y), averaged over the domain X x Y.
!  For a separable covariance function, the variance function is also
!  separable and can be written as the product of two 1-D variance
!  functions corresponding to the exponential correlation function
!
!         V(X,Y) = V(X)*V(Y)
!
!  where individual functions are based on the 1-D analytical model
!
!                     dthx^2
!         V(X) = 2 * ------- [(2|X|/dthx) + exp{-2|X|/dthx} - 1]
!                     4*X^2
!
!  The variance of the local average is then given by var*V(X,Y). See page
!  186 of "Random Fields" by E. Vanmarcke and G.A. Fenton's thesis for more
!  details.
!
!  The arguments to this routine, X and Y, are just the components of the
!  lag vector separating the two points (or the size of the averaging
!  domain, if var < 0).
!
!  The data value `eund' is designed to avoid underflow errors arising from
!  exp(-a) for a > eund. We quite happily accept a zero in this case.
!
!  REVISION HISTORY:
!  1.2	added lvarfn flag to end of /dparam/ common block. Set it to
!	true if the variance function is to be returned. Otherwise this
!	function returns the covariance. (May 1/00)
!  1.3	eliminated lvarfn -- now return covariances only if var < 0 (Mar 27/01)
!  1.4	properly handled sign on var for covariances (Apr 5/01)
!  1.5	reversed default - now return covariances if var > 0. (Apr 11/01)
!  1.51	revised above documentation to reflect revision 1.5 (May 9/01)
!---------------------------------------------------------------------------
      real*8 function dlsep2( X, Y )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/, half/0.5d0/, one/1.d0/, two/2.d0/, four/4.d0/
      data eund/300.d0/
      exp(s) = dexp(s)
      abs(s) = dabs(s)

      if( var .lt. zero ) then			! return variance function
!							x-direction var func
         a = two*abs(X)/dthx
         if( X .eq. zero ) then
            dsepa = half
         elseif( a .lt. eund ) then
            dsepa = (a + exp(-a) - one)/(a*a)
         else
            dsepa = (a - one)/(a*a)
         endif
!							y-direction var func
         b = two*abs(Y)/dthy
         if( Y .eq. zero ) then
            dsepb = half
         elseif( b .lt. eund ) then
            dsepb = (b + exp(-b) - one)/(b*b)
         else
            dsepb = (b - one)/(b*b)
         endif

         dlsep2 = -four*var*dsepa*dsepb
      else					! var > 0, return covariance
         if( dthx .eq. zero ) then
            if( X .eq. zero ) then
               a = zero
            else
               dlsep2 = zero
               return
            endif
         else
            a = two*abs(X)/dthx
         endif
         if( dthy .eq. zero ) then
            if( Y .eq. zero ) then
               b = zero
            else
               dlsep2 = zero
               return
            endif
         else
            b = two*abs(Y)/dthy
         endif
         ab = a + b
         if( ab .ge. eund ) then
            dlsep2 = zero
         else
            dlsep2 = var*exp(-ab)
         endif
      endif

      return
      end
