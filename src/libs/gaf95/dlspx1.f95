!  *******************************************************************
!  *                                                                 *
!  *                     Function dlspx1                             *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 1.41
!  Written by Gordon A. Fenton, Mar. 25, 1994
!  Latest Update: May 9, 2001
!
!  PURPOSE  returns the covariance between two points along a 1-D random
!           process having an exponential (Gaussian) type covariance function
!
!  This function returns the covariance between two points, separated by
!  distance T, along a random process having Gaussian type covariance function
!
!      B(T) = var * exp{ -pi*[(|T|/dthx)^2] }
!
!
!  If var < 0, then this function returns the variance of a local average
!  of the process, |var|*V(T), averaged over the distance T, where the
!  variance function is given by
!
!         V(T) = (1/a*a) * [pi*(|T|/dthx)*erf(a) + exp(-a*a) - 1]
!
!  where a = (|T|*sqrt(pi))/dthx, erf(.) is the error function,
!  dthx is the directional scale of fluctuation, and var is the point
!  variance. The variance of a local average is given by var*V(X).
!
!  See page 186 of "Random Fields" by E. Vanmarcke and G.A. Fenton's thesis
!  for more details.
!
!  The parameters var and dthx are brought in via the common block
!  /dparam/.
!
!  The argument to this routine, T, is the distance between the two points
!  for which the covariance is desired (or the distance over which the local
!  averaging is performed, if var < 0).
!
!  REVISION HISTORY:
!  1.1	fixed argument list and dthx = 0 case, improved limit handling
!	(Jun 19/96)
!  1.2	added lvarfn flag to end of /dparam/ common block. Set it to
!	true if the variance function is to be returned. Otherwise this
!	function returns the covariance. (May 1/00)
!  1.3	eliminated lvarfn -- now return covariances only if var < 0 (Mar 2/01)
!  1.4	reversed default - now return covariances if var > 0. (Apr 11/01)
!  1.41	revised above documentation to reflect revision 1.4 (May 9/01)
!---------------------------------------------------------------------------
      real*8 function dlspx1( T )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/, one/1.d0/, six/6.d0/
      data rtpi/1.77245385090551588d0/, pt1/0.1d0/
      data   pi/3.1415926535897932384d0/, big/1.d6/
      exp(s)   = dexp(s)
      abs(s)   = dabs(s)
      r_erf(s) = derf(s)

      if( var .lt. zero ) then			! return variance function
         aT = abs(T)
         if( dthx .ge. big*aT ) then
            if( dthx .eq. zero ) then
               dlspx1 = -var
            else
               rx = aT/dthx
               dlspx1 = -var*(one - pi*rx*rx/six)
            endif
         elseif( dthx .le. pt1*aT ) then
            rx = dthx/aT
            dlspx1 = -var*rx*(one - rx/pi)
         else
            rx = aT/dthx
            a  = rtpi*rx
            aa = a*a
            dlspx1 = -var*(pi*rx*r_erf(a) + exp(-aa) - one)/aa
         endif
      else					! return covariance
         if( dthx .eq. zero ) then
            if( T .eq. zero ) then
               dlspx1 = var
            else
               dlspx1 = zero
            endif
         else
            at = T/dthx
            dlspx1 = var * exp(-pi*at*at)
         endif
      endif

      return
      end

