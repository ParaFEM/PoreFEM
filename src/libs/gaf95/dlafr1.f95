!  *******************************************************************
!  *                                                                 *
!  *                     Real*8 Function dlafr1                      *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 1.41
!  Written by Gordon A. Fenton, Sept. 1989.
!  Latest Update: May 9, 2001
!
!  PURPOSE  returns the covariance between two points in a 1-D fractional
!           Gaussian noise process.
!
!  This function computes and returns the covariance between two points,
!  separated by the distance T, in a 1-D fractional Gaussian noise
!  (fGn) process. The continuous covariance function of the fGn process
!  is approximated by (for lag T)
!
!              var              2H       2H           2H
!     B(T) = -------- [ |T + pb|  -  2|T|  +  |T - pb|  ]		(1)
!                 2H
!             2*pb
!
!  where var is the point variance, H is the self-similarity parameter
!  as defined by Mandelbrot, also called the Hurst exponent, and pb is the
!  length over which the fractional Brownian noise is averaged in order to
!  define the derivative process which is this fGn. Typically pb should be
!  taken as equal to the minimum distance between field points in the
!  simulated process. var, H, and pb are brought into this routine
!  through the common block DPARAM (other parameters in this common are
!  ignored).
!
!  If var < 0, then this function returns the variance of a local average
!  of the process, |var|*V(T), averaged over the distance T.
!  The variance function for the fGn process is given by
!
!                        r        r      r            r
!                |T + pb| - 2( |T| + |pb| ) + |T - pb|
!       V(T) = -----------------------------------------		(4)
!                      2                  2H
!                     T * (2H+1)*(2H+2)*pb
!
!  where r = (2*H+2) and T is the averaging length. Notice that as T goes to
!  zero the above becomes 0/0. However, the limiting value is 1.
!  When var < 0, this function returns the variance of the local average,
!  that is it returns |var|*V(T).
!
!  REVISION HISTORY:
!  1.1	used the same notation in all the fGn process variance functions
!	(Mar 15/95)
!  1.2	added lvarfn flag to end of /dparam/ common block. Set it to
!	true if the variance function is to be returned. Otherwise this
!	function returns the covariance. (Apr 29/00)
!  1.3	eliminated lvarfn -- now return covariances only if var < 0 (Mar 2/01)
!  1.4	reversed default - now return covariances if var > 0 (Apr 11/01)
!  1.41	revised above documentation to reflect revision 1.4 (May 9/01)
! =======================================================================
      real*8 function dlafr1( X )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, H, da, db
      data zero/0.d0/, one/1.d0/, two/2.d0/
      abs(x) = dabs(x)

      e0 = two*H
      if( var .lt. zero ) then			! return variance function
         if( X .eq. zero ) then
            dlafr1 = -var
         else
            e1 = e0 + one
            e2 = e0 + two
            t1 = abs(X+pb)**e2-two*(abs(X)**e2+pb**e2)+abs(X-pb)**e2
            dlafr1 = -var*t1/(X*X*e1*e2*pb**e0)
         endif
      else					! return covariance
           t1 = abs(X+pb)**e0 - two*(abs(X)**e0) + abs(X-pb)**e0
           dlafr1 = var*t1/(two*(pb**e0))
      endif

      return
      end
