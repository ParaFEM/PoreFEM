!  *******************************************************************
!  *                                                                 *
!  *                     Real*8 Function dlavx1                      *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 1.62
!  Written by Gordon A. Fenton, 1989
!  Latest Update: Jul 10, 2002
!
!  PURPOSE  returns the covariance between two points along a 1-D Markov
!           process.
!
!  Returns the covariance between two points, separated by distance T, in
!  a 1-D Markov process. The covariance function is given by
!
!       B(T) = var * exp( -2|T|/dthx )
!
!  where var is the point variance and dthx is the scale of fluctuation.
!  var and dthx are brought into this routine through the common
!  block DPARAM.
!
!  If var < 0, then this function returns the variance of a local average
!  of the process, |var|*V(T), averaged over the distance T, where V(T) is
!  the variance function (see E.H. Vanmarcke, "Random Fields: Analysis and
!  Synthesis", MIT Press, 1984, Chapter 5).
!                                  T
!       V(T) = (2/T) * INT  (1-s/T)*r(s) ds
!                                  0
!  where r(s) is the correlation function: r(s) = B(s)/var. The integration
!  gives us
!
!       V(T) = (2/a*a) * [ a + exp(-a) - 1]
!
!  where a = 2*T/dthx. When dthx >> T, a third order Taylor's series expansion
!  gives the numerically improved approximation
!
!	V(T) = [1 - 2*|T|/(3*dthx) ],	dthx >> T
!
!  and when dthx << T, the exp(-a) term disappears because `a' becomes very
!  large, leaving us with
!
!	V(T) = (dthx/T) * [ 1 - (dthx/2T) ]
!
!
!  REVISION HISTORY:
!  1.1	explicitly accomodate the case when dthx = 0.
!  1.2	corrected bug introduced in 1.1 (missed multiplying by var)
!  1.3	improved limit handling
!  1.4	added lvarfn flag to end of /dparam/ common block. Set it to
!	true if the variance function is to be returned. Otherwise this
!	function returns the covariance. (Apr 29/00)
!  1.5	eliminated lvarfn -- now return covariances only if var < 0 (Mar 2/01)
!  1.6	reversed default - now return covariances if var > 0 (Apr 11/01)
!  1.61	revised above documentation to reflect revision 1.6 (May 9/01)
!  1.62	corrected erroneous documentation above re variance func (Jul 10/02)
! =======================================================================
      real*8 function dlavx1( T )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/, half/0.5d0/, one/1.d0/, two/2.d0/, three/3.d0/
      data small/0.0067d0/, big/1.d6/
      abs(x) = dabs(x)
      exp(x) = dexp(x)

      aT = abs(T)
      if( var .lt. zero ) then			! return variance function
         if( dthx .ge. big*aT ) then
            if( dthx .eq. zero ) then
               dlavx1 = -var
            else
               dlavx1 = -var*(one - two*aT/(three*dthx))
            endif
         elseif( dthx .le. small*aT ) then
            rx = dthx/aT
            dlavx1 = -var*rx*(one - half*rx)
         else
            a = dthx/abs(T)
            b = half*a
            c = exp(-one/b) - one
            dlavx1 = -var*a*(one + b*c)
         endif
      else					! return covariance
         if( dthx .eq. zero ) then
            if( T .eq. zero ) then
               dlavx1 = var
            else
               dlavx1 = zero
            endif
         else
            dlavx1 = var*exp( -two*aT/dthx )
         endif
      endif

      return
      end
