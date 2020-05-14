!  *********************************************************************
!  *                                                                   *
!  *                       real*8 function derf                        *
!  *                                                                   *
!  *********************************************************************
!  Double Precision Version 1.01
!  Written by Gordon A. Fenton, TUNS, Mar  7, 1999
!  Latest Update: Mar  7, 1999
!
!  PURPOSE  returns the error function
!
!  DESCRIPTION
!  This function returns the error function defined by
!
!                       x
!     erf(x) = (2/pi) int exp( -z^2 ) dz
!                       0
!
!  It does so by employing the cumulative standard normal distribution.
!  See dphi.f
!
!  ARGUMENTS
!
!	x	real value giving the upper bound on the integrand. If x
!		is negative, -erf(-x) is returned. (input)
!
!  REVISION HISTORY:
!
!-------------------------------------------------------------------------
      real*8 function derf(x)
      implicit real*8 (a-h,o-z)
      data rttwo/1.4142135623730950488d0/
      data zero/0.d0/, one/1.d0/, two/2.d0/

      derf = two*dphi(rttwo*x) - one

      return
      end
