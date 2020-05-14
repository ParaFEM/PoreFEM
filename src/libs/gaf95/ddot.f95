!  *******************************************************************
!  *                                                                 *
!  *                       subroutine daxpy                          *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 03/11/78 1.01
!  Written by Jack Dongarra, Linpack
!  Modified by Gordon A. Fenton, Aug. 1993
!  Latest Update: Jun 9, 1999
!
!  PURPOSE    forms the dot product of two vectors.
!
!
!  REVISION HISTORY:
!  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!-------------------------------------------------------------------------
      real*8 function ddot(n,dx,dy)
      real*8 dx(*),dy(*), zero
      integer i, n
      data zero/0.d0/

      ddot = zero
      if( n .le. 0 ) return

      do 10 i = 1, n
         ddot = ddot + dx(i)*dy(i)
  10  continue

      return
      end
