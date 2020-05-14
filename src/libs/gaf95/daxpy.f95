!  *******************************************************************
!  *                                                                 *
!  *                       subroutine daxpy                          *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 03/11/78 1.01
!  Written by Jack Dongarra, Linpack.
!  Modified by Gordon A. Fenton, Aug. 24, 1993
!  Latest Update: Jun 9, 1999
!
!  PURPOSE   constant times a vector plus a vector.
!
!  This routine computes the output vector
!
!         dy(1)        = dy(1)        + da*dx(1)
!         dy(2)        = dy(2)        + da*dx(2)
!         dy(3)        = dy(3)        + da*dx(3)
!          .              .               .
!          .              .               .
!
!  Arguments are
!
!   n      the number of elements in the output vector
!
!   da     the constant
!
!   dx     the vector being multiplied by da
!
!   dy     on input, dy contains the additive vector. On output, dy contains
!          the result.
!
!  REVISION HISTORY:
!  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!----------------------------------------------------------------------------
      subroutine daxpy(n,da,dx,dy)
      real*8 dx(*), dy(*), da, zero
      integer i, n
      data zero/0.d0/

      if( n .le. 0 .or. da .eq. zero ) return

      do 10 i = 1, n
         dy(i) = dy(i) + da*dx(i)
  10  continue

      return
      end
