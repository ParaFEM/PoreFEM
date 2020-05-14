!  *******************************************************************
!  *                                                                 *
!  *                       function idamax                           *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 06/17/77 1.01
!  Written by Jack Dongarra, Linpack
!  Modified by Gordon A. Fenton, Aug. 24, 1993
!  Latest Update: Jun 9, 1999
!
!  PURPOSE   find the index of the maximum absolute element in a double
!            precision vector of length n
!
!
!  REVISION HISTORY:
!  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!------------------------------------------------------------------------
      integer function idamax( n, dx )
      real*8 dx(*), dmax
      integer i, n

      if( n .le. 0 ) then
         idamax = 0
         return
      endif
      idamax = 1
      if( n .eq. 1 )return

      dmax = dabs(dx(1))
      do 10 i = 2, n
         if( dabs(dx(i)) .gt. dmax ) then
            idamax = i
            dmax   = dabs(dx(i))
         endif
  10  continue

      return
      end
