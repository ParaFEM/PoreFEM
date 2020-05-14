!   ********************************************************************
!   *                                                                  *
!   *                       function ntwom                             *
!   *                                                                  *
!   ********************************************************************
!   Integer Version 1.0
!   Written by G. A. Fenton, Princeton, June 7, 1988.
!
!   PURPOSE  return M such that a given N = 2**M
!
!   Returns the power that two must be raised to to get the value of the
!   argument, or in other words, returns m such that
!
!        N = 2**m
!
!   where N is the input argument integer and N > 0. If N is not a power
!   of 2, then the next smaller power is returned. The maximum power
!   returned is 40.
!-------------------------------------------------------------------------
      integer function ntwom( N )
      integer N, k

      k = 1
      do 20 ntwom = 0, 40
         if( N .le. k ) go to 30
         k = 2*k
  20  continue

  30  if( N .lt. k ) ntwom = ntwom - 1

      return
      end
