!  *********************************************************************
!  *                                                                   *
!  *                          Function second                          *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 3.0
!  Written by Gordon A. Fenton, TUNS, 1990
!
!  PURPOSE returns elapsed user execution time in seconds.
!
!  Returns elapsed user execution time of the calling process in seconds.
!  This routine tends to be sysem specific and you may need to customize it
!  for your environment. This routine works for Sun's and VaX's running
!  Ultrix f77.
!---------------------------------------------------------------------------
      real function second()
      real tyme(2)

      second = etime(tyme)
!                               and for the HP-9000
!     second = 1.e-6*float(clock())

      return
      end
