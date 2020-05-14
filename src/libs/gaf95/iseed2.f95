!  **********************************************************************
!  *                                                                    *
!  *                   Integer Function iseed                           *
!  *                                                                    *
!  **********************************************************************
!  Single Precision Version 1.1
!  Written by Gordon A. Fenton, Princeton, Dec. 8, 1988.
!  Latest Update: Oct 14, 1996
!
!  PURPOSE  initializes the system pseudo-random number generated using
!           process ID as default seed.
!
!  Initializes the local random number generator RANDF. If the argument integer
!  seed (kseed) is zero, a random seed is calculated based on the process ID
!  of the parent process. The function returns the actual seed used. This
!  routine is system specific.
!
!  Notes:
!   1) any function which is reasonably certain to return a different integer
!      on each invocation of the calling process can be used to generate a
!      pseudo-random seed here. Some possibilities might be wall clock time,
!      Process ID number (UNIX based)...
!
!  Requires:
!   1) from Fortran lib: GETPID
!
!  REVISION HISTORY:
!  1.1	now using new randu function (RAN2 from Numerical Recipes, 2nd Ed)
!	(Oct 14/96)
!--------------------------------------------------------------------------
      integer function iseed2( kseed )
      integer kseed, getpid
!
      iseed2 = kseed
      if( kseed .eq. 0 ) then
!                                avoid seeds of 0
	  iseed2 = getpid() + 2
      endif
!                                initialize the generator
      rfirst2 = randu2( iseed2 )
      return
      end
