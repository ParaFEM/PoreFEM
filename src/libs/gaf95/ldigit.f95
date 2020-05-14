!  *********************************************************************
!  *                                                                   *
!  *                  logical function ldigit                          *
!  *                                                                   *
!  *********************************************************************
!  Logical Version 1.0
!  Written by Gordon A. Fenton, Dalhousie University, Oct  1, 2001
!  Latest Update: Oct  1, 2001
!
!  PURPOSE  returns true if the provided character is a digit in the range
!           0-9 or a decimal '.'
!
!  DESCRIPTION
!
!  This function checks to see if the provided character, c, is a digit in
!  the range 0 to 9 or a decimal point (.). It does so by checking the
!  integer index of the character, which assumes that the character set
!  follows the ASCII indexing.
!
!  ARGUMENTS
!
!  REVISION HISTORY:
!
!-------------------------------------------------------------------------
      logical function ldigit(c)
      character*1 c
      integer ic

      ic     = ichar(c)
      ldigit = ((ic .ge. 48) .and. (ic .le. 57)) .or. (ic .eq. 46)

      return
      end
