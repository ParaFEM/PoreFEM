!  *********************************************************************
!  *                                                                   *
!  *                    Integer Function Lnblnk                        *
!  *                                                                   *
!  *********************************************************************
!  Integer Version 1.1
!  Written by Gordon A. Fenton, Princeton, 1989.
!  Latest Update: May 16, 1996
!
!  PURPOSE  returns the index of the last non-blank character in a string.
!
!  Returns the index of the last non-blank character in the argument string.
!  Returns 0 if they are all blank.
!  This routine is for use if your version of Fortran does not have lnblnk.
!  It does, however, depend on the function `len' which returns the
!  dimensioned length of the string `str'.
!  Note the char(0) = ^@ character is the same as '\0' in C (ie ASCII 0).
!
!  REVISION HISTORY:
!  1.1	replaced ^@ character with explicit call to char(0)
!----------------------------------------------------------------------------
      integer function lnblnk( str )
      character*(*) str
      character*1 space, tab, null, char
      data space/' '/, tab/'	'/

      null = char(0)
      i = len( str )
      do 10 j = i, 1, -1
         if(       str(j:j) .ne. space                                    &
             .and. str(j:j) .ne. null                                     &
             .and. str(j:j) .ne. tab ) then
            lnblnk = j
            return
         endif
  10  continue
      lnblnk = 0

      return
      end
