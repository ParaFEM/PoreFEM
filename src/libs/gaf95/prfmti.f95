!  *********************************************************************
!  *                                                                   *
!  *                         subroutine prfmti                         *
!  *                                                                   *
!  *********************************************************************
!  Integer Version 1.1
!  Written by Gordon A. Fenton, TUNS, Sun Jun  8 21:34:21 1997
!  Latest Update: Oct 2, 2001
!
!  PURPOSE  prints an integer to unit k according to an optional field
!           width
!
!  This routine writes the integer `ival' to an output file
!  using an optional field width specification.
!
!  Arguments to this routine are as follows;
!
!    ival	integer containing the number to print. (input)
!
!      iw	integer containing the desired field width. If iw < 0
!		then the minimum field width is computed internally (note
!		that integers larger than 256 digits cannot be printed).
!		If iw is larger than required, the number is right
!		justified. (input)
!
!       k	output unit number to which the number is printed. (input)
!
!  REVISION HISTORY:
!  1.1	calling routine now provides the iw format (Oct 2/01)
!-------------------------------------------------------------------------
      subroutine prfmti(ival,iw,k)
      character fstr*256, d(10)*1
!					basic digits
      data d/'0','1','2','3','4','5','6','7','8','9'/

   1  format(a,$)
!					get absolute value of val
      iaval = iabs(ival)

!					build character string of number
      m = 0
      n = 257				! work backwards in fstr
  20  m = m + 1
      n = n - 1
      i = iaval/10
      ii = 10*i
      j = iaval - ii
      fstr(n:n) = d(j+1)
      iaval = i
      if( iaval .gt. 0 ) go to 20
!					add `-' if negative
      if( ival .lt. 0 ) then
         m = m + 1
         n = n - 1
         fstr(n:n) = '-'
      endif
!					now print it
      do 30 i = m+1, iw
         write(k,1) ' '
  30  continue
      write(k,1) fstr(n:256)

      return
      end
