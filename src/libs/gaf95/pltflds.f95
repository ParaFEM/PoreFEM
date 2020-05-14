!  ********************************************************************
!  *                                                                  *
!  *                       subroutine pltfld                          *
!  *                                                                  *
!  ********************************************************************
!  Single Precision Version 1.4
!  Written by Gordon A. Fenton, TUNS, 1991
!  Latest Update: Jun 25, 2006
!
!  PURPOSE  to output an array of data to a file in DISPLAY format (for
!           display on a PostScript printer).
!
!  This routine takes an array (2-D) of data and dumps it to a file
!  having a format readable by DISPLAY. Arguments are as follows;
!
!    job    character string containing the title of the run. (input)
!
!   sub1    character string containing the subtitle of the run. (input)
!
!   sub2    character string containing the sub-subtitle of the run. (input)
!
!      Z    real array of size at least N1 x N2 containing the data to
!           be displayed graphically. (input)
!
!     iz    leading dimension of Z exactly as specified in the calling
!           routine. (input)
!
!     N1    column (1st index) dimension of Z. (input)
!
!     N2    row (2nd index) dimension of Z. (input)
!
!     D1    physical dimension of the data in the X (1st index) direction.
!           (input)
!
!     D2    physical dimension of the data in the Y (2nd index) direction.
!           (input)
!
!  title    character string containing the title of the display. (input)
!
!   ifld    unit number to which output is sent. (input)
!
! NOTE: since RFLOW has inverted the direction of Y, ie the second index
!       of Z increases in the downwards direction, this routine inverts
!       the data with respect to Y (ie Z(1,N2) is printed first so that
!       the DISPLAY'ed image is correct).
!
!  REVISION HISTORY:
!  1.3	removed extra '0' between number of points and field data (Feb 17/05)
!  1.4	added lrev flag to tell display to map low to black (Jun 25/06)
!-----------------------------------------------------------------------------
      subroutine pltflds( job, sub1, sub2, Z, iz, N1, N2, D1, D2,title, ifld )
      real*4 z
	  real*8 d1,d2
	  dimension Z(iz)
      character*(*) job, sub1, sub2, title

   1  format(a)
   2  format(2e13.5)
   3  format()
   4  format(2i8)
   5  format(e12.4)

      lj = lnblnk(job)
      l1 = lnblnk(sub1)
      l2 = lnblnk(sub2)
      lt = lnblnk(title)
      write(ifld,1) job(1:lj)
      write(ifld,1) title(1:lt)
      write(ifld,1) sub1(1:l1)
      write(ifld,1) sub2(1:l2)
      write(ifld,1) '2'
      write(ifld,1) 'x'
      write(ifld,1) 'y'
      write(ifld,1) 'random field'
      write(ifld,2) 0., D1
      write(ifld,2) 0., D2
      write(ifld,4) N1, N2
      write(ifld,5) (Z(i), i = 1, iz)

      return
      end
