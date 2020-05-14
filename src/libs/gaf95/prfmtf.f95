!  *********************************************************************
!  *                                                                   *
!  *                         subroutine prfmtf                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.1
!  Written by Gordon A. Fenton, TUNS, Sun Jun  1 17:46:55 1997
!  Latest Update: Jul 1, 1997
!
!  PURPOSE  prints a real number using f format
!
!  This routine writes the real value `val' to an output file
!  using the format
!
!	l.r
!
!  where `l' expands into as many digits as necessary, unless otherwise
!  specified by the iw and id arguments, and `r' represents
!  the digits to the right of the decimal place. Since this is a single
!  precision subroutine, the maximum number of significant digits
!  displayed will be 7 and by default no more than 7 digits will appear
!  to the right of the decimal point. These defaults may be overridden
!  by specifying non-negative values of iw and id, as discussed below.
!
!  The floating point format will only be used if val is between 1.e+7
!  and 1.e-6 or zero. If val is outside this range, then a scientific format
!  is used (see prfmte), and the values of iw and id are ignored.
!
!  If both iw and id are negative, then the minimum width format required to
!  represent the number is determined internally.
!
!  If either iw or id are non-negative, then the format used to represent the
!  number is one of the following;
!
!	iw.id	(iw and id non-negative) In this case, iw is the minimum total
!		field width (including the decimal point) and id is the
!		maximum number of digits to show to the right of the decimal
!		place. For example, the format %7.2f implies that the number
!		12.345600 would be represented as '  12.35'. If iw is too
!		small to contain the number, it is increased accordingly.
!
!	iw	(id is negative) In this case, iw is the minimum field width
!		(including the decimal point) used to represent the number.
!		For example, the format %12f implies that the number 12.345600
!		would be represented as '     12.3456'. If iw is less than
!		the number of digits to the left of the decimal place, then
!		iw is increased accordingly (except that if more than 7 digits
!		are required, scientific notation is used).
!
!	.id	(iw is negative) In this case, id is the maximum number of
!		digits to show to the right of the decimal place. For example,
!		the format %.2f implies that the number 12.345600 would be
!		represented as '12.35'. id must be less than or equal to
!		seven.
!
!	.	(iw and id negative, but decimal point provided) In this case,
!		the decimal point just implies that the value zero is
!		represented as '0.', rather than just '0'. For example, the
!		format %.f would represent the number zero as '0.'.
!
!   Arguments to this routine are as follows;
! 
!     val	real value containing the number to print. (input)
!
!      iw	integer containing the minimum total field width of the number
!		to be printed. If iw is insufficient to contain the number, it
!		is increased accordingly. If iw < 0, then it is assumed to not
!		be set and a minimum value is computed internally. (input)
!
!      id	integer containing the maximum number of digits to show to the
!		right of the decimal point. If id < 0, then it is assumed to
!		not be set and the minimum number of digits to the right of
!		the decimal point required to show the number is computed
!		internally. (input)
!
!    ldot	logical flag which is true if a decimal point has been
!		specified in the format descriptor. This is used only to
!		decide whether to represent zero as '0.' (ldot true) or
!		'0' (ldot false). (input)
!
!	k	output unit number to which the number is printed (without
!		concluding newline character). (input)
!
!  REVISION HISTORY:
!  1.1	corrected overrun on 7 sig digit numbers. (Jul 1/97)
!  1.2	calling routine now provides the iw.id format (Oct 1/01)
!-------------------------------------------------------------------------
      subroutine prfmtf(val,iw,id,ldot,k)
      character fstr*256, d(10)*1
      logical ldot, lround
!					basic digits
      data d/'0','1','2','3','4','5','6','7','8','9'/

   1  format(a,$)
!					get absolute value of val
      aval = abs(val)
!					check 1.e-6 < |val| < 1.e+7
      if( aval .ne. 0.0 ) then
         if( aval .lt. 1.e-6 .or. aval .gt. 1.e+7 ) then
            call prfmte(val,iw,id,k)
            return
         endif
      endif
!					transfer iw and id to temp vars
      jw = iw
      jd = id
!					check specification; iw > id
      md = 2
      if( val .lt. 0.0 ) md = 3
      if( (jw .ge. 0) .and. (jd .ge. 0) ) then
         if( jd .ge. jw ) jd = jw - md	! we like to keep at least 0.
         if( jd .lt.  0 ) jd = 0
      endif
!					pointer into fstr
      m = 1
!					special case: val = 0 ---------------
      if( val .eq. 0.0 ) then
         if( (jw .lt. 0) .and. (jd .lt. 0) ) then
            if( ldot ) then
               write(k,1) '0.'
            else
               write(k,1) '0'
            endif
         elseif( (jw .ge. 0) .and. (jd .lt. 0) ) then
            if( ldot ) then
               do 30 i = 1, jw-2
                  fstr(m:m) = ' '
                  m = m + 1
  30           continue
               fstr(m:m+1) = '0.'
               write(k,1) fstr(1:m+1)
            else
               do 40 i = 1, jw-1
                  fstr(m:m) = ' '
                  m = m + 1
  40           continue
               fstr(m:m) = '0'
               write(k,1) fstr(1:m)
            endif
         elseif( (jw .lt. 0) .and. (jd .ge. 0) ) then
            write(k,1) '0'
         else				! both iw and id are specified
            ib = jw - (jd + 2)
            do 50 i = 1, ib
               fstr(m:m) = ' '
               m = m + 1
  50        continue
            fstr(m:m+1) = '0.'
            m = m + 1
            do 60 i = 1, jd
               m = m + 1
               fstr(m:m) = '0'
  60        continue
            write(k,1) fstr(1:m)
         endif
         return
      endif
!					set sign ---------------------------
      if( val .lt. 0.0 ) then
         fstr(1:1) = '-'
         m = 2
      endif
      lround = .false.
!					derive magnitude of number
  70  al = alog10(aval)
      im = int( al )
      if( (aval .lt. 1.0) .and. (aval .ne. 10.**im) ) im = im - 1

!					get number of digits to left of decimal
      il = max0( md-1, im+md-1 )	! we like to keep at least the 0.

!					adjust iw if necessary and
!					derive roundoff factor
      if( jw .ge. 0 ) then
         if( jd .ge. 0 ) then
            iwr = jw - (jd + 1)
            if( iwr .lt. il ) then
               jd = jd - (il - iwr)
               if( jd .lt. 0 ) then
                  jw = jw - jd
                  jd = 0
               endif
            endif
            round = 0.5*10.**(-jd)	! show id digits to right of .
         else
            if( jw .lt. il ) jw = il
            if( aval .lt. 1.0 ) then
               iwr = min0( jw-md, 6-im ) ! for large iw, show 7 sig digits
               round = 0.5*10.**(-iwr)
            else
               iwr = min0( 6, jw-md ) - im
               round = 0.5*10.**(-iwr)
            endif
         endif
      else
         if( jd .lt. 0 ) then
            if( aval .lt. 1.0 ) then
               round = 0.5e-07		! show only 7 digits to right of .
            else
               iwr = 6 - im		! show 7 sig digits
               round = 0.5*10.**(-iwr)
            endif
         else
            round = 0.5*10.**(-jd)
         endif
      endif
!					apply round-off and go back to check im
      if( .not. lround ) then
         aval   = aval + round
         lround = .true.
         go to 70
      endif
!					set up the basic number string
      idot = 0
!					leading zeros if val < 1
      if( aval .lt. 1.0 ) then
         fstr(m:m+1) = '0.'
         idot = m+1
         m = m + 2
         do 80 i = 1, -im-1
            fstr(m:m) = '0'
            m = m + 1
            n = n + 1
  80     continue
      endif
!					the rest of the number, 7 digits max
!					unless >= 1e+7 (which shouldn't happen)
      n = 0				! the number of sig digits written
      t = aval
!					normalize t
      i1 = im
      t  = t/(10.**i1)
  90  if( (n .lt. 7) .or. (idot .eq. 0) ) then
         ij = int( t )
         fstr(m:m) = d(ij+1)
         n = n + 1
         m = m + 1
         if( i1 .eq. 0 ) then
            fstr(m:m) = '.'
            idot = m
            m = m + 1
         endif
         t  = (t - float( ij ))*10.
         i1 = i1 - 1
         go to 90
      endif
      m = m - 1				! points at last digit in fstr
!					eliminate trailing zeroes after .
      if( jw .lt. 0 ) then
         mm = 7 + idot
         if( jd .ge. 0 ) mm = jd+idot
         if( m .gt. mm ) m  = mm
         do 100 i = m, idot, -1
            if( fstr(i:i) .ne. '0' ) go to 110
            m = m - 1			! points at last non-zero digit
 100     continue
 110     if( m .eq. idot ) m = m - 1	! don't show decimal if nothing after
      endif

!					now write it out...

      if(     (jw .lt. 0) .and. (jd .lt. 0) ) then	! no specification
         write(k,1) fstr(1:m)
      elseif( (jw .ge. 0) .and. (jd .lt. 0) ) then
         if( aval .lt. 1.0 ) then
            nb = jw - md - 6 + im
            mm = md + min0( jw-md, 6-im )
         else
            nb = jw - md - 6
            mm = min0( m, jw )
         endif
         do 120 i = 1, nb
            write(k,1) ' '
 120     continue
         write(k,1) fstr(1:mm)
      elseif( (jw .lt. 0) .and. (jd .ge. 0) ) then
         ii = min0( idot+jd, m )
         write(k,1) fstr(1:ii)
      else				! both iw and id are specified
         do 130 i = m+1, jd + idot
            fstr(i:i) = ' '
 130     continue
         mm = jd + idot
         do 140 i = mm, jw-1
            write(k,1) ' '
 140     continue
         write(k,1) fstr(1:mm)
      endif
!					all done
      return
      end
