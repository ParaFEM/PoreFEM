!  *********************************************************************
!  *                                                                   *
!  *                         subroutine prfmte                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.2
!  Written by Gordon A. Fenton, TUNS, Jun  4, 1997
!  Latest Update: Oct 2, 2001
!
!  PURPOSE  prints a number using e format
!
!  This routine writes the real value `val' to an internal character string
!  using the format
!
!	d.ddde+ee
!
!  where `ddd' expands into as many digits as necessary, or as specified
!  by the iw and id arguments, to the left of the decimal place, and `ee'
!  represents a power of ten. For example, if the number 123.4567 is to be
!  printed using the format e12.5 (see the iw.id description below), then the
!  printed value will appear as
!
!	1.23456e+02
!
!  in scientific notation.
!
!  If both iw and id are negative, then the actual format used is e11.5 if
!  val is positive, and e12.5 if val is negative.
!
!  If either iw or id are non-negative, then the format used to represent the
!  number is one of the following, where the variable y takes values
!
!			y = 6 if val is non-negative
!			y = 7 if val is negative
!
!
!	iw.id	(both iw and id are non-negative) This specifies that a
!		field width of iw is to be used with id digits to the
!		right of the decimal place. iw must be at least (id+y)
!		to accomodate the leading digit, the decimal, and the
!		trailing 'e+dd' characters. If iw < (id+y) it is internally
!		set to (id+y). If iw is greater than (id+y), the number is
!		right justified with leading blanks.
!
!	iw	(id is negative) This specifies the field width only. The
!		number of decimal digits is determined as id = (iw-y), so
!		iw must be at least y (if it isn't, it is set to y). Since
!		this is a single precision routine, values of iw greater than
!		12 result in a 12 character right justified number
!		(that is, at most 7 significant digits are shown).
!
!	.id	(iw is negative) This specifies the number of digits to
!		display to the right of the decimal place. iw is computed
!		as (id+y). If id is greater than 6, it is set internally
!		to 6 to provide at most 7 significant digits.
!
!
!  Arguments to this routine are as follows;
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
!	k	output unit number to which the number is printed (without
!		concluding newline character). (input)
!
!  REVISION HISTORY:
!  1.1	corrected incorrect rendering of powers greater than 9. (Jul 28/97)
!  1.2	calling routine now provides iw.id format (Oct 2/01)
!-------------------------------------------------------------------------
      subroutine prfmte(val,iw,id,k)
      character fstr*256
      character d(10)*1
      data d/'0','1','2','3','4','5','6','7','8','9'/

   1  format(a,$)
!					transfer iw and id to temp vars
      jw = iw
      jd = id
!					get absolute value of val
      aval = abs(val)
!					check specification; iw > id + 6
      idm = jd
      iwm = jw
      iy  = 6
      if( val .lt. 0 ) iy = 7
      if( (jw .gt. 0) .and. (jd .gt. 0) ) then
         if( jw .lt. jd+iy ) jw = jd + iy
      elseif( jw .gt. 0 ) then
         if( jw .lt. iy ) jw = iy
         jd = jw - iy
      elseif( jd .gt. 0 ) then
         if( jd .gt. 6 ) jd = 6
         jw = jd + iy
      else
         jw = 5 + iy
         jd = 5
      endif
!					special case: val = 0 ---------------
      if( val .eq. 0.0 ) then
         if( (idm .lt. 0) .and. (iwm .lt. 0 ) ) then
            write(k,1) '0.e+00'
            return
         endif
         do 10 i = jd+7, jw
            write(k,1) ' '
  10     continue
         write(k,1) '0.'
         do 20 i = 1, jd
            write(k,1) '0'
  20     continue
         write(k,1) 'e+00'
         return
      endif
!					derive magnitude of val
  30  al = alog10(aval)
      im = int( al )
      if( (aval .lt. 1.0) .and. (aval .ne. 10.**im) ) im = im - 1
!					normalize
      aval = aval/(10.**im)
!					round val
      aval = aval + 0.5*10.**(-jd)
!					recheck val to see if it has changed
      ial = int( alog10(aval) )
      if( ial .ne. 0 ) then
         im = im + 1
         aval = aval/10.
      endif
!					pointer into fstr
      m = 1
!					set sign
      if( val .lt. 0.0 ) then
         fstr(m:m) = '-'
         m = m + 1
      endif
      t = aval
      ij = int( t )
      fstr(m:m) = d(1+ij)
      fstr(m+1:m+1) = '.'
      m = m + 2
      t = t - float( ij )
      mm = jd + m
  40  if( m .lt. mm ) then
         t  = t*10.
         ij = int( t )
         fstr(m:m) = d(1+ij)
         m  = m + 1
         t  = t - float( ij )
         go to 40
      endif
!					optionally eliminate trailing zeroes
      if( idm .lt. 0 ) then
         do 50 i = m-1, 1, -1
            if( fstr(i:i) .ne. '0' ) go to 60
            m  = m - 1
            jd = jd - 1
  50     continue
      endif
  60  fstr(m:m) = 'e'
      m = m + 1
      if( im .lt. 0 ) then
         fstr(m:m) = '-'
      else
         fstr(m:m) = '+'
      endif
      iam = iabs(im)
      if( iam .gt. 99 ) then		! this shouldn't happen in real*4
         if( iam .gt. 999 ) then
            fstr(m+1:m+3) = '***'
            m = m + 3
         else
            i1 = iam/100
            i2 = (iam - 100*i1)/10
            i3 = (iam - 100*i1 - 10*i2)
            m = m + 1
            fstr(m:m) = d(1+i1)
            m = m + 1
            fstr(m:m) = d(1+i2)
            m = m + 1
            fstr(m:m) = d(1+i3)
         endif
      elseif( iam .gt. 9 ) then
         i1 = iam/10
         i2 = iam - 10*i1
         m = m + 1
         fstr(m:m) = d(1+i1)
         m = m + 1
         fstr(m:m) = d(1+i2)
      else
         m = m + 1
         fstr(m:m) = '0'
         m = m + 1
         fstr(m:m) = d(1+iam)
      endif
!					output val
      if( (idm .lt. 0) .and. (iwm .lt. 0) ) jw = jd + iy
      do 70 i = jd+iy+1, jw
         write(k,1) ' '
  70  continue
      write(k,1) fstr(1:m)
!					all done
      return
      end
