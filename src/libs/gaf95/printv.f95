!  *********************************************************************
!  *                                                                   *
!  *                         subroutine printv                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.61
!  Written by Gordon A. Fenton, TUNS, Sat Sep  9 13:21:15 1995
!  Latest Update: Jun 25, 2003
!
!  PURPOSE  prints real/integer values according to a C style format
!
!  This routine takes as input a format, which is a character string
!  containing C style format descriptors, along with a vector of additional
!  arguments and prints them to the given unit. Normally, this routine is
!  called by one of print1, print2, ... See their headers for usage
!  instructions.
!
!  Arguments to this routine are as follows;
!
!       k	unit number to which output is directed. (input)
!
!     fmt	format string. (input)
!
!     val	real vector of length at least MP containing the elements
!		to be printed. An element of val can be either real or integer,
!		it's type is determined by the presence of a corresponding
!	 	%f or %e (for real) or %i (for integer) descriptor. See notes
!		below on the optional format specification which may follow
!		any of these format descriptors.
!		(input)
!
!       n	the number of elements in val. Currently this must be less
!		than 256. (input)
!
!   The following % formats are recognized by this routine;
!
!	   %<space> to produce a space
!
!	   %%	to produce a % character
!
!	   %b	to produce a '\' character
!
!	   %t	to produce a tab character
!
!	   %n	to produce a newline character. Note that the newline
!		character is NOT added by default to the end of an ouput
!		string, as it is in regular fortran. It must be explicitly
!		asked for via the %n format (as in C).
!
!	   %f	to format a real value in the form '12.236'. If the value
!		is outside the range (1.e-6, 1.e+7) then the %e format is
!		automatically used (except if the value is zero). The %f
!		format can be optionally written in one of the following
!		forms;
!
!		%w.df	specifies that a field width of `w' is to be used
!			with d digits to the right of the decimal place.
!			The decimal place takes up one character in the
!			total field width. `w' is actually a minimum width
!			and will be increased if more digits are required to
!			the left of the decimal point at the expense of
!			digits to the right of the decimal, if any. If the
!			number to be printed occupies less than w digits,
!			it is right justified. Digits to the right of the
!			7'th significant digit are left blank. For example,
!			%6.2f implies that the number 12.345 will be
!			printed as ' 12.35'.
!
!		%wf	specifies a minimum field width. The number of
!			decimal places is determined as necessary. Smaller
!			numbers are right justified. At most 7 significant
!			digits will be displayed. For example, %8f implies
!			that the number 12.34500 will be printed as
!			'  12.345'.
!
!		%.df	specifies a maximum number of digits to display
!			to the right of the decimal place.
!			Trailing zeroes are not shown. For example, %.2f
!			implies that the number 12.345 will be printed
!			as '12.35'.
!
!		%.f	at the moment, this is only specifies that zero is
!			printed as '0.'. The default is '0'.
!
!	   %e	to format a real value in scientific notion. For example
!		the number 12.345 is written as 1.2345e+01 in scientific
!		notation (where the e+01 means 'times 10 to the power 1').
!		The %e format can be optionally written in the following
!		forms, where;
!
!			y = 6 if the value being printed is non-negative
!			y = 7 if the value being printed is negative
!
!		%w.de	specifies that a field width of `w' is to be used
!			with d digits to the right of the decimal place.
!			w must be at least (d+y) to accomodate the leading
!			digit, the decimal, and the trailing 'e+dd' characters.
!			If w < (d+y) it is internally set to (d+y). If w
!			is greater than (d+y), the number is right justified
!			with leading blanks. For example, %10.2e implies that
!			the number 12.345 will be printed as '  1.23e+01'.
!
!		%we	specifies the field width only. The number of decimal
!			digits is determined as d = (w-y), so w must be at
!			least y (if it isn't, it is set to y). Since this is
!			a single precision routine, values of w greater than
!			12 result in a 12 character right justified number
!			(at most 7 significant digits are shown). For example,
!			%10e implies that the number 12.345 will be printed
!			as '1.2340e+01'.
!
!		%.de	specifies the number of digits to display to the right
!			of the decimal place. w is computed as (d+y). If d is
!			greater than 6, it is set internally to 6 to provide
!			at most 7 significant digits. For example, %.2e implies
!			that the number 12.345 will be printed as '1.23e+01'.
!
!	   %i	to format an integer value. Normally, the minimum field width
!		required to represent the number is used. However, the %i
!		format can be optionally written in the following forms;
!
!		%wi	specifies that a field width of at least 'w' is used.
!			For example, %10i implies that the number 12 will be
!			printed as '        12'.
!
!
!  REVISION HISTORY:
!  1.1	turned fw.0 into iw
!  1.2	added %fw specification, where w is interpreted as the desired max.
!	number of significant digits to show on real output (Oct. 9, 1995)
!  1.21	added comment about trailing blanks in fmt not being printed (10/14/96)
!  1.3	revised %i, %f, and %e output formatting (Jun 9/97)
!  1.4	replaced \\ with char(92) for portability (Dec 1/99)
!  1.5	replaced '\n' with char(10) for portability (Sep 24/01)
!  1.6	replaced all \ escapes with % escapes (for compatibility with
!	Win98's Visual Fortran). '\' is now treated as an ordinary character
!	(Oct 1/01)
!  1.61	now properly handle '% ' at end of format string (in fact '%' at
!	end of format string will always produce a space now). (Jun 25/03)
!-------------------------------------------------------------------------
      subroutine printv( k, fmt, val, n )
      parameter (MPX = 256)
      character*(*) fmt
      character*1 bslash, spc, tab
      real*4 val(*), tmp(MPX)
      integer*4 itmp(MPX)
      logical ldot, ldigit
      equivalence (itmp(1),tmp(1))
      external lnblnk

   1  format(a,$)
   3  format()

      tab    = char(9)
      spc    = char(32)
      bslash = char(92)

      MP     = n
      if( MP .gt. MPX ) MP = MPX
      do 5 i = 1, MP
         tmp(i) = val(i)
   5  continue

      jt = 1
      lf = lnblnk(fmt)
      i  = 0
  10  i = i + 1
      if( i .gt. lf ) return
      if( fmt(i:i) .eq. '%' ) then				! DESCRIPTOR
         i = i + 1
         if( i .gt. lf ) then
            write(k,1) spc
            return
         endif
         if( fmt(i:i) .eq. 'b' ) then
            write(k,1) bslash
            go to 10
         elseif( fmt(i:i) .eq. 'n' ) then
            write(k,3)
            go to 10
         elseif( fmt(i:i) .eq. 't' ) then
            write(k,1) tab
            go to 10
         elseif( fmt(i:i) .eq. '%' ) then
            write(k,1) '%'
            go to 10
         elseif( fmt(i:i) .eq. spc ) then
            write(k,1) spc
            go to 10
         elseif( ldigit(fmt(i:i)) ) then
            call getfsp(fmt,lf,i,j,iw,id,ldot)
            if( fmt(j:j) .eq. 'f' ) then
               call prfmtf(tmp(jt),iw,id,ldot,k)
            elseif( fmt(j:j) .eq. 'e' ) then
               call prfmte(tmp(jt),iw,id,k)
            elseif( fmt(j:j) .eq. 'i' ) then
               call prfmti(tmp(jt),iw,k)
            else
               write(k,1) fmt(i-1:j)
               i = j
               go to 10
            endif
            i = j
         elseif( fmt(i:i) .eq. 'f' ) then
            call prfmtf(tmp(jt),-1,-1,.false.,k)
         elseif( fmt(i:i) .eq. 'e' ) then
            call prfmte(tmp(jt),-1,-1,k)
         elseif( fmt(i:i) .eq. 'i' ) then
            call prfmti(tmp(jt),-1,k)
         else
            write(k,1) fmt(i-1:i)
            go to 10
         endif
         jt = jt + 1
         if( jt .gt. MP ) jt = MP
         go to 10
      else							! CHARACTER
         write(k,1) fmt(i:i)
      endif
      go to 10

      end
