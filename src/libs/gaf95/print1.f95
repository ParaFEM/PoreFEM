!  *********************************************************************
!  *                                                                   *
!  *                         subroutine print1                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.2
!  Written by Gordon A. Fenton, TUNS, Sat Sep  9 13:21:15 1995
!  Latest Update: Oct 4, 2001
!
!  PURPOSE  prints one value according to a C style format
!
!  This routine takes as input a format, which is a character string
!  containing C style format descriptors along with one additional
!  argument and prints it to the given unit. For example, say you
!  wish to print the variable d to unit 6 using the following format;
!
!	d = 16.49
!
!  (assuming d's value is 16.49), then this routine would be called as
!  follows;
!
!	call print1( 6, 'd = %f%n', d )
!
!  Note that if the trailing '%n' characters are absent from the format string,
!  then the carriage return is inhibited in the output.
!
!  Arguments to this routine are as follows;
!
!       k	unit number to which output is directed. (input)
!
!     fmt	format string. (input)
!
!     val	the input argument. val can be either real or integer,
!		it's type is determined by the presence of a corresponding
!	 	%f or %e (for real) or %i (for integer) descriptor. See notes
!		below on the optional format specification which may follow
!		any of these format descriptors.
!		(input)
!
!   The following % formats are recognized by this routine; only the %f, %e,
!   and %i formats actually consume an argument (ie, make use of one of
!   the val argument).
!
!	   %<space> to produce a space
!
!	   %%	to produce a % character
!
!	   %b	to produce a '\' character
!
!	   %t	to produce a tab character
!
!	   %n	to produce a newline character. This is NOT added to the
!		end of the output string by default. It must be explicitly
!		added (as in C).
!
!	   %f	to format a real value in the form '12.236'. If the value
!		is outside the range (1.e-6, 1.e+7) then the %e format is
!		automatically used. The %f format can be optionally written
!		in the following forms;
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
!		format can be optionally written in the following form;
!
!		%wi	specifies that a field width of at least 'w' is used.
!			For example, %10i implies that the number 12 will be
!			printed as '        12'.
!
!  NOTES:
!   1)	if more than 1 occurrence of the formats %f, %e or %i appear in
!	the format string, the argument value is repeatedly used (this may
!	lead to problems if it's type doesn't match the repeated format
!	descriptor).
!
!   2)	trailing white space in the format string is not printed in the output.
!	To actually get it in the output, it must be escaped, as in % % % ...
!
!  REVISION HISTORY
!  1.1	added %fw specification, where w is interpreted as the desired max.
!	number of significant digits to show on real output (Oct 9, 1995)
!  1.11	added comment about trailing blanks in fmt not being printed (10/14/96)
!  1.12	modified above writeup for new prints routines (Jun 9/97)
!  1.2	modified above writeup to reflect new formats and escapes (Oct 4/01)
!-------------------------------------------------------------------------
      subroutine print1( k, fmt, val )
      character*(*) fmt
      real val

      call printv( k, fmt, val, 1 )

      return
      end
