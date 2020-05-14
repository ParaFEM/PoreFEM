!  *******************************************************************
!  *                                                                 *
!  *                    subroutine vnorm                             *
!  *                                                                 *
!  *******************************************************************
!  Single Precision Version 2.21
!  Written by Gordon A. Fenton, TUNS, Mar. 21, 1993.
!  Latest Update: Jun 9, 1999
!
!  PURPOSE  return a vector of N(0,1) random distributed realizations
!
!  Returns a normally distributed, zero mean, unit variance random vector
!  in the argument `g'. `n' is the length of the desired vector.
!
!  Notes:
!   1) Ensure that the random number generator is initialized prior to calling
!      this routine.
!   2) since randu returns random numbers in the range (0,1), excluding
!      the endpoints, we need not guard against taking the log of zero
!      below.
!
!  REVISION HISTORY:
!  2.1	now using new randu function (RAN2 from Numerical Recipes, 2nd Ed)
!	(Oct 14/96)
!  2.2	eliminated unused local variables `zero' and `big' (Dec 5/96)
!  2.21	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!---------------------------------------------------------------------------
      subroutine vnorm2( g, n )
      real g(*)
      data twopi/6.2831853071795864769/
      data two/2.0/
!
      nn = n/2
      nn = 2*nn
      do 10 i = 1, nn, 2
         a = twopi*randu2(0)
         r = sqrt(-two*alog(randu2(0)))
         g(i)   = r*cos(a)
         g(i+1) = r*sin(a)
  10  continue

      if( n .eq. nn ) return
!                                     for n odd, set the last value
      a = twopi*randu2(0)
      r = sqrt(-two*alog(randu2(0)))
      g(n) = r*cos(a)

      return
      end
