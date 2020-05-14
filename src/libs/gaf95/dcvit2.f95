!  *********************************************************************
!  *                                                                   *
!  *                         subroutine dcvit2                         *
!  *                                                                   *
!  *********************************************************************
!  Double Precision Version 2.11
!  Written by Gordon A. Fenton, TUNS, July 18, 1992
!  Latest Update: Jun 9, 1999
!
!  PURPOSE  computes the initial k1*k2 x k1*k2 covariance matrix (and a 9x9
!           submatrix) between local averages. Used by LAS2G.
!          
!
!  This routine forms the covariance matrix Q between local averages arranged
!  on a k1 x k2 grid using a user provided variance function "dvarfn" that
!  returns the variance of a local average (note that this is the point
!  variance sigma^2 times the traditionally "variance" function, as defined
!  by Vanmarcke in "Random Fields", MIT Press, 1984). The grid and element
!  numbering schemes appear as follows (for k1 = k2 = 5)
!
!          | Dx |
!     ---  --------------------------
!      Dy  | 21 | 22 | 23 | 24 | 25 |
!     ---  --------------------------
!          | 16 | 17 | 18 | 19 | 20 |
!          --------------------------           ----------------
!          | 11 | 12 | 13 | 14 | 15 |           |  7 |  8 |  9 |
!          --------------------------           ----------------
!          |  6 |  7 |  8 |  9 | 10 |           |  4 |  5 |  6 |
!          --------------------------           ----------------
!          |  1 |  2 |  3 |  4 |  5 |           |  1 |  2 |  3 |
!          --------------------------           ----------------
!
!                  Q Array                          R Array
!
!  with the k1 x k2 array numbering on the left and a basic 3 x 3 subarray
!  shown on the right. If we call Z_i the local average of a random process
!  over the i'th cell, then for E[Z] = 0, the elements of Q are defined by
!
!           Q(i,j) = E[Z_i*Z_j]
!
!  Note that the elements of R are simply extracted directly from the
!  appropriate elements of Q (for example, R(1,5) = Q(1,7)) if k1 and k2 are
!  greater than 2. Note also that since both Q and R are symmetric, ONLY THE
!  UPPER TRIANGULAR VALUES are stored herein. Finally note that the random
! process is assumed to be quadrant symmetric (that is the covariance
! function is even in each individual lag component) to reduce the total
! number of covariance calculations.
!
!  Arguments to this routine are as follows;
!
!    dvarfn     external double precision function which returns the
!               variance of a local average of size X x Y. This function
!               is referenced using the call "V = dvarfn(X,Y)" (other
!               parameters to the function must be passed via common blocks).
!               On each invocation, this routine calls DVARFN approximately
!               225 times (via DCVAA2).
!
!    Q          double precision array of size at least k1*k2 x k1*k2 which
!               on output will contain the covariance matrix between the local
!               average elements shown above on the left. Note that since Q is
!               symmetric, only the upper triangular values are actually
!               stored (compatible with DCHOL2). (output)
!
!    iq         leading dimension of the array Q as specified in the calling
!               routine. (input)
!
!    R          double precision array of size at least 9 x 9 which on output
!               will contain the covariance matrix between the local average
!               elements shown above on the right. Note that since R is
!               symmetric, only the upper triangular values are actually
!               stored (compatible with DSIFA). (output)
!
!    ir         leading dimension of the array R as specified in the calling
!               routine. (input)
!
!    k1         number of local average cells in the x direction
!               (horizontally). (input)
!
!    k2         number of local average cells in the y direction
!               (vertically). (input)
!
!    Dx         x-dimension of the local average cells. (input)
!
!    Dy         y-dimension of the local average cells. (input)
!
!  Required:
!   1) from libGAFsim: DCVAA2
!
!  REVISION HISTORY:
!  2.1	eliminated unused local variables (Dec 5/96)
!  2.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!----------------------------------------------------------------------------
      subroutine dcvit2( dvarfn, Q, iq, R, ir, k1, k2, Dx, Dy )
      implicit real*8 (a-h,o-z)
      dimension Q(iq,*), R(ir,*)
      integer kx(9), ky(9)
      external dvarfn
      data kx/0,1,2,0,1,2,0,1,2/, ky/0,0,0,1,1,1,2,2,2/

!					first form the essential elements of Q
      ii = 0
      do 20 j = 1, k2
         ty = dble(j-1)
         do 10 i = 1, k1
            ii = ii + 1
            tx = dble(i-1)
            Q(1,ii) = dcvaa2( dvarfn, Dx, Dy, tx, ty )
  10     continue
  20  continue
!					now distribute into the upper triangle
      kk = k1*k2
      do 40 j = 2, kk
         mxj = mod(j-1,k1)
         myj = (j-1)/k1
         do 30 i = 2, j
            mxi    = mod(i-1,k1)
            myi    = (i-1)/k1
            m      = 1 + iabs(mxj - mxi) + k1*iabs(myj - myi)
            Q(i,j) = Q(1,m)
  30     continue
  40  continue
!					form R (9 x 9)
      if( k1 .lt. 3 .or. k2 .lt. 3 ) then
!						find essential elements of R
         ii = 0
         do 60 j = 1, 3
            ty = dble(j-1)
            do 50 i = 1, 3
               ii = ii + 1
               tx = dble(i-1)
               R(1,ii) = dcvaa2( dvarfn, Dx, Dy, tx, ty )
  50        continue
  60     continue
!						and distribute into upper tri.
         do 70 j = 2, 9
         do 70 i = 2, j
            m  = 1 + iabs(kx(j)-kx(i)) + 3*iabs(ky(j)-ky(i))
            R(i,j) = R(1,m)
  70     continue
      else
!						or extract elements of R from Q
         do 80 j = 1, 9
         do 80 i = 1, j
            m      = 1 + iabs(kx(j)-kx(i)) + k1*iabs(ky(j)-ky(i))
            R(i,j) = Q(1,m)
  80     continue
      endif

      return
      end
