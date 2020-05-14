!  *********************************************************************
!  *                                                                   *
!  *                         subroutine las1i                          *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 3.33
!  Written by Gordon A. Fenton, TUNS, Feb. 6, 1993
!  Latest Update: Aug 10, 2005
!
!  PURPOSE  sets up initial intercell covariance parameters for LAS1G.
!
!  This routine sets up the matrices required by LAS1G to construct
!  realizations of the random field. The covariances between local averages
!  at each subdivision stage are computed in double precision for accuracy,
!  but the final construction matrices are stored in single precision and
!  returned to the calling routine via the argument list. The general
!  recursive field construction follows the relationship
!
!         {Z^j} = [A^T]{Z^(j-1)} + [C]{U}
!
!  where {Z^j} is a vector of length 2 representing the values assigned to
!  the 2-cell subdivision and {Z^(j-1)} are the parent cell values in
!  some neighbourhood. The following figure illustrates the subdivision of
!  cell 3, and its neighbourhood,
!
!    ------------------------------------------------------------------
!    |            |            |     |      |            |            |
!    |            |            |  1' |  2'  |            |            |
!    |     1      |     2      |     3      |     4      |     5      |
!    |            |            |     |      |            |            |
!    |            |            |     |      |            |            |
!    ------------------------------------------------------------------
!
!  If NBH = 3, then the neighbourhood of subdivided cells (1' and 2') are the
!  parent cells 2, 3, and 4. If NBH = 5 then the neighbourhood consists of
!  parent cells 1 through 5.
!
!  The first stage of the simulation involves the direct generation of k1
!  cells, where k1 is an integer in the range [1,MXK] which satisfies the
!  relationship N = k1*2**m. k1 is chosen to be as large as possible in
!  the given range. The maximum value of m is set by the calling routine
!  (LAS1G) by the argument MXM.
!  Arguments to this routine are as follows;
!
!  dvarfn  external real*8 function which returns the covariance of the random
!          process between two points separated by a distance V1.
!          DVARFN is referenced as follows
!
!                var = dvarfn( V1 )
!
!          where (V1) is the separation distance. Any other
!          parameters to the function must be passed by common block from the
!          calling routine. The current version of LAS uses the sign on `var'
!          to tell dvarfn to return either the covariance, or the variance of
!          a local average over distance V1. The latter is returned if var < 0,
!          although this feature is no longer used by LAS.
!
!       N  number of cells to discretize the field. N must satisfy the
!          relationship N = k1*2**m where k1 is an integer in the range
!          [1,MXK] and m is an integer less than or equal to MXM. (input)
!
!      XL  physical dimensions of the process. (input)
!
!     MXM  integer giving the largest value that M can take. An error is
!          generated if the process size is such that M > MXM. (input)
!
!     MXK  integer giving the largest value that k1 can take. (input)
!
!     NBH  integer giving the desired neighborhood size. NBH can be either
!          3 or 5. (input)
!
!       A  real array of size at least 9 x m containing parameters required
!          for the field construction. (output)
!
!       C  real array of size at least 3 x m containing parameters required
!          for the field construction. (output)
!
!      C0  real array of size at least MXK x MXK containing the upper
!          triangular values of the Cholesky decomposition of the covariance
!          matrix for the initial stage of k1 cells. (output)
!
!       M  the number of subdivisions to perform. It is an error for
!          M to be greater than MXM. (output)
!
!      k1  integers giving the size of the initial field (see C0). It is an
!          error for k1 to exceed MXK. (output)
!
!    iout  unit number to which error and warning messages are to be logged.
!          (input)
!
!     tol  the maximum relative error allowed in the Cholesky decomposition
!          of any covariance matrices before a warning message is issued.
!          (input)
!---------------------------------------------------------------------------
!  PARAMETERS:
!   MXmxk the maximum value of MXK that can be accomodated.
!
!   MXnbh the maximum value of NBH that can be accomodated. (Note that if NBH
!         greater than 5 is to be used, this routine must be modified
!         extensively).
!
!  Requires:
!    1) from libGAFsim:	DCVIT1, DCHOL2, DSIFA, DSISL, DAXPY, DSWAP,
!			IDAMAX, DDOT
!    2) external user defined variance function (see DVARFN).
!
!  Notes:
!    1) the treatment of the Cholesky decomposition is not as well done here
!       as in LAS2G and LAS3G. Some day I'll have to follow the algorithm
!       derived in DCHOL2. (GAF, Feb. 22, 1994)
!
!  REVISION HISTORY:
!  3.3	attempt recovery on Cholesky decomposition failure. (Jul 23, 1996)
!  3.31	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!  3.32	modified writeup above for new dvarfn return value. (Apr 17/01)
!  3.33	increased MXK from 16 to 256 in las1g, so MXmxk increased to 256
!	here as well. (Aug 10/05)
!----------------------------------------------------------------------------
      subroutine las1i( dvarfn, n, xl, MXM, MXK, NBH, A, C, c0,m, k1, iout, tol )
      implicit real*8 (d)
      parameter( MXmxk = 256, MXnbh = 5 )
      real A(9,*), C(3,*), c0(MXK,*)
      real*8 dg(MXnbh*MXnbh), dc0(MXmxk,MXmxk)
      real*8 daf(MXnbh), dr(MXmxk), ds(MXnbh)
      integer indx(MXnbh)
      external dvarfn

   1  format(a)
   2  format(a,e13.4)
   3  format(a,i4,a,i4,a)

      if( NBH .ne. 3 .and. NBH .ne. 5 ) then
          write(iout,1)'Warning: only neighborhood sizes of 3 or 5 are supported by LAS1G.'
          write(iout,1)'         Using a neighborhood size of 3.'
          NBH = 3
      endif
!                                       determine number of subdivisions
      do 10 k1 = MXK, 2, -1
         nbk = n/k1
         if( nbk*k1 .eq. n ) go to 20
  10  continue
      k1  = 1
      nbk = n
  20  m   = ntwom(nbk)
      if( (2**m) .ne. nbk .or. m .gt. MXM ) then
         write(iout,1)'Error: unable to determine an acceptable combination of k1 and m'
         write(iout,1)'       such that k1*2**m = N'
         write(iout,3)'       k1 = ',k1,' (must be less than ',MXK,')'
         write(iout,3)'       m  = ',m, ' (must be less than ',MXM,')'
         write(iout,1)'       Try changing N.'
         stop
      endif
!				form initial covariance matrix
      d1 = dble(xl)/dble(k1)
      call dcvit1( dvarfn, dr, ds, k1, NBH, d1 )
!					distribute into upper triangle
      do 30 j = 1, k1
      do 30 i = 1, j
         dc0(i,j) = dr(j-i+1)
  30  continue
!					and compute its Cholesky decomposition
      call dchol2( dc0, MXmxk, k1, rerr )
      if( rerr .gt. tol ) then
         write(iout,1)'Warning: Cholesky decomposition of stage 0 covariance matrix'
         write(iout,2)'         has maximum relative error of ',rerr
      endif
!					transfer to real*4
      do 40 j = 1, k1
      do 40 i = 1, j
         c0(i,j) = dc0(i,j)
  40  continue
!					compute parameters for each stage
      if( nbh .eq. 3 ) then
!						for a 3-neighborhood
         do 50 i = 1, m
            dd1  = dr(1)*dr(1)
            dd2  = dr(2)*dr(2)
            dd3  =(dr(1)-dr(3))*(dd1 + dr(1)*dr(3) - 2.d0*dd2)
            if( dd3 .eq. 0.d0 ) then
               write(iout,1)'Error: unable to factorize covariance matrix in LAS1I.'
               stop
            endif

            df1 = (dd1 - dr(2)*ds(3))/(dd1-dd2)
            df2 = (dr(1)*ds(3) - dr(1)*dr(2))/(dd1-dd2)
            df3 = (dd1*(ds(1)-dr(2)) + dd2*(ds(3)-ds(1)) &
                   + dr(1)*dr(3)*(dr(2)-ds(3)))/dd3

            db1 = df1*dr(1) + df2*ds(3)
            db2 = df3*(ds(1)-ds(3)) + ds(2)

            A(1,i) = df1
            A(2,i) = df2
            A(3,i) = df3

            d1 = 0.5d0*d1
            call dcvit1( dvarfn, dr, ds, NBH, NBH, d1 )

            if( dr(1) .lt. db1 ) then
               write(iout,3)'Warning: stage ',i, &
              ' covariance matrix not positive definite (db1).'
               write(iout,1)'         attempting recovery.'
               dr(1) = 1.1d0*db1
            endif
            if( dr(1) .lt. db2 ) then
               write(iout,3)'Warning: stage ',i, &
              ' covariance matrix not positive definite (db2).'
               write(iout,1)'         attempting recovery.'
               dr(1) = 1.1d0*db2
            endif
            C(1,i) = dsqrt(dr(1) - db1)
            C(2,i) = dsqrt(dr(1) - db2)
  50     continue
      else
!						for a neighborhood of 5
         do 60 i = 1, m
!							edge case
            dg(1) = dr(1)
            dg(4) = dr(2)
            dg(5) = dr(1)
            dg(7) = dr(3)
            dg(8) = dr(2)
            dg(9) = dr(1)
            call dsifa( dg, 3, 3, indx, ierr )
            if( ierr .ne. 0 ) then
               write(iout,1)'Error: unable to factorize edge covariance matrix in LAS1I.'
               stop
            endif
            daf(1) = ds(3)
            daf(2) = ds(4)
            daf(3) = ds(5)
            call dsisl( dg, 3, 3, indx, daf )
            A(1,i) = daf(1)
            A(2,i) = daf(2)
            A(3,i) = daf(3)
            db1 = daf(1)*ds(3) + daf(2)*ds(4) + daf(3)*ds(5)

!							edge case + 1
            dg( 1) = dr(1)
            dg( 5) = dr(2)
            dg( 6) = dr(1)
            dg( 9) = dr(3)
            dg(10) = dr(2)
            dg(11) = dr(1)
            dg(13) = dr(4)
            dg(14) = dr(3)
            dg(15) = dr(2)
            dg(16) = dr(1)
            call dsifa( dg, 4, 4, indx, ierr )
            if( ierr .ne. 0 ) then
               write(iout,1)'Error: unable to factorize near edge covariance matrix in LAS1I.'
               stop
            endif
            daf(1) = ds(2)
            daf(2) = ds(3)
            daf(3) = ds(4)
            daf(4) = ds(5)
            call dsisl( dg, 4, 4, indx, daf )
            A(4,i) = daf(1)
            A(5,i) = daf(2)
            A(6,i) = daf(3)
            A(7,i) = daf(4)
            db2 = daf(1)*ds(2)+daf(2)*ds(3)+daf(3)*ds(4)+daf(4)*ds(5)

!							interior case
            dg( 1) = dr(1)
            dg( 6) = dr(2)
            dg( 7) = dr(1)
            dg(11) = dr(3)
            dg(12) = dr(2)
            dg(13) = dr(1)
            dg(16) = dr(4)
            dg(17) = dr(3)
            dg(18) = dr(2)
            dg(19) = dr(1)
            dg(21) = dr(5)
            dg(22) = dr(4)
            dg(23) = dr(3)
            dg(24) = dr(2)
            dg(25) = dr(1)
            call dsifa( dg, 5, 5, indx, ierr )
            if( ierr .ne. 0 ) then
               write(iout,1)'Error: unable to factorize interior covariance matrix in LAS1I.'
               stop
            endif
            daf(1) = ds(1)
            daf(2) = ds(2)
            daf(3) = ds(3)
            daf(4) = ds(4)
            daf(5) = ds(5)
            call dsisl( dg, 5, 5, indx, daf )
            A(8,i) = daf(1)
            A(9,i) = daf(2)
            db3 = daf(1)*(ds(1)-ds(5)) + daf(2)*(ds(2)-ds(4)) + ds(3)

!						compute C elements
            d1 = 0.5d0*d1
            call dcvit1( dvarfn, dr, ds, NBH, NBH, d1 )
            if( dr(1).lt.db1 .or. dr(1).lt.db2 .or. dr(1).lt.db3) then
               write(iout,1)'Error: Cholesky decomposition of covariance matrix in LAS1I'
               write(iout,1)'         has imaginary components.'
               stop
            endif
	    C(1,i) = dsqrt(dr(1) - db1)
	    C(2,i) = dsqrt(dr(1) - db2)
	    C(3,i) = dsqrt(dr(1) - db3)
  60     continue
      endif

      return
      end
