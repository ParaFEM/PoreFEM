!  *********************************************************************
!  *                                                                   *
!  *                         subroutine las2i                          *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 3.41
!  Written by Gordon A. Fenton, TUNS, Aug. 26, 1992
!  Latest Update: Jun 9, 1999
!
!  PURPOSE   initializes parameters for LAS2G
!
!  This routine sets up the matrices required by LAS2G to construct
!  realizations of the random field. The covariances between local averages
!  at each subdivision stage are computed in double precision for accuracy,
!  but the final construction matrices are stored in single precision and
!  returned to the calling routine via the argument list. The general
!  recursive field construction follows the relationship
!
!         {Z^j} = [A^T]{Z^(j-1)} + [C]{U}
!
!  where {Z^j} is a vector of length 4 representing the values assigned to
!  the 2 x 2 cell subdivision and {Z^(j-1)} are the parent cell values in
!  some neighbourhood. The following figure illustrates the cell subdivision,
!  the neighbourhood, and the numbering scheme for an interior cell (special
!  subsets of the neighbourhood are used for the corners and sides):
!
!                ----------------------------------------
!                |            |            |            |
!                |            |            |            |
!                |     7      |     8      |      9     |
!                |            |            |            |
!                |            |            |            |
!                |------------|------------|------------|
!                |            | 3   |    4 |            |
!                |            |     |      |            |
!                |     4      |-----5------|      6     |
!                |            |     |      |            |
!                |            | 1   |    2 |            |
!                |------------|------------|------------|
!                |            |            |            |
!                |            |            |            |
!                |     1      |     2      |     3      |
!                |            |            |            |
!                |            |            |            |
!                ----------------------------------------
!
!  We see that for the upper left corner, the parent cell neighbourhood used
!  consists of just cells {4,5,7,8} and similarly for the other corners and
!  sides.
!
!  The first stage of the simulation involves the direct generation of a
!  k1 x k2 cell array, where k1 and k2 are integers which satisfying the
!  decomposition N1 = k1*2**m, N2 = k2*2**m for a common factor 2**m. The
!  integers k1 and k2 are chosen to be as large as possible while requiring
!  the product k1*k2 to be less than or equal to MXK. Note that the direct
!  simulation involves the inversion of a MXK x MXK matrix (at the upper
!  limit) and so MXK should not be overly large.
!  This formulation is somewhat less restrictive than simply requiring
!  N1 and N2 to be powers of 2. Also N1 and N2 do not have to be equal.
!  However N1 and N2 cannot still be chosen arbitrarily, for example the
!  set (N1,N2) = (336,256) results in k1 = 21, k2 = 16, m = 4 which is
!  not acceptable here (since k1*k2 > MXK for MXK = 256), while the set
!  (N1,N2) = (160,256) is acceptable since k1 = 10, k2 = 16, m = 4. In general
!  it may be easier to choose k1, k2, and m before specifying N1 and N2. In
!  the event that an unacceptable (k1,k2,m) combination is selected, IERR
!  is set to -1 and control returned to the calling routine.
!  The maximum value of m is set by the calling routine in the argument MXM.
!
!  Arguments to this routine are as follows;
!
! dvarfn  external real*8 function which returns the variance of the random
!         process averaged over a given area. dvarfn is referenced as follows
!
!                var = dvarfn( V1, V2 )
!
!         where (V1,V2) are the side dimensions of the rectangular averaging
!         domain. Any other parameters to the function may be passed by
!         common block from the calling routine. Note that the variance of
!         the random process averaged over the area (V1 x V2) is the product
!         of the point variance and the traditionally defined "variance"
!         function, as discussed by Vanmarcke (pg 186).
!
!  N1,N2  number of cells to discretize the field in the x and y directions
!         respectively (corresponding to the first and second indices of Z
!         respectively). Both N1 and N2 must have the form N1 = k1 * 2**m
!         and N2 = k2 * 2**m where m is common to both and k1 and k2 are
!         positive integers satisfying k1*k2 <= MXK. Generally k1 and k2
!         are chosen to be as large as possible and still satisfy the above
!         requirements so the the first stage involves directly simulating
!         a k1 x k2 field by inversion of a covariance matrix.
!         A potential example is (N1,N2) = (160,256) which gives k1 = 5,
!         k2 = 8, and m = 5. Note that in general N1 and N2 cannot be chosen
!         arbitrarily - it is usually best to choose m first then
!         k1 and k2 so as to satisfy or exceed the problem requirements. Note
!         that because of the requirements on k1*k2, N1 cannot be more than
!         MXK times as big (or smaller) than N2. (input)
!
!  XL,YL  physical dimensions of the process. (input)
!
!  KSEED  integer seed to be used for the pseudo-random number generator.
!         If KSEED = 0, then a random seed will be used (based on the
!         process ID of the current program -- see iseed.f).
!         On output, KSEED is set to the value of the actual seed used.
!
!    MXM  integer giving the largest value that M can take. An error is
!         generated if the process size is such that M > MXM. (input)
!
!    C0   real vector containing the upper triangular values of the Cholesky
!         decomposition of the covariance matrix for the initial stage (0) of
!         k1 x k2 cells. (output)
!
!    CT   real vector containing the upper triangular values of the Cholesky
!         decomposition of the covariance matrix for the k1 or k2 = 1 special
!         case. (output)
!
!    CC   real vector containing the upper triangular values of the Cholesky
!         decomposition of the covariance matrix for the corner cell 2 x 2
!         subdivisions. (output)
!
!    CS   real vector containing the upper triangular values of the Cholesky
!         decomposition of the covariance matrix for the side cell 2 x 2
!         subdivisions. (output)
!
!    CI   real vector containing the upper triangular values of the Cholesky
!         decomposition of the covariance matrix for the interior cell 2 x 2
!         subdivisions. (output)
!
!    AT   real array containing the best linear estimation coefficients for
!         the k1 or k2 = 1 special case. (output)
!
!    AC   real array containing the best linear estimation coefficients for
!         the corner cell subdivisions. (output)
!
!    AS   real array containing the best linear estimation coefficients for
!         the side cell subdivisions. (output)
!
!    AI   real array containing the best linear estimation coefficients for
!         the interior cell subdivisions. (output)
!
!     M   the number of 2 x 2 subdivisions to perform. It is an error for
!         M to be greater than MXM. (output)
!
! k1,k2   integers giving the size of the initial field (see C0). It is an
!         error for the product k1*k2 to exceed MXK. (output)
!
!    kk   integer giving the size of the initial field covariance matrix
!         (kk = k1*k2). (output)
!
!  iout   unit number to which error and warning messages are logged. (input)
!
!   tol   maximum relative error allowed in the Cholesky decomposition of
!         covariance matrices before a warning message is issued. (input)
!---------------------------------------------------------------------------
!  PARAMETERS:
!   MXK   represents the maximum number of cells (k1 x k2) in the initial
!         field, k1*k2 <= MXK. If the value of MXK is changed here, it must
!         also be changed in LAS2G.
!
!  Requires:
!    1) from libGAFsim:	ISEED, DCVIT2, DCVMT2, DCHOL2, CORN2D, SIDE2D,
!			INTR2D, DCVAA2, DCVAB2, DSIFA, DSISL, DAXPY, DSWAP,
!			IDAMAX, DDOT
!    2) user defined external variance function (see DVARFN)
!
!  REVISION HISTORY:
!  3.41	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!----------------------------------------------------------------------------
      subroutine las2i( dvarfn, N1, N2, XL, YL, kseed, MXM,C0, CT, CC, CS, CI, AT, AC, AS, AI, M,k1, k2, kk, iout, tol )
      parameter ( MXK = 256 )
      real C0(*), CT(6,*), CC(6,4,*), CS(6,4,*), CI(6,*)
      real AT(3,3,*), AC(4,3,4,*), AS(6,3,4,*), AI(9,3,*)
      real*8 R0(MXK*MXK)
      real*8 R(9,9,2), B(4,4), S(9,4)
      real*8 T1, T2, dvarfn, dble
      logical lformR
      integer mc(4,4), ms(6,4), mi(9)
      external dvarfn
      data mc/5,6,8,9, 4,5,7,8, 2,3,5,6, 1,2,4,5/
      data ms/4,5,6,7,8,9, 2,3,5,6,8,9, 1,2,4,5,7,8, 1,2,3,4,5,6/
      data mi/1,2,3,4,5,6,7,8,9/

   1  format(a,a,a)
   2  format(a,e13.4)
   3  format(a,i4,a,i4,a,i4,a)
!						decompose N1 and N2
      k1 = N1
      k2 = N2
      do 10 m = 0, MXM
         kk = k1*k2
         if( kk .le. MXK ) go to 20
         j1 = k1/2
         j2 = k2/2
         if( 2*j1 .ne. k1 .or. 2*j2 .ne. k2 ) go to 15
         k1 = j1
         k2 = j2
  10  continue
  15  write(iout,1)'Error: unable to determine an acceptable combination of k1, k2 and m'
      write(iout,1)'       such that k1*2**m = N1, k2*2**m = N2.'
      write(iout,3)'       k1 = ',k1,', k2 = ',k2,', m = ',m
      write(iout,3)'       (k1*k2 must be less than ',MXK,' and m must be less than ',MXM,')'
      write(iout,1)'       Try changing N1 and N2.'
      stop
!						initialize internal generator
  20  kseed = iseed( kseed )
!						form initial covariance matrix
      T1 = dble(XL)/dble(k1)
      T2 = dble(YL)/dble(k2)
      call dcvit2( dvarfn, R0, kk, R, 9, k1, k2, T1, T2 )

!						and compute its cholesky decomp
      call dchol2( R0, kk, kk, rerr )
      if( rerr .gt. tol ) then
         write(iout,1)'Warning: Cholesky decomposition of stage 0 covariance matrix'
         write(iout,2)'         has maximum relative error of ',rerr
      endif
!						save in real*4 for LAS2G
      L = 0
      do 30 j = 1, kk
      do 30 i = 1, j
         L = L + 1
         C0(L) = R0(i+(j-1)*kk)
  30  continue
!						setup for subsequent stages
      in = 2
      io = 1
      nm = 1
      if( (k1 .eq. 1 .or. k2 .eq. 1) .and. M .gt. 0 ) then
!							special k1,k2 = 1 case
         T1 = 0.5d0*T1
         T2 = 0.5d0*T2
!							get basic cov matrices
         lformR = (1 .lt. M)
         call dcvmt2(dvarfn,R(1,1,in),9,B,4,S,9,T1,T2,lformR)
         i2 = 5
         if( k1 .eq. 1 ) then
            i1 = 2
            i3 = 8
         else
            i1 = 4
            i3 = 6
         endif
         call thin1d(R(1,1,io),9,B,4,S,9,AT,3,AT(1,1,2),3,CT,CT(1,2),i1,i2,i3,3,iout,tol)
         in = 1
         io = 2
         nm = 2
      endif

      do 40 k = nm, M
         T1 = 0.5d0*T1
         T2 = 0.5d0*T2
!							get basic cov matrices
         lformR = (k .lt. M)
         call dcvmt2(dvarfn,R(1,1,in),9,B,4,S,9,T1,T2,lformR)

!							corner parameters

         call corn2d(R(1,1,io),9,B,4,S,9,CC(1,1,k),6,3,AC(1,1,1,k),mc,iout,tol)
!							side parameters

         call side2d(R(1,1,io),9,B,4,S,9,CS(1,1,k),6,3,AS(1,1,1,k),ms,iout,tol)
!							interior parameters

         call intr2d(R(1,1,io),9,B,4,S,9,CI(1,k),3,AI(1,1,k),mi,iout,tol)
!							swap old/new indices
         ii = in
         in = io
         io = ii
  40  continue

      return
      end
