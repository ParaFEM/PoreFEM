!  *******************************************************************
!  *                                                                 *
!  *                       Subroutine las1g                          *
!  *                                                                 *
!  *******************************************************************
!  Single Precision Version 3.34
!  Written by Gordon A. Fenton, TUNS, Feb. 4, 1993
!  Latest Update: Jun 3, 2006
!
!  PURPOSE  generates a realization of a 1-D Gaussian stationary random
!           local average process
!
!  This routine creates a realization of a zero-mean 1-D random process given
!  its variance function (see DVARFN below). Each discrete
!  value generated represents the local average of a realization of
!  the process over the cell length XL/N, where XL is the physical length
!  of the process and N is the number of cells desired. N must be such that
!  N = k1*2**m, for some integer k1 in the range [1,256] and for some integer
!  m in the range [0,16]. Thus the largest value of N is 2^(24) = 16,777,216.
!  The values of k1 and m are computed internally. The construction of the
!  realization proceeds in a recursive fashion as follows;
!    1) generate the first k1 cells directly,
!    2) subdivide each of the k1 cells into two equal parts and generate
!       random values for each part preserving the cell average and closely
!       approximating the covariance structure (the latter is exact within
!       the cell),
!    3) subdivide the 2*k1 cells obtained in step (2) into two equal parts
!       and again generate random values for each part,
!    4) continue the subdivision process until the domain is divided into
!       N equal cells.
!
!  Note that this routine sets up a number of parameters on the first
!  call and thus the time required to produce the first realization
!  is substantially greater than on subsequent calls (see INIT). For
!  more details on this Local Average Subdivision algorithm, see
!
!    1) Fenton, G.A., and Vanmarcke, E.H., "Simulation of Random Fields
!          via Local Average Subdivision", ASCE Journal of Engineering
!          Mechanics, 116(8), 1733-1749, 1990.
!    2) Fenton, G.A., Simulation and Analysis of Random Fields, Ph.D. Thesis,
!          Dept. of Civil Engineering and Operations Research, Princeton
!          University, Princeton, NJ, 1990.
!    3) Vanmarcke, E.H., Random Fields: Analysis and Synthesis, MIT Press,
!          Boston, Massachusetts, 1984.
!    4) Fenton, G.A., Error evaluation of three random field generators,
!          ASCE Journal of Engineering Mechanics (to appear), 1993.
!
!  Arguments to this routine are as follows;
!
!      Z  real vector of length at least (3*N/2) which on output will contain
!         the desired realization in the first N locations. The remainder is
!         used for workspace. (output)
!
!      N  desired number of cells discretizing the process. N must be such
!         that it can be written N = k1*2**m, k1 and m are integers with
!         1 <= k1 <= MXK and 0 <= m <= MXM (see parameters below). (input)
!
!     XL  physical length of the process. (input)
!
! DVARFN  external real*8 function which returns the covariance of the random
!         process between two points separated by a distance V1.
!         DVARFN is referenced as follows
!
!                var = dvarfn( V1 )
!
!         where (V1) is the separation distance. Any other
!         parameters to the function must be passed by common block from the
!         calling routine. The current version of LAS uses the sign on `var'
!         to tell dvarfn to return either the covariance, or the variance of
!         a local average over distance V1. The latter is returned if var < 0,
!         although this feature is no longer used by LAS.
!
!  NBH    integer giving the desired size of the neighborhood to include
!         in the determination of the best linear estimate of the mean value
!         of a sub-cell. At present only two neighborhood sizes are supported:
!         NBH = 3 or 5. Note that the code will run significantly slower using
!         NBH = 5, but it yields improved results for some types of processes
!         (such as damped oscillatory noise). (input)
!
!  KSEED  integer seed to be used to initialize the pseudo-random number
!         generator (which is assumed to be performed by calling rand(kseed)).
!         If KSEED = 0, then a random seed will be used (based on the
!         clock time when this routine is called for the first time).
!         On output, KSEED is set to the value of the actual seed used.
!         (input/output)
!
!  INIT   integer flag which must be set to +/- 1 when parameters of the
!         process are to be calculated or recalculated. If multiple
!         realizations of the same process are desired, then subsequent
!         calls should use INIT not equal to 1 and N1 the same as used
!         initially. If INIT = +1, then both the parameters are calculated
!         and a random realization is returned. If INIT = -1, then just
!         the initial parameters are calculated. (input)
!
!  IOUT   unit number to which error and warning messages are to be logged.
!         (input)
!  ---------------------------------------------------------------------------
!
!  PARAMETERS:
!
!    MXK  the maximum value of k1. Currently 256.
!
!    MXM  represents the maximum number of subdivisions (ie the maximum value
!         of m) that the routine can carry out. Currently 16. Note that this
!         means that N cannot exceed MXK*2**MXM = 16*2**16 = 1,048,576.
!
!    NGS  this routine will generate NGS independent Gaussian random variates
!         at a time (this is because some machines slow down drastically when
!         making multiple calls to an external routine - that is we attempt to
!         minimize the number of calls to VNORM herein).
!
!  Notes:
!    1) Simulation timing is available through common block LASTYM where
!       TI is the time required to set up the parameters of the process
!       and TS is the cumulative simulation time. These timings are
!       obtained through the function SECOND which returns elapsed user
!       time in seconds since the start of the program.
!    2) all variables that start with "d" are double precision - this is
!       used to perform the covariance calculations in extended precision.
!    3) The parameter `tol' is the maximum relative error allowed on the
!       Cholesky decomposition of covariance matrices before a warning
!       message is emitted. The relative error is estimated by computing
!       the lower-right most element of L*L' and comparing to the original
!       element of A (where L*L' = A is the decomposition). This give some
!       measure of the roundoff errors accumulated in the computation of L.
!
!  Requires:
!    1) from libGAFsim:	ISEED, DCVIT1, DCHOL2, DSIFA, DSISL, DAXPY, DSWAP,
!			IDAMAX, DDOT, LAS1I, VNORM, RANDF, SECOND
!    2) external user defined variance function (see DVARFN).
!
!  REVISION HISTORY:
!  3.31	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!  3.32	modified writeup above for new dvarfn return value. (Apr 17/01)
!  3.33	increased MXK from 16 to 256. (Aug 10/05)
!  3.34	just initialize if init = -1. (Jun 3/06)
!--------------------------------------------------------------------
      subroutine las1g( z, n, xl, dvarfn, nbh, kseed, init, iout )
      implicit real*8 (d)
      parameter( MXM = 16, MXK = 256, NGS = 4096 )
      real z(*), A(9,MXM), C(3,MXM), C0(MXK,MXK), g(NGS), ge(4*MXM)
      external dvarfn
      save A, C, C0, m, k1
      common/LASTYM/ ti, ts
      data ifirst/1/, two/2.0/
      data tol/1.e-3/
!------------------------------ initialization ---------------------------

      if( ifirst .eq. 1 .or. iabs(init) .eq. 1 ) then
         ti = second()
         ifirst = 0
!					compute parameters for this routine
         call las1i( dvarfn, n, xl, MXM, MXK, nbh, A, C, C0,              &
                      m, k1, iout, tol )
!					initialize the random number generator
         kseed = iseed( kseed )
!					set timers
         ts = 0.
         ti = second() - ti
         if( init .eq. -1 ) return
      endif
!
!------------------------------ create the realization ---------------------
!
      tt = second()
!					create first k1 cells directly
      if( mod(m,2) .eq. 0 ) then
         iz = 0
         jz = n
      else
         iz = n
         jz = 0
      endif
!					generate stage 0 field
      call vnorm( g, k1 )
      do 20 i = 1, k1
         Z(iz+i) = C0(1,i)*g(1)
         do 10 j = 2, i
            Z(iz+i) = Z(iz+i) + C0(j,i)*g(j)
  10     continue
  20  continue
!					generate stage 1, 2, ... M fields
      jx = k1
      if ( nbh .eq. 3 ) then
!					for a neighborhood of 3 (NBH = 3)...
         call vnorm( ge, 2*M )
         do 40 i = 1, M
            ii = 2*i
            it = jz
            jz = iz
            iz = it
            i0 = iz + 1
            j0 = jz + 1
            Z(i0)   = A(1,i)*Z(j0) + A(2,i)*Z(j0+1) + C(1,i)*ge(ii-1)
            Z(i0+1) = two*Z(j0) - Z(i0)
            do 30 jj = 1, jx-2, NGS
               kk = min0( jx-jj-1, NGS )
               call vnorm( g, kk )
               do 30 j = 1, kk
                  i0      = i0 + 2
                  Z(i0)   = A(3,i)*(Z(j0)-Z(j0+2))+Z(j0+1) + C(2,i)*g(j)
                  Z(i0+1) = two*Z(j0+1) - Z(i0)
                  j0      = j0 + 1
  30        continue
            i0      = i0 + 2
            Z(i0+1) = A(1,i)*Z(j0+1) + A(2,i)*Z(j0) + C(1,i)*ge(ii)
            Z(i0)   = two*Z(j0+1) - Z(i0+1)
            jx      = 2*jx
  40     continue
      else
!					for a neighborhood of 5 (NBH = 5)...
         call vnorm( ge, 4*M )
         do 60 i = 1, M
            ii = 4*i
            it = jz
            jz = iz
            iz = it
            i0 = iz + 1
            j0 = jz + 1
            Z(i0)   = A(1,i)*Z(j0) + A(2,i)*Z(j0+1) + A(3,i)*Z(j0+2)      &
                                  + C(1,i)*ge(ii-3)
            Z(i0+1) = two*Z(j0) - Z(i0)
            i0 = i0 + 2
            Z(i0)   = A(4,i)*Z(j0) + A(5,i)*Z(j0+1) + A(6,i)*Z(j0+2)      &
                                  + A(7,i)*Z(j0+3) + C(2,i)*ge(ii-2)
            Z(i0+1) = two*Z(j0+1) - Z(i0)
            do 50 jj = 1, jx-4, NGS
               kk = min0( jx-jj-3, NGS )
               call vnorm( g, kk )
               do 50 j = 1, kk
                  i0      = i0 + 2
                  Z(i0)   = A(8,i)*(Z(j0  ) - Z(j0+4)) + Z(j0+2)          &
                          + A(9,i)*(Z(j0+1) - Z(j0+3)) + C(3,i)*g(j)
                  Z(i0+1) = two*Z(j0+2) - Z(i0)
                  j0      = j0 + 1
  50        continue
            i0      = i0 + 2
            Z(i0+1) = A(4,i)*Z(j0+3) + A(5,i)*Z(j0+2) + A(6,i)*Z(j0+1)    &
                                    + A(7,i)*Z(j0) + C(2,i)*ge(ii-1)
            Z(i0)   = two*Z(j0+2) - Z(i0+1)
            i0      = i0 + 2
            Z(i0+1) = A(1,i)*Z(j0+3) + A(2,i)*Z(j0+2) + A(3,i)*Z(j0+1)    &
                                    + C(1,i)*ge(ii)
            Z(i0)   = two*Z(j0+3) - Z(i0+1)
            jx      = 2*jx
  60     continue
      endif
!					all done, update timer
      ts = ts + (second() - tt)

      return
      end
