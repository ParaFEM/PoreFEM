!  *********************************************************************
!  *                                                                   *
!  *                        subroutine las2g                           *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 3.42
!  Written by Gordon A. Fenton, TUNS, July 19, 1992
!  Latest Update: Jun 9, 1999
!
!  PURPOSE  produces a 2-D quadrant symmetric stationary Gaussian random field
!           using the Local Average Subdivision algorithm.
!
!  This routine creates a zero-mean realization of a 2-D random process given
!  its variance function (as defined by E.H. Vanmarcke in "Random Fields:
!  Analysis and Synthesis", MIT Press, 1984). Each discrete value generated
!  herein represents the local average of a realization of the process over
!  the area Dx x Dy, where (Dx,Dy) is the grid spacing of the desired field.
!  The construction of the realization proceeds recursively as follows;
!
!    1) generate a low resolution field of size k1 x k2. If (N1 x N2) is
!       the desired field resolution, then k1 and k2 are determined so
!       that N1 = k1*2**m, N2 = k2*2**m, where 2**m is a common factor and
!       k1*k2 <= MXK. Generally k1 and k2 are maximized where possible.
!       This is refered to subsequently as the Stage 0 generation.
!
!    2) subdivide the domain m times by dividing each cell into 4 equal
!       parts (2 x 2). In each subdivision, new random values are generated
!       for each new cell. The parent cells of the previous stage are used
!       to obtain best linear estimates of the mean of each new cell so that
!       the spatial correlation is approximated. Also upwards averaging is
!       preserved (ie the average of the 4 new values is the same as the
!       parent cell value). Only parent cells in a neighbourhood of 3 x 3
!       are considered (thus the approximation to the spatial correlation).
!
!  The linear estimation of the mean is accomplished by using the covariance
!  between local averages over each cell, consistent with the goal of
!  producing a local average field. Note that this conditioning process
!  implies that the construction of cells near the edge of the boundary will
!  require the use of values which are, strictly speaking, outside the
!  boundary of the field in question. This is handled by using special
!  reduced neighbourhoods along the boundaries (equivalent to saying that
!  what goes on beyond the boundary has no effect on the process within the
!  boundary).
!
!  Note that this routine sets up a number of parameters on the
!  first call and thus the time required to produce the first realization
!  is substantially greater than on subsequent calls (see INIT). For
!  more information on local average processes, see
!
!   1) G.A. Fenton, "Simulation and Analysis of Random Fields", Ph.D. thesis,
!      Dept. of Civil Eng. and Op. Research, Princeton University, Princeton,
!      New Jersey, 1990.
!   2) G.A. Fenton and E.H. Vanmarcke "Simulation of Random Fields
!      via Local Average Subdivision", ASCE Journal of Engineering Mechanics,
!      Vol. 116, No. 8, August 1990.
!   3) E. H. Vanmarcke, Random Fields: Analysis and Synthesis, MIT Press,
!      1984
!   4) G.A. Fenton, Error evaluation of three random field generators,
!      ASCE Journal of Engineering Mechanics (to appear), 1993.
!
!  Arguments to this routine are as follows;
!
!      Z  real array of size at least N1 x (5*N2/4) which on output will
!         contain the realization of the 2-D process in the first N1 x N2
!         locations. When dimensioned as Z(N1,5*N2/4) in the calling routine
!         and indexed as Z(i,j), Z(1,1) is the lower left cell, Z(2,1) is
!         the next cell in the X direction (to the right), etc.,
!         while Z(1,2) is the next cell in the Y direction (upwards), etc.
!         The extra space is required for workspace (specifically to store
!         the previous stage in the subdivision). (output)
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
!  dvarfn  external real*8 function which returns the variance of the random
!         process averaged over a given area. dvarfn is referenced as follows
!
!                var = dvarfn( V1, V2 )
!
!         where (V1,V2) are the side dimensions of the rectangular averaging
!         domain. Any other parameters to the function must be passed by
!         common block from the calling routine. Note that the variance of
!         the random process averaged over the area (V1 x V2) is the product
!         of the point variance and the traditionally defined "variance"
!         function, as discussed by Vanmarcke (pg 186).
!
!  KSEED  integer seed to be used to initialize the pseudo-random number
!         generator. The generator is only initialized on the first call
!         to this routine or when abs(INIT) = 1.
!         If KSEED = 0, then a random seed will be used (based on the
!         process ID of the current program invocation - see iseed.f)
!         On output, KSEED is set to the value of the actual seed used.
!
!   INIT  integer flag which must be 1 when parameters of the process
!         are to be calculated or recalculated. If multiple realizations
!         of the same process are desired, then subsequent calls should
!         use INIT equal to 0 and M less than or equal to the value used 
!         initially. Note that if INIT = -1 is used, then only the
!         initialization is performed.
!
!   IOUT  unit number to which error and warning messages are to be logged.
!         (input)
!  -------------------------------------------------------------------------
!  PARAMETERS:
!   MXM   represents the maximum value that m can have in the decomposition
!         N = k * 2**m.
!
!   MXK   represents the maximum number of cells (k1 x k2) in the initial
!         field, k1*k2 <= MXK. If the value of MXK is changed here, it must
!         also be changed in LAS2I.
!
!   NGS   represents the maximum number of random Gaussian variates that can
!         be produced on each call to VNORM. NGS should be 3*[9*2**(MXM-1)],
!         but not less than MXK.
!
!  Notes: 
!    1) Simulation timing is available through common block LASTYM where
!       TI is the time required to initialize the parameters of the
!       process and TS is the cumulative simulation time. These timings
!       are obtained through the function SECOND which returns elapsed
!       user time in seconds since the start of the program.
!    2) In 2 and higher dimensions, the LAS method yeilds a slight pattern
!       in the point variance field. This is due to the neighborhood
!       approximations. This error lessens for larger scales of fluctuation
!       and smaller numbers of subdivisions.
!    3) The parameter `tol' is the maximum relative error allowed on the
!       Cholesky decomposition of covariance matrices before a warning
!       message is emitted. The relative error is estimated by computing
!       the lower-right most element of L*L' and comparing to the original
!       element of A (where L*L' = A is the decomposition). This give some
!       measure of the roundoff errors accumulated in the computation of L.
!
!  Requires:
!    1) from libGAFsim:	ISEED, LAS2I, DCVIT2, DCVMT2, DCHOL2, CORN2D, SIDE2D,
!			INTR2D, DCVAA2, DCVAB2, DSIFA, DSISL, DAXPY, DSWAP,
!			IDAMAX, DDOT, VNORM, RANDF, SECOND
!    2) user defined external variance function (see DVARFN)
!
!  REVISION HISTORY:
!  3.41	export saved parameters via /las2gb/ (previously just static)
!  3.42	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!  ==========================================================================
      subroutine las2g2( Z, N1, N2, xl2, YL2, dvarfn, kseed, init, iout )
      parameter( MXM = 9, MXK = 256, NGS = 6912 )
      real Z(*)
      real C0((MXK*(MXK + 1))/2), U(NGS)
      real CT(6,2), CC(6,4,MXM), CS(6,4,MXM), CI(6,MXM)
      real AT(3,3,2), AC(4,3,4,MXM), AS(6,3,4,MXM), AI(9,3,MXM)
      external dvarfn
      common/las2gb2/C0,CT,CC,CS,CI,AT,AC,AS,AI,M,K1,K2,KK,NN
      common/LASTYM/ ti, ts
      data ifirst2/1/, zero/0.0/, four/4.0/
      data tol/1.e-3/
!-------------------------------------- initialize LAS generator -------------

      if( ifirst2 .eq. 1 .or. iabs(init) .eq. 1 ) then
!						start timer
         ti = second()
         call las2i2( dvarfn, N1, N2, xl2, YL2, kseed, MXM,C0, CT, CC, CS,   &
		            CI, AT, AC, AS, AI, M,k1, k2, kk, iout, tol )
         NN = N1*N2
         ifirst2 = 0
!						all done, set timers
         ts   = zero
         ti   = second() - ti
         if( init .eq. -1 ) return
      endif
!-------------------------------------- generate realization -----------------
      tt = second()
      if( mod(M,2) .eq. 0 ) then
         iz = 0
         jz = NN
      else
         iz = NN
         jz = 0
      endif
!					generate stage 0 field
      call vnorm2( U, kk )
      L = 1
      do 20 i = 1, kk
         Z(iz+i) = C0(L)*U(1)
         do 10 j = 2, i
            Z(iz+i) = Z(iz+i) + C0(L+j-1)*U(j)
  10     continue
         L = L + i
  20  continue
!					generate stage 1, 2, ..., M sub-fields
      jx = k1
      jy = k2
      nm = 1
      if( (k1.eq.1 .or. k2.eq.1) .and. M .gt. 0 ) then
         jq = max0(k1,k2)
         call vnorm2( U, 3*jq )
         it = jz
         jz = iz
         iz = it
         ix = 2*k1
         iy = 2
         if( k1 .eq. 1 ) iy = 4
         i0 = iz + 1
         i1 = i0 + ix
         j0 = jz + 1
         Z(i0)   = AT(1,1,1)*Z(j0) + AT(2,1,1)*Z(j0+1) + CT(1,1)*U(1)
         Z(i0+1) = AT(1,2,1)*Z(j0) + AT(2,2,1)*Z(j0+1) + CT(2,1)*U(1) + CT(3,1)*U(2)
         Z(i1)   = AT(1,3,1)*Z(j0) + AT(2,3,1)*Z(j0+1) + CT(4,1)*U(1) + CT(5,1)*U(2) + CT(6,1)*U(3)
         Z(i1+1) = four*Z(j0) - Z(i0) - Z(i0+1) - Z(i1)
         L = 4
         do 30 js = 1, jq-2
            i0 = i0 + iy
            i1 = i1 + iy
            Z(i0)   = AT(1,1,2)*Z(j0  ) + AT(2,1,2)*Z(j0+1) + AT(3,1,2)*Z(j0+2) + CT(1,2)*U(L)
            Z(i0+1) = AT(1,2,2)*Z(j0  ) + AT(2,2,2)*Z(j0+1) + AT(3,2,2)*Z(j0+2) + CT(2,2)*U(L) + CT(3,2)*U(L+1)
            Z(i1)   = AT(1,3,2)*Z(j0  ) + AT(2,3,2)*Z(j0+1) + AT(3,3,2)*Z(j0+2) + CT(4,2)*U(L) + CT(5,2)*U(L+1) + CT(6,2)*U(L+2)
            Z(i1+1) = four*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)
            j0 = j0 + 1
            L  = L  + 3
  30     continue
         i0 = i0 + iy
         i1 = i1 + iy
         Z(i0)   = AT(2,1,1)*Z(j0) + AT(1,1,1)*Z(j0+1) + CT(1,1)*U(1)
         Z(i0+1) = AT(2,2,1)*Z(j0) + AT(1,2,1)*Z(j0+1) + CT(2,1)*U(1) + CT(3,1)*U(2)
         Z(i1)   = AT(2,3,1)*Z(j0) + AT(1,3,1)*Z(j0+1) + CT(4,1)*U(1) + CT(5,1)*U(2) + CT(6,1)*U(3)
         Z(i1+1) = four*Z(j0) - Z(i0) - Z(i0+1) - Z(i1)
         jx = 2*k1
         jy = 2*k2
         nm = 2
      endif

      do 80 i = nm, M
!						swap current and prev fields
         it = jz
         jz = iz
         iz = it
!						new field dimensions
         ix = 2*jx
         iy = 2*jy
!						pointers into Z
         j0 = jz + 1
         j1 = j0 + jx
         i0 = iz + 1
         i1 = i0 + ix
         call vnorm2( U, 3*jx )
!								corner #1
         Z(i0)   = AC(1,1,1,i)*Z(j0) + AC(2,1,1,i)*Z(j0+1)                &
		           + AC(3,1,1,i)*Z(j1) + AC(4,1,1,i)*Z(j1+1)              &
				   + CC(1,  1,i)*U(1)
!
         Z(i0+1) = AC(1,2,1,i)*Z(j0) + AC(2,2,1,i)*Z(j0+1)                &
		           + AC(3,2,1,i)*Z(j1) + AC(4,2,1,i)*Z(j1+1)              &
				   + CC(2,  1,i)*U(1) +  CC(3,  1,i)*U(2)
         Z(i1)   = AC(1,3,1,i)*Z(j0) + AC(2,3,1,i)*Z(j0+1)                &
		           + AC(3,3,1,i)*Z(j1) + AC(4,3,1,i)*Z(j1+1)              &
				   + CC(4,  1,i)*U(1) +  CC(5,  1,i)*U(2)                 &
				   + CC(6,  1,i)*U(3)
         Z(i1+1) = four*Z(j0) - Z(i0) - Z(i0+1) - Z(i1)
!								side #1
         L = 4
         do 40 js = 1, jx-2
            i0 = i0 + 2
            i1 = i1 + 2
            Z(i0)   = AS(1,1,1,i)*Z(j0  ) + AS(2,1,1,i)*Z(j0+1)           &
			          + AS(3,1,1,i)*Z(j0+2) + AS(4,1,1,i)*Z(j1  )         &
					  + AS(5,1,1,i)*Z(j1+1) + AS(6,1,1,i)*Z(j1+2)         &
					  + CS(1,  1,i)*U(L)
            Z(i0+1) = AS(1,2,1,i)*Z(j0  ) + AS(2,2,1,i)*Z(j0+1)           &
			          + AS(3,2,1,i)*Z(j0+2) + AS(4,2,1,i)*Z(j1  )         &
					  + AS(5,2,1,i)*Z(j1+1) + AS(6,2,1,i)*Z(j1+2)         &
					  + CS(2,  1,i)*U(L)    + CS(3,  1,i)*U(L+1)
            Z(i1)   = AS(1,3,1,i)*Z(j0  ) + AS(2,3,1,i)*Z(j0+1)           &
			          + AS(3,3,1,i)*Z(j0+2) + AS(4,3,1,i)*Z(j1  )         &
					  + AS(5,3,1,i)*Z(j1+1) + AS(6,3,1,i)*Z(j1+2)         &
					  + CS(4,  1,i)*U(L)    + CS(5,  1,i)*U(L+1)          &
					  + CS(6,  1,i)*U(L+2)
            Z(i1+1) = four*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)
            j0 = j0 + 1
            j1 = j1 + 1
            L  = L  + 3
  40     continue
!								corner #2
         i0 = i0 + 2
         i1 = i1 + 2
         Z(i0)   = AC(1,1,2,i)*Z(j0) + AC(2,1,2,i)*Z(j0+1)                &
		           + AC(3,1,2,i)*Z(j1) + AC(4,1,2,i)*Z(j1+1)              &
				   + CC(1,  2,i)*U(L)
         Z(i0+1) = AC(1,2,2,i)*Z(j0) + AC(2,2,2,i)*Z(j0+1)                &
		           + AC(3,2,2,i)*Z(j1) + AC(4,2,2,i)*Z(j1+1)              &
				   + CC(2,  2,i)*U(L)  + CC(3,  2,i)*U(L+1)
         Z(i1)   = AC(1,3,2,i)*Z(j0) + AC(2,3,2,i)*Z(j0+1)                &
		           + AC(3,3,2,i)*Z(j1) + AC(4,3,2,i)*Z(j1+1)              &
				   + CC(4,  2,i)*U(L)  + CC(5,  2,i)*U(L+1)               &
				   + CC(6,  2,i)*U(L+2)
         Z(i1+1) = four*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)

         j0 = jz + 1
         do 60 ks = 1, jy-2
            j1 = j0 + jx
            j2 = j1 + jx
            i0 = i1 + 2
            i1 = i0 + ix
            call vnorm2( U, 3*jx )
!								side #2
            Z(i0)   = AS(1,1,2,i)*Z(j0) + AS(2,1,2,i)*Z(j0+1)             &
			          + AS(3,1,2,i)*Z(j1) + AS(4,1,2,i)*Z(j1+1)           &
					  + AS(5,1,2,i)*Z(j2) + AS(6,1,2,i)*Z(j2+1)           &
					  + CS(1,  2,i)*U(1)
            Z(i0+1) = AS(1,2,2,i)*Z(j0) + AS(2,2,2,i)*Z(j0+1)             &
			          + AS(3,2,2,i)*Z(j1) + AS(4,2,2,i)*Z(j1+1)           &
					  + AS(5,2,2,i)*Z(j2) + AS(6,2,2,i)*Z(j2+1)           &
					  + CS(2,  2,i)*U(1)  + CS(3,  2,i)*U(2)
            Z(i1)   = AS(1,3,2,i)*Z(j0) + AS(2,3,2,i)*Z(j0+1)             &
			          + AS(3,3,2,i)*Z(j1) + AS(4,3,2,i)*Z(j1+1)           &
					  + AS(5,3,2,i)*Z(j2) + AS(6,3,2,i)*Z(j2+1)           &
					  + CS(4,  2,i)*U(1)  + CS(5,  2,i)*U(2)              &
					  + CS(6,  2,i)*U(3)
            Z(i1+1) = four*Z(j1) - Z(i0) - Z(i0+1) - Z(i1)
!								interior
            L = 4
            do 50 js = 1, jx-2
               i0 = i0 + 2
               it = i0 + 1
               i1 = i1 + 2
               Z(i0)=AI(1,1,i)*Z(j0)+AI(2,1,i)*Z(j0+1)+AI(3,1,i)*Z(j0+2)  &
			         +AI(4,1,i)*Z(j1)+AI(5,1,i)*Z(j1+1)+AI(6,1,i)*Z(j1+2) &
					 +AI(7,1,i)*Z(j2)+AI(8,1,i)*Z(j2+1)+AI(9,1,i)*Z(j2+2) &
					 +CI(1,  i)*U(L)
               Z(it)=AI(1,2,i)*Z(j0)+AI(2,2,i)*Z(j0+1)+AI(3,2,i)*Z(j0+2)  &
			         +AI(4,2,i)*Z(j1)+AI(5,2,i)*Z(j1+1)+AI(6,2,i)*Z(j1+2) &
					 +AI(7,2,i)*Z(j2)+AI(8,2,i)*Z(j2+1)+AI(9,2,i)*Z(j2+2) &
					 +CI(2,  i)*U(L) +CI(3,  i)*U(L+1)
               Z(i1)=AI(1,3,i)*Z(j0)+AI(2,3,i)*Z(j0+1)+AI(3,3,i)*Z(j0+2)  &
			         +AI(4,3,i)*Z(j1)+AI(5,3,i)*Z(j1+1)+AI(6,3,i)*Z(j1+2) &
					 +AI(7,3,i)*Z(j2)+AI(8,3,i)*Z(j2+1)+AI(9,3,i)*Z(j2+2) &
					 +CI(4,  i)*U(L) +CI(5,  i)*U(L+1) +CI(6,  i)*U(L+2)
               Z(i1+1) = four*Z(j1+1) - Z(i0) - Z(it) - Z(i1)
               j0 = j0 + 1
               j1 = j1 + 1
               j2 = j2 + 1
               L  = L  + 3
  50        continue
!								side #3
            i0 = i0 + 2
            i1 = i1 + 2
            Z(i0)   = AS(1,1,3,i)*Z(j0) + AS(2,1,3,i)*Z(j0+1)             &
			          + AS(3,1,3,i)*Z(j1) + AS(4,1,3,i)*Z(j1+1)           &
					  + AS(5,1,3,i)*Z(j2) + AS(6,1,3,i)*Z(j2+1)           &
					  + CS(1,  3,i)*U(L)
            Z(i0+1) = AS(1,2,3,i)*Z(j0) + AS(2,2,3,i)*Z(j0+1)             &
			          + AS(3,2,3,i)*Z(j1) + AS(4,2,3,i)*Z(j1+1)           &
					  + AS(5,2,3,i)*Z(j2) + AS(6,2,3,i)*Z(j2+1)           &
					  + CS(2,  3,i)*U(L)  + CS(3,  3,i)*U(L+1)
            Z(i1)   = AS(1,3,3,i)*Z(j0) + AS(2,3,3,i)*Z(j0+1)             &
			          + AS(3,3,3,i)*Z(j1) + AS(4,3,3,i)*Z(j1+1)           &
					  + AS(5,3,3,i)*Z(j2) + AS(6,3,3,i)*Z(j2+1)           &
					  + CS(4,  3,i)*U(L)  + CS(5,  3,i)*U(L+1)            &
					  + CS(6,  3,i)*U(L+2)
            Z(i1+1) = four*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
            j0 = j0 + 2
  60     continue

         j1 = j0 + jx
         i0 = i1 + 2
         i1 = i0 + ix
         call vnorm2( U, 3*jx )
!								corner #3
         Z(i0)   = AC(1,1,3,i)*Z(j0) + AC(2,1,3,i)*Z(j0+1)                 &
		           + AC(3,1,3,i)*Z(j1) + AC(4,1,3,i)*Z(j1+1)               &
				   + CC(1,  3,i)*U(1)
		 Z(i0+1) = AC(1,2,3,i)*Z(j0) + AC(2,2,3,i)*Z(j0+1)                 &
		           + AC(3,2,3,i)*Z(j1) + AC(4,2,3,i)*Z(j1+1)               &
				   + CC(2,  3,i)*U(1)  + CC(3,  3,i)*U(2)
         Z(i1)   = AC(1,3,3,i)*Z(j0) + AC(2,3,3,i)*Z(j0+1)                 &
		           + AC(3,3,3,i)*Z(j1) + AC(4,3,3,i)*Z(j1+1)               &
				   + CC(4,  3,i)*U(1)  + CC(5,  3,i)*U(2)                  &
				   + CC(6,  3,i)*U(3)
         Z(i1+1) = four*Z(j1) - Z(i0) - Z(i0+1) - Z(i1)
!								side #4
         L = 4
         do 70 js = 1, jx-2
            i0 = i0 + 2
            i1 = i1 + 2
            Z(i0)   = AS(1,1,4,i)*Z(j0  ) + AS(2,1,4,i)*Z(j0+1)            &
			         + AS(3,1,4,i)*Z(j0+2) + AS(4,1,4,i)*Z(j1  )           &
					 + AS(5,1,4,i)*Z(j1+1) + AS(6,1,4,i)*Z(j1+2)           &
					 + CS(1,  4,i)*U(L)
            Z(i0+1) = AS(1,2,4,i)*Z(j0  ) + AS(2,2,4,i)*Z(j0+1)            &
			          + AS(3,2,4,i)*Z(j0+2) + AS(4,2,4,i)*Z(j1  )          &
					  + AS(5,2,4,i)*Z(j1+1) + AS(6,2,4,i)*Z(j1+2)          &
					  + CS(2,  4,i)*U(L)    + CS(3,  4,i)*U(L+1)
            Z(i1)   = AS(1,3,4,i)*Z(j0  ) + AS(2,3,4,i)*Z(j0+1)            &
			          + AS(3,3,4,i)*Z(j0+2) + AS(4,3,4,i)*Z(j1  )          &
					  + AS(5,3,4,i)*Z(j1+1) + AS(6,3,4,i)*Z(j1+2)          &
					  + CS(4,  4,i)*U(L)    + CS(5,  4,i)*U(L+1)           &
					  + CS(6,  4,i)*U(L+2)
            Z(i1+1) = four*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
            j0 = j0 + 1
            j1 = j1 + 1
            L  = L  + 3
  70     continue
!								corner #4
         i0 = i0 + 2
         i1 = i1 + 2
         Z(i0)   = AC(1,1,4,i)*Z(j0) + AC(2,1,4,i)*Z(j0+1)                 &
		           + AC(3,1,4,i)*Z(j1) + AC(4,1,4,i)*Z(j1+1)               &
				   + CC(1,  4,i)*U(L)
         Z(i0+1) = AC(1,2,4,i)*Z(j0) + AC(2,2,4,i)*Z(j0+1)                 &
		           + AC(3,2,4,i)*Z(j1) + AC(4,2,4,i)*Z(j1+1)               &
				   + CC(2,  4,i)*U(L)  + CC(3,  4,i)*U(L+1)
         Z(i1)   = AC(1,3,4,i)*Z(j0) + AC(2,3,4,i)*Z(j0+1)                 &
		           + AC(3,3,4,i)*Z(j1) + AC(4,3,4,i)*Z(j1+1)               &
				   + CC(4,  4,i)*U(L)  + CC(5,  4,i)*U(L+1)                &
				   + CC(6,  4,i)*U(L+2)
         Z(i1+1) = four*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)

         jx = ix
         jy = iy
  80  continue
!						all done, compute elapsed time
      ts = ts + (second() - tt)

      return
      end
