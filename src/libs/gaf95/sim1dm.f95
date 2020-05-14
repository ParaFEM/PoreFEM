!  *********************************************************************
!  *                                                                   *
!  *                         subroutine sim1dm                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.0
!  Written by Gordon A. Fenton, Dalhousie University, Mar 12, 2007
!  Latest Update: Mar 12, 2007
!
!  PURPOSE  simulates random processes for MSIM1M
!
!  DESCRIPTION
!  This routine simulates m random processes, which are
!  possibly intercorrelated. Individual fields are generated as standard
!  Gaussian random fields using the 1-D Local Average Subdivision (LAS)
!  algorithm. These fields are then combined to produce correlated standard
!  Gaussian fields using the Covariance Matrix Decomposition approach.
!  Finally the individual fields are transformed so that they have the
!  desired marginal distributions. These transformations are as follows;
!
!	P(x) = mean + sd*G(x)		       if normally distributed
!
!	P(x) = exp{ log-mean + log-sd*G(x) }   if lognormally distributed
!
!	P(x) = a + 0.5*(b-a)*[ 1 + tanh((m + s*G(x))/2*pi) ]
!						if bounded
!
!  where P(x) is the desired 1-D random process,
!  G(x) is one of the standard correlated Gaussian fields, and (mean,sd)
!  and (a,b,m,s) are parameters of the distributions.
!
!  If the process is deterministic, the entire process is simply set to its
!  mean.
!
!  ARGUMENTS
!
!   iterm	unit number to which error and warning messages are sent.
!		(input)
!
!	m	integer containing the number of random processes produced
!		in each realization. (input)
!
!      xl	real value containing the physical length of the process.
!		(input)
!
!     nxe	integer containing the number of cells that each process
!		is subdivided into. (input)
!
!     thx	real value containing the desired scale of fluctuation of
!		the processes. (input)
!
!  varfnc	character string of length 6 containing the name of the
!		1-D covariance function to use in the simulation. Possible
!		values of varfnc are;
!
!			dlace1	- damped oscillatory noise process
!			dlafr1	- fractional Gaussian noise process
!			dlavx1	- Markov process
!			dlsmp1	- simple polynomial decaying covariance
!			dlspx1	- Gaussian covariance function
!		(input)
!
!   kseed	integer containing the seed used to initialize the pseudo-
!		random number generator. (input)
!
!   lxfld	logical flag which is true if more than one of the processes
!		are cross-correlated. That is, if all processes are
!		independent, then lxfld is false. (input)
!
!	a	real array of size at least 7 x m which contains
!		the parameters of each of the i = 1, 2, ..., m processes.
!		Notably,
!		   a(1,i) = mean,
!		   a(2,i) = standard deviation,
!		   a(3,i) = distribution type;
!			  = 0.0 if process is deterministic (at mean value)
!			  = 1.0 if process is normally distributed
!			  = 2.0 if process is lognormally distributed (logn)
!			  = 3.0 if process is bounded
!		   a(4,i) = lower bound (bounded), or mean of log-process(logn)
!		   a(5,i) = upper bound (bounded), or sd of log-process (logn)
!		   a(6,i) = m parameter (if bounded)
!		   a(7,i) = s parameter (if bounded)
!		If process i (i = 1, 2, ..., m) is bounded, then a(1,i) and
!		a(2,i) are ignored and the parameters a(4,i) through a(7,i)
!		completely describe the distribution. (input)
!
!	R	real array of size at least m x m which contains
!		the correlation matrix between the m
!		random processes. For example R(i,j) will be the
!		correlation coefficient, assumed to act at all points
!		in the field, acting between process i and process j.
!		(input)
!
!       G	real array of size at least (3*nxe/2) x m which, on output,
!		will contain the m process realizations, each of length nxe.
!		The extra nxe/2, in each of the m columns in the array, is
!		used as workspace. (output)
!
!     MXM	integer containing the maximum allowable number of random
!		processes. (input)
!
!     MXP	integer containing the maximum allowable number of cells
!		in each process. (input)
!
! lsquare	logical flag which is true if the cells in each process are
!		actually two-dimensional squares. What this does is to
!		provide for local averaging over the square (assuming a
!		separable variance function) rather than just along the
!		1-D process. (input)
!
!   debug	logical flag which is true if debug info is to be written
!		to unit iterm. (input)
!
!  REVISION HISTORY:
!-------------------------------------------------------------------------
      subroutine sim1dm(iterm,m,xl,nxe,thx,varfnc,kseed,lxfld,            &
                       a,R,G,MXM,MXQ,lsquare,debug)
      real xl, thx, a(7,*), R(MXM,*), G(MXQ,*)
      integer m, nxe, kseed, MXM, MXQ
      logical lxfld, lsquare, debug
      character*(*) varfnc
      real*4 rerr
      real*8 dvar, dpb, dthx, dthy, dthz, ddx
      real*8 dlace1, dlafr1, dlavx1, dlsmp1, dlspx1
      real*8 dmin1, dble
      save ienter
      external dlace1, dlafr1, dlavx1, dlsmp1, dlspx1
      common/dparam/ dvar, dpb, dthx, dthy, dthz
      data zero/0.0/, half/0.5/, one/1.0/, onept5/1.5/, twopt5/2.5/, d4/4.0/
      data tolr/0.001/
      data twopi/6.2831853071795864769/
      data ienter/0/

   1  format(a,a)

!-------------------------------------- initialize -------------------------
!					compute required field size (once)
      if( debug ) then
         write(iterm,1)'sim1dm: just entered...'
         call flush(iterm)
      endif
      ienter = ienter + 1
      if( ienter .eq. 1 ) then
         if( lxfld ) then
            if( debug ) then
               write(iterm,1)'sim1dm: Decomposing R...'
               call flush(iterm)
            endif
            call chol2(R,MXM,m,rerr)		! decompose correlation matrix
            if( rerr .gt. tolr ) then
              call print1(istat,                                     &
		   'Warning: relative error on Cholesky decomposition of R is %f%n'&
		   ,rerr)
            endif
            if( debug ) then
               write(iterm,1)'sim1dm: R decomposed...'
               call flush(iterm)
            endif
         endif

         if( thx .gt. zero ) then
            ddx = dble(xl)/dble(nxe)
!						set covariance fnc parameters
            dthx = thx
            if( lsquare ) then
               dvar = -1.d0
               if(     varfnc .eq. 'dlace1' ) then
                  dpb = dlace1(ddx)
               elseif( varfnc .eq. 'dlafr1' ) then
                  dpb = dlafr1(ddx)
               elseif( varfnc .eq. 'dlavx1' ) then
                  dpb = dlavx1(ddx)
               elseif( varfnc .eq. 'dlsmp1' ) then
                  dpb = dlsmp1(ddx)
               elseif( varfnc .eq. 'dlspx1' ) then
                  dpb = dlspx1(ddx)
               endif
               dvar = dpb
            else
               dvar = 1.d0
            endif

            if( varfnc .eq. 'dlafr1' ) dpb = ddx

!						initialize LAS1G
            if( debug ) then
               write(iterm,1)'sim1dm: Initializing LAS1G...'
               call flush(iterm)
            endif
            if( varfnc .eq. 'dlace1' ) then
               call las1g(G,nxe,xl,dlace1,3,kseed,-1,iterm)
            elseif( varfnc .eq. 'dlafr1' ) then
               call las1g(G,nxe,xl,dlafr1,3,kseed,-1,iterm)
            elseif( varfnc .eq. 'dlavx1' ) then
               call las1g(G,nxe,xl,dlavx1,3,kseed,-1,iterm)
            elseif( varfnc .eq. 'dlsmp1' ) then
               call las1g(G,nxe,xl,dlsmp1,3,kseed,-1,iterm)
            elseif( varfnc .eq. 'dlspx1' ) then
               call las1g(G,nxe,xl,dlspx1,3,kseed,-1,iterm)
            else
               write(iterm,1)                                             &
                 '*** Error: unknown variance function name: ', varfnc
               stop
            endif
            if( debug ) then
               write(iterm,1)'sim1dm: LAS1G initialized...'
               call flush(iterm)
            endif
         endif
      endif
!					now produce m standard normal fields
      do 30 j = 1, m
         if( a(3,j) .lt. half ) then	! deterministic process
            do 10 i = 1, nxe
               G(i,j) = zero
  10        continue
         elseif( thx .eq. zero ) then	! white noise process
            do 20 i = 1, nxe
               G(i,j) = gausv(one)
  20        continue
         else
            if( debug ) then
               call print1(iterm,'sim1dm: Generating field #%i%n',rerr)
               call flush(iterm)
            endif
            if( varfnc .eq. 'dlace1' ) then
               call las1g(G(1,j),nxe,xl,dlace1,3,kseed,0,iterm)
            elseif( varfnc .eq. 'dlafr1' ) then
               call las1g(G(1,j),nxe,xl,dlace1,3,kseed,0,iterm)
            elseif( varfnc .eq. 'dlavx1' ) then
               call las1g(G(1,j),nxe,xl,dlavx1,3,kseed,0,iterm)
            elseif( varfnc .eq. 'dlsmp1' ) then
               call las1g(G(1,j),nxe,xl,dlsmp1,3,kseed,0,iterm)
            elseif( varfnc .eq. 'dlspx1' ) then
               call las1g(G(1,j),nxe,xl,dlspx1,3,kseed,0,iterm)
            endif
            if( debug ) then
               call print1(iterm,'sim1dm: Field #%i generated.%n',rerr)
               call flush(iterm)
            endif
         endif
  30  continue
!					combine for correlated fields...
      if( debug ) then
         write(iterm,1)'sim1dm: Correlating fields...'
         call flush(iterm)
      endif
      if( lxfld ) then
         do 60 i = 1, nxe
            do 50 j = m, 2, -1
               G(i,j) = R(j,j)*G(i,j)
               do 40 k = 1, j-1
                  G(i,j) = G(i,j) + R(k,j)*G(i,k)
  40           continue
  50        continue
  60     continue
      endif
      if( debug ) then
         write(iterm,1)'sim1dm: Fields correlated...'
         call flush(iterm)
      endif
!					convert to final processes
      do 110 j = 1, m
         if( debug ) then
            call print1(iterm,'sim1dm: Converting field #%i.%n',rerr)
            call flush(iterm)
         endif
         if( a(3,j) .lt. half ) then			! deterministic
            do 70 i = 1, nxe
               G(i,j) = a(1,j)
  70        continue
         elseif( a(3,j) .lt. onept5 ) then			! normal
            do 80 i = 1, nxe
               G(i,j) = a(1,j) + a(2,j)*G(i,j)
  80        continue
         elseif( a(3,j) .lt. twopt5 ) then			! lognormal
            do 90 i = 1, nxe
               G(i,j) = exp( a(4,j) + a(5,j)*G(i,j) )
  90        continue
         else						! bounded
            do 100 i = 1, nxe
			 G(i,j) = a(4,j) + half*(a(5,j)-a(4,j))                     &
                  *( one + tanh((a(6,j)+a(7,j)*G(i,j))/twopi) )
 100        continue
         endif
         if( debug ) then
            call print1(iterm,'sim1dm: Field #%i converted.%n',rerr)
            call flush(iterm)
         endif
 110  continue

      if( debug ) then
         write(iterm,1)'sim1dm: all done, just leaving...'
         call flush(iterm)
      endif

      return
      end
