!  *********************************************************************
!  *                                                                   *
!  *                         subroutine sim2sd                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.4
!  Written by Gordon A. Fenton, DalTech, Nov 17, 1999
!  Latest Update: Oct 18, 2002
!
!  PURPOSE  simulates soil property fields for RSLOPE2D
!
!  DESCRIPTION
!  This routine simulates up to 6 soil property fields (cohesion, friction
!  angle, dilation angle, unit weight, elastic modulus, and Poisson's ratio
!  which are possibly intercorrelated. Individual fields are generated as
!  standard Gaussian random fields using the 2-D Local Average Subdivision
!  (LAS) algorithm. These fields are then combined to produce correlated
!  standard Gaussian fields using the Covariance Matrix Decomposition
!  approach. Finally the individual fields are transformed so that they
!  have the desired marginal distributions. These transformations are as
!  follows;
!
!	P(x,y) = mean + sd*G(x,y)		if normally distributed
!
!	P(x,y) = exp{ log-mean + log-sd*G(x,y) }if lognormally distributed
!
!	P(x,y) = a + 0.5*(b-a)*[ 1 + tanh((m + s*G(x,y))/2*pi) ]
!						if bounded
!
!  where P(x,y) is the desired random field (one of the soil properties)
!  G(x,y) is one of the standard correlated Gaussian fields, and (mean,sd)
!  and (a,b,m,s) are parameters of the distributions.
!
!  If the property is deterministic, the entire field is simply set to the
!  mean.
!
!  ARGUMENTS
!
!   istat	unit number connected to a file to which the input data is to
!		echoed if echo is true. (input)
!
!   iterm	unit number connected to the screen. If verbos is true, then
!		error and warning messages are also sent to the screen. (input)
!
!  verbos	logical flag which is true if error, warning, and progress
!		messages are allowed to be sent to the screen. (input)
!
!    cfld	real array of size at least (nrfx x nrfy) which will contain
!		the (optionally) random (locally averaged) cohesion field.
!		(output)
!
!  phifld	real array of size at least (nrfx x nrfy) which will contain
!		the (optionally) random (locally averaged) friction angle
!		field. Note that phifld may optionally contain the tan(phi)
!		field, as specified in the data file read by readsd. (output)
!
!  psifld	real array of size at least (nrfx x nrfy) which will contain
!		the (optionally) random (locally averaged) dilation angle
!		field. (output)
!
!  gamfld	real array of size at least (nrfx x nrfy) which will contain
!		the (optionally) random (locally averaged) unit weight field.
!		(output)
!
!    efld	real array of size at least (nrfx x nrfy) which will contain
!		the (optionally) random (locally averaged) elastic modulus
!		field. (output)
!
!    vfld	real array of size at least (nrfx x nrfy) which will contain
!		the (optionally) random (locally averaged) poisson's ratio
!		field. (output)
!
!	c	real vector of length at least 7 containing the parameters of
!		the soil cohesion. Notably,
!		   c(1) = mean cohesion,
!		   c(2) = cohesion standard deviation,
!		   c(3) = cohesion distribution type;
!			  = 0.0 if cohesion is deterministic (at mean value)
!			  = 1.0 if cohesion is normally distributed
!			  = 2.0 if cohesion is lognormally distributed (logn)
!			  = 3.0 if cohesion is bounded
!		   c(4) = lower bound (bounded), or mean of log-cohesion (logn)
!		   c(5) = upper bound (bounded), or sd of log-cohesion (logn)
!		   c(6) = m parameter (if bounded)
!		   c(7) = s parameter (if bounded)
!		If c is bounded, then c(1) and c(2) are ignored and the
!		parameters c(4) through c(7) completely describe the
!		distribution. (input)
!
!     phi	real vector of length at least 7 containing the parameters of
!		soil friction angle. See `c' for what the various elements of
!		phi contain. (input)
!
!     psi	real vector of length at least 7 containing the parameters of
!		soil dilation angle. See `c' for what the various elements of
!		psi contain. (input)
!
!     gam	real vector of length at least 7 containing the parameters of
!		soil unit weight. See `c' for what the various elements of
!		gam contain. (input)
!
!	e	real vector of length at least 7 containing the parameters of
!		soil elastic modulus. See `c' for what the various elements of
!		e contain. (input)
!
!	v	real vector of length at least 7 containing the parameters of
!		soil Poisson ratio. See `c' for what the various elements of
!		v contain. (input)
!
!	R	real array of size at least 6 x 6 which, on output, will
!		contain the correlation matrix between the 6 (possibly)
!		random soil properties.  Indexing into R is as follows;
!		  1 = cohesion
!		  2 = friction angle
!		  3 = dilation angle
!		  4 = unit weight
!		  5 = elastic modulus
!		  6 = Poisson's ratio
!		(input)
!
!   lxfld	logical flag which is true if more than one soil property
!		are cross-correlated. That is, if all soil properties are
!		independent, then lxfld is false. (input)
!
!     thx	real value giving the x-direction scale of fluctuation
!		(or, at least, this is a parameter of the covariance function).
!		(input)
!
!     thy	real value giving the y-direction scale of fluctuation
!		(or, at least, this is another parameter of the covariance
!		function). (input)
!
!     nrfx	integer giving the number of elements describing the soil
!		mass in the x-direction (horizontally). (input)
!
!     nrfy	integer giving the number of elements describing the soil
!		mass in the y-direction (vertically). (input)
!
!      dx	real value giving the physical size of an element in the
!		x-direction. (input)
!
!      dy	real value giving the physical size of an element in the
!		y-direction. (input)
!
!  dmpfld	logical flag which is true if a random field is to be
!		plotted to *.fld. (input)
!
!    nfld	integer giving the realization number of the random field
!		which is to be plotted to *.fld. (input)
!
!    jfld	integer denoting which random field is to be plotted;
!		 = 1 for the cohesion (c) field,
!		 = 2 for the friction angle (phi) field,
!		 = 3 for the dilation angle (psi) field,
!		 = 4 for the unit weight (gam) field,
!		 = 5 for the elastic modulus (psi) field,
!		 = 6 for the poisson ratio (psi) field,
!		(input)
!
!    ifld	unit number connected to the file to which the random field
!		plot is to be written in the event that dmpfld is true.
!		(input)
!
!     job	character string containing the main title of the run.
!		(input)
!
!    sub1	character string, which on output will contain the first
!		subtitle for this run. (output)
!
!    sub2	character string, which on output will contain the second
!		subtitle for this run. (output)
!
!  varfnc	character string containing the name of the covariance
!		function controlling the random fields. Possible covariance
!		functions are
!		`dlavx2' - 2-D exponentially decaying (Markov) model
!			   requires X- and Y-direction scales of fluctuation
!		`dlafr2' - 2-D isotropic fractional Gaussian noise model
!			   requires (H,delta) as parameters. In this case,
!			   thx is H, and delta is the minimum element
!			   dimension.
!		`dlsep2' - 2-D separable (1D x 1D) Markov model
!			   requires X- and Y-direction scales of fluctuation
!		`dlsfr2' - 2-D separable fractional Gaussian noise model
!			   requires (H_x,H_y,delta) as parameters. In this
!			   case, thx is H_x, thy is H_y, and delta is the
!			   minimum element dimension.
!		`dlspx2' - 2-D separable Gaussian decaying model
!			   requires X- and Y-direction scales of fluctuation
!		(input)
!
!   kseed	integer giving the seed to be used to initialize the
!		pseudo-random number generator. Subsequent runs using
!		the same seed will result in the same sequence of random
!		numbers. (input)
!
!   debug	logical flag which is true if debug information is to be
!		sent to the *.stt file. (input)
!
!  ltanfi	logical flag which is true if the phifld array contains the
!		tan(phi) field, rather than the phi field, where phi is the
!		friction angle. (input)
!
!  REVISION HISTORY:
!  1.1	the LAS generator now integrates the covariance function directly,
!	rather than using the variance function, see lvarfn (May 26/00)
!  1.2	eliminated lvarfn from common /dparam/ (May 26/01)
!  1.21	fixed truncated documentation above (Mar 1/02)
!  1.3	added ltanfi flag (used for title to pltfld) (Jun 13/02)
!  1.31	output error message if covariance function is unknown (Sep 9/02)
!  1.4	added data definition for twopi! (Oct 18/02)
!-----------------------------------------------------------------------
      subroutine sim2sd1_2(istat,iterm,verbos,cfld,phifld,     &
                       c,phi,R,lxfld,thx,thy,       &
                       nrfx,nrfy,dx,dy,dmpfld,nfld,jfld,ifld,job,sub1,    &
                       sub2,varfnc,kseed,debug,ltanfi)

      real cfld(nrfx,*), phifld(nrfx,*)
      real c(*), phi(*),R(2,*), thx, thy
      integer nrfx, nrfy, nfld, ifld
      character*(*) job, sub1, sub2, varfnc
      logical verbos, dmpfld, debug, liid, shofld, lxfld, ltanfi
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dmin1, dble
      save XL, YL, ienter, liid
      external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      common/dparam/ dvar, dpb, dthx, dthy, dthz
      data zero/0.0/, half/0.5/, one/1.0/, onept5/1.5/, twopt5/2.5/
      data twopi/6.2831853071795864769/
      data ienter/0/

   1  format(a,a)

!-------------------------------------- initialize ---------------------
!					compute required field size (once)
      ienter = ienter + 1
      if( ienter .eq. 1 ) then
         liid   = ((thx .eq. zero) .and. (thy .eq. zero))
         XL     = float(nrfx)*dx
         YL     = float(nrfy)*dy

         if( debug ) write(istat,1)'SIM2SD: setting field parameters...'

         if( .not. liid ) then
!						set variance fnc parameters
            dvar = 1.d0
            dthx = thx
            dthy = thy

            if( varfnc .eq. 'dlafr2' .or. varfnc .eq. 'dlsfr2' ) then
               dpb = dx
               if( dy .lt. dx ) dpb = dy
            endif
!						initialize LAS2G
            if( varfnc .eq. 'dlavx2' ) then
               call las2g(cfld,nrfx,nrfy,XL,YL,dlavx2,kseed,-1,istat)
            elseif( varfnc .eq. 'dlsep2' ) then
               call las2g(cfld,nrfx,nrfy,XL,YL,dlsep2,kseed,-1,istat)
            elseif( varfnc .eq. 'dlspx2' ) then
               call las2g(cfld,nrfx,nrfy,XL,YL,dlspx2,kseed,-1,istat)
            elseif( varfnc .eq. 'dlafr2' ) then
               call las2g(cfld,nrfx,nrfy,XL,YL,dlafr2,kseed,-1,istat)
            elseif( varfnc .eq. 'dlsfr2' ) then
               call las2g(cfld,nrfx,nrfy,XL,YL,dlsfr2,kseed,-1,istat)
            else
               if( verbos ) then
                  write(iterm,1)                                          &
                 '*** Error: unknown variance function name: ', varfnc
               endif
               write(istat,1)                                             &
                 '*** Error: unknown variance function name: ', varfnc
               stop
            endif
         endif
      endif
!					are we going to plot anything?
      shofld = dmpfld .and. (nfld .eq. ienter)

!					now produce 5 standard normal fields


!					for cohesion...
      if( c(3) .lt. half ) then			! cohesion is deterministic
         do 10 j = 1, nrfy
         do 10 i = 1, nrfx
            cfld(i,j) = zero			! we'll add mean back on later
  10     continue
      elseif( liid ) then			! white noise field
         do 20 j = 1, nrfy
         do 20 i = 1, nrfx
            cfld(i,j) = gausv(one)
  20     continue
      else					! generate standard normal RF
         if( varfnc .eq. 'dlavx2' ) then
            call las2g(cfld,nrfx,nrfy,XL,YL,dlavx2,kseed,0,istat)
         elseif( varfnc .eq. 'dlsep2' ) then
            call las2g(cfld,nrfx,nrfy,XL,YL,dlsep2,kseed,0,istat)
         elseif( varfnc .eq. 'dlspx2' ) then
            call las2g(cfld,nrfx,nrfy,XL,YL,dlspx2,kseed,0,istat)
         elseif( varfnc .eq. 'dlafr2' ) then
            call las2g(cfld,nrfx,nrfy,XL,YL,dlafr2,kseed,0,istat)
         elseif( varfnc .eq. 'dlsfr2' ) then
            call las2g(cfld,nrfx,nrfy,XL,YL,dlsfr2,kseed,0,istat)
         endif
      endif


!					for friction angle...
      if( phi(3) .lt. half ) then		! friction is deterministic
         do 30 j = 1, nrfy
         do 30 i = 1, nrfx
            phifld(i,j) = zero			! we'll add mean back on later
  30     continue
      elseif( liid ) then			! white noise field
         do 40 j = 1, nrfy
         do 40 i = 1, nrfx
            phifld(i,j) = gausv(one)
  40     continue
      else					! generate standard normal RF
         if( varfnc .eq. 'dlavx2' ) then
            call las2g(phifld,nrfx,nrfy,XL,YL,dlavx2,kseed,0,istat)
         elseif( varfnc .eq. 'dlsep2' ) then
            call las2g(phifld,nrfx,nrfy,XL,YL,dlsep2,kseed,0,istat)
         elseif( varfnc .eq. 'dlspx2' ) then
            call las2g(phifld,nrfx,nrfy,XL,YL,dlspx2,kseed,0,istat)
         elseif( varfnc .eq. 'dlafr2' ) then
            call las2g(phifld,nrfx,nrfy,XL,YL,dlafr2,kseed,0,istat)
         elseif( varfnc .eq. 'dlsfr2' ) then
            call las2g(phifld,nrfx,nrfy,XL,YL,dlsfr2,kseed,0,istat)
         endif
      endif


!					combine for correlated fields...
      if( lxfld ) then
         do 110 j = 1, nrfy
         do 110 i = 1, nrfx
            phifld(i,j) = R(1,2)*cfld(i,j)   + R(2,2)*phifld(i,j)
 110     continue
      endif
!					plot the random field?
      if( shofld ) then
         if( jfld .eq. 1 ) then
            call pltfld( job, sub1, sub2, cfld, nrfx,nrfx,nrfy,XL,YL,     &
                     '(Underlying) Cohesion Field', ifld )
         elseif( jfld .eq. 2 ) then
            if( ltanfi ) then
               call pltfld( job,sub1,sub2,phifld,nrfx,nrfx,nrfy,XL,YL,    &
                        '(Underlying) Friction Angle Field', ifld )
            else
               call pltfld( job,sub1,sub2,phifld,nrfx,nrfx,nrfy,XL,YL,    &
                        '(Underlying) tan(Friction Angle) Field',ifld)
            endif
         endif
      endif
!					convert to final fields
!					for cohesion...
      if( c(3) .lt. half ) then				! deterministic
         do 120 j = 1, nrfy
         do 120 i = 1, nrfx
            cfld(i,j) = c(1)
 120     continue
      elseif( c(3) .lt. onept5 ) then			! normal
         do 130 j = 1, nrfy
         do 130 i = 1, nrfx
            cfld(i,j) = c(1) + c(2)*cfld(i,j)
 130     continue
      elseif( c(3) .lt. twopt5 ) then			! lognormal
         do 140 j = 1, nrfy
         do 140 i = 1, nrfx
            cfld(i,j) = exp( c(4) + c(5)*cfld(i,j) )
 140     continue
      else						! bounded
         do 150 j = 1, nrfy
         do 150 i = 1, nrfx
            cfld(i,j) = c(4) + half*(c(5)-c(4))                           &
               *( one + tanh((c(6)+c(7)*cfld(i,j))/twopi) )
 150     continue
      endif
!					for friction angle...
      if( phi(3) .lt. half ) then			! deterministic
         do 160 j = 1, nrfy
         do 160 i = 1, nrfx
            phifld(i,j) = phi(1)
 160     continue
      elseif( phi(3) .lt. onept5 ) then			! normal
         do 170 j = 1, nrfy
         do 170 i = 1, nrfx
            phifld(i,j) = phi(1) + phi(2)*phifld(i,j)
 170     continue
      elseif( phi(3) .lt. twopt5 ) then			! lognormal
         do 180 j = 1, nrfy
         do 180 i = 1, nrfx
            phifld(i,j) = exp( phi(4) + phi(5)*phifld(i,j) )
 180     continue
      else						! bounded
         do 190 j = 1, nrfy
         do 190 i = 1, nrfx
            phifld(i,j) = phi(4) + half*(phi(5)-phi(4))                   &
                     *( one + tanh((phi(6)+phi(7)*phifld(i,j))/twopi) )
 190     continue
      endif
      return
      end
