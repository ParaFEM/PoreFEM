!  *********************************************************************
!  *                                                                   *
!  *                         subroutine simpl3                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.3
!  Written by Gordon A. Fenton, TUNS, Jul 15, 2000
!  Latest Update: Feb 28, 2005
!
!  PURPOSE  simulates soil/rock property fields for MRPILL3D
!
!  DESCRIPTION
!  This routine simulates up to 5 soil property fields (cohesion, friction
!  angle, dilation angle, elastic modulus, and Poisson's ratio), which are
!  possibly intercorrelated. Individual fields are generated as standard
!  Gaussian random fields using the 2-D Local Average Subdivision (LAS)
!  algorithm. These fields are then combined to produce correlated standard
!  Gaussian fields using the Covariance Matrix Decomposition approach.
!  Finally the individual fields are transformed so that they have the
!  desired marginal distributions. These transformations are as follows;
!
!	P(x,y) = mean + sd*G(x,y)		if normally distributed
!
!	P(x,y) = exp{ log-mean + log-sd*G(x,y) }if lognormally distributed
!
!	P(x,y) = a + 0.5*(b-a)*[ 1 + tanh((m + s*G(x,y))/2*pi) ]
!						if bounded
!
!  where P(x,y) is the desired random field (one of the soil properties),
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
!    cfld	real array of size at least (nrfx x nrfy x nrfz) which on
!		output will contain the (optionally random) cohesion field.
!		(output)
!
! phifld	real array of size at least (nrfx x nrfy x nrfz) which on
!		output will contain the (optionally random) friction angle or
!		tan(friction angle) field. (output)
!
! psifld	real array of size at least (nrfx x nrfy x nrfz) which on
!		output will contain the (optionally random) dilation angle
!		field. (output)
!
!   efld	real array of size at least (nrfx x nrfy x nrfz) which on
!		output will contain the (optionally random) elastic modulus
!		field. (output)
!
!   vfld	real array of size at least (nrfx x nrfy x nrfz) which on
!		output will contain the (optionally random) Poisson's ratio
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
!    phi	real vector of length at least 7 containing the parameters of
!		soil friction angle. See `c' for what the various elements of
!		phi contain. (input)
!
!    psi	real vector of length at least 7 containing the parameters of
!		soil dilation angle. See `c' for what the various elements of
!		psi contain. (input)
!
!	e	real vector of length at least 7 containing the parameters of
!		soil elastic modulus. See `c' for what the various elements of
!		e contain. (input)
!
!	v	real vector of length at least 7 containing the parameters of
!		soil Poisson ratio. See `c' for what the various elements of
!		v contain. (input)
!
!	R	real array of size at least 5 x 5 which, on output, will
!		contain the correlation matrix between the 5 (possibly)
!		random soil properties.  Indexing into R is as follows;
!		  1 = cohesion
!		  2 = friction angle
!		  3 = dilation angle
!		  4 = elastic modulus
!		  5 = Poisson's ratio
!		(input)
!
!   lxfld	logical flag which is true if more than one soil property
!		are cross-correlated. That is, if all soil properties are
!		independent, then lxfld is false. (input)
!
!     thx	real value giving the x-direction scale of fluctuation
!		(or, at least, this is the 1st parameter of the variance
!		function). (input)
!
!     thy	real value giving the y-direction scale of fluctuation
!		(or, at least, this is the 2nd parameter of the variance
!		function). (input)
!
!     thz	real value giving the z-direction scale of fluctuation
!		(or, at least, this is the 3rd parameter of the variance
!		function). (input)
!
!    nrfx	integer giving the number of elements describing the random
!		field in the x-direction (width). (input)
!
!    nrfy	integer giving the number of elements describing the random
!		field in the y-direction (depth). (input)
!
!    nrfz	integer giving the number of elements describing the random
!		field in the z-direction (height). (input)
!
!      dx	real value giving the physical size of an element in the
!		x-direction. (input)
!
!      dy	real value giving the physical size of an element in the
!		y-direction. (input)
!
!      dz	real value giving the physical size of an element in the
!		z-direction. (input)
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
!		 = 4 for the elastic modulus (psi) field,
!		 = 5 for the poisson ratio (psi) field,
!		(input)
!
!   ieplt	integer vector of length at least 3 containing the element
!		indices corresponding to the plane over which the random field
!		plot is to be drawn. The elements of ieplt correspond
!		to the (x,y,z) element indices and zero values in this vector
!		correspond to all elements in that direction. Only one of the
!		components is expected to be non-zero. So, for example, the
!		vector (0,26,0) means that the x-z plane at y-element = 26 is
!		to be drawn. If a 50 x 50 x 30 element domain has two
!		footings, one centered at x-node = 26 and y-node = 15 and
!		the other centered at x-node = 26 and y-node = 35, then
!		normally the y-z plane through x-node = 26 would be drawn
!		to see the random field under the two footings. Thus, in this
!		case, idplt might be (25,0,0) or (26,0,0) (depending on which
!		side of the nodal plane you want to be on). The z-indices are
!		measured from the top surface of the soil mass, thus z = 1 is
!		the top layer of the soil mass and z = nze is the bottom
!		layer, etc. (input)
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
!  varfnc	character string containing the name of the variance
!		function controlling the random fields. Possible variance
!		functions are (input)
!		`dlavx3' - 3-D exponentially decaying (Markov) model
!			   requires X-, Y-, and Z-direction scales of
!			   fluctuation
!		`dlafs3' - 3-D partially separable fractional Gaussian noise
!			   model. Requires (H,G,delta) as parameters. In this
!			   case, thx is H, thy is G, and delta is the minimum
!			   element dimension. G is assumed to apply in the
!			   Y-Z plane.
!		`dlsep3' - 3-D separable (1D x 1D x 1D) Markov model
!			   requires X-, Y-, and Z-direction scales of
!			   fluctuation
!		`dlsfr3' - 3-D separable fractional Gaussian noise model
!			   requires (H_x,H_y,H_z,delta) as parameters. In this
!			   case, thx is H_x, thy is H_y, thz is H_z, and
!			   delta is the minimum element dimension.
!		`dlspx3' - 3-D separable Gaussian decaying model
!			   requires X-, Y-, and Z-direction scales of
!			   fluctuation.
!
!   kseed	integer giving the seed to be used to initialize the
!		pseudo-random number generator. Subsequent runs using
!		the same seed will result in the same sequence of random
!		numbers. (input)
!
!   debug	logical flag which is true if debug information is to be
!		sent to the *.stt file. (input)
!
!  REVISION HISTORY:
!  1.1	eliminate lvarfn, corrected computation of pb, ZL corrected (Feb 13/02)
!  1.2	added ieplt(*) to specify RF output plot plane (Feb 28/05)
!-------------------------------------------------------------------------
      subroutine simpl3(istat,iterm,verbos,cfld,phifld,psifld,gamfld,efld,&
	                   vfld,c,phi,psi,gam,e,v,R,lxfld,thx,thy,thz,nrfx,   &
                       nrfy,nrfz,dx,dy,dz,dmpfld,nfld,jfld,ieplt,         &
                       ifld,job,sub1,sub2,varfnc,kseed,debug)

      real cfld(nrfx,nrfy,*), phifld(nrfx,nrfy,*), psifld(nrfx,nrfy,*),   &
	       gamfld(nrfx,nrfy,*)
      real efld(nrfx,nrfy,*), vfld(nrfx,nrfy,*)
      real c(*), phi(*), psi(*), gam(*),e(*), v(*), R(6,*), thx, thy, thz
      real dx, dy, dz
      integer nrfx, nrfy, nrfz, nfld, ifld, ieplt(*)
      character*(*) job, sub1, sub2, varfnc
      logical verbos, dmpfld, debug, liid, shofld, lxfld
      real*8 dvar, dpb, dthx, dthy, dthz
      real*8 dlavx3, dlsep3, dlspx3, dlafs3, dlsfr3
      save XL, YL, ZL, ienter, liid
      external dlavx3, dlsep3, dlspx3, dlafs3, dlsfr3
      common/dparam/ dvar, dpb, dthx, dthy, dthz
      data zero/0.0/, half/0.5/, one/1.0/, onept5/1.5/, twopt5/2.5/
      data twopi/6.2831853071795864769d0/
      data ienter/0/

   1  format(a,a)

!-------------------------------------- initialize -------------------------
!					compute required field size (once)
      ienter = ienter + 1
      if( ienter .eq. 1 ) then
         liid   = (thx .eq. zero) .and. (thy .eq. zero) .and.             &
                 (thz .eq. zero)
         XL     = float(nrfx)*dx
         YL     = float(nrfy)*dy
         ZL     = float(nrfz)*dz

         if( debug ) then
            write(istat,1)'simpl3: initializing field parameters...'
            call flush(istat)
         endif

         if( .not. liid ) then
!						set covariance fnc parameters
            dvar   = 1.d0
            dthx   = thx
            dthy   = thy
            dthz   = thz
            if( varfnc .eq. 'dlafs3' .or. varfnc .eq. 'dlsfr3' ) then
               pb = dx
               if( dy .lt. pb ) pb = dy
               if( dz .lt. pb ) pb = dz
               dpb = pb
            endif
!						initialize LAS2G
            if( varfnc .eq. 'dlavx3' ) then
               call las3g(cfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlavx3,            &
                         kseed,-1,istat)
            elseif( varfnc .eq. 'dlsep3' ) then
               call las3g(cfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsep3,            &
                         kseed,-1,istat)
            elseif( varfnc .eq. 'dlspx3' ) then
               call las3g(cfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlspx3,            &
                         kseed,-1,istat)
            elseif( varfnc .eq. 'dlafs3' ) then
               call las3g(cfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlafs3,            &
                         kseed,-1,istat)
            elseif( varfnc .eq. 'dlsfr3' ) then
               call las3g(cfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsfr3,            &
                         kseed,-1,istat)
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
         if( debug ) then
            write(istat,1)'simpl3: finished initialization...'
            call flush(istat)
         endif
      endif
!					are we going to plot anything?
      shofld = dmpfld .and. (nfld .eq. ienter)

!					now produce 6 standard normal fields
!					for cohesion...
      if( c(3) .lt. half ) then			! cohesion is deterministic
         do 10 k = 1, nrfz
         do 10 j = 1, nrfy
         do 10 i = 1, nrfx
            cfld(i,j,k) = zero			! we'll add mean back on later
  10     continue
      elseif( liid ) then			! white noise field
         do 20 k = 1, nrfz
         do 20 j = 1, nrfy
         do 20 i = 1, nrfx
            cfld(i,j,k) = gausv(one)
  20     continue
      else					! generate standard normal
         if( varfnc .eq. 'dlavx3' ) then
            call las3g(cfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlavx3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsep3' ) then
            call las3g(cfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsep3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlspx3' ) then
            call las3g(cfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlspx3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlafs3' ) then
            call las3g(cfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlafs3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsfr3' ) then                               
            call las3g(cfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsfr3,               &
                      kseed,0,istat)
         endif
      endif
!					for friction angle...
      if( phi(3) .lt. half ) then		! friction is deterministic
         do 30 k = 1, nrfz
         do 30 j = 1, nrfy
         do 30 i = 1, nrfx
            phifld(i,j,k) = zero		! we'll add mean back on later
  30     continue
      elseif( liid ) then			! white noise field
         do 40 k = 1, nrfz
         do 40 j = 1, nrfy
         do 40 i = 1, nrfx
            phifld(i,j,k) = gausv(one)
  40     continue
      else					! generate standard normal
         if( varfnc .eq. 'dlavx3' ) then
            call las3g(phifld,nrfx,nrfy,nrfz,XL,YL,ZL,dlavx3,             &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsep3' ) then
            call las3g(phifld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsep3,             &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlspx3' ) then
            call las3g(phifld,nrfx,nrfy,nrfz,XL,YL,ZL,dlspx3,             &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlafs3' ) then
            call las3g(phifld,nrfx,nrfy,nrfz,XL,YL,ZL,dlafs3,             &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsfr3' ) then
            call las3g(phifld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsfr3,             &
                      kseed,0,istat)
         endif
      endif
!					for dilation angle...
      if( psi(3) .lt. half ) then		! dilation is deterministic
         do 50 k = 1, nrfz
         do 50 j = 1, nrfy
         do 50 i = 1, nrfx
            psifld(i,j,k) = zero		! we'll add mean back on later
  50     continue
      elseif( liid ) then			! white noise field
         do 60 k = 1, nrfz
         do 60 j = 1, nrfy
         do 60 i = 1, nrfx
            psifld(i,j,k) = gausv(one)
  60     continue
      else					! generate standard normal
         if( varfnc .eq. 'dlavx3' ) then
            call las3g(psifld,nrfx,nrfy,nrfz,XL,YL,ZL,dlavx3,             &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsep3' ) then
            call las3g(psifld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsep3,             &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlspx3' ) then
            call las3g(psifld,nrfx,nrfy,nrfz,XL,YL,ZL,dlspx3,             &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlafs3' ) then
            call las3g(psifld,nrfx,nrfy,nrfz,XL,YL,ZL,dlafs3,             &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsfr3' ) then
            call las3g(psifld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsfr3,             &
                      kseed,0,istat)
         endif
      endif
!					for elastic modulus...
      if( e(3) .lt. half ) then		! modulus is deterministic
         do 70 k = 1, nrfz
         do 70 j = 1, nrfy
         do 70 i = 1, nrfx
            efld(i,j,k) = zero			! we'll add mean back on later
  70     continue
      elseif( liid ) then			! white noise field
         do 80 k = 1, nrfz
         do 80 j = 1, nrfy
         do 80 i = 1, nrfx
            efld(i,j,k) = gausv(one)
  80     continue
      else					! generate standard normal
         if( varfnc .eq. 'dlavx3' ) then
            call las3g(efld,nrfx,nrfy,nrfz,XL,YL,ZL,dlavx3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsep3' ) then
            call las3g(efld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsep3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlspx3' ) then
            call las3g(efld,nrfx,nrfy,nrfz,XL,YL,ZL,dlspx3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlafs3' ) then
            call las3g(efld,nrfx,nrfy,nrfz,XL,YL,ZL,dlafs3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsfr3' ) then
            call las3g(efld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsfr3,               &
                      kseed,0,istat)
         endif
      endif
!					for Poisson's ratio...
      if( v(3) .lt. half ) then		! ratio is deterministic
         do 90 k = 1, nrfz
         do 90 j = 1, nrfy
         do 90 i = 1, nrfx
            vfld(i,j,k) = zero			! we'll add mean back on later
  90     continue
      elseif( liid ) then			! white noise field
         do 100 k = 1, nrfz
         do 100 j = 1, nrfy
         do 100 i = 1, nrfx
            vfld(i,j,k) = gausv(one)
 100     continue
      else					! generate standard normal
         if( varfnc .eq. 'dlavx3' ) then
            call las3g(vfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlavx3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsep3' ) then
            call las3g(vfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsep3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlspx3' ) then
            call las3g(vfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlspx3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlafs3' ) then
            call las3g(vfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlafs3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsfr3' ) then
            call las3g(vfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsfr3,               &
                      kseed,0,istat)
         endif
      endif
!					for Unit weight...
      if( gam(3) .lt. half ) then		! ratio is deterministic
         do 901 k = 1, nrfz
         do 901 j = 1, nrfy
         do 901 i = 1, nrfx
            gamfld(i,j,k) = zero			! we'll add mean back on later
  901     continue
      elseif( liid ) then			! white noise field
         do 1001 k = 1, nrfz
         do 1001 j = 1, nrfy
         do 1001 i = 1, nrfx
            gamfld(i,j,k) = gausv(one)
 1001     continue
      else					! generate standard normal
         if( varfnc .eq. 'dlavx3' ) then
            call las3g(gamfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlavx3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsep3' ) then
            call las3g(gamfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsep3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlspx3' ) then
            call las3g(gamfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlspx3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlafs3' ) then
            call las3g(gamfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlafs3,               &
                      kseed,0,istat)
         elseif( varfnc .eq. 'dlsfr3' ) then
            call las3g(gamfld,nrfx,nrfy,nrfz,XL,YL,ZL,dlsfr3,               &
                      kseed,0,istat)
         endif
      endif
!					combine for correlated fields...
      if( lxfld ) then
         do 110 k = 1, nrfz
         do 110 j = 1, nrfy
         do 110 i = 1, nrfx
            vfld(i,j,k)   = R(1,6)*cfld(i,j,k)   + R(2,6)*phifld(i,j,k)   &
                         + R(3,6)*psifld(i,j,k) + R(4,6)*gamfld(i,j,k)    &
						 + R(5,6)*efld(i,j,k)   + R(6,6)*vfld(i,j,k)
            efld(i,j,k)   = R(1,5)*cfld(i,j,k)   + R(2,5)*phifld(i,j,k)   &
                         + R(3,5)*psifld(i,j,k) + R(4,5)*gamfld(i,j,k)    &
						 + R(5,5)*efld(i,j,k)
            gamfld(i,j,k)  = R(1,4)*cfld(i,j,k)   + R(2,4)*phifld(i,j,k)   &
                         + R(3,4)*psifld(i,j,k) + R(4,4)*gamfld(i,j,k)
            psifld(i,j,k) = R(1,3)*cfld(i,j,k)   + R(2,3)*phifld(i,j,k)   &
                         + R(3,3)*psifld(i,j,k)
            phifld(i,j,k) = R(1,2)*cfld(i,j,k)   + R(2,2)*phifld(i,j,k)
 110     continue
      endif
!					plot the random field?
      if( shofld ) then
         if( jfld .eq. 1 ) then
	    call pltfld( job, sub1, sub2, cfld, nrfx, nrfy,                   &
     		         nrfx, nrfy, nrfz, XL, YL, ZL, ieplt,                 &
                     '(Underlying) Cohesion Field', ifld )
         elseif( jfld .eq. 2 ) then
	    call pltfld( job, sub1, sub2, phifld, nrfx, nrfy,                 &
     		         nrfx, nrfy, nrfz, XL, YL, ZL, ieplt,                 &
     		     '(Underlying) Friction Angle Field', ifld )
         elseif( jfld .eq. 3 ) then
	    call pltfld( job, sub1, sub2, psifld, nrfx, nrfy,                 &
     		         nrfx, nrfy, nrfz, XL, YL, ZL, ieplt,                 &
     		     '(Underlying) Dilation Angle Field', ifld )
           elseif( jfld .eq. 4 ) then
	    call pltfld( job, sub1, sub2, gamfld, nrfx, nrfy,                   &
     		         nrfx, nrfy, nrfz, XL, YL, ZL, ieplt,                 &
     		     '(Underlying) Elastic Modulus Field', ifld )
         elseif( jfld .eq. 5 ) then
	    call pltfld( job, sub1, sub2, efld, nrfx, nrfy,                   &
     		         nrfx, nrfy, nrfz, XL, YL, ZL, ieplt,                 &
     		     '(Underlying) Elastic Modulus Field', ifld )
         elseif( jfld .eq. 6 ) then
	    call pltfld( job, sub1, sub2, vfld, nrfx, nrfy,                   &
     		         nrfx, nrfy, nrfz, XL, YL, ZL, ieplt,                 &
     		     '(Underlying) Poisson''s Ratio Field', ifld )
         endif
      endif
!					convert to final fields
!					for cohesion...
      if( c(3) .lt. half ) then				! deterministic
         do 120 k = 1, nrfz
         do 120 j = 1, nrfy
         do 120 i = 1, nrfx
            cfld(i,j,k) = c(1)
 120     continue
      elseif( c(3) .lt. onept5 ) then			! normal
         do 130 k = 1, nrfz
         do 130 j = 1, nrfy
         do 130 i = 1, nrfx
            cfld(i,j,k) = c(1) + c(2)*cfld(i,j,k)
 130     continue
      elseif( c(3) .lt. twopt5 ) then			! lognormal
         do 140 k = 1, nrfz
         do 140 j = 1, nrfy
         do 140 i = 1, nrfx
            cfld(i,j,k) = exp( c(4) + c(5)*cfld(i,j,k) )
 140     continue
      else						! bounded
         do 150 k = 1, nrfz
         do 150 j = 1, nrfy
         do 150 i = 1, nrfx
            cfld(i,j,k) = c(4) + half*(c(5)-c(4))                         &
               *( one + tanh((c(6)+c(7)*cfld(i,j,k))/twopi) )
 150     continue
      endif
!					for friction angle...
      if( phi(3) .lt. half ) then			! deterministic
         do 160 k = 1, nrfz
         do 160 j = 1, nrfy
         do 160 i = 1, nrfx
            phifld(i,j,k) = phi(1)
 160     continue
      elseif( phi(3) .lt. onept5 ) then			! normal
         do 170 k = 1, nrfz
         do 170 j = 1, nrfy
         do 170 i = 1, nrfx
            phifld(i,j,k) = phi(1) + phi(2)*phifld(i,j,k)
 170     continue
      elseif( phi(3) .lt. twopt5 ) then			! lognormal
         do 180 k = 1, nrfz
         do 180 j = 1, nrfy
         do 180 i = 1, nrfx
            phifld(i,j,k) = exp( phi(4) + phi(5)*phifld(i,j,k) )
 180     continue
      else						! bounded
         do 190 k = 1, nrfz
         do 190 j = 1, nrfy
         do 190 i = 1, nrfx
            phifld(i,j,k) = phi(4) + half*(phi(5)-phi(4))                 &
                     *(one + tanh((phi(6)+phi(7)*phifld(i,j,k))/twopi))
 190     continue
      endif
!					for dilation angle...
      if( psi(3) .lt. half ) then			! deterministic
         do 200 k = 1, nrfz
         do 200 j = 1, nrfy
         do 200 i = 1, nrfx
            psifld(i,j,k) = psi(1)
 200     continue
      elseif( psi(3) .lt. onept5 ) then			! normal
         do 210 k = 1, nrfz
         do 210 j = 1, nrfy
         do 210 i = 1, nrfx
            psifld(i,j,k) = psi(1) + psi(2)*psifld(i,j,k)
 210     continue
      elseif( psi(3) .lt. twopt5 ) then			! lognormal
         do 220 k = 1, nrfz
         do 220 j = 1, nrfy
         do 220 i = 1, nrfx
            psifld(i,j,k) = exp( psi(4) + psi(5)*psifld(i,j,k) )
 220     continue
      else						! bounded
         do 230 k = 1, nrfz
         do 230 j = 1, nrfy
         do 230 i = 1, nrfx
            psifld(i,j,k) = psi(4) + half*(psi(5)-psi(4))                 &
                     *(one + tanh((psi(6)+psi(7)*psifld(i,j,k))/twopi))
 230     continue
      endif
!					for elastic modulus...
      if( e(3) .lt. half ) then				! deterministic
         do 240 k = 1, nrfz
         do 240 j = 1, nrfy
         do 240 i = 1, nrfx
            efld(i,j,k) = e(1)
 240     continue
      elseif( e(3) .lt. onept5 ) then			! normal
         do 250 k = 1, nrfz
         do 250 j = 1, nrfy
         do 250 i = 1, nrfx
            efld(i,j,k) = e(1) + e(2)*efld(i,j,k)
 250     continue
      elseif( e(3) .lt. twopt5 ) then			! lognormal
         do 260 k = 1, nrfz
         do 260 j = 1, nrfy
         do 260 i = 1, nrfx
            efld(i,j,k) = exp( e(4) + e(5)*efld(i,j,k) )
 260     continue
      else						! bounded
         do 270 k = 1, nrfz
         do 270 j = 1, nrfy
         do 270 i = 1, nrfx
            efld(i,j,k) = e(4) + half*(e(5)-e(4))                         &
                     *( one + tanh((e(6)+e(7)*efld(i,j,k))/twopi) )
 270     continue
      endif
!					for Poisson's ratio...
      if( v(3) .lt. half ) then				! deterministic
         do 280 k = 1, nrfz
         do 280 j = 1, nrfy
         do 280 i = 1, nrfx
            vfld(i,j,k) = v(1)
 280     continue
      elseif( v(3) .lt. onept5 ) then			! normal
         do 290 k = 1, nrfz
         do 290 j = 1, nrfy
         do 290 i = 1, nrfx
            vfld(i,j,k) = v(1) + v(2)*vfld(i,j,k)
 290     continue
      elseif( v(3) .lt. twopt5 ) then			! lognormal
         do 300 k = 1, nrfz
         do 300 j = 1, nrfy
         do 300 i = 1, nrfx
            vfld(i,j,k) = exp( v(4) + v(5)*vfld(i,j,k) )
 300     continue
      else						! bounded
         do 310 k = 1, nrfz
         do 310 j = 1, nrfy
         do 310 i = 1, nrfx
            vfld(i,j,k) = v(4) + half*(v(5)-v(4))                         &
                     *( one + tanh((v(6)+v(7)*vfld(i,j,k))/twopi) )
 310     continue
      endif
!					for Unit weight...
      if( gam(3) .lt. half ) then				! deterministic
         do 281 k = 1, nrfz
         do 281 j = 1, nrfy
         do 281 i = 1, nrfx
            gamfld(i,j,k) = gam(1)
 281     continue
      elseif( gam(3) .lt. onept5 ) then			! normal
         do 291 k = 1, nrfz
         do 291 j = 1, nrfy
         do 291 i = 1, nrfx
            gamfld(i,j,k) = gam(1) + gam(2)*gamfld(i,j,k)
 291     continue
      elseif( gam(3) .lt. twopt5 ) then			! lognormal
         do 301 k = 1, nrfz
         do 301 j = 1, nrfy
         do 301 i = 1, nrfx
            gamfld(i,j,k) = exp( gam(4) + gam(5)*gamfld(i,j,k) )
 301     continue
      else						! bounded
         do 311 k = 1, nrfz
         do 311 j = 1, nrfy
         do 311 i = 1, nrfx
            gamfld(i,j,k) = gam(4) + half*(gam(5)-gam(4))                         &
                     *( one + tanh((gam(6)+gam(7)*gamfld(i,j,k))/twopi) )
 311     continue
      endif

      if( debug ) then
         write(istat,1)'simpl3: finished everything...'
         call flush(istat)
      endif
      return
      end
