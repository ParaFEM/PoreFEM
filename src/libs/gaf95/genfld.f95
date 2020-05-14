!  *********************************************************************
!  *                                                                   *
!  *                         subroutine genfld                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.1
!  Written by Gordon A. Fenton, DalTech, Jan 31, 2001
!  Latest Update: Feb 22, 2005
!
!  PURPOSE  generates a soil property field
!
!  DESCRIPTION
!  This routine generates a soil property field. The field may be
!  deterministic, white noise, or a correlated random field.
!
!  ARGUMENTS
!
!   istat	unit number connected to a file to which the input data is to
!		echoed if echo is true. (input)
!
!    pfld	real array of size at least (nxe x nye) which, on output,
!		will contain the underlying soil property field (ie, the
!		normally distributed field, if the property is random).
!		(output)
!
!       p	real vector of length at least 7 containing the statistics
!		of the soil property field. In particular,
!		   p(1) = soil property mean
!		   p(2) = soil property standard deviation,
!		   p(3) = soil property distribution type;
!			  = 0.0 if property is deterministic (at mean value)
!			  = 1.0 if property is normally distributed
!			  = 2.0 if property is lognormally distributed (logn)
!			  = 3.0 if property is bounded
!			  = 4.0 if property is a function of phi
!		   p(4) = lower bound (bounded), or mean log-property (logn),
!			  or constant in functional relationship (f(phi))
!		   p(5) = upper bound (bounded), or sd of log-property (logn),
!			  or slope in functional relationship (f(phi))
!		   p(6) = m parameter (if bounded), or function type (f(phi))
!			  where in the last case,
!				p(6) = 0.0 for 1*phi
!				p(6) = 1.0 for sin(phi)
!				p(6) = 2.0 for tan(phi)
!		   p(7) = s parameter (if bounded)
!		If the soil property is bounded or functional, then p(1) and
!		p(2) are ignored and the parameters c(4) through c(7)
!		completely describe the distribution (or function). (input)
!
!     nxe	integer giving the number of elements describing the soil
!		mass in the x-direction (horizontally). (input)
!
!     nye	integer giving the number of elements describing the soil
!		mass in the y-direction (vertically). (input)
!
!      XL	real value containing the physical size of the soil mass in
!		the x-direction (horizontally). (input)
!
!      YL	real value containing the physical size of the soil mass in
!		the y-direction (vertically). (input)
!
!  varfnc	character string containing the name of the variance
!		function controlling the random fields. Possible variance
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
!
!  REVISION HISTORY:
!  1.1	reordered numbering of f(phi) function types (0=1,1=sin,2=tan)
!	(this only changed the documentation above) (Feb 22/05)
!------------------------------------------------------------------------
      subroutine genfld(istat,pfld,p,nxe,nye,XL,YL,varfnc,kseed,liid)
      real pfld(nxe,*), p(*)
      real XL, YL
      integer istat, nxe, nye, kseed
      character*(*) varfnc
      logical liid
      real*8   dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      data zero/0.0/, one/1.0/
      data half/0.5/, thrpt5/3.5/

      if( p(3) .lt. half ) then			! property is deterministic
         do 10 j = 1, nye
         do 10 i = 1, nxe
            pfld(i,j) = zero			! we'll add mean back on later
  10     continue
      elseif( liid ) then			! white noise field
         do 20 j = 1, nye
         do 20 i = 1, nxe
            pfld(i,j) = gausv(one)
  20     continue
      elseif( p(3) .lt. thrpt5 ) then		! generate standard normal
         if( varfnc .eq. 'dlavx2' ) then
            call las2g(pfld,nxe,nye,XL,YL,dlavx2,kseed,0,istat)
         elseif( varfnc .eq. 'dlsep2' ) then
            call las2g(pfld,nxe,nye,XL,YL,dlsep2,kseed,0,istat)
         elseif( varfnc .eq. 'dlspx2' ) then
            call las2g(pfld,nxe,nye,XL,YL,dlspx2,kseed,0,istat)
         elseif( varfnc .eq. 'dlafr2' ) then
            call las2g(pfld,nxe,nye,XL,YL,dlafr2,kseed,0,istat)
         elseif( varfnc .eq. 'dlsfr2' ) then
            call las2g(pfld,nxe,nye,XL,YL,dlsfr2,kseed,0,istat)
         endif
      endif
!						all done
      return
      end
