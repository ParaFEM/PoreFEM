!  ********************************************************************
!  *                                                                  *
!  *                        subroutine sim2dk                         *
!  *                                                                  *
!  ********************************************************************
!  Single Precision Version 3.0
!  Written by Gordon A. Fenton, TUNS, 1991
!  Latest Update: May 10, 2001
!
!  PURPOSE  to generate a realization of a soil conductivity field
!
!  This routine creates a realization of a soil conductivity field in 2-D.
!  The generation proceeds as follows;
!
!   1) generate a zero-mean, unit-variance, Gaussian random field in 2-D
!      with prescribed covariance structure.
!
!   2) transform the field into a lognormally distributed conductivity field
!      on a point-wise basis by setting Y = exp(ksdl + kmnl*G) at each field
!      coordinate (kmnl and ksdl being the mean and variance of the
!      log-conductivity field).
!
!   3) map discrete conductivity values into the finite element mesh, so that
!      each element takes on a particular conductivity value
!
!  Arguments are described as follows;
!
!    perm    real array of size at least NXE x NYE which on output will
!            contain the desired realization of conductivity.
!            The leading dimension of perm is assumed to be NXE.
!c            (output)
!c
!c     nxe    number of elements in the X direction. (input)
!c
!     nye    number of elements in the Y direction. (input)
!
!      dx    physical size of an element in the X direction. (input)
!
!      dy    physical size of an element in the Y direction. (input)
!
!     thx    scale of fluctuation of log-conductivity in the X direction.
!            Note that if either thx or thy is zero, the entire field is
!            assumed to be made up of independent values (ie a white
!            noise process). In this case, the mean and variance of each
!            cell is determined directly from kmn? and ksd?, that is the
!            effects of local averaging is not considered. (input)
!
!     thy    scale of fluctuation of log-conductivity in the Y direction.
!            See note on thx. (input)
!
!     kmn    (real) target mean of conductivity at a point. (input)
!
!     ksd    (real) target standard deviation of conductivity at a point.
!            (input)
!
!   lunif    logical flag which is true if the conductivity field is
!            uniform, corresponding to an infinite scale of fluctuation
!            in both directions. In this case, the conductivity field
!            is set equal to a single random value having mean kmn and
!            standard deviation ksd. (input)
!
!   kseed    integer seed to be used for the pseudo-random number generator.
!            If KSEED = 0, then a random seed will be used (based on the
!            clock time when LAS2D is called for the first time).
!            On output, KSEED is set to the value of the actual seed used.
!            (input/output)
!
!  dmpfld    logical flag which is true if the first generated random field
!            is to be dumped to a DISPLAY format file. Only one field is
!            dumped, so dmpfld is set to false after dumping. (input/output)
!
!   iperm    unit number to which the first field of log-conductivity is
!            dumped. (input)
!
!     job    character string containing the title of the run. (input)
!
!    sub1    character string containing the subtitle of the run. (input)
!
!    sub2    character string containing the sub-subtitle of the run. (input)
!
!  varfnc    character*6 string containing the name of the var/covar function
!            to be used by las2g. See las2g(3f) for details. (input)
!
!   debug    logical flag which is true if debugging information is to be
!            output to unit istat. (input)
!
!   istat    unit number to which error and warning messages are issued.
!            (input)
!
!  dmpagh    logical flag which is true if arithmetic, geometric and harmonic
!            means are to be calculated for the conductivity field. (input)
!
!  kstats    real array of size at least 3 x 4 containing conductivity stats;
!            (input/output)
!	       (1,1) = current block conductivity value
!	       (2,1) = sum of block conductivities (over the realizations)
!	       (3,1) = sum of (block conductivities)**2
!	       (1,2) = current arithmetic mean of conductivity field
!	       (2,2) = sum of arithmetic means (over the realizations)
!	       (3,2) = sum of (arithmetic means)**2
!	       (1,3) = current geometric mean of conductivity field
!	       (2,3) = sum of geometric means (over the realizations)
!	       (3,3) = sum of (geometric means)**2
!	       (1,4) = current harmonic mean of conductivity field
!	       (2,4) = sum of harmonic means (over the realizations)
!	       (3,4) = sum of (harmonic means)**2
!
!    kmne    (real) target mean of conductivity over a finite element (ie.
!            after locally averaging over (dx x dy)). (output)
!
!    ksde    (real) target standard deviation of conductivity over a finite
!            element (ie. after locally averaging over (dx x dy)). (output)
!
!    akmn    real value which which contains a running sum of
!            conductivities over each field and over all realizations.
!            (input/output)
!
!    aksd    real value which which contains a running sum of squared
!            conductivities over each field and over all realizations.
!            (input/output)
!
!  dcheck    logical flag which is true if this is a data-check run only, in
!            which case only the required field size and target statistics
!            are actually computed herein. (input)
!
!  Revision History:
!
!  2.5	now using libGAFsim variance functions dlavx2, dlsep2, and dlspx2,
!	rather than vfrtxx, vfsepx, and vfspxx (these are equivalent, they
!	are just part of the library now and renamed).
!  2.6	added the 2-D fractal process variance functions dlvfr2 and dlsfr2
!	(these come from libGAFsim).
!  2.7	made conductivity isotropic, assumes that kmn and ksd refer to the
!	point statistics (now kmne, ksde return the cell statistics).
!  2.71	corrected header.
!  2.8	improved the handling of the uniform field case (lunif = true)
!  2.9	LAS2G now uses covariances to produce simulations (via lvarfn)
!	(May 23/00)
!  3.0	removed lvarfn from common/dparam/ (LAS2G now uses Gauss quadrature
!	for all covariance calculations) (May 10/01)
!---------------------------------------------------------------------------
      subroutine sim2dk(perm,nxe,nye,dx,dy,thx,thy,kmn,ksd,               &
	               lunif,kseed,dmpfld,iperm,job,sub1,                     &
                       sub2,varfnc,debug,istat,dmpagh,                    &
                       kstats,kmne,ksde,akmn,aksd,dcheck)
      real perm(nxe,*), kstats(3,*)
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy, oned
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2, uservf
      real*8 dmin1, dble
      real kmn, kmne, ksd, ksde
      logical debug, dmpfld, liid, lunif, dmpagh, dcheck
      character*(*) job, sub1, sub2, varfnc
      external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2, uservf
      save pmn, psd, liid, XL, YL, nn
!				export parameters to variance function
      common/dparam/ dvar, dpb, dthx, dthy, dthz
      data ifirst/1/
      data zero/0.0/, half/0.5/, one/1.0/, two/2.0/, oned/1.d0/

   1  format(a,a)
   2  format(a,i4,a,i4,a)
   3  format(a,e12.5)
   4  format()

!-------------------------------------- initialize -------------------------
!					compute required field size (once)
      if( ifirst .eq. 1 ) then
         ifirst = 0
         liid   = (thx .eq. zero .or. thy .eq. zero)
         nn     = nxe*nye
         XL     = float(nxe)*dx
         YL     = float(nye)*dy

         ierr = 0
         if( debug ) write(istat,1)'SIM2DK: setting field parameters...'
         pvr  = alog(one + ksd*ksd/(kmn*kmn))
         pmn  = alog(kmn) - half*pvr
         psd  = sqrt(pvr)
         if( .not. liid .and. .not. lunif ) then
!						set variance fnc parameters
            dvar = oned
            dthx = thx
            dthy = thy
!					find mean and s.d. of log-conductivity
!					averaged over the element
            ddx  = dble(dx)
            ddy  = dble(dy)
            dvar = -dvar		! tell cov functions to return varfunc
            if( varfnc .eq. 'dlavx2' ) then
               vf  = dlavx2(ddx,ddy)
            elseif( varfnc .eq. 'dlsep2' ) then
               vf = dlsep2(ddx,ddy)
            elseif( varfnc .eq. 'dlspx2' ) then
               vf = dlspx2(ddx,ddy)
            elseif( varfnc .eq. 'dlafr2' ) then
               dpb = dmin1(ddx,ddy)
               vf  = dlafr2(ddx,ddy)
            elseif( varfnc .eq. 'dlsfr2' ) then
               dpb = dmin1(ddx,ddy)
               vf  = dlsfr2(ddx,ddy)
            elseif( varfnc .eq. 'uservf' ) then
               vf = uservf(ddx,ddy)
            else
               write(istat,1)'*** Error: unknown variance function name: ', varfnc
               stop
            endif
            dvar = -dvar		! tell cov functions to return cov
            kmne = exp(pmn + half*vf*pvr)
            ksde = kmne*sqrt((exp(vf*pvr) - one))
!						print out some information

            write(istat,2)'Specified Field size = (',nxe,' x ',nye,')'
            write(istat,3)'Variance reduction   = ',vf
            write(istat,3)'Point Log-Cond. mean = ',pmn
            write(istat,3)'Point Log-Cond. S.D. = ',psd
            write(istat,3)'Cell Cond. mean      = ',kmne
            write(istat,3)'Cell Cond. S.D.      = ',ksde
            write(istat,'()')
!						initialize LAS2G
	    if( varfnc .eq. 'dlavx2' ) then
	       call las2g(perm,nxe,nye,XL,YL,dlavx2,kseed,-1,istat)
	    elseif( varfnc .eq. 'dlsep2' ) then
	       call las2g(perm,nxe,nye,XL,YL,dlsep2,kseed,-1,istat)
	    elseif( varfnc .eq. 'dlspx2' ) then
	       call las2g(perm,nxe,nye,XL,YL,dlspx2,kseed,-1,istat)
	    elseif( varfnc .eq. 'dlafr2' ) then
	       call las2g(perm,nxe,nye,XL,YL,dlafr2,kseed,-1,istat)
	    elseif( varfnc .eq. 'dlsfr2' ) then
	       call las2g(perm,nxe,nye,XL,YL,dlsfr2,kseed,-1,istat)
	    elseif( varfnc .eq. 'uservf' ) then
	       call las2g(perm,nxe,nye,XL,YL,uservf,kseed,-1,istat)
	    endif
         else
!					for a white noise or uniform field
            kmne  = kmn
            ksde  = ksd
            kseed = iseed(kseed)
         endif
      endif

!------------------------------------ generate the random field(s) -----------
      if( dcheck ) return
      aa = zero
      va = zero
!						white noise field
      if( liid ) then
         if( dmpagh ) then
            ag = zero
            ah = zero
            do 10 iq = 1, nye
            do 10 ip = 1, nxe
               gk = pmn + psd*gausv(one)
               pk = exp( gk )
               ph = one/pk
               aa = aa + pk
               va = va + pk*pk
               ag = ag + gk
               ah = ah + ph
               perm(ip,iq) = pk
  10        continue
         else
            do 20 iq = 1, nye
            do 20 ip = 1, nxe
               perm(ip,iq)   = exp( pmn + psd*gausv(one) )
               aa = aa + perm(ip,iq)
               va = va + perm(ip,iq)*perm(ip,iq)
  20        continue
         endif
         if( dmpfld ) then
            call pltfld( job, sub1, sub2, perm, nxe, nxe, nye, XL, YL,'White Noise Conductivity Field', iperm )
            dmpfld = .false.
         endif
!						uniform field
      elseif( lunif ) then
         gk = pmn + psd*gausv(one)
         pk = exp( gk )
         aa = float(nn)*pk
         va = float(nn)*pk*pk
         if( dmpagh ) then
            ag = float(nn)*gk
            ah = float(nn)/pk
         endif
         do 30 iq = 1, nye
         do 30 ip = 1, nxe
            perm(ip,iq) = pk
  30     continue
      else
!						random field
	 if( varfnc .eq. 'dlavx2' ) then
	    call las2g(perm,nxe,nye,XL,YL,dlavx2,kseed,0,istat)
	 elseif( varfnc .eq. 'dlsep2' ) then
	    call las2g(perm,nxe,nye,XL,YL,dlsep2,kseed,0,istat)
	 elseif( varfnc .eq. 'dlspx2' ) then
	    call las2g(perm,nxe,nye,XL,YL,dlspx2,kseed,0,istat)
	 elseif( varfnc .eq. 'dlafr2' ) then
	    call las2g(perm,nxe,nye,XL,YL,dlafr2,kseed,0,istat)
	 elseif( varfnc .eq. 'dlsfr2' ) then
	    call las2g(perm,nxe,nye,XL,YL,dlsfr2,kseed,0,istat)
	 elseif( varfnc .eq. 'uservf' ) then
	    call las2g(perm,nxe,nye,XL,YL,uservf,kseed,0,istat)
	 endif
!						plot the first field?
	 if( dmpfld ) then
	    call pltfld( job, sub1, sub2, perm, nxe, nxe, nye, XL, YL,'Log-Conductivity Field', iperm )
	    dmpfld = .false.
	 endif
!						log-field --> field
         if( dmpagh ) then
            ag = zero
            ah = zero
            do 40 iq = 1, nye
            do 40 ip = 1, nxe
               gk = pmn + psd*perm(ip,iq)
               ag = ag + gk
               perm(ip,iq) = exp( gk )
               aa = aa + perm(ip,iq)
               va = va + perm(ip,iq)*perm(ip,iq)
               ah = ah + one/perm(ip,iq)
  40        continue
         else
            do 50 iq = 1, nye
            do 50 ip = 1, nxe
               perm(ip,iq) = exp( pmn + psd*perm(ip,iq) )
               aa = aa + perm(ip,iq)
               va = va + perm(ip,iq)*perm(ip,iq)
  50        continue
         endif
      endif
!						conductivity statistics
      akmn = akmn + aa
      aksd = aksd + va
      kstats(1,2) = aa/float(nn)
      if( dmpagh ) then
         kstats(2,2) = kstats(2,2) + kstats(1,2)
         kstats(3,2) = kstats(3,2) + kstats(1,2)*kstats(1,2)
         kstats(1,3) = exp(ag/float(nn))
         kstats(2,3) = kstats(2,3) + kstats(1,3)
         kstats(3,3) = kstats(3,3) + kstats(1,3)*kstats(1,3)
         kstats(1,4) = float(nn)/ah
         kstats(2,4) = kstats(2,4) + kstats(1,4)
         kstats(3,4) = kstats(3,4) + kstats(1,4)*kstats(1,4)
      endif
   do i=1,16
   enddo
      return
      end
