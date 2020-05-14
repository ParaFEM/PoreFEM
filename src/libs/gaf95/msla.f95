! Comopute locally averaged statistics of lognormal variable
!                                          JS May 15 2007
      subroutine msla(dx,dy,thx,thy,kmn,ksd,varfnc,istat,kmne,ksde)
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy, oned
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2, uservf
      real*8 dmin1,dble
      real kmn, kmne, ksd, ksde
      character*(*) varfnc
      external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2, uservf
      common/dparam/ dvar, dpb, dthx, dthy, dthz
      data zero/0.0/, half/0.5/, one/1.0/, two/2.0/, oned/1.d0/
   1  format(a,a)
   3  format(a,e12.5)
         pvr  = alog(one + ksd*ksd/(kmn*kmn))
         pmn  = alog(kmn) - half*pvr
         psd  = sqrt(pvr)
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
            write(istat,3)'Variance reduction   = ',vf
            write(istat,3)'Point Log-Cond. mean = ',pmn
            write(istat,3)'Point Log-Cond. S.D. = ',psd
            write(istat,3)'Cell Cond. mean      = ',kmne
            write(istat,3)'Cell Cond. S.D.      = ',ksde
            write(istat,'()')
      return
      end
