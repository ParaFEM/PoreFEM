! If var < 0, then this function returns the variance of a local average
!  of the process, |var|*V(X,Y), averaged over the domain X x Y. This variance
!  is obtained by 16-pt Gauss quadrature integration of the covariance function.
!  The 4-fold integral is reduced to a 2-fold integration by taking advantage
!  of the quadrant symmetry of the covariance function.
!
! return variance reduction factor vf for normal distribution  JS May 15 2007
!
      subroutine exca(dx,dy,thx,thy,varfnc,vf)
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy, oned
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2, uservf
      real*8 dmin1,dble
      character*(*) varfnc
      external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2, uservf
      common/dparam/ dvar, dpb, dthx, dthy, dthz
      data zero/0.0/, half/0.5/, one/1.0/, two/2.0/, oned/1.d0/
            dvar = oned
            dthx = thx
            dthy = thy
            ddx  = dble(dx)
            ddy  = dble(dy)
            dvar = -dvar		! tell cov functions to return varfunc
!            if( varfnc .eq. 'dlavx2' ) then
               vf  = dlavx2(ddx,ddy)
!            elseif( varfnc .eq. 'dlsep2' ) then
!               vf = dlsep2(ddx,ddy)
!            elseif( varfnc .eq. 'dlspx2' ) then
!               vf = dlspx2(ddx,ddy)
!            elseif( varfnc .eq. 'dlafr2' ) then
!               dpb = dmin1(ddx,ddy)
!               vf  = dlafr2(ddx,ddy)
!            elseif( varfnc .eq. 'dlsfr2' ) then
!               dpb = dmin1(ddx,ddy)
!               vf  = dlsfr2(ddx,ddy)
!            elseif( varfnc .eq. 'uservf' ) then
!               vf = uservf(ddx,ddy)
!            endif
            dvar = -dvar		! tell cov functions to return cov
      return
      end
