!  *******************************************************************
!  *                                                                 *
!  *                     Function dlavx2                             *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 1.6
!  Written by Gordon A. Fenton, TUNS, 1991
!  Latest Update: Apr 8, 2003
!
!  PURPOSE  returns the covariance between two points in a 2-D Markovian
!           random field.
!
!  This function returns the covariance between two points in a 2-D
!  Markov random field, separated by lag vector {X,Y}.
!  The covariance function for this process is given by
!
!         B(x,y) = var * exp{ -sqrt[(2*x/dthx)^2 + (2*y/dthy)^2] }
!
!  where var is the point variance, and dthx and dthy are the scales of
!  fluctuation in the x- and y-directions, respectively. The parameters
!  var, dthx, and dthy, are brought in via the common block /dparam/.
!
!  If var < 0, then this function returns the variance of a local average
!  of the process, |var|*V(X,Y), averaged over the domain X x Y. This variance
!  is obtained by 16-pt Gauss quadrature integration of the covariance function.
!  The 4-fold integral is reduced to a 2-fold integration by taking advantage
!  of the quadrant symmetry of the covariance function.
!
!  See "Random Fields" by E. Vanmarcke and G.A. Fenton's thesis
!  for more details.
!
!  The arguments to this routine are just the X and Y components of the
!  lag vector between the points (or the dimensions of the averaging area,
!  if var < 0).
!
!  NOTE: when the element is non-square, extremely small scales of fluctuation
!	 tend to yield inaccurate local average variances (either directly
!	 here using var < 0, or through dcvaa2). More Gauss points are
!	 needed in the integration. See F77/sim/utils for a version of
!	 dlavx2 with 20-pt Gaussian quadrature.
!
!  REVISION HISTORY:
!  1.1	brought the pi/2 factor into the function directly, rather than
!	making it a parameter (was `dum1').
!  1.2	added lvarfn flag to end of /dparam/ common block. Set it to
!	true if the variance function is to be returned. Otherwise this
!	function returns the covariance. (Apr 29/00)
!  1.3	eliminated lvarfn -- now return covariances only if var < 0 (Mar 27/01)
!  1.4	properly handled sign on var for covariances (Apr 5/01)
!  1.5	reversed default - now return covariances if var > 0.
!	Variance function values now computed by Gauss quadrature (Apr 11/01)
!  1.51	revised above documentation to reflect revision 1.5 (May 9/01)
!  1.52	added comment on accuracy above (Jun 7/01)
!  1.6	now use 16-pt Gauss quadrature for increased accuracy (Apr 8/02)
!---------------------------------------------------------------------------
      real*8 function dlavv2( X, Y )
      parameter (NG = 16)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)
      common/dparam2/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/, one/1.d0/, two/2.d0/, three/3.d0/
      data quart/0.25d0/, half/0.5d0/
!     data pb/1.570796326794896619231322d0/, onept5/1.5d0/
!			Gauss weights and points for NG = 5
!     data w/ .236926885056189d0,.478628670499366d0,.568888888888889d0,
!    >        .478628670499366d0,.236926885056189d0/
!     data z/-.906179845938664d0,-.538469310105683d0,.000000000000000d0,
!    >        .538469310105683d0, .906179845938664d0/
!			these are for NG = 16
      data w/0.027152459411754094852d0, 0.062253523938647892863d0,&
             0.095158511682492784810d0, 0.124628971255533872052d0,&
             0.149595988816576732081d0, 0.169156519395002538189d0,&
             0.182603415044923588867d0, 0.189450610455068496285d0,&
             0.189450610455068496285d0, 0.182603415044923588867d0,&
             0.169156519395002538189d0, 0.149595988816576732081d0,&
             0.124628971255533872052d0, 0.095158511682492784810d0,&
             0.062253523938647892863d0, 0.027152459411754094852d0/
      data z/-.989400934991649932596d0, -.944575023073232576078d0,&
             -.865631202387831743880d0, -.755404408355003033895d0,&
             -.617876244402643748447d0, -.458016777657227386342d0,&
             -.281603550779258913230d0, -.095012509837637440185d0,&
             0.095012509837637440185d0, 0.281603550779258913230d0,&
             0.458016777657227386342d0, 0.617876244402643748447d0,&
             0.755404408355003033895d0, 0.865631202387831743880d0,&
             0.944575023073232576078d0, 0.989400934991649932596d0/

      exp(zz)  = dexp(zz)
      abs(zz)  = dabs(zz)
      sqrt(zz) = dsqrt(zz)

      aY = abs(Y)
      aX = abs(X)
      if( var .lt. zero ) then			! return variance function
         if( (dthx .eq. zero) .and. (dthy .eq. zero) ) then
            if( (X .eq. zero) .and. (Y .eq. zero) ) then
               dlavv2 = -var
            else
               dlavv2 = zero
            endif
         elseif( dthx .eq. zero ) then
            if( X .eq. zero ) then
               r2 = half*aY
               ty = two/dthy
               d1 = zero
               do 10 j = 1, NG
                  yj = r2*(one + z(j))
                  a2 = ty*yj
                  d1 = d1 + w(j)*(one-z(j))*exp(-a2)
  10           continue
               dlavv2 = -half*var*d1
            else
               dlavv2 = zero
            endif
         elseif( dthy .eq. zero ) then
            if( Y .eq. zero ) then
               r1 = half*aX
               tx = two/dthx
               d1 = zero
               do 20 i = 1, NG
                  xi = r1*(one + z(i))
                  a1 = tx*xi
                  d1 = d1 + w(i)*(one-z(i))*exp(-a1)
  20           continue
               dlavv2 = -half*var*d1
            else
               dlavv2 = zero
            endif
         else
            r1 = half*aX
            r2 = half*aY
            tx = two/dthx
            ty = two/dthy
            d1 = zero
            do 40 i = 1, NG
               xi = r1*(one + z(i))
               a1 = tx*xi
               d2 = zero
               do 30 j = 1, NG
                  yj = r2*(one + z(j))
                  a2 = ty*yj
                  T  = sqrt(a1*a1 + a2*a2)
                  d2 = d2 + w(j)*(one-z(j))*exp(-T)
  30           continue
               d1 = d1 + w(i)*(one-z(i))*d2
  40        continue
            dlavv2 = -quart*var*d1
         endif
!				this is the old var func approximation
!				which isn't a particularly great approx.
!        twothd = -two/three
!        r  = abs(X)
!        s  = abs(Y)
!        xx = r/dthx
!        yy = s/dthy
!        pp = pb*pb
!        a1 = dthy*(pb + (one - pb)*exp(-xx*xx/pp))
!        a2 = dthx*(pb + (one - pb)*exp(-yy*yy/pp))
!        t1 = ( (one+yy**onept5)*(one+(r/a2)**onept5) )**twothd
!        t2 = ( (one+xx**onept5)*(one+(s/a1)**onept5) )**twothd
!        dlavv2 = half*var*(t1 + t2)
      else					! var > 0, return covariance
         if( (dthx .eq. zero) .and. (dthy .eq. zero) ) then
            if( (X .eq. zero) .and. (Y .eq. zero) ) then
               dlavv2 = var
            else
               dlavv2 = zero
            endif
         elseif( dthx .eq. zero ) then
            if( X .eq. zero ) then
               dlavv2 = var*exp(-two*aY/dthy)
            else
               dlavv2 = zero
            endif
         elseif( dthy .eq. zero ) then
            if( Y .eq. zero ) then
               dlavv2 = var*exp(-two*aX/dthx)
            else
               dlavv2 = zero
            endif
         else
            a1 = two*aX/dthx
            a2 = two*aY/dthy
            T  = sqrt(a1*a1 + a2*a2)
            dlavv2 = var*exp( -T )
         endif
      endif

      return
      end
