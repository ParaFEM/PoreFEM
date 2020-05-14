!  *********************************************************************
!  *                                                                   *
!  *                           function dcvab3                         *
!  *                                                                   *
!  *********************************************************************
!  Double Precision Version 1.5
!  Written by Gordon A. Fenton, TUNS, July 17, 1992
!  Latest Update: Jul 2, 2003
!
!  PURPOSE  returns the covariance between two 3-D local averages, one
!           having eight times the volume of the other. Used by LAS3G.
!
!  This function evaluates the covariance between two local averages in
!  3-dimensional space. One local average is assumed to be of size Dx x Dy x Dz
!  and the other of size 2*Dx x 2*Dy x 2*Dz and they are separated in space by
!  a central distance having components Tx = C1*Dx, Ty = C2*Dy, and
!  Tz = C3*Dz.
!
!  The covariance is obtained by a 16-pt Gaussian quadrature of a
!  6-D integral collapsed (by assuming the covariance function to be
!  quadrant symmetric) to a 3-D integral.
!
!  NOTE: if the covariance function is not quadrant symmetric, this
!        routine will return erroneous results. Quadrant symmetry means
!        that cov(X,Y,Z) = cov(-X,Y,Z) = cov(X,-Y,Z), etc, where
!        cov(X,Y,Z) is a function returning the covariance between
!        points separated by lag (X,Y,Z), as discussed next.
!
!  The variance/covariance function is referenced herein using a call of the
!  form
!
!          V = cov( X, Y, Z )
!
!  where X, Y, and Z are the lag distances between the points in the
!  field.
!  Parameters, such as var, pb, dthx, dthy, and dthz are passed to the
!  variance/covariance function and this routine via the common block
!  'dparam'.
!
!  Arguments to this function are as follows
!
!    cov       external function provided by the user. On each invocation,
!              this routine calls cov up to 27*(5^3) = 3375 times.
!
!    Dx        x-dimension of the smaller local average cell. (input)
!
!    Dy        y-dimension of the smaller local average cell. (input)
!
!    Dz        z-dimension of the smaller local average cell. (input)
!
!    C1        x-direction distance between local average centers is C1*Dx.
!              (input)
!
!    C2        y-direction distance between local average centers is C2*Dy.
!              (input)
!
!    C3        z-direction distance between local average centers is C3*Dz.
!              (input)
!
!  REVISION HISTORY:
!  1.1	now including a Gaussian Quadrature integration option as an
!	alternative to the variance function approach to evaluate
!	covariances between local averages. (Jun 16/00)
!  1.2	use Gauss quadrature unless intervals overlap, eliminated lvarfn
!	(Apr 5/01)
!  1.3	now use Gauss quadrature for everything. (May 8/01)
!  1.4	now use 16 point Gauss quadrature (Mar 20/03)
!  1.5	simplification below only correct for volume completely enclosed.
!	Now using 20-point Gauss quadrature (Jul 2/03)
!---------------------------------------------------------------------------
      real*8 function dcvab3( cov, Dx, Dy, Dz, C1, C2, C3 )
      parameter (NG = 20)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)

      external cov
      common/dparam/ var, pb, dthx, dthy, dthz

      data zero/0.d0/, half/0.5d0/, one/1.d0/, onept5/1.5d0/, two/2.d0/
      data eighth/0.125d0/, fv12th/0.001953125d0/	!(= 1/512)
      data halfe/0.500000000001d0/			! half plus a wee bit
      data w/0.017614007139152118312d0, 0.040601429800386941331d0,        &
             0.062672048334109063570d0, 0.083276741576704748725d0,        &
             0.101930119817240435037d0, 0.118194531961518417312d0,        &
             0.131688638449176626898d0, 0.142096109318382051329d0,        &
             0.149172986472603746788d0, 0.152753387130725850698d0,        &
             0.152753387130725850698d0, 0.149172986472603746788d0,        &
             0.142096109318382051329d0, 0.131688638449176626898d0,        &
             0.118194531961518417312d0, 0.101930119817240435037d0,        &
             0.083276741576704748725d0, 0.062672048334109063570d0,        &
             0.040601429800386941331d0, 0.017614007139152118312d0/
      data z/-.993128599185094924786d0, -.963971927277913791268d0,        &
             -.912234428251325905868d0, -.839116971822218823395d0,        &
             -.746331906460150792614d0, -.636053680726515025453d0,        &
             -.510867001950827098004d0, -.373706088715419560673d0,        &
             -.227785851141645078080d0, -.076526521133497333755d0,        &
             0.076526521133497333755d0, 0.227785851141645078080d0,        &
             0.373706088715419560673d0, 0.510867001950827098004d0,        &
             0.636053680726515025453d0, 0.746331906460150792614d0,        &
             0.839116971822218823395d0, 0.912234428251325905868d0,        &
             0.963971927277913791268d0, 0.993128599185094924786d0/

!					enclosed volume: GQ simplifies
!					to variance of larger volume

      if( (C1.lt.halfe).and.(C2.lt.halfe).and.(C3.lt.halfe) ) then
         d1 = zero
         do 30 i = 1, NG
            xi = Dx*(one + z(i))
            d2 = zero
            do 20 j = 1, NG
               yj = Dy*(one + z(j))
               d3 = zero
               do 10 k = 1, NG
                  zk = Dz*(one + z(k))
                  d3 = d3 + w(k)*(one - z(k))*cov(xi,yj,zk)
  10           continue
               d2 = d2 + w(j)*(one-z(j))*d3
  20        continue
            d1 = d1 + w(i)*(one-z(i))*d2
  30     continue
         dcvab3 = eighth*d1
         return
      endif
!					otherwise, non-enclosed volume
      r1 = half*Dx
      r2 = half*Dy
      r3 = half*Dz

      s2 = two*C1
      s1 = s2 - two
      s3 = s2 + two
      v2 = two*C2
      v1 = v2 - two
      v3 = v2 + two
      u2 = two*C3
      u1 = u2 - two
      u3 = u2 + two

      d1 = zero
      do 60 i = 1, NG
         x1  = r1*(z(i) + s1)
         x2  = r1*(z(i) + s2)
         x3  = r1*(z(i) + s3)
         d21 = zero
         d22 = zero
         d23 = zero
         do 50 j = 1, NG
            y1 = r2*(z(j) + v1)
            y2 = r2*(z(j) + v2)
            y3 = r2*(z(j) + v3)
            d31 = zero
            d32 = zero
            d33 = zero
            d34 = zero
            d35 = zero
            d36 = zero
            d37 = zero
            d38 = zero
            d39 = zero
            do 40 k = 1, NG
               z1 = r3*(z(k) + u1)
               z2 = r3*(z(k) + u2)
               z3 = r3*(z(k) + u3)
               dp = one + z(k)
               dm = one - z(k)
               d31 = d31 + w(k)*(dp*cov(x1,y1,z1) + two*cov(x1,y1,z2)     &
                                +dm*cov(x1,y1,z3))
               d32 = d32 + w(k)*(dp*cov(x1,y2,z1) + two*cov(x1,y2,z2)     &
                                +dm*cov(x1,y2,z3))
               d33 = d33 + w(k)*(dp*cov(x1,y3,z1) + two*cov(x1,y3,z2)     &
                                +dm*cov(x1,y3,z3))
               d34 = d34 + w(k)*(dp*cov(x2,y1,z1) + two*cov(x2,y1,z2)     &
                                +dm*cov(x2,y1,z3))
               d35 = d35 + w(k)*(dp*cov(x2,y2,z1) + two*cov(x2,y2,z2)     &
                                +dm*cov(x2,y2,z3)) 
               d36 = d36 + w(k)*(dp*cov(x2,y3,z1) + two*cov(x2,y3,z2)     &
                                +dm*cov(x2,y3,z3))
               d37 = d37 + w(k)*(dp*cov(x3,y1,z1) + two*cov(x3,y1,z2)     &
                                +dm*cov(x3,y1,z3))
               d38 = d38 + w(k)*(dp*cov(x3,y2,z1) + two*cov(x3,y2,z2)     &
                                +dm*cov(x3,y2,z3))
               d39 = d39 + w(k)*(dp*cov(x3,y3,z1) + two*cov(x3,y3,z2)     &
                                +dm*cov(x3,y3,z3))
  40        continue
            dp = one + z(j)
            dm = one - z(j)
            d21 = d21 + w(j)*(dp*d31 + two*d32 + dm*d33)
            d22 = d22 + w(j)*(dp*d34 + two*d35 + dm*d36)
            d23 = d23 + w(j)*(dp*d37 + two*d38 + dm*d39)
  50     continue
         d1 = d1 + w(i)*((one+z(i))*d21 + two*d22 + (one-z(i))*d23)
  60  continue
      dcvab3 = fv12th*d1

      return
      end
