MODULE gaf95
!
INTERFACE
!
SUBROUTINE meshWT(g_coord,surf,g_num,argv,nlen,nosurf,ips)
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::g_coord(:,:),surf(:,:)
 INTEGER,INTENT(IN)::g_num(:,:),ips,nlen,nosurf
 CHARACTER(*),INTENT(IN)::argv
END SUBROUTINE meshWT
!
      subroutine las2g( Z, N1, N2, XL, YL, dvarfn, kseed, init, iout )
      parameter( MXM = 9, MXK = 256, NGS = 6912 )
      real Z(*)
      real C0((MXK*(MXK + 1))/2), U(NGS)
      real CT(6,2), CC(6,4,MXM), CS(6,4,MXM), CI(6,MXM)
      real AT(3,3,2), AC(4,3,4,MXM), AS(6,3,4,MXM), AI(9,3,MXM)
      end
!
      subroutine las2g2( Z, N1, N2, XL2, YL2, dvarfn, kseed, init, iout )
      parameter( MXM = 9, MXK = 256, NGS = 6912 )
      real Z(*)
      real C0((MXK*(MXK + 1))/2), U(NGS)
      real CT(6,2), CC(6,4,MXM), CS(6,4,MXM), CI(6,MXM)
      real AT(3,3,2), AC(4,3,4,MXM), AS(6,3,4,MXM), AI(9,3,MXM)
      end
!
      subroutine pltfld( job, sub1, sub2, Z, iz, N1, N2, D1, D2,title, ifld )
      real*8 z
      dimension Z(iz,1)
      character*(*) job, sub1, sub2, title
      end
!
      subroutine pltflds( job, sub1, sub2, Z, iz, N1, N2, D1, D2,title, ifld )
      real*4 z
      real*8 d1,d2
      dimension Z(iz)
      character*(*) job, sub1, sub2, title
      end
!
      subroutine las2i( dvarfn, N1, N2, XL2, YL2, kseed, MXM,C0, CT, CC, CS, CI, AT, AC, AS, AI, M,k1, k2, kk, iout, tol )
      parameter ( MXK = 256 )
      real C0(*), CT(6,*), CC(6,4,*), CS(6,4,*), CI(6,*)
      real AT(3,3,*), AC(4,3,4,*), AS(6,3,4,*), AI(9,3,*)
      real*8 R0(MXK*MXK)
      real*8 R(9,9,2), B(4,4), S(9,4)
      real*8 T1, T2, dvarfn, dble
      logical lformR
      integer mc(4,4), ms(6,4), mi(9)
      end
!
      subroutine las2i2( dvarfn, N1, N2, XL, YL, kseed, MXM,C0, CT, CC, CS, CI, AT, AC, AS, AI, M,k1, k2, kk, iout, tol )
      parameter ( MXK = 256 )
      real C0(*), CT(6,*), CC(6,4,*), CS(6,4,*), CI(6,*)
      real AT(3,3,*), AC(4,3,4,*), AS(6,3,4,*), AI(9,3,*)
      real*8 R0(MXK*MXK)
      real*8 R(9,9,2), B(4,4), S(9,4)
      real*8 T1, T2, dvarfn, dble
      logical lformR
      integer mc(4,4), ms(6,4), mi(9)
      end

!
      subroutine dcvit2( dvarfn, Q, iq, R, ir, k1, k2, Dx, Dy )
      implicit real*8 (a-h,o-z)
      dimension Q(iq,*), R(ir,*)
      integer kx(9), ky(9)
      end
!
      subroutine dcvmt2( dvarfn, R, ir, B, ib, S, is, Dx, Dy, lformR )
      implicit real*8 (a-h,o-z)
      dimension R(ir,*), B(ib,*), S(is,*)
      logical lformR
      end
!
      subroutine thin1d(R,ir,B,ib,S,is,AS,ias,AI,iai,CS,CI,i1,i2,i3,n,iout,tol)
      real AS(ias,n,*), AI(iai,*), CS(*), CI(*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RT(2,2), RI(3,3), DA(3), B1(7,7), B2(7,7)
      integer indT(2), indI(3)
      end
!
      subroutine corn2d(R,ir,B,ib,S,is,CC,ic,n,AC,mc,iout,tol)
      real CC(ic,*), AC(4,n,*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RC(4,4), DA(4), BB(7,7)
      integer mc(4,*), indx(4)
      end
!
      subroutine side2d(R,ir,B,ib,S,is,CS,ic,n,AS,ms,iout,tol)
      real CS(ic,*), AS(6,n,*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RS(6,6), DA(6), BB(7,7)
      integer ms(6,*), indx(6)
      end
!
      subroutine intr2d( R, ir, B, ib, S, is, CI, n, AI, mi, iout, tol )
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RI(9,9), DA(9), BB(7,7)
      real CI(*), AI(9,*)
      integer mi(*), indx(9)
      end
!
      subroutine daxpy(n,da,dx,dy)
      real*8 dx(*), dy(*), da, zero
      integer i, n
      end
!
      subroutine dchol2( A, ia, n, rerr )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*)
      real*4 rerr
      end
!
      subroutine dsifa(a,lda,n,kpvt,ierr)
      implicit real*8 (a-h,o-z)
      dimension a(lda,*)
      integer kpvt(*)
      logical swap
      end
!
      subroutine dsisl(a,lda,n,kpvt,b)
      real*8 a(lda,*),b(*)
      real*8 ak,akm1,bk,bkm1,ddot,denom,temp
      integer lda,n,kpvt(*)
      integer k,kp
      end
!
      subroutine  dswap (n,dx,dy)
      real*8 dx(*), dy(*), dtemp
      integer i, m, mp1, n
      end
!
      integer function iseed( kseed )
      integer kseed, getpid
      end
!
      integer function iseed2( kseed )
      integer kseed, getpid
      end
!
      subroutine vnorm( g, n )
      real g(*)
      end
!
      subroutine vnorm2( g, n )
      real g(*)
      end
!
      real*8 function ddot(n,dx,dy)
      real*8 dx(*),dy(*), zero
      integer i, n
      end
!
      integer function idamax( n, dx )
      real*8 dx(*), dmax
      integer i, n
      end
!
      integer function ntwom( N )
      integer N, k
      end
!
      real function randu(jseed)
      integer jseed, idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      PARAMETER (IM1  = 2147483563, &
                IM2  = 2147483399,  &
                AM   = 1./IM1,      &
                IMM1 = IM1-1,       &
                IA1  = 40014,       &
                IA2  = 40692,       &
                IQ1  = 53668,       &
                IQ2  = 52774,       &
                IR1  = 12211,       &
                IR2  = 3791,        & 
                NTAB = 32,          &
                NDIV = 1+IMM1/NTAB, &
                EPS  = 1.2e-7,      &
                RNMX = 1.-EPS)      
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum,idum2,ifirst
      end
!
      real function randu2(jseed)
      integer jseed, idumm,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      PARAMETER (IM1  = 2147483563, &
                IM2  = 2147483399,  &
                AM   = 1./IM1,      &
                IMM1 = IM1-1,       &
                IA1  = 40014,       &
                IA2  = 40692,       &
                IQ1  = 53668,       &
                IQ2  = 52774,       &
                IR1  = 12211,       &
                IR2  = 3791,        & 
                NTAB = 32,          &
                NDIV = 1+IMM1/NTAB, &
                EPS  = 1.2e-7,      &
                RNMX = 1.-EPS)      
      INTEGER idum22,j,k,ivv(NTAB),iyy
      SAVE ivv,iyy,idumm,idum22,ifirst2
      end
!
      real function second()
      real tyme(2)
      end
!
      real*8 function dcvaa2( cov, Dx, Dy, C1, C2 )
      parameter (NG = 16)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)
      end
!
      real*8 function dcvab2( cov, Dx, Dy, C1, C2 )
      parameter (NG = 16)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)
      end
!
      real*8 function dlafr2( X, Y )
      implicit real*8 (a-h,o-z)
      end
!
      real*8 function dlavx2( X, Y )
      parameter (NG = 16)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)
      end
!
      real*8 function dlavv2( X, Y )
      parameter (NG = 16)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)
      end
!
      real*8 function dlsep2( X, Y )
      implicit real*8 (a-h,o-z)
      end
!
      real*8 function dlsfr2( X, Y )
      implicit real*8 (a-h,o-z)
      end
!
      real*8 function dlspx2( X, Y )
      implicit real*8 (a-h,o-z)
      end
!
      real function gausv( var )
      real var, a, r
      logical getnxt
      save a, r
      end
!
      real*8 function uservf( x, y )
      implicit real*8 (a-h,o-z)
      end
!
      subroutine genfld(istat,pfld,p,nxe,nye,XL,YL,varfnc,kseed,liid)
      real pfld(nxe,*), p(*)
      real XL, YL
      integer istat, nxe, nye, kseed
      character*(*) varfnc
      logical liid
      real*8   dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      end
!
      subroutine genfld2(istat,pfld,p,nxe,nye,XL,YL,varfnc,kseed,liid)
      real pfld(nxe,*), p(*)
      real XL, YL
      integer istat, nxe, nye, kseed
      character*(*) varfnc
      logical liid
      real*8   dlavv2, dlsep2, dlspx2, dlafr2, dlsfr2
      end
!
      subroutine genprp(pfld,p,phifld,nrfx,nrfy,ltanfi)
      real pfld(nrfx,*), p(*), phifld(nrfx,*)
      integer nrfx, nrfy
      logical ltanfi
      end
!
      subroutine sim2ep(istat,iterm,verbos,cfld,phifld,psifld,efld, &
                       vfld,gamfld,k0fld,nrfx,nrfy,c,phi,psi,e,v,   &
                       gam,k0,R,ltanfi,lxfld,thx,thy,dx,dy,         &
                       dmpfld,nfld,jfld,ifld,job,sub1,sub2,varfnc,  &
                       kseed,debug)

      real cfld(nrfx,*), phifld(nrfx,*), psifld(nrfx,*)
      real efld(nrfx,*), vfld(nrfx,*), gamfld(nrfx,*), k0fld(nrfx,*)
      real c(*), phi(*), psi(*), e(*), v(*), gam(*), k0(*)
      real R(7,*), thx, thy
      integer istat, iterm, nfld, jfld, ifld, kseed
      integer nrfx, nrfy
      character*(*) job, sub1, sub2, varfnc
      logical verbos, dmpfld, debug, liid, shofld, lxfld, ltanfi
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy
      end
!
      subroutine propst(p)
      real p(*)
      end
!
      subroutine chol2( A, ia, n, rerr )
      dimension A(ia,*)
      real*4 rerr
      end
!
      subroutine chknxe(istat,iterm,verbos,nxe,nye,nrfx,nrfy)
      parameter (MXK = 256)
      logical verbos
      end
!
      subroutine print2( k, fmt, val1, val2 )
      parameter (MP = 2)
      character*(*) fmt
      real val1, val2, tmp(MP)
      end
!
      subroutine printv( k, fmt, val, n )
      parameter (MPX = 256)
      character*(*) fmt
      character*1 bslash, spc, tab
      real*4 val(*), tmp(MPX)
      integer*4 itmp(MPX)
      logical ldot, ldigit
      equivalence (itmp(1),tmp(1))
      end
!
      subroutine getfsp(fmt,lf,i,j,iw,id,ldot)
      character*(*) fmt
      character*1 c
      integer iq(2), il(2)
      logical ldot
      end
!
      subroutine prfmtf(val,iw,id,ldot,k)
      character fstr*256, d(10)*1
      logical ldot, lround
      end
!
      subroutine prfmte(val,iw,id,k)
      character fstr*256
      character d(10)*1
      end
!
      subroutine prfmti(ival,iw,k)
      character fstr*256, d(10)*1
      end
!
      integer function lnblnk( str )
      character*(*) str
      character*1 space, tab, null, char
      end
!
      logical function ldigit(c)
      character*1 c
      integer ic
      end
!
      subroutine permprop (perm,nrfx,nxe,nye,ndim,dir,prop)
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      real perm(nrfx,*)
      character*(*) dir
      integer nxe,nye,iel,iq,ip,ndim,nels,nrfx
      REAL(iwp) prop(nxe*nye)
      end
!
      subroutine msla(dx,dy,thx,thy,kmn,ksd,varfnc,istat,kmne,ksde)
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy, oned
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2, uservf
      real*8 dmin1,dble
      real kmn, kmne, ksd, ksde
      character*(*) varfnc
      end
!
      subroutine exca(dx,dy,thx,thy,varfnc,vf)
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy, oned
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2, uservf
      real*8 dmin1,dble
      character*(*) varfnc
      end
!
      subroutine exca2(dx,dy,thx,thy,varfnc,vf)
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy, oned
      real*8 dlavv2, dlsep2, dlspx2, dlafr2, dlsfr2, uservf
      real*8 dmin1,dble
      character*(*) varfnc
      end
!
      real*8 function dcvaa1( cov, Dx, C1 )
      parameter (NG = 16)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)
      end
!
      real*8 function dcvab1( cov, Dx, C1 )
      parameter (NG = 16)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)
      end
!
      subroutine dcvit1( vfn, R, S, k, NBH, T )
      implicit real*8 (a-h,o-z)
      dimension R(*), S(*)
      end
!
      real*8 function dlace1( T )
      implicit real*8 (a-h,o-z)
      end
!
      real*8 function dlafr1( X )
      implicit real*8 (a-h,o-z)
      end
!
      real*8 function dlavx1( T )
      implicit real*8 (a-h,o-z)
!      common/dparam/ var, dpb, dthx, dthy, dthz
      end
!
      real*8 function dlsmp1( T )
      implicit real*8 (a-h,o-z)
      end
!
      real*8 function dlspx1( T )
      implicit real*8 (a-h,o-z)
      end
!
      subroutine las1g( z, n, xl, dvarfn, nbh, kseed, init, iout )
      implicit real*8 (d)
      parameter( MXM = 16, MXK = 256, NGS = 4096 )
      real z(*), A(9,MXM), C(3,MXM), C0(MXK,MXK), g(NGS), ge(4*MXM)
      external dvarfn
      save A, C, C0, m, k1
!      common/LASTYM/ ti, ts
      end
!
     subroutine las1i( dvarfn, n, xl, MXM, MXK, NBH, A, C, c0,m, k1, iout, tol )
      implicit real*8 (d)
      parameter( MXmxk = 256, MXnbh = 5 )
      real A(9,*), C(3,*), c0(MXK,*)
      real*8 dg(MXnbh*MXnbh), dc0(MXmxk,MXmxk)
      real*8 daf(MXnbh), dr(MXmxk), ds(MXnbh)
      integer indx(MXnbh)
      external dvarfn
      end
!
      subroutine print1( k, fmt, val )
      character*(*) fmt
      real val
      end
!
      real*8 function derf(x)
      implicit real*8 (a-h,o-z)
      end
!
      real*8 function dphi(z)
      implicit real*8 (a-h,o-z)
      parameter (NPROBS = 768)
      dimension p(NPROBS)
      real*8 nine
      end
!      
      real*8 function dphinv(q)
      implicit real*8 (a-h,o-z)
      parameter (NPROBS = 768)
      dimension p(NPROBS)
      end
!
      subroutine sim1dm1(iterm,m,xl,nxe,thx,varfnc,kseed,lxfld,            &
                       a,ac,R,G,MXM,MXQ,lsquare,debug)
      real xl, thx, a(7,*), ac(2,nxe,*),R(MXM,*), G(MXQ,*)
      integer m, nxe, kseed, MXM, MXQ
      logical lxfld, lsquare, debug
      character*(*) varfnc
      real*4 rerr
      real*8 dvar, dpb, dthx, dthy, dthz, ddx
      real*8 dlace1, dlafr1, dlavx1, dlsmp1, dlspx1
      real*8 dmin1, dble
      end
!
      subroutine sim1dm2(iterm,m,xl,nxe,thx,varfnc,kseed,lxfld,            &
                       a,ac,R,G,MXM,MXQ,lsquare,debug)
      real xl, thx, a(7,*), ac(2,nxe,*),R(MXM,*), G(MXQ,*)
      integer m, nxe, kseed, MXM, MXQ
      logical lxfld, lsquare, debug
      character*(*) varfnc
      real*4 rerr
      real*8 dvar, dpb, dthx, dthy, dthz, ddx
      real*8 dlace1, dlafr1, dlavx1, dlsmp1, dlspx1
      real*8 dmin1, dble
      end
!
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
      END
!
      subroutine sim2bc1(istat,iterm,verbos,kxfld,kyfld,efld,  &
                        nrfx,nrfy,kx,ky,e,R,lxfld,thx,thy,dx,dy,    &
                        dmpfld,nfld,jfld,ifld,job,sub1,sub2,              &
                        varfnc,kseed,debug)

      real kxfld(nrfx,*), kyfld(nrfx,*), efld(nrfx,*)
      real kx(*), ky(*), e(*), R(3,*), thx, thy
      integer nfld, ifld, nrfx, nrfy
      character*(*) job, sub1, sub2, varfnc
      logical verbos, dmpfld, debug, liid, shofld, lxfld
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dmin1, dble
      save XL, YL, ienter, liid
      external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
!      common/dparam/ dvar, dpb, dthx, dthy, dthz
      end
!
      subroutine sim2bc2(istat,iterm,verbos,kxfld,efld,  &
                        nrfx,nrfy,kx,e,R,lxfld,thx,thy,dx,dy,    &
                        dmpfld,nfld,jfld,ifld,job,sub1,sub2,              &
                        varfnc,kseed,debug)

      real kxfld(nrfx,*),  efld(nrfx,*)
      real kx(*), e(*), R(2,*), thx, thy
      integer nfld, ifld, nrfx, nrfy
      character*(*) job, sub1, sub2, varfnc
      logical verbos, dmpfld, debug, liid, shofld, lxfld
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dmin1, dble
      save XL, YL, ienter, liid
      external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
!      common/dparam/ dvar, dpb, dthx, dthy, dthz
      end
!
      subroutine sim2ep1(istat,iterm,verbos,phifld,  &
                       nrfx,nrfy,phi,   &
                       ltanfi,lxfld,thx,thy,dx,dy,         &
                       dmpfld,nfld,jfld,ifld,job,sub1,sub2,varfnc,  &
                       kseed,debug)

      real phifld(nrfx,*)
      real phi(*)
      real thx, thy
      integer istat, iterm, nfld, jfld, ifld, kseed
      integer nrfx, nrfy
      character*(*) job, sub1, sub2, varfnc
      logical verbos, dmpfld, debug, liid, shofld, lxfld, ltanfi
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy
      save XL, YL, ienter, liid
      external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
!      common/dparam/ dvar, dpb, dthx, dthy, dthz
      end
!
      subroutine sim2ep2(istat,iterm,verbos,phifld,  &
                       nrfx,nrfy,phi,   &
                       ltanfi,lxfld,thx,thy,dx,dy,         &
                       dmpfld,nfld,jfld,ifld,job,sub1,sub2,varfnc,  &
                       kseed,debug)

      real phifld(nrfx,*)
      real phi(*)
      real thx, thy
      integer istat, iterm, nfld, jfld, ifld, kseed
      integer nrfx, nrfy
      character*(*) job, sub1, sub2, varfnc
      logical verbos, dmpfld, debug, liid2, shofld, lxfld, ltanfi
      real*8 dlavv2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy
      save XL2, YL2, ienter2, liid2
      external dlavv2, dlsep2, dlspx2, dlafr2, dlsfr2
!      common/dparam2/ dvar, dpb, dthx, dthy, dthz
      end
!
      subroutine fem2rf(istat,iterm,verbos,nx1,nx2,ny1,ny2,ngrad,nels,    &
                       nrfx,nrfy,imap)
      integer imap(*)
      logical verbos
      end
!
      subroutine fem2rf2(istat,iterm,verbos,nx1,nx2,nx3,ny1,ny2,ny3,ngrad,ngrad1,nels,&
                       nrfx,nrfy,imap)
      integer imap(*),ny4
      logical verbos
      end
!
      subroutine rect(nx1,nx2,ny1,ny2,ngrad,dx,g_coord,nf,inf,            &
                     g_num,maxel,nels,nn,n)
      real g_coord(2,inf)
      integer nf(inf,2),g_num(8,maxel)
      real coord(8,2)
      integer num(8)
      end
!
      subroutine rect2(nx1,nx2,nx3,ny1,ny2,ny3,ngrad,ngrad1,dx,g_coord,nf,inf,&
                     g_num,maxel,nels,nn,n)
      real g_coord(2,inf)
      integer nf(inf,2),g_num(8,maxel),ny4
      real coord(8,2)
      integer num(8)
      end
!
      subroutine sim2sd1(istat,iterm,verbos,cfld,phifld,psifld,gamfld,     &
                       efld,vfld,c,phi,psi,gam,e,v,R,lxfld,thx,thy,       &
                       nrfx,nrfy,dx,dy,dmpfld,nfld,jfld,ifld,job,sub1,    &
                       sub2,varfnc,kseed,debug,ltanfi)

      real cfld(nrfx,*), phifld(nrfx,*), psifld(nrfx,*)
      real gamfld(nrfx,*), efld(nrfx,*), vfld(nrfx,*)
      real c(*), phi(*), psi(*), gam(*), e(*), v(*), R(6,*), thx, thy
      integer nrfx, nrfy, nfld, ifld
      character*(*) job, sub1, sub2, varfnc
      logical verbos, dmpfld, debug, liid, shofld, lxfld, ltanfi
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dmin1, dble
      save XL, YL, ienter, liid
      external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      common/dparam/ dvar, dpb, dthx, dthy, dthz
      end
!
      subroutine sim2sd2(istat,iterm,verbos,cfld,phifld,psifld,gamfld,     &
                       efld,vfld,c,phi,psi,gam,e,v,R,lxfld,thx,thy,       &
                       nrfx,nrfy,dx,dy,dmpfld,nfld,jfld,ifld,job,sub1,    &
                       sub2,varfnc,kseed,debug,ltanfi)

      real cfld(nrfx,*), phifld(nrfx,*), psifld(nrfx,*)
      real gamfld(nrfx,*), efld(nrfx,*), vfld(nrfx,*)
      real c(*), phi(*), psi(*), gam(*), e(*), v(*), R(6,*), thx, thy
      integer nrfx, nrfy, nfld, ifld
      character*(*) job, sub1, sub2, varfnc
      logical verbos, dmpfld, debug, liid2, shofld, lxfld, ltanfi
      real*8 dvar, dpb, dthx, dthy, dthz, ddx, ddy
      real*8 dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dmin1, dble
      save XL2, YL2, ienter2, liid2
      external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      common/dparam/ dvar, dpb, dthx, dthy, dthz
      end
!
      subroutine dismsh2f(loads,nf,inf,g_coord,g_num,n,nn,nod,nels,       &
                        base,ib,ips,prop,imap,lmesh,lgrey,llog,           &
                        xsize,ysize,xoff,yoff,lfail,propb,nrfx,ny1,ny2)
      real prop(*),xload,propb(*),g_coord(2,*)
      real*8 loads(*)
      integer g_num(8,*),nf(nn,*), imap(*)
      integer nrfx,ny1,ny2
      character base*(*), ofile*256
      logical lmesh, lgrey, llog, lfail, ldet, ldetb
      end
!
      subroutine dismsh3f(loads,nf,inf,g_coord,g_num,n,nn,nod,nels,       &
                        base,ib,ips,prop,imap,lmesh,lgrey,llog,           &
                        xsize,ysize,xoff,yoff,lfail,propb,nrfx,ny1,ny2)
      real prop(*),xload,propb(*),g_coord(2,*)
      real*8 loads(*)
      integer g_num(8,*),nf(nn,*), imap(*)
      integer nrfx,ny1,ny2
      character base*(*), ofile*256
      logical lmesh, lgrey, llog, lfail, ldet, ldetb
      end
!
      subroutine dismsh4f(loads,nf,inf,g_coord,g_num,n,nn,nod,nels,       &
                        base,ib,ips,vf,c,prop,imap,lmesh,lgrey,llog,           &
                        xsize,ysize,xoff,yoff,lfail,propb,nrfx,ny1,ny2)
      real prop(*),xload,propb(*),g_coord(2,*),c(*),vf
      real*8 loads(*)
      integer g_num(8,*),nf(nn,*), imap(*)
      integer nrfx,ny1,ny2
      character base*(*), ofile*256
      logical lmesh, lgrey, llog, lfail, ldet, ldetb
      end
!
      subroutine print4( k, fmt, val1, val2, val3, val4 )
      parameter (MP = 4)
      character*(*) fmt
      real val1, val2, val3, val4, tmp(MP)
      end
!
      subroutine chkne3(istat,iterm,verbos,nxe,nye,nze,nrfx,nrfy,nrfz)
      parameter (MXK = 512)
      logical verbos
      end
!
      subroutine corn3d(R,ir,B,ib,S,is,CC,AC,iout,tol)
      real CC(28,*), AC(8,7,*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RC(8,8), DA(8), BB(7,7)
      integer ic(8,8), indx(8)
      end
!
      real*8 function dcvaa3( cov, Dx, Dy, Dz, C1, C2, C3 )
      parameter (NG = 20)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)
      external cov
      common/dparam/ var, pb, dthx, dthy, dthz
      end
!
      real*8 function dcvab3( cov, Dx, Dy, Dz, C1, C2, C3 )
      parameter (NG = 20)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)
      external cov
      common/dparam/ var, pb, dthx, dthy, dthz
      end
!     
      subroutine dcvit3( dvarfn, R0, iq, R, ir, k1, k2, k3, Dx, Dy, Dz )
      implicit real*8 (a-h,o-z)
      dimension R0(iq,*), R(ir,*)
      external dvarfn
      end
!
      subroutine dcvmt3( dvfn, R, ir, B, ib, S, is, Dx, Dy, Dz, lformR )
      implicit real*8 (a-h,o-z)
      dimension R(ir,*), B(ib,*), S(is,*)
      dimension T(3)
      integer map(5)
      logical lformR
      external dvfn
      end
!
      real*8 function dlafs3( X, Y, Z )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, H, G, da
      end
!
      real*8 function dlavx3( X, Y, Z )
      parameter (NG = 20)
      implicit real*8 (a-h,o-z)
      dimension w(NG), q(NG)
      common/dparam/ var, dpb, dthx, dthy, dthz
      end
!
      real*8 function dlsep3( X, Y, Z )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, dpb, dthx, dthy, dthz
      end
!
      real*8 function dlsfr3( X, Y, Z )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, H, G, F
      end
!
      real*8 function dlspx3( X, Y, Z )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, thx, thy, thz
      end
!
      subroutine edge3d(R,ir,B,ib,S,is,CE,AE,iout,tol)
      real CE(28,*), AE(12,7,*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RS(12,12), DA(12), BB(7,7)
      integer ic(12,12), indx(12)
      end
!
      subroutine fem3rf(istat,iterm,verbos,nx1,nx2,ny1,ny2,nze,ngrad,nels,    &
                       nrfx,nrfy,nrfz,imap)
      integer imap(*)
      integer iii
      logical verbos
      end
!
      subroutine intr3d( R, ir, B, ib, S, is, CI, AI, iout, tol )
      real*8 R(ir,*), B(ib,*), S(is,*), DA(27)
      real CI(*), AI(27,*)
      integer indx(27)
      end
!
      subroutine las3g( Z, N1, N2, N3, XL, YL, ZL, dvfn, kseed, init,     &
                       iout )
      parameter( MXM = 6, MXK = 512, NGS = 32000 )
      real Z(*)
      real C0(MXK*(MXK + 1)/2), U(NGS)
      real CC(28,8,MXM), CE(28,12,MXM), CS(28,6,MXM), CI(28,MXM)
      real AC(8,7,8,MXM), AE(12,7,12,MXM), AS(18,7,6,MXM), AI(27,7,MXM)
      real ATC(4,7,4), ATS(6,7,4), ATI(9,7)
      real CTC(28,4), CTS(28,4), CTI(28)
      save C0, CC, CS, CE, CI, AC, AE, AS, AI
      save ATC, ATS, ATI, CTC, CTS, CTI
      save M, k1, k2, k3, kk, NN
      external dvfn, dot1, dot2, dot3, dot4, dot5, dot6, dot7
      common/LASTYM/ ti, ts
      end
!
      subroutine las3i( dvfn, N1, N2, N3, XL, YL, ZL, kseed, MXM,         &
                       C0, CC, CE, CS, CI, AC, AE, AS, AI,                &
                       ATC, ATS, ATI, CTC, CTS, CTI,                      & 
                       M, k1, k2, k3, kk, iout, tol )
      parameter ( MXK = 512 )
      real C0(*), CC(28,8,*), CE(28,12,*), CS(28,6,*), CI(28,*)
      real AC(8,7,8,*), AE(12,7,12,*), AS(18,7,6,*), AI(27,7,*)
      real ATC(4,7,*), ATS(6,7,*), ATI(9,*)
      real CTC(28,*), CTS(28,*), CTI(*)
      real*8 R0(MXK*MXK)
      real*8 R(27,27,2), B(8,8), S(27,8)
      real*8 T1, T2, T3, dvfn, dble
      logical lformR, lk1, lk2, lk3
      integer mc(4,4,3), ms(6,4,3), mi(9,3)
      external dvfn
      end
!
      subroutine plan3d(Z,iq,jq,k1,k2,k3,AC,AS,AI,CC,CS,CI,U,iout)
      real Z(*), U(*)
      real AC(4,7,*), AS(6,7,*), AI(9,7)
      real CC(28,*), CS(28,*), CI(*)
      logical lk1, lk2, lk3
      real     dot1, dot2, dot3, dot4, dot5, dot6, dot7
      external dot1, dot2, dot3, dot4, dot5, dot6, dot7
      end
!
      subroutine side3d(R,ir,B,ib,S,is,CS,AS,iout,tol)
      real CS(28,*), AS(18,7,*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RS(18,18), DA(18), BB(7,7)
      integer ic(18,6), indx(18)
      end
!
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
      end
!
      real function dot8c( Z, A, C, U, dotn, j0, j1, j2, j3 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn
      end

      real function dot12h( Z, A, C, U, dotn, j0, j1, j2, j3 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn
      end
!
      real function dot12v( Z, A, C, U, dotn, j0, j1, j2, j3, j4, j5 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn
      end
!
      real function dot18a( Z, A, C, U, dotn, j0,j1,j2,j3,j4,j5 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn
      end
!
      real function dot18s( Z, A, C, U, dotn,j0,j1,j2,j3,j4,j5,j6,j7,j8)
      dimension Z(*), A(*), C(*), U(*)
      external dotn
      end
!
      real function dot27i( Z, A, C, U, dotn,j0,j1,j2,j3,j4,j5,j6,j7,j8)
      dimension Z(*), A(*), C(*), U(*)
      external dotn
      end
!
      real function dot1( X, Y )
      dimension X(*), Y(*)
      end
!
      real function dot2( X, Y )
      dimension X(*), Y(*)
      end
!
      real function dot3( X, Y )
      dimension X(*), Y(*)
      end
!
      real function dot4( X, Y )
      dimension X(*), Y(*)
      end
!
      real function dot5( X, Y )
      dimension X(*), Y(*)
      end
!
      real function dot6( X, Y )
      dimension X(*), Y(*)
      end
!
      real function dot7( X, Y )
      dimension X(*), Y(*)
      end
!
      real function dot4c( Z, A, C, U, dotn, j0, j1 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn
      end
!
      real function dot6h( Z, A, C, U, dotn, j0, j1 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn
      end
!
      real function dot6v( Z, A, C, U, dotn, j0, j1, j2 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn
      end
!
      real function dot9i( Z, A, C, U, dotn, j0, j1, j2 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn
      end
!
      subroutine print3( k, fmt, val1, val2, val3 )
      parameter (MP = 3)
      character*(*) fmt
      real val1, val2, val3, tmp(MP)
     end
!
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
      end
!
      real function betaf(x,a,b)
      parameter (IDEBUG = 1, IMAX = 200, TOL = 1.e-6 )
      end
!
      real function gamln(z)
      parameter (IDEBUG = 1, MXFAC = 34)
      dimension cof(6), fac(MXFAC)
      end
!
      real function gmdst( x, a, b )
      parameter (ITMAX = 400)
	  end
!
      real function gammaf( zz )
      parameter (IDEBUG = 1)
      dimension c(26)
      end
!
      subroutine rect2mod(nx1,nx2,nx3,ny1,ny2,ny3,rgrad,rgrad1,dx,g_coord,nf,inf,            &
                     g_num,maxel,nels,nn,n)
      real g_coord(2,inf)
      integer nf(inf,2),g_num(8,maxel),ny4
      real coord(8,2)
      integer num(8)
	  end
!
END INTERFACE
!
END MODULE gaf95
