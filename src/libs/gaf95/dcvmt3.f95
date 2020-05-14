!  *********************************************************************
!  *                                                                   *
!  *                       subroutine dcvmt3                           *
!  *                                                                   *
!  *********************************************************************
!  Double Precision Version 1.0
!  Written by Gordon A. Fenton, TUNS, July 17, 1992
!
!  PURPOSE   computes the covariance matrices between local averages needed
!            by LAS3G for all stages of the subdivision.
!
!  This routine computes the three covariance matrices, R (27x27), B (8x8),
!  and S (27x8), required by LAS3G. The matrices R and B represent the
!  covariances between a 3 x 3 x 3 array of equal sized local average elements
!  and a 2 x 2 x 2 subarray respectively, numbered as follows for the first
!  plane only (in the depth or Z direction, the numbering is similar; ie the
!  element behind element 1 of the R array is element number 10, etc.);
!
!              |-- Dx--|   R Array                         B Array
!        ---   -------------------------           -------------------------
!         |    |       |       |       |           |       |       |       |
!         Dy   |   7   |   8   |   9   |           |       |       |       |
!         |    |       |       |       |           |       |       |       |
!        ---   -------------------------           -------------------------
!              |       |       |       |           |       |       |       |
!              |   4   |   5   |   6   |           |   3   |   4   |       |
!              |       |       |       |           |       |       |       |
!              -------------------------           -------------------------
!              |       |       |       |           |       |       |       |
!              |   1   |   2   |   3   |           |   1   |   2   |       |
!              |       |       |       |           |       |       |       |
!              -------------------------           -------------------------
!
!                                        Figure 1
!
!  Notice that the array B is just equal to selected elements of the array R
!  (ie. B(1,1) = R(1,1), B(1,4) = R(1,5), etc.). S represents the 27x8
!  covariance matrix between a doubly large array and an interior subdivided
!  2 x 2 x 2 array as shown, where the numbering of the larger array is the
!  same as on the left in Figure 1.
!
!              |--2Dx--|
!        ---   -------------------------
!         |    |       |       |       |
!        2Dy   |       |       |       |
!         |    |       |       |       |
!        ---   -------------------------
!              |       | 3 | 4 |       |
!              |       |-------|       |                     Figure 2
!              |       | 1 | 2 |       |
!              -------------------------
!              |       |       |       |
!              |       |       |       |
!              |       |       |       |
!              -------------------------
!
!  Finally note that the random process is assumed to be quadrant symmetric
!  (that is the covariance function is even in each individual lag component)
!  to reduce the total number of covariance calculations.
!
!  Arguments to this routine are as follows;
!
!    dvfn       external double precision function which returns the
!               variance of a local average of size X x Y x Z. This function
!               is referenced using the call "V = dvfn(X,Y,Z)" (other
!               parameters to the function must be passed via common blocks).
!
!    R          double precision array of size at least 27 x 27 which on output
!               will contain the covariance matrix between the local average
!               elements shown above in Figure 1, left. Note that since R is
!               symmetric, only the upper triangular values are actually
!               stored (this is to be compatible with DFACLD). (output)
!
!    ir         leading dimension of the array R as specified in the calling
!               routine. (input)
!
!    B          double precision vector of length at least 8 which on output
!               will contain the essential elements of the covariance matrix
!               between the local average elements shown above in Figure 1,
!               right (upper triangle only). (output)
!
!    S          double precision array of size at least 27 x 8 which on output
!               will contain the covariances between the set of subdivided
!               central elements and the parent 3 x 3 x 3 set of doubly
!               large elements (see Figure 2). (output)
!
!    is         leading dimension of the array S as specified in the calling
!               routine. (input)
!
!    Dx         x-dimension of the subdivided elements. The dimension of the
!               parent elements is assumed to be 2*Dx. (input)
!
!    Dy         y-dimension of the subdivided elements. The dimension of the
!               parent elements is assumed to be 2*Dy. (input)
!
!    Dz         z-dimension of the subdivided elements. The dimension of the
!               parent elements is assumed to be 2*Dz. (input)
!
!    lformR     logical flag which is true if the matrix R is to be formed
!               (in some cases, only the smaller matrix B is required). (input)
!----------------------------------------------------------------------------
      subroutine dcvmt3( dvfn, R, ir, B, ib, S, is, Dx, Dy, Dz, lformR )
      implicit real*8 (a-h,o-z)
      dimension R(ir,*), B(ib,*), S(is,*)
      dimension T(3)
      integer map(5)
      logical lformR
      external dvfn
      data zero/0.d0/, half/0.5d0/, one/1.d0/, onept5/1.5d0/
      data two/2.d0/, twopt5/2.5d0/
      data T/ 1.5d0, 0.5d0, 2.5d0 /
      data map/1, -1, 0, -1, 2/
!					first find B's essential elements
      B(1,1) = dcvaa3(dvfn,Dx,Dy,Dz,zero,zero,zero)
      B(1,2) = dcvaa3(dvfn,Dx,Dy,Dz, one,zero,zero)
      B(1,3) = dcvaa3(dvfn,Dx,Dy,Dz,zero, one,zero)
      B(1,4) = dcvaa3(dvfn,Dx,Dy,Dz, one, one,zero)
      B(1,5) = dcvaa3(dvfn,Dx,Dy,Dz,zero,zero, one)
      B(1,6) = dcvaa3(dvfn,Dx,Dy,Dz, one,zero, one)
      B(1,7) = dcvaa3(dvfn,Dx,Dy,Dz,zero, one, one)
      B(1,8) = dcvaa3(dvfn,Dx,Dy,Dz, one, one, one)
!						fill in the upper triangle
      do 20 n = 2, 8
         n3 = (n-1)/4
         n2 = (n - 4*n3 - 1)/2
         do 10 m = 2, n
            m3 = (m-1)/4
            m2 = (m - 4*m3 - 1)/2
            i3 = m3 - n3
            i2 = m2 - n2
            i1 = (m-n) - 2*i2 - 4*i3
            k  = 1 + iabs(i1) + 2*iabs(i2) + 4*iabs(i3)
            B(m,n) = B(1,k)
  10     continue
  20  continue

!					now find R optionally (upper triangle)
      if( lformR ) then
!						essential elements first ...
         R(1, 1) = B(1,1)
         R(1, 2) = B(1,2)
         R(1, 3) = dcvaa3(dvfn,Dx,Dy,Dz,  two, zero, zero )
         R(1, 4) = B(1,3)
         R(1, 5) = B(1,4)
         R(1, 6) = dcvaa3(dvfn,Dx,Dy,Dz,  two,  one, zero )
         R(1, 7) = dcvaa3(dvfn,Dx,Dy,Dz, zero,  two, zero )
         R(1, 8) = dcvaa3(dvfn,Dx,Dy,Dz,  one,  two, zero )
         R(1, 9) = dcvaa3(dvfn,Dx,Dy,Dz,  two,  two, zero )
         R(1,10) = B(1,5)
         R(1,11) = B(1,6)
         R(1,12) = dcvaa3(dvfn,Dx,Dy,Dz,  two, zero,  one )
         R(1,13) = B(1,7)
         R(1,14) = B(1,8)
         R(1,15) = dcvaa3(dvfn,Dx,Dy,Dz,  two,  one,  one )
         R(1,16) = dcvaa3(dvfn,Dx,Dy,Dz, zero,  two,  one )
         R(1,17) = dcvaa3(dvfn,Dx,Dy,Dz,  one,  two,  one )
         R(1,18) = dcvaa3(dvfn,Dx,Dy,Dz,  two,  two,  one )
         R(1,19) = dcvaa3(dvfn,Dx,Dy,Dz, zero, zero,  two )
         R(1,20) = dcvaa3(dvfn,Dx,Dy,Dz,  one, zero,  two )
         R(1,21) = dcvaa3(dvfn,Dx,Dy,Dz,  two, zero,  two )
         R(1,22) = dcvaa3(dvfn,Dx,Dy,Dz, zero,  one,  two )
         R(1,23) = dcvaa3(dvfn,Dx,Dy,Dz,  one,  one,  two )
         R(1,24) = dcvaa3(dvfn,Dx,Dy,Dz,  two,  one,  two )
         R(1,25) = dcvaa3(dvfn,Dx,Dy,Dz, zero,  two,  two )
         R(1,26) = dcvaa3(dvfn,Dx,Dy,Dz,  one,  two,  two )
         R(1,27) = dcvaa3(dvfn,Dx,Dy,Dz,  two,  two,  two )
!						distribute into upper triangle
         do 40 n = 2, 27
            n3 = (n-1)/9
            n2 = (n - 9*n3 - 1)/3
            do 30 m = 2, n
               m3 = (m-1)/9
               m2 = (m - 9*m3 - 1)/3
               i3 = m3 - n3
               i2 = m2 - n2
               i1 = (m-n) - 3*i2 - 9*i3
               k  = 1 + iabs(i1) + 3*iabs(i2) + 9*iabs(i3)
               R(m,n) = R(1,k)
  30        continue
  40     continue
      endif
!					now find S (essential elements first)
      ii = 1
      do 50 k = 1, 3
      do 50 j = 1, 3
         S(ii,  1) = dcvab3(dvfn,Dx,Dy,Dz,onept5,T(j),T(k))
         S(ii+1,1) = dcvab3(dvfn,Dx,Dy,Dz,half,  T(j),T(k))
         S(ii+2,1) = dcvab3(dvfn,Dx,Dy,Dz,twopt5,T(j),T(k))
         ii = ii + 3
  50  continue
!					and distribute to remainder of S
      do 70 n = 2, 8
         n3 = (n-1)/4
         n2 = (n - 4*n3 - 1)/2
         do 60 m = 1, 27
            m3 = (m-1)/9
            m2 = (m - 9*m3 - 1)/3
            i3 = iabs(4*m3 - 2*n3 - 3)
            i2 = iabs(4*m2 - 2*n2 - 3)
            i1 = iabs(4*m - 12*m2 - 36*m3 - 2*n + 4*n2 + 8*n3 - 5)
            k  = 1 + map(i1) + 3*map(i2) + 9*map(i3)
            S(m,n) = S(k,1)
  60     continue
  70  continue

      return
      end
