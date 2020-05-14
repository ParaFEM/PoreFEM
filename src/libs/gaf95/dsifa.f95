!  *******************************************************************
!  *                                                                 *
!  *                       subroutine dsifa                          *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 08/14/78 1.01
!  Written by James Bunch, Univ. Calif. San Diego, Argonne Nat. Lab.
!  Modified by Gordon A. Fenton, Aug. 24, 1993.
!  Latest Update: Jun 9, 1999
!
!  PURPOSE  factors a double precision symmetric matrix by elimination
!           with symmetric pivoting
!
!     To solve  A*X = B ,          follow dsifa by dsisl.
!
!     On entry
!
!        A       double precision(lda,n)
!                the symmetric matrix to be factored.
!                only the diagonal and upper triangle are used.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        A       a block diagonal matrix and the multipliers which
!                were used to obtain it.
!                the factorization can be written  a = u*d*trans(u)
!                where  u  is a product of permutation and unit
!                upper triangular matrices , trans(u) is the
!               transpose of  u , and  d  is block diagonal
!               with 1 by 1 and 2 by 2 blocks.
!
!       kpvt    integer(n)
!               an integer vector of pivot indices.
!
!       ierr    integer
!               = 0  if all goes well.
!               = k  if the k-th pivot block is singular. this is
!                    not an error condition for this subroutine,
!                    but it does indicate that dsisl or dsidi may
!                    divide by zero if called.
!
! Requires:
!   1) from libGAFsim:	DAXPY, DSWAP, IDAMAX
!
! REVISION HISTORY:
! 1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!------------------------------------------------------------------------
      subroutine dsifa(a,lda,n,kpvt,ierr)
      implicit real*8 (a-h,o-z)
      dimension a(lda,*)
      integer kpvt(*)
      logical swap
      data zero/0.d0/, one/1.d0/, eight/8.d0/, svtn/17.d0/

!						choose pivot block size.
      alpha = (one + dsqrt(svtn))/eight
      ierr  = 0
!						loop on k, from n to 1.
      k = n
!						leave the loop if k=0 or k=1.
  10  if (k .eq. 0) return
      if (k .eq. 1) then
         kpvt(1) = 1
         if (a(1,1) .eq. zero) ierr = 1
         return
      endif
      km1 = k - 1
      absakk = dabs(a(k,k))
!						largest element in col k
      imax = idamax(k-1,a(1,k))
      colmax = dabs(a(imax,k))
      if (absakk .ge. alpha*colmax) then
         kstep = 1
         swap = .false.
      else
!						largest element in row imax
         rowmax = zero
         imaxp1 = imax + 1
         do 20 j = imaxp1, k
            rowmax = dmax1(rowmax,dabs(a(imax,j)))
  20     continue
         if (imax .ne. 1) then
            jmax = idamax(imax-1,a(1,imax))
            rowmax = dmax1(rowmax,dabs(a(jmax,imax)))
         endif
         if (dabs(a(imax,imax)) .ge. alpha*rowmax) then
            kstep = 1
            swap = .true.
         else
            if (absakk .ge. alpha*colmax*(colmax/rowmax)) then
               kstep = 1
               swap = .false.
            else
               kstep = 2
               swap = imax .ne. km1
            endif
         endif
      endif
      if (dmax1(absakk,colmax) .eq. zero) then
!						column k is zero.  set ierr
         kpvt(k) = k
         ierr = k
      elseif (kstep .eq. 1) then
!						1 x 1 pivot block.
         if (swap) then
            call dswap(imax,a(1,imax),a(1,k))
            do 30 jj = imax, k
               j = k + imax - jj
               t = a(j,k)
               a(j,k) = a(imax,j)
               a(imax,j) = t
  30        continue
         endif
!						perform the elimination.
         do 40 jj = 1, km1
            j = k - jj
            amk = -a(j,k)/a(k,k)
            t = amk
            call daxpy(j,t,a(1,k),a(1,j))
            a(j,k) = amk
  40     continue
!	     				set the pivot array.
         kpvt(k) = k
         if (swap) kpvt(k) = imax
      else
!						2 x 2 pivot block.
         if (swap) then
            call dswap(imax,a(1,imax),a(1,k-1))
            do 50 jj = imax, km1
               j = km1 + imax - jj
               t = a(j,k-1)
               a(j,k-1) = a(imax,j)
               a(imax,j) = t
  50        continue
            t = a(k-1,k)
            a(k-1,k) = a(imax,k)
            a(imax,k) = t
         endif
!						perform the elimination.
         km2 = k - 2
         if (km2 .ne. 0) then
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            denom = one - ak*akm1
            do 60 jj = 1, km2
               j = km1 - jj
               bk = a(j,k)/a(k-1,k)
               bkm1 = a(j,k-1)/a(k-1,k)
               amk = (akm1*bk - bkm1)/denom
               amkm1 = (ak*bkm1 - bk)/denom
               t = amk
               call daxpy(j,t,a(1,k),a(1,j))
               t = amkm1
               call daxpy(j,t,a(1,k-1),a(1,j))
               a(j,k) = amk
               a(j,k-1) = amkm1
  60        continue
         endif
!            				set the pivot array.
         kpvt(k) = 1 - k
         if (swap) kpvt(k) = -imax
         kpvt(k-1) = kpvt(k)
      endif
      k = k - kstep
      go to 10

      end
