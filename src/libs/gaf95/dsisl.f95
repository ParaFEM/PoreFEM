!  *******************************************************************
!  *                                                                 *
!  *                       subroutine dsisl                          *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 08/14/78 1.01
!  Written by James Bunch, Univ. Calif. San Diego, Argonne Nat. Lab.
!  Modified by Gordon A. Fenton, Aug. 24, 1993
!  Latest Update: Jun 9, 1999
!
!  PURPOSE   solves the double precision symmetric system [A]{X} = {B}
!            using the factors computed by DSIFA.
!
!     On entry
!
!        A       double precision(lda,n)
!                the output from dsifa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        kpvt    integer(n)
!                the pivot vector from dsifa.
!
!        B       double precision(n)
!                the right hand side vector.
!
!     on return
!
!        B       the solution vector  x .
!
!  Notes:
!   1) a division by zero may occur if DSIFA  has set ierr .ne. 0.
!
!  Requires:
!   1) from libGAFsim:	DAXPY, DDOT
!
!  REVISION HISTORY:
!  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!-------------------------------------------------------------------------
      subroutine dsisl(a,lda,n,kpvt,b)
      real*8 a(lda,*),b(*)
      real*8 ak,akm1,bk,bkm1,ddot,denom,temp
      integer lda,n,kpvt(*)
      integer k,kp
!							work backwards
      k = n
   10 if (k .ne. 0) then
         if (kpvt(k) .ge. 0) then
!							1 x 1 pivot block.
            if (k .ne. 1) then
               kp = kpvt(k)
               if (kp .ne. k) then
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
               endif
               call daxpy(k-1,b(k),a(1,k),b(1))
            endif
            b(k) = b(k)/a(k,k)
            k = k - 1
         else
!							2 x 2 pivot block.
            if (k .ne. 2) then
               kp = iabs(kpvt(k))
               if (kp .ne. k - 1) then
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
               endif
               call daxpy(k-2,b(k),a(1,k),b(1))
               call daxpy(k-2,b(k-1),a(1,k-1),b(1))
            endif
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/a(k-1,k)
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0d0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
         endif
         go to 10
      endif
!							now work forwards
      k = 1
   90 if (k .le. n) then
         if (kpvt(k) .ge. 0) then
!							1 x 1 pivot block.
            if (k .ne. 1) then
               b(k) = b(k) + ddot(k-1,a(1,k),b(1))
               kp = kpvt(k)
               if (kp .ne. k) then
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
               endif
            endif
            k = k + 1
         else
!							2 x 2 pivot block.
            if (k .ne. 1) then
               b(k) = b(k) + ddot(k-1,a(1,k),b(1))
               b(k+1) = b(k+1) + ddot(k-1,a(1,k+1),b(1))
               kp = iabs(kpvt(k))
               if (kp .ne. k) then
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
               endif
            endif
            k = k + 2
         endif
         go to 90
      endif

      return
      end
