!  *********************************************************************
!  *                                                                   *
!  *                       subroutine dcvit                            *
!  *                                                                   *
!  *********************************************************************
!  Double Precision Version 1.31
!  Written by Gordon A. Fenton, TUNS, Feb. 11, 1993
!  Latest Update: Apr 17, 2001
!
!  PURPOSE  computes the covariance vector between a sequence of
!           local averages in 1-D. Used by LAS1G.
!
!  This routine computes the covariances between a sequence of adjacent
!  equal-sized local averages of a 1-D random function. As well, the
!  covariances between the half cell 1' and the full cells 1, 2, ..., NBH
!  are computed. That is, for local average cells arranged as follows,
!
!    |--T--|
!     ------------------------------------
!    |  1  |  2  |  3  |  4  | ... |   k  |
!     ------------------------------------
!
!     ------------------------------------
!    |  1  |  2  |1'|2'|  4  | ... |   k  | (subdivided cells 1' and 2')
!     ------------------------------------
!
!  and where Z_i represents the local average of cell i, this routine
!  finds the vector {R} such that
!
!    R_j = E[ Z_1 * Z_j ]  ( = Cov[Z_1,Z_j] for Z zero-mean)
!
!  and the vector {S} such that
!
!    S_j = E[ Z_1' * Z_j]  ( = Cov[Z_1',Z_j] for Z zero-mean);
!
!  The actual location of the subdivided cell is given by the index (NBH+1)/2
!  (ie, NBH = 3 implies that the cells 1' and 2' occur in the full cell
!  number 2). Thus NBH must be odd.
!
!  Arguments to this routine are as follows;
!
!    vfn      external real*8 function which returns the covariance of the
!             random process between two points separated by a distance T.
!             VFN is referenced as follows
!
!                C = vfn(T)
!
!             where (T) is the separation distance. Any other parameters
!             to the function must be passed by common block from the
!             calling routine. The current version of LAS uses the sign on
!             `var' to tell vfn to return either the covariance, or the
!             variance of a local average over distance V1. The latter is
!             returned if var < 0, although this feature is no longer used
!             by LAS. VFN must be an even function, ie vfn(-T) = vfn(T).
!
!    R        real vector of length at least k which on output will
!             contain the covariances discussed above. (output)
!
!    S        real vector of length at least NBH which on output will contain
!             the covariances between the half cell 1' and the full cells
!             1, 2, ..., NBH. (output)
!
!    k        the number of adjacent cells considered in the calculation of
!             R. (input)
!
!    NBH      the cell which is subdivided in the calculation of {S} is
!             (NBH+1)/2. NBH must be odd. (input)
!
!    T        the cell dimension. (input)
!
!  REVISION HISTORY:
!  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!  1.2	now including a Gaussian Quadrature integration option as an
!	alternative to the variance function approach to evaluate
!	covariances between local averages. (Jun 16/00)
!  1.3	revised above docs to reflect elimination of lvarfn (Mar 27/01)
!  1.31	modified writeup above for new vfn return value. (Apr 17/01)
!---------------------------------------------------------------------------
      subroutine dcvit1( vfn, R, S, k, NBH, T )
      implicit real*8 (a-h,o-z)
      dimension R(*), S(*)
      external vfn
      data half/0.5d0/

!					compute cell 1 - cell i covariances
      do 10 i = 1, k
         C1   = dble(i-1)
         R(i) = dcvaa1(vfn,T,C1)
  10  continue
!					compute cell 1' - cell i covariances
      n1 = -1 - 2*NBH
      t2 = half*T
      do 20 i = 1, NBH
         C1   = half*dble(n1+4*i)
         S(i) = dcvab1(vfn,t2,C1)
  20  continue

      return
      end
