!  *********************************************************************
!  *                                                                   *
!  *                          function randu                           *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 2.0
!  Adapted from Press etal. "Numerical Recipes in Fortran", 2nd Ed.,
!  RAN1 program, pg 271
!  by Gordon A. Fenton, TUNS, Oct 14, 1996
!
!  PURPOSE  returns a pseudo-random number uniformly distributed on the
!           interval (0,1)
!
!  This routine uses a pair of multiplicative congruential generators to
!  produce a sequence of random numbers uniformly distributed on the interval
!  (0,1) which exceeds Park and Miller's ``Minimal Standard.'' The generator
!  has a period of about 2.3E+18 and uses a shuffling algorithm due to Bays
!  and Durham to remove low-order serial correlations. This routine is
!  about 2.0 times slower than the minimal routine `RAN0' presented in
!  Numerical Recipes. The endpoints of the interval (0,1) are excluded.
!
!  Arguments to the function are as follows;
!
!      jseed    integer seed. If jseed is zero, the next pseudo random number
!               in the sequence is returned. Otherwise if jseed is positive,
!               it is used to initialize the generator and the first value
!               in the sequence is returned. (input)
!
!  Notes:
!    1) if you always call this routine without initializing it (using
!       a positive seed) then you'll always get the same sequence of
!       pseudo-random numbers.
!
!  REVISION HISTORY:
!  2.0	implemented RAN2 from 2nd Edition of Numerical Recipes
!-----------------------------------------------------------------------------
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
      DATA idum22/123456789/, ivv/NTAB*0/, iyy/0/, ifirst2/0/

      if( (jseed .gt. 0) .or. (ifirst2 .eq. 0) ) then
         ifirst2 = 1
         idumm   = max0(jseed,1)
         idum22  = idumm
         do 10 j = NTAB+8, 1, -1
            k    = idumm/IQ1
            idumm = IA1*(idumm-k*IQ1) - k*IR1
            if( idumm .lt. 0 ) idumm  = idumm + IM1
            if( j .le. NTAB ) ivv(j) = idumm
  10     continue
         iyy = ivv(1)
      endif

      k    = idumm/IQ1
      idumm = IA1*(idumm-k*IQ1) - k*IR1
      if( idumm .lt. 0 ) idumm = idumm + IM1
      k     = idum22/IQ2
      idum22 = IA2*(idum22-k*IQ2) - k*IR2
      if( idum22 .lt. 0 ) idum22 = idum22 + IM2
      j     = 1 + iyy/NDIV
      iyy    = ivv(j) - idum22
      ivv(j) = idumm
      if( iyy .lt. 1 ) iyy = iyy + IMM1
      randu2 = amin1( AM*iyy, RNMX )
      return
      end
