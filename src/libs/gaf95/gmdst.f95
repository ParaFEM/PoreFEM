!  ********************************************************************
!  *                                                                  *
!  *                          function gmdst                          *
!  *                                                                  *
!  ********************************************************************
!  Single Precision Version 2.0
!  Written by Press etal., Numerical Recipes with some modifications
!  by Gordon A. Fenton, TUNS, 1994
!  Latest Update: May 9, 1994
!
!  PURPOSE  returns the Gamma distribution function for given x and parameters
!           alpha and beta.
!
!  This routine employs a series or partial fraction solution to compute
!  the cumulative probability associated with a variable following the
!  Gamma distribution,
!
!                   1       x
!     F(x) =   ---------  INT exp{-t/b}*t**(a-1) dt
!              G(a)*b**a    0
!
!  where G(.) is the Gamma function (G(a) = (a-1)! if `a' is an integer).
!  When b = 1, this is the incomplete Gamma function. When b = 2 and a = k/2,
!  this is the Chi-squared distribution with k degrees of freedom.
!  This routine is derived from Press etal. in "Numerical Recipes",
!  2nd Edition, pg. 218.
!  Arguments to the routine are as follows;
!
!     x   real value giving the position at which the cumulative distribution
!         is desired. (input)
!
!     a   the alpha parameter of the Gamma distribution (a > 0). (input)
!
!     b   the beta parameter of the Gamma distribution (b > 0). (input)
!---------------------------------------------------------------------------
      real function gmdst( x, a, b )
      parameter (ITMAX = 400)
      data zero/0.0/, one/1.0/, two/2.0/, tol/1.e-6/, small/1.e-24/

   1  format(a)
!-------------------------------- start executable statements ------------

      if( (x .le. zero) .or. (a .le. zero) .or. (b .lt. zero) ) then
         gmdst = zero
         return
      endif

      xb = x/b
      gl = gamln(a)
      ta = exp(a*alog(xb) - xb - gl)
      if( xb .lt. (a + one) ) then
!						use series representation
         ap = a
         d  = one/a
         s  = d
         do 10 i = 1, ITMAX
            ap = ap + one
            d  = d*xb/ap
            s  = s + d
            if( d .lt. s*tol ) then
               gmdst = s*ta
               return
            endif
  10     continue
      else
!						use continued fraction rep.
         ap = xb + one - a
         c  = one/small
         d  = one/ap
         h  = d
         do 20 i = 1, ITMAX
            an = -float(i)*(float(i) - a)
            ap = ap + two
            d  = an*d + ap
            if( abs(d) .lt. small ) d = small
            c  = ap + an/c
            if( abs(c) .lt. small ) c = small
            d  = one/d
            e  = d*c
            h  = h*e
            if( abs(e-one) .lt. tol ) then
               gmdst = one - h*ta
               return
            endif
  20     continue
      endif

      write(6,1)'Gamma Distribution (gmdst) unable to converge.'
      gmdst = zero
      return
      end
