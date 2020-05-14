!  *********************************************************************
!  *                                                                   *
!  *                  More Dot Product Functions                       *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.0
!  Written by Gordon A. Fenton, TUNS, Dec. 6, 1993
!
!  PURPOSE  another set of specialized dot product functions designed for
!           use by LAS3G.
!
!  This file includes a set of functions which compute dot products for
!  LAS3G of the form
!
!             dot = {A^T}{Z} + {C^T}{U}
!
!  as well as some simple fixed length dot products.
!  In general arguments are as follows;
!
!    Z    real vector containing the parent cell values. Z(1) is assumed
!         to be the first element included in the dot product. (input)
!
!    A    real vector containing the BLUE coefficients. (input)
!
!    C    real vector containing the Covariance coefficients. (input)
!
!    U    real vector containing the random noise inputs. (input)
!
!    dotn external simple dot product function of fixed length. This is
!         used to compute {C^T}{U}.
!
!    j0-j2
!         pointers to cells in left edge of parent cell neighborhood.
!---------------------------------------------------------------------------
!-------------------------------------------------- 4 element corner -------
      real function dot4c( Z, A, C, U, dotn, j0, j1 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot4c = A(1)*Z(j0) + A(2)*Z(1+j0)                                   &
           + A(3)*Z(j1) + A(4)*Z(1+j1)                                    &
           + dotn( C, U )

      return
      end
!-------------------------------------------------- 6 element side (hor) ---
      real function dot6h( Z, A, C, U, dotn, j0, j1 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot6h = A( 1)*Z(j0) + A( 2)*Z(1+j0) + A( 3)*Z(2+j0)                 &
            + A( 4)*Z(j1) + A( 5)*Z(1+j1) + A( 6)*Z(2+j1)                 &
            + dotn( C, U )

      return
      end
!-------------------------------------------------- 6 element side (vert) --
      real function dot6v( Z, A, C, U, dotn, j0, j1, j2 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot6v = A( 1)*Z(j0) + A( 2)*Z(1+j0) + A( 3)*Z(j1) + A( 4)*Z(1+j1)   &
            + A( 5)*Z(j2) + A( 6)*Z(1+j2)                                 &
            + dotn( C, U )

      return
      end
!-------------------------------------------------- 9 element interior -----
      real function dot9i( Z, A, C, U, dotn, j0, j1, j2 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot9i = A( 1)*Z(j0) + A( 2)*Z(1+j0) + A( 3)*Z(2+j0)                 &
            + A( 4)*Z(j1) + A( 5)*Z(1+j1) + A( 6)*Z(2+j1)                 &
            + A( 7)*Z(j2) + A( 8)*Z(1+j2) + A( 9)*Z(2+j2)                 &
            + dotn( C, U )

      return
      end
