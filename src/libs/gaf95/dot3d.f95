!  *********************************************************************
!  *                                                                   *
!  *                       Dot Product Functions                       *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.0
!  Written by Gordon A. Fenton, TUNS, Aug. 2, 1993
!
!  PURPOSE  a set of specialized dot product functions designed for use by
!           LAS3G.
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
!    j0-j8
!         pointers to cells in left edge of parent cell neighborhood.
!---------------------------------------------------------------------------
!-------------------------------------------------- 8 element corner -------
      real function dot8c( Z, A, C, U, dotn, j0, j1, j2, j3 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn
!
      dot8c = A(1)*Z(j0) + A(2)*Z(1+j0)                                  &
           + A(3)*Z(j1) + A(4)*Z(1+j1)                                   &
           + A(5)*Z(j2) + A(6)*Z(1+j2)                                   &
           + A(7)*Z(j3) + A(8)*Z(1+j3)                                   
!          + dotn( C, U )
      dot8c = dot8c + dotn( C, U )

      return
      end
!-------------------------------------- 12 element edge (horizontal) -------
      real function dot12h( Z, A, C, U, dotn, j0, j1, j2, j3 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot12h = A( 1)*Z(j0) + A( 2)*Z(1+j0) + A( 3)*Z(2+j0)                &
            + A( 4)*Z(j1) + A( 5)*Z(1+j1) + A( 6)*Z(2+j1)                &
            + A( 7)*Z(j2) + A( 8)*Z(1+j2) + A( 9)*Z(2+j2)                &
            + A(10)*Z(j3) + A(11)*Z(1+j3) + A(12)*Z(2+j3)                &
            + dotn( C, U )

      return
      end
!-------------------------------------- 12 element edge (vertical) -------
      real function dot12v( Z, A, C, U, dotn, j0, j1, j2, j3, j4, j5 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot12v = A( 1)*Z(j0) + A( 2)*Z(1+j0) + A( 3)*Z(j1) + A( 4)*Z(1+j1)  &
            + A( 5)*Z(j2) + A( 6)*Z(1+j2) + A( 7)*Z(j3) + A( 8)*Z(1+j3)  &
            + A( 9)*Z(j4) + A(10)*Z(1+j4) + A(11)*Z(j5) + A(12)*Z(1+j5)  &
            + dotn( C, U )

      return
      end
!-------------------------------------- 18 element side (front,back,top,bot) -
      real function dot18a( Z, A, C, U, dotn, j0,j1,j2,j3,j4,j5 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot18a = A( 1)*Z(j0) + A( 2)*Z(1+j0) + A( 3)*Z(2+j0)                &
             + A( 4)*Z(j1) + A( 5)*Z(1+j1) + A( 6)*Z(2+j1)                &
             + A( 7)*Z(j2) + A( 8)*Z(1+j2) + A( 9)*Z(2+j2)                &
             + A(10)*Z(j3) + A(11)*Z(1+j3) + A(12)*Z(2+j3)                &
             + A(13)*Z(j4) + A(14)*Z(1+j4) + A(15)*Z(2+j4)                &
             + A(16)*Z(j5) + A(17)*Z(1+j5) + A(18)*Z(2+j5)                &
             + dotn( C, U )

      return
      end
!-------------------------------------- 18 element side (left and right) -
      real function dot18s( Z, A, C, U, dotn,j0,j1,j2,j3,j4,j5,j6,j7,j8)
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot18s = A( 1)*Z(j0) + A( 2)*Z(1+j0)                                &
             + A( 3)*Z(j1) + A( 4)*Z(1+j1)                                &
             + A( 5)*Z(j2) + A( 6)*Z(1+j2)                                &
             + A( 7)*Z(j3) + A( 8)*Z(1+j3)                                &
             + A( 9)*Z(j4) + A(10)*Z(1+j4)                                &
             + A(11)*Z(j5) + A(12)*Z(1+j5)                                &
             + A(13)*Z(j6) + A(14)*Z(1+j6)                                &
             + A(15)*Z(j7) + A(16)*Z(1+j7)                                &
             + A(17)*Z(j8) + A(18)*Z(1+j8)                                &
             + dotn( C, U )

      return
      end
!-------------------------------------- 27 element interior --------------
      real function dot27i( Z, A, C, U, dotn,j0,j1,j2,j3,j4,j5,j6,j7,j8)
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot27i = A( 1)*Z(j0) + A( 2)*Z(1+j0) + A( 3)*Z(2+j0)                &
             + A( 4)*Z(j1) + A( 5)*Z(1+j1) + A( 6)*Z(2+j1)                &
             + A( 7)*Z(j2) + A( 8)*Z(1+j2) + A( 9)*Z(2+j2)                &
             + A(10)*Z(j3) + A(11)*Z(1+j3) + A(12)*Z(2+j3)                &
             + A(13)*Z(j4) + A(14)*Z(1+j4) + A(15)*Z(2+j4)                &
             + A(16)*Z(j5) + A(17)*Z(1+j5) + A(18)*Z(2+j5)                &
             + A(19)*Z(j6) + A(20)*Z(1+j6) + A(21)*Z(2+j6)                &
             + A(22)*Z(j7) + A(23)*Z(1+j7) + A(24)*Z(2+j7)                &
             + A(25)*Z(j8) + A(26)*Z(1+j8) + A(27)*Z(2+j8)                &
             + dotn( C, U )

      return
      end
!===========================================================================
!  Here are some general purpose dot product functions (fixed length)
!---------------------------------------------------------------------------
!  DOT1: vectors X and Y are of length 1
!---------------------------------------------------------------------------
      real function dot1( X, Y )
      dimension X(*), Y(*)

      dot1 = X(1)*Y(1)

      return
      end
!---------------------------------------------------------------------------
!  DOT2: vectors X and Y are of length 2
!---------------------------------------------------------------------------
      real function dot2( X, Y )
      dimension X(*), Y(*)

      dot2 = X(1)*Y(1) + X(2)*Y(2)

      return
      end
!---------------------------------------------------------------------------
!  DOT3: vectors X and Y are of length 3
!---------------------------------------------------------------------------
      real function dot3( X, Y )
      dimension X(*), Y(*)

      dot3 = X(1)*Y(1) + X(2)*Y(2) + X(3)*Y(3)

      return
      end
!---------------------------------------------------------------------------
!  DOT4: vectors X and Y are of length 4
!---------------------------------------------------------------------------
      real function dot4( X, Y )
      dimension X(*), Y(*)

      dot4 = X(1)*Y(1) + X(2)*Y(2) + X(3)*Y(3) + X(4)*Y(4)

      return
      end
!---------------------------------------------------------------------------
!  DOT5: vectors X and Y are of length 5
!---------------------------------------------------------------------------
      real function dot5( X, Y )
      dimension X(*), Y(*)

      dot5 = X(1)*Y(1) + X(2)*Y(2) + X(3)*Y(3) + X(4)*Y(4) + X(5)*Y(5)

      return
      end
!---------------------------------------------------------------------------
!  DOT6: vectors X and Y are of length 6
!---------------------------------------------------------------------------
      real function dot6( X, Y )
      dimension X(*), Y(*)

      dot6 = X(1)*Y(1) + X(2)*Y(2) + X(3)*Y(3) + X(4)*Y(4) + X(5)*Y(5)    &
           + X(6)*Y(6)

      return
      end
!---------------------------------------------------------------------------
!  DOT7: vectors X and Y are of length 7
!---------------------------------------------------------------------------
      real function dot7( X, Y )
      dimension X(*), Y(*)

      dot7 = X(1)*Y(1) + X(2)*Y(2) + X(3)*Y(3) + X(4)*Y(4) + X(5)*Y(5)    &
           + X(6)*Y(6) + X(7)*Y(7)

      return
      end

