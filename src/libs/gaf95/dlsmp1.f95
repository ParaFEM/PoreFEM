!  *******************************************************************
!  *                                                                 *
!  *                     Real*8 Function dlsmp1                      *
!  *                                                                 *
!  *******************************************************************
!  Double Precision Version 1.31
!  Written by Gordon A. Fenton, May 13, 1996
!  Latest Update: May 9, 2001
!
!  PURPOSE  returns the covariance between two points along a 1-D simple
!           variance function process 
!
!  This routine returns the covariance between two points, separated by
!  the distance T, along a process having the simple covariance function
!
!                  var*dthx^3
!	  B(T) = --------------
!                (dthx + |T|)^3
!
!  where var is the point variance, and dthx is the scale of fluctuation
!  of the process. Both parameters, var and dthx, are brought into this
!  function via the common block /dparam/.
!
!  If var < 0, then this function returns the variance of a local average
!  of the process, |var|*V(T), averaged over the distance T.
!  The random process considered has a particularly simple variance
!  function,
!
!                     dthx
!	    V(T) = ----------
!                  dthx + |T|
!
!  where `dthx' is the scale of fluctuation of the process. This function
!  returns var*V(T), where `var' is the point variance of the process.
!  For more details, see page 204 of E.H. Vanmarcke's ``Random Fields:
!  Analysis and Synthesis.''
!
!  REVISION HISTORY:
!  1.1	added lvarfn flag to end of /dparam/ common block. Set it to
!	true if the variance function is to be returned. Otherwise this
!	function returns the covariance. (May 1/00)
!  1.2	eliminated lvarfn -- now return covariances only if var < 0 (Mar 2/01)
!  1.3	reversed default - now return covariances if var > 0. (Apr 11/01)
!  1.31	revised above documentation to reflect revision 1.3 (May 9/01)
! =======================================================================
      real*8 function dlsmp1( T )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/
      abs(x) = dabs(x)

      if( var .lt. zero ) then			! return variance function
         dlsmp1 = -var*dthx/(dthx + abs(T))
      else					! return covariance
         at = dthx + abs(T)
         if( at .eq. zero ) then
            dlsmp1 = var
         else
            bt = dthx/at
            dlsmp1 = var*bt*bt*bt
         endif
      endif

      return
      end
