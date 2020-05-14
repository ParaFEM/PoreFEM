!  ********************************************************************
!  *                                                                  *
!  *                        function uservf                           *
!  *                                                                  *
!  ********************************************************************
!  Double Precision Version 0.0
!  Written by ...
!
!  PURPOSE  returns the variance function corresponding to the user's model
!
!  This is currently a dummy function which can be modified by the user to
!  return his/her variance function. The name 'uservf' is already known to
!  RFLOW2D. Simply edit this module and then remake RFLOW2D. See 'vftrxx'
!  for an example function.
!
!  The parameters of the process are brought in through the common DPARAM
!  in which var is the point variance, PB is the ratio of the characteristic
!  area to the product of the directional scales of
!  fluctuation (see Vanmarcke, pg 240) (PB = PI/2 for exponential processes),
!  and thx, thy are the directional scales of fluctuation. Other parameters
!  are not needed and are included purely for consistency amongst my
!  covariance/variance functions. The function returns the variance function
!  value multiplied by the point variance var to obtain the variance of the
!  local average. Arguments to the function are as follows
!
!     x, y      dimensions of the rectangular averaging domain. (input)
!
!--------------------------------------------------------------------------
      real*8 function uservf( x, y )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, thx, thy, thz
!						replace below this line...
   1  format(a)
      write(6,1)'The variance function USERVF is not currently implemented'
      write(6,1)'Execution stopped in USERVF.'
      uservf = 0.d0
      stop
      end













