!  *********************************************************************
!  *                                                                   *
!  *                         subroutine genprp                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.3
!  Written by Gordon A. Fenton, DalTech, Jan 31, 2001
!  Latest Update: Feb 22, 2005
!
!  PURPOSE  generates final soil property fields
!
!  DESCRIPTION
!  This routine sets the final soil property fields. In the event that
!  the soil property field follows a normal, lognormal, or bounded
!  distribution, this involves a transformation from the underlying
!  normally distributed field. If the soil property is deterministic,
!  this routine sets the entire field to the provided `mean' value. If
!  the soil property is a function of phi, this routine sets the field
!  to the required function.
!
!  ARGUMENTS
!
!    pfld	real array of size at least (nxe x nye) which, on output,
!		will contain the final soil property field. (output)
!
!       p	real vector of length at least 7 containing the statistics
!		of the soil property field. In particular,
!		   p(1) = soil property mean
!		   p(2) = soil property standard deviation,
!		   p(3) = soil property distribution type;
!			  = 0.0 if property is deterministic (at mean value)
!			  = 1.0 if property is normally distributed
!			  = 2.0 if property is lognormally distributed (logn)
!			  = 3.0 if property is bounded
!			  = 4.0 if property is a function of phi
!		   p(4) = lower bound (bounded), or mean log-property (logn),
!			  or constant in functional relationship (f(phi))
!		   p(5) = upper bound (bounded), or sd of log-property (logn),
!			  or slope in functional relationship (f(phi))
!		   p(6) = m parameter (if bounded), or function type (f(phi))
!			  where in the last case,
!				p(6) = 0.0 for 1*phi
!				p(6) = 1.0 for sin(phi)
!				p(6) = 2.0 for tan(phi)
!		   p(7) = s parameter (if bounded)
!         If the soil property is bounded or functional, then p(1) and
!		p(2) are ignored and the parameters c(4) through c(7)
!		completely describe the distribution (or function). (input)
!
!  phifld	real array of size at least (nxe x nye) which is assumed to
!		contain the friction angle soil property field. Ensure that
!		this is set before calling this routine since it is used
!		here in the event that the soil property is a function
!		of phi. (input)
!
!    nrfx	integer giving the number of elements describing the random
!		fields pfld and phifld in the x-direction (horizontally).
!		(input)
!
!    nrfy	integer giving the number of elements describing the random
!		fields pfld and phifld in the y-direction (vertically).
!		(input)
!
!  ltanfi	logical flag which is true if phifld contains the tan(phi)
!		random field, rather than phi directly. (input)
!
!  REVISION HISTORY:
!  1.1	correctly set twopt5 = 2.5 (Sep 12/01)
!  1.2	allow for phifld to contain tan(phi) field (Aug 14/02)
!  1.3	reordered numbering of f(phi) function types (0=1,1=sin,2=tan)
!	(Feb 22/05)
!-------------------------------------------------------------------------
      subroutine genprp(pfld,p,phifld,nrfx,nrfy,ltanfi)
      real pfld(nrfx,*), p(*), phifld(nrfx,*)
      integer nrfx, nrfy
      logical ltanfi
      data half/0.5/, one/1.0/, onept5/1.5/, twopt5/2.5/, thrpt5/3.5/
      data twopi/6.2831853071795864769/		! 2*pi
      data deg2rad/0.0174532925199/		! pi/180

      if( p(3) .lt. half ) then			! deterministic
         do 10 j = 1, nrfy
         do 10 i = 1, nrfx
            pfld(i,j) = p(1)
  10     continue
      elseif( p(3) .lt. onept5 ) then		! normal
         do 20 j = 1, nrfy
         do 20 i = 1, nrfx
            pfld(i,j) = p(1) + p(2)*pfld(i,j)
  20     continue
      elseif( p(3) .lt. twopt5 ) then		! lognormal
         do 30 j = 1, nrfy
         do 30 i = 1, nrfx
            pfld(i,j) = exp( p(4) + p(5)*pfld(i,j) )
  30     continue
      elseif( p(3) .lt. thrpt5 ) then		! bounded
         do 40 j = 1, nrfy
         do 40 i = 1, nrfx
            pfld(i,j) = p(4) + half*(p(5)-p(4))                           &
                *( one + tanh((p(6)+p(7)*pfld(i,j))/twopi) )
  40     continue
      else					! f(phi)
         if( p(6) .lt. half ) then		! prop = p(4)+p(5)*phi
            do 50 j = 1, nrfy
            do 50 i = 1, nrfx
               if( ltanfi ) then
                  ang = atan(phifld(i,j))/deg2rad ! angle in degrees
               else
                  ang = phifld(i,j)
               endif
               pfld(i,j) = p(4) + p(5)*ang
  50        continue
         elseif( p(6) .lt. onept5 ) then	! prop = p(4)+p(5)*sin(phi)
            do 60 j = 1, nrfy
            do 60 i = 1, nrfx
               if( ltanfi ) then
                  ang = atan(phifld(i,j))	! angle in radians
               else
                  ang = deg2rad*phifld(i,j)
               endif
               pfld(i,j) = p(4) + p(5)*sin(ang)
  60        continue
         elseif( p(6) .lt. twopt5 ) then	! prop = p(4)+p(5)*tan(phi)
            do 70 j = 1, nrfy
            do 70 i = 1, nrfx
               if( ltanfi ) then
                  ang = phifld(i,j)
               else
                  ang = tan(deg2rad*phifld(i,j))
               endif
               pfld(i,j) = p(4) + p(5)*ang
  70        continue
         endif
      endif

      return
      end
