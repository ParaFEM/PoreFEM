!  ********************************************************************
!  *                                                                  *
!  *                          subroutine dismsh                       *
!  *                                                                  *
!  ********************************************************************
!  Single Precision Version 1.7
!  Written by Gordon A. Fenton, Geoffrey Paice, and D. Vaughan Griffiths
!  Latest Update: Apr 8, 2005
!
!  PURPOSE  displays the displaced finite element mesh with random property
!           field shown (optionally) as grey-scale
!
!  DESCRIPTION
!  This routine produces a PostScript plot of the displaced finite element
!  mesh. A grey-scale representation of the material property random field
!  is, optionally, overlain onto the mesh.
!
!  ARGUMENTS
!
!   loads	real vector of length at least equal to the number of free
!		degrees-of-freedom in the finite element mesh which contains
!		the computed nodal displacements. (input)
!
!      nf	integer array of size at least the number of nodes (nn) by 2
!		which contains one of
!		  i) a zero if the dof is constrained (ie. prescribed)
!		 ii) the index into `loads' for a (node,dir) if the dof
!		     is free (where node = 1, 2, ... nn, and dir = 1 for
!		     x and 2 for y directions).
!		(input)
!
!     inf	integer containing the leading dimension of the array nf.
!		(input)
!
! g_coord	real array of size at least 2 x nn containing the undisplaced
!		coordinates (x,y = g_coord(1,i),g_coord(2,i)) of each node.
!		(input)
!
!   g_num	integer array of size at least 8 x nels containing the global
!		dof numbers for each (node,dir) at the corners of each of the
!		nels elements. Nodes are numbered clockwise and each node has
!		2 dofs corresponding to the x and y directions. (input)
!
!	n	integer containing the number of unconstrained dof's (the
!		number of equations in the matrix problem). (input)
!
!      nn	integer containing the total number of nodes in the finite
!		element mesh. (input)
!
!     nod	integer containing the number of nodes per element. (input)
!
!    nels	integer containing the total number of elements in the problem
!		(input)
!
!    base	character string containing the basename of the output
!		PostScript file in its first ib characters. The output
!		file will have extension `.dis'. (input)
!
!      ib	integer such that base(1:ib) contains the basename (that is,
!		the name before the extension) of the output PostScript
!		file. (input)
!
!     ips	unit number to use to open the output PostScript file. This
!		file is both opened and closed in this routine. (input)
!
!    prop	real array of size at least nxe x nye which contains the
!		field of properties to (optionally) show as a grey-scale
!		on the plot. (input)
!
!    imap	integer vector of length at least equal to the number of
!		elements in the mesh containing the mapping between the
!		element numbers and the corresponding element of the
!		property array (the latter treated as a one-dimensional
!		vector stored by columns). (input)
!
!   lmesh	logical flag which is true if the finite element mesh is
!		to be shown on the output plot. (It doesn't make sense to
!		make both lmesh and lgrey false, so if lgrey is false, the
!		mesh is automatically shown.) (input)
!
!   lgrey	logical flag which is true if a grey-scale of the material
!		property field is to be shown on the mesh. (input)
!
!    llog	logical flag which is true if the natural logarithm of the
!		property field is to be shown in its grey-scale representation.
!		(Ignored if lgrey is false.) (input)
!
!   xsize	physical width of the soil mass on the printed page in inches.
!		(input)
!
!   ysize	physical height of the soil mass (undeformed) on the printed
!		page in inches. (input)
!
!    xoff	offset in inches from the right edge of the printed page to
!		the near edge of the soil mass. (input)
!
!    yoff	offset in inches from the bottom edge of the page to the
!		bottom edge of the soil mass. (input)
!
!   lfail	logical flag which is true if the slope is deemed to have
!		failed. (input)
!
!  REVISION HISTORY:
!  1.01	corrected truncated documentation above (Mar 1/02)
!  1.1	added plot size and offset parameters (Jul 22/02)
!  1.2	added Bounding box info to output PS file (Nov 9/02)
!  1.4	avoid division by 0 when computing dismag for undeformed mesh
!	(Oct 14/04)
!  1.5	added title to plot which indicates if the slope failed or not.
!	Accommodate case where greyscale prop is deterministic.  (Oct 17/04)
!  1.6	now map low to white, high to black, in greyscale. (Feb 17/05)
!  1.7	doubled displacement magnification (Apr 8/05)
!-----------------------------------------------------------------------
      subroutine dismsh2f(loads,nf,inf,g_coord,g_num,n,nn,nod,nels,       &
	                    base,ib,ips,prop,imap,lmesh,lgrey,llog,           &
                        xsize,ysize,xoff,yoff,lfail,propb,nrfx,ny1,ny2)
	  real prop(*),xload,propb(*),g_coord(2,*)
      real*8 loads(*)
	  integer g_num(8,*),nf(nn,*), imap(*)
	  integer nrfx,ny1,ny2
      character base*(*), ofile*256
      logical lmesh, lgrey, llog, lfail, ldet, ldetb
      data zero/0.0/, one/1.0/, scale/72.0/

   1  format(5a)
   4  format(f9.2,a)
   5  format(2f9.2,a)
   6  format(a,f5.3,a)
!                                      open output file
      if( ib+4 .gt. 256 ) then                                            
         write(6,1)                                                       &
           'Base filename too long to produce displaced mesh diagram.'
         return
      endif
      ofile = base(1:ib)
      ofile(ib+1:) = '.dis'
      open(ips, file = ofile, status = 'unknown' )
!					get range in greyscales (if lgrey)
      if( lgrey ) then
         if( llog ) then
			pmin = alog(prop(1))
            pmax = alog(prop(1))
			pminb = alog(propb(1))
            pmaxb = alog(propb(1))
			do 10 i = 1, nels
               ii = imap(i)
                if(ii>nrfx*ny2)then
                ii=ii-nrfx*ny2
                ap = alog(prop(ii))
                if( ap .lt. pmin ) pmin = ap
                if( ap .gt. pmax ) pmax = ap
				else
                ap = alog(propb(ii))
                if( ap .lt. pminb ) pminb = ap
                if( ap .gt. pmaxb ) pmaxb = ap
                endif
  10        continue
         else
			pmin = prop(1)
            pmax = prop(1)
			pminb =propb(1)
            pmaxb =propb(1)
            do 20 i = 1, nels
               ii = imap(i)
                if(ii>nrfx*ny2)then
                ii=ii-nrfx*ny2
                if( prop(ii) .lt. pmin ) pmin = prop(ii)
                if( prop(ii) .gt. pmax ) pmax = prop(ii)
				else
                if( propb(ii) .lt. pminb ) pminb = propb(ii)
                if( propb(ii) .gt. pmaxb ) pmaxb = propb(ii)
                endif
  20        continue
         endif
         if( pmax .eq. pmin ) then
            ldet = .true.
         else
            scgr = one/(pmax-pmin)
            ldet = .false.
         endif
         if( pmaxb .eq. pminb ) then
            ldetb = .true.
         else
            scgrb = one/(pmaxb-pminb)
            ldetb = .false.
         endif
      endif
!				find the range in the physical soil regime
      xmin = g_coord(1,1)		! g_coord must contain the (x,y)
      xmax = g_coord(1,1)		! coordinates of the displaced mesh
      ymin = g_coord(2,1)
      ymax = g_coord(2,1)
      do 30 i=2,nn
         if(g_coord(1,i) .lt. xmin) xmin = g_coord(1,i)
         if(g_coord(1,i) .gt. xmax) xmax = g_coord(1,i)
         if(g_coord(2,i) .lt. ymin) ymin = g_coord(2,i)
         if(g_coord(2,i) .gt. ymax) ymax = g_coord(2,i)
   30 continue
      width  = xmax - xmin
      height = ymax - ymin
!				find max displacement to physical size scale
      xload=zero
      do 40 i=1,n
         if( abs(loads(i)) .gt. xload ) xload = abs(loads(i))
   40 continue
      if( xload .eq. zero ) then
         dismag = one
      else
         dismag=0.4*height/xload
      endif
!				plot offset and scale factors
      xo = scale*xoff
      yo = scale*yoff
      sx = scale*xsize/width
      sy = scale*ysize/height
!				bounding box dimensions
      bxo = xo - 10.
      byo = yo - 10.
      bx1 = scale*(xoff + xsize) + 10.
      by1 = scale*(yoff + ysize) + 30.
!						start PostScript output
      write(ips,1)'%!PS-Adobe-1.0'
      write(ips,1)'%%DocumentFonts: Times-Roman'
      call print4(ips,'%%%%BoundingBox: %f %f %f %f%n',bxo,byo,bx1,by1)
      write(ips,1)'%%Pages: 1'
      write(ips,1)'%%EndComments'
      write(ips,1)'/m {moveto} def'
      write(ips,1)'/l {lineto} def'
      write(ips,1)'/s {stroke} def'
      write(ips,1)'/c {closepath} def'
      write(ips,1)'%%EndProlog'
      write(ips,1)'%%Page: 0 1'
      write(ips,1)'gsave'
!						move to plot origin
      write(ips,5) xo, yo, ' translate'
!						write out title
      write(ips,1)'12 /Times-Roman findfont exch scalefont setfont'
      write(ips,5) 10., scale*ysize + 15., ' m'
      if( lfail ) then
         write(ips,1)'(Slope failed.) show'
      else
         write(ips,1)'(Slope didn''t fail.) show'
      endif
!						draw the deformed mesh
      write(ips,4) 0.5, ' setlinewidth'
      do 60 i=1,nels
         ig = g_num(1,i)
         if(nf(ig,1).eq.0)then
            xdisp=0.
         else
            xdisp=loads(nf(ig,1))
         end if
         if(nf(ig,2).eq.0)then
            ydisp=0.
         else
            ydisp=loads(nf(ig,2))
         end if
         x=sx*(g_coord(1,ig)+dismag*xdisp- xmin)
         y=sy*(g_coord(2,ig)+dismag*ydisp- ymin)
         write(ips,5) x, y,' m'
         do 70 j=2,nod
            jj = g_num(j,i)
            if(nf(jj,1).eq.0)then
               xdisp=0.
            else
               xdisp=loads(nf(jj,1))
            end if
            if(nf(jj,2).eq.0)then
               ydisp=0.
            else
               ydisp=loads(nf(jj,2))
            end if
            x  = sx*(g_coord(1,jj)+ dismag*xdisp - xmin)
            y  = sy*(g_coord(2,jj)+ dismag*ydisp - ymin)
            write(ips,5) x, y,' l'
   70    continue
!					overlay greyscale of random property?
         if( lgrey ) then
            ii = imap(i)
                if(ii>nrfx*ny2)then
                ii=ii-nrfx*ny2
                 if( ldet ) then
                  gr = 0.7
                 else
                  if( llog ) then
                   gr = scgr*(pmax - alog(prop(ii)))
                  else
                   gr = scgr*(pmax - prop(ii))
                  endif
                 endif
                write(ips,6)'gsave c ',gr,' setgray fill grestore'
                else
                 if( ldetb ) then
                  gr = 0.7
                 else
                  if( llog ) then
                   gr = scgrb*(pmaxb - alog(propb(ii)))
                  else
                   gr = scgrb*(pmaxb - propb(ii))
                  endif
                 endif
!                write(ips,6)'gsave c ',gr,' setgray fill grestore'
 WRITE(ips, '(A,2F8.3,A)')'gsave c  ',gr,1.0-gr,' 0.0 setrgbcolor fill grestore'
				endif
         endif
         if( lmesh .or. .not. lgrey ) then
            write(ips,1)'c s'
         else
            write(ips,1)'newpath'
         endif
   60 continue
      write(ips,1)'grestore'
      write(ips,1)'showpage'
      close(ips)
!
      return
      end
