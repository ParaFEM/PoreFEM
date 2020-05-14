!  *********************************************************************
!  *                                                                   *
!  *                         subroutine fem2rf                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.01
!  Written by Gordon A. Fenton, TUNS, Nov 17, 1999
!  Latest Update: Mar 1, 2002
!
!  PURPOSE  maps the finite element mesh to the random field simulation
!           for RSLOPE2D
!
!  DESCRIPTION
!  This routine maps global element numbers from RSLOPE2D (see below) into
!  the two-dimensional random field described by the array rf(nxe x nye).
!  The axis orientations and element numbering are as follows;
!
!                                        y
!     --  +----------------              ^
!     ^   | 1 |             \            |
!     |   |---               \           |
!     |   | 2 |               \          |
!    ny1  |---                 \         +-----------> x
!     |   | 3 |                 \
!    _|_  |---   block     block \__________________
!     ^   | .      1         2                     |
!     |   | .                           block      |
!    ny2  | .                             3        |
!     |   |                                        |
!    _|_  +========================================+
!
!         |<---- nx1 ---->|      |<--- nx2 ------->|
!
!  The first index of rf proceeds in the positive x direction, while the
!  second proceeds in the positive y direction. The element numbers for the
!  finite element analysis of RSLOPE2D proceeds starting at the top left
!  corner of the mesh and counts top to bottom, then left to right. So the
!  last element is at the bottom right corner of the mesh.
!  Since the element number counting proceeds in the negative y direction
!  first, and then the positive x direction the correspondence between the
!  finite elements and the random field elements is as follows:
!
!	element 1 corresponds to rf(1,nye)   = rf(1 + (nye-1)*nxe)
!	element 2 corresponds to rf(1,nye-1) = rf(1 + (nye-2)*nxe)
!	   .
!	   .
!	   .
!	element (1+ny) corresponds to rf(2,nye) = rf(2 + (nye-1)*nxe)
!	   .
!	   .
!  where ny = ny1 + ny2.
!
!  ARGUMENTS
!
!   istat	unit number to which error and warning messages are to be
!		written. (input)
!
!   iterm	standard output unit number to which error and warning
!		messages are to be written if verbos is true. (input)
!
!  verbos	logical flag which is true if error and warning messages
!		are to be written to standard output. (input)
!
!     nx1	integer giving the number of elements in the x-direction
!		(horizontal) to the left of the embankment. (input)
!
!     nx2	integer giving the number of elements in the x-direction
!		(horizontal) to the right of the embankment. (input)
!
!     ny1	integer giving the number of elements in the y-direction
!		(vertically) in the embankment. (input)
!
!     ny2	integer giving the number of elements in the y-direction
!		(vertically) in the foundation, underlying the embankment.
!		(input)
!
!   ngrad	integer giving the slope of the embankment as the ratio of
!		the x-distance traversed for every unit fall in the
!		y-direction. (input)
!
!    nels	integer giving the number of elements in the problem. (input)
!
!    nrfx	integer giving the number of cells in the x-direction in the
!		simulated random field. This is normally equal to nxe, but if
!		nxe is not an integer of the form nxe = k1*(2**m), where m is
!		a non-negative integer, and k1 is a positive integer with
!		k1*k2 < MXK (see chknxe) then nrfx is the next larger
!		possible integer that does satisfy this requirement of LAS2G.
!		(input)
!
!    nrfy	integer giving the number of cells in the y-direction in the
!		simulated random field. This is normally equal to nye, but if
!		nye is not an integer of the form nye = k2*(2**m), where m is
!		a non-negative integer, and k2 is a positive integer with
!		k1*k2 < MXK (see chknxe) then nrfy is the next larger
!		possible integer that does satisfy this requirement of LAS2G.
!		(input)
!
!    imap	integer vector of length at least equal to the number of
!		elements in the problem which is used to store the FEM to
!		random field element mapping. This vector is formed upon
!		initialization of this routine and is thereafter used to
!		extract each element's appropriate material properties from
!		the random fields, cfld, phifld, etc. (output)
!
!
!  REVISION HISTORY:
!  1.01	corrected truncated documentation above, replaced \n with %n (Mar 1/02)
!-----------------------------------------------------------------------
      subroutine fem2rf(istat,iterm,verbos,nx1,nx2,ny1,ny2,ngrad,nels,    &
                       nrfx,nrfy,imap)
      integer imap(*)
      logical verbos

!					some constants we'll need later
      nxs = ny1*ngrad
      nx  = nx1 + nxs + nx2
      ny  = ny1 + ny2
!					block 1
      k = 0					! k counts the elements
      do 10 j = 1, nx1
      do 10 i = 1, ny
         k = k + 1
         imap(k) = j + (nrfy-i)*nrfx
  10  continue
!					block 2
      do 30 j = nx1 + 1, nx1 + nxs
         jj = j - nx1
         ii = (jj-1)/ngrad
         do 20 i = 1+ii, ny
            k = k + 1
            imap(k) = j + (nrfy-i)*nrfx
  20     continue
  30  continue
!					block 3
      do 40 j = nx1 + nxs + 1, nx
      do 40 i = ny1+1, ny
         k = k + 1
         imap(k) = j + (nrfy-i)*nrfx
  40  continue
!					check that k == nels
      if( k .ne. nels ) then
         call print1(istat,                                               &
           '*** ERROR: computed number of elements in fem2rf (%i)%n',     &
           k)
         call print1(istat,                                               &
        '           does not agree with that computed by rect (%i)%n',    &
           nels)
         if( verbos ) then
         call print1(iterm,                                               &
           '*** ERROR: computed number of elements in fem2rf (%i)%n',     &
           k)
         call print1(iterm,                                               &
        '           does not agree with that computed by rect (%i)%n',    &
           nels)
         endif
         close(istat)
         stop
      endif
!					all done, return
      return
      end
