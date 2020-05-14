!  *********************************************************************
!  *                                                                   *
!  *                         subroutine chknxe                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.2
!  Written by Gordon A. Fenton, TUNS, Jul 16, 1999
!  Latest Update: Oct 17, 2002
!
!  PURPOSE	checks nxe and nye to ensure that they are compatible with
!		random fields produced by LAS2G
!
!  DESCRIPTION
!  This routine checks nxe and nye to ensure that they can be written in the
!  form;
!		nxe = k1*(2**m)
!		nye = k2*(2**m)
!
!  If not, nxe and nye are adjusted upwards until they do satisfy the above
!  equations for some non-negative integer m and positive integers (k1,k2)
!  with k1*k2 <= MXK.
!
!  ARGUMENTS
!
!   istat	unit number connected to a file to which error and warning
!		messages can be issued. (input)
!
!   iterm	unit number connected to the screen. If verbos is true, then
!		error and warning messages are also sent to the screen. (input)
!
!  verbos	logical flag which is true if error, warning, and progress
!		messages are allowed to be sent to the screen. (input)
!
!     nxe	integer giving the number of elements describing the soil
!		mass in the x-direction (horizontally). (input)
!
!     nye	integer giving the number of elements describing the soil
!		mass in the y-direction (vertically). (input)
!
!    nrfx	integer giving the number of cells in the x-direction in the
!		simulated random field. This is normally equal to nxe, but if
!		nxe is not an integer of the form nxe = k1*(2**m), where m is
!		a non-negative integer, and k1 is a positive integer with
!		k1*k2 < MXK (see below) then nrfx is set to the next larger
!		possible integer that does satisfy this requirement of LAS2G.
!		(output)
!
!    nrfy	integer giving the number of cells in the y-direction in the
!		simulated random field. This is normally equal to nye, but if
!		nye is not an integer of the form nye = k2*(2**m), where m is
!		a non-negative integer, and k2 is a positive integer with
!		k1*k2 < MXK (see below) then nrfy is set to the next larger
!		possible integer that does satisfy this requirement of LAS2G.
!		(output)
!
!  PARAMETERS:
!  MXK	maximum value that k1*k2 can be. This should be identical to that
!	prescribed in GAF77 library routine last2g.f.
!
!  REVISION HISTORY:
!  1.01	corrected truncated documentation above (Mar 1/02)
!  1.1	separated RF size from FEM size (Jul 23/02)
!  1.2	added data statement for zero (Oct 17/02)
!-----------------------------------------------------------------------
      subroutine chknxe(istat,iterm,verbos,nxe,nye,nrfx,nrfy)
      parameter (MXK = 256)
      logical verbos
      data zero/0.0/, half/0.5/, pt9/0.99999/, two/2.0/

   1  format(a)
!					compute minimum m
      aa = float(nxe*nye)/float(MXK)
      am = alog(aa)/(two*alog(two))
      if( am .lt. zero ) am = zero
      m  = int(am + pt9)		! round m up
      tm = two**m
!					now compute minimum k's
      ak1 = float(nxe)/tm
      k1  = int(ak1 + pt9)
      ak2 = float(nye)/tm
      k2  = int(ak2 + pt9)
!					new random field dimensions
      nrfx = k1*int(tm + half)
      nrfy = k2*int(tm + half)
!					note if we're changing anything

      if( (nrfx .ne. nxe) .or. (nrfy .ne. nye) ) then
         if( verbos ) then
            call print2(iterm,                                            &
              'Note: provided number of elements (%i, %i)%n',nxe,nye)
            write(iterm,1)'         incompatible with LAS2G simulator.'
            call print2(iterm,                                            &
              '      Adjusting random field size to (%i, %i)%n',          &
                        nrfx,nrfy)
         endif
         call print2(istat,                                               &
           'Note: provided number of elements (%i, %i)%n',nxe,nye)
         write(istat,1)'         incompatible with LAS2G simulator.'
         call print2(istat,                                               &
           '      Adjusting random field size to (%i, %i)%n',             &
                     nrfx,nrfy)
      endif

      return
      end
