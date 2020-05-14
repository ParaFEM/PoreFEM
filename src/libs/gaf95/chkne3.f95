!  *********************************************************************
!  *                                                                   *
!  *                         subroutine chkne3                         *
!  *                                                                   *
!  *********************************************************************
!  Single Precision Version 1.11
!  Written by Gordon A. Fenton, TUNS, Jul 14, 2000
!  Latest Update: Oct 18, 2002
!
!  PURPOSE	checks nxe, nye, and nze to ensure that they are compatible
!		with random fields produced by LAS3G
!
!  DESCRIPTION
!  This routine checks nxe, nye, and nze to ensure that they can be written
!  in the form;
!		nxe = k1*(2**m)
!		nye = k2*(2**m)
!		nze = k3*(2**m)
!
!  If not, the random field size given by (nrfx,nrfy,nrfz) is adjusted
!  upwards until the integers nrfx, nrfy, and nrfz do satisfy the
!  above equations for some non-negative integer m and positive integers
!  (k1,k2,k3) with k1*k2*k3 <= MXK.
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
!		messages are allowed to be sent to the screen. (output)
!
!     nxe	integer giving the number of elements describing the pillar
!		in the x-direction (width). (input)
!
!     nye	integer giving the number of elements describing the pillar
!		in the y-direction (depth). (input)
!
!     nze	integer giving the number of elements describing the pillar
!		in the z-direction (height). (input)
!
!    nrfx	integer giving the number of cells in the x-direction in the
!		simulated random field. This is normally equal to nxe, but if
!		nxe is not an integer of the form nxe = k1*(2**m), where m is
!		a non-negative integer, and k1 is a positive integer with
!		k1*k2*k3 < MXK then nrfx is set to the next larger
!		possible integer that does satisfy this requirement of LAS3G.
!		(input)
!
!    nrfy	integer giving the number of cells in the y-direction in the
!		simulated random field. This is normally equal to nye, but if
!		nye is not an integer of the form nye = k2*(2**m), where m is
!		a non-negative integer, and k2 is a positive integer with
!		k1*k2*k3 < MXK then nrfy is set to the next larger
!		possible integer that does satisfy this requirement of LAS3G.
!		(input)
!
!    nrfz	integer giving the number of cells in the z-direction in the
!		simulated random field. This is normally equal to nze, but if
!		nze is not an integer of the form nze = k3*(2**m), where m is
!		a non-negative integer, and k3 is a positive integer with
!		k1*k2*k3 < MXK then nrfz is set to the next larger
!		possible integer that does satisfy this requirement of LAS3G.
!		(input)
!
!  PARAMETERS:
!  MXK	maximum value that k1*k2*k3 can be. This should be identical to that
!	prescribed in GAF77 library routine las3g.f.
!
!  REVISION HISTORY:
!  1.11	added data definition for zero (Oct 18/02)
!-------------------------------------------------------------------------
      subroutine chkne3(istat,iterm,verbos,nxe,nye,nze,nrfx,nrfy,nrfz)
      parameter (MXK = 512)
      logical verbos
      data zero/0.0/, half/0.5/, pt9/0.99999/, two/2.0/

   1  format(a)
!					compute minimum m
      aa = float(nxe*nye*nze)/float(MXK)
      am = alog(aa)/(two*alog(two))
      if( am .lt. zero ) am = zero
      m  = int(am + pt9)		! round m up
      tm = two**m
!					now compute minimum k's
      ak1 = float(nxe)/tm
      k1  = int(ak1 + pt9)
      ak2 = float(nye)/tm
      k2  = int(ak2 + pt9)
      ak3 = float(nze)/tm
      k3  = int(ak3 + pt9)
!					new nxe, nye, and nze
      nrfx = k1*int(tm + half)
      nrfy = k2*int(tm + half)
      nrfz = k3*int(tm + half)
!					warn if we're changing anything

      if( (nrfx .ne. nxe) .or.                                            &
         (nrfy .ne. nye) .or.                                             &
         (nrfz .ne. nze) ) then
         if( verbos ) then
            call print3(iterm,                                            &
             'Note: desired number of elements (%i, %i, %i)%n',           &
              nxe,nye,nze)
            write(iterm,1)'         incompatible with LAS3G simulator.'
            call print3(iterm,                                            &
             '         Adjusting random field size to (%i, %i, %i)%n',    &
                       nrfx,nrfy,nrfz)
         endif
         call print3(istat,                                               &
          'Note: desired number of elements (%i, %i, %i)%n',              &
           nxe,nye,nze)
         write(istat,1)'         incompatible with LAS3G simulator.'
         call print3(istat,                                               &
          '         Adjusting random field size to (%i, %i, %i)%n',       &
                    nrfx,nrfy,nrfz)
      endif

      return
      end
