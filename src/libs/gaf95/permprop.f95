      subroutine permprop (perm,nrfx,nxe,nye,ndim,dir,prop)
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      real perm(nrfx,*)
      character*(*) dir
      integer nxe,nye,iel,iq,ip,ndim,nels,nrfx
      REAL(iwp) prop(nxe*nye)
	  nels=nxe*nye
	  DO iel=1,nels
	  IF(dir=='x'.OR.dir=='r')THEN
      iq=(iel-1)/nxe+1
      ip=iel-(iq-1)*nxe
      ELSE
      ip=(iel-1)/nye+1
      iq=iel-(ip-1)*nye
      END IF
      prop(iel)=perm(ip,iq)
      enddo
      return
      end
