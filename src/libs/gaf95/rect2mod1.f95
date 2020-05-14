      subroutine rect2mod(nx1,nx2,nx3,ny1,ny2,ny3,rgrad,rgrad1,dx,g_coord,nf,inf,            &
                     g_num,maxel,nels,nn,n)
      real g_coord(2,inf)
      integer nf(inf,2),g_num(8,maxel),ny4
      real coord(8,2)
      integer num(8)
!
      ngrad=NINT(rgrad)
      ngrad1=NINT(rgrad1)
      IF(ngrad==0)ngrad=1
      IF(ngrad1==0)ngrad1=1
      dy=dx
      nye=ny1+ny2
      ny4=nye-ny3
      w1=dx*nx1
      w2=dx*nx2
      w3=dx*nx3
      h1=dy*ny1
      h2=dy*ny2
      h3=dy*ny3
      h4=dy*ny4
      
   IF(nx3>0)THEN
!
!      down the left hand side
!
     nm=0
     ic=0
     DO i=1,2*ny4
       nm=nm+1
       nf(nm,1)=0
       ic=ic+1
       nf(nm,2)=ic
     END DO
!
!      bottom left corner
!
     nm=nm+1
     nf(nm,1)=0
     nf(nm,2)=0
!
!      internal nodes in left rectangular block
!
     DO j=1,nx3-1
       DO i=1,ny4
         nm=nm+1
         ic=ic+1
         nf(nm,1)=ic
         ic=ic+1
         nf(nm,2)=ic
       END DO
       nm=nm+1
       nf(nm,1)=0
       nf(nm,2)=0
      
       DO i=1,2*ny4
         nm=nm+1
         ic=ic+1
         nf(nm,1)=ic
         ic=ic+1
         nf(nm,2)=ic
       END DO
       nm=nm+1
       nf(nm,1)=0
       nf(nm,2)=0
     END DO
      
     DO i=1,ny4
       nm=nm+1
       ic=ic+1
       nf(nm,1)=ic
       ic=ic+1
       nf(nm,2)=ic
     END DO
     nm=nm+1
     nf(nm,1)=0
     nf(nm,2)=0
   END IF
!
!      left sloping bit
!

   nxs=ngrad1*ny3

   IF(nx3==0)THEN
     nm=0
     ic=0
     DO i=1,2
       nm=nm+1
       ic=ic+1
       nf(nm,1)=ic
       ic=ic+1
       nf(nm,2)=ic
     END DO
     DO i=1,2*ny4
       nm=nm+1
       nf(nm,1)=0
       ic=ic+1
       nf(nm,2)=ic
     END DO 
   ELSE
     DO i=1,2*(1+ny4)
       nm=nm+1
       ic=ic+1
       nf(nm,1)=ic
       ic=ic+1
       nf(nm,2)=ic
     END DO
   END IF
  
   nm=nm+1
   nf(nm,1)=0
   nf(nm,2)=0  
   
   DO i=1,ny4+1
     nm=nm+1
     ic=ic+1
     nf(nm,1)=ic
     ic=ic+1
     nf(nm,2)=ic
   END DO

   nm=nm+1
   nf(nm,1)=0
   nf(nm,2)=0

   DO i=1,2*nye
     nm=nm+1
     ic=ic+1
     nf(nm,1)=ic
     ic=ic+1
     nf(nm,2)=ic
   END DO
   nm=nm+1
   nf(nm,1)=0
   nf(nm,2)=0
 
      DO j=2,nxs
   
       DO i=1,2*(ny4+1+(j-1)/ngrad1)
         nm=nm+1
         ic=ic+1
         nf(nm,1)=ic
         ic=ic+1
         nf(nm,2)=ic
       END DO
       nm=nm+1
       nf(nm,1)=0
       nf(nm,2)=0

       DO i=1,ny4+1+(j-1)/ngrad1
         nm=nm+1
         ic=ic+1
         nf(nm,1)=ic
         ic=ic+1
         nf(nm,2)=ic
       END DO
       nm=nm+1
       nf(nm,1)=0
       nf(nm,2)=0

       DO i=1,2*nye
         nm=nm+1
         ic=ic+1
         nf(nm,1)=ic
         ic=ic+1
         nf(nm,2)=ic
       END DO
       nm=nm+1
       nf(nm,1)=0
       nf(nm,2)=0
     END DO

!
!      internal nodes in middle rectangular block
!
	  if(nx3==0)then
      nm=0
      ic=0
      do 8111 i=1,2*(ny1+ny2)
      nm=nm+1
      nf(nm,1)=0
      ic=ic+1
 8111 nf(nm,2)=ic
!
!      bottom left corner
!
      nm=nm+1
      nf(nm,1)=0
      nf(nm,2)=0
	  endif
      do 82 j=1,nx1
      do 83 i=1,nye
      nm=nm+1
      ic=ic+1
      nf(nm,1)=ic
      ic=ic+1
   83 nf(nm,2)=ic
      nm=nm+1
      nf(nm,1)=0
      nf(nm,2)=0

      do 84 i=1,2*nye
      nm=nm+1
      ic=ic+1
      nf(nm,1)=ic
      ic=ic+1
   84 nf(nm,2)=ic
      nm=nm+1
      nf(nm,1)=0
      nf(nm,2)=0
   82 continue
!
      if(ngrad.gt.0)then
!
!      right sloping bit
!
      nxs=ngrad*ny1
      do 85 j=1,nxs
      do 86 i=1,nye-(j-1)/ngrad
      nm=nm+1
      ic=ic+1
      nf(nm,1)=ic
      ic=ic+1
 86   nf(nm,2)=ic
      nm=nm+1
      nf(nm,1)=0
      nf(nm,2)=0

      do 87 i=1,2*(nye-(j-1)/ngrad)
      nm=nm+1
      ic=ic+1
      nf(nm,1)=ic
      ic=ic+1
 87   nf(nm,2)=ic
      nm=nm+1
      nf(nm,1)=0
      nf(nm,2)=0
 85   continue
      end if
!
!      right hand side
!
      do 88 j=1,nx2
      do 90 i=1,ny2
      nm=nm+1
      ic=ic+1
      nf(nm,1)=ic
      ic=ic+1
 90   nf(nm,2)=ic
      nm=nm+1
      nf(nm,1)=0
      nf(nm,2)=0
!
      do 91 i=1,2*ny2
      nm=nm+1
      if(j.eq.nx2)then
      nf(nm,1)=0
      else
      ic=ic+1
      nf(nm,1)=ic
      end if
      ic=ic+1
 91   nf(nm,2)=ic
      nm=nm+1
      nf(nm,1)=0
      nf(nm,2)=0
 88   continue
      n=ic
      nn=nm
!
!      now for node numbering and coordinates
!
!      left hand bit
!
!
   If(nx3>=0)then
      yshift=dy*ny3
      nm=0
      do 18 ip=1,nx3
      do 18 iq=1,ny4
      nm=nm+1
      num(1)=(ip-1)*(3*ny4+2)+2*iq+1
      num(2)=num(1)-1
      num(3)=num(1)-2
      num(4)=(ip-1)*(3*ny4+2)+2*ny4+iq+1
      if(ip<nx3)then
      num(5)=ip*(3*ny4+2)+2*iq-1
      else
      num(5)=ip*(3*ny4+2)+2*iq+1
      endif
      num(6)=num(5)+1
      num(7)=num(5)+2
      num(8)=num(4)+1
      n1=num(8)
      coord(1,1)=(ip-1)*dx
      coord(2,1)=(ip-1)*dx
      coord(3,1)=(ip-1)*dx
      coord(5,1)=ip*dx
      coord(6,1)=ip*dx
      coord(7,1)=ip*dx
      coord(4,1)=(coord(3,1)+coord(5,1))*.5
      coord(8,1)=coord(4,1)
      coord(3,2)=-yshift-(iq-1)*dy
      coord(4,2)=-yshift-(iq-1)*dy
      coord(5,2)=-yshift-(iq-1)*dy
      coord(1,2)=-yshift-iq*dy
      coord(8,2)=-yshift-iq*dy
      coord(7,2)=-yshift-iq*dy
      coord(2,2)=(coord(1,2)+coord(3,2))*.5
      coord(6,2)=coord(2,2)
      do 18 i=1,8
      g_num(i,nm)=num(i)
      g_coord(1,num(i))=coord(i,1)
      g_coord(2,num(i))=coord(i,2)
18    continue
!
      if(ngrad1.gt.0)then
!
!      left sloping bit
!

      cx1=nx3*dx
      do 19 ii=1,ny3
      nys=ny4+ii
      nshift=n1
      xshift=cx1
      do 19 ip=1,ngrad1
      do 19 iq=1,nys
      nm=nm+1
      num(1)=nshift+(ip-1)*(3*nys+2)+2*iq+1
      num(2)=num(1)-1
      num(3)=num(1)-2
      num(4)=nshift+(ip-1)*(3*nys+2)+2*nys+iq+1
      if(ip==ngrad1.and.ii.lt.ny3)then
       num(5)=nshift+ip*(3*nys+2)+2*iq+1
       else
       num(5)=nshift+ip*(3*nys+2)+2*iq-1
      endif 
      num(6)=num(5)+1
      num(7)=num(5)+2
      num(8)=num(4)+1
      n1=num(8)
!
!      x-coordinates
!
      dxl=dx*rgrad1/dble(ngrad1)
      coord(1,1)=xshift+(ip-1)*dxl
      coord(3,1)=xshift+(ip-1)*dxl
      coord(5,1)=xshift+ip*dxl
      coord(7,1)=coord(5,1)
      if(iq.eq.1.and.ip.eq.1)coord(3,1)=coord(3,1)+0.5*dxl
      coord(2,1)=.5*(coord(1,1)+coord(3,1))
      coord(4,1)=.5*(coord(3,1)+coord(5,1))
      coord(6,1)=.5*(coord(5,1)+coord(7,1))
      coord(8,1)=.5*(coord(7,1)+coord(1,1))
      if(ip.eq.ngrad1)cx1=coord(7,1)
!
!      y-coordinates
!
      coord(1,2)=-(nye-nys)*dy-iq*dy
      coord(3,2)=-(nye-nys)*dy-(iq-1)*dy
      coord(5,2)=coord(3,2)
      coord(7,2)=coord(1,2)
      if(iq.eq.1.and.ip.gt.1)then
        coord(5,2)=coord(5,2)-float(ngrad1-ip)/float(ngrad1)*dy
        coord(3,2)=coord(3,2)-float(ngrad1-ip+1)/float(ngrad1)*dy
      end if
      if(iq.eq.1.and.ip.eq.1)then
        coord(5,2)=coord(7,2)+dy/float(ngrad1)
        coord(3,2)=coord(1,2)+dy/float(2*ngrad1)
      endif
      coord(2,2)=.5*(coord(1,2)+coord(3,2))
      coord(4,2)=.5*(coord(3,2)+coord(5,2))
      coord(6,2)=.5*(coord(5,2)+coord(7,2))
      coord(8,2)=.5*(coord(7,2)+coord(1,2))
      cx1=coord(7,1)
      do 19 i=1,8
      g_num(i,nm)=num(i)
      g_coord(1,num(i))=coord(i,1)
      g_coord(2,num(i))=coord(i,2)
19    continue
      end if
  endif
!
!      middle bit
!
!
	  if(nx3==0)then
       xshift=0
       nshift=0
	  else
       xshift=cx1
       nshift=n1
      endif
      do 8 ip=1,nx1
      do 8 iq=1,nye
      nm=nm+1
      num(1)=nshift+(ip-1)*(3*nye+2)+2*iq+1
      num(2)=num(1)-1
      num(3)=num(1)-2
      num(4)=nshift+(ip-1)*(3*nye+2)+2*nye+iq+1
      num(5)=nshift+ip*(3*nye+2)+2*iq-1
      num(6)=num(5)+1
      num(7)=num(5)+2
      num(8)=num(4)+1
      n1=num(8)
      if(ngrad.eq.0)n1=num(6)
      coord(1,1)=xshift+(ip-1)*dx
      coord(2,1)=xshift+(ip-1)*dx
      coord(3,1)=xshift+(ip-1)*dx
      coord(5,1)=xshift+ip*dx
      coord(6,1)=xshift+ip*dx
      coord(7,1)=xshift+ip*dx
      coord(4,1)=(coord(3,1)+coord(5,1))*.5
      coord(8,1)=coord(4,1)
      coord(3,2)=-(iq-1)*dy
      coord(4,2)=-(iq-1)*dy
      coord(5,2)=-(iq-1)*dy
      coord(1,2)=-iq*dy
      coord(8,2)=-iq*dy
      coord(7,2)=-iq*dy
      coord(2,2)=(coord(1,2)+coord(3,2))*.5
      coord(6,2)=coord(2,2)
      cx1=coord(7,1)
      do 8 i=1,8
      g_num(i,nm)=num(i)
      g_coord(1,num(i))=coord(i,1)
      g_coord(2,num(i))=coord(i,2)
 8    continue
!
      if(ngrad.gt.0)then
!
!      sloping bit
!
      dxr=dx*rgrad/dble(ngrad)
      do 9 ii=1,ny1
      nys=nye-(ii-1)
      nshift=n1
      xshift=cx1
      do 9 ip=1,ngrad
      do 9 iq=1,nys
      nm=nm+1
      num(1)=nshift+(ip-1)*(3*nys+2)+2*iq+1
      num(2)=num(1)-1
      num(3)=num(1)-2
      num(4)=nshift+(ip-1)*(3*nys+2)+2*nys+iq+1
      num(5)=nshift+ip*(3*nys+2)+2*iq-1
      num(6)=num(5)+1
      num(7)=num(5)+2
      num(8)=num(4)+1
      if(iq.eq.1)n1=num(6)
!
!      x-coordinates
!
      coord(1,1)=xshift+(ip-1)*dxr
      coord(3,1)=xshift+(ip-1)*dxr
      coord(5,1)=xshift+ip*dxr
      coord(7,1)=coord(5,1)
      if(iq.eq.1.and.ip.eq.ngrad)coord(5,1)=coord(5,1)-0.5*dxr
      coord(2,1)=.5*(coord(1,1)+coord(3,1))
      coord(4,1)=.5*(coord(3,1)+coord(5,1))
      coord(6,1)=.5*(coord(5,1)+coord(7,1))
      coord(8,1)=.5*(coord(7,1)+coord(1,1))
      if(ip.eq.ngrad)cx1=coord(7,1)
!
!      y-coordinates
!
      coord(1,2)=-(nye-nys)*dy-iq*dy
      coord(3,2)=-(nye-nys)*dy-(iq-1)*dy
      coord(5,2)=coord(3,2)
      coord(7,2)=coord(1,2)
      if(iq.eq.1.and.ip.lt.ngrad)then
        coord(3,2)=coord(3,2)-float(ip-1)/float(ngrad)*dy
        coord(5,2)=coord(5,2)-float(ip)/float(ngrad)*dy
      end if
      if(iq.eq.1.and.ip.eq.ngrad)then
        coord(3,2)=coord(3,2)-float(ip-1)/float(ngrad)*dy
        coord(5,2)=coord(1,2)+dy/float(2*ngrad)
      endif
      coord(2,2)=.5*(coord(1,2)+coord(3,2))
      coord(4,2)=.5*(coord(3,2)+coord(5,2))
      coord(6,2)=.5*(coord(5,2)+coord(7,2))
      coord(8,2)=.5*(coord(7,2)+coord(1,2))
      cx1=coord(7,1)
      do 9 i=1,8
      g_num(i,nm)=num(i)
      g_coord(1,num(i))=coord(i,1)
      g_coord(2,num(i))=coord(i,2)
 9    continue
      end if
!
!      right hand bit
!
!      if(ngrad.eq.0)n1=g_num(6,(nx1-1)*nye+ny1+)
!      xshift=(nx1+nx3+ny3*ngrad1+ny1*ngrad)*dx
      xshift=cx1
      yshift=dy*ny1
      nshift=n1
      do 10 ip=1,nx2
      do 10 iq=1,ny2
      nm=nm+1
        num(1)=nshift+(ip-1)*(3*ny2+2)+2*iq+1
        num(2)=num(1)-1
        num(3)=num(1)-2
        num(4)=nshift+(ip-1)*(3*ny2+2)+2*ny2+iq+1
        num(5)=nshift+ip*(3*ny2+2)+2*iq-1
        num(6)=num(5)+1
        num(7)=num(5)+2
        num(8)=num(4)+1
      coord(1,1)=xshift+(ip-1)*dx
      coord(2,1)=xshift+(ip-1)*dx
      coord(3,1)=xshift+(ip-1)*dx
      coord(5,1)=xshift+ip*dx
      coord(6,1)=xshift+ip*dx
      coord(7,1)=xshift+ip*dx
      coord(4,1)=(coord(3,1)+coord(5,1))*.5
      coord(8,1)=coord(4,1)
      coord(3,2)=-yshift-(iq-1)*dy
      coord(4,2)=-yshift-(iq-1)*dy
      coord(5,2)=-yshift-(iq-1)*dy
      coord(1,2)=-yshift-iq*dy
      coord(8,2)=-yshift-iq*dy
      coord(7,2)=-yshift-iq*dy
      coord(2,2)=(coord(1,2)+coord(3,2))*.5
      coord(6,2)=coord(2,2)
      do 10 i=1,8
      g_num(i,nm)=num(i)
      g_coord(1,num(i))=coord(i,1)
      g_coord(2,num(i))=coord(i,2)
 10   continue
      nels=nm
      return
      end
