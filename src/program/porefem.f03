PROGRAM porefem
!-------------------------------------------------------------------------
! Program porefem
!             
!       Three-dimensional strain of an elastic solid using
!       8-, 14- or 20-node brick hexahedra. Mesh numbered in x-y
!       planes then in the z-direction. No global stiffness matrix
!       assembly. Diagonally preconditioned conjugate gradient solver.
!-------------------------------------------------------------------------
 USE main
 USE geom
 USE gaf95
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),npri=1
 INTEGER::cg_iters,cg_limit,fixed_freedoms,i,iel,k,loaded_nodes,ndim=3,   &
   ndof,nels,neq,nip,nn,nprops=2,np_types,nod,nodof=3,nr,nst=6,nxe,nye,   &
   nze,nlen,nlen_out,ii,iii,s
 INTEGER::ip,iq,maxfld,j,ipore
 REAL(iwp)::alpha,beta,cg_tol,det,one=1.0_iwp,penalty=1.0e20_iwp,up,      &
   zero=0.0_iwp,xl,yl,zl,dx1,dy1,dz1,small=1.0e-6_iwp,dtim=1.0_iwp
 REAL(iwp)::ee,ei,vi,meanqe,sdqe,meanqv,meanqv1,sdqv,sdqv1,pernom,porosity,ll,ul, &
            three=3.0_iwp,two=2.0_iwp,meanpo,sdpo,start_time,end_time,    &
            pcg_start_time,pcg_end_time
 CHARACTER(LEN=15)   :: element='hexahedron',argv,trash
 CHARACTER(LEN=15)   :: filename_out, format_string
 LOGICAL::cg_converged,solid=.true.
 LOGICAL::debug2 = .false.
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),no(:),    &
   node(:),num(:),sense(:),imap(:)
 REAL(iwp),ALLOCATABLE::bee(:,:),coord(:,:),d(:),dee(:,:),der(:,:),       &
   deriv(:,:),diag_precon(:),eld(:),fun(:),gc(:),g_coord(:,:),jac(:,:),   &
   km(:,:),loads(:),p(:),points(:,:),prop(:,:),sigma(:),store(:),         &
   storkm(:,:,:),u(:),value(:),weights(:),x(:),xnew(:),x_coords(:),       &
   y_coords(:),z_coords(:),val(:,:),temp1(:),temp2(:)
 REAL(iwp),ALLOCATABLE::enom(:),qe(:),qv(:),perc(:),qv1(:),kmi(:,:)
!-----------------------input and initialisation--------------------------
 integer ieplt(3)
 real, allocatable:: cfld(:),phifld(:),psifld(:),gamfld(:),efld(:),vfld(:)
 real  c(7),phi(7),psi(7),gam(7),e(7),v(7),R(6,6)
 integer kseed,istat,nsim,ns,nseed,iterm,nrfx,nrfy,nrfz,nfld,jfld,ifld
 real dx,dy,dz,thx,thy,thz,ti,ts
 logical verbos,debug,lxfld,dmpfld
 character sub2*128, job*80, sub1*80,varfnc*6
 common/dbgrfl/ istat, debug
 common/lastym/ ti, ts
 istat = 13
 iterm = 14
 dmpfld=.false.
 debug=.false.
 verbos=.false.
 lxfld=.false.
 sub1='none'
 sub2='none'
 job='equivalent_3D'
 varfnc='dlavx3'
 nfld=1
 jfld=2
 ifld=19
 call flush(6)
 CALL getname(argv,nlen)
 OPEN(10,FILE=argv(1:nlen)//'.dat')
 OPEN(11,FILE=argv(1:nlen)//'.res')
 CALL cpu_time(start_time)
 READ(10,*)nod,nip
 READ(10,*)xl,yl,zl
 READ(10,*)nxe,nye,nze
 READ(10,*)ei,vi
 READ(10,*)cg_tol,cg_limit
 READ(10,*)trash
 read(10,*)c
 read(10,*)phi
 read(10,*)psi
 read(10,*)gam
 read(10,*)e
 read(10,*)v
 READ(10,*)trash
 READ(10,*)r
 READ(10,*)thx,thy,thz
 READ(10,*)kseed,nsim
 READ(10,*)porosity
 kseed = iseed(kseed)
 WRITE(11,'(A)')job
 WRITE(11,'(A,I2,A)')'Analysis uses', nod,' noded elements'
 WRITE(11,*)
 WRITE(11,'(A,I2,A)')'Overall element dimensions'
 WRITE(11,'(A,F5.2)')'xl =',xl
 WRITE(11,'(A,F5.2)')'yl =',yl
 WRITE(11,'(A,F5.2)')'zl =',zl
 WRITE(11,*)
 WRITE(11,'(A,I2,A)')'Number of elements used'
 WRITE(11,'(A,I3)')'x-direction =',nxe
 WRITE(11,'(A,I3)')'y-direction =',nye
 WRITE(11,'(A,I3)')'z-direction =',nze
 WRITE(11,*)
 WRITE(11,'(A,I2,A)')'Stiffness parameters of intact rock'
 WRITE(11,'(A,F8.2)')"Young's modulus =",ei
 WRITE(11,'(A,F5.2)')"Poisson's ratio =",vi
 WRITE(11,*)
! WRITE(11,'(A)')"c phi psi gam e v prop erties of embankment (mean, SD, dist, L, U, m, s)"
 WRITE(11,'(A)')"Normal distribution (voids)"
 write(11,*)c
! write(11,*)phi
! write(11,*)psi
! write(11,*)gam
! write(11,*)e
! write(11,*)v
! WRITE(11,*)
! WRITE(11,'(A)')"Correlation matrix"
! WRITE(11,'(6f8.2)')r
! WRITE(11,*)
 WRITE(11,'(A,I2,A)')'Spatial correlation length of voids'
 WRITE(11,'(A,F8.2)')"x-direction =",thx
 WRITE(11,'(A,F8.2)')"y-direction =",thy
 WRITE(11,'(A,F8.2)')"z-direction =",thz
 WRITE(11,*)
 WRITE(11,'(A,I2,A)')'Starting seed'
 WRITE(11,'(I5)')kseed
 WRITE(11,*)
 WRITE(11,'(A,I2,A)')'Number of simulations'
 WRITE(11,'(I5)')nsim
 WRITE(11,*)
 WRITE(11,'(A,I2,A)')'Target porosity'
 WRITE(11,'(F5.2)')porosity
 ll=dphinv(0.5_iwp-porosity*0.5_iwp)
 ul=-ll
 !  dphinv(x) gives the abscissa that has the area x to its left in a normal distibution
 WRITE(11,*)
 WRITE(11,'(A,I2,A)')'Limits of porosity zone in standard space'
 WRITE(11,'(A,F7.4)')"lower =",ll
 WRITE(11,'(A,F7.4)')"upper =",ul
 CALL chkne3(istat,iterm,verbos,nxe,nye,nze,nrfx,nrfy,nrfz)
 maxfld=9*nrfx*nrfy*nrfz/8
 WRITE(11,'(A,I2,A)')'nrfx nrfy nrfz'
 WRITE(11,*)nrfx,nrfy,nrfz
 WRITE(11,*)
 CALL mesh_size(element,nod,nels,nn,nxe,nye,nze)
 ndof=nod*nodof

! Allocation of large arrays that dominate memory usage

 ALLOCATE(nf(nodof,nn))
 ALLOCATE(g_coord(ndim,nn))
 ALLOCATE(g_num(nod,nels),g_g(ndof,nels))
 ALLOCATE(imap(nels),etype(nels))
 
 ! ALLOCATE(enom(nxe*nye*nze)) ! - LM - allocated but not used

! Allocation of small arrays
 
 ALLOCATE(weights(nip),num(nod))
 ALLOCATE(points(nip,ndim),dee(nst,nst),coord(nod,ndim))
 ALLOCATE(jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),fun(nod),gc(ndim))
 ALLOCATE(bee(nst,ndof),km(ndof,ndof),eld(ndof),sigma(nst))
 ALLOCATE(g(ndof),x_coords(nxe+1),y_coords(nye+1),z_coords(nze+1),kmi(ndof,ndof))
 ALLOCATE(storkm(ndof,ndof,2))

! ALLOCATE(g(ndof),y_coords(nye+1),z_coords(nze+1),storkm(ndof,ndof,nels))
! ALLOCATE(nf(nodof,nn),points(nip,ndim),dee(nst,nst),coord(nod,ndim),     &
!   jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),fun(nod),gc(ndim),        &
!   bee(nst,ndof),km(ndof,ndof),eld(ndof),sigma(nst),g_coord(ndim,nn),     &
!   g_num(nod,nels),weights(nip),num(nod),g_g(ndof,nels),x_coords(nxe+1),  &
!   g(ndof),y_coords(nye+1),z_coords(nze+1),storkm(ndof,ndof,nels),        &
!   cfld(maxfld),phifld(maxfld),psifld(maxfld),gamfld(maxfld),efld(maxfld),&
!   vfld(maxfld),enom(nxe*nye*nze),imap(nels),etype(nels))

 x_coords=zero
 y_coords=zero
 z_coords=zero
 dx1=xl/dble(nxe)
 do i=2,nxe+1
 x_coords(i)=x_coords(i-1)+dx1
 enddo
 dy1=yl/dble(nye)
 do i=2,nye+1
 y_coords(i)=y_coords(i-1)+dy1
 enddo
 dz1=yl/dble(nze)
 do i=2,nze+1
 z_coords(i)=z_coords(i-1)-dz1
 enddo
 dx=dx1*dble(nxe)/dble(nrfx)
 dy=dy1*dble(nye)/dble(nrfy)
 dz=dz1*dble(nze)/dble(nrfz)
! call exca(dx,dy,thx,thy,varfnc,vf)
!
! prop(1)=c, prop(2)=phi, prop(3)=psi, prop(4)=gamma, prop(5)=e, prop(6)=v
!
 call propst(c)
 call propst(phi)
 call propst(psi)
 call propst(gam)
 call propst(e)
 call propst(v)
 IF(SUM(ABS(r))>6)lxfld=.TRUE.
 nf=1
  !-----------------------loop the elements to find global arrays sizes-----
 DO iel=1,nels
   CALL hexahedron_xz(iel,x_coords,y_coords,z_coords,coord,num)
   g_coord(:,num)=TRANSPOSE(coord)
 END DO
do i=1,nn
  if(abs(zl+g_coord(3,i))<=small)nf(3,i)=0
  if(abs(g_coord(1,i))<=small)nf(1,i)=0
  if(abs(g_coord(2,i)-yl)<=small)nf(2,i)=0
  if(abs(g_coord(3,i))<=small)nf(3,i)=0
  if(abs(g_coord(2,i))<=small)nf(2,i)=0
  if(abs(g_coord(1,i)-xl)<=small)nf(1,i)=0
enddo
 nf(2,1)=1
 nf(3,1)=1
 nf(1,nxe+1)=1
CALL formnf(nf)
do i=1,nn
 if(abs(g_coord(2,i))<=small)nf(2,i)=nf(2,1)
 if(abs(g_coord(3,i))<=small)nf(3,i)=nf(3,1)
 if(abs(g_coord(1,i)-xl)<=small)nf(1,i)=nf(1,nxe+1)
enddo
 neq=MAXVAL(nf)
 loaded_nodes=(nxe+1)*(nye+1)
 ALLOCATE(no(loaded_nodes),val(loaded_nodes,ndim))
 ALLOCATE(loads(0:neq),qe(nsim),qv(nsim),qv1(nsim),perc(nsim))
 WRITE(11,'(A,I5,A)')" There are",neq," equations"

!----------element stiffness integration, storage and preconditioner------
elements_1: DO iel=1,nels
  CALL hexahedron_xz(iel,x_coords,y_coords,z_coords,coord,num)
  CALL num_to_g(num,nf,g)
  g_num(:,iel)=num
  g_coord(:,num)=TRANSPOSE(coord)
  g_g(:,iel)=g
 END DO elements_1
 loads=zero
!-----------------------element stiffness integration and assembly--------
 CALL sample(element,points,weights)

   CALL deemat(dee,ei,vi)
   num=g_num(:,1)
   coord=TRANSPOSE(g_coord(:,num))
   kmi=zero
   eld=zero
   gauss_pts_1: DO i=1,nip
     CALL shape_der(der,points,i)
     jac=MATMUL(der,coord)
     det=determinant(jac)
     CALL invert(jac)
     deriv=MATMUL(jac,der)
     CALL beemat(bee,deriv)
     kmi=kmi+MATMUL(matmul(transpose(bee),dee),bee)*det*weights(i)
   END DO gauss_pts_1

   storkm(:,:,2) = kmi          ! medium
   storkm(:,:,1) = kmi/100._iwp ! void
   
   DEALLOCATE(g_num,g_coord)
   
!---------------------------------------------------------------------------

 do ns =1, nsim

 ALLOCATE(cfld(maxfld))
!ALLOCATE(phifld(maxfld),psifld(maxfld),gamfld(maxfld))
!ALLOCATE(efld(maxfld),vfld(maxfld))

! write(*,*)ns
 if( ns/=1)nseed = iseed(kseed + ns - 1)
!CALL simpl3(istat,iterm,verbos,cfld,phifld,psifld,gamfld,efld,           &
!                   vfld,c,phi,psi,gam,e,v,R,lxfld,thx,thy,thz,nrfx,      &
!                      nrfy,nrfz,dx,dy,dz,dmpfld,nfld,jfld,ieplt,         &
!                      ifld,job,sub1,sub2,varfnc,kseed,debug)
                       
 CALL simpl3cfld(istat,iterm,verbos,cfld,c,phi,psi,gam,e,v,R,lxfld,thx,   &
                       thy,thz,nrfx,                                      &
                       nrfy,nrfz,dx,dy,dz,dmpfld,nfld,jfld,ieplt,         &
                       ifld,job,sub1,sub2,varfnc,kseed,debug)                       
 
!DEALLOCATE(phifld,psifld,gamfld,efld,vfld)
 
 ! phifld, psifld, gamfld, efld & vfld are not used in the main program
 ! perhaps these should be allocated and deallocated in subroutine simpl3
                       
      k = 0
      do iii=1, nze
        do s = 1, nxe
          do i = 1, nye
            k = k + 1
            imap(k) = s + (nrfy-i)*nrfx+(nrfz-iii)*nrfx*nrfy
          enddo
        enddo
      enddo

 ALLOCATE(diag_precon(0:neq))
 diag_precon=zero
 
 j=0
 elements_2: DO iel=1,nels
   g=g_g(:,iel)
   ii=imap(iel)
   if(cfld(ii)>=ll.AND.cfld(ii)<=ul)then
    etype(iel)=1    !  type 1   (voids)
    km=kmi/100.0_iwp
    j=j+1
   else
    etype(iel)=2    !  type 2   (rock)
    km=kmi
   endif
   DO k=1,ndof
     diag_precon(g(k))=diag_precon(g(k))+km(k,k)
   END DO
 END DO elements_2
 
 IF(debug2) CALL sortvoid(etype)
 
 DEALLOCATE(cfld) ! not used again in the program until the next iteration
 
 ALLOCATE(p(0:neq))
 ALLOCATE(x(0:neq))
 ALLOCATE(xnew(0:neq))
 ALLOCATE(u(0:neq))
 ALLOCATE(d(0:neq))
 ALLOCATE(temp1(ndof),temp2(ndof))
          

 !-----------------------invert the preconditioner and get starting loads--
 perc(ns)= dble(j)/dble(nxe*nye*nze)
 loads=zero
 loads(nf(3,1))=xl*yl
 loads(0)=0.0
 diag_precon(1:)=one/diag_precon(1:)
 diag_precon(0)=zero
 d=diag_precon*loads
 p=d
!-----------------------pcg equation solution-----------------------------
 x=zero
 cg_iters=0
 CALL cpu_time(pcg_start_time) 
 pcg: DO
   cg_iters=cg_iters+1
   WRITE(*,*)"cg_iters = ", cg_iters
   u=zero
   elements_3: DO iel=1,nels
     g    = g_g(:,iel)
     u(g) = u(g)+MATMUL(storkm(:,:,etype(iel)),p(g))
   END DO elements_3
   up    = DOT_PRODUCT(loads,d)
   alpha = up/DOT_PRODUCT(p,u)
   xnew  = x+p*alpha
   loads = loads-u*alpha
   d     = diag_precon*loads
   beta  = DOT_PRODUCT(loads,d)/up
   p     = d+p*beta
   CALL checon(xnew,x,cg_tol,cg_converged)
   IF(cg_converged.OR.cg_iters==cg_limit)EXIT
 END DO pcg
 CALL cpu_time(pcg_end_time) 
 xnew(0)=zero
 !write(11,*)xnew(nf(3,1)),loads(nf(3,1)),xl
 qe(ns)=(xl*yl)/(xl*abs(xnew(nf(3,1))))
 qv(ns)=abs(xnew(nf(2,1)))/abs(xnew(nf(3,1)))
 qv1(ns)=abs(xnew(nf(1,nxe+1)))/abs(xnew(nf(3,1)))
 write(*,*)ns,cg_iters
 write(*,*) "Time in PCG = ", pcg_end_time - pcg_start_time
 write(11,'(2I4)')ns,cg_iters


   DEALLOCATE(p,x,xnew,u,temp1,temp2,diag_precon,d)          

!  IF(ns==1) THEN

     IF(ns < 10) THEN
       format_string = "(2A,I1)"
     ELSE IF(ns < 100) THEN
       format_string = "(2A,I2)"
     ELSE IF(ns < 1000) THEN
       format_string = "(2A,I3)"
     ELSE IF(ns < 10000) THEN
       format_string = "(2A,I4)"
     ELSE
       WRITE(*,*) "Limit of 9999 realisations in format_string"
       STOP
     END IF
     
     WRITE(filename_out,format_string) trim(argv),"_r",ns
     nlen_out = len(trim(filename_out))
     
     
     ALLOCATE(g_coord(ndim,nn),g_num(nod,nels))
     DO iel=1,nels
       CALL hexahedron_xz(iel,x_coords,y_coords,z_coords,coord,num)
       ! CALL num_to_g(num,nf,g)
       g_num(:,iel)=num
       g_coord(:,num)=TRANSPOSE(coord)
     END DO
!    CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,                    &
!                   loads(1:),1,npri,dtim,solid)
     CALL mesh_ensi(filename_out,nlen_out,g_coord,g_num,element,etype,nf,     &
                    loads(1:),1,npri,dtim,solid)
     DEALLOCATE(g_coord,g_num)
!  END IF


  ENDDO
 !!End simulation
 meanqe=sum(qe)/nsim
 sdqe=sqrt(dot_product(qe-meanqe,qe-meanqe)/(nsim-1))
 meanqv=sum(qv)/nsim
 sdqv=sqrt(dot_product(qv-meanqv,qv-meanqv)/(nsim-1))
 meanqv1=sum(qv1)/nsim
 sdqv1=sqrt(dot_product(qv1-meanqv1,qv1-meanqv1)/(nsim-1))
 meanpo=sum(perc)/nsim
 sdpo=sqrt(dot_product(perc-meanpo,perc-meanpo)/(nsim-1))
 write(11,*)
 WRITE(11,'(A)')"Simulation   Young's Modulus  Poisson's ratio(y)   Poisson's ratio(x)  Porosity"
 write(11,'(I5,8X,F10.4,8X,F10.3,10X,F10.3,5X,f10.3)')(i,qe(i),qv(i),qv1(i),perc(i),i=1,nsim)
 WRITE(11,*)
 WRITE(11,'(A,8X,F10.4,8X,F10.3,10X,F10.3,5X,f10.3)')'Means',meanqe,meanqv,meanqv1,meanpo
 WRITE(11,'(A,8X,F10.4,8X,F10.3,10X,F10.3,5X,f10.3)')'SDs  ',sdqe,sdqv,sdqv1,sdpo
! write(11,*)
! WRITE(11,'(A,I5)')" Number of cg iterations to convergence was",cg_iters
! WRITE(11,'(/A)')"  Node   x-disp      y-disp      z-disp"
! DO k=1,nn
!   WRITE(11,'(I5,3f12.4)')k,xnew(nf(:,k))
! END DO
CALL cpu_time(end_time)
WRITE(11,'(/,A,f12.4)')"time taken=",end_time-start_time
WRITE(*,*)"Time taken = ", end_time-start_time
STOP

CONTAINS

SUBROUTINE sortvoid(etype)

INTEGER,INTENT(INOUT) :: etype(:)
INTEGER               :: iel, nels, npores, count
REAL(iwp)             :: start_sort, end_sort
LOGICAL               :: debug = .false.
LOGICAL               :: found = .false.

CALL cpu_time(start_sort)

nels = UBOUND(etype,1)

IF(debug) THEN
  DO iel = 1,nels
    WRITE(*,*) "Element ", iel, " = type ", etype(iel)
  END DO
END IF

ipore = 0
DO iel = 1,nels
  IF(etype(iel)==1) npores = npores + 1
END DO

DO j = 1, npores
  count = 0
  IF(found) THEN
    found = .false.
  END IF
  DO iel = 1, nels
    count = count + 1
    IF(etype(iel)==1) THEN 
      found = .true.
      IF(debug) WRITE(*,*) "Position = ",count
      CYCLE
    END IF
  END DO
END DO

IF(debug) THEN
  WRITE(*,*) "Number of pore elements  = ",npores
  WRITE(*,*) "Number of solid elements = ",nels-npores
END IF

CALL cpu_time(end_sort)
WRITE(*,*)"Time to sort voids = ", end_sort - start_sort

END SUBROUTINE sortvoid


END PROGRAM porefem
