


! amber vdw terms (using r0/2 !) and coulomb
subroutine nonbonded_amber(nfrag,fmol,evdw,ec)
use moldata
implicit none
integer i,j,k,l,nfrag
real(8) evdw,eij,rij,rij6,rab2,ec
real(8) r0ij,a,b,r0ij6
type(molecule) fmol(nfrag)
real(8) dx,dy,dz,irij

evdw=0d0
ec=0d0

! loop over pairs of fragments
!$omp parallel do default(none) shared(fmol,nfrag) private(k,l,i,j,r0ij,eij,dx,dy,dz,rij,rij6,irij,r0ij6,a,b) reduction(+:ec,evdw)
do k=1,nfrag-1
 do l=k+1,nfrag

  ! loop over all intermolecular atom pairs
  do i=1,fmol(k)%nat
   do j=1,fmol(l)%nat
      r0ij=fmol(k)%LJrad(i)+fmol(l)%LJrad(j)
      eij=sqrt(fmol(k)%LJe(i)*fmol(l)%LJe(j))

      dx=fmol(k)%xyz(1,i)-fmol(l)%xyz(1,j)
      dy=fmol(k)%xyz(2,i)-fmol(l)%xyz(2,j)
      dz=fmol(k)%xyz(3,i)-fmol(l)%xyz(3,j)

      rij=(dx*dx+dy*dy+dz*dz) ! rij^2
      irij=1d0/sqrt(rij)
      rij6=rij**3                  ! rij^6
      r0ij6=r0ij**6                  ! rij^6
      a=eij*(r0ij6*r0ij6)
      b=2d0*eij*r0ij6
      evdw=evdw + a/(rij6*rij6) - b/rij6
      ec=ec+(fmol(k)%chrg(i)*fmol(l)%chrg(j))*irij
   enddo
  enddo
enddo
enddo
!$omp end parallel do


write(*,'(a,F12.4,a)') 'E(coulomb) :',ec,' [kcal/mol]'
write(*,'(a,F12.4,a)') 'E(LJ 6-12) :',evdw,' [kcal/mol]'

end subroutine


! amber vdw terms (using r0/2 !) and coulomb
subroutine nonbonded_amber_engrad(nfrag,fmol,evdw,ec)
use moldata
implicit none
integer i,j,k,l,nfrag
real(8) evdw,eij,rij2,rij6,rab2,ec
real(8) r0ij,a,b,r0ij6
type(molecule) fmol(nfrag)
real(8) dx,dy,dz,irij,irij2
real(8) rij14,rij8,gvdw(3),tmp,qq,tmp2

evdw=0d0
ec=0d0

! loop over pairs of fragments
!$omp parallel do default(none) shared(fmol,nfrag) private(k,l,i,j,r0ij,eij,dx,dy,dz,rij2,rij6,irij,irij2,r0ij6,a,b,gvdw,rij14,rij8,tmp,qq,tmp2) reduction(+:ec,evdw)
do k=1,nfrag-1
 do l=k+1,nfrag

  ! loop over all intermolecular atom pairs
  do i=1,fmol(k)%nat
   do j=1,fmol(l)%nat
      ! combination rules
      r0ij=fmol(k)%LJrad(i)+fmol(l)%LJrad(j)
      eij=sqrt(fmol(k)%LJe(i)*fmol(l)%LJe(j))
      ! intermediates
      dx=fmol(k)%xyz(1,i)-fmol(l)%xyz(1,j)
      dy=fmol(k)%xyz(2,i)-fmol(l)%xyz(2,j)
      dz=fmol(k)%xyz(3,i)-fmol(l)%xyz(3,j)
      rij2=(dx*dx+dy*dy+dz*dz) ! rij^2
      irij2=1d0/rij2           ! 1/rij^2
      irij=1d0/sqrt(rij2)      ! 1/rij
      rij6=rij2**3             ! rij^6
      r0ij6=r0ij**6            ! r0ij^6
      a=eij*(r0ij6*r0ij6)      ! A
      b=2d0*eij*r0ij6          ! B
      evdw=evdw + a/(rij6*rij6) - b/rij6 ! LJ
      qq=(fmol(k)%chrg(i)*fmol(l)%chrg(j))
      ! print*,fmol(k)%chrg(i),fmol(l)%chrg(j),fmol(k)%iat(i),fmol(k)%iat(i)
      ec=ec+qq*irij ! Coulomb

      ! cart. vdw gradient intermediate
      tmp=(-12d0*a*(irij2**7) + 6d0*b*(irij2**4))
      ! cart. coulomb gradient intermediate
      tmp2=-1d0*qq*irij2*irij
      ! local non-bonded pair-gradient
      gvdw(1)=(tmp + tmp2)*dx
      gvdw(2)=(tmp + tmp2)*dy
      gvdw(3)=(tmp + tmp2)*dz
      ! add to fragment arrays
      fmol(k)%g(1,i)=fmol(k)%g(1,i)+gvdw(1)
      fmol(k)%g(2,i)=fmol(k)%g(2,i)+gvdw(2)
      fmol(k)%g(3,i)=fmol(k)%g(3,i)+gvdw(3)
      fmol(l)%g(1,j)=fmol(l)%g(1,j)-gvdw(1)
      fmol(l)%g(2,j)=fmol(l)%g(2,j)-gvdw(2)
      fmol(l)%g(3,j)=fmol(l)%g(3,j)-gvdw(3)
   enddo
  enddo
enddo
enddo
!$omp end parallel do

write(*,'(a,F12.4,a)') 'E(coulomb) :',ec,' [kcal/mol]'
write(*,'(a,F12.4,a)') 'E(LJ 6-12) :',evdw,' [kcal/mol]'

end subroutine
!
! subroutine get_internal_numgrad(nat,iat,xyz,g)
! implicit none
! integer i,j,k,l,nat,iat(nat)
! real(8) xyz0(3,nat),xyz(3,nat)
! real(8) el,er,step,e,g(3,nat)
!
! step=0.005d0
! xyz0=xyz
! do i=1,nat
! write(*,'(a,I2,a,I2,a)') 'gradient of atom [',i, ']/[', nat,']'
!  do j=1,3
!    xyz(j,i)=xyz(j,i)+step
!    call nonbonded_amber(nfrag,fmol,evdw,ec)
!    er=energy
!    xyz(j,i)=xyz(j,i)-2d0*step
!    call nonbonded_amber(nfrag,fmol,evdw,ec)
!    el=energy
!    xyz(j,i)=xyz(j,i)+step
!    g(j,i)=(er-el)/(step*2d0)
!  enddo
! enddo
! end subroutine


! coulomb energy + cart. gradient for all atom pairs
subroutine eg_coulomb(nat,xyz,iat,atype,chrg,ec,g)
use atomdata, only: rcov
use constant, only: au2ang
implicit none
integer i,j,k,l,nat,iat,atype(nat),ii,jj
real(8) ec,dx,dy,dz,qq,chrg(nat)
real(8) irij,irij2,rij2
real(8) rab,gvdw(3),tmp2
real(8) g(3,nat),xyz(3,nat)
real(8) thr,gnorm

ec=0d0
g=0d0
thr=(15d0**2)

! atom-pair loop
do i=1,nat-1
 do j=i,nat
      if(i==j) cycle
      ii=atype(i)
      jj=atype(j)
      dx=xyz(1,i)-xyz(1,j)
      dy=xyz(2,i)-xyz(2,j)
      dz=xyz(3,i)-xyz(3,j)
      rij2=(dx*dx+dy*dy+dz*dz) ! rij^2
      ! if(rij2>thr) cycle
      irij2=1d0/rij2           ! 1/rij^2
      irij=1d0/sqrt(rij2)      ! 1/rij
      qq=chrg(i)*chrg(j)
      ! print'(I3,x,I3,2F10.4)',ii,jj,qq,irij
      ec=ec+qq*irij ! Coulomb

      ! cart. coulomb gradient intermediate
      tmp2=-1d0*qq*irij2*irij
      ! local non-bonded pair-gradient
      gvdw(1)=(tmp2)*dx
      gvdw(2)=(tmp2)*dy
      gvdw(3)=(tmp2)*dz
      ! add to gradient array
      g(1,i)=g(1,i)+gvdw(1)
      g(2,i)=g(2,i)+gvdw(2)
      g(3,i)=g(3,i)+gvdw(3)
      g(1,j)=g(1,j)-gvdw(1)
      g(2,j)=g(2,j)-gvdw(2)
      g(3,j)=g(3,j)-gvdw(3)

      ! hessian
 enddo
enddo

! call calc_gnorm(nat,g,gnorm)
! write(*,'(a,F12.4,a)') 'E(coulomb) :',ec,' [kcal/mol]'
! write(*,'(a,ES10.3)') 'Gnorm(coulomb) :',gnorm

end subroutine
