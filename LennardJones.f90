


! amber vdw terms (using r0/2 !) and coulomb
subroutine nonbonded_amber(evdw,ec)
use fragment, only : fmol, nfrag, flen, fiat
use FFparm, only: nbfix_npair,nbfix_ipair,nbfix_shift
implicit none
integer i,j,k,l !,nfrag
real(8) evdw,eij,rij,rij6,ec
real(8) r0ij,a,b,r0ij6
! type(molecule) fmol(nfrag)
real(8) dx,dy,dz,irij,er6,er12
integer ifix

integer, allocatable :: ifrag(:,:)
integer io
allocate(ifrag(sum(flen),9999))
open(newunit=io,file='bff_ifrag',form='unformatted')
read(io) ifrag
close(io)


evdw=0d0
ec=0d0
er6=0d0
er12=0d0

! loop over pairs of fragments
! !$omp parallel do default(none) shared(fmol,nfrag) private(k,l,i,j,r0ij,eij,dx,dy,dz,rij,rij6,irij,r0ij6,a,b) reduction(+:ec,er12,er6)
! turned off for now. Looping over fragments it not really worth it
!acc parallel loop reduction(+:ec,er12,er6) private(k,l,i,j,r0ij,eij,dx,dy,dz,rij,rij6,irij,r0ij6,a,b) copy(fmol,nfrag)
!$acc kernels  
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

      ! NFbix (as we dont pre-compute LJ params we need to do it here). 
      ! I REALLY would like to avoid ifrag matrix here
      if (any(nbfix_shift > 0.0d0)) then
        do ifix=1,nbfix_npair
          if ( (ifrag(i,k) == nbfix_ipair(ifix,1)) .and. (ifrag(j,l) == nbfix_ipair(ifix,2)) ) then
            r0ij=fmol(k)%LJrad(i)+fmol(l)%LJrad(j)-nbfix_shift(ifix)
              print*, nbfix_ipair(ifix,1),  nbfix_ipair(ifix,2)
            ! print*,'hit',r0ij
          endif
       enddo
     endif


      rij=(dx*dx+dy*dy+dz*dz) ! rij^2
      irij=1d0/sqrt(rij)
      rij6=rij**3         ! rij^6
      r0ij6=r0ij**6       ! rij^6
      a=eij*(r0ij6*r0ij6)
      b=2.0d0*eij*r0ij6
!      evdw=evdw + a/(rij6*rij6) - b/rij6
      er12=er12 +a/(rij6*rij6) 
      er6=er6 -b/rij6
      print*,'LJ',ifrag(i,k),ifrag(j,l),a/(rij6*rij6)-b/rij6,r0ij,fiat(i,k),fiat(j,l)
      ec=ec+(fmol(k)%chrg(i)*fmol(l)%chrg(j))*irij
   enddo
  enddo
enddo
enddo
!$acc end kernels  
!$omp end parallel do

evdw=evdw + er12 + er6

write(*,'(a,F12.4,a)') 'E(coulomb) :',ec,' [kcal/mol]'
write(*,'(a,F12.4,a)') 'E(LJ 6-12) :',evdw,' [kcal/mol]'
write(*,'(a,F12.4,a)') 'E(LJ 12) :',er12,' [kcal/mol]'
write(*,'(a,F12.4,a)') 'E(LJ 6) :',er6,' [kcal/mol]'

end subroutine


! amber vdw terms (using r0/2 !) and coulomb
subroutine nonbonded_amber_engrad(nfrag,fmol,evdw,ec)
use moldata
implicit none
integer i,j,k,l,nfrag
real(8) evdw,eij,rij2,rij6,ec
real(8) r0ij,a,b,r0ij6
type(molecule) fmol(nfrag)
real(8) dx,dy,dz,irij,irij2
real(8) rij14,rij8,gvdw(3),tmp,qq,tmp2

evdw=0d0
ec=0d0

! loop over pairs of fragments
! turned off for now. Looping over fragments it not really worth it
! !$omp parallel do default(none) shared(fmol,nfrag) private(k,l,i,j,r0ij,eij,dx,dy,dz,rij2,rij6,irij,irij2,r0ij6,a,b,gvdw,rij14,rij8,tmp,qq,tmp2) reduction(+:ec,evdw)
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
      ec=ec+qq*irij ! Coulomb
      ! print*,'LJ2', a/(rij6*rij6)-b/rij6
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

! coulomb energy with charge-penetration correction (e only)
! ALPHA
! Q. Wang et al JCTC 2015,11,2609; 10.1021/acs.jctc.5b00267
subroutine ScreenedCoulomb(nfrag,fmol)
use atomdata, only: rcov
use constant, only: au2ang,AmberElec
use moldata
implicit none
integer i,j,k,l,nat,ii,jj,nfrag
type(molecule) fmol(nfrag)
real(8) ec,dx,dy,dz,qq
real(8) irij,irij2,rij2
real(8) rab,gvdw(3),tmp2
real(8) thr,gnorm

real(8) alphai,alphaj,betai,betaj,r,zi,zj,val


ec=0d0

print*, 'screened   vs   conventional'
do k=1,nfrag-1
 do l=k+1,nfrag
  do i=1,fmol(k)%nat
   do j=1,fmol(l)%nat
      zi=val(fmol(k)%iat(i))
      zj=val(fmol(l)%iat(j))
      call setCP(fmol(k)%iff(i),alphai,betai)
      call setCP(fmol(l)%iff(j),alphaj,betaj)
      alphai=zi
      alphaj=zj
      zi=zi*AmberElec
      zj=zj*AmberElec
      dx=fmol(k)%xyz(1,i)-fmol(l)%xyz(1,j)
      dy=fmol(k)%xyz(2,i)-fmol(l)%xyz(2,j)
      dz=fmol(k)%xyz(3,i)-fmol(l)%xyz(3,j)
      r=sqrt(dx*dx+dy*dy+dz*dz) ! rij^2
      qq=zi*zj                                       &
             - zi*(zj-fmol(l)%chrg(j))*(1d0-exp(-alphaj*r))  &
             - zj*(zi-fmol(k)%chrg(i))*(1d0-exp(-alphai*r))  &
       + (zi-fmol(k)%chrg(i))*(zj-fmol(l)%chrg(j))*(1d0-exp(-betai*r))*(1d0-exp(-betaJ*r))
!      print*,zi*(zj-fmol(l)%chrg(j))*(1d0-exp(-alphaj*r)),zj*(zi-fmol(k)%chrg(i))*(1d0-exp(-alphai*r)),(zi-fmol(k)%chrg(i))*(zj-fmol(l)%chrg(j))*(1d0-exp(-betai*r))*(1d0-exp(-betaJ*r))
      print*, qq/r,(fmol(k)%chrg(i)*fmol(l)%chrg(j))/r
      ec=ec+qq/r
    enddo
  enddo
 enddo
enddo


write(*,'(a,F12.4,a)') '[test] E(screened coulomb) :',ec,' [kcal/mol]'

end subroutine


!!! below is quick & dirty and not fully correct !!!

real(8) function val(i)
implicit none
integer i
real(8) dat(5)
! valence electrons
!        H, C, N, O , P
data dat/2, 4, 5, 6, 5/

val=-1d0
if(i==1) val=dat(1)
if(i==6) val=dat(2)
if(i==7) val=dat(3)
if(i==8) val=dat(4)
if(i==15) val=dat(5)
if(val<0d0) stop 'error at the valence assignment!'
end function


subroutine setCP(atype,alpha,beta)
implicit none
integer atype
real(8) alpha,beta
real(8) adat(13),bdat(13)

data bdat/ &
!H(nonpol), H(arom), H(water),C(sp3), C(arom),C(sp2)
! 1         2           3         4       5     6
 1.999,    2.010,     2.004,    2.646, 2.708, 2.685 &
! N(sp3),N(arom),N(sp2)
!   7        8      9
 ,3.097   ,3.072  ,3.054 &
! O(sp3,oh,water),O(arom),O(sp2,carbonyl),P(phosphate)
! 10               11        12           13
,3.661,         4.282,      4.469,       2.360/

select case(atype)
 case(1) 
  beta=bdat(1)
 case(2)
  beta=bdat(4) ! Csp3
 case(3)
  beta=bdat(6) ! Csp2
 case(5)
  beta=bdat(12) ! O(arom)
 case(6,10,11)
  beta=bdat(11) ! O(arom)
 case(12)
  beta=bdat(10)
 case(7:8)
  beta=bdat(8) ! N(arom)
 case(13)
  beta=bdat(13) !phosphate
 case default
 stop 'unknown beta'
end select

end subroutine
