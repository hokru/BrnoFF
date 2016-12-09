


! amber vdw terms (using r0/2 !) and coulomb
subroutine nonbonded_amber(nfrag,fmol,evdw,ec)
use moldata
implicit none
integer i,j,k,l,nfrag
real(8) evdw,eij,rij,rij6,rab2,ec
real(8) r0ij,a,b,r0ij6
type(molecule) fmol(nfrag)

evdw=0d0
ec=0d0

! loop over pairs of fragments
do k=1,nfrag-1
 do l=k+1,nfrag

  ! loop over all intermolecular atom pairs
  do i=1,fmol(k)%nat
   do j=1,fmol(l)%nat
  !  print*, k,l,i,j
      r0ij=fmol(k)%LJrad(i)+fmol(l)%LJrad(j)
      eij=sqrt(fmol(k)%LJe(i)*fmol(l)%LJe(j))
      rij=rab2(fmol(k)%xyz(1,i),fmol(l)%xyz(1,j))  ! rij^2
      rij6=rij**3                  ! rij^6
      r0ij6=r0ij**6                  ! rij^6
      a=eij*(r0ij6*r0ij6)
      b=2d0*eij*r0ij6
      evdw=evdw + a/(rij6*rij6) - b/rij6
      ec=ec+(fmol(k)%chrg(i)*fmol(l)%chrg(j))/sqrt(rij)
    !  print*, fmol(l)%chrg(j),fmol(l)%chrg(j)
   enddo
  enddo
enddo
enddo


write(*,'(a,F12.4,a)') 'E(coulomb) :',ec,' [kcal/mol]'
write(*,'(a,F12.4,a)') 'E(LJ 6-12) :',evdw,' [kcal/mol]'

end subroutine


