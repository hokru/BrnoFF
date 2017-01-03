! bonded potentials


! reminder
! r' -->  dx*irij

subroutine cov_bond_harm(mol,e)
use moldata
implicit none
type(molecule) :: mol
integer :: i,j,io
integer :: bond(mol%nat,mol%nat)
real(8) :: rij,e,r0,rab,rij2,irij
real(8) dx,dy,dz,tmp

open(newunit=io,file='bondmat.bin')
  read(io,*) bond
close(io)

e=0d0
do i=1,mol%nat-1
 do j=i,mol%nat
    if(i==j) cycle
      dx=mol%xyz(1,i)-mol%xyz(1,j)
      dy=mol%xyz(2,i)-mol%xyz(2,j)
      dz=mol%xyz(3,i)-mol%xyz(3,j)
      rij2=(dx*dx+dy*dy+dz*dz)
      rij=sqrt(rij2)
      irij=1d0/sqrt(rij2)      ! 1/rij
    ! r0=mol%r0(i,j)
      r0=1.0
      e=e+99d0*(rij-r0)**2

      ! gradient, updates mol%g
      tmp=2d0*99d0*(rij-r0)*irij
      mol%g(1,i)=mol%g(1,i)+tmp*dx
      mol%g(2,i)=mol%g(2,i)+tmp*dy
      mol%g(3,i)=mol%g(3,i)+tmp*dz

      mol%g(1,j)=mol%g(1,j)-tmp*dx
      mol%g(2,j)=mol%g(2,j)-tmp*dy
      mol%g(3,j)=mol%g(3,j)-tmp*dz
 enddo
enddo


end subroutine
