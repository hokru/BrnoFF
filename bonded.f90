! bonded potentials


! reminder
! r' -->  dx*irij


subroutine cov_bond_harm(mol,e)
use moldata
implicit none
type(molecule) :: mol
integer :: i,j,io,a,b
real(8) :: rij,e,r0,rab,rij2,irij,kb
real(8) dx,dy,dz,tmp

print*,'Covalent bond energies',mol%nbonds
e=0d0
do i=1,mol%nbonds
!      dx=mol%xyz(1,i)-mol%xyz(1,j)
!      dy=mol%xyz(2,i)-mol%xyz(2,j)
!      dz=mol%xyz(3,i)-mol%xyz(3,j)
      a=mol%ibond(i,1)
      b=mol%ibond(i,2)
      dx=mol%xyz(1,a)-mol%xyz(1,b)
      dy=mol%xyz(2,a)-mol%xyz(2,b)
      dz=mol%xyz(3,a)-mol%xyz(3,b)
      rij2=(dx*dx+dy*dy+dz*dz)
      rij=sqrt(rij2)
      irij=1d0/sqrt(rij2)      ! 1/rij
      r0=mol%r0(i)
      kb=mol%rk(i)
!      r0=99
!      kb=1
      e=e+kb*(rij-r0)**2
      print*,'ebond',i,a,b,e

      ! gradient, updates mol%g
      tmp=2d0*kb*(rij-r0)*irij
      mol%g(1,a)=mol%g(1,a)+tmp*dx
      mol%g(2,a)=mol%g(2,a)+tmp*dy
      mol%g(3,a)=mol%g(3,a)+tmp*dz

      mol%g(1,b)=mol%g(1,b)-tmp*dx
      mol%g(2,b)=mol%g(2,b)-tmp*dy
      mol%g(3,b)=mol%g(3,b)-tmp*dz
enddo

end subroutine


subroutine cov_angle_harm(mol,e)
use moldata
implicit none
type(molecule) :: mol
integer :: i,j,k,io
integer :: bond(mol%nat,mol%nat)
real(8) :: rij,rij2,e
real(8) :: ra(3),rb(3),rc(3)

open(newunit=io,file='bondmat.bin')
  read(io,*) bond
close(io)


do i=1,mol%nat-1
 do j=i,mol%nat
  if(i==j) cycle 
  do k=i+1,mol%nat
   if(k==j.or.k==i) cycle
   if(bond(i,j)==1.and.bond(j,k)==1) then
    ra=mol%xyz(1,i)
    rb=mol%xyz(1,j)
    rc=mol%xyz(1,k)
   endif
  enddo
 enddo
enddo

end subroutine
