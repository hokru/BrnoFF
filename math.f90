! distance between vectors a(3),b(3)
real(8) pure function rab(a,b)
implicit none
real(8), intent(in) :: a(3),b(3)
real(8) dx,dy,dz
 dx=a(1)-b(1)
 dy=a(2)-b(2)
 dz=a(3)-b(3)
 rab=sqrt(dx*dx+dy*dy+dz*dz)
end function

!real(8) pure function rab(a,b)
!implicit none
!real(8), intent(in) :: a(3),b(3)
!real(8) d(3)
! d=a-b
! rab= sqrt(dot_product(d,d))
!end function


! distance^2
real(8) pure function rab2(a,b)
implicit none
real(8), intent(in) :: a(3),b(3)
real(8) dx,dy,dz
 dx=a(1)-b(1)
 dy=a(2)-b(2)
 dz=a(3)-b(3)
 rab2=(dx*dx+dy*dy+dz*dz)
end function


subroutine DiagSM(xvar,mat,eig)
! lapack diag
implicit none
integer i,j,k
real(8), allocatable :: aux(:)
integer info,lwork,xvar
real(8) ,intent(inout) :: mat(xvar,xvar)
real(8) xx
real(8), intent(out) :: eig(xvar)

eig=0
call dsyev ('V','U',xvar,mat,xvar,eig,xx,-1,info)
lwork=int(xx)
allocate(aux(lwork))
call dsyev ('V','U',xvar,mat,xvar,eig,aux,lwork,info)
if(info/=0) print*,'Diagonalization failed !!'
end subroutine

