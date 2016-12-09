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


