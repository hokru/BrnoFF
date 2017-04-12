!subroutine get_rbc





subroutine rigid_coord(mol,rbc)
! return rigid-body coordinates for a molecule/fragment
use moldata
implicit none
type(molecule) mol
real(8) rotmat(3,3),tmp,axis(3,3),com(3)
real(8) phi,theta,psi,rbc(6,mol%nat)
integer i,j,k

! center of mass
call getcom(mol,com,.false.)
call getIntertia(mol,axis)
!

! axis direction ?
tmp=0d0
do j=1,2
 do i=1,mol%nat
  do k=1,3
   tmp=tmp+axis(k,j)*(mol%xyz(k,i)-com(k))
  enddo
  if(tmp<0d0) then
   axis(1:3,j)=-axis(1:3,j)
  endif
enddo
enddo

! principle axis are the row of rotation matrix
do i=1,3
 do j=1,3
  rotmat(i,j)=axis(j,i)
 enddo
enddo

! Euler angles from rotation matrix
! call roteuler(rotmat,phi,theta,psi)

! set rigid body coordinates


end subroutine

subroutine getCOM(mol,com,orient)
use moldata
use atomdata, only: ams
use constant, only: au2ang
implicit none
type(molecule) mol
logical orient
real(8) com(3), mmass
integer i,j,nat


!xyz=xyz*au2ang
mmass=0.0d0
com=0.0d0
do i=1,nat
 mmass=mmass+ams(mol%iat(i))
 do j=1,3
  com(j)=com(j)+mol%xyz(j,i)*ams(mol%iat(i))
 enddo
enddo
com=com/mmass

write(*,'(3x,a,3F12.5)') ' molecular mass      : ',mmass
write(*,'(3x,a,3F12.5)') ' center of mass (A)  : ',com(1:3)*au2ang

if(orient) then
! move molecule to COM
do i=1,nat
 do j=1,3
    mol%xyz(j,i)=mol%xyz(j,i)-com(j)
 enddo
enddo
endif
end

subroutine getIntertia(mol,mom)
! calculate moment of inertia and principle axis
use moldata
use atomdata, only: ams
implicit none
type(molecule) mol
integer i,nat
real(8) mom(3,3),e(3),rot(3)
real(8) x,y,z,m
real(8) conv,s,ddot

mom=0.0d0
!xyz=xyz*au2ang
do i=1,nat
x=mol%xyz(1,i)
y=mol%xyz(2,i)
z=mol%xyz(3,i)
m=ams(mol%iat(i))
mom(1,1)=mom(1,1)+(z**2+y**2)*m
mom(1,2)=mom(1,2)-x*y*m
mom(1,3)=mom(1,3)-x*z*m
mom(2,2)=mom(2,2)+(z**2+x**2)*m
mom(2,3)=mom(2,3)-y*z*m
mom(3,3)=mom(3,3)+(x**2+y**2)*m
enddo

!print*,xyz(1:3,1)
mom(2,1)=mom(1,2)
mom(3,1)=mom(1,3)
mom(3,2)=mom(2,3)


! diag, mom contains now eigenvectors
 call DiagSM(3,mom,e)

! handedness
s=mom(1,1)*(mom(2,2)*mom(3,3)-mom(3,2)*mom(2,3)) +  &
  mom(1,2)*(mom(2,3)*mom(3,1)-mom(2,1)*mom(3,3)) +  &
  mom(1,3)*(mom(2,1)*mom(3,2)-mom(2,2)*mom(3,1))
! invert if left-handed
if(s<0) then
  do i=1,3
  mom(i,1)=-mom(i,1)
  enddo
endif


end subroutine




real(8) function minmax(val)
implicit none
real(8) val
 minmax=min(1.0d0,max(-1.0d0,val))
end function


logical function do_switch(val)
implicit none
real(8) val,thr
 thr=1d-7
 do_switch=.false.
 if (abs(val) > thr) do_switch=.true.
end function

subroutine roteuler (rot,phi,theta,psi)
! Euler angles from rotation matrix
! following stackexchange and some lecture notes .. o_O
use constant, only: pi
implicit none
integer i
real(8) phi,theta,psi,thr,minmax
real(8) co_phi,co_theta,co_psi
real(8) s_phi,s_theta,s_psi
real(8) rot(3,3),diag(3)
logical do_switch


!threshold tolerance, what is a good value?
thr = 5d-8

!temp theta
theta = asin(minmax(-rot(1,3))) 
co_theta = cos(theta)
s_theta = -rot(1,3)

! theta : 90 or -90
if (abs(co_theta) <= thr) then
   phi = 0.0d0
   if (abs(rot(3,1)) <  thr) then
      psi = asin(minmax(-rot(2,1)/rot(1,3)))
   else if (abs(rot(2,1)) <  thr) then
      psi = acos(minmax(-rot(3,1)/rot(1,3)))
   else
      psi = atan(rot(2,1)/rot(3,1))
   end if

! theta /= +-90 
else
   if (abs(rot(1,1)) <  thr) then
      phi = asin(minmax(rot(1,2)/co_theta))
   else if (abs(rot(1,2)) <  thr) then
      phi = acos(minmax(rot(1,1)/co_theta))
   else
      phi = atan(rot(1,2)/rot(1,1))
   end if
   if (abs(rot(3,3)) <  thr) then
      psi = asin(minmax(rot(2,3)/co_theta))
   else if (abs(rot(2,3)) <  thr) then
      psi = acos(minmax(rot(3,3)/co_theta))
   else
      psi = atan(rot(2,3)/rot(3,3))
   end if
end if

! temp  phi and psi values
co_phi = cos(phi)
co_psi = cos(psi)
s_phi = sin(phi)
s_psi = sin(psi)

!  new diagonal of rotation matrix
diag(1) = co_theta * co_phi
diag(2) = s_psi*s_theta*s_phi + co_psi*co_phi
diag(3) = co_theta * co_psi

!  correct if needed 
if( do_switch(rot(1,1)-diag(1)).and.do_switch(rot(2,2)-diag(2))) phi=phi-sign(pi,phi)
if( do_switch(rot(1,1)-diag(1)).and.do_switch(rot(3,3)-diag(3))) theta=-theta+sign(pi,theta)
if( do_switch(rot(2,2)-diag(2)).and.do_switch(rot(3,3)-diag(3))) psi=psi-sign(pi,psi)

!     want positive values
if (phi   <= -pi) phi = pi
if (theta <= -pi) theta = pi
if (psi   <= -pi) psi = pi

end



subroutine rigid_to_xyz(mol,rbc)
! rigid-body coordinates rbc to cartesian xyz for a molecule object (updates xyz)
use moldata
implicit none
type(molecule) mol
integer i,j,l
real(8) rbc(6,mol%nat),com(3)
real(8) phi,theta,psi
real(8) coord(3,mol%nat)
real(8) rot(3,3)

! we need ref coords in the rb frame
!coord=?

do i=1,mol%nat
 com(1) = rbc(1,i)
 com(2) = rbc(2,i)
 com(3) = rbc(3,i)
 phi = rbc(4,i)
 theta = rbc(5,i)
 psi = rbc(6,i)

call get_rotmat(rot,theta,phi,psi)


! rotate+translate  atom
 mol%xyz(1:3,i)=matmul(rot,coord(:,i))+com 

enddo


end subroutine


subroutine get_rotmat(rot,theta,phi,psi)
! rot matrix from euler angles (xyz convention)
implicit none
real(8), intent(out) :: rot(3,3)
real(8), intent(in) :: theta, phi,psi
real(8) co_phi,co_theta,co_psi
real(8) s_phi,s_theta,s_psi

co_phi   = cos(phi)
s_phi    = sin(phi)
co_theta = cos(theta)
s_theta  = sin(theta)
co_psi   = cos(psi)
s_psi    = cos(psi)

rot=0d0
! xyz convention 
rot(1,1)=co_theta*co_phi
rot(1,2)=co_theta*s_phi
rot(1,3)=-s_theta
rot(2,1)=s_psi*s_theta*co_phi-co_psi*s_phi
rot(2,2)=s_psi*s_theta*s_psi+co_psi*co_phi
rot(2,3)=co_theta*s_psi
rot(3,1)=co_psi*s_theta*co_phi+s_psi*s_phi
rot(3,2)=co_psi*s_theta*s_phi-s_psi*co_phi
rot(3,3)=co_theta*co_psi


end subroutine



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

