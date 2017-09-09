! read params from file
subroutine read_parm(pfile,FF)
use FFparm
implicit none
character(*) pfile
character(120) aa
integer io,i,len
logical, external :: fstr
type(FFdata) FF
!integer, dimension(1:2) :: bondM_dim !my. dimensions of a bond matrix of a molecule
! real(8), parameter :: amber2kcal=18.2223d0

FF%chrg=0d0
FF%LJe=0d0
FF%LJrad=0d0
FF%id='undefined'
FF%aname=''
FF%qname=''
FF%nbondpar=0
FF%bond_name=''
FF%r0=0d0
FF%rk=0d0
!

io=22
open(io,file=pfile)
do
 read(io,'(a)',end=100) aa !each line read into aa string
 if(fstr(aa,'#FFid'))  read(io,*) FF%id
 if(fstr(aa,'#npar')) then
    read(io,*) FF%npar
    if(FF%npar>maxpar) stop 'increase maxpar !!'
    endif
 if(fstr(aa,'#vdw')) then
   print*,'reading vdw parameters'
   do i=1,FF%npar
     read(io,*,end=100) FF%aname(i),FF%LJrad(i),FF%LJe(i)
   enddo
  endif
  if(fstr(aa,'#chrg')) then
   print*,'reading charges '
     do i=1,FF%npar
       read(io,*,end=100) FF%qname(i),FF%chrg(i)
       FF%chrg(i)=FF%chrg(i)
     enddo
  endif
  if(fstr(aa,'#bond')) then
   print*,'reading bond parameters '
    do i=1,FF%npar**2
      read(io,*,end=100) FF%bond_name(i), FF%rk(i), FF%r0(i) !this supposes that there are NO spaces in definition of bond name, i.e C-C iv valid but C -C as in amber parm.dat will raise error!
      call rmspace(FF%bond_name(i),FF%bond_name(i),len) ! should fix spaces?
      FF%nbondpar = FF%nbondpar + 1
    enddo
  endif
  !
enddo

100 continue
close(io)

print*,''
print*,'FF name      : ',FF%id
print*,'non-bonded FF parameters:',FF%npar
print*,'bond parameters: ',FF%nbondpar ! my
print*,''

end subroutine


subroutine assign_parm(FF,mol)
use FFparm
use moldata
implicit none
integer :: i,j,k,l,numb 
type(FFdata) FF
type(molecule) mol
logical assigned
character(2) esym, a1, a2   !my: a1 and a2 - atom names, max 2 chars

if(skip_charge) print*,'skipping charge assignment'

numb=0  ! my  to keep track of number of bonds for which params are assigned

do i=1,mol%nat

assigned=.false.
 do j=1,FF%npar
   if( adjustl(trim(mol%aname(i))) == adjustl(trim(FF%aname(j))) ) then
!   print*,'Assigning.. vdw ',trim(mol%aname(i)),i,esym(mol%iat(i))
   assigned=.true.
   mol%LJe(i)=FF%LJe(j)
   mol%LJrad(i)=FF%LJrad(j)
   exit
  endif
 enddo

 if(.not.assigned) then
  print*,'Error in vdW assignment!'
  print*,'unknown atom:',i, esym(mol%iat(i)),mol%aname(i)
 endif

if(skip_charge) cycle
assigned=.false.
 do j=1,FF%npar
   if( trim(mol%aname(i)) == trim(FF%qname(j)) ) then
!    print*,'Assigning..chrg ',trim(mol%aname(i)),i,esym(mol%iat(i))
   assigned=.true.
   mol%chrg(i)=FF%chrg(j)
   exit
  endif
 enddo

 if(.not.assigned) then
  print*,'Error in charge assignment!'
  print*,'unknown atom:',i, esym(mol%iat(i)),mol%aname(i)
 endif

 ! bond assignement - r0 and rk:
 l=1
 do while (l < i) ! loop over all atom pairs in bond matrix but not diagonal and each pair only once
  assigned=.false.
  if (mol%bond(l,i)==1) then
    do k=1, FF%nbondpar
      call split_string(FF%bond_name(k),a1,a2,'-')  !for CX-HY from parm file makes: a1=CX, a2=HY
      if( (trim(mol%aname(l)) == trim(a1) .AND. trim(mol%aname(i)) == trim(a2)) .OR. (trim(mol%aname(i)) == trim(a1) .AND. trim(mol%aname(l)) == trim(a2)) ) then
        numb = numb + 1
        mol%nbonds = mol%nbonds + 1
        !print*,'nbonds' , mol%nbonds
        mol%r0(numb) = FF%r0(k) 
        mol%rk(numb) = FF%rk(k)
        assigned=.true.
        exit
      endif
    enddo
    if(.not.assigned) then
      print*,'Error in bond assignment!'
      print*,'unknown bond: ', trim(mol%aname(i)), '-' ,mol%aname(l)
    endif
  endif
  l=l+1
 enddo

enddo
end subroutine


subroutine print_info_FFnb(mol)
use moldata
type(molecule), intent(in):: mol
character(2) esym
real(8) s,d

print*, 'non-bonded parameter info'
write(*,'(2x,a6,1x,a2,1x,a4,1x,a10,1x,a10,1x,a10)') 'index','el','type','charge','LJe','LJrad'
       !xxFFFFFFxAAxAAAAxFFFFFFFFFFxFFFFFFFFFFxFFFFFFFFFF
do i=1,mol%nat
 write(*,'(2x,I6,1x,a2,1x,a4,1x,F10.4,1x,F10.4,1x,F10.4)') &
    i,esym(mol%iat(i)),mol%aname(i),mol%chrg(i),mol%LJe(i),mol%LJrad(i)
enddo

s=sum(mol%chrg)
write(*,'(2x,a,2x,F10.4)') 'molecular charge: ',s
! d=abs(int(s))-abs(s)
! if(abs(d)>epsilon(1d0)) then
! write(*,'(2x,a,ES10.3,2x,ES10.3)') 'deviation from integer value: ', d
! endif

end subroutine

! my
subroutine print_info_FFb(mol)
use moldata
type(molecule), intent(in):: mol
integer :: i,j,k
character(2) esym

print*, 'bond parameter info'
print*, '  bonds: ',mol%nbonds
write(*,'(2x,a6,1x,a5,1x,a5,1x,a10,1x,a10)') 'index','el','type','r0','rk'
k=1
do i=1,mol%nat
j=1
do while (j < i)
  if (mol%bond(j,i)==1) then
    write(*,'(2x,I6,1x,a2,a1,a2,1x,a2,a1,a2,1x,F10.4,1x,F10.4)') &
    k, trim(esym(mol%iat(j))), '-', trim(esym(mol%iat(i))), trim(mol%aname(j)), '-', trim(mol%aname(i)), mol%r0(k), mol%rk(k)
    k=k+1
  endif
  j=j+1
enddo
enddo

end subroutine
!


subroutine load_internalFF(FF)
use FFparm
implicit none
integer :: i,io,s,q
type(FFdata) FF

FF%id='model1'
FF%npar=25

print*,'type   charge'
open(newunit=io,file='params')
do i=1,FF%npar
 read(io,*) FF%itype(i), FF%chrg(i)
 print'(I3,x,F10.4)',FF%itype(i), FF%chrg(i)
enddo
close(io)

end

subroutine assign_internal(FF,mol)
use FFparm
use moldata
implicit none
integer i,j,s
real(8) q
logical assigned
type(FFdata) FF
type(molecule) mol

s=0
do i=1,mol%nat
assigned=.false.
do j=1,FF%npar
  if(FF%itype(j)==mol%iff(i)) then
    mol%chrg(i)=FF%chrg(j)
    assigned=.true.
    if(FF%itype(j)==1) s=s+1
  endif

enddo
if(.not.assigned) stop ' error in param assignment (assign_internal)'
enddo

print*,''
print'(a,I6)', ' # H atoms', s
print*,''
print*,' initial charge: ', sum(mol%chrg)
q=(0d0-sum(mol%chrg))/s
print*, 'H charge correction', q

print*,''
print*,'index  type  charge'
do i=1,mol%nat
  if(mol%iff(i)==1) then
    mol%chrg(i)=mol%chrg(i)+q
  endif
print'(I6,x,I3x,F8.4)',i,mol%iff(i),mol%chrg(i)
enddo

print*,''
print*,' final charge: ', sum(mol%chrg)
print*,''


end
