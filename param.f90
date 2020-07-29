! MEMO
!  parm.dat
!  #FFid <string> : name of force field
!  #field : marks a new field
!  !      : comment/skipped
!  

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
integer :: n_c,n_vdw,n_bond,n_angle,n_torsion,n_chrg
real(8) :: par(10)
logical block(10)

block=.false.
n_c=0; n_vdw=0; n_angle=0; n_torsion=0;n_chrg=0

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
FF%npar=0

nbfix_ipair=0
nbfix_shift=0

io=22
open(io,file=pfile)
do
 read(io,'(a)',end=100) aa !each line read into aa string
 !-----------------------------------------------------------------
 if(fstr(aa,'#nbfix')) then
   ! <at_nr>  <at_nr> <shift>
   nbfix_npair=0
   i=1
   do
     read(io,'(a)',end=100) aa
     if(fstr(aa,'#')) exit
     if(fstr(aa,'!')) cycle
     backspace(io)
     read(io,*,end=100) nbfix_ipair(i,1), nbfix_ipair(i,2), nbfix_shift(i)
     i=i+1
   enddo
   nbfix_npair=i-1
   if(nbfix_npair>0) then 
    print*,'Found ',nbfix_npair, ' NBfix pairs!'
    do i=1,nbfix_npair
      print*,i, nbfix_ipair(i,:), nbfix_shift(i)
    enddo
   endif
 endif
 !-----------------------------------------------------------------
 if(fstr(aa,'#FFid'))  read(io,*) FF%id
 !-----------------------------------------------------------------
 if(fstr(aa,'#vdw')) then
!   if(n_vdw>0) exit
   i=0
   print*,'reading vdw parameters'
   do 
     read(io,'(a)',end=100) aa !each line read into aa string
       if(fstr(aa,'#')) then
          exit
          backspace(io)
       endif
     if(fstr(aa,'!')) cycle
     i=i+1
!     backspace(io)
     if(i>maxpar) stop 'maxpar reached!'
       call read_paramline(aa,FF%aname(i),par)
       FF%LJrad(i)=par(1)
       FF%LJe(i)=par(2)
       ! FF%npar=FF%npar+1
       ! write(*,*) FF%aname(i),FF%LJrad(i),FF%LJe(i)
       n_vdw=i
    ! read(io,*,end=100) FF%aname(i),FF%LJrad(i),FF%LJe(i)
   enddo
!   backspace(io)
  endif
  !-----------------------------------------------------------------
  if(fstr(aa,'#chrg')) then
!   if(n_chrg>0) exit
   i=0
   print*,'reading charges '
     do 
       read(io,'(a)',end=100) aa !each line read into aa string
       if(fstr(aa,'#')) then
          exit
          backspace(io)
       endif
       if(fstr(aa,'!')) cycle
       i=i+1
!       backspace(io)
       if(i>maxpar) stop 'maxpar reached!'
       call read_paramline(aa,FF%qname(i),par)
       FF%chrg(i)=par(1)
       n_chrg=i
!       read(io,*,end=100) FF%qname(i),FF%chrg(i)
     enddo
!     backspace(io)
  endif
  !-----------------------------------------------------------------
  if(fstr(aa,'#bond')) then
!   if(n_bond>0) exit
   i=0
   print*,'reading bond parameters '
    do 
       read(io,'(a)',end=100) aa !each line read into aa string
       if(fstr(aa,'#')) then
          exit
          backspace(io)
       endif
       if(fstr(aa,'!')) cycle
       i=i+1
       if(i>maxpar2) stop 'maxpar^2 reached!'
       call read_paramline(aa,FF%bond_name(i),par)
       FF%rk(i)=par(1)
       FF%r0(i)=par(2)
      FF%nbondpar = FF%nbondpar + 1
    enddo
    n_bond=i
  endif
  !
enddo

100 continue
close(io)

FF%npar=n_chrg+n_vdw
print*,''
print*,'FF name      : ',FF%id
print*,'non-bonded FF parameters:',FF%npar
print*,'vdW FF parameters:',n_vdw
print*,'charge FF parameters:',n_chrg
print*,'bond parameters: ',FF%nbondpar 
print*,''

end subroutine


subroutine read_paramline(str,namepar,par)
use string_parse
implicit none
character(*), intent(in) :: str
character(*), intent(out) :: namepar
real(8) :: par(10)
integer :: idx,len,n,i
character(255) :: rest

par=0

! find first occurance of 2 spaces
idx=index(str,'  ')

namepar=adjustl(trim(str(1:idx)))
call rmspace(namepar,namepar,len) 
! print* , '<',trim(namepar),'>'

rest=str(idx:)
call cstring(rest,n)
if(n>10) stop '>10 params in line'

do i=1,n
 call str_parse(rest,i,par(i))
enddo
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
   print*,'Assigning.. vdW ',trim(mol%aname(i)),i,esym(mol%iat(i))
   assigned=.true.
   mol%LJe(i)=FF%LJe(j)
   mol%LJrad(i)=FF%LJrad(j)
   exit
  endif
 enddo

 if(.not.assigned) then
  print*,'unknown atom:',i, esym(mol%iat(i)),mol%aname(i)
  stop 'Error in vdW assignment!'
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
  print*,'unknown atom:',i, esym(mol%iat(i)),mol%aname(i)
  stop 'Error in charge assignment!'
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
!        mol%nbonds = mol%nbonds + 1
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
