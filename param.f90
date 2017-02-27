! read params from file
subroutine read_parm(pfile,FF)
use FFparm
implicit none
character(*) pfile
character(120) aa
integer io,i
logical, external :: fstr
type(FFdata) FF
real(8), parameter :: amber2kcal=18.2223d0

FF%chrg=0d0
FF%LJe=0d0
FF%LJrad=0d0
FF%id='undefined'
FF%aname=''
FF%qname=''

io=22
open(io,file=pfile)
do
 read(io,'(a)',end=100) aa
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
       FF%chrg(i)=FF%chrg(i)*amber2kcal
     enddo
 endif
enddo

100 continue
close(io)

print*,''
print*,'FF name      : ',FF%id
print*,'FF parameters:',FF%npar
print*,''

end subroutine



subroutine assign_parm(FF,mol)
use FFparm
use moldata
implicit none
integer :: i,j
type(FFdata) FF
type(molecule) mol
logical assigned
character(2) esym

if(skip_charge) print*,'skipping charge assignment'
do i=1,mol%nat

assigned=.false.
 do j=1,FF%npar
   if( adjustl(trim(mol%aname(i))) == adjustl(trim(FF%aname(j))) ) then
   print*,'Assigning.. vdw ',trim(mol%aname(i)),i,esym(mol%iat(i))
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
    print*,'Assigning..chrg ',trim(mol%aname(i)),i,esym(mol%iat(i))
   assigned=.true.
   mol%chrg(i)=FF%chrg(j)
   exit
  endif
 enddo

 if(.not.assigned) then
  print*,'Error in charge assignment!'
  print*,'unknown atom:',i, esym(mol%iat(i)),mol%aname(i)
 endif

enddo




end subroutine


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
