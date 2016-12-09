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



! print*,'',FF%npar

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
