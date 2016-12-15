! intermolecular FF program
!
!
program iff
use moldata
use FFparm
use logic
use constant, only : au2ang
use fragment, only: fmol,nfrag
implicit none
integer nat,i,n
character(120) infile,outfile
real(8) ep,el,ec,evdw,rab,morse_OH
type(molecule) :: mol1, mol2
type(FFdata)   :: FF

real(8) :: start_time, end_time
real(8) :: time_nb,time_frag, time_ff


print*,''
print*,'|^^^^^^^^^^^^^^^^^^^^^^^|'
print*,'|                       |'
print*,'|   intermolecular FF   |'
print*,'|   (fragment based)    |'
print*,'|                       |'
print*,'| V: ALPHA              |'
print*,'| H.Kruse               |'
print*,'| mail2holger@gmail.com |'
print*,'|                       |'
print*,'|-----------------------|'
print*,''

print*,''
print*,'usage:'
print*,'  imff <structure>'
print*,''

call version

! defaults
echo=.true.
grad=.false.



call eval_options()


print*,''
!process molecules
print*,'structure: ',trim(filevec(1))
  call tmolrd(trim(filevec(1)),.true.,1,1,nat)
  mol1%nat=nat
  allocate(mol1%xyz(3,nat))
  allocate(mol1%iat(nat))
  allocate(mol1%iff(nat))
  allocate(mol1%mass(nat))
  allocate(mol1%chrg(nat))
  allocate(mol1%bond(nat,nat))
  allocate(mol1%cn(nat))
  allocate(mol1%aname(nat))
  allocate(mol1%LJe(nat),mol1%LJrad(nat))

  ! call tmolrd(trim(filevec(1)),.false.,mol1%xyz,mol1%iat,nat)
  call read_xyz(filevec(1),nat,mol1%iat,mol1%xyz,mol1%aname)
print*,''
! print*,'mol 2:'
!   call tmolrd(trim(filevec(2)),.true.,1,1,nat)
!   mol2%nat=nat
!   allocate(mol2%xyz(3,nat))
!   allocate(mol2%iat(nat))
!   allocate(mol2%iff(nat))
!   allocate(mol2%mass(nat))
!   allocate(mol2%chrg(nat))
!   allocate(mol2%bond(nat,nat))
!   allocate(mol2%cn(nat))
!   call tmolrd(trim(filevec(2)),.false.,mol2%xyz,mol2%iat,nat)
! print*,''

print*,'Reading parameter file..'
call read_parm('parm.dat',FF)
print*,''
! call assign_parm(FF,mol1)

print*,'Determining bonding info and fragments ...'
 call cpu_time(start_time)
call bondmat(mol1%nat,mol1%iat,mol1%xyz,mol1%aname,mol1%bond,mol1%cn)
 call cpu_time(end_time)
time_frag=end_time-start_time
print*,''

print*,'FF setup..'
 call cpu_time(start_time)
do i=1,nfrag
n=fmol(i)%nat
allocate(fmol(i)%LJe(n),fmol(i)%LJrad(n),fmol(i)%chrg(n))
print*,' Assigning FF parameters fragment ',i
call assign_parm(FF,fmol(i))
enddo
 call cpu_time(end_time)
time_FF=end_time-start_time
print*,''

print*,' running non-bonded part..'
call cpu_time(start_time)
call nonbonded_amber(nfrag,fmol,evdw,ec)

write(*,'(a)') ''
write(*,'(a)') ''
write(*,'(2x,a,F18.5,a)') 'E(tot)',evdw+ec,' [kcal/mol]'
write(*,'(a)') ''

print*,'done '
call cpu_time(end_time)
time_nb=end_time-start_time

print*,' timings: '
call prtime(6,time_frag,'fragment setup   ')
call prtime(6,time_FF,  'force field setup')
call prtime(6,time_nb,  'non-bonded part  ')
call prtime(6,time_nb+time_FF+time_frag,  'total            ')

end





subroutine prtime(io,t,string)
! prints elapsed time
implicit none
real(8) x,t,d,h,m,s
integer io
character(*) string

d=int(t/86400.0)!days
t=t-d*86400
h=int(t/3600.0) !hours
t=t-h*3600
m=int(t/60.0)     !mins
t=t-m*60
s=t          !secs

write(io,'(a,": ",x,i3," d",x,i3," h",x,i3," m",f5.1," s")') string,int(d),int(h),int(m),s

end subroutine prtime

