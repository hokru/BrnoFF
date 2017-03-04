! intermolecular FF program
!
!
program iff
use moldata
use FFparm
use logic
use constant, only : au2ang,amberelec
use fragment, only: fmol,nfrag
implicit none
integer nat,i,n,j,k,io,ixyz
character(120) infile,outfile
real(8) etotal,ec,evdw,rab,morse_OH,gnorm,ebond,s
type(molecule) :: mol1
type(FFdata)   :: FF

real(8) :: start_time, end_time
real(8) :: time_nb,time_frag, time_ff

! handle this better!??
integer, allocatable :: ifrag(:,:)

#ifdef OPENMP
integer omp_get_num_threads,nomp
#endif


call print_header()

! defaults
ixyz=1
echo=.true.
grad=.false.
nchess=.false.
s=0.0d0

#ifdef OPENMP
    print*,''
    !$omp parallel
    nomp=omp_get_num_threads()
    !$omp end parallel
    print'(a,I4)',' OpenMP threads: ',nomp
    print*,''
#endif

call eval_options()

if(nchess) print*,'NON-COVALENT MODEL HESSIAN'

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
  allocate(mol1%g(3,nat))

! read molecule data
if(nchess) then
    call tmolrd(trim(filevec(1)),.false.,mol1%xyz,mol1%iat,nat)
    print*,''
    print*,'loading internal parameters'
    call load_internalFF(FF)
    print*,''
    ! allocate(mol1%chrg(nat))
    print*,' Assigning FF parameters  '
    call get_atype(mol1%nat,mol1%iat,mol1%xyz,mol1%iff)
    call assign_internal(FF,mol1)
    print*,''
    call eg_coulomb(mol1%nat,mol1%xyz,mol1%iat,mol1%iff,mol1%chrg,ec,mol1%g)
    call calc_gnorm(mol1%nat,mol1%g,gnorm)
    write(*,'(a,F12.4,a)') 'E(coulomb) :',ec,' [kcal/mol]'
    write(*,'(a,ES10.3)') 'Gnorm(coulomb) :',gnorm
    call cov_bond_harm(mol1,ebond)
    write(*,'(a,F12.4,a)') 'E(bond) :',ebond,' [kcal/mol]'
    ! call numhess(mol1)
else ! FRAGMENT-BASED AMBER-LIKE FF
    call det_xyz_type(filevec(1),ixyz)
    select case (ixyz)
      case(1)
      case default
       call error(6,'no enhanced xyz file found')
      case(2)
      print*,'Found atom types in xyz file'
       call read_xyz_with_type(filevec(1),nat,mol1%iat,mol1%xyz,mol1%aname)
      case(3)
       skip_charge=.true.
       print*,'Found atom types and charges in xyz file'
       call read_xyz_with_type_charge(filevec(1),nat,mol1%iat,mol1%xyz,mol1%aname,mol1%chrg)
      end select
    print*,''

    print*,'Reading parameter file:  <parm.dat>'
    call read_parm('parm.dat',FF)
    print*,''

    print*,'Determining bonding info and fragments ...'
     call cpu_time(start_time)
    call bondmat(mol1%nat,mol1%iat,mol1%xyz,mol1%aname,mol1%bond,mol1%cn,mol1%chrg)
     call cpu_time(end_time)
    time_frag=end_time-start_time
    print*,''

    print*,'FF setup..'
     call cpu_time(start_time)
    do i=1,nfrag
    n=fmol(i)%nat
    print*,' Assigning FF parameters fragment ',i
    ! if we already have the charges we must set the fragment charges now
    call assign_parm(FF,fmol(i))
    call print_info_FFnb(fmol(i))
    print*,''
    ! need to scale the charge for units used in amber
    ! (units of the electron charge and kcal)
    s=s+sum(fmol(i)%chrg)
    fmol(i)%chrg=fmol(i)%chrg*AmberElec
    enddo
    write(*,'(2x,a,F10.4)') 'global charge= ',s
     call cpu_time(end_time)
    time_FF=end_time-start_time
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


!    print*,' running non-bonded engrad ..'
!    call cpu_time(start_time)
!    call nonbonded_amber_engrad(nfrag,fmol,evdw,ec)
!    call cpu_time(end_time)
!    time_nb=end_time-start_time
!    call prtime(6,time_nb,  'non-bonded engrad  ')
!
    allocate(ifrag(mol1%nat,9999))
    open(newunit=io,file='imff_ifrag',form='unformatted')
    read(io) ifrag
    close(io,status='delete')


    ! DEBUG DEBUG DEBUG DEBUG DEBUG
    ! gradient
!    do i=1,nfrag
!     do j=1,fmol(i)%nat
!      do k=1,3
!       mol1%g(k,ifrag(j,i))=fmol(i)%g(k,j)*au2ang ! in bohrs!
!       enddo
!     enddo
!    enddo


endif

!etotal=evdw+ec+ebond

! do i=1,mol1%nat
!   write(6,'(3E22.13)'), mol1%g(1:3,i)  ! ??
! enddo

!open(newunit=io,file='imff_gradient')
!  write(io,'(2x,F20.12)') etotal
!  do i=1,mol1%nat
!    write(io,'(3E22.13)'), mol1%g(1:3,i)*au2ang
!  enddo
!close(io)

!write(*,'(2x,F20.12)') etotal
! DEBUG DEBUG DEBUG DEBUG DEBUG

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



subroutine calc_gnorm(nat,g,gnorm)
implicit none
real(8) g(3,nat),gnorm,dnrm2
integer nat

gnorm=dnrm2(nat,g(1,1),1)
gnorm=gnorm+dnrm2(nat,g(2,1),1)
gnorm=gnorm+dnrm2(nat,g(3,1),1)

end
