subroutine numhess(mol)
! numerical cartesian hessian
use moldata
use constant, only: au2ang
implicit none
integer i,j,ii,jj,hi,hj,yy,nat,nat3
real(8)  :: step,ec
real(8)  :: t1,t0,xx,thrE
real(8), allocatable ::  hess(:,:), gr(:,:),gl(:,:)
type(molecule) mol

! machine precision epsilon
thrE=epsilon(1d0)
nat3=mol%nat*3
nat=mol%nat
step=0.005d0
  print*, 'Doing Hessian numerically ..'
  ! print*, '#displacements: ',6*nat
  allocate(hess(nat3,nat3),gl(3,nat),gr(3,nat))
  hess=0.0d0
  do i=1,nat
  if(int(mod(i,nint(nat/10.0))).eq.0) write(6,'(I3,A)') nint(100d0*i/dble(nat)),' % done'
   do j=1,3
    hi=(i-1)*3+j
    mol%xyz(j,i)=mol%xyz(j,i)+step
    call eg_coulomb(mol%nat,mol%xyz,mol%iat,mol%iff,mol%chrg,ec,gr)
    ! call getgrad
    mol%xyz(j,i)=mol%xyz(j,i)-2d0*step
    ! call getgrad
    call eg_coulomb(mol%nat,mol%xyz,mol%iat,mol%iff,mol%chrg,ec,gl)
    mol%xyz(j,i)=mol%xyz(j,i)+step
  do ii=1,nat
   do jj=1,3
    hj=(ii-1)*3+jj
    xx=(gr(jj,ii)-gl(jj,ii))/(2d0*step)
    if(abs(xx).gt.thrE)  hess(hi,hj)=xx
   enddo ! jj-loop
  enddo  ! ii-loop

 enddo ! j-loop
 enddo  ! i-loop
  print*, ''
deallocate(gl,gr)

! symmetrize
  ! ! print*, 'Symmetrizing Hessian ..'
  ! do i=1,nat3
  !    do j=1,nat3
  !       chess(j,i)=(hess(i,j)+hess(j,i))*0.5d0
  !    enddo
  ! enddo

! write in bohrs
hess=hess*au2ang**2
call wrhess(nat3,hess,'FF.hess')

! call g98fake(nat,mol%iat,mol%xyz,1,1,1)

! do mass weighting+projection
! call Hmass(chess,'project')
!~ call DiagSM(nat3,chess,eig)
! and printing of Freq+ZPVE, and g98fake output


deallocate(hess)
end subroutine


! TM style hessian
subroutine wrhess(nat3,h,fname)
implicit none
integer nat3
real*8 h(nat3,nat3)
character*(*) fname
integer iunit,i,j,mincol,maxcol
character*5 adum
character*80 a80

adum='   '
iunit=11
open(unit=iunit,file=fname)
a80='$hessian'
write(iunit,'(a)')a80
do i=1,nat3
  maxcol = 0
  do while (maxcol.lt.nat3)
    mincol = maxcol + 1
    maxcol = min(maxcol+5,nat3)
    write(iunit,'(a5,5f15.10)')adum,(h(i,j),j=mincol,maxcol)
  enddo !if (maxcol.lt.nat3) goto 200
enddo
write(iunit,'(''$end'')')
close(iunit)
end




  subroutine g98fake(n,at,xyz,freq,u2,u)
  implicit none
  integer n,at(n)
  real*8 u(3*n,3*n),freq(3*n),xyz(3,n),u2(3*n,3*n)

  integer gu,i,j,ka,kb,kc,la,lb,k
  character*2 irrep
  real*8 red_mass(3*n)
  real*8 force   (3*n)
  real*8 ir_int  (3*n)
  real*8 f2      (3*n)
  real*8 zero

  irrep='a'
  red_mass=99.0
  force   =99.0
  ir_int  =99.0
  zero    =0.0

  k=0
  do i=1,3*n
     if(abs(freq(i)).gt.1.d-4)then
        k=k+1
        u(1:3*n,k)=u2(1:3*n,i)
        f2(k)=freq(i)
     endif
  enddo

  gu=55
  open(unit=gu,file='g98.out',status='unknown')
  write (gu,'('' Entering Gaussian System'')')
  write (gu,'('' *********************************************'')')
  write (gu,'('' Gaussian 98:'')')
  write (gu,'('' frequency fake output'')')
  write (gu,'('' *********************************************'')')

  write (gu,*) '                        Standard orientation:'
  write (gu,*) '---------------------------------------------', '-----------------------'
  write (gu,*) ' Center     Atomic     Atomic', '              Coordinates (Angstroms)'
  write (gu,*) ' Number     Number      Type ', '             X           Y           Z'
  write (gu,*) '-----------------------','---------------------------------------------'
  j=0
  do i=1,n
   write(gu,111) i,at(i),j,xyz(1:3,i)*0.52917
  enddo
  write (gu,*) '----------------------','----------------------------------------------'
  write (gu,*) '    1 basis functions        1 primitive gaussians'
  write (gu,*) '    1 alpha electrons        1 beta electrons'
  write (gu,*)
  111   format(i5,i11,i14,4x,3f12.6)

  write (gu,*) 'Harmonic frequencies (cm**-1), IR intensities',' (KM/Mole),'
  write (gu,*) 'Raman scattering activities (A**4/amu),',' Raman depolarization ratios,'
  write (gu,*) 'reduced masses (AMU), force constants (mDyne/A)',' and normal coordinates:'

  ka=1
  kc=3
  60    kb=min0(kc,k)
  write (gu,100) (j,j=ka,kb)
  write (gu,105) (irrep,j=ka,kb)
  write (gu,110) ' Frequencies --',(f2(j),j=ka,kb)
  write (gu,110) ' Red. masses --',(red_mass(j),j=ka,kb)
  write (gu,110) ' Frc consts  --',(force(j),j=ka,kb)
  write (gu,110) ' IR Inten    --',(ir_int(j),j=ka,kb)
  write (gu,110) ' Raman Activ --',(zero,j=ka,kb)
  write (gu,110) ' Depolar     --',(zero,j=ka,kb)
  write (gu,*)'Atom AN      X      Y      Z        X      Y','      Z        X      Y      Z'
  la=1
  70    lb=n
  do  i=la,lb
    write (gu,130) i,at(i), &
               (u(i*3-2,j),&
                u(i*3-1,j),&
                u(i*3  ,j),j=ka,kb)
  enddo
  if (lb.eq.n) go to 90
  go to 70
  90    if (kb.eq.k) then
   return
  endif
  ka=kc+1
  kc=kc+3
  go to 60

  100   format (3(20x,i3))
  105   format (3x,3(18x,a5))
  110   format (a15,f11.4,12x,f11.4,12x,f11.4)
  130   format (2i4,3(2x,3f7.2))

  close(gu)

  end
