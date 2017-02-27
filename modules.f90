
module moldata
 ! holds data for a molecule
 type molecule
  integer, allocatable :: iat(:) ! atom type (H=1, C=6, etc)
  integer, allocatable :: iff(:) ! FF type type integer-index
  integer, allocatable :: bond(:,:) ! bond matrix
  integer, allocatable :: cn(:) ! CNs
  integer              :: nat    ! # atom
  character(4), allocatable :: aname(:) ! atom name array
  real(8), allocatable :: xyz(:,:) ! cartesian coordinates(3,nat)
  real(8), allocatable :: g(:,:) ! cartesian gradient(3,nat)
  real(8), allocatable :: dist(:) ! distance between all atoms
  real(8), allocatable :: chrg(:) ! charge array
  real(8), allocatable :: mass(:) ! masses array

  real(8), allocatable :: LJe(:)  !  LJ well-depth
  real(8), allocatable :: LJrad(:) ! LJ radius

  ! pair arrays 
  real(8), allocatable :: r0(:)
  real(8), allocatable :: rk(:)
 end type molecule
logical :: skip_charge
end module

module fragment


use moldata
integer, parameter :: maxfrag=9999

integer :: flen(maxfrag)       ! # atoms per fragment
integer :: iidn(maxfrag)       ! # atom index per fragment
integer :: nfrag  ! frag_atom_index, frag_index, number of fragments
real(8), allocatable :: fxyz(:,:,:) ! fragment coordinates size= (3,max(flen)),nfrag)
integer, allocatable :: fiat(:,:)
type(molecule) fmol(9999)

end module


module FFparm
integer, parameter :: maxpar = 50
 type FFdata
  integer        :: npar ! # parameters
  character(120) :: id   ! FF name
  character(4) :: aname(maxpar) ! atom name array for vdw
  character(4) :: qname(maxpar) ! atom name array for charges
  real(8) :: chrg(maxpar) ! charge array
  real(8) :: LJe(maxpar)  !  LJ well-depth
  real(8) :: LJrad(maxpar) ! LJ radius
  integer :: itype(maxpar) ! integer atom type (if needed)
 end type
end module


module logic
logical echo ! printout
logical grad ! do gradient
character(120) filevec(5)
logical nchess
end module logic


module radii
real(kind=8) rvdw(94),rcov(94)
! atomic radii from Mantina, Valero, Cramer, Truhlar "Atomic radii of elements"
! Copyed by hand,  maybe contain typos...
!            H       He
data rvdw /1.10d0,1.40d0, &
    ! Li     Be     B     C       N      O     F      Ne
    1.82d0,1.53d0,1.92d0,1.70d0,1.54d0,1.52d0,1.47d0,1.54d0, &
    ! Na    Mg     Al     Si     P      S       Cl     Ar
    2.27d0,1.73d0,1.84d0,2.10d0,1.80d0,1.80d0,1.75d0,1.88d0, &
    ! K      Ca     Sc     Ti     V      Cr      Mn     Fe     Co    Ni     Cu     Zn     Ga     Ge     As     Se     Br    Kr
    2.75d0,2.31d0,2.15d0,2.11d0,2.07d0,2.06d0,2.05d0,2.04d0,2.00d0,1.97d0,1.96d0,2.01d0,1.87d0,2.11d0,1.85d0,1.90d0,1.85d0,2.02d0, &
    !  Rb    Sr     Y      Zr      Nb     Mo    Tc     Ru     Rh     Pd     Ag     Cd     In    Sn      Sb      Te     I     Xe
    3.03d0,2.49d0,2.26d0,2.23d0,2.18d0,2.17d0,2.16d0,2.13d0,2.10d0,2.10d0,2.11d0,2.18d0,1.93d0,2.17d0,2.06d0,2.06d0,1.98d0,2.16d0, &
    ! Cs Ba
    3.32d0,2.68d0, &
    ! La-Lu
    2.43d0,2.42d0,2.40d0,2.46d0,2.38d0,2.36d0,2.35d0,2.34d0,2.33d0,2.31d0,2.30d0,2.29d0,2.27d0,2.26d0,2.24d0, &
    ! Hf     Ta     W       Re     Os    Ir     Pt     Au     Hg     Ti     Pb     Bi     Po     At     Rn
    2.23d0,2.22d0,2.18d0,2.16d0,2.16d0,2.13d0,2.13d0,2.23d0,2.23d0,2.11d0,2.02d0,2.07d0,1.97d0,2.02d0,2.20d0, &
    ! Fr-Pu
    3.48d0,2.83d0,2.47d0,2.45d0,2.43d0,2.41d0,2.39d0,2.43d0/


!            H       He
data rcov /0.32d0,0.37d0, &
    ! Li     Be     B     C       N      O     F      Ne
    1.30d0,0.99d0,0.84d0,0.75d0,0.71d0,0.64d0,0.60d0,0.62d0,  &
    ! Na    Mg     Al     Si     P      S       Cl     Ar
    1.60d0,1.40d0,1.24d0,1.14d0,1.09d0,1.04d0,1.00d0,1.01d0,  &
    ! K      Ca     Sc     Ti     V      Cr      Mn     Fe     Co    Ni     Cu     Zn     Ga     Ge     As     Se     Br    Kr
    2.00d0,1.74d0,1.59d0,1.48d0,1.44d0,1.30d0,1.29d0,1.24d0,1.18d0,1.17d0,1.22d0,1.20d0,1.23d0,1.20d0,1.20d0,1.18d0,1.17d0,1.24d0, &
    !  Rb    Sr     Y      Zr      Nb     Mo    Tc     Ru     Rh     Pd     Ag     Cd     In    Sn      Sb      Te     I     Xe
    2.15d0,1.90d0,1.78d0,1.64d0,1.56d0,1.46d0,1.38d0,1.36d0,1.34d0,1.30d0,1.36d0,1.40d0,1.42d0,1.40d0,1.40d0,1.37d0,1.32d0,1.36d0, &
    ! Cs Ba
    2.38d0,2.06d0,  &
    ! La-Lu
     1.94d0,1.84d0,1.90d0,1.73d0,1.86d0,1.85d0,1.83d0,1.82d0,1.81d0,1.80d0,1.79d0,1.77d0,1.77d0,1.78d0,1.74d0,  &
    ! Hf     Ta     W       Re     Os    Ir     Pt     Au     Hg     Ti     Pb     Bi     Po     At     Rn
    1.64d0,1.58d0,1.50d0,1.41d0,1.36d0,1.32d0,1.30d0,1.64d0,1.88d0,1.48d0,1.45d0,1.50d0,1.42d0,1.47d0,1.46d0,  &
    ! Fr-Pu
    2.42d0,2.11d0,2.01d0,1.90d0,1.84d0,1.83d0,1.80d0,1.80d0/



end module radii


module constant
real(8), parameter:: pi = 3.141592653589793d0
real(8), parameter:: au2ang = 0.5291770d0
real(8), parameter:: amu2au=1.66053886E-27/9.10938215E-31
real(8), parameter:: au2cm =219474.63067d0
end module
