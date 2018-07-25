
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

  integer, allocatable :: ipair14(:,:,:,:)

  ! pair arrays 
  real(8), allocatable :: r0(:)
  real(8), allocatable :: rk(:)
  integer :: nbonds ! number of bonds
  integer, allocatable :: ibond(:,:)  ! atom identifier
  integer :: nangles ! number of valence angles
  integer, allocatable :: iangles  ! atom identifier
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
integer, parameter :: maxpar = 50, maxpar2 = maxpar**2
 type FFdata
  integer        :: npar ! # parameters
  character(120) :: id   ! FF name
  character(4) :: aname(maxpar) ! atom name array for vdw
  character(4) :: qname(maxpar) ! atom name array for charges
  real(8) :: chrg(maxpar) ! charge array
  real(8) :: LJe(maxpar)  !  LJ well-depth
  real(8) :: LJrad(maxpar) ! LJ radius
  integer :: itype(maxpar) ! integer atom type (if needed)
  integer :: nbondpar   ! number of bond parameters
  character(8) :: bond_name(maxpar2) ! atom pair name array for bonds
  real(8) :: r0(maxpar2) ! r0 for bonds
  real(8) :: rk(maxpar2) ! force constants for bonds
 end type
end module


module logic
logical echo ! printout
logical grad ! do gradient
character(120) filevec(5)
logical nchess
end module logic



module atomdata
real(8) ams(118)
data ams /1.00782503207,4.00260325415,7.016004548,9.012182201,11.009305406, &
12,14.00307400478,15.99491461956,18.998403224,19.99244017542, &
22.98976928087,23.985041699,26.981538627,27.97692653246,30.973761629, &
31.972070999,34.968852682,39.96238312251,38.963706679,39.962590983, &
44.955911909,47.947946281,50.943959507,51.940507472,54.938045141,  &
55.934937475,58.933195048,57.935342907,62.929597474,63.929142222,  &
68.925573587,73.921177767,74.921596478,79.916521271,78.918337087,  &
85.910610729,84.911789737,87.905612124,88.905848295,89.904704416,  &
92.906378058,97.905408169,98.906254747,101.904349312,102.905504292 , &
105.903485715,106.90509682,113.90335854,114.903878484,119.902194676, &
120.903815686,129.906224399,126.904472681,131.904153457,132.905451932, &
137.905247237,138.906353267,139.905438706,140.907652769,141.907723297, &
144.912749023,151.919732425,152.921230339,157.924103912,158.925346757, &
163.929174751,164.93032207,165.930293061,168.93421325,173.938862089, &
174.940771819,179.946549953,180.947995763,183.950931188,186.955753109, &
191.96148069,192.96292643,194.964791134,196.966568662,201.970643011, &
204.974427541,207.976652071,208.980398734,208.982430435,210.987496271, &
222.017577738,222.01755173,228.031070292,227.027752127,232.038055325, &
231.03588399,238.050788247,237.048173444,242.058742611,243.06138108, &
247.07035354,247.07030708,251.079586788,252.082978512,257.095104724, &
258.098431319,255.093241131,260.105504,263.112547,255.107398,259.114500, &
262.122892,263.128558,265.136151,281.162061,272.153615,283.171792,283.176451, &
285.183698,287.191186,292.199786,291.206564,293.214670/

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

end module



module constant
real(8), parameter:: pi = 3.141592653589793d0
real(8), parameter:: au2ang = 0.5291770d0
real(8), parameter:: amu2au=1.66053886E-27/9.10938215E-31
real(8), parameter:: au2cm =219474.63067d0
real(8), parameter :: AmberElec=18.2223d0 ! 
real(8), parameter:: au2kcal = 627.5095d0
end module

! If you measure charge in units of the electron charge, and distance in Angstroms (as is done in amber), then an electrostatic energy looks like:
!     E (kcal/mol) =  332 * q1*q2/r
! where q1 and q2 are charges and r is a distance. The square root of 332 is 18.2; hence, to save the multiplication by 332 all of the time, the charges are modified by 18.2, so that
!      E (kcal/mol) =  (18.2*q1) * (18.2*q2) / r
