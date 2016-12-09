! intermolecular FF 

module parm
 ! holds data for a molecule

 type molecule
  integer, allocatable :: iat(:) ! atom type (H=1, C=6, etc)
  integer, allocatable :: iff(:) ! FF type type
  integer,             :: nat    ! # atom

  character(4), allocatable :: aname(:) ! atom name array

  real(8), allocatable :: xyz(:,:) ! cartesian coordinates(3,nat)
  real(8), allocatable :: dist(:) ! distance between all atoms
  real(8), allocatable :: chrg(:) ! charge array
  real(8), allocatable :: mass(:) ! masses array
  
 end type molecule

end module



subroutine interFF
implicit none





end subroutine
