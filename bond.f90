! determine bond matrix, coordination numbers and fragments
subroutine bondmat(nat,iat,xyz,aname,bond,cn,charge)
use atomdata, only: rcov
use moldata
use fragment
implicit none

integer nat,iat(nat),i,j,l,k,c
integer bond(nat,nat),cn(nat),atype(nat),cm(nat,nat),io
real(8) r1,r2,rab
real(8) xyz(3,nat),charge(nat)
character(2) esym
character(120) fmt,name
character(4) aname(nat)

logical exclude(nat),assigned
integer tmp_flen, tmp_nfrag,idx
logical write_xyz
logical debug

integer :: ifrag(nat,maxfrag)  !frag info array (can become LARGE!)
integer :: idn,idf  ! frag_atom_index, frag_index, number of fragments

debug=.false.
! debug=.true.
write_xyz=.false.

! full symmetric bond matrix
! somewhat expensive
!bond=0
!do i=1,nat
! do j=1,nat
!   if(i==j) cycle
!   r1=rcov(iat(i))+rcov(iat(j))*0.5
!   r2=rab(xyz(1,i),xyz(1,j))
!   if(abs(r1-r2) < 0.5) bond(i,j)=1
! enddo
!enddo


!  bond matrix (triangular -> full)
bond=0
do i=1,nat-1
 do j=i+1,nat
   if(i==j) cycle
   r1=rcov(iat(i))+rcov(iat(j))*0.5
   r2=rab(xyz(1,i),xyz(1,j))
   if(abs(r1-r2) < 0.5) then 
      bond(i,j)=1
      bond(j,i)=1
   endif
 enddo
enddo




!if(nat < 100) then
! call printimat(6,nat,nat,bond,'bond matrix')
!else
! print*,'print of bondmatrix omitted'
!endif

! coordination number (integer)
cn=0
do i=1,nat
 do j=1,nat
 if(bond(i,j)==1) then
  cn(i)=cn(i)+1
 endif
 enddo
enddo



! fragments
exclude=.false.
ifrag=0
flen=0
! 1st atom = 1st fragment
!idf=1
!idn=1
iidn=0
nfrag=0
!ifrag(idn,idf)=1
!flen(nfrag)=1
!exclude(idn)=.true.
!fragment infor array: ifrag(idn,idf)
! frag_atom_index(idn), frag_index(idf) = global atom indext (=stored value)
! 1 1 = 1
! 2 1 = 2
! 3 1 = 4
! 1 2 = 3

! count fragments and allocate dynamically...
! to-be-done


! assignt fragments
exclude(1)=.true.

! becomes expensive for many fragments (solvents!)
atomloop: do i=1,nat

 if(debug) then
  print*,' ATOM: ', i
  print*,'   flen',flen(1:nfrag)
 endif

  ! atom bound to a previous fragment?
  assigned=.false.
  fragloop: do idf=1,nfrag
	      if(debug)      print*,'    frag:',idf,flen(idf)
	       do k=1,flen(idf)
		  idx=ifrag(k,idf)
		  if(debug) print*,'      idf/idn/atom ',idf,k,idx
		  if(bond(i,idx)==1.and..not.assigned) then
		    if(debug) print*,'               --> bound to atom ',idx
		     assigned=.true.
		     iidn(idf)=iidn(idf)+1
		     idn=iidn(idf)
		     flen(idf)=flen(idf)+1
		     ifrag(idn,idf)=i
		     exclude(i)=.true.
		  endif
	       enddo
  enddo fragloop

  if(.not.assigned) then
  if(debug) print*,'--> nfrag+1 for atom ',i
  !idn=1
  nfrag=nfrag+1
  iidn(nfrag)=1
  idn=iidn(idf)
  if(nfrag>maxfrag) stop ' increase maxfrag!'
  flen(nfrag)=1
  ifrag(idn,nfrag)=i
  endif

  if(debug) call printimat(6,6,2,ifrag,'ifrag')
enddo atomloop


print*,'# fragments found',nfrag
! print*,flen(1 :nfrag)
do i=1,nfrag

! count bonds in fragments
do k=1,flen(i)-1
do j=k+1,flen(i)
    if (bond(ifrag(j,i),ifrag(k,i))==1) then
       fmol(i)%nbonds = fmol(i)%nbonds+1  !need to get num of bonds in each frag to allocate space for them
    endif
enddo
enddo

! print fragment info
print'(2x,a,I4)','Fragment: ',i
print'(2x,a,I4)','    size: ',flen(i)
print'(2x,a,I4)','    nbonds: ',fmol(i)%nbonds
   print'(2x,a)','    atoms : '
write(fmt,'("(3x",I5,x,"I5)")') flen(i)
write(*,fmt) ifrag(1:flen(i),i)
 

enddo

! assign fragment coordinates and fmol
print*,' '
print*, 'Largest fragment size: ',maxval(flen)

!allocate(fxyz(3,maxval(flen),nfrag))
!allocate(fgrad(3,maxval(flen),nfrag))
allocate(fiat(maxval(flen),nfrag))

! BUSY loop over all fragments
do i=1,nfrag
allocate(fmol(i)%xyz(3,flen(i)),fmol(i)%iat(flen(i)),fmol(i)%aname(flen(i)),fmol(i)%bond(flen(i),flen(i)))
allocate(fmol(i)%g(3,flen(i)))
allocate(fmol(i)%LJe(flen(i)),fmol(i)%LJrad(flen(i)),fmol(i)%chrg(flen(i)))
allocate(fmol(i)%r0(fmol(i)%nbonds),fmol(i)%rk(fmol(i)%nbonds))
allocate(fmol(i)%ibond(fmol(i)%nbonds,2))
 fmol(i)%nat=flen(i)
 fmol(i)%g=0d0
! ?? fmol(i)%nbonds = 0 
 do j=1,flen(i)
  ! use pre-assigned charges
  if(skip_charge) then
   fmol(i)%chrg(j)=charge(ifrag(j,i))
  endif
!   fxyz(1,j,i)=xyz(1,ifrag(j,i))
!   fxyz(2,j,i)=xyz(2,ifrag(j,i))
!   fxyz(3,j,i)=xyz(3,ifrag(j,i))
  fiat(j,i)=iat(ifrag(j,i))

   fmol(i)%xyz(1,j)=xyz(1,ifrag(j,i))
   fmol(i)%xyz(2,j)=xyz(2,ifrag(j,i))
   fmol(i)%xyz(3,j)=xyz(3,ifrag(j,i))
   fmol(i)%iat(j)=iat(ifrag(j,i))
   fmol(i)%aname(j)=aname(ifrag(j,i))

   do k=1,flen(i)
    fmol(i)%bond(j,k)=bond(ifrag(j,i),ifrag(k,i)) ! fill bond matrix for fragment
   enddo
 enddo
 call get_frag_primitives(fmol(i))

enddo

open(newunit=io,file='bff_ifrag',form='unformatted')
write(io) ifrag
close(io)

end


subroutine get_frag_primitives(mol)
use moldata
implicit none
integer, parameter :: maxcon=12
type(molecule) mol
integer :: c,i,j,k,l
integer, allocatable :: inat(:) ! number of attached atoms
integer, allocatable :: ipair12(:,:) ! id of attached atoms, max 8

allocate(inat(mol%nat),ipair12(mol%nat,maxcon))
! assign bonds for each fragment
! this uses local xyz-indices
mol%ibond=0
c=0
do i=1,mol%nat-1
 do j=i+1,mol%nat
   if(i==j) cycle
  if (mol%bond(j,i)==1) then
     c=c+1
     mol%ibond(c,1)=i
     mol%ibond(c,2)=j
  endif
 enddo
enddo

! 1-2 pair list for each fragment
ipair12=0
inat=0
do i=1,mol%nat
 do j=1,mol%nat
  if(i==j) cycle
  if(mol%bond(i,j)==1) then
    inat(i)=inat(i)+1
    ipair12(i,inat(i))=j
    if(inat(i)>maxcon) call error('too many connected atoms!')
  endif
 enddo
enddo

! angles for each fragment
! torsions for each fragment


end subroutine
