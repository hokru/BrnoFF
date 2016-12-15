! determine bond matrix, coordination numbers and fragments
subroutine bondmat(nat,iat,xyz,aname,bond,cn)
use radii
use moldata
use fragment
implicit none

integer nat,iat(nat),i,j,l,k
integer bond(nat,nat),cn(nat),atype(nat),cm(nat,nat),io
real(8) r1,r2,rab
real(8) xyz(3,nat)
character(2) esym
character(120) fmt,name
character(4) aname(nat)

logical exclude(nat),assigned
integer tmp_flen, tmp_nfrag,idx
logical write_xyz
logical debug

integer :: ifrag(nat,maxfrag)  !frag info array
integer :: idn,idf  ! frag_atom_index, frag_index, number of fragments

debug=.false.
! debug=.true.
write_xyz=.false.

! full symmetric bond matrix
bond=0
do i=1,nat
 do j=1,nat
   if(i==j) cycle
   r1=rcov(iat(i))+rcov(iat(j))*0.5
   r2=rab(xyz(1,i),xyz(1,j))
   if(abs(r1-r2) < 0.5) bond(i,j)=1
 enddo
enddo

if(nat < 100) then
 call prmati(6,bond,nat,nat,'bond matrix')
else
 print*,'print of bondmatrix omitted'
endif

! coordination number (integer)
cn=0
do i=1,nat
 do j=1,nat
 if(bond(i,j)==1) then
  cn(i)=cn(i)+1
 endif
 enddo
enddo


! write(6,'(a)') ' Coordination Numbers '
! write(6,'(a)') '   index   CN'
! do i=1,nat
!  write(6,'(2x,I8,x,I3)') i,cn(i)
! enddo

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
  !                   idn=idn+1
		     iidn(idf)=iidn(idf)+1
		     idn=iidn(idf)
		     flen(idf)=flen(idf)+1
  !                   print*,'X',idn,flen(idf)
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

  if(debug) call prmati(6,ifrag,6,2,'ifrag')
enddo atomloop


print*,'# fragments found',nfrag
! print*,flen(1:nfrag)
do i=1,nfrag
print'(2x,a,I4)','Fragment: ',i
print'(2x,a,I4)','    size: ',flen(i)
print'(2x,a)','    atoms: '
write(fmt,'("(3x",I5,x,"I5)")') flen(i)
write(*,fmt) ifrag(1:flen(i),i)
enddo

! assign fragment coordinates and fmol
print*,' '
print*, 'Largest fragment size: ',maxval(flen)

allocate(fxyz(3,maxval(flen),nfrag))
allocate(fiat(maxval(flen),nfrag))


do i=1,nfrag
allocate(fmol(i)%xyz(3,flen(i)),fmol(i)%iat(flen(i)),fmol(i)%aname(flen(i)))
 fmol(i)%nat=flen(i)
 do j=1,flen(i)
   fxyz(1,j,i)=xyz(1,ifrag(j,i))
   fxyz(2,j,i)=xyz(2,ifrag(j,i))
   fxyz(3,j,i)=xyz(3,ifrag(j,i))
   fiat(j,i)=iat(ifrag(j,i))

   fmol(i)%xyz(1,j)=xyz(1,ifrag(j,i))
   fmol(i)%xyz(2,j)=xyz(2,ifrag(j,i))
   fmol(i)%xyz(3,j)=xyz(3,ifrag(j,i))
   fmol(i)%iat(j)=iat(ifrag(j,i))
   fmol(i)%aname(j)=aname(ifrag(j,i))
 enddo
write(name,'(I4,a)') i,'_frag.xyz'
if(write_xyz)call wrxyz(fiat(1,i),flen(i),fxyz(1,1,i),adjustl(trim(name)))
enddo



end
