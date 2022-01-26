! determine atom type from bond matrix and CN
subroutine get_atype(nat,iat,xyz,atype)
use atomdata, only: rcov
implicit none
integer nat,iat(nat),atype(nat),i,j
real(8) xyz(3,nat)
integer bond(nat,nat),cn(nat),io,ii,jj
real(8) r1,r2,rab
character(2) esym
character(30) atstring(25)




! full symmetric bond matrix
bond=0
do i=1,nat
 do j=1,nat
   if(i==j) cycle
   r1=rcov(iat(i))+rcov(iat(j))*0.5
   r2=rab(xyz(1,i),xyz(1,j))
!   write(*,'(3F7.2,2I2)') r1,r2,abs(r1-r2),i,j
   if(abs(r1-r2) < 0.5) bond(i,j)=1
 enddo
enddo

! print bond matrix to file
open(newunit=io,file='bondmat.bin')
  write(io,*) bond
close(io)


! call printimat(6,nat,nat,bond,'bond matrix')

! coordination number (integer)
cn=0
do i=1,nat
 do j=1,nat
 if(bond(j,i)==1) then
  cn(i)=cn(i)+1
 endif
 enddo
enddo

!write(*,'(a)')'CN:'
!do i=1,nat
! write(*,'(2x,I2,''['',a2,'']'',2x,I2)') i,esym(iat(i)),cn(i)
!enddo




! HYDROGEN
! 1 -->  H-R (general)
!     H-O(water)
!     H-O(hydroxy)
!     H-C(aliphatic)
!     H-N (+heavy analogues)
!     H-P (+heavy analogues)
!     H-S (+ "       "     )
!     H-M(alkali)
!     H-M(earth alkali)
!     H-TM
!     H-F
!     H-Cl (+ heavy analogues)
! CARBON
!2:   C(sp3) ; CN=4
!3:   C(sp2  ; CN=3
!4:   C(sp)  ; CN=2
! OXYGEN
!5:   O*-R2 ; CN=2
!6:   O*=C  ; CN=1
!10   O*=P  ; CN=1
!11   O*-P  ; CN=2
!12:  O*-H  ; CN=2
! NITROGEN
!7:   N*-R  ; CN=3
!8:   N*=R  ; CN=2
! PHOSPHATE
!9:   P*O4  ; CN=4
! 13,14 empty
!20 Na, 21=K 22=Mg 23=Ca 24=F 25=Cl
!15: Transition metal
!16: alkali
!17: earth alkali
!18: rare gases
!19: halogens
data atstring/ &
!      1           2        3          4     5       6       7      8
'H ','C(sp3) ','C(sp2)','C(sp)','O*-C2 ','O*=C','N*-R','N*=R', &
! 9      10      11    12     13     14      15     16       17
'P*O4','O*=P','O*-P','O*-H','XX','XX','TM','alkali (others)','earth alkali (others)', &
!  18          19
'rare gas','heavy halogens', &
! 20    21    22     23    24   25
 'Na+','K+','Mg2+','Ca2+','F','Cl' &
/


atype=0
write(*,'(a)')' Atom types:'
! assign atom type
do i=1,nat

if(esym(iat(i))=='h ') then
atype(i)=1
endif


if(esym(iat(i))=='c ') then
 if(cn(i)==2) atype(i)=4
 if(cn(i)==3) atype(i)=3
 if(cn(i)==4) atype(i)=2
endif

if(esym(iat(i))=='o ') then
 do j=1,nat
   if(bond(i,j)==1.and.iat(j)==15.and.cn(i)==1) atype=10
   if(bond(i,j)==1.and.iat(j)==15.and.cn(i)==2) atype=11
   if(bond(i,j)==1.and.iat(j)==6.and.cn(i)==1) atype=6
   if(bond(i,j)==1.and.iat(j)==6.and.cn(i)==2) atype=5
   if(bond(i,j)==1.and.iat(j)==1.and.cn(i)==2) atype=12
 enddo
endif

if(esym(iat(i))=='n ') then
 if(cn(i)==2) atype(i)=8
 if(cn(i)==3) atype(i)=7
endif


if(esym(iat(i))=='p ') then
 if(cn(i)==4) atype(i)=9
endif

if(esym(iat(i))=='k ') then
  atype(i)=17
 endif

write(*,'(i4,2x,a2,2x,i3,2x,a20)') i,esym(iat(i)),atype(i),atstring(atype(i))
enddo

end
