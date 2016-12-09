subroutine eval_options()
use logic
integer i,maxarg
character(80), allocatable :: arg(:)
character(120) infile,outfile
character(80) ftmp

maxarg=iargc()
if(maxarg.gt.0) then

 allocate(arg(maxarg))
   do i=1,maxarg
     call getarg(i,arg(i))
   enddo

 filevec(1)=arg(1)
 ! filevec(2)=arg(2)

 do i=1,maxarg
  ftmp=arg(i)
  if(index(ftmp,'-h ').ne.0) then
   print*,' help!..Help? HELP!!! ...h e l p ? !...  .'
  ! call help
  stop
  endif
  if(index(ftmp,'-noprint ').ne.0) echo=.false.
  if(index(ftmp,'-grad ').ne.0) grad=.true.
 enddo


endif


end subroutine
