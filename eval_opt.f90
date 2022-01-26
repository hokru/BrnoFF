subroutine eval_options()
use logic
integer i,maxarg
character(80), allocatable :: arg(:)
character(120) infile,outfile
character(80) ftmp
logical fstr

maxarg=iargc()


if(maxarg.gt.0) then

allocate(arg(maxarg))
do i=1,maxarg
 call getarg(i,arg(i))
enddo
filevec(1)=arg(1)


 do i=1,maxarg
  ftmp=arg(i)
  if(index(ftmp,'-h ').ne.0) then
   print*,' help!..Help? HELP!!! ...h e l p ? !...  .'
  ! call help
  stop
  endif
  if(fstr(ftmp,'-nchess')) nchess=.true.
  if(fstr(ftmp,'-manual')) user_fragments=.true.
 enddo


endif


end subroutine
