module M_getvals
private
public getvals
contains
subroutine getvals(line,values,icount,ierr)
implicit none
character(len=*),parameter  :: ident='@(#)getvals: read arbitrary number of values from a character variable up to size of values'
! JSU 20170831
! http://fortranwiki.org/fortran/show/getvals
! fortranwiki.com

character(len=*),intent(in)  :: line
class(*),intent(in)          :: values(:)
integer,intent(out)          :: icount
integer,intent(out),optional :: ierr

   character(len=:),allocatable :: buffer
   character(len=len(line))     :: words(size(values))
   integer                      :: ios, i, ierr_local,isize

   select type(values)
   type is (integer);          isize=size(values)
   type is (real);             isize=size(values)
   type is (doubleprecision);  isize=size(values)
   type is (character(len=*)); isize=size(values)
   end select

   ierr_local=0

   words=' '                            ! make sure words() is initialized to null+blanks
   buffer=trim(line)//"/"               ! add a slash to the end so how the read behaves with missing values is clearly defined
   read(buffer,*,iostat=ios) words      ! undelimited strings are read into an array
   icount=0
   do i=1,isize                         ! loop thru array and convert non-blank words to numbers
      if(words(i).eq.' ')cycle

      select type(values)
      type is (integer);          read(words(i),*,iostat=ios)values(icount+1)
      type is (real);             read(words(i),*,iostat=ios)values(icount+1)
      type is (doubleprecision);  read(words(i),*,iostat=ios)values(icount+1)
      type is (character(len=*)); values(icount+1)=words(i)
      end select

      if(ios.eq.0)then
         icount=icount+1
      else
         ierr_local=ios
         write(*,*)'*getvals* WARNING:['//trim(words(i))//'] is not a number of specified type'
      endif
   enddo

   if(present(ierr))then
      ierr=ierr_local
   elseif(ierr_local.ne.0)then        ! error occurred and not returning error to main program to print message and stop program
      write(*,*)'*getval* error reading line ['//trim(line)//']'
      stop 2
   endif

end subroutine getvals
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
end module M_getvals
