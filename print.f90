
subroutine print_header()

print*,''
print*,'|---------------------|'
print*,'| * Brno FF program * |'
print*,'|---------------------|'
!print*,'| model hessians +    |'
!print*,'| intermolecular      |'
!print*,'| rigid-body          |'
!print*,'| force fields        |'
!print*,'|---------------------|'
print*,''
print*,' ! ALPHA VERSION !'
print*,' H.Kruse               '
print*,' mail2holger@gmail.com '
print*,''
print*,'usage:'
print*,'  bff <structure>'
print*,''

call version

end subroutine
