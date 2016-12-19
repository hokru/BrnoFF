
subroutine print_header()

print*,''
print*,'|^^^^^^^^^^^^^^^^^^^^^|'
print*,'|      I M F F        |'
print*,'| intermolecular      |'
print*,'| rigid-body          |'
print*,'| force field         |'
print*,'|---------------------|'
print*,''
print*,' ! ALPHA VERSION !'
print*,' H.Kruse               '
print*,' mail2holger@gmail.com '
print*,''
print*,'usage:'
print*,'  imff <structure>'
print*,''

call version

end subroutine
