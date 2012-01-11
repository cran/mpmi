! This module is just to get the variable kind
! for double precision floats. 
module iface
    implicit none

    ! Double precision in R is real(kind=rdble)
    ! Used in all other functions
    integer, parameter :: rdble = kind(0.d0)
end module
