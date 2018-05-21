program ltest
    implicit none
    integer, parameter :: DP = selected_real_kind(14)
    real(kind = DP) :: q = 22.0, r = 22.0, s = 15.0
    if (q < r) then
        print*,'q is less than r'
    else if (q > r) then
        print*,'q is greater than r'
    else 
        print*,'q is equal to r'
    end if
end program ltest