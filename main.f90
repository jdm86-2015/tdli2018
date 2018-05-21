program hello
    implicit none
    integer i
    real, dimension(3) :: A
    real, dimension(3) :: B
    real, dimension(3) :: q

    A = (/ 1, 2, 3 /)
    B = (/ 0, 2, 7 /)

    q = max(A,B)

    do i=1,3
        print*,q(i)
    end do
end program hello
