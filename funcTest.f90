program funcTest
    implicit none

    interface
        function square(x, bins)
            integer,  intent(in) :: bins
            real, dimension(bins) :: square
            real, dimension(bins), intent(in) :: x 
        end function square
        function cube(x,bins)
            integer,  intent(in) :: bins
            real, dimension(bins) :: cube
            real, dimension(bins), intent(in) :: x 
        end function cube
        function power(x,bins,N)
            integer,  intent(in) :: bins
            integer,  intent(in) :: N
            real, dimension(bins) :: power
            real, dimension(bins), intent(in) :: x 
        end function power
    end interface

    integer, parameter :: b = 3
    integer :: i
    integer :: k

    real, dimension(b) :: r 
    real, dimension(b) :: y 
    real, dimension(b) :: z 
    real, dimension(b) :: q

    k = 3

    r = (/ 1,2,3 /)

    q = square(r,b)
    do i=1,b
        print*,q(i)
    end do
    print*,'***************************'
    q = cube(r,b)
    do i=1,b
        print*,q(i)
    end do
    print*,'***************************'
    q = power(r,b,k)
    do i=1,b
        print*,q(i)
    end do
end program funcTest

function square(x, bins)
    integer :: bins
    real, dimension(bins) :: x
    real, dimension(bins) :: square
    square = x*x
end function square

function cube(x, bins)
    interface
        function square(q, b)
            integer,  intent(in) :: b
            real, dimension(b) :: square
            real, dimension(b), intent(in) :: q 
        end function square
    end interface
    integer :: bins
    real, dimension(bins) :: x
    real, dimension(bins) :: cube
    cube = square(x,bins)
    cube = cube*x
end function cube

function power(x, bins, N)
    integer :: bins
    integer :: N
    real, dimension(bins) :: x
    real, dimension(bins) :: power
    power = x
    do n=1,N-1,1
        power = power*x
    end do
end function power
