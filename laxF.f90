program laxF
    ! Configuration
    implicit none
    integer, parameter :: DP = SELECTED_REAL_KIND(14)

    interface
        function uNext(u,deltaT,deltaX,xBins)
            integer, intent(in) :: xBins
            real, intent(in) :: deltaT
            real, intent(in) :: deltaX
            real, dimension(xBins), intent(in) :: u
            real, dimension(xBins) :: uNext
        end function uNext
    end interface

    ! x axis objects
    integer :: J
    integer :: i
    integer, parameter :: xB = 100
    real :: xLength = 1.0
    real :: dX
    real, dimension(xB) :: xVec

    ! u object
    real, dimension(xB) :: u
    real, dimension(xB) :: up

    ! t objects
    real :: dT = 1e-2
    real :: tFinal = 0.02
    integer :: tSteps
    
    ! Define the delta X parameter
    dX = xLength/xB

    ! Build the initial condition
    xVec = (/ (J, J=0, xB-1) /)
    where (xVec < xB/2)
        u = 1.
    elsewhere
        u = 0.
    end where 

    ! figure out how many steps we need
    tSteps = ceiling(tFinal/dT)

    do i=1,tSteps,1
        u = uNext(u,dT,dX,xB)
    end do
    
    do i =1,xB,1
        print*,u(i)
    end do

end program laxF

function uNext(u, deltaT, deltaX, xBins)
    integer :: xBins
    real :: deltaT
    real :: deltaX
    real, dimension(xBins) :: u
    real, dimension(xBins) :: um1
    real, dimension(xBins) :: up1
    real, dimension(xBins) :: uNext

    up1 = cshift(u,1)
    um1 = cshift(u,-1)
    uNext = 0.5*(up1+um1) - 0.5*(deltaT/deltaX)*(up1 - um1)
end