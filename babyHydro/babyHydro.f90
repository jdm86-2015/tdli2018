program babyHydro
    implicit none
    interface
        function flux(u, uBins, xBins, gamma)
            integer,  intent(in) :: xBins
            integer,  intent(in) :: uBins
            real, intent(in) :: gamma
            real, dimension(uBins,xBins), intent(in) :: u
            real, dimension(uBins,xBins) :: flux
        end function flux
    end interface

    integer, parameter :: uB = 3
    integer, parameter :: xB = 1000

    integer :: i
    integer :: j

    real, dimension(uB,xB) :: uVector
    real, dimension(uB,xB) :: fVector

    real, parameter :: gam = 5.0/3.0

    do i = 1,uB
        do j = 1,xB
            uVector(i,j) = 10*(i-1) + j
        end do
    end do

    ! Compute the fluxes along the x-axis
    fVector = flux(uVector, uB, xB, gam)
    ! interpolate the fluxes at the midpoints with fVector(1) starting at 1/2
    fVector = 0.5*(fVector + cshift(fVector,1))

end program babyHydro


function flux(u, uBins, xBins, gamma)
    ! computes the flux object 
    integer :: xBins
    integer :: uBins
    real :: gamma
    real, dimension(uBins,xBins) :: u
    real, dimension(xBins) :: vSqr
    real, dimension(uBins, xBins) :: flux

    real :: pvCoeff

    pvCoeff = gamma/(gamma - 1.0)

    vSqr = u(2,:)*u(2,:)
    flux(1,:) = u(1,:)*u(2,:)
    flux(2,:) = u(1,:)*u(2,:)
    flux(3,:) = pvCoeff*(u(2,:)*u(3,:)) + 0.5*vSqr*u(1,:)
end function flux

    
