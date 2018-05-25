module setupGrid
    implicit none
    contains
        subroutine gridValues(xAxis,xLength,xDim,deltaX)
            integer, intent(in) :: xDim
            integer :: i

            real, intent(inout), dimension(xDim) :: xAxis
            real, intent(inout) :: deltaX
            real, intent(in) :: xLength

            deltaX = xLength/real(xDim)

            do i=1,xDim
                xAxis(i) = (real(i)-1.0)*deltaX
            end do
        end subroutine
end module setupGrid