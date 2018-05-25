module timeStepCalc
    implicit none
    contains
        subroutine timeStep(uPrim,uDim,xDim,deltaX,deltaT,safety)
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim
            
            real, intent(in), dimension(uDim,-1:xDim+2) :: uPrim
            real, dimension(-1:xDim+2) :: v

            real, intent(in) :: safety
            real, intent(in) :: deltaX
            real, intent(out) :: deltaT
            real :: lambda

            v = sqrt(uPrim(2,:)*uPrim(2,:) + uPrim(3,:)*uPrim(3,:) + uPrim(4,:)*uPrim(4,:))

            lambda = maxval(v + uPrim(7,:))

            deltaT = safety*(deltaX/lambda)
        end subroutine
end module timeStepCalc
            