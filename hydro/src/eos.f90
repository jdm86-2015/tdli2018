module equationOfState
    implicit none
    contains
        subroutine eos(uPrim,gamma,uDim,xDim)
            integer, intent(in) :: xDim
            integer, intent(in) :: uDim
            
            real, intent(in) :: gamma
            real, intent(inout), dimension(uDim,-1:xDim+2) :: uPrim

            uPrim(6,:) = (gamma-1.0)*uPrim(1,:)*uPrim(5,:)
            uPrim(7,:) = sqrt(gamma*(gamma-1.0)*uPrim(5,:))
        end subroutine
end module equationOfState