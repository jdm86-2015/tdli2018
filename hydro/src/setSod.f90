module setSod
    implicit none
    contains
        subroutine sodIC(uPrim,uDim,xDim,set)
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim
            integer, intent(in) :: set

            real, intent(out), dimension(uDim,-1:xDim+2) :: uPrim

            select case (set)
                case (1)
                    uPrim(1,1:xDim/2) = 1.0
                    uPrim(1,((xDim/2)+1):xDim) = 0.125
                    uPrim(5,1:xDim/2) = 1.0/((5.0/3.0) - 1.0)
                    uPrim(5,((xDim/2)+1):xDim) = 0.1/((-1.0 + 5.0/3.0)*0.125)
                case default
                    print*,"No initial condition selected.  Initializing to zero."
                    uPrim(:,:) = 0.0
            end select
        end subroutine
end module setSod
