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
                    ! First initial condition in the worksheet
                    uPrim(1,1:xDim/2) = 1.0
                    uPrim(1,((xDim/2)+1):xDim) = 0.125
                    uPrim(5,1:xDim/2) = 1.0/((5.0/3.0) - 1.0)
                    uPrim(5,((xDim/2)+1):xDim) = 0.1/((-1.0 + 5.0/3.0)*0.125)
                case (2)
                    ! A double shock tube discontinuity
                    uPrim(2,:) = 0.0
                    uPrim(3,:) = 0.0
                    uPrim(4,:) = 0.0
                    uPrim(5,:) = 0.5
                    uPrim(6,:) = 0.5
                    uPrim(1,1:xDim/3) = 0.5
                    uPrim(1,xDim/3+1:2*xDim/3) = 1.0
                    uPrim(1,(2*xDim/3+1):) = 0.5
                case default
                    ! A contact discontinuity
                    uPrim(2,:) = 0.5
                    uPrim(3,:) = 0.0
                    uPrim(4,:) = 0.0
                    uPrim(6,:) = 0.5
                    uPrim(1,1:xDim/2) = 1.0
                    uPrim(1,(xDim/2 +1):) = 0.5
                    uPrim(5,1:xDim/2) = 1.0
                    uPrim(5,(xDim/2 +1):) = 2.0
            end select
        end subroutine
end module setSod
