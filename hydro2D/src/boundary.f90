module boundary
    implicit none
    contains
        subroutine boundaries(u,uDim,xDim,bounds)
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim
            integer, intent(in) :: bounds
            
            real, intent(inout), dimension(uDim,-1:xDim+2) :: u

            select case (bounds)
                case (1)
                    ! Periodic boundary condition
                    u(:,-1) = u(:,xDim-1)
                    u(:,0) = u(:,xDim)

                    u(:,xDim+1) = u(:,1)
                    u(:,xDim+2) = u(:,2)

                case default
                    print*,"NO BOUNDARY CONDITION IMPLEMENTED." 
            end select
        end subroutine
end module boundary