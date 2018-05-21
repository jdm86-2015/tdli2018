module boundary
    implicit none
    contains
        subroutine boundaries(uPrim,uCons,uDim,xDim,bounds)
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim
            integer, intent(in) :: bounds
            
            real, intent(inout), dimension(uDim,-1:xDim+2) :: uPrim
            real, intent(inout), dimension(uDim-2,-1:xDim+2) :: uCons

            select case (bounds)
                case (1)
                    ! Periodic boundary condition
                    uPrim(:,-1) = uPrim(:,xDim-1)
                    uPrim(:,0) = uPrim(:,xDim)

                    uPrim(:,xDim+1) = uPrim(:,1)
                    uPrim(:,xDim+2) = uPrim(:,2)

                    uCons(:,-1) = uCons(:,xDim-1)
                    uCons(:,0) = uCons(:,xDim)

                    uCons(:,xDim+1) = uCons(:,1)
                    uCons(:,xDim+2) = uCons(:,2)
                case default
                    print*,"NO BOUNDARY CONDITION IMPLEMENTED." 
            end select
        end subroutine
end module boundary