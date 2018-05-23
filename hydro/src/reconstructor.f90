module reconstructor
    implicit none
    contains
        subroutine reconstruct(uPrim, uEdgeL, uEdgeR, uDim, xDim, method)
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim
            integer, optional, intent(in) :: method
            
            real, intent(in), dimension(uDim,-1:xDim+2) :: uPrim
            real, intent(out), dimension(uDim,-1:xDim+2) :: uEdgeL
            real, intent(out), dimension(uDim,-1:xDim+2) :: uEdgeR

            real, dimension(uDim) :: sL
            real, dimension(uDim) :: sR

            integer :: method_
            real :: dUdX

            integer :: i
            integer :: j

            if ( present(method) ) then
                method_ = method
            else
                method_ = 1
            end if

            select case (method_)
                case (1)
                    ! Linear reconstruction
                    do i = 1,xDim
                        sL(:) = (uPrim(:,i) - uPrim(:,i-1))
                        sR(:) = (uPrim(:,i+1) - uPrim(:,i))
                        do j = 1,uDim
                            if ((sL(j) > 0.0) .AND. (sR(j) > 0.0)) then
                                dUdX = min(abs(sL(j)),abs(sR(j)))
                            else if ((sL(j) < 0.0) .AND. (sR(j) < 0.0)) then
                                dUdX = -min(abs(sL(j)),abs(sR(j)))
                            else
                                dUdX = 0.0
                            end if
                            uEdgeL(j,i) = uPrim(j,i) - 0.5*dUdX
                            uEdgeR(j,i) = uPrim(j,i) + 0.5*dUdX
                        end do
                    end do
                    uEdgeL(:,-1) = uEdgeL(:,xDim-1)
                    uEdgeR(:,-1) = uEdgeR(:,xDim-1)
                    uEdgeL(:,0) = uEdgeL(:,xDim)
                    uEdgeR(:,0) = uEdgeR(:,xDim)
                    uEdgeL(:,xDim+1) = uEdgeL(:,1)
                    uEdgeR(:,xDim+1) = uEdgeR(:,1)
                    uEdgeL(:,xDim+2) = uEdgeL(:,2)
                    uEdgeR(:,xDim+2) = uEdgeR(:,2)
                case default
                    ! No reconstruction at all.
                    uEdgeL = uPrim
                    uEdgeR = uPrim
            end select
        end subroutine
end module reconstructor


