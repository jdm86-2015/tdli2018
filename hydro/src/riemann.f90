module riemann
    use reconstructor, only:reconstruct
    implicit none
    contains
        subroutine riemannSolver(uPrim, uFlux, uDim, xDim)
            ! calculates the fluxes at each cell interface.
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim
            integer :: iii

            real, intent(in), dimension(uDim,-1:xDim+2) :: uPrim
            ! uFlux matches dimension with uCons
            real, intent(out), dimension(uDim-2,-1:xDim+2) :: uFlux

            real, dimension(uDim,-1:xDim+2) :: uEdgeL
            real, dimension(uDim,-1:xDim+2) :: uEdgeR

            real :: fl
            real :: fr
            real :: denom
            real :: prod

            real, dimension(-1:xDim+2) :: vSqrR
            real, dimension(-1:xDim+2) :: vSqrL

            real :: lambdaP
            real :: lambdaM

            ! Begin Analysis

            ! Interpolate the value of the primary variables at cell edges.
            call reconstruct(uPrim,uEdgeL,uEdgeR,uDim,xDim,-1)

            ! get the squared velocities
            vSqrR = uEdgeR(2,:)*uEdgeR(2,:) + uEdgeR(3,:)*uEdgeR(3,:) + uEdgeR(4,:)*uEdgeR(4,:)
            vSqrL = uEdgeL(2,:)*uEdgeL(2,:) + uEdgeL(3,:)*uEdgeL(3,:) + uEdgeL(4,:)*uEdgeL(4,:)

            do iii=1,xDim
                ! calculate the characteristic speeds
                lambdaM = min(uEdgeL(2,iii) - uEdgeL(7,iii),uEdgeR(2,iii-1) - uEdgeR(7,iii-1))
                lambdaP = max(uEdgeL(2,iii) + uEdgeL(7,iii),uEdgeR(2,iii-1) + uEdgeR(7,iii-1))
                ! print*,lambdaM
                ! Implement fluxes given by the hydrodynamic equations
                if (lambdaM > 0) then
                    uFlux(1,iii) = uEdgeR(1,iii-1)*uEdgeR(2,iii-1)
                    uFlux(2,iii) = uEdgeR(1,iii-1)*vSqrR(iii-1) + uEdgeR(6,iii-1)
                    uFlux(5,iii) = (uEdgeR(5,iii-1) + & 
                        0.5*vSqrR(iii-1))*uEdgeR(1,iii-1)*uEdgeR(2,iii-1) + &
                        uEdgeR(6,iii-1)*uEdgeR(2,iii-1) 
                else if (lambdaP < 0) then
                    uFlux(1,iii) = uEdgeL(1,iii)*uEdgeL(2,iii)
                    uFlux(2,iii) = uEdgeL(1,iii)*vSqrL(iii) + uEdgeL(6,iii)
                    uFlux(5,iii) = (uEdgeL(5,iii) + &
                        0.5*vSqrL(iii))*uEdgeL(1,iii)*uEdgeL(2,iii) + &
                        uEdgeL(6,iii)*uEdgeL(2,iii)               
                else
                    denom = 1.0/(lambdaP - lambdaM)
                    prod = lambdaP*lambdaM
                    ! mass
                    fl = uEdgeR(1,iii-1)*uEdgeR(2,iii-1)
                    fr = uEdgeL(1,iii)*uEdgeL(2,iii)
                    uFlux(1,iii) = denom*(lambdaP*fl - lambdaM*fr + &
                        prod*(uEdgeL(1,iii) - uEdgeR(1,iii-1)))
                    ! momentum
                    fl = uEdgeR(1,iii-1)*vSqrR(iii-1) + uEdgeR(6,iii-1)
                    fr = uEdgeL(1,iii)*vSqrL(iii) + uEdgeL(6,iii)
                    uFlux(2,iii) = denom*(lambdaP*fl - lambdaM*fr + &
                        prod*(uEdgeL(2,iii) - uEdgeR(2,iii-1)))
                    ! energy
                    fl = (uEdgeR(5,iii-1) + &
                            0.5*vSqrR(iii-1))*uEdgeR(1,iii-1)*uEdgeR(2,iii-1) + &
                            uEdgeR(6,iii-1)*uEdgeR(2,iii-1)        
                    fr = (uEdgeL(5,iii) + 0.5*vSqrL(iii))*uEdgeL(1,iii)*uEdgeL(2,iii) + &
                        uEdgeL(6,iii)*uEdgeL(2,iii) 
                    uFlux(5,iii) = denom*(lambdaP*fl - lambdaM*fr + &
                        prod*(uEdgeL(1,iii) - uEdgeR(1,iii-1)))
                end if
            end do
            uFlux(:,-1) = uFlux(:,xDim-1)
            uFlux(:,0) = uFlux(:,xDim)
            uFlux(:,xDim+1) = uFlux(:,1)
            uFlux(:,xDim+2) = uFlux(:,2)
        end subroutine
end module riemann
        