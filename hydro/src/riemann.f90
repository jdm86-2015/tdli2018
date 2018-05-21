module riemann
    implicit none
    contains
        subroutine riemannSolver(uPrim, uFlux, uDim, xDim)
            ! calculates the fluxes at each cell interface.
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim
            integer :: iii

            real, intent(in), dimension(uDim,-1:xDim+2) :: uPrim
            real, intent(inout), dimension(uDim-2,-1:xDim+2) :: uFlux

            real :: fl
            real :: fr
            real :: denom
            real :: prod

            ! uFlux matches dimension with uCons
            
            real, dimension(uDim,-1:xDim+2) :: uL

            real, dimension(-1:xDim+2) :: vSqr
            real, dimension(-1:xDim+2) :: vSqrL

            real, dimension(-1:xDim+2) :: lambdaP
            real, dimension(-1:xDim+2) :: lambdaM

            uL = cshift(uPrim,1)

            vSqr = uPrim(2,:)*uPrim(2,:) + uPrim(3,:)*uPrim(3,:) + uPrim(4,:)*uPrim(4,:)
            vSqrL = cshift(vSqr,1)

            ! Begin Analysis
            lambdaM = min(uL(2,:) - uL(7,:),uPrim(2,:) - uPrim(7,:))
            lambdaP = max(uL(2,:) + uL(7,:),uPrim(2,:) + uPrim(7,:))

            do iii=1,xDim
                ! Implement fluxes given by the hydrodynamic equations
                if (lambdaM(iii) < 0) then
                    uFlux(1,iii) = uL(1,iii)*uL(2,iii)
                    uFlux(2,iii) = uL(1,iii)*vSqrL(iii) + uL(6,iii)
                    uFlux(5,iii) = (uL(5,iii) + 0.5*vSqrL(iii))*uL(1,iii)*uL(2,iii) + &
                        uL(6,iii)*uL(2,iii)                       
                else if (lambdaP(iii) > 0) then
                    uFlux(1,iii) = uPrim(1,iii)*uPrim(2,iii)
                    uFlux(2,iii) = uPrim(1,iii)*vSqr(iii) + uPrim(6,iii)
                    uFlux(5,iii) = (uPrim(5,iii) + 0.5*vSqr(iii))*uPrim(1,iii)*uPrim(2,iii) + &
                        uPrim(6,iii)*uPrim(2,iii)                    
                else
                    denom = 1.0/(lambdaP(iii) - lambdaM(iii))
                    prod = lambdaP(iii)*lambdaM(iii)
                    ! mass
                    fl = uL(1,iii)*uL(2,iii)
                    fr = uPrim(1,iii)*uPrim(2,iii)
                    uFlux(1,iii) = denom*(lambdaP(iii)*fl - lambdaM(iii)*fr + &
                        prod*(uPrim(1,iii) - uL(1,iii)))
                    ! momentum
                    fl = uL(1,iii)*vSqrL(iii) + uL(6,iii)
                    fr = uPrim(1,iii)*vSqr(iii) + uPrim(6,iii)
                    uFlux(2,iii) = denom*(lambdaP(iii)*fl - lambdaM(iii)*fr + &
                        prod*(uPrim(1,iii) - uL(1,iii)))
                    ! energy
                    fl = (uL(5,iii) + 0.5*vSqrL(iii))*uL(1,iii)*uL(2,iii) + &
                        uL(6,iii)*uL(2,iii)        
                    fr = (uPrim(5,iii) + 0.5*vSqr(iii))*uPrim(1,iii)*uPrim(2,iii) + &
                        uPrim(6,iii)*uPrim(2,iii) 
                    uFlux(5,iii) = denom*(lambdaP(iii)*fl - lambdaM(iii)*fr + &
                        prod*(uPrim(1,iii) - uL(1,iii)))
                end if
            end do
        end subroutine
end module
        