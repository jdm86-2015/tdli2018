module riemann
    use reconstructor, only:reconstruct
    implicit none
    contains
        subroutine riemannSolver(uPrim, uFlux, uDim, xDim, method)
            ! calculates the fluxes at each cell interface.
            ! fluxes are defined with flux(i) being the left edge of cell i
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim
            integer, intent(in) :: method
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

            real :: ps
            real :: rl
            real :: rr
            real :: rul
            real :: rur
            real :: rvl
            real :: rvr
            real :: rwl
            real :: rwr
            real :: ul
            real :: ur
            real :: vl
            real :: vr
            real :: wl
            real :: wr
            real :: rel
            real :: rer

            real, dimension(-1:xDim+2) :: vSqrR
            real, dimension(-1:xDim+2) :: vSqrL

            real :: lambdaP
            real :: lambdaM
            real :: lambdaS

            lambdaP = 0.0
            lambdaM = 0.0
            lambdaS = 0.0

            uEdgeL(:,:) = 0.0
            uEdgeR(:,:) = 0.0

            vSqrR(:) = 0.0
            vSqrL(:) = 0.0
            
            ps = 0.0
            rl = 0.0
            rr = 0.0
            rul = 0.0
            rur = 0.0
            rvl = 0.0
            rvr = 0.0
            rwl = 0.0
            rwr = 0.0
            ul = 0.0
            ur = 0.0
            vl = 0.0
            vr = 0.0
            wl = 0.0
            wr = 0.0
            rel = 0.0
            rer = 0.0


            ! Begin Analysis

            ! Interpolate the value of the primary variables at cell edges.
            call reconstruct(uPrim,uEdgeL,uEdgeR,uDim,xDim,1)

            ! get the squared velocities
            vSqrR = uEdgeR(2,:)*uEdgeR(2,:) + uEdgeR(3,:)*uEdgeR(3,:) + uEdgeR(4,:)*uEdgeR(4,:)
            vSqrL = uEdgeL(2,:)*uEdgeL(2,:) + uEdgeL(3,:)*uEdgeL(3,:) + uEdgeL(4,:)*uEdgeL(4,:)

            select case (method)
            case (1)
                ! HLLC Solver
#ifdef NOTHING
                do iii=1,xDim
                    ! calculate the extreme characteristic speeds
                    lambdaM = min(uEdgeL(2,iii) - uEdgeL(7,iii),uEdgeR(2,iii-1) - uEdgeR(7,iii-1))
                    lambdaP = max(uEdgeL(2,iii) + uEdgeL(7,iii),uEdgeR(2,iii-1) + uEdgeR(7,iii-1))
                    ! Implement fluxes given by the hydrodynamic equations
                    denom =  uEdgeR(1,iii-1)*(lambdaM - uEdgeR(2,iii-1)) - &
                                uEdgeL(1,iii)*(lambdaP - uEdgeL(2,iii))
                    lambdaS = (1./denom)*(uEdgeR(6,iii-1) - uEdgeL(6,iii) + &
                                uEdgeR(2,iii-1)*uEdgeR(1,iii-1)*(lambdaM - uEdgeR(2,iii-1)) - &
                                uEdgeL(2,iii)*uEdgeL(1,iii)*(lambdaM - uEdgeL(2,iii)))

                    ps = uEdgeL(6,iii) + uEdgeL(1,iii)*(lambdaS-uEdge(2,iii))*(lambdaP - uEdge(2,iii))

                    rl = uEdgeR(1,iii-1)*(lambdaP-uEdgeR(2,iii-1))/(lambdaP - lambdaS)

                    rr = uEdgeL(1,iii)*(lambdaP-uEdgeL(2,iii))/(lambdaP - lambdaS)

                    rul = (uEdgeR(1,iii-1)*uEdgeR(2,iii-1)*(lambdaM - uEdgeR(2,iii-1)) + &
                        ps - uEdgeR(6,iii-1))/(lambdaM - lambdaS)
                    rur = (uEdgeL(1,iii)*uEdgeL(2,iii)*(lambdaP - uEdgeL(2,iii)) + &
                        ps - uEdgeL(6,iii))/(lambdaP - lambdaS)

                    rvl = (uEdgeR(1,iii-1)*uEdgeR(3,iii-1))*((lambdaM - uEdgeR(3,iii-1)) / &
                            (lambdaM - lambdaS))
                    rvr = (uEdgeL(1,iii)*uEdgeL(3,iii))*((lambdaP - uEdgeL(3,iii)) / &
                            (lambdaP - lambdaS))

                    rwl = (uEdgeR(1,iii-1)*uEdgeR(4,iii-1))*((lambdaM - uEdgeR(4,iii-1)) / &
                            (lambdaM - lambdaS))
                    rvr = (uEdgeL(1,iii)*uEdgeL(3,iii))*((lambdaP - uEdgeL(3,iii)) / &
                            (lambdaP - lambdaS))

                end do
#endif
            case default
                ! HLLE solver
                do iii=1,xDim
                    ! calculate the extreme characteristic speeds
                    lambdaM = min(uEdgeL(2,iii) - uEdgeL(7,iii),uEdgeR(2,iii-1) - uEdgeR(7,iii-1))
                    lambdaP = max(uEdgeL(2,iii) + uEdgeL(7,iii),uEdgeR(2,iii-1) + uEdgeR(7,iii-1))
                    ! Implement fluxes given by the hydrodynamic equations
                    if (lambdaM > 0) then
                        uFlux(1,iii) = uEdgeR(1,iii-1)*uEdgeR(2,iii-1)
                        uFlux(2,iii) = uEdgeR(1,iii-1)*uEdgeR(2,iii-1)*uEdgeR(2,iii-1) + uEdgeR(6,iii-1)
                        uFlux(5,iii) = ((uEdgeR(5,iii-1) + & 
                            0.5*vSqrR(iii-1))*uEdgeR(1,iii-1)*uEdgeR(2,iii-1) + &
                            uEdgeR(6,iii-1)*uEdgeR(2,iii-1))
                    else if (lambdaP < 0) then
                        uFlux(1,iii) = uEdgeL(1,iii)*uEdgeL(2,iii)
                        uFlux(2,iii) = uEdgeL(1,iii)*uEdgeL(2,iii)*uEdgeL(2,iii) + uEdgeL(6,iii)
                        uFlux(5,iii) = ((uEdgeL(5,iii) + &
                            0.5*vSqrL(iii))*uEdgeL(1,iii)*uEdgeL(2,iii) + &
                            uEdgeL(6,iii)*uEdgeL(2,iii))            
                    else
                        denom = 1.0/(lambdaP - lambdaM)
                        prod = lambdaP*lambdaM
                        ! mass
                        fl = uEdgeR(1,iii-1)*uEdgeR(2,iii-1)
                        fr = uEdgeL(1,iii)*uEdgeL(2,iii)
                        uFlux(1,iii) = denom*(lambdaP*fl - lambdaM*fr + &
                            prod*(uEdgeL(1,iii) - uEdgeR(1,iii-1)))
                        ! momentum
                        fl = uEdgeR(1,iii-1)*uEdgeR(2,iii-1)*uEdgeR(2,iii-1) + uEdgeR(6,iii-1)
                        fr = uEdgeL(1,iii)*uEdgeL(2,iii)*uEdgeL(2,iii) + uEdgeL(6,iii)
                        uFlux(2,iii) = denom*(lambdaP*fl - lambdaM*fr + &
                            prod*(uEdgeL(1,iii)*uEdgeL(2,iii) - uEdgeR(1,iii-1)*uEdgeR(2,iii-1)))
                        ! energy
                        fl = ((uEdgeR(5,iii-1) + &
                                0.5*vSqrR(iii-1))*uEdgeR(1,iii-1)*uEdgeR(2,iii-1) + &
                                uEdgeR(6,iii-1)*uEdgeR(2,iii-1))
                        fr = ((uEdgeL(5,iii) + 0.5*vSqrL(iii))*uEdgeL(1,iii)*uEdgeL(2,iii) + &
                            uEdgeL(6,iii)*uEdgeL(2,iii))
                        uFlux(5,iii) = denom*(lambdaP*fl - lambdaM*fr + &
                            prod*(uEdgeL(1,iii)*uEdgeL(5,iii) + 0.5*uEdgeL(1,iii)*vSqrL(iii) - &
                            (uEdgeR(1,iii-1)*uEdgeR(5,iii-1) + 0.5*uEdgeR(1,iii-1)*vSqrR(iii-1))))
                    end if
                end do
            end select
            ! uFlux(:,-1) = uFlux(:,xDim-1)
            ! uFlux(:,0) = uFlux(:,xDim)
            ! uFlux(:,xDim+1) = uFlux(:,1)
            ! uFlux(:,xDim+2) = uFlux(:,2)
            
        end subroutine
end module riemann
        