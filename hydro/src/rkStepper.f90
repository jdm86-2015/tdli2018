module rkStepper
    use equationOfState, only:eos
    use conservedVars, only:cons_calc
    use primatives, only:prim_calc
    use boundary, only:boundaries
    use sweepFunc, only:sweep
    implicit none
    contains
        subroutine rkStep(uPrim,uDim,xDim,deltaT,deltaX)
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim

            real, intent(in) :: deltaT
            real, intent(in) :: deltaX
            real, intent(inout), dimension(uDim,-1:xDim+2) :: uPrim

            real, dimension(uDim-2,-1:xDim+2) :: uCons
            real, dimension(uDim-2,-1:xDim+2) :: uConsP
            real, dimension(uDim-2,-1:xDim+2) :: uConsNext
            real, dimension(uDim-2,-1:xDim+2) :: dFlux

            real, parameter :: gamma = 5.0/3.0

            ! Calculate the equation of state
            call eos(uPrim,gamma,uDim,xDim)

            ! Calculate the conserved variables
            call cons_calc(uPrim,uCons,uDim,xDim)
            
            ! Populate the boundary condition for the conserved variables
            call boundaries(uCons,uDim-2,xDim,1)

            ! Calculate the fluxes using sweep
            call sweep(uPrim,dFlux,uDim,xDim,deltaX)

            ! Take a time step
            uConsP = uCons + deltaT*dFlux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THIS WILL OVERWRITE UPRIM AT THIS POINT.
            ! Populate the primary variables for this time step
            call prim_calc(uPrim,uConsP,uDim,xDim)

            ! calculate the equation of state for this time step
            call eos(uPrim,gamma,uDim,xDim)

            ! calculate the fluxes using the new uPrim and populate dFlux
            call sweep(uPrim,dFlux,uDim,xDim,deltaX)

            ! Take the final time step.
            uConsNext = 0.5*(uCons + uConsP + deltaT*dFlux)

            ! Write the solution into uPrim
            call prim_calc(uPrim,uConsNext,uDim,xDim)

            ! calculate the equation of state for the output primative variables.
            call eos(uPrim,gamma,uDim,xDim)

        end subroutine
end module rkStepper