module sweepFunc
    ! Input the value of the primative variabled and the flux derivative and the grid spacing.
    use boundary, only:boundaries
    use riemann, only:riemannSolver
    implicit none
    contains
        subroutine sweep(uPrim,dFlux,uDim,xDim,deltaX)
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim

            real, intent(inout), dimension(uDim,-1:xDim+2) :: uPrim
            real, intent(out), dimension(uDim-2,-1:xDim+2) :: dFlux

            real, dimension(uDim-2,-1:xDim+2) :: uFlux

            real, intent(in) :: deltaX

            ! ensure boundary conditions for uPrim
            call boundaries(uPrim,uDim,xDim,1)
            ! calculate the fluxes by solving the riemann problem.
            call riemannSolver(uPrim,uFlux,uDim,xDim)

            ! Should do flux(i + 1/2) - flux(i-1/2)
            dFlux = (1./deltaX)*(cshift(uFlux,1,2) - uFlux)
        end subroutine
end module sweepFunc


