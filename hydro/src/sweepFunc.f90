module sweepFunc
    use boundary, only:boundaries
    use riemann, only:riemannSolver
    implicit none
    contains
        subroutine sweepFunc(uPrim,uDim,xDim)
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim

            real, intent(inout), dimension(uDim,-1:xDim+2) :: uPrim

            real, dimension(uDim-2,-1:xDim+2) :: uCons
            real, dimension(uDim-2,-1:xDim+2) :: uFlux
            real, dimension(uDim-2,-1:xDim+2) :: dFlux


            call boundaries(uPrim,uCons,uDim,xDim,1)
            call riemann(uPrim,uFlux,uDim,xDim)

            ! Should do flux(i + 1) - flux(i)
            dFlux = (cshift(uFlux,1,2) - uFlux)


