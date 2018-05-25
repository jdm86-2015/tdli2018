module primatives
    implicit none
    contains
        subroutine prim_calc(uPrim,uCons,uDim,xDim)
            ! uDim is the dimensionality of the primitives
            ! uCDim is the dimension of the conserved quantities, it should be uDim - 2
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim

            real, intent(out), dimension(uDim,-1:xDim+2) :: uPrim
            real, intent(in), dimension(uDim-2,-1:xDim+2) :: uCons

            real, dimension(-1:xDim+2) :: vSqr

            ! uPrim == (rho, vx, vy, vz, e, p, cs)
            ! uCons == (rho, rho*vx, rho*vy, rho*vz, rho*e+0.5*rho*v^2)

            uPrim(1,1:xDim) = uCons(1,1:xDim)
            uPrim(2,1:xDim) = uCons(2,1:xDim)/uCons(1,1:xDim)
            uPrim(3,1:xDim) = uCons(3,1:xDim)/uCons(1,1:xDim)
            uPrim(4,1:xDim) = uCons(4,1:xDim)/uCons(1,1:xDim)

            vSqr = uPrim(2,:)*uPrim(2,:) + uPrim(3,:)*uPrim(3,:) + uPrim(4,:)*uPrim(4,:)

            uPrim(5,1:xDim) = (uCons(5,1:xDim) - 0.5*uPrim(1,1:xDim)*vSqr(1:xDim))/uPrim(1,1:xDim)
        end subroutine
end module primatives

