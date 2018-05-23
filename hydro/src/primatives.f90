module primatives
    implicit none
    contains
        subroutine prim_calc(uPrim,uCons,uDim,xDim)
            ! uDim is the dimensionality of the primitives
            ! uCDim is the dimension of the conserved quantities, it should be uDim - 2
            integer, intent(in) :: uDim
            integer, intent(in) :: xDim

            real, intent(inout), dimension(uDim,-1:xDim+2) :: uPrim
            real, intent(in), dimension(uDim-2,-1:xDim+2) :: uCons

            real, dimension(-1:xDim+2) :: vSqr

            ! uPrim == (rho, vx, vy, vz, e, p, cs)
            ! uCons == (rho, rho*vx, rho*vy, rho*vz, e+0.5*rho*v^2)

            uPrim(1,:) = uCons(1,:)
            uPrim(2,:) = uCons(2,:)/uCons(1,:)
            uPrim(3,:) = uCons(3,:)/uCons(1,:)
            uPrim(4,:) = uCons(4,:)/uCons(1,:)

            vSqr = uPrim(2,:)*uPrim(2,:) + uPrim(3,:)*uPrim(3,:) + uPrim(4,:)*uPrim(4,:)

            uPrim(5,:) = (uCons(5,:) - 0.5*uPrim(1,:)*vSqr)/uPrim(1,:)
        end subroutine
end module primatives

