program test
    use primatives, only:prim_calc
    use conservedVars, only:cons_calc
    use boundary, only:boundaries
    use equationOfState, only:eos
    implicit none
    integer, parameter :: uDim = 7
    integer, parameter :: xDim = 3
    integer, parameter :: uCDim = 5
    integer :: testFlag

    real, dimension(uDim,-1:xDim+2)  :: uPrim
    real, dimension(uCDim,-1:xDim+2) :: uCons
    real :: dummy

    real, parameter :: gam = 5.0/3.0

    print*,'Test program running...'

    testFlag = 0
    uPrim(:,:) = 0.0
    uCons(:,:) = 0.0

    uPrim(1,1:3) = (/ 0.25, 0.5, 0.75 /)
    uPrim(2,1:3) = (/ -0.25, 0.1, 1.2 /)
    uPrim(5,1:3) = (/ 0.5, 1.5, 2.0 /)

    print*,'Testing conserved quantity calculation....'
    call cons_calc(uPrim,uCons,uDim,xDim)
! 
    dummy = 1.5 + 0.5*0.1*0.1*0.5
    if(abs(uCons(5,2) - dummy) < 1.0e-14) then
        print*,'OK'
    else
        print*,'Conserved quantity calculation failed.'
    end if

    print*,'Testing primative quantity calculation....'
    uPrim(:,:) = 0.0
    call prim_calc(uPrim,uCons,uDim,xDim)
    if(abs(uPrim(5,2) - 1.5) < 1.0e-14) then
        print*,'OK'
    else
        print*,'Primative quantity calculation failed.'
    end if

    print*,'Testing boundary implementation....'
    call boundaries(uPrim,uCons,uDim,xDim,1)
    if(abs(uPrim(5,-1)-uPrim(5,xDim-1)) < 1.0e-14) then
        print*,'OK'
    else
        print*,'Boundary implementation failed.'
    end if

    dummy = (gam - 1.0)*1.5*0.5
    print*,'Testing equation of state....'
    call eos(uPrim,gam,uDim,xDim)
    if(abs(uPrim(6,2) - dummy) < 1.0e-14) then
        print*,'OK'
    else    
        print*,'Equation of state calculation failed'
    end if

end program test