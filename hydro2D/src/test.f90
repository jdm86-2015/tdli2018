program test
    ! Imports
    use primatives, only:prim_calc
    use conservedVars, only:cons_calc
    use boundary, only:boundaries
    use equationOfState, only:eos
    use reconstructor, only:reconstruct
    use riemann, only:riemannSolver
    use sweepFunc, only:sweep
    use timeStepCalc, only:timeStep
    use rkStepper, only:rkStep
    use setupGrid, only:gridValues
    use setSod, only:sodIC

    implicit none

    integer, parameter :: uDim = 7
    integer, parameter :: xDim = 3
    integer, parameter :: uCDim = 5
    integer :: testFlag

    real, dimension(uDim,-1:xDim+2) :: uPrim
    real, dimension(uDim,-1:xDim+2) :: uEdgeL
    real, dimension(uDim,-1:xDim+2) :: uEdgeR

    real, dimension(uCDim,-1:xDim+2) :: uCons
    real, dimension(uCDim,-1:xDim+2) :: uFlux

    real, dimension(xDim) :: xAxis

    real :: dummy
    real :: deltaT
    real :: deltaX

    real, parameter :: gam = 5.0/3.0

    print*,'Test program running...'

    deltaX = 0.1

    ! ********************** CONSERVED TEST ********************** 
    testFlag = 0
    uPrim(:,:) = 0.0
    uCons(:,:) = 0.0

    uPrim(1,1:3) = (/ 0.25, 0.5, 0.75 /)
    uPrim(2,1:3) = (/ -0.25, 0.1, 1.2 /)
    uPrim(5,1:3) = (/ 0.5, 1.5, 2.0 /)

    print*,'Testing conserved quantity calculation....'
    call cons_calc(uPrim,uCons,uDim,xDim)
! 
    dummy = (0.5*1.5 + 0.5*0.1*0.1*0.5)
    if(abs(uCons(5,2) - dummy) < 1.0e-14) then
        print*,'OK'
    else
        print*,'Conserved quantity calculation failed.'
    end if

    ! ********************** PRIMATIVE TEST ********************** 

    print*,'Testing primative quantity calculation....'
    uPrim(:,:) = 0.0
    call prim_calc(uPrim,uCons,uDim,xDim)
    if(abs(uPrim(5,2) - 1.5) < 1.0e-14) then
        print*,'OK'
    else
        print*,'Primative quantity calculation failed.'
    end if

    ! ********************** BOUNDARY TEST ********************** 

    print*,'Testing boundary implementation....'
    call boundaries(uPrim,uDim,xDim,1)
    if(abs(uPrim(5,-1)-uPrim(5,xDim-1)) < 1.0e-14) then
        print*,'OK'
    else
        print*,'Boundary implementation failed.'
    end if

    ! ********************** EOS TEST ********************** 
    dummy = (gam - 1.0)*1.5*0.5
    print*,'Testing equation of state....'
    call eos(uPrim,gam,uDim,xDim)
    if(abs(uPrim(6,2) - dummy) < 1.0e-14) then
        print*,'OK'
    else    
        print*,'Equation of state calculation failed'
    end if

    ! ********************** RECONSTRUCTOR TEST ********************** 

    testFlag = 0
    uPrim(:,:) = 0.0
    uCons(:,:) = 0.0
    uEdgeL(:,:) = 0.0
    uEdgeR(:,:) = 0.0

    uPrim(1,1:3) = (/ 2.0, 4.0, 6.0 /)

    print*,'Testing reconstructor....'
    call boundaries(uPrim,uDim,xDim,1)
    call reconstruct(uPrim,uEdgeL,uEdgeR,uDim,xDim,1)

    dummy = 5
    print*,'uPrim(1,:): '
    print*,uPrim(1,:)

    print*,'uEdgeR(1,:): '
    print*,uEdgeR(1,:)

    print*,'uEdgeL(1,:): '
    print*,uEdgeL(1,:)

    if(abs(uEdgeR(1,2) - dummy) < 1.0e-14) then
        print*,'Right Edge OK'
    else    
        print*,'Right edge calculation failed'
    end if
    dummy = 3
    if(abs(uEdgeL(1,2) - dummy) < 1.0e-14) then
        print*,'Left Edge OK'
    else    
        print*,'Left edge calculation failed'
    end if

    ! ********************** RECONSTRUCTOR TEST ********************** 

    testFlag = 0
    uPrim(:,:) = 0.0
    uCons(:,:) = 0.0
    uEdgeL(:,:) = 0.0
    uEdgeR(:,:) = 0.0

    uPrim(1,1:3) = (/ 2.0, 4.0, 6.0 /)

    print*,'Testing reconstructor with no reconstruction enabled....'
    call boundaries(uPrim,uDim,xDim,1)
    call reconstruct(uPrim,uEdgeL,uEdgeR,uDim,xDim,-1)

    print*,'uPrim(1,:): '
    print*,uPrim(1,:)

    print*,'uEdgeR(1,:): '
    print*,uEdgeR(1,:)

    print*,'uEdgeL(1,:): '
    print*,uEdgeL(1,:)

    dummy = uPrim(1,2)
    if(abs(uEdgeR(1,2) - dummy) < 1.0e-14) then
        print*,'Right Edge OK'
    else    
        print*,'Right edge calculation failed'
    end if
    dummy = uPrim(1,2)
    if(abs(uEdgeL(1,2) - dummy) < 1.0e-14) then
        print*,'Left Edge OK'
    else    
        print*,'Left edge calculation failed'
    end if


    ! ********************** RIEMANN SOLVER TEST ********************** 

    testFlag = 0
    uPrim(:,:) = 0.0
    uCons(:,:) = 0.0
    uFlux(:,:) = 0.0

    uEdgeL(:,:) = 0.0
    uEdgeR(:,:) = 0.0

    uPrim(1,1:3) = (/ 3.0, 3.0, 3.0 /)
    uPrim(2,1:3) = (/ 0.0, 0.0, 0.0 /)
    uPrim(5,1:3) = (/ 0.1, 0.1, 0.1 /)
    uPrim(6,1:3) = (/ 0.9, 0.9, 0.9 /)
    uPrim(7,1:3) = (/ 7.0, 7.0, 7.0 /)
    print*,'Testing Riemann Solver with constant input....'
    call boundaries(uPrim,uDim,xDim,1)
    ! call eos(uPrim,gam,uDim,xDim)
    call reconstruct(uPrim,uEdgeL,uEdgeR,uDim,xDim,1)
    call riemannSolver(uPrim,uFlux,uDim,xDim)

    print*,'uPrim(1,:): '
    print*,uPrim(1,:)
    print*,'uEdgeR(1,:): '
    print*,uEdgeR(1,:)
    print*,'uEdgeL(1,:): '
    print*,uEdgeL(1,:)
    print*,'uFlux(1,:): '
    print*,uFlux(1,:)
    print*,'uFlux(2,:): '
    print*,uFlux(2,:)
    print*,'uFlux(5,:): '
    print*,uFlux(5,:)

    ! ********************** SWEEP TEST ********************** 

    print*,'Testing Sweep with constant input....'

    dummy = 0.0
    call sweep(uPrim,uFlux,uDim,xDim,deltaX)    
    if(abs(uFlux(1,2) - dummy) < 1.0e-14) then
        print*,'OK'
    else    
        print*,'Sweep failed'
    end if

    ! ********************** SWEEP TEST ********************** 

    print*,'Testing time step....'

    dummy =  0.1/7.0
    call timeStep(uPrim, uDim, xDim, deltaX, deltaT, 1.0)
    if(abs(deltaT - dummy) < 1.0e-14) then
        print*,'OK'
    else    
        print*,'Time step failed'
    end if

    ! ********************** SWEEP TEST ********************** 

    print*,'Testing rk step with constant variables....'

    dummy = 3.0 
    call rkstep(uPrim,uDim,xDim,deltaT,deltaX)
    if(abs(uPrim(1,3) - dummy) < 1.0e-14) then
        print*,'OK'
    else    
        print*,'Time step failed'
    end if

    ! ********************** GRID TEST ********************** 

    print*,'Testing grid values....'

    dummy = 1.0/3.0
    call gridValues(xAxis,1.0,xDim,deltaX)
    print*,xAxis
    if(abs(deltaX - dummy) < 1.0e-14) then
        print*,'OK'
    else    
        print*,'Time step failed'
    end if

    ! ********************** SOD TEST ********************** 

    print*,'Testing Sod initial condition setup....'

    uPrim(:,:) = 0.0
    call sodIC(uPrim,uDim,xDim,1)
    call eos(uPrim,5.0/3.0,uDim,xDim)

    dummy = 0.1
    if(abs(uPrim(6,xDim) - dummy) < 1.0e-14) then
        print*,'OK'
    else    
        print*,'Sod initial condition failed'
    end if

   




end program test