module hydroDriverRoutine
    use equationOfState, only:eos
    use setupGrid, only:gridValues
    use setSod, only:sodIC
    use rkStepper, only:rkStep
    use timeStepCalc, only:timeStep
    use output, only:write
    implicit none
    contains
        subroutine hydroDriver(xLength, xDim, tMax)
            integer, intent(in) :: xDim
            integer, parameter :: uDim = 7
            integer :: stepNumber
            
            real, intent(in) :: tMax
            real, intent(in) :: xLength

            real, dimension(xDim) :: xAxis
            real, dimension(uDim,-1:xDim+2) :: uPrim

            real :: deltaX
            real :: deltaT
            real :: time
            real :: timeOut
            real :: safety
            real :: outFreq
            real :: gamma

            safety = 0.7
            outFreq = 0.1
            timeOut = time + outFreq
            stepNumber = 0
            gamma = 5.0/3.0

            time = 0.0

            ! get the x coordinate values and the spacing
            call gridValues(xAxis,xLength,xDim,deltaX)
            ! populate the initial condition
            call sodIC(uPrim,uDim,xDim,-1)
            ! populate the pressure and sound speed
            call eos(uPrim,gamma,uDim,xDim)
            ! determine the initial time step
            call timeStep(uPrim,uDim,xDim,deltaX,deltaT,safety)
            ! write the initial condition to disk
            call write(uPrim,xAxis,uDim,xDim,stepNumber,time)

            do while ((time < tMax) .OR. deltaT < 1e-10)

                stepNumber = stepNumber + 1
                call timeStep(uPrim,uDim,xDim,deltaX,deltaT,safety)
                call rkStep(uPrim,uDim,xDim,deltaT,deltaX)
                time = time + deltaT
                call write(uPrim,xAxis,uDim,xDim,stepNumber,time)
                if(time > timeOut) then
                    call write(uPrim,xAxis,uDim,xDim,stepNumber,time)
                    timeOut = time + outFreq
                end if
                if(deltaT < 1e-10) then
                    print*,"Underflow...."
                end if
                ! print*,time
            end do
        end subroutine hydroDriver
end module hydroDriverRoutine