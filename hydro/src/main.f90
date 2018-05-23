program hydro
    use hydroDriverRoutine, only:hydroDriver
    implicit none

    real, parameter :: xLength = 1.0
    integer, parameter :: xDim = 256
    real, parameter :: tMax = 0.01

    call hydroDriver(xLength, xDim, tMax)
end program hydro