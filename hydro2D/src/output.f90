module output
    implicit none
    contains
        subroutine write(uPrim,xAxis,uDim,xDim,stepNumber,time)
            integer, intent(in) :: xDim
            integer, intent(in) :: uDim
            integer, intent(in) :: stepNumber

            real, intent(in), dimension(1:xDim) :: xAxis
            real, intent(in), dimension(uDim,-1:xDim+2) :: uPrim

            real, intent(in) :: time
  
            integer :: i
            character(len=20) :: filename
  
            write(filename, "('output/out',I6.6,'.dat')") stepNumber
  
            open(80, file=filename, status='new', form = 'formatted')
  
            do i = 1,xDim
                write(80,*) time, xAxis(i), uPrim(:,i)
            end do
  
            close(80)
        end subroutine write
end module output