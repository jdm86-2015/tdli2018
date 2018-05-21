program arrayTest
    implicit none    
    integer, parameter :: bins = 10
    real, dimension(bins) :: A
    real, dimension(bins) :: B
    real, dimension(bins,bins) :: C
    integer :: i
    integer :: J

    A = (/ (J, J=1, 10) /)
    do i=1,10
        print*,'A(',i,'): ',A(i)
    end do
    B = cshift(A,-1)
    ! B = A/A
    print*,'*****************************'
    do i=1,10
        print*,'B(',i,'): ',B(i)
    end do

    ! do i = 1, bins
        ! do j = 1, bins
            ! C(i,j) = 10*(i-1) + j
        ! end do
    ! end do

    ! do i=1,bins
    !     do j = 1,bins
    !         print*,'C(',i,',',j,'): ',C(i,j)
    !     end do
    ! end do

    ! C(10,:) = C(10,:)*C(10,:)

    ! do i = 1, bins
        ! print*,'C(10,',i,'): ',C(10,i)
    ! end do

end program arrayTest