program hydro
    use conservedVars, only:conVar
    implicit none
    print*,'hello world!'
#ifdef PPT
    print*,'this tests the preprocessor.'
#endif
end program hydro