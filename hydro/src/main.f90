program hydro
    use primatives, only:prim_calc
    use conservedVars, only:cons_calc
    use boundary, only:boundaries
    use equationOfState, only:eos
    implicit none
    print*,'hello world!'
#ifdef PPT
    print*,'this tests the preprocessor.'
#endif
end program hydro