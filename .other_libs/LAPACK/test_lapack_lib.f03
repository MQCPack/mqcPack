      Program Test_LAPACK
        implicit none
! This program is built in install.  It is not meant to be run 
! It was created to test that the user supplied a LAPACK lib that contains 
! something, adding a call to dgegv that does nothing.  This 
! is just a link time test to see that symbols are resolved. 

        integer i
        integer j
        integer k
        integer lrow
        integer info
        integer ivec(20)
        real    amat(46,20)
        real    vec(20)

        i = 20
        j = 4
        k = 5
        lrow = 2 * j + k + 1

        call dgbtrs ( 'n', i, j, k, 1, amat, lrow, ivec, vec, i, info )

end program Test_LAPACK

