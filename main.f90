program QRDecomposition
    implicit none

    ! Define Variables
    integer :: i, j, n, iter, o,nxn
    real(kind=8) :: a_red1, gam, beta, beta_sq
    real(kind=8), dimension(:), allocatable :: a_red, h_small, e
    real(kind=8), dimension(:,:), allocatable :: A, R, Q, H_sub, H, identity, eXe

    character(64):: ofmt !creating rounded formate for matrix

    ! Read matrix A from input file (assuming a text file with space-separated values)
    open(unit=10, file='input.inp', status='old', action='read')
    open(unit=20, file='output.log', status='replace', action='write')
    read(10, *) n    ! Reading matrix dimention

    ! Error Handling
    if (n /= 3) then
        write (*,*) 'Program restricted to 3x3 Matrices only'
        stop
    end if

    allocate(A(n,n))
    do i = 1, n
        read(10, *, iostat=nxn) (A(i,j), j=1,n)
        if (nxn < 0) then
           write(*,*) "Matrix not", n ,"x", n
           stop
        end if
    end do
    close(unit=10)

    write(ofmt, '(a,i2,a)') "(", n, "(1x,f8.3))"  ! defining the formate

    ! Allocating and Initialize Q and R
    allocate(Q(n,n))
    allocate(R(n,n))
    Q = 0.0
    do i = 1,n
        Q(i,i) = 1.0
    end do
    R = A

    ! Allocate arrays for Householder transformations
    allocate(a_red(n))
    allocate(h_small(n))
    allocate(e(n))
    allocate(identity(n,n))
    allocate(eXe(n,n))

    ! QR Decomposition with Householder transformations
    do iter = 1, n

        ! Step 1: Get the reduced column vector a_red for the current matrix
        a_red = R(iter:n, iter)
        a_red1 = a_red(1)   ! assigning value of a_red1 from first element of a_red vector

        ! Step 2: Determine coefficient gamma
        if (a_red1 >= 0.0) then
            gam = sqrt(sum(a_red**2))
        else
            gam = -sqrt(sum(a_red**2))
        end if

        ! Step 3: Compute vector h
        h_small = a_red
        h_small(1) = h_small(1) + gam

        ! Step 4: Determine coefficient beta
        beta_sq = 2.0 * (gam**2 + gam * a_red1)
        beta = sqrt(beta_sq)

        ! Step 5: Determine unit vector e and Householder submatrix H_sub
        e = 1.0
        e = h_small / beta

        ! Creating identity matrix
        identity = 0.0
        do i = 1, n
            identity(i,i) = 1.0
        end do

        ! creating e dyadic e matrix
        eXe = 0.0
        do i = 1,n-iter+1
            do j = 1, n-iter+1
                eXe(i,j) = e(i) * e(j)
            end do
        end do

        H_sub = 0.0
        H_sub = identity - 2.0 * eXe

        ! Step 6: Assemble Householder matrix H with identity matrix
        H = 0.0
        H = identity
        H(iter:n, iter:n) = H_sub

        ! Step 7: Update R via Ri+1 = H · Ri
        R = matmul(H, R)

        ! Step 8: Update Q via Qi+1 = Qi · H^T
        Q = matmul(Q, transpose(H))

        ! Print results for this iteration(iter)
        write (20,*) "Iteration: ", iter
        write (20,*) "a_red: "
        do o = 1, n-iter+1
            write (20,'(F8.3)') a_red(o)
        end do

        write (20,*) "gamma: ", gam

        write (20,*) "h_small: "
        do o = 1, n-iter+1
            write (20,'(F8.3)') h_small(o)
        end do

        write (20,*) "beta: ", beta
        write (20,*) "e: "
        do o = 1, n-iter+1
            write (20,'(F8.3)') e(o)
        end do

        write (20,*) "H_sub: "
        do o = 1, n-iter+1
            write (20,*) (H_sub(o,i), i=1,n-iter+1)
        end do

        write (20,*) "H: "
        do o = 1, n
            write (20,*) (H(o,i), i=1,n)
        end do

        write (20,*) "R: "
        do o = 1, n
            write (20,*) (R(o,i), i=1,n)
        end do

        write (20,*) "Q:"
        do o = 1, n
            write (20,*) (Q(o,i), i=1,n)
        end do
    end do

    ! Final result
    ! Matrix multiplication of Q and R for...
    !...comparing input values with our result values of A
    A = matmul(Q, R)

    write (20,*) "Final Q matrix:"
    do j = 1,n
        write(20,ofmt) (Q(j,i), i= 1,n)
    end do
    write (20,*) "Final R matrix:"
    do j = 1,n
        write(20,ofmt) (R(j,i), i= 1,n)
    end do

    write (20,*) "Final A matrix: Q . R"
    do j = 1,n
        write(20,ofmt) (A(j,i), i= 1,n)
    end do

    ! Deallocate memory
    deallocate(A)
    deallocate(Q)
    deallocate(R)
    deallocate(a_red)
    deallocate(h_small)
    deallocate(e)
    deallocate(identity)
    deallocate(eXe)

end program QRDecomposition
