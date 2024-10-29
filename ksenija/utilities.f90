subroutine amux (n, nloc, x, y, a,ja,ia, myid) 
    !----------------------------------------------------------------------
    !PARALLEL Matrix Vector Multiplication (from BLASSM matvec.f)
    !----------------------------------------------------------------------
          complex*16  x(n), y(nloc), a(*) 
          integer n, ja(*), ia(*)
          complex*16 t
          integer i, k, l, nmin, nmax, myid
          nmin=1+myid*nloc
          nmax=(myid+1)*nloc
          l = 0
          do i=nmin,min(nmax,n)
                t = 0.0d0
                do k=ia(i), ia(i+1)-1 
                   t = t + a(k)*x(ja(k))
                enddo
                l = l + 1
                y(l) = t
          enddo
          return
    end subroutine amux

subroutine amux2 (n, x, y, a,ja,ia)
    !----------------------------------------------------------------------
    !Matrix Vector Multiplication (from BLASSM matvec.f)
    !----------------------------------------------------------------------
          complex*16  x(n), y(n), a(*)
          integer n, ja(*), ia(*)
          complex*16 t
          integer i, k
          do i=1,n
                t = 0.0d0
                do k=ia(i), ia(i+1)-1
                   t = t + a(k)*x(ja(k))
                enddo
                y(i) = t
          enddo
          return
    end subroutine amux2

    subroutine sparse_matrix(a, ir, jc, n)
        implicit none
        integer :: i, j, ix, iy, iz, orbital_index, seed(8), k     ! counters
        integer n                                  ! number of atomic sites along one direstion
        integer N_sites                            ! total number of atomic sites
        integer :: nnz, nzmax                                           ! for matrix construction in coo format
        complex*16, parameter :: t = (1.0d0, 0.0d0)    ! hopping (nn)
        real*8 :: eps, eps2                                  ! onsite energy
        real*8, parameter :: W = 16.5d0                ! disorder strength
!        real*8, parameter :: W = 0.000001d0                ! noise strength
        complex*16 :: a(*)       ! value arrays for coo and crs formats
        integer :: ir(*), jc(*) ! index arrays for coo and crs formats
    
        N_sites = n**3
        !---Estimate the number of non-zero elements
        nzmax = 6 * N_sites - 6 * n*n + N_sites 
        ! 6 neighbors per site
        ! minus the edges
        ! plus diagonal entries
    
        !---Construct the Hamiltonian matrix with nearest-neighbor hoppings
        nnz = 0  ! Initialize the number of non-zero entries
        i = 0
        DO iz = 1, n
            DO iy = 1, n
                DO ix = 1, n
                    !---count all the sites in the lattice
                    i = i + 1
                    ! Add nearest neighbors in 3D, avoiding going over the edge
                    IF (ix < n) THEN
                        ! in the x direction we just jump to the next index
                        j = i + 1
                        nnz = nnz + 1
                        a(nnz) = t
                        ir(nnz) = i
                        jc(nnz) = j
                        nnz = nnz + 1
                        a(nnz) = conjg(t)
                        ir(nnz) = j
                        jc(nnz) = i
                    END IF
                    IF (iy < n) THEN
                        ! in the y direction we just jump over a row
                        j = i + n
                        nnz = nnz + 1
                        a(nnz) = t
                        ir(nnz) = i
                        jc(nnz) = j
                        nnz = nnz + 1
                        a(nnz) = conjg(t)
                        ir(nnz) = j
                        jc(nnz) = i
                    END IF
                    IF (iz < n) THEN
                        ! in the z direction we just jump over a row
                        j = i + n*n
                        nnz = nnz + 1
                        a(nnz) = t 
                        ir(nnz) = i
                        jc(nnz) = j
                        nnz = nnz + 1
                        a(nnz) = conjg(t)
                        ir(nnz) = j
                        jc(nnz) = i
                    END IF
                END DO
            END DO
        END DO
    
    
        !---add random disorder
        data seed / 12345, 67890, 13579, 24680, 98765, 43210, 11111, 22222 /
        CALL RANDOM_SEED(PUT=seed) ! to have reproducible results, remove seed for new disorder realisation each time
        !open(333, file='onsites.dat', status='replace')
        i = 0
        DO ix = 1, N
            DO iy = 1, N
                DO iz = 1, N
                    i = i + 1
                    CALL RANDOM_NUMBER(eps)
                    nnz = nnz + 1
                    !write(333, '(F12.8, 1X)', advance="no") eps * W - W/2.0
                    !a(nnz) = eps/100000000
                    a(nnz) = eps * W - W/2.0
                    ! a(nnz) = cmplx(eps * W - W/2.0, eps * W - W/2.0)  !for random noise if has degenerate values
                    ! diagonal terms -> same index
                    ir(nnz) = i
                    jc(nnz) = i
                END DO
            END DO
        END DO
        !close(333)

    end subroutine sparse_matrix
    
    
    subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
    ! nrow dimension of the (dense) matrix; nnz - size of the coo arrays 
    ! nrow,nnz - coo; ao,jao,iao - csr
        complex*16 a(*),ao(*),x
        integer ir(*), jc(*), jao(*), iao(nrow+1), nrow, nnz, k0, iad, i, j, k
        
    ! Initialize the row pointer array to zero
        iao = 0
    
    ! determine row-lengths.
        do k=1, nnz
            iao(ir(k)) = iao(ir(k)) + 1
        enddo
    
    ! starting position of each row
        k = 1
        do j=1,nrow+1
                k0 = iao(j)
                iao(j) = k
                k = k+k0
        enddo
    
    ! go through the structure  once more. Fill in output matrix.
        do k=1, nnz
            i = ir(k)
            iad = iao(i)
            ao(iad) =  a(k)
            jao(iad) = jc(k)
            iao(i) = iad+1
        enddo
    
    ! shift back iao
        do j=nrow,1,-1
            iao(j+1) = iao(j)
        enddo
    
        iao(1) = 1
    end subroutine coocsr
     
subroutine dnscsr(nrow, ncol, nzmax, dns, ndns, a, ja, ia, ierr)
            complex*16 dns(ndns,*), a(*)
            integer ia(*), ja(*), nrow, ncol, ierr, nzmax, ndns, next, j, k
        !----------------------------------------------------------------------
        !Dense        to    Compressed Row Sparse (from FORMATS/formats.f)
        !---------------------------------------------------------------------- 
        !Note: this routine does not check whether an element 
        !is zero or not. i.e. you should only call it when you are 
        !sure that all the zero elements are not stored in dns. 
            ierr = 0
            next = 1
            ia(1) = 1
            do i=1, nrow
                do j=1, ncol
                    if (dns(i,j) .ne. 0.0d0) then
                        a(next) = dns(i,j)
                        a(next) = dns(i,j)
                        ja(next) = j
                        next = next + 1
                    endif
                end do
                ia(i+1) = next
            end do
            return
        end subroutine dnscsr

   subroutine generate_hamiltonian(n, H)
!----simple trigonal matrix to check against LAPACK
        integer, intent(in) :: n
        complex*16 , intent(out) :: H(n, n)
        integer :: i
        ! Initialize the matrix to zero
        H = (0.0d0, 0.0d0)

        ! Fill the diagonal with 2
        do i = 1, n
            H(i, i) = (2.0d0, 0.0d0)
        end do

        ! Fill the subdiagonal and superdiagonal with -1 - i
        do i = 1, n-1
            H(i, i+1) = (-1.0d0, 1.0d0)
            H(i+1, i) = (-1.0d0, -1.0d0)
        end do
    end subroutine generate_hamiltonian

    subroutine csrdns(nrow, ncol, a, ja, ia, dns, ndns, ierr)
        complex*16 dns(ndns,*), a(*)
        integer ja(*), ia(*), nrow, ncol, ndns, ierr, j, k, i
    !----------------------------------------------------------------------
    !Compressed Sparse Row    to    Dense  (from FORMATS/formats.f)
    !----------------------------------------------------------------------
        ierr = 0
        do i=1, nrow
            do j=1, ncol
                dns(i,j) = 0.0d0
            end do
        end do

        do i=1, nrow
            do k=ia(i), ia(i+1)-1
                j = ja(k)
                if (j .gt. ncol) then
                    ierr = i
                    return
                endif
                dns(i,j) = a(k)
            end do
        end do
        return

    end subroutine csrdns


subroutine haldane_sparse(a, ir, jc, n, W)
    implicit none
    integer :: i, j, l, m, ix, iy, seed(8)    ! counters
    integer n                                  ! number of atomic sites along one direstion
    integer :: nnz                                         ! for matrix construction in coo format
    complex*16, parameter :: t = (1.0d0, 0.0d0)    ! hopping (nn)
    real*8 :: eps, eps2                                  ! onsite energy
!    real*8, parameter :: W = 5.0d0                ! disorder strength
    real*8 :: W               ! disorder strength
    complex*16 :: a(*), t1, t2, phi       ! value arrays for coo and crs formats
    integer :: ir(*), jc(*)          ! index arrays for coo and crs formats
    real*8 :: delta
    real(16), parameter :: PI = 4 * atan (1.0d0)

    delta = 2.0d0
!    delta = 0.2d0
    phi = cmplx(0.0d0, PI / 2.0d0)
    t1 = -1.0d0
    t2 = real(t1) / 3.0d0 * exp(phi)

    !---Construct the Hamiltonian matrix with nearest-neighbor hoppings
    nnz = 0  ! Initialize the number of non-zero entries
    i = -1 ! to start with 1, 3, ...
    DO iy = 1, n
        DO ix = 1, n
            !---count all the unit cells (site A) in the lattice
            i = i + 2  ! site A
            ! Add nn hopping within the cell (nerby site B)
            j = i + 1  ! site B
            nnz = nnz + 1
            a(nnz) = t1
            ir(nnz) = i
            jc(nnz) = j
            nnz = nnz + 1
            a(nnz) = conjg(t1)
            ir(nnz) = j
            jc(nnz) = i

            ! hoppings to (1, 0) until we reach the edge
            IF (ix < n) THEN
                ! in the x direction we just jump to the next index 
                ! nn hopping
                l = i + 2   ! site A
                m = j + 2   ! site B
                nnz = nnz + 1
                a(nnz) = t1
                ir(nnz) = j
                jc(nnz) = l
                nnz = nnz + 1
                a(nnz) = conjg(t1)
                ir(nnz) = l
                jc(nnz) = j
                ! nnn hoppings
                nnz = nnz + 1
                a(nnz) = conjg(t2)
                ir(nnz) = i
                jc(nnz) = l
                nnz = nnz + 1
                a(nnz) = t2
                ir(nnz) = j
                jc(nnz) = m
                nnz = nnz + 1
                a(nnz) = t2
                ir(nnz) = l
                jc(nnz) = i
                nnz = nnz + 1
                a(nnz) = conjg(t2)
                ir(nnz) = m
                jc(nnz) = j
            END IF
            ! hoppings to (0, 1) until we reach the edge
            IF (iy < n) THEN
                ! in the y direction we just jump to the next row index
                ! nn hopping
                l = i + 2*n   ! site A
                m = j + 2*n   ! site B
                nnz = nnz + 1
                a(nnz) = t1
                ir(nnz) = j
                jc(nnz) = l
                nnz = nnz + 1
                a(nnz) = conjg(t1)
                ir(nnz) = l
                jc(nnz) = j
                ! nnn hoppings
                nnz = nnz + 1
                a(nnz) = t2
                ir(nnz) = i
                jc(nnz) = l
                nnz = nnz + 1
                a(nnz) = conjg(t2)
                ir(nnz) = j
                jc(nnz) = m
                nnz = nnz + 1
                a(nnz) = conjg(t2)
                ir(nnz) = l
                jc(nnz) = i
                nnz = nnz + 1
                a(nnz) = t2
                ir(nnz) = m
                jc(nnz) = j
                END IF
            ! hoppings to (1, -1) after we passe difrst row and until we reach the edge
            IF (ix < n .and. iy > 1) THEN
                l = i + 2 - 2*n   ! site A
                m = j + 2 - 2*n   ! site B
                ! nnn hoppings
                nnz = nnz + 1
                a(nnz) = t2
                ir(nnz) = i
                jc(nnz) = l
                nnz = nnz + 1
                a(nnz) = conjg(t2)
                ir(nnz) = j
                jc(nnz) = m
                nnz = nnz + 1
                a(nnz) = conjg(t2)
                ir(nnz) = l
                jc(nnz) = i
                nnz = nnz + 1
                a(nnz) = t2
                ir(nnz) = m
                jc(nnz) = j
            ENDIF
        END DO
    END DO

    !---add ONSITE ENERGIES
    !data seed / 12345, 67890, 13579, 24680, 98765, 43210, 11111, 22222 /
    CALL RANDOM_SEED() ! to have reproducible results, remove seed for new disorder realisation each time
    i = 0
    DO iy = 1, N
        DO ix = 1, N
            i = i + 1
            CALL RANDOM_NUMBER(eps)
            CALL RANDOM_NUMBER(eps2)
            nnz = nnz + 1
            a(nnz) = delta + (eps * W - W/2.0)
            ir(nnz) = i
            jc(nnz) = i
            i = i + 1
            nnz = nnz + 1
            a(nnz) = - delta + (eps2 * W - W/2.0)
            ir(nnz) = i
            jc(nnz) = i
        END DO
    END DO

end subroutine haldane_sparse
