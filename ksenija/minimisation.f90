 subroutine pzndrv2
    use DEBUG_MODULE
    use STAT_MODULE
    IMPLICIT NONE
    include 'mpif.h'
    
! 
!     %---------------%
!     | MPI INTERFACE |
!     %---------------%
 
      INTEGER ::comm, myid, nprocs, rc, nloc
!     %-----------------------------%
!     | Define maximum dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      INTEGER ::maxn, maxnev, maxncv, ldv, nsites, nzmax
      parameter (maxn=10000, maxnev=100, maxncv=300, ldv=maxn)

!  MY VARIABLES FOR THE DISORDERED MATRIX
   logical, parameter :: disorder = .false., haldane = .true.
   logical :: continue_flag
   complex*16,allocatable::cooa(:), cooir(:), coojc(:)
   complex*16,allocatable:: AMATRIX(:)
   integer,allocatable :: IA(:),JA(:)
   complex*16,allocatable:: v_global(:), w_global(:)         ! SHARED VECTOR FOR MULTIPLICATION
   
!  MY VARIABLES FOR CHECK AGAINST LAPACK
   Complex*16, allocatable :: dense_check(:,:), dns(:,:)


!     %--------------%
!     | Local Arrays |
!     %--------------%

      INTEGER ::iparam(11), ipntr(14)
      logical::select(maxncv)
      Complex*16 ::ax(maxn), d(maxncv),             &
                        &  v(ldv,maxncv), workd(3*maxn),    &
                        &  workev(3*maxncv), resid(maxn),   &
                        &  workl(3*maxncv*maxncv+5*maxncv)
      REAL*8::rwork(maxncv), rd(maxncv,3)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      INTEGER ::ido, n, nx, nev, ncv, lworkl, info, i, j
      INTEGER ::ierr, nconv, maxitr, ishfts, mode, nmin, nmax
      Complex*16 ::sigma
      REAL*8 ::tol, W
      real*8 :: t_pzaupd, t_gather, t_amux, time1, time2, time3, time4
      logical ::rvec
      integer :: run_number
      character(len=20) :: evecfilename
      character(len=20) :: input_file
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      REAL*8 ::pdznorm2
      external          pdznorm2, zaxpy 
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%

      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )

!
!     %--------------------------------------------------%
!     | The number NX is the number of interior points   |
!     | in the discretization of the 2-dimensional       |
!     | convection-diffusion operator on the unit        |
!     | square with zero Dirichlet boundary condition.   | 
!     | The number N(=NX*NX) is the dimension of the     |
!     | matrix.  A standard eigenvalue problem is        |
!     | solved (BMAT = 'I').  NEV is the number of       |
!     | eigenvalues to be approximated.  The user can    |
!     | modify NX, NEV, NCV, WHICH to solve problems of  |
!     | different sizes, and to get different parts of   |
!     | the spectrum.  However, The following            |
!     | conditions must be satisfied:                    |
!     |                   N <= MAXN                      |
!     |                 NEV <= MAXNEV                    |
!     |           NEV + 2 <= NCV <= MAXNCV               | 
!     %--------------------------------------------------% 
!

      ! Get the input file from the command line
      call get_command_argument(1, input_file)
      

      open(100, file=trim(input_file), status="old")
      read(100, *) nx, run_number
      close(100)

      W = 4.0d0

      if (haldane) then 
         n = 2*nx*nx
      else
         n = nx*nx*nx
      endif
     
      nev   = 4
      ncv   = 20 

! global v vector should have the dimensions of the matrix
      allocate(v_global(n), w_global(n))
!      allocate(dns(2*n,2*n))

!
!     %--------------------------------------%
!     | Set up distribution of data to nodes |
!     %--------------------------------------%
!
      nloc = (n / nprocs)
      if ( mod(n, nprocs) .gt. myid ) nloc = nloc + 1
     
      write(*,*)'myid:',myid,'nloc:',nloc
      
      if ( nloc .gt. maxn ) then
         print *, ' ERROR with _NDRV1: NLOC is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'SM'
!
!     %---------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD as         | 
!     | workspace.  Its dimension LWORKL is set as        |
!     | illustrated below.  The parameter TOL determines  |
!     | the stopping criterion. If TOL<=0, machine        |
!     | precision is used.  The variable IDO is used for  |
!     | reverse communication, and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is  |
!     | generated to start the ARNOLDI iteration.         | 
!     %---------------------------------------------------%
!
      lworkl  = 3*ncv**2+5*ncv 
      tol    = 0.001 
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shift with respect to     |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of ZNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 10000000
      mode   = 1

      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode 
      
      if (disorder) then

        ! create the disordered matrix 
        nzmax=6*n-6*nx**2 + n

        allocate(AMATRIX(nzmax),IA(nzmax),JA(nzmax))
        allocate(cooa(nzmax),cooir(nzmax),coojc(nzmax))

        call sparse_matrix(cooa,cooir,coojc,nx)
        call coocsr(n,nzmax,cooa,cooir,coojc,amatrix,ja,ia)

      else if (haldane) then

        ! create the disordered matrix 
        nzmax= (9*(nx-2)**2 + 26*(nx-2) + 18)*2 + 2*(nx**2)

        allocate(AMATRIX(nzmax),IA(nzmax),JA(nzmax))
        allocate(cooa(nzmax),cooir(nzmax),coojc(nzmax))

        call haldane_sparse(cooa,cooir,coojc,nx,W)
        call coocsr(n,nzmax,cooa,cooir,coojc,amatrix,ja,ia)

      else

         allocate(dense_check(n,n), AMATRIX(3*n), IA(3*n), JA(3*n))
         ! create a matrix and convert to csr format
         call generate_hamiltonian(n, dense_check)
         call dnscsr(n, n, 3*n, dense_check, n, AMATRIX, JA, IA, ierr)

      endif

!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
!
    print*, myid, "starting znAupd"
    t_pzaupd = 0d0
    t_gather = 0d0
    t_amux = 0d0
    do while(ido .ge. -1 .and. ido .le. 1)
!
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         time1 = MPI_WTIME()
         call pznaupd ( comm, ido, bmat, nloc, which,      &
     &        nev, tol, resid, ncv, v, ldv, iparam, ipntr, & 
     &        workd, workl, lworkl, rwork,info )
         time2 = MPI_WTIME()
         

!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- OP*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%

             CALL MPI_ALLGATHER(workd(ipntr(1):ipntr(1)+nloc-1),nloc,MPI_DOUBLE_COMPLEX,   &
             &     v_global,nloc,MPI_DOUBLE_COMPLEX,   &
             &     comm,IERR)
             call MPI_Barrier(comm,IERR) 
             time3 = MPI_WTIME()
             call amux(n, nloc, v_global, workd(ipntr(2):ipntr(2)+nloc-1), AMATRIX, JA, IA, myid)
             call MPI_Barrier(comm,IERR)
             time4 = MPI_WTIME()
             t_pzaupd = t_pzaupd + (time2 - time1)/60
             t_gather = t_gather + (time3 - time2)/60
             t_amux = t_amux + (time4 -time3)/60

      enddo
      if (myid .eq. 0) print*, 'time:', t_pzaupd, t_gather, t_amux
! 
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD  |
!        %--------------------------%
!
         if ( myid .eq. 0 ) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation of _naupd'
            print *, ' '
         endif

      else 
      print*, myid, "finished znAupd"

!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.

         call pzneupd (comm, rvec, 'A', select, d, v, ldv, sigma, &
     &        workev, bmat, nloc, which, nev, tol, resid, ncv,    &
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,        &
     &        rwork, ierr)
!
!       %----------------------------------------------%
!       | Eigenvalues are returned in the one          |
!       | dimensional array D.  The corresponding      |
!       | eigenvectors are returned in the first NCONV |
!       | (=IPARAM(5)) columns of the two dimensional  | 
!       | array V if requested.  Otherwise, an         |
!       | orthogonal basis for the invariant subspace  |
!       | corresponding to the eigenvalues in D is     |
!       | returned in V.                               |
!       %----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!          %------------------------------------%
!          | Error condition:                   |
!          | Check the documentation of ZNEUPD. |
!          %------------------------------------%
!
            if ( myid .eq. 0 ) then
                print *, ' '
                print *, ' Error with _neupd, info = ', ierr
                print *, ' Check the documentation of _neupd. '
                print *, ' '
            endif

         else
             print*, myid, "finished znEupd"

             nconv = iparam(5)
             write(evecfilename, '(I0, "_evecs_", I0, ".dat")') nx, run_number
             open(222, file=trim(evecfilename), status='new', action='write')
             do j=1, nconv
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |   ||  A*x - lambda*x ||   |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (iparam(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
                CALL MPI_ALLGATHER(v(1,j),nloc,MPI_DOUBLE_COMPLEX,   &
                &     v_global,nloc,MPI_DOUBLE_COMPLEX,   &
                &     comm, IERR)

                if (myid .eq. 0) then

                   !print*, j, v_global(1)
                   !print*, j, v_global(4)

                   do i = 1, size(v_global)
                     write(222, '(F20.16, 1X, F20.16)') real(v_global(i)), aimag(v_global(i))
                   end do
                   write(222, *)

                endif
                call amux(n, nloc, v_global, ax, AMATRIX, JA, IA, myid)
                call zaxpy(nloc, -d(j), v(1,j), 1, ax, 1)
                rd(j,1) = dble(d(j))
                rd(j,2) = dimag(d(j))
                rd(j,3) = pdznorm2(comm, nloc, ax, 1)
            enddo
            close(222)
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
          open(999, file="outputlog.txt", status='unknown') 
             call pdmout(comm, 999, nconv, 3, rd, maxncv, -6, 'Ritz values (Real, Imag) and direct residuals')
          end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
         if (myid .eq. 0)then
         if ( info .eq. 1) then
             WRITE(999, *)  ' '
             WRITE(999, *)  ' Maximum number of iterations reached.'
             WRITE(999, *)  ' '
         else if ( info .eq. 3) then
             WRITE(999, *) ' ' 
             WRITE(999, *) ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
             WRITE(999, *) ' '
         end if      

          WRITE(999, *) ' '
          WRITE(999, *) '_NDRV1'
          WRITE(999, *) '====== '
          WRITE(999, *) ' '
          WRITE(999, *) ' Size of the matrix is ', n
          WRITE(999, *) ' The number of processors is ', nprocs
          WRITE(999, *) ' The number of Ritz values requested is ', nev
          WRITE(999, *) ' The number of Arnoldi vectors generated (NCV) is ', ncv
          WRITE(999, *) ' What portion of the spectrum: ', which
          WRITE(999, *) ' The number of converged Ritz values is ', nconv 
          WRITE(999, *) ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
          WRITE(999, *) ' The number of OP*x is ', iparam(9)
          WRITE(999, *) ' The convergence criterion is ', tol
          WRITE(999, *) ' '

          close(999)
         endif
      end if
!
!     %----------------------------%
!     | Done with program pzndrv1. |
!     %----------------------------%
!
 9000 continue


    END SUBROUTINE pzndrv2
