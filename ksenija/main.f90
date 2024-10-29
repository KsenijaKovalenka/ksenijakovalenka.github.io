
    PROGRAM myparpack
    IMPLICIT NONE
    include 'mpif.h'
    INTEGER ::ierr,rc,myid
    real :: start_time, end_time, execution_time

    call cpu_time(start_time)
    
    CALL MPI_INIT( ierr )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    
    CALL pzndrv2
    
    CALL MPI_FINALIZE(rc)

    call cpu_time(end_time)
    execution_time = end_time - start_time
    print *, 'Execution time: ', execution_time, ' seconds'
    
    END PROGRAM myparpack
