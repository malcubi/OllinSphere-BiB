
  subroutine die

  use mpi
  use procinfo

! Subroutine for a clean code abort.

! ************************
! ***   FINALIZE MPI   ***
! ************************

! When compilig without MPI this call does nothing at all.

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)


! ****************
! ***   STOP   ***
! ****************

  stop


! ***************
! ***   END   ***
! ***************

  end subroutine die

