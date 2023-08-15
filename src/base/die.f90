!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/die.f90,v 1.3 2021/11/04 18:44:20 malcubi Exp $

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

