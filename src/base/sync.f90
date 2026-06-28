
  subroutine sync

! ********************************
! ***   SYNCHRONIZE VARIABLE   ***
! ********************************

! Load modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer :: left,right
  integer :: status(MPI_STATUS_SIZE)

  real(8) :: wsend(ghost),wrecv(ghost)


! *************************************
! ***   EXCHANGE LEFT GHOST ZONES   ***
! *************************************

  left  = rank-1

  if (left >= 0) then

     wsend = syncvar(1:ghost)

     call MPI_Sendrecv( &
        wsend, ghost, MPI_REAL8, left, 1, &
        wrecv, ghost, MPI_REAL8, left, 2, &
        MPI_COMM_WORLD, status, ierr)

     syncvar(1-ghost:0) = wrecv

  end if


! **************************************
! ***   EXCHANGE RIGHT GHOST ZONES   ***
! **************************************

  right = rank+1

  if (right < size) then

     wsend = syncvar(Nr-2*ghost+1:Nr-ghost)

     call MPI_Sendrecv( &
        wsend, ghost, MPI_REAL8, right, 2, &
        wrecv, ghost, MPI_REAL8, right, 1, &
        MPI_COMM_WORLD, status, ierr)

     syncvar(Nr-ghost+1:Nr) = wrecv

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sync

