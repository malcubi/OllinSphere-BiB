
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

  integer i,p
  integer status(MPI_STATUS_SIZE)

  real(8), dimension(1:ghost) :: w     ! Auxiliary array for efficient communications.


! ***************************
! ***   SYNCHRONIZATION   ***
! ***************************

! Send information to processor on the left
! (only if we are not on axis).

  if (rank/=0) then
     p = rank-1
     do i=1,ghost
        w(i) = syncvar(i)
     end do
     call MPI_SEND(w,ghost,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

! Send information to processor on the right
! (only if we are not at outer boundary).

  if (rank<size-1) then
     p = rank+1
     do i=1,ghost
        w(i) = syncvar(Nr-2*ghost+i)
     end do
     call MPI_SEND(w,ghost,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

! Receive information from processor on the right
! (only if we are not at outer boundary).

  if (rank<size-1) then
     p = rank+1
     call MPI_RECV(w,ghost,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=1,ghost
        syncvar(Nr-ghost+i) = w(i)
     end do
  end if

! Receive information from processor on the left
! (only if we are not on axis).

  if (rank/=0) then
     p = rank-1
     call MPI_RECV(w,ghost,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=1,ghost
        syncvar(i-ghost) = w(i)
     end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sync

