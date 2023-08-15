!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/sync.f90,v 1.6 2019/05/22 18:11:09 malcubi Exp $

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
  integer i0
  integer status(MPI_STATUS_SIZE)

  real(8), dimension(1:ghost) :: w     ! Auxiliary array for efficient communications.


! ***************************
! ***   SYNCHRONIZATION   ***
! ***************************

! Since we are using a pointer "syncvar" we have a problem with
! the lower and upper bounds.  If you have a pointer with
! identical dimensions to an array it copies its bounds, but
! if the pointer has lower dimensions the lower bound is
! always set to 1.  So here we need to find out the real
! lower bound of the pointer array and compensate.
 
  i0 = lbound(syncvar,dim=1) - 1 + ghost

! We only need to synchronize if we have
! more than one processor.

  if (size>1) then

!    Send information to processor on the left
!    (only if we are not on axis).

     if (rank/=0) then
        p = rank-1
        do i=1,ghost
           w(i) = syncvar(i+i0)
        end do
        call MPI_SEND(w,ghost,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
     end if

!    Send information to processor on the right
!    (only if we are not at outer boundary).

     if (rank<size-1) then
        p = rank+1
        do i=1,ghost
           w(i) = syncvar(Nr-2*ghost+i+i0)
        end do
        call MPI_SEND(w,ghost,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
     end if

!    Receive information from processor on the right
!    (only if we are not at outer boundary).

     if (rank<size-1) then
        p = rank+1
        call MPI_RECV(w,ghost,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        do i=1,ghost
           syncvar(Nr-ghost+i+i0) = w(i)
        end do
     end if

!    Receive information from processor on the left
!    (only if we are not on axis).

     if (rank/=0) then
        p = rank-1
        call MPI_RECV(w,ghost,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        do i=1,ghost
           syncvar(i-ghost+i0) = w(i)
        end do
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sync

