!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/procinfo.f90,v 1.5 2018/09/04 15:57:23 malcubi Exp $

! *********************************
! ***   PROCESSOR INFORMATION   ***
! *********************************

! size              Number of processors.
! rank              Local processor number.
! ierr              MPI error code.
! errorcode         MPI error code.
!
! Nrmax:            Maximum local number of grid points in r direction.
!
! Nrl:              Array with local number of grid points in r direction.
! Nmin:             Local first grid point.
! Nmax:             Local last grid point.

  module procinfo

  integer :: rank
  integer :: size
  integer :: ierr,errorcode
  integer :: Nrmax

! These arrays have one index that identifies
! the processor number for MPI runs.

  integer, allocatable, dimension(:) :: Nrl,Nmin,Nmax

! These arrays have two indices, one that identifies
! the processor number, and the second one the grid
! refinement level.

  real(8), allocatable, dimension(:,:) :: rleft,rright

  end module procinfo

