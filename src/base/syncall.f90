
  subroutine syncall(l)

! *****************************************
! ***   SYNCHRONIZE ACROSS PROCESSORS   ***
! *****************************************

! This routine should only be called for parallel runs.

  use param

! Extra variables.

  implicit none

  integer l

! Synchronize geometric variables.

  call syncgeo(l)

! Synchronize matter variables.
  
  if (mattertype/="vacuum") then
     call syncmatt(l)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine syncall


