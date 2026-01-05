
  subroutine idata_chargeddiracstar

! ***************************************************************
! ***   CHARGED DIRAC STAR INITIAL DATA IN POLAR AREAL GAUGE  ***
! ***************************************************************

! This subroutine calculates initial data for a charged
! Dirac star using a shooting method in the polar areal gauge.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays
  use derivatives

! Extra variables.

  implicit none

  integer i,l,iter                      ! Counters.

  real(8) half,smallpi                  ! Numbers.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                        ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g,alpha_g               ! Radial metric and lapse global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: F_g,G_g                   ! Global arrays for Dirac functions (f,g).
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: maxwellF_g,maxwellPhi_g   ! Global arrays for Maxwell fields (F,Phi)
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: maxwellE_g                ! Global array for Maxwell field (E)


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  smallpi = acos(-1.d0)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for a charged Dirac star in polar-areal gauge using shooting method'
     print *, 'This is not yet implemented ...'
     print *
  end if

! Sanity check.

  if (spacetime=="minkowski") then
     print *, 'Charged Dirac star initial data is not compatible with a Minkowski background ...'
     print *, 'Aborting! (subroutine idata_chargeddiracstar)'
     print *
     call die
  end if


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

  do l=0,Nl-1
     do i=1-ghost,Nrtotal
        rr(l,i) = (dble(Nmin(rank) + i) - half)*dr(l)
     end do
  end do


! ***********************
! ***   INTEGRATION   ***
! ***********************

! Here we use a shooting method to get the correct
! asymptotic behaviour.  On each iteration, the
! integration is done with fourth order Runge-Kutta.

! For the moment only processor 0 solves for the
! initial data.  The solution is later distributed
! among processors.

! Initialize arrays.

  A_g     = 1.d0
  alpha_g = 1.d0

  F_g = 0.d0
  G_g = 0.d0

  maxwellF_g = 0.d0
  maxwellE_g = 0.d0
  maxwellPhi_g = 0.d0


! ***************
! ***   END   ***
! ***************

  end subroutine idata_chargeddiracstar
