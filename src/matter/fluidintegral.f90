!$Header:$

  subroutine fluidintegral

! *****************************************************
! ***   CALCULATION OF INTEGRATED FLUID REST MASS   ***
! *****************************************************

! This is the integrated Dirac charge:
!
!              /                       /              1/2      6   2
! fluid_Nb  =  | fluid_rho dV  =  4 pi | fluid_rho [ A   B  psi ] r dr
!              /                       /
!
! Notice that this assumes that the spacetime is REGULAR at the
! origin, so it will not work for eternal black holes such as
! Schwarzschild or Reissner-Nordstrom.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use integrals

! Extra variables.

  implicit none

  integer i,l,i0

  real(8) delta
  real(8) NTOT
  real(8) aux,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  smallpi = acos(-1.d0)


! **********************
! ***   INITIALIZE   ***
! **********************

  fluid_Nb = 0.d0


! *********************
! ***   INTEGRATE   ***
! *********************

! Remember to subtract the atmosphere since we don't
! want toi add its contribution to the total rest mass.

  auxarray = 4.d0*smallpi*(fluid_rho-fluid_atmos)*dsqrt(A)*B*psi**6*r**2

  intvar => auxarray

  do l=0,Nl-1
     fluid_Nb(l,:) = integral(l)
  end do


! ********************
! ***   RESTRICT   ***
! ********************

! For the case of several grid levels we need to restrict from the
! fine to the coarse grids.

  restrictvar => fluid_Nb

  do l=Nl-1,1,-1

!    Substract constant difference at edge of fine grid. Since the
!    integrals are done independently at each grid level, we find that
!    after restriction there will be small jumps due to accumulated
!    numerical error when we pass from the end on a fine grid to the
!    next coarse grid. Here we correct for that.

     i0 = Nr/2

     if (size==1) then

!       Find difference between the value at the edge of the fine grid, 
!       and the value in the coarse grid.

        delta = fluid_Nb(l-1,i0) - 0.5d0*(fluid_Nb(l,2*i0) + fluid_Nb(l,2*i0-1))

     else

!       NOT YET IMPLEMENTED.

     end if

!    Restrict.

     call restrict(l,.false.)

!    Fix symmetries.

     if (rank==0) then
        do i=1,ghost
           fluid_Nb(l-1,1-i) = fluid_Nb(l-1,i)
        end do
     end if

  end do


! ****************************************
! ***   OUTPUT TOTAL FLUID REST MASS   ***
! ****************************************

! Only at t=0.

  if (t(0)==0.d0) then

     NTOT = fluid_Nb(0,Nr)

!    Processor 0 writes result to screen.

     if (rank==0) then
        write(*,'(A,E19.12)') ' Total baryon rest mass (Nb) = ',NTOT
        print *
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine fluidintegral


