!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/fluidintegral.f90,v 1.3 2025/03/04 19:38:25 malcubi Exp $

  subroutine fluidintegral

! *****************************************************
! ***   CALCULATION OF INTEGRATED FLUID REST MASS   ***
! *****************************************************

! This is the integrated Dirac charge:
!
!              /                      /             1/2      6   2
! fluid_Nb  =  | fluid_cD dV  =  4 pi | fluid_cD [ A   B  psi ] r dr
!              /                      /
!
! Notice that this assumes that the spacetime is REGULAR at the
! origin, so it will not work for eternal black holes such as
! Schwarzschild or Reissner-Nordstrom.
!
! This routine also does the same intregral for dust in case
! it is required.

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


! *************************************
! ***   INTEGRATE FLUID REST MASS   ***
! *************************************

! The total rest mass is the integral of fluid_cD (not fluid_rho!).

  if (allocated(fluid_Nb)) then

     fluid_Nb = 0.d0

!    Remember not to take into account the atmosphere.

     do l=0,Nl-1
        do i=1-ghost,Nr
           if (fluid_rho(l,i)>fluid_atmos) then 
              auxarray(l,i) = 4.d0*smallpi*fluid_cD(l,i)*dsqrt(A(l,i))*B(l,i)*psi(l,i)**6*r(l,i)**2
           else
              auxarray(l,i) = 0.d0
           end if
        end do
     end do

     intvar => auxarray

     do l=0,Nl-1
        fluid_Nb(l,:) = integral(l)
     end do


!    ********************
!    ***   RESTRICT   ***
!    ********************

!    For the case of several grid levels we need to restrict from the
!    fine to the coarse grids.

     restrictvar => fluid_Nb

     do l=Nl-1,1,-1

!       Substract constant difference at edge of fine grid. Since the
!       integrals are done independently at each grid level, we find that
!       after restriction there will be small jumps due to accumulated
!       numerical error when we pass from the end on a fine grid to the
!       next coarse grid. Here we correct for that.

        i0 = Nr/2

        if (size==1) then

!          Find difference between the value at the edge of the fine grid, 
!          and the value in the coarse grid.

           delta = fluid_Nb(l-1,i0) - 0.5d0*(fluid_Nb(l,2*i0) + fluid_Nb(l,2*i0-1))

        else

!          NOT YET IMPLEMENTED.

        end if

!       Restrict.

        call restrict(l,.false.)

!       Fix symmetries.

        if (rank==0) then
           do i=1,ghost
              fluid_Nb(l-1,1-i) = fluid_Nb(l-1,i)
           end do
        end if

     end do


!    ****************************************
!    ***   OUTPUT TOTAL FLUID REST MASS   ***
!    ****************************************

!    Only at t=0.

     if (t(0)==0.d0) then

        NTOT = fluid_Nb(0,Nr)

!       Processor 0 writes result to screen.

        if (rank==0) then
           write(*,'(A,E19.12)') ' Total fluid rest mass (Nb) = ',NTOT
        end if

     end if

  end if


! *************************************
! ***   INTEGRATE DUST REST MASS   ***
! *************************************

! The total rest mass is the integral of dust_cD (not dust_rho!).

  if (allocated(dust_Nb)) then

     dust_Nb = 0.d0

!    Remember not to take into account the atmosphere.

     do l=0,Nl-1
        do i=1-ghost,Nr
           if (dust_rho(l,i)>dust_atmos) then 
              auxarray(l,i) = 4.d0*smallpi*dust_cD(l,i)*dsqrt(A(l,i))*B(l,i)*psi(l,i)**6*r(l,i)**2
           else
              auxarray(l,i) = 0.d0
           end if
        end do
     end do

     intvar => auxarray

     do l=0,Nl-1
        dust_Nb(l,:) = integral(l)
     end do


!    ********************
!    ***   RESTRICT   ***
!    ********************

!    For the case of several grid levels we need to restrict from the
!    fine to the coarse grids.

     restrictvar => dust_Nb

!    NOT YET IMPLEMENTED.


!    ***************************************
!    ***   OUTPUT TOTAL DUST REST MASS   ***
!    ***************************************

!    Only at t=0.

     if (t(0)==0.d0) then

        NTOT = dust_Nb(0,Nr)

!       Processor 0 writes result to screen.

        if (rank==0) then
           write(*,'(A,E19.12)') ' Total dust rest mass (Nb) = ',NTOT
        end if

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine fluidintegral


