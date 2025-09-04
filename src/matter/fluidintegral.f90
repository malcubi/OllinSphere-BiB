!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/fluidintegral.f90,v 1.4 2025/09/04 16:06:41 malcubi Exp $

  subroutine fluidintegral

! *****************************************************
! ***   CALCULATION OF INTEGRATED FLUID REST MASS   ***
! *****************************************************

! This is the integrated fluid rest mass:
!
!                /                      /             1/2      6   2
! fluid_mass  =  | fluid_cD dV  =  4 pi | fluid_cD [ A   B  psi ] r dr
!                /                      /
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

! The total rest mass is the integral of fluid_cD (not fluid_rho).
! Notice that fluid_cD and fluid_rho coincide for a fluid at rest,
! but if the fluid is moving then fluid_rho is the rest mass density
! in the fluid's rest frame, while fluid_cD is the rest mass denisty
! as measured by the Eulerian observers.  The conserved quantity is
! then the integral of fluid_cD.

  if (allocated(fluid_mass)) then

     fluid_mass = 0.d0

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
        fluid_mass(l,:) = integral(l)
     end do


!    ********************
!    ***   RESTRICT   ***
!    ********************

!    For the case of several grid levels we need to restrict from the
!    fine to the coarse grids.

     restrictvar => fluid_mass

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

           delta = fluid_mass(l-1,i0) - 0.5d0*(fluid_mass(l,2*i0) + fluid_mass(l,2*i0-1))

        else

!          NOT YET IMPLEMENTED.

        end if

!       Restrict.

        call restrict(l,.false.)

!       Fix symmetries.

        if (rank==0) then
           do i=1,ghost
              fluid_mass(l-1,1-i) = fluid_mass(l-1,i)
           end do
        end if

     end do


!    ****************************************
!    ***   OUTPUT TOTAL FLUID REST MASS   ***
!    ****************************************

!    Only at t=0.

     if (t(0)==0.d0) then

        NTOT = fluid_mass(0,Nr)

!       Processor 0 writes result to screen.

        if (rank==0) then
           write(*,*)
           write(*,'(A,E19.12)') ' Total fluid rest mass = ',NTOT
        end if

     end if

  end if


! *************************************
! ***   INTEGRATE DUST REST MASS   ***
! *************************************

! The total rest mass is the integral of dust_cD.
! (see comment above for fluid_mass).

  if (allocated(dust_mass)) then

     dust_mass = 0.d0

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
        dust_mass(l,:) = integral(l)
     end do


!    ********************
!    ***   RESTRICT   ***
!    ********************

!    For the case of several grid levels we need to restrict from the
!    fine to the coarse grids.

     restrictvar => dust_mass

!    NOT YET IMPLEMENTED.


!    ***************************************
!    ***   OUTPUT TOTAL DUST REST MASS   ***
!    ***************************************

!    Only at t=0.

     if (t(0)==0.d0) then

        NTOT = dust_mass(0,Nr)

!       Processor 0 writes result to screen.

        if (rank==0) then
           write (*,*)
           write(*,'(A,E19.12)') ' Total dust rest mass = ',NTOT
        end if

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine fluidintegral


