!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/fluidintegral.f90,v 1.7 2025/10/01 17:52:44 malcubi Exp $

  subroutine fluidintegral

! ************************************************
! ***   CALCULATION OF INTEGRATED FLUID MASS   ***
! ************************************************

! Integrated fluid rest mass:
!
!                    /                      /             1/2      6   2
! fluid_restmass  =  | fluid_cD dV  =  4 pi | fluid_cD [ A   B  psi ] r dr
!                    /                      /
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

  logical contains

  integer i,l,i0,p
  integer status(MPI_STATUS_SIZE)

  real(8) r0,interp
  real(8) delta,m0,m1
  real(8) NTOT
  real(8) aux,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  smallpi = acos(-1.d0)


! ********************************
! ***   INTEGRATE FLUID MASS   ***
! ********************************

! The total rest mass is the integral of fluid_cD (not fluid_rho).
! Notice that fluid_cD and fluid_rho coincide for a fluid at rest,
! but if the fluid is moving then fluid_rho is the rest mass density
! in the fluid's rest frame, while fluid_cD is the rest mass denisty
! as measured by the Eulerian observers.  The conserved quantity is
! then the integral of fluid_cD.

  if (contains(mattertype,"fluid")) then

     fluid_restmass = 0.d0

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
        fluid_restmass(l,:) = integral(l)
     end do

!    Restrict integral.

     if (Nl>1) then
        intvar => fluid_restmass
        call restrictintegral
     end if


!    ****************************************
!    ***   OUTPUT TOTAL FLUID REST MASS   ***
!    ****************************************

!    Only at t=0.

     if (t(0)==0.d0) then

        NTOT = fluid_restmass(0,Nr)

        if (size==1) then
           write(*,'(A,E19.12)') ' Total fluid rest mass M0 = ',NTOT
        else
           if (rank==0) then
              p = size-1
              call MPI_RECV(NTOT,1,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              write(*,'(A,E19.12)') ' Total fluid rest mass = ',NTOT
           else if (rank==size-1) then
              call MPI_SEND(NTOT,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           end if
        end if

     end if

  end if


! *************************************
! ***   INTEGRATE DUST REST MASS   ***
! *************************************

! The total rest mass is the integral of dust_cD.
! (see comment above for fluid_mass).

  if (contains(mattertype,"dust")) then

     dust_restmass = 0.d0

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
        dust_restmass(l,:) = integral(l)
     end do

!    Restrict integral.

     if (Nl>1) then
        intvar => dust_restmass
        call restrictintegral
     end if


!    ***************************************
!    ***   OUTPUT TOTAL DUST REST MASS   ***
!    ***************************************

!    Only at t=0.

     if (t(0)==0.d0) then

        NTOT = dust_restmass(0,Nr)

!       Processor 0 writes result to screen.

        if (size==1) then
           write(*,'(A,E19.12)') ' Total dust rest mass = ',NTOT
        else
           if (rank==0) then
              p = size-1
              call MPI_RECV(NTOT,1,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              write(*,'(A,E19.12)') ' Total dust rest mass = ',NTOT
           else if (rank==size-1) then
              call MPI_SEND(NTOT,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           end if
        end if

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine fluidintegral


