!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/fluidintegral.f90,v 1.5 2025/09/12 18:34:37 malcubi Exp $

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


!    ********************
!    ***   RESTRICT   ***
!    ********************

!    For the case of several grid levels we need to restrict from the
!    fine to the coarse grids.

     restrictvar => fluid_restmass

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

           delta = fluid_restmass(l-1,i0) - 0.5d0*(fluid_restmass(l,2*i0) + fluid_restmass(l,2*i0-1))

        else

!          Find values of radius (r0) and fluid_restmass (m0) at edge of fine grid,
!          and send them to all processors.  Notice that the edge of the grid
!          belongs to the processor with rank=size-1.

           r0 = 0.5d0*(r(l,2*i0)+r(l,2*i0-1))
           m0 = 0.5d0*(fluid_restmass(l,2*i0)+fluid_restmass(l,2*i0-1))

           call MPI_BCAST(r0,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(m0,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)

!          Now interpolate value of fluid_restmass at r0 from coarse grid,
!          and send it to all processors.  Notice that the point r0
!          is in the middle of the coarse grid, and belongs to the
!          processor with rank=(size-1)/2.

           interpvar => fluid_restmass
           m1 = interp(l-1,r0,.false.)
           call MPI_BCAST(m1,1,MPI_DOUBLE_PRECISION,(size-1)/2,MPI_COMM_WORLD,ierr)

!          Find difference bewteen m1 and m0.

           delta = m1 - m0

        end if

!       Subtract delta.

        fluid_restmass(l-1,:) = fluid_restmass(l-1,:) - delta

!       Restrict.

        call restrict(l,.false.)

!       Fix symmetries.

        if (rank==0) then
           do i=1,ghost
              fluid_restmass(l-1,1-i) = fluid_restmass(l-1,i)
           end do
        end if

     end do


!    ****************************************
!    ***   OUTPUT TOTAL FLUID REST MASS   ***
!    ****************************************

!    Only at t=0.

     if (t(0)==0.d0) then

        NTOT = fluid_restmass(0,Nr)

        if (size==1) then
           write(*,'(A,E19.12)') ' Total fluid rest mass = ',NTOT
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


!    ********************
!    ***   RESTRICT   ***
!    ********************

!    For the case of several grid levels we need to restrict from the
!    fine to the coarse grids.

     restrictvar => dust_restmass

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

           delta = dust_restmass(l-1,i0) - 0.5d0*(dust_restmass(l,2*i0) + dust_restmass(l,2*i0-1))

        else

!          Find values of radius (r0) and dust_restmass (m0) at edge of fine grid,
!          and send them to all processors.  Notice that the edge of the grid
!          belongs to the processor with rank=size-1.

           r0 = 0.5d0*(r(l,2*i0)+r(l,2*i0-1))
           m0 = 0.5d0*(dust_restmass(l,2*i0)+dust_restmass(l,2*i0-1))

           call MPI_BCAST(r0,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(m0,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)

!          Now interpolate value of dust_restmass at r0 from coarse grid,
!          and send it to all processors.  Notice that the point r0
!          is in the middle of the coarse grid, and belongs to the
!          processor with rank=(size-1)/2.

           interpvar => DUST_restmass
           m1 = interp(l-1,r0,.false.)
           call MPI_BCAST(m1,1,MPI_DOUBLE_PRECISION,(size-1)/2,MPI_COMM_WORLD,ierr)

!          Find difference bewteen m1 and m0.

           delta = m1 - m0

        end if

!       Subtract delta.

        dust_restmass(l-1,:) = dust_restmass(l-1,:) - delta

!       Restrict.

        call restrict(l,.false.)

!       Fix symmetries.

        if (rank==0) then
           do i=1,ghost
              dust_restmass(l-1,1-i) = dust_restmass(l-1,i)
           end do
        end if

     end do


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


