!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/diracintegral.f90,v 1.2 2024/09/23 18:50:17 malcubi Exp $

  subroutine diracintegral

! **************************************************
! ***   CALCULATION OF INTEGRATED DIRAC CHARGE   ***
! **************************************************

! This is the integrated Dirac charge:
!
!                /                        /               1/2      6   2
! dirac_Nint  =  | dirac_dens dV  =  4 pi | dirac_dens [ A   B  psi ] r dr
!                /                        /
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
  real(8) R99,NTOT
  real(8) aux,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  smallpi = acos(-1.d0)


! **********************
! ***   INITIALIZE   ***
! **********************

  dirac_Nint = 0.d0


! *********************
! ***   INTEGRATE   ***
! *********************

  auxarray = 4.d0*smallpi*dirac_dens*dsqrt(A)*B*psi**6*r**2

  intvar => auxarray

  do l=0,Nl-1
     dirac_Nint(l,:) = integral(l)
  end do


! ********************
! ***   RESTRICT   ***
! ********************

! For the case of several grid levels we need to restrict from the
! fine to the coarse grids.

  restrictvar => dirac_Nint

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

        delta = dirac_Nint(l-1,i0) - 0.5d0*(dirac_Nint(l,2*i0) + dirac_Nint(l,2*i0-1))

     else

!       NOT YET IMPLEMENTED.

     end if

!    Restrict.

     call restrict(l,.false.)

!    Fix symmetries.

     if (rank==0) then
        do i=1,ghost
           dirac_Nint(l-1,1-i) = dirac_Nint(l-1,i)
        end do
     end if

  end do


! *********************************************
! ***   OUTPUT TOTAL DIRAC CHARGE AND R99   ***
! *********************************************

! Calculate the radius for which we have 99% of the
! total charge. We do this only on the coarse grid,
! and use linear interpolation.
!
! Notice that at the moment we only do this calculation at
! t=0 (for the initial data).  I might change this later.

  if (t(0)==0.d0) then

     R99 = 0.d0
     NTOT = dirac_Nint(0,Nr)

!    Single processor run.

     if (size==1) then

        do i=1,Nr
           if ((abs(dirac_Nint(0,i-1))<0.99d0*abs(NTOT)).and.(abs(dirac_Nint(0,i))>=0.99d0*abs(NTOT))) then
              R99 = r(0,i-1) + dr(0)*(0.99d0*NTOT-dirac_Nint(0,i-1))/(dirac_Nint(0,i)-dirac_Nint(0,i-1))
           end if
        end do

!    Multiple processor run. Broadcast the boson number from
!    processor size-1 (the boundary) to all other processors.

     else

     end if

!    Processor 0 writes result to screen.

     if (rank==0) then
        write(*,'(A,E19.12)') ' Total Dirac charge    = ',NTOT
        write(*,'(A,E19.12)') ' R99 from Dirac charge = ',R99
        print *
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine diracintegral


