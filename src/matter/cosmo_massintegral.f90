!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/cosmo_massintegral.f90,v 1.9 2025/09/24 17:25:02 malcubi Exp $

  subroutine cosmo_massintegral

! ***************************************************************
! ***   CALCULATION OF INTEGRATED MASS AND DENSITY CONTRAST   ***
! ***************************************************************

! This routine calculates the integral of perturbed density
! and density contrast for cosmological spacetimes.
!
! The calculation is identical to that of the routine massintegral.f90,
! except for the fact that we subtract the integral of the background
! density.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use integrals

! Extra variables.

  implicit none

  integer i,l,i0

  real(8) r0,interp
  real(8) delta,m0,m1
  real(8) INTOT
  real(8) half,one,two,third,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  one = 1.d0
  two = 2.d0
  third = 1.d0/3.d0

  smallpi = acos(-1.d0)


! **********************
! ***   INITIALIZE   ***
! **********************

  rho_contrast_int = 0.d0


! ********************************************
! ***   INTEGRAL OF DENSITY PERTURBATION   ***
! ********************************************

! This is the integral of the denisty perturbation
! with respect to the areal radius:
!
!            /             2
! m  =  4 pi | rho_pert r_a  dr_a
!            /

! Since the code does not work in areal coordinates,
! we have in general the following transformation:
!
! r_a  =  r psi**2 sqrt(B)
!
! with "r" the standard radial coordinate for the code. This implies:
!
!    2          2    6  3/2
! r_a dr_a  =  r  psi  B   ( 1 + r d B / (2 B)  +  2 r d psi / psi ) dr
!                                   r                   r

  auxarray = 4.d0*smallpi*rho_pert*r**2*psi**6*B**1.5d0*(one + r*(half*D1_B/B + two*D1_phi))

  intvar => auxarray

  do l=0,Nl-1
     rho_pert_int(l,:) = integral(l)
  end do

! Restrict integral.

  if (Nl>1) then
     intvar => rho_pert_int
     call restrictintegral
  end if


! ****************************************
! ***   INTEGRAL OF DENSITY CONTRAST   ***
! ****************************************

! The integral of the density contrast is just the
! integral of the density perturbation divided by
! the background density.

  do l=0,Nl-1
     rho_contrast_int(l,:) = rho_pert_int(l,:)/cosmobg_rho(l)
  end do


! ***********************
! ***   COMPACTNESS   ***
! ***********************

! The compactness of the perturbation is defined as:
!
! C(r)  :=  [ M(r) - M_b(r) ] / r_area  =  rho_pert_int(r) / r_area
!
! where M(r) is the integrated total mass, and M_b(r)
! is the integrated background mass.

  cosmo_compactness = rho_pert_int/r_area


! *********************************
! ***   OUTPUT TOTAL INTEGRAL   ***
! *********************************

! Only at t=0.

  if (s(0)==0) then

     INTOT = rho_pert_int(0,Nr)

!    Multiple processor run. Broadcast the boson number from
!    processor size-1 (the boundary) to all other processors.

     if (size>1) then
        call MPI_BCAST(INTOT,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)
     end if

!    Processor 0 writes result to screen.

     if (rank==0) then
        print *
        write(*,'(A,E15.8)') ' Total integrated density perturbation: ',INTOT
        write(*,'(A,E15.8)') ' Total integrated density contrast:     ',INTOT/cosmobg_rho(0)
        print *
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine cosmo_massintegral

