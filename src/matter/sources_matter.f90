!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/sources_matter.f90,v 1.19 2023/08/17 20:19:23 malcubi Exp $

  subroutine sources_matter(l)

! ******************************
! ***   SOURCES FOR MATTER   ***
! ******************************

! This routine calculates the sources for
! the different matter types.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  logical contains

  integer l


! *****************************
! ***   REAL SCALAR FIELD   ***
! *****************************

  if (contains(mattertype,"scalar")) then
     call sources_scalar(l)
  end if


! ***********************
! ***   GHOST FIELD   ***
! ***********************

  if (contains(mattertype,"ghost")) then
     call sources_ghost(l)
  end if


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

  if (contains(mattertype,"complex")) then
     call sources_complex(l)
  end if


! *******************************
! ***   COMPLEX GHOST FIELD   ***
! *******************************

  if (contains(mattertype,"complexghost")) then
     call sources_complexghost(l)
  end if


! ************************
! ***   NONMIN FIELD   ***
! ************************

  if (contains(mattertype,"nonmin")) then
     call sources_nonmin(l)
  end if


! **************************
! ***   ELECTRIC FIELD   ***
! **************************

  if (contains(mattertype,"electric")) then
     call sources_electric(l)
  end if


! ***********************
! ***   PROCA FIELD   ***
! ***********************

  if (contains(mattertype,"proca")) then
     call sources_proca(l)
  end if


! *******************************
! ***   COMPLEX PROCA FIELD   ***
! *******************************

  if (contains(mattertype,"complexproca")) then
     call sources_complexproca(l)
  end if


! ***********************
! ***   DIRAC FIELD   ***
! ***********************

  if (contains(mattertype,"dirac")) then
     call sources_dirac(l)
  end if


! ****************
! ***   DUST   ***
! ****************

  if (contains(mattertype,"dust")) then
     call sources_dust(l)
  end if


! *************************
! ***   PERFECT FLUID   ***
! *************************

  if (contains(mattertype,"fluid")) then
     call sources_fluid(l)
  end if


! ******************************
! ***   KODAMA ENERGY FLUX   ***
! ******************************

  if (allocated(P_Kodama)) then
     call KEnergy_flux(l)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_matter

