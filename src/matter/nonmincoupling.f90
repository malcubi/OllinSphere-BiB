!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/nonmincoupling.f90,v 1.2 2021/05/07 00:18:11 malcubi Exp $

  subroutine nonmincoupling

! ************************************************
! ***   NON-MINIMAL COUPLING FUNCTION f(phi)   ***
! ************************************************

! This routine calculates the non-minimal coupling
! function f(phi) and its derivatives.
!
! Notice that the normalization is such that the
! effective gravitational constant is Geff=1/(8 pi f),
! so that f=1/(8 pi) corresponds to G=1 (the standard
! choice of units taken in the code).

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  real(8) smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  smallpi = acos(-1.d0)


! *******************
! ***   TRIVIAL   ***
! *******************

! By default set f = 1/(8 pi).

  nonmin_f   = 0.125d0/smallpi
  nonmin_fp  = 0.d0
  nonmin_fpp = 0.d0


! *********************
! ***   QUADRATIC   ***
! *********************

! Quadratic coupling:
!
! nonmin_f = 1/(8 pi) + xi*(phi-phi0)^2

  if (nonminf=="quadratic") then
     nonmin_f   = 0.125d0/smallpi + nonmin_fxi*(nonmin_phi-nonmin_phi0)**2
     nonmin_fp  = 2.d0*nonmin_fxi*(nonmin_phi-nonmin_phi0)
     nonmin_fpp = 2.d0*nonmin_fxi
  end if


! ***********************
! ***   BRANS-DICKE   ***
! ***********************

  if (nonminf=="BransDicke") then

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine nonmincoupling

