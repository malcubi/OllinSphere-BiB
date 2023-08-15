!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_complexpulse_charged.f90,v 1.3 2022/05/13 17:21:52 malcubi Exp $

  subroutine idata_complexpulse_charged

! **********************************************
! ***   CHARGED COMPLEX PULSE INITIAL DATA   ***
! **********************************************

! Initial data for a charged complex scalar field pulse
! with an initial non-zero charge density.
!
! NOTE: You can in fact also have initial data with zero
! initial charge density and non-zero initial current
! density, in which case you don't need this routine (that
! case is handeled by the routine idata_complexpulse.f90).
!
! This initial data leaves the lapse, shift, conformal
! spatial metric and extrinsic curvature as in Minkowski
! (which have already been set up in "initial".)
!
! We then give a simple initial profile in the scalar field
! and solve the Hamiltonian constraint for the conformal
! factor, coupled to the Gauss constraint for the electric
! field.
!
! This initial data is described in Section 3.2 of the paper:
! "Gravitational collapse of charged scalar field", J.M. Torres
! and  M. Alcubierre, Gen.Rel.Grav. 46:1773 (2014).

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays
  use derivatives
  use radialfunctions

! Extra variables.

  implicit none


! *************************************************
! ***   SET UP GAUSSIAN PULSE IN SCALAR FIELD   ***
! *************************************************

! Initial pulse.  Remember that the complex scalar field must be even.

! *******************
! ***   NUMBERS   ***
! *******************


! ***************
! ***   END   ***
! ***************

  if (rank==0) then
     print *, 'Charged complex pulse initial data not yet implemented ...'
     print *, 'Aborting! (subroutine idata_complexpulse_charged)'
     print *
     call die
  end if

  end subroutine idata_complexpulse_charged

