!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_diracpulse.f90,v 1.1 2023/08/17 20:18:06 malcubi Exp $

  subroutine idata_diracpulse

! ************************************
! ***   DIRAC PULSE INITIAL DATA   ***
! ************************************

! This initial data leaves the lapse, shift, conformal
! spatial metric and extrinsic curvature as in Minkowski
! (which have already been set up in "initial").
!
! We then give a simple initial profile in the Dirac field
! and solve the Hamiltonian constraint for the conformal
! factor.
!
! Notice that the initial pulse is assumed to be time-symmetric,
! as otherwise we would have a non-zero energy flux and we would
! need to solve the coupled momentum and Hamiltonian constraints.
!
! The Hamiltonian constraint that we need to solve has the form:
!
!  __2                5
!  \/ psi  +  2 pi psi rho  =  0
!
!
! with rho the energy density of the Dirac field:
!
!
! rho  =  1/(2 pi) [ m ( FR**2 + FI**2 - GR**2 - GI**2 )
!
!      +  2/(r psi**2 sqrt(B)) ( FR*GI - FI*GR )
!
!      +  1/(psi**2 sqrt(A)) ( FR*dGI/dr - FI*dGR/dr + GR*dFI/dr - DI*dFR/dr ) ]
!
!
! (assuming that the initial pulse is time-symmetric), and
! where the Laplacian of psi is given in general by:
!
! __2          2
! \/ psi  = [ d psi  +  (2/r - d A/2A + d B/B) d psi ] / A
!              r                r        r      r
!
! In most cases we will have A=B=1, so the Laplacian simplifies.
! But this won't be the case for a transformed radial coordinate.
! Notice also that due to the form of the energy density, the
! Hamiltonian constraint is always non-linear in psi.
!
! NOTE FOR PARALLEL RUNS:  The initial data is not really
! solved in parallel.  It is in fact solved only on processor
! zero on a full size array, and then it is distributed
! among the other processors.  This is slow, but works.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays
  use derivatives
  use radialfunctions

! Extra variables.

  implicit none

  integer i,l


! *******************
! ***   NUMBERS   ***
! *******************


! ***********************************************
! ***   SET UP INITIAL PULSE IN DIRAC FIELD   ***
! ***********************************************

! Initial pulse.  Remember that F must be even and G must be odd.

  if (diracGR_r0*diracGI_r0==0.d0) then
     print *, 'For "diracpulse" initial data diracGR_r0 and diracGI_r0 must both be non-zero ...'
     print *, 'Aborting! (subroutine idata_diracpulse)'
     print *
     call die
  end if

! Gaussian profile.

  if (diracprofile=="gaussian") then

!    FR and FI.

     if (diracFR_r0==0.d0) then
        dirac_FR = gaussian(diracFR_a0,0.d0,diracFR_s0)
     else
        dirac_FR = gaussian(diracFR_a0,+diracFR_r0,diracFR_s0) &
                 + gaussian(diracFR_a0,-diracFR_r0,diracFR_s0)
     end if

     if (diracFI_r0==0.d0) then
        dirac_FI = gaussian(diracFI_a0,0.d0,diracFI_s0)
     else
        dirac_FI = gaussian(diracFI_a0,+diracFI_r0,diracFI_s0) &
                 + gaussian(diracFI_a0,-diracFI_r0,diracFI_s0)
     end if

!    GR and GI.

     dirac_GR = gaussian(diracGR_a0,+diracGR_r0,diracGR_s0) &
              - gaussian(diracGR_a0,-diracGR_r0,diracGR_s0)

     dirac_GI = gaussian(diracGI_a0,+diracGI_r0,diracGI_s0) &
              - gaussian(diracGI_a0,-diracGI_r0,diracGI_s0)

! Smooth top-hat profile.

  else if (diracprofile=="tophat") then

!   FR and FI.

     if (diracFR_r0==0.d0) then
        dirac_FR = tophat(diracFR_a0,0.d0,diracFR_s0,diracFR_t0)
     else
        dirac_FR = tophat(diracFR_a0,+diracFR_r0,diracFR_s0,diracFR_t0) &
                 + tophat(diracFR_a0,-diracFR_r0,diracFR_s0,diracFR_t0)
     end if

     if (diracFI_r0==0.d0) then
        dirac_FI = tophat(diracFI_a0,0.d0,diracFI_s0,diracFI_t0)
     else
        dirac_FI = tophat(diracFI_a0,+diracFI_r0,diracFI_s0,diracFI_t0) &
                 + tophat(diracFI_a0,-diracFI_r0,diracFI_s0,diracFI_t0)
     end if

!   GR and GI.

    dirac_GR = tophat(diracGR_a0,+diracGR_r0,diracGR_s0,diracGR_t0) &
             - tophat(diracGR_a0,-diracGR_r0,diracGR_s0,diracGR_t0)

    dirac_GI = tophat(diracGI_a0,+diracGI_r0,diracGI_s0,diracGI_t0) &
             - tophat(diracGI_a0,-diracGI_r0,diracGI_s0,diracGI_t0)

  end if

! Spatial derivatives with finite differences.

  diffvar => dirac_FR
  do l=0,Nl-1
     D1_dirac_FR(l,:) = diff1(l,+1)
  end do

  diffvar => dirac_FI
  do l=0,Nl-1
     D1_dirac_FI(l,:) = diff1(l,+1)
  end do

  diffvar => dirac_GR
  do l=0,Nl-1
     D1_dirac_GR(l,:) = diff1(l,-1)
  end do

  diffvar => dirac_GI
  do l=0,Nl-1
     D1_dirac_GI(l,:) = diff1(l,-1)
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine idata_diracpulse
