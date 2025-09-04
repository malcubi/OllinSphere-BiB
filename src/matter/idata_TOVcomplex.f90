!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_TOVcomplex.f90,v 1.1 2025/09/04 16:23:06 malcubi Exp $

  subroutine idata_TOVcomplex

! *******************************************
! ***   TOV STAR + COMPLEX SCALAR FIELD   ***
! ***        (FERMION-BOSON SYSTEM)       ***
! *******************************************

! This initial data is for a TOV star plus a complex
! scalar field pulse.  It is a fermion-boson system,
! but NOT a fermion-boosn star since the complex
! scalar field is not assumed to be stationary.
!
! The idea is to first solve for the TOV star and
! later add the scalar field and solve the 
! Hamiltonian constraint and lapse condition
! again.
!
! We the first call the idata_TOVstar routine.
! Once we come out from this we need to solve
! the Hamiltonian constraint for the conformal
! factor psi, which takes the form:
!
!  __2                5
!  \/ psi  +  2 pi psi rho  =  0
!
! where rho is the total energy denisty that includes
! the fluid and complex scalar field contributions:
!
! rho  =  rho_fluid + rho_complex
!
! with:
!
! rho_fluid    =  fluid_cE(l,i) + fluid_cD(l,i)
!
! and:
!                      2        2        2       2          4
! rho_complex  =  [ (pi_r  +  pi_i) + (xi_r  + xi_i) / A psi  ] / 2  +  V
!
! ! The Laplacian of psi is given in general by:
!
! __2           2
! \/ psi  =  [ d psi  +  (2/r - d A/2A + d B/B) d psi ] / A
!               r                r        r      r
!
! In most cases we will have B=1, but A will not be 1 since
! the TOV star initial data solved for A.
!
! Notice that the initial pulse is assumed to be time-symmetric,
! (i.e. pi=0), or else purely real with purely imaginary time
! derivative, as otherwise we would have a non-zero energy
! flux and we would need to solve the coupled momentum and
! hamiltonian constraints.
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

  integer :: l

  real(8) :: one,two,smallpi

! Local arrays.

  real(8), dimension (0:Nl-1,1-ghost:Nrmax)   :: rhoF,rhoS1,rhoS2
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: u,u_old

! u        Global solution of ODE.
! CL0,C0   Local and global coefficient of first derivative in ODE.
! CL1,C1   Local and global coefficient of linear term in ODE.
! CL2,C2   Local and global source term in ODE.


! *******************
! ***   NUMBERS   ***
! *******************

  one = 1.d0
  two = 2.d0

  smallpi = acos(-1.d0)


! **************************************
! ***   FIRST SOLVE FOR A TOV STAR   ***
! **************************************

  call idata_TOVstar


! ***********************************************
! ***  NOW ADD A COMPLEX SCALAR FIELD PULSE   ***
! ***********************************************

! Remember that the complex scalar field must be even.

! Gaussian profle.

  if (complexprofile=="gaussian") then

     if (complexR_r0==0.d0) then
        complex_phiR = gaussian(complexR_a0,0.d0,complexR_s0)
     else
        complex_phiR = gaussian(complexR_a0,+complexR_r0,complexR_s0) &
                     + gaussian(complexR_a0,-complexR_r0,complexR_s0)
     end if

     if (complexI_r0==0.d0) then
        complex_phiI = gaussian(complexI_a0,0.d0,complexI_s0)
     else
        complex_phiI = gaussian(complexI_a0,+complexI_r0,complexI_s0) &
                     + gaussian(complexI_a0,-complexI_r0,complexI_s0)
     end if

! Smooth top-hat profile.

  else if (complexprofile=="tophat") then

     if (complexR_r0==0.d0) then
        complex_phiR = tophat(complexR_a0,0.d0,complexR_s0,complexR_t0)
     else
        complex_phiR = tophat(complexR_a0,+complexR_r0,complexR_s0,complexR_t0) &
                     + tophat(complexR_a0,-complexR_r0,complexR_s0,complexR_t0)
     end if

     if (complexI_r0==0.d0) then
        complex_phiI = tophat(complexI_a0,0.d0,complexI_s0,complexI_t0)
     else
        complex_phiI = tophat(complexI_a0,+complexI_r0,complexI_s0,complexI_t0) &
                     + tophat(complexI_a0,-complexI_r0,complexI_s0,complexI_t0)
     end if

  end if

! Time derivatives.

  if (k_parameter==0.d0) then

!    The initial time derivatives are set to zero.

     complex_piR = 0.d0
     complex_piI = 0.d0

  else if ((k_parameter/=0.d0).and.(complex_mass>0.d0)) then

!    When k_parameter is not zero we set up an initially
!    purely real field.

     complex_phiI = 0.d0

!    The initial time derivative is harmonic with frequency k_parameter.
!    This guarantees that the initial data still satisfies the momentum
!    constraint trivially.

     complex_piR = 0.d0
     complex_piI = k_parameter*complex_phiR

  end if

! The spatial derivatives are calculated with finite differences.

  diffvar => complex_phiR

  do l=0,Nl-1
     complex_xiR(l,:) = diff1(l,+1)
  end do

  diffvar => complex_phiI

  do l=0,Nl-1
     complex_xiI(l,:) = diff1(l,+1)
  end do


! ******************************************
! ***   COMPLEX SCALAR FIELD POTENTIAL   ***
! ******************************************

! Find complex field potential.

  do l=0,Nl-1
     call potential(l)
  end do


! ********************************
! ***   TOTAL ENERGY DENSITY   ***
! ********************************

! Fluid contribution:
!
! rho_fluid  =  fluid_cE + flud_cD

  rhoF = fluid_cE + fluid_cD

! Complex scalar field contribution. We separate
! this in two contributions depending on the power
! of psi, so that we have the full scalar contribution
! given by:
!
! rho_scalar = rhoS1/psi**4 + rhoS2
!
! Notice that in the Hamiltonian constraint we need
! in fact to multiply this with A*psi^5, so that we
! will have:
!
! 2 pi A psi^5 rho  =  2 pi ( A rhoS2 psi^5 + rhoS1 psi )

  rhoS1 = 0.5d0*(complex_xiR**2 + complex_xiI**2)/A
  rhoS2 = 0.5d0*(complex_piR**2 + complex_piI**2) + complex_V


! ****************************************
! ***   SOLVE HAMILTONIAN CONSTRAINT   ***
! ****************************************


! ************************************************
! ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
! ************************************************

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! array with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  if (size==1) then
     !psi = u
  else
     call distribute(0,Nl-1,psi,u)
  end if

! Derivative of psi.

  diffvar => psi

  do l=0,Nl-1
     D1_psi(l,:) = diff1(l,+1)
  end do


! ***************************************
! ***   FIND CONFORMAL FUNCTION phi   ***
! ***************************************

! Find phi.

  phi  = dlog(psi)
  D1_phi = D1_psi/psi

  if (chimethod) then
     chi  = one/psi**chipower
     D1_chi = - dble(chipower)*D1_psi/psi**3
  end if

! Find psi2 and psi4.

  psi2 = psi**2
  psi4 = psi**4


! **********************
! ***   FIND LAPSE   ***
! **********************

! Now solve for the lapse using maximal slicing.
! But we first need to calculate all auxiliary
! variables on all grid levels, which include
! the stress-energy tensor.

  do l=0,Nl-1
     call auxiliary(l)
  end do

  call alphamaximal(0,Nl-1,"robin",1.d0)


! ***************
! ***   END   ***
! ***************

  end subroutine idata_TOVcomplex

