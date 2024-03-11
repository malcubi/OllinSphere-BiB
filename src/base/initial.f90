!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/initial.f90,v 1.89 2024/02/27 18:30:52 malcubi Exp $

  subroutine initial

! ************************
! ***   INITIAL DATA   ***
! ************************

! Include modules.

  use procinfo
  use param
  use arrays
  use derivatives
  use radialfunctions

! Extra variables.

  implicit none

  logical contains

  integer i,l

  real(8) zero,half,one,two
  real(8) smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  half = 0.5d0
  one  = 1.d0
  two  = 2.d0

  smallpi = acos(-one)


! *********************
! ***   MINKOWSKI   ***
! *********************

! By default, initial data is always set to Minkowski first.
! This is in order to avoid the code crashing because the
! metric componets are all zero.

! Lapse.

  alpha = one
  D1_alpha = zero

! Shift.

  if (shift/="none") then

     beta = zero
     D1_beta = zero

     dtbeta = zero
     D1_dtbeta = zero

     fdriver = one
     D1_fdriver = zero

  end if

! Conformal factor.

  psi = one
  phi = zero
  chi = one

  D1_psi = zero
  D1_phi = zero
  D1_chi = zero

  D2_psi = zero
  D2_phi = zero
  D2_chi = zero

! Spatial metric.  For standard runs we just set A=B=1.
! When we use the transformed radial coordinate we need
! to call "rtransform".

  if (.not.newr) then

    A = one
    B = one
  
    D1_A = zero
    D1_B = zero

  else

    call rtransform

  end if

! Extrinsic curvature.

  trK = zero
  KTA = zero
  KTB = zero

! Z4c theta variable (for BSSN it does not evolve).

  if (formulation=="z4c") then
     z4theta = zero
  end if

! Cosmological background.

  if (cosmic_run) then

     cosmobg_a = one
     cosmobg_H = zero

     cosmobg_alpha = one
     cosmobg_phi   = zero
     cosmobg_trK   = zero

  end if


! ****************************************
! ***   RESTART FROM CHECKPOINT FILE   ***
! ****************************************

! If we are restarting from a checkpoint file, call
! the restart routine and then return.

  if (idata=="checkpoint") then
     call checkpointrestart
     return
  end if


! *************************
! ***   SCHWARZSCHILD   ***
! *************************

! Initial data for a Schwarzschild black hole
! in isotropic (conformally flat) coordinates.
! In this case the conformal factor is:
!
! psi = 1 + M/(2r)
!
! Notice that for this metric the horizon is
! located (initially) at r=M/2.
!
! The isotropic lapse (for which the solution
! is static) is given by:
!
! alpha  =  (1 - M/(2r)) / (1 + M/(2r))  =  (2 - psi) / psi
!
! But notice that this lapse will cause numerical
! problems at the origin.

  if (idata=="schwarzschild") then

        print *, 'Schwarzschild initial data in isotropic cordinates with:'
        write(*,'(A,E14.8)') ' M = ',BHmass
        print *
        print *, 'The initial horizon position should be at:'
        write(*,'(A,E14.8)') ' r_H = M/2 = ',BHmass/2.d0
        print *

     if (newr) then

!       Conformal factor.

        psi = one + half*BHmass/abs(r_trans)
        phi = log(psi)

        psi2 = psi**2
        psi4 = psi**4

        D1_psi = - half*BHmass/abs(r_trans)**2*sqrt(A)
        D1_phi = D1_psi/psi*sqrt(A)

!       Lapse.

        if (ilapse=="isotropic") then
           alpha = (one - half*BHmass/abs(r_trans))/(one + half*BHmass/abs(r_trans))
        end if

     else

!       Conformal factor.
  
        psi = one + half*BHmass/abs(r)
        phi = log(psi)

        psi2 = psi**2
        psi4 = psi**4

        D1_psi = - half*BHmass/abs(r)**2
        D1_phi = D1_psi/psi

!       Lapse.

        if (ilapse=="isotropic") then
           alpha = (one - half*BHmass/abs(r))/(one + half*BHmass/abs(r))
        end if

     end if

  end if


! ******************************
! ***   TRUMPET BLACK HOLE   ***
! ******************************

! Initial data for a Schwarzschild black hole
! using the so-called "trumpet" gauge.  This
! corresponds to a static black hole with
! maximal slicing and non-trivial shift.

  if (idata=="trumpetBH") then
     call idata_trumpetBH
     phi = log(psi)
  end if


! ******************************
! ***   REISSNER-NORDSTROM   ***
! ******************************

! Charged (Reissner-Nordstrom) black hole in isotropic
! (conformally flat) coordinates. In this case the
! conformal factor is:
!
! psi  =  sqrt [ (1 + M/(2r))**2 - Q**2 / (4 r**2) ]
!      =  sqrt [ 1 + (M+Q)/2r ] [ 1 + (M-Q)/2r ]
!
! with the areal coordinate R related to the radial
! coordinate as:
!
! R  =  r psi^2  =  r [ 1 + (M+Q)/2r ] [ 1 + (M-Q)/2r ]
!
! Notice that in these coordinates the outer horizon
! is located (initially) at r=sqrt(M**2 - Q**2)/2,
! and the coordinates do not penetrate the horizon
! (smaller values of r correspond to the other side
! of the wormhole).
!
! On the other hand, the horizon mass is not given
! by M, but rather by M_H = ( sqrt(M**2 - Q**2) + M ) /2.
!
! Also, in this case there is a non-trival electric
! field that also needs to be initialized:
!
! E^r  =  Q /(r**2 psi**6)
!
! The isotropic lapse (for which the solution is
! static) is given by:
!
! alpha  =  [ ( 1 + M/(2r)) (1 - M/(2r)) + Q**2 / (4 r**2) ]
!        
!        /  [ ( 1 + M/(2r))**2 - Q**2 / (4 r**2) ]
!
! For the initial potentials, we can always take the
! vector potential equal to zero. The scalar potential Phi
! in fact depends on the initial lapse. This is because
! for Ar=0 the electric field and scalar potential
! are related through:
!
! E_r  =  - d (alpha Phi) / alpha
!            r
!
! One can show that for Reissner-Nordstrom this generally
! implies that:
!
! Phi  =  Q / (alpha R)
!
! with R the areal radius.  In isotropic coordinates
! this becomes:
!
! Phi  =  Q / (alpha r psi**2)
!
! But notice that for the isotropic lapse Phi becomes
! singular at the horizon since alpha vanishes there!

  if (idata=="reissnernordstrom") then

     if (contains(mattertype,"electric")) then

        print *, 'Reissner-Nordstrom initial data in isotropic coordinates with:'
        write(*,'(A,E14.8,A,E14.8)') ' M = ',BHmass,'     Q = ',BHcharge
        print *
        print *, 'The horizon mass should be:'
        write(*,'(A,E14.8)') ' M_H = ',(sqrt(BHmass**2-BHcharge**2)+BHmass)/2.d0
        print *, 'The initial horizon position should be at:'
        write(*,'(A,E14.8)') ' r_H = ',sqrt(BHmass**2-BHcharge**2)/2.d0
        print *, 'The areal radius for the inner and outer horizons should be:'
        write(*,'(A,E14.8,A,E14.8)') ' R_IN  = ',BHmass-sqrt(BHmass**2-BHcharge**2), &
                                 '     R_OUT = ',BHmass+sqrt(BHmass**2-BHcharge**2)
        print *

!       Conformal factor.

        psi = sqrt((one + half*BHmass/abs(r))**2 - 0.25d0*BHcharge**2/r**2)
        phi = log(psi)

        psi2 = psi**2
        psi4 = psi**4

!       Lapse.

        if (ilapse=="isotropic") then
           alpha = ((1.d0 - 0.25d0*BHmass**2/r**2) + 0.25d0*BHcharge**2/r**2)/psi2
        end if

!       Electric field (index up).

        electric = BHcharge/r**2/psi**6

!       Scalar and vector potentials.

        ePhi = BHcharge/alpha/psi2/abs(r)
        eAr  = 0.d0

     else

        print *, 'Reissner-Nordstrom initial data needs "electric" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! ************************
! ***   SCALAR PULSE   ***
! ************************

! This initial data corresponds to a small initial pulse
! in the real scalar field.

  if (idata=="scalarpulse") then

     if (contains(mattertype,"scalar")) then

        call idata_scalarpulse

     else

        print *, 'Scalar pulse initial data needs "scalar" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! ***********************
! ***   GHOST PULSE   ***
! ***********************

! This initial data corresponds to a small initial pulse
! in the ghost scalar field.

  if (idata=="ghostpulse") then

     if (contains(mattertype,"ghost")) then

        call idata_ghostpulse

     else

        print *, 'Ghost pulse initial data needs "ghost" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! **************************
! ***   GHOST WORMHOLE   ***
! **************************

! Simple wormhole sustained by a ghost scalar field.
!
! It turns out that a ghost scalar field can sustain
! static wormholes (though they are unstable). The
! example here is a simple symmetric wormhole in
! isotropic coordinates with a conformal factor
! given by:
!
! psi = sqrt(1 + 1/(4 r^2))
!
! and a scalar field given by:
!
! ghost_phi = [ 1/(2*sqrt(pi)) ] atan(r - 1/(4r))
!
! This is an exact static solution with alpha=1.
! The throat is located at r=1/2, and the asymptotic ADM
! mass is equal to 0 on both sides of the wormhole.
!
! The metric has an isometry if we take r -> 1/(4r),
! and the ghost scalar field changes sign. This does
! not change the energy density as it depends on
! the square of the derivatve of phi.

  if (idata=="ghostwormhole") then

     if (contains(mattertype,"ghost")) then

        psi = sqrt(1.d0 + 0.25d0/r**2)
        phi = log(psi)

        ghost_phi = 0.5d0/sqrt(smallpi)*atan(abs(r) - 0.25d0/abs(r))
        ghost_xi  = 2.d0/sqrt(smallpi)/(1 + 4.d0*r**2)

     else

        print *, 'Ghost wormhole initial data needs "ghost" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! *************************
! ***   COMPLEX PULSE   ***
! *************************

! This initial data corresponds to a small initial pulse
! in the complex scalar field.

  if (idata=="complexpulse") then

     if (contains(mattertype,"complex")) then

!       Non-charged initial data.  Notice that this initial data
!       can also be used for the case of a charged complex field
!       with zero initial charge denisty.  As long as the real and
!       imaginary parts are different and overlap in some region
!       we will have a non-zero initial current density.
!       (see J. Torres and M. Alcubierre, Gen. Rel. Grav. (2014) 46:1773)

        if (.not.complex_initialcharged) then

           call idata_complexpulse
 
!       Initially charged complex scalar field.  For an initially
!       charged complex scalar field we need to solve both the
!       Hamiltonian and Gauss constraints.

        else

           if (contains(mattertype,"electric")) then
              call idata_complexpulse_charged
           else
              print *, 'Charged complex scalar field initial data needs electric type matter ...'
              print *, 'Aborting! (subroutine initial)'
              print *
              call die
           end if

        end if

!    Warning.

     else

        print *, 'Complex pulse initial data needs "complex" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! *******************************
! ***   COMPLEX GHOST PULSE   ***
! *******************************

! Not yet implemented.

  if (idata=="complexghostpulse") then

     if (contains(mattertype,"complexghost")) then

        call idata_complexghostpulse

     else

        print *, 'Complex ghost pulse initial data needs "complexghost" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! ************************
! ***   NONMIN PULSE   ***
! ************************

! This initial data corresponds to a small initial pulse
! in the non-minimaly coupled scalar field.

  if (idata=="nonminpulse") then

     if (contains(mattertype,"nonmin")) then

        call idata_nonminpulse

     else

        print *, 'Nonmin pulse initial data needs "nonmin" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! ***********************
! ***   PROCA PULSE   ***
! ***********************

! This initial data corresponds to a small initial pulse
! in the proca field.

  if (idata=="procapulse") then

     if (contains(mattertype,"proca")) then

        call idata_procapulse

     else

        print *, 'Proca pulse initial data needs "proca" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! **********************
! ***   BOSON STAR   ***
! **********************

! Boson star initial data.

  if (idata=="bosonstar") then

     if (contains(mattertype,"complex")) then

        if (complexpotential=="none") then
           print *, 'For boson star initial data we need a massive scalar field.'
           print *, 'Aborting! (subroutine initial)'
           print *
           call die
        else
           if (boson_gauge=="PA") then
              call idata_BosonstarPA
           else if (boson_gauge=="CF") then
              call idata_BosonstarCF
           end if
        end if

     else

        print *, 'Boson star initial data needs "complex" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! ******************************
! ***   CHARGED BOSON STAR   ***
! ******************************

! Charged boson star initial data.

  if (idata=="chargedboson") then

     if (contains(mattertype,"complex").and.contains(mattertype,"electric")) then

       call idata_chargedboson

     else

        print *, 'Charged boson star initial data needs both "complex" and "electric" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! **********************
! ***   PROCA STAR   ***
! **********************

! Proca star initial data.

  if (idata=="procastar") then

     if (contains(mattertype,"complexproca")) then

       call idata_Procastar

     else

        print *, 'Proca star initial data needs "complexproca" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! ************************
! ***   L-PROCA STAR   ***
! ************************

! Proca star initial data.

  if (idata=="l-procastar") then

     if (contains(mattertype,"complexproca")) then

        call idata_LProcastar

     else

        print *, 'Proca star initial data needs "complexproca" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! **********************
! ***  DIRAC PULSE   ***
! **********************

  if (idata=="diracpulse") then

     if (contains(mattertype,"dirac")) then

        call idata_diracpulse

     else

        print *, 'Dirac pulse initial data needs "dirac" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! **********************
! ***  DIRAC PULSE   ***
! **********************

  if (idata=="diracstar") then

     if (contains(mattertype,"dirac")) then

        call idata_diracstar

     else

        print *, 'Dirac star initial data needs "dirac" type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! *****************************
! ***   DUST INITIAL DATA   ***
! *****************************

! Initial data for a dust spherical shell.

  if (idata=="dustshell") then

     if (contains(mattertype,"dust")) then

        call idata_dustshell

     else

        print *, 'dustshell initial data needs fluid or dust type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if



! ******************************
! ***   FLUID INITIAL DATA   ***
! ******************************

! Initial data for a fluid spherical shell.

  if (idata=="fluidshell") then

     if (contains(mattertype,"fluid")) then

        call idata_fluidshell

     else

        print *, 'Fluidshell initial data needs fluid or dust type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

! Initial data for a TOV star.

  else if (idata=="TOVstar") then

     if (contains(mattertype,"fluid")) then

        call idata_TOVstar

     else

        print *, 'TOV initial data needs fluid type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

! Blast wave initial data.

  else if (idata=="blastwave") then

     if (contains(mattertype,"fluid")) then

        call idata_blastwave

     else

        print *, 'Blastwave initial data needs fluid type matter ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! *********************
! ***   DE SITTER   ***
! *********************

! Initial data for deSitter spacetime (seen as an exponentially
! expanding flat space). Notice that this only makes sense for a
! positive cosmological constant.
!
! The initial spatial metric for deSitter is flat.  The initial
! lapse is 1, and the initial shift zero. The only difference with
! Minkowski is in the initial extrinsic curvature which is not zero.

  if (idata=="desitter") then

     if (contains(mattertype,"cosmo")) then

!       Initial lapse.

        alpha = one

!       Initial spatial metric.

        A = one
        B = one

!       Check sign of cosmological constant and calculate trK.
!       Notice that for a cosmological spacetime we always have
!       trK = - 3H, with H the Hubble constant. For deSitter
!       in particular we have H = sqrt(Lambda/3).

        if (lambda_cosmo>=zero) then
           trK = - 3.d0*sqrt(lambda_cosmo/3.d0)
        else
           print *, 'deSitter initial data needs a positive cosmological constant ...'
           print *, 'Aborting! (subroutine initial)'
           print *
           call die
        end if

!       Background 0D quantities.

        if (cosmic_run) then
           cosmobg_H   = sqrt(lambda_cosmo/3.d0)
           cosmobg_trK = - 3.d0*cosmobg_H
        end if

     else

        print *, 'deSitter initial data needs a positive cosmological constant ...'
        print *, 'Aborting! (subroutine initial)'
        print *
        call die

     end if

  end if


! ****************************************
! ***   SCALAR DARK MATTER COSMOLOGY   ***
! ****************************************

! This initial data corresponds to a flat cosmology
! with a scalar dark matter field.

  if (idata=="scalarDM") then

     if (contains(mattertype,"scalar")) then

        call idata_scalarDM

     else

        print *, 'Scalar field Dark Matter initial data needs scalar field type matter ...'
        print *, 'Aborting!  (subroutine initial)'
        print *
        call die

     end if

  end if


! ************************************************
! ***   COMPLEX SCALAR DARK MATTER COSMOLOGY   ***
! ************************************************

! This initial data corresponds to a flat cosmology
! with a complex scalar dark matter field. The initial
! data routine allows for a perturbation to be added.

  if (idata=="complexDM") then

!    Sanity check.

     if (.not.contains(mattertype,"complex")) then

        print *, 'Complex field Dark Matter initial data needs complex field type matter ...'
        print *, 'Aborting!  (subroutine initial)'
        print *
        call die

     end if

!    Time symmetric initial data with harmonic time dependence.

     if (complexDM_type=='harmonicTS') then
 
        call idata_complexDM

!    In this case we assume that the initial data has non-zero
!    momentum denisty, so we need to solve the coupled hamiltonian
!    and momentum constraints.

     else

        call idata_complexDMmom

     end if

  end if


! ********************************
! ***   AUXILIARY QUANTITIES   ***
! ********************************

! WHAT FOLLOWS MUST BE DONE AFTER ALL INITIAL DATA ROUTINES,
! SO PLEASE DO NOT PUT ANY CALLS TO INITIAL DATA AFTER THIS LINE!


! ***************************
! ***   PSI2, PSI4, CHI   ***
! ***************************

  psi2 = psi**2
  psi4 = psi**4

  chi = 1.d0/psi**chipower

! Derivatives of phi.

  diffvar => phi

  do l=0,Nl-1
     D1_phi(l,:) = diff1(l,+1)
     D2_phi(l,:) = diff2(l,+1)
  end do


! **********************************************
! ***   DERIVATIVES OF METRIC COEFFICIENTS   ***
! **********************************************

! We need to find the derivatives of the metric
! coefficients (A,B) since they are used  to
! calculate some quantities below.

  diffvar => A

  do l=0,Nl-1
     D1_A(l,:) = diff1(l,+1)
     D2_A(l,:) = diff2(l,+1)
  end do

  diffvar => B

  do l=0,Nl-1
     D1_B(l,:) = diff1(l,+1)
     D2_B(l,:) = diff2(l,+1)
  end do


! ************************************
! ***   REGULARIZATION VARIABLES   ***
! ************************************

! Klambda and lambda.

  lambda  = (one - A/B)/r**2
  Klambda = (KTA - KTB)/r**2

  D1_lambda  = - two*lambda/r  - A/B/r**2*(D1_A/A - D1_B/B)
  D1_Klambda = - two*Klambda/r + (D1_KTA - D1_KTB)/r**2

  if (regular2) then
     lambda2  =  lambda/exp(lambdapower*phi)
     Klambda2 = Klambda/exp(lambdapower*phi)
     D1_lambda2  =  D1_lambda/exp(lambdapower*phi) - lambdapower*lambda2*D1_phi
     D1_Klambda2 = D1_Klambda/exp(lambdapower*phi) - lambdapower*Klambda2*D1_phi
  end if


! **********************
! ***   BSSN Delta   ***
! **********************

! Deltar and derivative.

  Deltar = (half*D1_A/A - D1_B/B - two*r*lambda)/A

  diffvar => Deltar

  do l=0,Nl-1
     D1_Deltar(l,:) = diff1(l,-1)
  end do

! Save initial value.

  Deltar0 = Deltar


! *************************
! ***   INITIAL LAPSE   ***
! *************************

! This is a place for special types of initial lapse
! that are quite generic. Notice that by default we
! have ilapse="none", so that this is ignored and the
! lapse retains the value it was given above.

  if (ilapse/="none") then

     if (ilapse=="one") then

        alpha = one

     else if (ilapse=="psiminus2") then

        if (cosmic_run) then
           do l=0,Nl-1
              alpha(l,:) = one/exp(2.d0*(phi(l,:)-cosmobg_phi(l)))
           end do
        else
           alpha = one/psi2
        end if

     else if (ilapse=="psiminus4") then

        if (cosmic_run) then
           do l=0,Nl-1
              alpha(l,:) = one/exp(4.d0*(phi(l,:)-cosmobg_phi(l)))
           end do
        else
           alpha = psi4
        end if

     else if (ilapse=="psi2") then

        if (cosmic_run) then
           do l=0,Nl-1
              alpha(l,:) = exp(2.d0*(phi(l,:)-cosmobg_phi(l)))
           end do
        else
           alpha = psi2
        end if

     else if (ilapse=="psi4") then

        if (cosmic_run) then
           do l=0,Nl-1
              alpha(l,:) = exp(4.d0*(phi(l,:)-cosmobg_phi(l)))
           end do
        else
           alpha = psi4
        end if

     end if

!    Derivatives.

     diffvar => alpha

     do l=0,Nl-1
        D1_alpha(l,:) = diff1(l,+1)
        D2_alpha(l,:) = diff2(l,+1)
     end do

  end if


! ******************************
! ***   LAPSE PERTURBATION   ***
! ******************************

! Add perturbation to initial lapse. The lapse perturbation must be even!

! Gaussian perturbation.

  if (lapsepert/="none") then

!    Gaussian profile.

     if (lapsepert=="gaussian") then

        if (lapse_r0==0.d0) then
           alpha = alpha + gaussian(lapse_a0,0.d0,lapse_s0)
        else
           alpha = alpha + gaussian(lapse_a0,+lapse_r0,lapse_s0) &
                         + gaussian(lapse_a0,-lapse_r0,lapse_s0)
        end if

!    Smooth top-hat profile.

     else if (lapsepert=="tophat") then

        if (lapse_r0==0.d0) then
           alpha = alpha + tophat(lapse_a0,0.d0,lapse_s0,lapse_t0)
        else
           alpha = alpha + tophat(lapse_a0,+lapse_r0,lapse_s0,lapse_t0) &
                         + tophat(lapse_a0,-lapse_r0,lapse_s0,lapse_t0) 
        end if

     end if

!    The spatial derivative is calculated with finite differences.

     diffvar => alpha
 
     do l=0,Nl-1
        D1_alpha(l,:) = diff1(l,+1)
     end do

  end if


! ******************************
! ***   SHIFT PERTURBATION   ***
! ******************************

! Add perturbation to initial shift. The shift perturbation must be odd!

  if (shiftpert/="none") then

!    Sanity check.

     if ((shift=="none").or.(shift=="zero")) then

        print *, 'shift=(none,zero) is incompatible with adding a perturbation.'
        print *, 'Aborting!  (subroutine initial)'
        print *
        call die

!    For a shift perturbation shift_r0 must be non-zero.

     else if (shift_r0==0.d0) then

        print *, 'The shift perturbation must be odd, so shift_r0 must be non-zero.'
        print *, 'Aborting!  (subroutine initial)'
        print *
        call die

!    Gaussian perturbation.

     else if (shiftpert=="gaussian") then

        beta = beta + gaussian(shift_a0,+shift_r0,shift_s0) &
                    - gaussian(shift_a0,-shift_r0,shift_s0)

!    Smooth top-hat profile.

     else if (shiftpert=="tophat") then

        beta = beta + tophat(shift_a0,+shift_r0,shift_s0,shift_t0) &
                    - tophat(shift_a0,-shift_r0,shift_s0,shift_t0)

     end if

!    The spatial derivative is calculated with finite differences.

     diffvar => beta
 
     do l=0,Nl-1
        D1_beta(l,:) = diff1(l,-1)
     end do

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine initial
