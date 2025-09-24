!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/analysis_geometry.f90,v 1.34 2025/09/24 17:22:23 malcubi Exp $

  subroutine analysis_geometry

! ********************************************************
! ***   CALCULATION OF ANALYSIS VARIABLES FOR OUTPUT   ***
! ********************************************************

! Include modules.

  use arrays
  use param
  use derivatives
  use integrals

! Extra variables.

  implicit none

  integer i,l

  real(8) half,third,one,two,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  half  = 0.5d0
  third = 1.d0/3.d0

  one = 1.d0
  two = 2.d0

  smallpi = acos(-1.d0)


! ************************
! ***   AREAL RADIUS   ***
! ************************

! This is the Schwarzschild (or areal) radius,
! and is given by:
!
! r_area  =  r psi**2 sqrt(B)

  r_area = r*psi**2*sqrt(B)

! Derivative.

  diffvar => r_area

  do l=0,Nl-1
     D1_r_area(l,:) = diff1(l,-1)
  end do

! Radial metric for areal radius. This can be useful in
! many cases, and is given by:
!
!                                            2
! GRR = (A/B) / [ 1 + r ( 2 d phi + d B/2B) ]
!                            r       r

  if (allocated(GRR)) then
     GRR = (A/B)/(1.d0 + r*(2.d0*D1_phi + 0.5d0*D1_B/B))**2
  end if

! Embedding height function h. The derivative of h
! with respect to r is given by:
!
!                                                                                 1/2
! dh/dr  =  a_area [ (A/B-1)/r - (2 d phi + d B/(2B)) (2/r + 2 d phi + d B/(2B)) ]
!                                    r       r                  r       r
!
! Notice that for this to make sense we need the
! whole expression inside the square root to be
! positive. In particular, for conformally flat data
! (A=B=1), we need to have :
!
! -1/r  <  dphi/dr  <  0

  if (allocated(hembed)) then

     auxarray = r_area*sqrt(abs((A/B - 1.d0)/r - (2.d0*D1_phi + 0.5d0*D1_B/B) &
              *(2.d0/r + 2.d0*D1_phi + 0.d0*D1_B/B)))

     intvar => auxarray

     do l=Nl-1,0,-1
        hembed(l,:) = integral(l)
     end do

!    Restrict integral.

     if (Nl>1) then
        intvar => hembed
        call restrictintegral
     end if

  end if


! ******************************
! ***   SCHWARZSCHILD MASS   ***
! ******************************

! This is the Schwarzschild mass.  It is obtained by noticing
! that for the Schwarzschild metric we have:
!
! g   =  1 / (1 - 2M/R)
!  RR
!
! with R=r_area the Schwarzschild or areal radius (see above).
! This implies:
!
!                                             2       4
! M  =  R/2 ( 1 - 1/g  )  =  R/2 ( 1 - (dR/dr) / A psi )
!                    RR
!
! with:
!
!            1/2    2
! dR/dr  =  B    psi  ( 1 + r d B / (2 B)  +  2 r d psi / psi )
!                              r                   r
!
! So that finally:
!
!            2  1/2                                                  2
! M  =  r psi  B   / 2  [ 1 - (B/A) ( 1 + r d B / (2 B) + 2 r d phi )  ]
!                                            r                 r

  if (allocated(mass_sch)) then

     mass_sch = half*r_area*(one - B/A*(one + r*(half*D1_B/B + two*D1_phi))**2)

!    Fix boundary.

     mass_sch(:,Nr) = 2.d0*mass_sch(:,Nr-1) - mass_sch(:,Nr-2)

!    Restrict.

     restrictvar => mass_sch

     do l=Nl-1,1,-1
        call restrict(l,.false.)
     end do

!    Symmetries at origin.

     do l=0,Nl-1
        do i=1,ghost
           mass_sch(l,1-i) = mass_sch(l,i)
        end do
     end do

  end if


! ****************************************************
! ***   EXPANSION OF OUTGOING/INGOING NULL LINES   ***
! *****************************************************

! The expansion of outward null lines is given by:
!
! theta = (2/r + d B / B + 4 d phi) / (psi**2 sqrt(A)) - 2 K_theta^theta
!                 r           r
!
! We multiply it by r so that it is constant and equal to +2 for
! a flat spacetime.

  if (allocated(expansion)) then
     expansion = r*((2.d0/r + D1_B/B + 4.d0*D1_phi)/sqrt(A)/psi**2 - 2.d0*KBPHYS)
  end if

! Expansion of ingoing null lines.  This is the same as above, but changing
! the sign of the first term. For a flat spacetime it should be -2.

  if (allocated(expansion_in)) then
     expansion_in = - r*((2.d0/r + D1_B/B + 4.d0*D1_phi)/sqrt(A)/psi**2 + 2.d0*KBPHYS)
  end if


! *************************************
! ***   COORDINATE SPEED OF LIGHT   ***
! *************************************

! The coordinate speed of light is given by:  vl = alpha/sqrt(A)/psi**2.
! This does not take into account the shift.

  if (allocated(vlight)) then
     vlight = alpha/sqrt(abs(A))/psi**2
  end if


! **********************************
! ***   BONA-MASSO GAUGE SPEED   ***
! **********************************

  if (allocated(vgauge)) then
     if ((slicing=="harmonic").or.(slicing=="1+log").or. &
         (slicing=="shockavoid").or.(slicing=="alphaminus2")) then
        vgauge = sqrt(falpha/abs(A))/psi**2
     end if
  end if


! *********************************
! ***   GEOMETRIC EIGENFIELDS   ***
! *********************************

! Lapse gauge eigenfields.

  if (slicing/="maximal") then

     if (allocated(wp_gauge)) then
        wp_gauge = alpha*sqrt(falpha*A)*psi**2*trK + D1_alpha
     end if

     if (allocated(wm_gauge)) then
        wm_gauge = alpha*sqrt(falpha*A)*psi**2*trK - D1_alpha
     end if

  end if

! Eigenfields associated with KTA.

  if (allocated(wp_eta)) then
     wp_eta = sqrt(A)*psi**2*sqrt((2.d0*eta-1.d0)/3.d0)*(KTA - 2.d0*trK/3.d0) &
            + ((D1_A/A - D1_B/B) - 2.d0*A*Deltar + 4.d0*D1_phi)/3.d0
  end if

  if (allocated(wm_eta)) then
     wm_eta = sqrt(A)*psi**2*sqrt((2.d0*eta-1.d0)/3.d0)*(KTA - 2.d0*trK/3.d0) &
            - ((D1_A/A - D1_B/B) - 2.d0*A*Deltar + 4.d0*D1_phi)/3.d0
  end if

! Non-propagating eigenfield.

  if (allocated(wDelta)) then
     wDelta = 4.d0*eta*D1_phi - half*(eta - 2.d0)*D1_A/A - A*Deltar
  end if


! ***********************
! ***   WEYL TENSOR   ***
! ***********************

! The electric part of the Weyl tensor is defined as:
!
!  i      i           i      i  k
! E   =  R    +  trK K   -  K  K
!  j      j           j      k  j
!
!                i         i
!     -  4 pi [ S  +  delta  ( 4 rho - S) / 3 ]
!                j         j
!
! Notice that the Hamiltonian constraint implies that
! this tensor is traceless.

! Mixed radial component: EWEYLA = E^r_r.

  EWEYLA = RICA + trK*KAPHYS - KAPHYS**2

  if (mattertype/="vacuum") then
     EWEYLA = EWEYLA - 4.d0*smallpi*(SAA + third*(4.d0*rho - (SAA + 2.d0*SBB)))
  end if

! Mixed angular component: EWEYLB = E^theta_theta.
! Here we use the fact that the electric tensor is traceless.

  EWEYLB = - 0.5d0*EWEYLA


! ********************************
! ***   CURVATURE INVARIANTS   ***
! ********************************

  if (allocated(invariantI).or.allocated(invariantJ)) then
     call invariants
  end if


! ***************************
! ***   4D RICCI SCALAR   ***
! ***************************

  if (allocated(Ricci4D)) then
     if (mattertype/="vacuum") then
        Ricci4D = - 8.d0*smallpi*(-rho + trS)
     else
        Ricci4D = 0.d0
     end if
  end if

! 4D Ricci Scalar in the Einstein frame.

  if (allocated(Ricci4D_E)) then
     if (mattertype=="nonmin".and.shift=="none") then

        Ricci4D_E = Ricci4D*(1.d0/nonmin_f+1.5d0*nonmin_fp**2/nonmin_f**2) &
                  - 1.5d0*(2.d0*nonmin_f*nonmin_fpp-nonmin_fp**2)/nonmin_f**3 &
                   *(D1_nonmin_phi**2/(A*psi4)-nonmin_pi**2)

        Ricci4D_E = Ricci4D_E/(8.d0*smallpi)           

     end if

  end if


! *********************
! ***   COSMOLOGY   ***
! *********************

! For cosmological runs, we define some quantities
! subtracting the background for easier analysis.
! We also calculate the background speed of light.

  if (cosmic_run) then
     do l=0,Nl-1

!       Perturbed lapse.

        alpha_pert(l,:) = alpha(l,:) - cosmobg_alpha(l)

!       Perturbed conformal factor.

        phi_pert(l,:) = phi(l,:) - cosmobg_phi(l)

!       Perturbed trace of extrinsic curvature.

        trK_pert(l,:) = trK(l,:) - cosmobg_trK(l)

!       Speed of light.

        cosmobg_vlight(l) = cosmobg_alpha(l)/cosmobg_a(l)

     end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine analysis_geometry

