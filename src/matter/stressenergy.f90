
  subroutine stressenergy(l)

! *********************************************
! ***   EVALUATION OF STRESS-ENERGY TERMS   ***
! *********************************************

! This subroutine evaluates the stress-energy variables
! that appear in the Einstein equations.
!
! These quantities are defined in general as:
! 
!          mu  nu
! rho  =  n   n   T             (energy density)
!                  mu nu
!
!            mu  nu
! J    =  - n   P   T           (momentum density)
!  i             i   mu nu
!
!          mu  nu
! S    =  P   P    T            (stress tensor)
!  ij      i   j    mu nu
!
! with T_{mu,nu} the stress-energy tensor, n^mu the normal unit
! vector to the spatial hypersurfaces and P^mu_nu the projector
! operator onto the hypersurfaces.
!
! Notice that in spherical symmetry there is only one independent
! component of the momentum density J_i (JA), and two independent
! components of the stress tensor: SAA=S^r_r and SBB=S^theta_theta.
!
!
! IMPORTANT:  For all types of matter, we always ADD to the values
! of the stress-energy variables.  This is because we might want
! to have more than one type of matter present in a given simulation.

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  logical contains
  logical usediracPi

  integer i,l

  real(8) zero,half,one,two
  real(8) third,smallpi
  real(8) rhoatmos,Eatmos,patmos
  real(8) aux


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  half = 0.5d0
  one  = 1.d0
  two  = 2.d0
  third = 1.d0/3.d0

  smallpi = acos(-one)


! *************************************************
! ***   MAKE SURE CONFORMAL FACTOR IS UPDATED   ***
! *************************************************

! This routine is typically called before updating the
! auxiliary geometry variables, so psi should not be
! assumed to have the correct values.

  if (chimethod) then
     phi(l,:) = - log(abs(chi(l,:)))/dble(chipower)
  end if

  psi(l,:) = exp(phi(l,:))

  psi2(l,:) = psi(l,:)**2
  psi4(l,:) = psi(l,:)**4


! ********************************
! ***   INITIALIZE TO VACUUM   ***
! ********************************

! Energy density.

  rho(l,:) = zero

! Momentum density.

  JA(l,:) = zero

! Stress tensor.

  SAA(l,:) = zero
  SBB(l,:) = zero

! SLL = (SAA - SBB)/r**2.

  if (.not.nolambda) then
     SLL(l,:) = zero
  end if

! Electric charge and current.

  if (contains(mattertype,"electric")) then
     echarge  = zero
     ecurrent = zero
  end if

! Cosmological background.

  if (cosmic_run) then
     cosmobg_rho(l) = zero
     cosmobg_p(l)   = zero
  end if


! *********************************
! ***   COSMOLOGICAL CONSTANT   ***
! *********************************

! In this case we have:
!
! rho  = - SAA  =  - SBB  = Lambda / (8 pi)
!
! All other components of the stress-energy are zero.

  if (contains(mattertype,"cosmo")) then

!    Energy density (the momentum density is zero for a
!    cosmological constant, so we don't touch it).

     rho(l,:) = rho(l,:) + lambda_cosmo/8.d0/smallpi

!    Stress tensor.  Notice that since both components
!    are equal, there is no contribution to SLL so
!    we don't touch it.

     SAA(l,:) = SAA(l,:) - lambda_cosmo/8.d0/smallpi
     SBB(l,:) = SBB(l,:) - lambda_cosmo/8.d0/smallpi

!    SLL = (SAA - SBB)/r**2.

     if (.not.nolambda) then
        SLL(l,:) = SLL(l,:)
     end if

!    Cosmological background.

     if (cosmic_run) then

!       Energy density.

        cosmobg_rho(l) = cosmobg_rho(l) + lambda_cosmo/8.d0/smallpi

!       Pressure.

        cosmobg_p(l)   = cosmobg_p(l)   - lambda_cosmo/8.d0/smallpi

     end if

  end if


! ************************
! ***   SCALAR FIELD   ***
! ************************

! The stress-energy tensor for a scalar field "phi" with a
! self-interaction potential V has the form:
!
!                                     /  beta                     \
! T       =  d  phi d  phi  -  g      | d    phi d    phi  +  2 V | / 2
!  mu nu      mu     nu         mu nu \           beta            /
!
!
! From this one gets the following values of {rho,JA,SAA,SBB}:
!
!
!          mu  nu                   2      2       4
! rho  =  n   n   T      =  1/2 ( pi  +  xi / A psi  )  +  V
!                  mu nu
!
!
!                 mu
! JA   =  - n   T        =  - pi xi
!             mu    r
!
!
!                    4              2      2       4
! SAA  =  T   / A psi    =  1/2 ( pi  +  xi / A psi  )  -  V
!           rr
!
!                            2      4              2      2       4
! SBB  =  T             / ( r  B psi )  =  1/2 ( pi  -  xi / A psi  )  -  V
!           theta theta

  if (contains(mattertype,"scalar")) then

!    Energy density.

     rho(l,:) = rho(l,:) + half*(scalar_pi(l,:)**2 + scalar_xi(l,:)**2/(A(l,:)*psi4(l,:))) &
              + scalar_V(l,:)

!    Momentum density (index down).

     JA(l,:) = JA(l,:) - scalar_pi(l,:)*scalar_xi(l,:)

!    Stress tensor.

     SAA(l,:) = SAA(l,:) + half*(scalar_pi(l,:)**2 + scalar_xi(l,:)**2/(A(l,:)*psi4(l,:))) &
              - scalar_V(l,:)

     SBB(l,:) = SBB(l,:) + half*(scalar_pi(l,:)**2 - scalar_xi(l,:)**2/(A(l,:)*psi4(l,:))) &
              - scalar_V(l,:)

!    SLL = (SAA - SBB)/r**2.

     if (.not.nolambda) then
        SLL(l,:) = SLL(l,:) + scalar_xi(l,:)**2/(A(l,:)*psi4(l,:))/r(l,:)**2
     end if

!    Cosmological background.

     if (cosmic_run) then

!       Energy density.

        cosmobg_rho(l) = cosmobg_rho(l) + half*cosmobg_scalar_pi(l)**2 + cosmobg_scalar_V(l)

!       Pressure.

        cosmobg_p(l) = cosmobg_p(l) + half*cosmobg_scalar_pi(l)**2 - cosmobg_scalar_V(l)

     end if

  end if


! ******************************
! ***   GHOST SCALAR FIELD   ***
! ******************************

! The stress-energy tensor for a ghost field is identical
! to that of the scalar field, but with a kinetic term of
! opposite sign.
!
! Notice that the convention of the code is such that the
! potential terms have THE SAME sign as for a standard scalar
! field, but the potential is in fact negative, so the whole
! stress-energy tensor has opposite sign to the standard case.

  if (contains(mattertype,"ghost")) then

!    Energy density.

     rho(l,:) = rho(l,:) - half*(ghost_pi(l,:)**2 + ghost_xi(l,:)**2/(A(l,:)*psi4(l,:))) &
              + ghost_V(l,:)

!    Momentum density (index down).

     JA(l,:) = JA(l,:) + ghost_pi(l,:)*ghost_xi(l,:)

!    Stress tensor.

     SAA(l,:) = SAA(l,:) - half*(ghost_pi(l,:)**2 + ghost_xi(l,:)**2/(A(l,:)*psi4(l,:))) &
              - ghost_V(l,:)

     SBB(l,:) = SBB(l,:) - half*(ghost_pi(l,:)**2 - ghost_xi(l,:)**2/(A(l,:)*psi4(l,:))) &
              - ghost_V(l,:)

!    SLL = (SAA - SBB)/r**2.

     if (.not.nolambda) then
        SLL(l,:) = SLL(l,:) - ghost_xi(l,:)**2/(A(l,:)*psi4(l,:))/r(l,:)**2
     end if

  end if


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

! The stress-energy tensor for a complex scalar field "phi" with a
! self-interaction potential V has the form:
!
!                  *                  /  beta   *                 \
! T       =  d  phi d  phi  -  g      | d    phi d    phi  +  2 V | / 2
!  mu nu     (mu     nu)        mu nu \           beta            /
!
! From this one gets the same values of {rho,JA,SAA,SBB} as in the case
! of a real scalar field but twice (for real and immaginary parts).
!
!
! NOTE FOR A CHARGED SCALAR FIELD:
!
! For a charged scalar field we need to use the EM "covariant derivatives":

! d -> D = d + iqA  (ommiting indexes)
!
! it is pretty straightforward to check the relations:
!
!   mu   
!  n  D  phi =  Pi - iq ePhi phi  =: gPi
!      mu
!
!  D phi     =  Xi + iq A phi     =: gXi
!   r                    r
!
! So the matter terms are exactly the same after the substitution:
!
! Pi, Xi -> gPi, gXi

  if (contains(mattertype,"complex")) then

!    First we calculate the auxiliar gauge derivatives (gXi,gPi).
!    These are the gauge derivatives for the case of a charged
!    scalar field.
     
     complex_gxiR = complex_xiR
     complex_gxiI = complex_xiI
     
     complex_gpiR = complex_piR
     complex_gpiI = complex_piI

     if (contains(mattertype,"electric")) then

        complex_gxiR = complex_gxiR - complex_q*eAr*complex_phiI
        complex_gxiI = complex_gxiI + complex_q*eAr*complex_phiR

        complex_gpiR = complex_gpiR + complex_q*ePhi*complex_phiI
        complex_gpiI = complex_gpiI - complex_q*ePhi*complex_phiR

     end if

!    Energy density.

     complex_rhotot(l,:) = half*(complex_gpiR(l,:)**2 + complex_gxiR(l,:)**2/(A(l,:)*psi4(l,:))) &
                         + half*(complex_gpiI(l,:)**2 + complex_gxiI(l,:)**2/(A(l,:)*psi4(l,:))) &
                         + complex_V(l,:)

     rho(l,:) = rho(l,:) + complex_rhotot(l,:)

!    Momentum density (index down).

     JA(l,:) = JA(l,:) - complex_gpiR(l,:)*complex_gxiR(l,:) - complex_gpiI(l,:)*complex_gxiI(l,:)

!    Stress tensor.

     SAA(l,:) = SAA(l,:) + half*(complex_gpiR(l,:)**2 + complex_gxiR(l,:)**2/(A(l,:)*psi4(l,:))) &
                         + half*(complex_gpiI(l,:)**2 + complex_gxiI(l,:)**2/(A(l,:)*psi4(l,:))) &
                         - complex_V(l,:)

     SBB(l,:) = SBB(l,:) + half*(complex_gpiR(l,:)**2 - complex_gxiR(l,:)**2/(A(l,:)*psi4(l,:))) &
                         + half*(complex_gpiI(l,:)**2 - complex_gxiI(l,:)**2/(A(l,:)*psi4(l,:))) &
                         - complex_V(l,:)

!    SLL = (SAA - SBB)/r**2

     if (.not.nolambda) then
        SLL(l,:) = SLL(l,:) + (complex_gxiR(l,:)**2 + complex_gxiI(l,:)**2)/(A(l,:)*psi4(l,:))/r(l,:)**2
     end if

!    Angular momentum.  For non-zero angular momentum we need
!    to add the following terms:
!
!    rho  ->  rho + [l(l+1)/2] |phi|^2 / ( r^2 B psi^4)
!
!    SAA  ->  SAA - [l(l+1)/2] |phi|^2 / ( r^2 B psi^4)
!
!    Notice that there is no extra term for SBB.

     if (complex_l>0) then

        rho(l,:) = rho(l,:) + half*dble(complex_l*(complex_l+1))/(B(l,:)*psi4(l,:))/r(l,:)**2 &
                 *(complex_phiR(l,:)**2 + complex_phiI(l,:)**2)

        SAA(l,:) = SAA(l,:) - half*dble(complex_l*(complex_l+1))/(B(l,:)*psi4(l,:))/r(l,:)**2 &
                 *(complex_phiR(l,:)**2 + complex_phiI(l,:)**2)

!       SLL = (SAA - SBB)/r**2.

        if (.not.nolambda) then
           SLL(l,:) = SLL(l,:) - half*dble(complex_l*(complex_l+1))/(B(l,:)*psi4(l,:))/r(l,:)**4 &
                    *(complex_phiR(l,:)**2 + complex_phiI(l,:)**2)
        end if

     end if

!    Boson density and current.  These are calculated here and not
!    in the analysis routine since for charged fields we need the
!    current for the source of the electric field.
!
!    For a complex scalar field there is a conserved Noether
!    current given by:
!
!                         *                       *
!    j    =  (i/2)  (  phi  d  phi  -  phi d   phi )
!     mu                     mu             mu
!
!         =  ( phiI  d  phiR  -  phiR  d  phiI )
!                     mu                mu
!
!    From this we can define a boson density as:
!
!                   mu
!    rho      =  - n   j   =  phiR*piI  -  phiI*piR
!       boson           mu
!
!    and a boson "flux" as:
!
!     i              ik
!    j        =  gamma   j
!     boson               k
!
!    which implies that the radial component is:
!
!     r               rr
!    j        =  gamma   j   =  (phiI*xiR  -  phiR*xiI) / (A psi^4)
!     boson               r
!
!    Notice that the code calls the boson density rho_boson "complex_Bdens",
!    and the radial flux j^r_boson "complex_Bflux".
!
!    Also, for a charged field we need to replace (pi,xi) with (gpi,gxi).

     complex_Bdens = (complex_phiR*complex_gpiI - complex_phiI*complex_gpiR)
     complex_Bflux = (complex_phiI*complex_gxiR - complex_phiR*complex_gxiI)/(psi4*A)

!    Charge and current density (index up).

     if (contains(mattertype,"electric")) then 
        echarge  = echarge  + complex_q*complex_Bdens
        ecurrent = ecurrent + complex_q*complex_Bflux
     end if

!    Cosmological background.

     if (cosmic_run) then

!       Energy density.

        cosmobg_rho(l) = cosmobg_rho(l) + half*(cosmobg_complex_piR(l)**2 + cosmobg_complex_piI(l)**2) &
                       + cosmobg_complex_V(l)

!       Pressure.

        cosmobg_p(l) = cosmobg_p(l) + half*(cosmobg_complex_piR(l)**2 + cosmobg_complex_piI(l)**2) &
                     - cosmobg_complex_V(l)

     end if

  end if


! *******************************
! ***   COMPLEX GHOST FIELD   ***
! *******************************

! The stress-energy tensor for a ghost field is identical
! to that of the scalar field, but with a kinetic term of
! opposite sign.
!
! Notice that the convention of the code is such that the
! potential terms have THE SAME sign as for a standard scalar
! field, but the potential is in fact negative, so the whole
! stress-energy tensor has opposite sign to the standard case.
!
! At the moment the expressions the complex ghost field
! is assumed t have NO electric charge and NO angular
! momentum contributions.  I might change this latter.

  if (contains(mattertype,"complexghost")) then

!    Energy density.

     rho(l,:) = rho(l,:) - half*(complexghost_piR(l,:)**2 + complexghost_xiR(l,:)**2/(A(l,:)*psi4(l,:))) &
                         - half*(complexghost_piI(l,:)**2 + complexghost_xiI(l,:)**2/(A(l,:)*psi4(l,:))) &
                         + complexghost_V(l,:)

!    Momentum density (index down).

     JA(l,:) = JA(l,:) + complexghost_piR(l,:)*complexghost_xiR(l,:) + complexghost_piI(l,:)*complexghost_xiI(l,:)

!    Stress tensor.

     SAA(l,:) = SAA(l,:) - half*(complexghost_piR(l,:)**2 + complexghost_xiR(l,:)**2/(A(l,:)*psi4(l,:))) &
                         - half*(complexghost_piI(l,:)**2 + complexghost_xiI(l,:)**2/(A(l,:)*psi4(l,:))) &
                         - complexghost_V(l,:)

     SBB(l,:) = SBB(l,:) - half*(complexghost_piR(l,:)**2 - complexghost_xiR(l,:)**2/(A(l,:)*psi4(l,:))) &
                         - half*(complexghost_piI(l,:)**2 - complexghost_xiI(l,:)**2/(A(l,:)*psi4(l,:))) &
                         - complexghost_V(l,:)

!    SLL = (SAA - SBB)/r**2.

     if (.not.nolambda) then
        SLL(l,:) = SLL(l,:) - (complexghost_xiR(l,:)**2 + complexghost_xiI(l,:)**2)/(A(l,:)*psi4(l,:))/r(l,:)**2
     end if

  end if


! **************************
! ***   ELECTRIC FIELD   ***
! **************************

! Stress-energy quantities for an electric field:
!
!                         4         2
! rho  = + 1/(8 pi)  A psi  electric
!
!
! JA   =   0   (no magnetic field in spherical symmetry)
!
!                         4         2
! SAA  = - 1/(8 pi)  A psi  electric
!
!                         4         2
! SBB  = + 1/(8 pi)  A psi  electric
!
!
! Notice that:  trS = SAA + 2 SBB = rho.

  if (contains(mattertype,"electric")) then

!    Energy density.

     rho(l,:) = rho(l,:) + (0.125d0/smallpi)*A(l,:)*psi4(l,:)*electric(l,:)**2

!    Momentum density (index down).  There is no contribution from the electric
!    field to the momentum density (this is because in spherical symmetry there
!    is no magnetic field).

     JA(l,:) = JA(l,:)

!    Stress tensor.

     SAA(l,:) = SAA(l,:) - (0.125d0/smallpi)*A(l,:)*psi4(l,:)*electric(l,:)**2
     SBB(l,:) = SBB(l,:) + (0.125d0/smallpi)*A(l,:)*psi4(l,:)*electric(l,:)**2

!    SLL = (SAA - SBB)/r**2.

     if (.not.nolambda) then
        SLL(l,:) = SLL(l,:) - (0.25d0/smallpi)*A(l,:)*psi4(l,:)*electric(l,:)**2/r(l,:)**2
     end if

  end if


! ****************************
! ***   REAL PROCA FIELD   ***
! ****************************

! The Proca field is a massive EM field (remember that both JA and ProcaA
! have index down, while procaE has index up):
!
!                           4       2               2        2       4            2
! rho  = + 1/(8 pi)  [ A psi  procaE   +  proca_mass ( procaA / A psi  +  procaPhi  ) ]
!
!                             2
! JA   = + 1/(4 pi) proca_mass  proca_Phi  Proca_A
!
!                           4       2               2        2       4            2
! SAA  = - 1/(8 pi)  [ A psi  procaE   -  proca_mass ( procaA / A psi  +  procaPhi  ) ]
!
!                           4       2               2        2       4            2 
! SBB  = + 1/(8 pi)  [ A psi  procaE   -  proca_mass ( procaA / A psi  -  procaPhi  ) ]
!
!
! In the previous expressions we have:
!
! procaE:      Electric part of the Proca field (index up).
! procaA:      Vector potential (index down).
! procaPhi:    Scalar potential.
! proca_mass:  Mass parameter.

  if (contains(mattertype,"proca")) then

!    Energy density.

     rho(l,:) = rho(l,:) + (0.125d0/smallpi)*(A(l,:)*psi4(l,:)*procaE(l,:)**2 &
              + proca_mass**2*(procaA(l,:)**2/(A(l,:)*psi4(l,:)) + procaPhi(l,:)**2))

!    Momentum density (index down).

     JA(l,:) = JA(l,:) + (0.25d0/smallpi)*proca_mass**2*procaPhi(l,:)*procaA(l,:)

!    Stress tensor.

     SAA(l,:) = SAA(l,:) - (0.125d0/smallpi)*(A(l,:)*psi4(l,:)*procaE(l,:)**2 &
              - proca_mass**2*(procaA(l,:)**2/(A(l,:)*psi4(l,:)) + procaPhi(l,:)**2))

     SBB(l,:) = SBB(l,:) + (0.125d0/smallpi)*(A(l,:)*psi4(l,:)*procaE(l,:)**2 &
              - proca_mass**2*(procaA(l,:)**2/(A(l,:)*psi4(l,:)) - procaPhi(l,:)**2))

!    SLL = (SAA - SBB)/r**2.

     if (.not.nolambda) then
        SLL(l,:) = SLL(l,:) - (0.25d0/smallpi)*(A(l,:)*psi4(l,:)*procaE(l,:)**2 &
                 - proca_mass**2*procaA(l,:)**2/(A(l,:)*psi4(l,:)))/r(l,:)**2
     end if

  end if


! *******************************
! ***   COMPLEX PROCA FIELD   ***
! *******************************

! The Proca field is a massive EM field (remember that both JA and ProcaA
! have index down, while procaE has index up).
!
! The expressiones for the complex Proca field are the same as for the
! real field above, one just needs to add contributions from the real
! and imaginary parts.

  if (contains(mattertype,"complexproca")) then

!    Energy density.

     rho(l,:) = rho(l,:) + (0.125d0/smallpi) &
              *(A(l,:)*psi4(l,:)*(cprocaE_R(l,:)**2 + cprocaE_I(l,:)**2) &
              + cproca_mass**2*((cprocaA_R(l,:)**2 + cprocaA_I(l,:)**2)/(A(l,:)*psi4(l,:)) &
              + cprocaPhi_R(l,:)**2 + cprocaPhi_I(l,:)**2))

!    Radial momentum density (index down).

     JA(l,:) = JA(l,:) + (0.25d0/smallpi)*cproca_mass**2 &
             *(cprocaPhi_R(l,:)*cprocaA_R(l,:) + cprocaPhi_I(l,:)*cprocaA_I(l,:))

!    Stress tensor.

     SAA(l,:) = SAA(l,:) - (0.125d0/smallpi) &
              *(A(l,:)*psi4(l,:)*(cprocaE_R(l,:)**2 + cprocaE_I(l,:)**2) &
              - cproca_mass**2*((cprocaA_R(l,:)**2 + cprocaA_I(l,:)**2)/(A(l,:)*psi4(l,:)) &
              + (cprocaPhi_R(l,:)**2 + cprocaPhi_I(l,:)**2)))

     SBB(l,:) = SBB(l,:) + (0.125d0/smallpi) &
              *(A(l,:)*psi4(l,:)*(cprocaE_R(l,:)**2 + cprocaE_I(l,:)**2) &
              - cproca_mass**2*((cprocaA_R(l,:)**2 + cprocaA_I(l,:)**2)/(A(l,:)*psi4(l,:)) &
              - (cprocaPhi_R(l,:)**2 + cprocaPhi_I(l,:)**2)))

!    SLL = (SAA - SBB)/r**2.

     if (.not.nolambda) then
        SLL(l,:) = SLL(l,:) - (0.25d0/smallpi) &
                 *(A(l,:)*psi4(l,:)*(cprocaE_R(l,:)**2 + cprocaE_I(l,:)**2) &
                 - cproca_mass**2*(cprocaA_R(l,:)**2 + cprocaA_I(l,:)**2)/(A(l,:)*psi4(l,:)))/r(l,:)**2
     end if

!    Aditional terms added for non-zero angular momentum of constituent fields:
!
!                                                   2        2            2                    2
!    rho  ->  rho +  1/(8 pi) l(l+1) [ + cproca_mass |ProcaB|  + |procaXi|  + |procaA - procaG| /(A psi^4) ] / ( r^2 B psi^4)
!
!
!    JA   ->   JA -  1/(4 pi) l(l+1) [ procaXi_R ( procaA_R - procaG_R ) + procaXi_I (procaA_I - procaG_I) ] / ( r^2 B psi^4)
!
!                                                   2        2            2                    2
!    SAA  ->  SAA +  1/(8 pi) l(l+1) [ - cproca_mass |ProcaB|  + |procaXi|  + |procaA - procaG| /(A psi^4) ] / ( r^2 B psi^4)
!
!
!    Notice that there is no extra term for SBB.
!
!                                                   2        2            2                    2
!    SLL ->   SLL +  1/(8 pi) l(l+1) [ - cproca_mass |ProcaB|  + |procaXi|  + |procaA - procaG| /(A psi^4) ] / ( r^4 B psi^4)


     if (cproca_l/=0) then

        rho(l,:) = rho(l,:) + (0.125d0/smallpi)*dble(cproca_l*(cproca_l+1)) &
                 *((cprocaXi_R(l,:)**2 + cprocaXi_I(l,:)**2 + cproca_mass**2*(cprocaB_R(l,:)**2 &
                 + cprocaB_I(l,:)**2))/r(l,:)**2 &
                 + r(l,:)**2*(cprocaL_R(l,:)**2 + cprocaL_I(l,:)**2)/(A(l,:)*psi4(l,:))) &
                 /(B(l,:)*psi4(l,:))

        JA(l,:)  = JA(l,:) - (0.25d0/smallpi)*dble(cproca_l*(cproca_l+1)) &
                 *(cprocaXi_R(l,:)*cprocaL_R(l,:) + cprocaXi_I(l,:)*cprocaL_I(l,:)) &
                 /(B(l,:)*psi4(l,:))

        SAA(l,:) = SAA(l,:) + (0.125d0/smallpi)*dble(cproca_l*(cproca_l+1)) &
                 *((cprocaXi_R(l,:)**2 + cprocaXi_I(l,:)**2 - cproca_mass**2*(cprocaB_R(l,:)**2 &
                 + cprocaB_I(l,:)**2))/r(l,:)**2 &
                 + r(l,:)**2*(cprocaL_R(l,:)**2 + cprocaL_I(l,:)**2)/(A(l,:)*psi4(l,:))) &
                 /(B(l,:)*psi4(l,:))

        if (.not.nolambda) then
            SLL(l,:) = SLL(l,:) + (0.125d0/smallpi)*dble(cproca_l*(cproca_l+1)) &
                     *((cprocaXi_R(l,:)**2 + cprocaXi_I(l,:)**2 - cproca_mass**2*(cprocaB_R(l,:)**2 &
                     + cprocaB_I(l,:)**2))/r(l,:)**2 &
                     + r(l,:)**2*(cprocaL_R(l,:)**2 + cprocaL_I(l,:)**2)/(A(l,:)*psi4(l,:))) &
                     /(B(l,:)*psi4(l,:))
        end if

     end if 

!    Proca charge density and current. These are calculated here and
!    not in the analysis routine since for charged fields we need the
!    current for the source of the electric field.
!
!    For a complex Proca field there is a conserved Noether current given by:
!
!                        *      nu            * nu
!    j    =  i/(8 pi) [ W      X   -  W      X    ]
!     mu                 mu nu         mu nu
!
!    where X_mu is the potential 1-form, W_{mu nu} is the field (Faraday)
!    tensor, and where * denotes complex conjugate. The (8*pi) factor is
!    the one consistent with our normalization.
!
!    From this we can define a "Proca" charge density as:
!
!                   mu
!    rho  :=  - n  j   =  [ cprocaA cprocaE  -  cprocaA cprocaE ] / (4 pi)
!       Q        mu                R       I           I       R
!
!    The radial flux is defined with the projection operator 
!
!     r          r   mu
!    j  :=  gamma   j  =  [ cprocaPhi cprocaE  -  cprocaPhi cprocaE ] / (4 pi)
!      Q          mu                  R       I             I       R

     cproca_Qdens(l,:) = (cprocaA_R(l,:)  *cprocaE_I(l,:) - cprocaA_I(l,:)  *cprocaE_R(l,:))/(4.d0*smallpi)
     cproca_Qflux(l,:) = (cprocaPhi_R(l,:)*cprocaE_I(l,:) - cprocaPhi_I(l,:)*cprocaE_R(l,:))/(4.d0*smallpi)

!    Additional terms for angular momentum.

     if (cproca_l/=0) then

!         rho  = rho  +  1/(4 pi) l(l+1)[ cprocaB cprocaXi  -  cprocaB cprocaXi  ] / r^2
!            Q      Q                            R        I           I        R

          cproca_Qdens(l,:) = cproca_Qdens(l,:) + dble(cproca_l*(cproca_l+1))/(4.d0*smallpi)/r(l,:)**2 &
                            *(cprocaB_R(l,:)*cprocaXi_I(l,:) - cprocaB_I(l,:)*cprocaXi_R(l,:))

!          r    r
!         j  = j  -  1/(4 pi) l(l+1)[ cprocaB ( cprocaA - cprocaG )  -  cprocaB ( cprocaA - cprocaG ) ] / (A B psi^8 r^2)
!                                            R         I         I             I         R         R

          cproca_Qflux(l,:) = cproca_Qdens(l,:) - dble(cproca_l*(cproca_l+1)) &
                            /(4.d0*smallpi)/(A(l,:)*B(l,:)*psi4(l,:)**2)/r(l,:)**2 &
                            *(cprocaB_R(l,:)*(cprocaA_I(l,:) - cprocaG_I(l,:)) &
                            - cprocaB_I(l,:)*(cprocaA_R(l,:) - cprocaG_R(l,:)))

     end if

!    Charge and current density. Notice that here I multiply with -q.
!    The reason for the negative sign is that we used an opposite
!    convention of the charge sign in the definition of the gauge
!    covariant derivative.

     if (contains(mattertype,"electric")) then
        echarge(l,:)  = echarge(l,:)  - cproca_q*cproca_Qdens(l,:)
        ecurrent(l,:) = ecurrent(l,:) - cproca_q*cproca_Qflux(l,:)
     end if

  end if


! ***********************
! ***   DIRAC FIELD   ***
! ***********************

! Stress-energy quantities for the Dirac field:
!
! rho  =  1/(2 pi) [ FI PiFR - FR PiFI + GI PiGR - GR PiGI ]
!
!
! JA   =  1/(4 pi) [ ( FR dFI/dr - FI dFR/dr + GR dGI/dr - GI dGR/r )
!
!      -  sqrt(A) psi**2 ( FR PiGI - FI PiGR + GR PiFI - GI PiFR ) ]
!
!
! SAA  =  1/(2 pi sqrt(A) psi**2) [ FR dGI/dr - FI dGR/dr + GR dFI/dr - GI dFR/r ]
!
!
! SBB  =  1/(2 pi r sqrt(B) psi**2) [ FR GI - FI GR ]
!
! Notice that in order to calculate (rho,JA) we need to know before the
! values of PiF and PiG, which are calculated in auxiliary_matter.f90
!
! Alternative expressions for rho and JA obtained by substituting the
! Dirac equation are:
!
! rho  =  1/(2 pi) [ m ( FR**2 + FI**2 - GR**2 - GI**2 )
!
!      +  2/(r sqrt(B) psi**2) ( FR GI - FI GR )
!
!      +  1/(sqrt(A) psi**2) ( FR dGI/dr - FI dGR/dr + GR dFI/dr - GI dFR/dr ) ]
!
!
! JA   =  1/(2 pi) [ FR dFI/dr - FI dFR/dr + GR dGI/dr - GI dGR/r ]

  if (contains(mattertype,"dirac")) then

!    This flag is there to change the way rho and JA are calculated.
!    The are two different ways explained in the comment above.
!    Both ways should work fine, but it is a good idea to check.

     usediracPi = .false.

!    1/(2*pi).

     aux = 1.d0/(2.d0*smallpi)

!    Energy density.

     if (usediracPi) then
        rho(l,:) = rho(l,:) + aux &
                 *(dirac_FI(l,:)*dirac_PiFR(l,:) - dirac_FR(l,:)*dirac_PiFI(l,:) &
                 + dirac_GI(l,:)*dirac_PiGR(l,:) - dirac_GR(l,:)*dirac_PiGI(l,:))
     else
        rho(l,:) = rho(l,:) + aux &
                 *(dirac_mass*(dirac_FR(l,:)**2 + dirac_FI(l,:)**2 - dirac_GR(l,:)**2 - dirac_GI(l,:)**2) &
                 + 2.d0/(r(l,:)*sqrt(B(l,:))*psi2(l,:))*(dirac_FR(l,:)*dirac_GI(l,:) - dirac_FI(l,:)*dirac_GR(l,:)) &
                 +(dirac_FR(l,:)*D1_dirac_GI(l,:) - dirac_FI(l,:)*D1_dirac_GR(l,:) &
                 + dirac_GR(l,:)*D1_dirac_FI(l,:) - dirac_GI(l,:)*D1_dirac_FR(l,:))/(sqrt(A(l,:))*psi2(l,:)))
     end if

!    Radial momentum density (index down).

     if (usediracPi) then
        JA(l,:) = JA(l,:) + 0.5d0*aux &
                *(dirac_FR(l,:)*D1_dirac_FI(l,:) - dirac_FI(l,:)*D1_dirac_FR(l,:) &
                + dirac_GR(l,:)*D1_dirac_GI(l,:) - dirac_GI(l,:)*D1_dirac_GR(l,:) &
                - sqrt(A(l,:))*psi2(l,:) &
                *(dirac_FR(l,:)*dirac_PiGI(l,:) - dirac_FI(l,:)*dirac_PiGR(l,:) &
                + dirac_GR(l,:)*dirac_PiFI(l,:) - dirac_GI(l,:)*dirac_PiFR(l,:)))
     else
        JA(l,:) = JA(l,:) + aux &
                *(dirac_FR(l,:)*D1_dirac_FI(l,:) - dirac_FI(l,:)*D1_dirac_FR(l,:) &
                + dirac_GR(l,:)*D1_dirac_GI(l,:) - dirac_GI(l,:)*D1_dirac_GR(l,:))
     end if

!    Stress tensor.

     SAA(l,:) = SAA(l,:) + aux/(sqrt(A(l,:))*psi2(l,:)) &
              *(dirac_FR(l,:)*D1_dirac_GI(l,:) - dirac_FI(l,:)*D1_dirac_GR(l,:) &
              + dirac_GR(l,:)*D1_dirac_FI(l,:) - dirac_GI(l,:)*D1_dirac_FR(l,:))

     SBB(l,:) = SBB(l,:) + aux/(r(l,:)*sqrt(B(l,:))*psi2(l,:)) &
              *(dirac_FR(l,:)*dirac_GI(l,:) - dirac_FI(l,:)*dirac_GR(l,:))

!    SLL = (SAA - SBB)/r**2.  This term is somewhat complicated, and in order
!    to make it regular it must be written in terms of the auxiliary variables
!    H = G/r, their derivatives, and lambda.  We find after some algebra:
!
!    SLL  =  1/(2 pi sqrt(A) psi**2) [ (FR dHI/dr - FI dHR/dr + HR dF//dr - HI dFR/dr) / r
!
!         + lambda / (1 + sqrt(A/B)) ( FR HI - HR FI ) ]

     if (.not.nolambda) then
        SLL(l,:) = SLL(l,:) + aux/(sqrt(A(l,:))*psi2(l,:)) &
                 *((dirac_FR(l,:)*D1_dirac_HI(l,:) - dirac_FI(l,:)*D1_dirac_HR(l,:) &
                 +  dirac_HR(l,:)*D1_dirac_FI(l,:) - dirac_HI(l,:)*D1_dirac_FR(l,:))/r(l,:) &
                 + lambda(l,:)/(1.d0 + sqrt(A(l,:)/B(l,:)))*(dirac_FR(l,:)*dirac_HI(l,:) - dirac_HR(l,:)*dirac_FI(l,:)))
     end if

!    Dirac particle density:
!                             2       2                     2    2    2    2
!    dens  =  1 / (2 pi) [ |F|  +  |G| ]  =  1 / (2 pi) [ FR + FI + GR + GI ]

     dirac_dens(l,:) = aux*(dirac_FR(l,:)**2 + dirac_FI(l,:)**2 + dirac_GR(l,:)**2 + dirac_GI(l,:)**2)

!    Dirac particle flux (index up):
!                                  2       *       *                          2
!    flux  =  1 / (2 pi sqrt(A) psi ) [ F G  +  G F  ]  =  1 / (pi sqrt(A) psi ) [ FR GR + FI GI ]

     dirac_flux(l,:) = 2.d0*aux/(sqrt(A(l,:))*psi2(l,:))*(dirac_FR(l,:)*dirac_GR(l,:) + dirac_FI(l,:)*dirac_GI(l,:))

!    Corrections for charged Dirac fields.

     if (contains(mattertype,"electric")) then

!       Energy density. The term that must be added depends on how we wrote rho
!       above.  If we did not substitute the Dirac equation we must add:
!
!       - q ePhi dirac_dens
!
!       but if we did substitute it we must add instead:
!
!       - q eAr dirac_flux
!
!       This is because the Dirac equation itself has terms that depend on q.
!
!       Notice that both terms above are scalars.  The first one is obvious,
!       and the second one corresponds to the contraction of the potential
!       eAr (index down) with dirac_flux j^i (index up).

        if (usediracPi) then
           rho(l,:) = rho(l,:) - dirac_q*ePhi(l,:)*dirac_dens(l,:)
        else
           rho(l,:) = rho(l,:) - dirac_q*eAr(l,:)*dirac_flux(l,:)
        end if

!       Momentum density (index down).  Again, the term we must add depends
!       on how we wrote JA above.  If we did not substitute the Dirac
!       equation we must add:
!
!       - q/2 ( eAr dirac_dens + ePhi dirac_flux A psi4 )
!
!       but if we did substitute it we must add instead:
!
!       - q eAr dirac_dens
!
!       Notice that in both cases we add a 1-form since JA is written
!       with the index down, and so is eAr, with dirac_flux has the
!       index up but the metric terms lower it.

        if (usediracPi) then
           JA(l,:) = JA(l,:) - 0.5d0*dirac_q*(eAr(l,:)*dirac_dens(l,:) + ephi(l,:)*dirac_flux(l,:)*A(l,:)*psi4(l,:))
        else
           JA(l,:) = JA(l,:) - dirac_q*eAr(l,:)*dirac_dens(l,:)
        end if

!       Stress tensor.  For the radial component we must add the term: - q*eAr*dirac_flux
!       The angular component is unchanged, but we must correct SLL.

        SAA(l,:) = SAA(l,:) - dirac_q*eAr(l,:)*dirac_flux(l,:)

        if (.not.nolambda) then
           SLL(l,:) = SLL(l,:) - dirac_q*eAr(l,:)*dirac_flux(l,:)/r(l,:)**2
        end if

!       Charge and current density (index up).

        echarge(l,:)  = echarge(l,:)  + dirac_q*dirac_dens(l,:)
        ecurrent(l,:) = ecurrent(l,:) + dirac_q*dirac_flux(l,:)

     end if

  end if


! ****************
! ***   DUST   ***
! ****************

! This is basically a fluid with zero pressure and
! zero internal energy. The stress-energy tensor is:
!
! T       =  rho0  u  u
!  mu nu            mu nu
!
! with rho0 the rest-mass density in the the fluid's
! reference frame and u^mu the fluid's 4-velocity:
!
!  mu
! u    =  ( W , W v)
!
! where v is the fluid's 3-velocity and W is the
! Lorentz factor:
!
!                             4  2
! W  =  1 / sqrt[ 1  -  A  psi  v  ]
!
! From this we find:
!
!               2
! rho  =  rho0 W  =  E + D
!
!               2          4
! JA   =  rho0 W  v ( A psi )  =  S
!
!               2  2        4
! SAA  =  rho0 W  v  ( A psi ) =  v S
!
!
! SBB  =  0

  if (contains(mattertype,"dust")) then

!    Remember not to take into account the artificial atmosphere.

     do i=1-ghost,Nr
        if (dust_cD(l,i)>dust_atmos) then

!          Energy density.

           rho(l,i) = rho(l,i) + dust_cE(l,i) + dust_cD(l,i)

!          Momentum density.

           JA(l,i) = JA(l,i) + dust_cS(l,i)

!          Stress tensor.

           SAA(l,i) = SAA(l,i) + dust_v(l,i)*dust_cS(l,i)
           SBB(l,i) = SBB(l,i)

!          SLL = (SAA - SBB)/r**2.

           if (.not.nolambda) then
              SLL(l,i) = SLL(l,i) + dust_v(l,i)*dust_cS(l,i)/r(l,i)**2
           end if

        end if
     end do

  end if


! *************************
! ***   PERFECT FLUID   ***
! *************************

! The stress energy tensor for a perfect fluid is:
!
! T       =  rho0  u  u     +  g      p
!  mu nu            mu nu       mu nu
!
! with rho0 the rest-mass density in the the fluid's
! reference frame, p the pressure, and u^mu the fluid's
! 4-velocity:
!
!  mu
! u    =  ( W , W v)
!
! where v is the fluid's 3-velocity and W is the
! Lorentz factor:
!
!                             4  2
! W  =  1 / sqrt[ 1  -  A  psi  v  ]
!
! From this we find:
!
!                 2
! rho  =  rho0 h W  -  p  =  E  +  D
!
!                 2          4
! JA   =  rho0 h W  v ( A psi )  =  S
!
!                 2  2        4
! SAA  =  rho0 h W  v  ( A psi )  +  p  =  v S  +  p
!
!
! SBB  =  p
!
! Notice that when we use artificial viscosity we must add
! the corresponding contribution to the pressure "fluid_q".

  if (contains(mattertype,"fluid")) then

!    Remember not to take into account the artificial atmosphere.

     do i=1-ghost,Nr
        if (fluid_cD(l,i)>fluid_atmos) then

!          Energy density.

           fluid_rhotot(l,i) = fluid_cE(l,i) + fluid_cD(l,i)

           rho(l,i) = rho(l,i) + fluid_rhotot(l,i)

!          Momentum density.

           JA(l,i) = JA(l,i) + fluid_cS(l,i)

!          Stress tensor.

           SAA(l,i) = SAA(l,i) + fluid_v(l,i)*fluid_cS(l,i) + fluid_p(l,i) + fluid_q(l,i)
           SBB(l,i) = SBB(l,i) + fluid_p(l,i) + fluid_q(l,i)

!          SLL = (SAA - SBB)/r**2.

           if (.not.nolambda) then
              SLL(l,i) = SLL(l,i) + fluid_v(l,i)*fluid_cS(l,i)/r(l,i)**2
           end if

        end if
     end do

  end if


! ************************
! ***   NONMIN FIELD   ***
! ************************

! VERY IMPORTANT:  THE NONMIN FIELD MUST BE AT THE END OF THE ROUTINE SINCE
! IT NEEDS THE FULL MATTER TERMS FROM ALL OTHER FIELDS.
!
! PLEASE PLACE ALL NEW MATTER FIELDS ABOVE THIS!

! The effective stress-energy for the nonmin field has two contributions,
! one identical to the scalar field, and one coming from the function
! nonmin_f.  For a constant f (fp=fpp=0), everything should reduce to the
! expressions for a real scalar field.
!
! The final expressions are (the superindex "Phi" means the standard
! expression for a scalar field):
!
!            Phi          2       4       /                            
! rho  =  rho    +  fpp xi / A psi  +  fp | pi trK
!                                         \
!
!                                                                      4  \
!      + ( d xi  -   xi ( d A / (2A) - d B / B - 2 d phi - 2/r )) / A psi |
!           r              r            r           r                     /
!
!           Phi                   r
! JA   =  JA    -  fp ( d pi  +  K  xi )  -  fpp pi xi
!                        r        r
!
!
!            Phi         2       /
! SAA  =  SAA   +  fpp pi  +  fp | pi (KTA + trK/3) - RHS
!                                \
!
!                                                      4 \
!      + ( d xi  -  xi ( d A / (2A) + 2 d phi ) / A psi  |
!           r             r              r               /
!
!
!            Phi           2     2       4         /
! SBB  =  SBB   +  fpp ( pi  - xi / A psi )  +  fp |  pi (KTB + trK/3) - RHS
!                                                  \
!
!                                                  4 \
!      +  xi ( d B / (2B) + 2 d phi + 1/r ) / A psi  |
!               r              r                     /
!
! where the quantities with superindex "Phi" are the standard expressions for
! a scalar field, and where RHS is a function of the nonmininal coupling f(phi)
! that has the form:
!
!         /             2        \ -1
! RHS  =  | f ( 1 + 3 fp / 2 f ) |
!         \                      /
!
!         /                                        2        4      2
!         | f VP - 2 fp V - (fp/2) (1 + 3 fpp) ( xi / (A psi ) - pi )
!         \
!
!                                             \
!         + (fp/2) (SAA + 2 SBB - rho)        |
!                                     Matter  /
!
! with (SAA + 2*SBB - rho) the trace of the stress-energy tensor of ALL
! matter fields OTHER than the non-minimally coupled field itself.
! (See the paper by M. Salgado: Class. Quant. Grav. 23, 4719-4741, 2006.)

  if (contains(mattertype,"nonmin")) then

!    Right hand side of field equation: Box(phi) = RHS.
!    Notice that here I need to use the values of (rho,SAA,SBB)
!    of ALL matter fields OTHER than the nonmin field itself.
!    This is the reason why the nonmin field must be at the
!    end of the routine.

     !$OMP PARALLEL DO SCHEDULE(GUIDED)
     do i=1-ghost,Nrmax

     nonmin_rhs(l,i) = (nonmin_f(l,i)*nonmin_VP(l,i) - 2.d0*nonmin_fp(l,i)*nonmin_V(l,i) &
          - 0.5d0*nonmin_fp(l,i)*(1.d0 + 3.d0*nonmin_fpp(l,i))*(nonmin_xi(l,i)**2/(A(l,i)*psi4(l,i)) &
          - nonmin_pi(l,i)**2) + 0.5d0*nonmin_fp(l,i)*(SAA(l,i) + two*SBB(l,i) - rho(l,i))) &
          /(nonmin_f(l,i)*(1.d0 + 1.5d0*nonmin_fp(l,i)**2/nonmin_f(l,i)))

!    Energy density.

     rho(l,i) = rho(l,i) + half*(nonmin_pi(l,i)**2 + nonmin_xi(l,i)**2*(1.d0 + two*nonmin_fpp(l,i))/(A(l,i)*psi4(l,i))) &
         + nonmin_fp(l,i)*((D1_nonmin_xi(l,i) - nonmin_xi(l,i)*(half*D1_A(l,i)/A(l,i) &
         - D1_B(l,i)/B(l,i) - two*D1_phi(l,i) - two/r(l,i)))/(A(l,i)*psi4(l,i)) + nonmin_pi(l,i)*trK(l,i)) &
         + nonmin_V(l,i)

!    Momentum density (lower index).

     JA(l,i) = JA(l,i) - nonmin_pi(l,i)*nonmin_xi(l,i)*(1.d0 + nonmin_fpp(l,i)) &
             - nonmin_fp(l,i)*(D1_nonmin_pi(l,i) + (KTA(l,i) + third*trK(l,i))*nonmin_xi(l,i))

!    Stress tensor.

     SAA(l,i) = SAA(l,i) + half*(nonmin_pi(l,i)**2 + nonmin_xi(l,i)**2/(A(l,i)*psi4(l,i))) &
         + nonmin_fp(l,i)*(nonmin_pi(l,i)*(KTA(l,i) + third*trK(l,i)) &
         + (D1_nonmin_xi(l,i) - nonmin_xi(l,i)*(half*D1_A(l,i)/A(l,i) + two*D1_phi(l,i)))/(A(l,i)*psi4(l,i)) &
         - nonmin_rhs(l,i)) + nonmin_fpp(l,i)*nonmin_pi(l,i)**2 &
         - nonmin_V(l,i)

     SBB(l,i) = SBB(l,i) + half*(nonmin_pi(l,i)**2 - nonmin_xi(l,i)**2*(1.d0 + two*nonmin_fpp(l,i))/(A(l,i)*psi4(l,i))) &
         + nonmin_fp(l,i)*(nonmin_pi(l,i)*(KTB(l,i) + third*trK(l,i)) &
         + nonmin_xi(l,i)*(half*D1_B(l,i)/B(l,i) + two*D1_phi(l,i) + 1.d0/r(l,i))/(A(l,i)*psi4(l,i)) &
         - nonmin_rhs(l,i)) + nonmin_fpp(l,i)*nonmin_pi(l,i)**2 &
         - nonmin_V(l,i)

!    SLL = (SAA - SBB)/r**2.

     if (.not.nolambda) then
        SLL(l,i) = SLL(l,i) + nonmin_xi(l,i)**2*(1.d0 + nonmin_fpp(l,i))/(A(l,i)*psi4(l,i))/r(l,i)**2 &
            + nonmin_fp(l,i)*(nonmin_pi(l,i)*Klambda(l,i) + (DD_nonmin_xir(l,i)/r(l,i) &
            - nonmin_xi(l,i)/r(l,i)**2*(half*(D1_A(l,i)/A(l,i) + D1_B(l,i)/B(l,i)) + 4.d0*D1_phi(l,i)))/(A(l,i)*psi4(l,i)))
     end if

     end do
     !$OMP END PARALLEL DO

  end if


! **************************************
! ***   NO MORE MATTER FIELDS HERE   ***
! **************************************

! STOP:  DON'T ADD NEW MATTER FIELDS HERE, THE NONMIN FIELD MUST BE THE LAST.
! PLEASE PLACE ALL NEW MATTER FIELDS ABOVE THE NONMIN FIELD!


! **********************************
! ***   TRACE OF STRESS TENSOR   ***
! **********************************

! After the components of the stress tensor have been calculated
! we just sum over the components.

  trS(l,:) = SAA(l,:) + 2.d0*SBB(l,:)


! ***************
! ***   END   ***
! ***************

  end subroutine stressenergy

