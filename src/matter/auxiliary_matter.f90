
  subroutine auxiliary_matter(l)

! ******************************************
! ***   AUXILIARY VARIABLES FOR MATTER   ***
! ******************************************

! Include modules.

  use param
  use arrays
  use derivatives
  use derivadvect

! Extra variables.

  implicit none

  logical contains

  integer i,l,sym


! *****************************
! ***   REAL SCALAR FIELD   ***
! *****************************

  if (contains(mattertype,"scalar")) then

!    Derivatives of phi.

     diffvar => scalar_phi
     D1_scalar_phi(l,:) = diff1(l,+1)
     D2_scalar_phi(l,:) = diff2(l,+1)

!    Derivatives of xi.

     diffvar => scalar_xi
     D1_scalar_xi(l,:) = diff1(l,-1)

!    Derivatives of pi.

     diffvar => scalar_pi
     D1_scalar_pi(l,:) = diff1(l,+1)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of phi.

        diffvar => scalar_phi
        DA_scalar_phi(l,:) = diffadv(l,+1)

!       Derivartives of xi.

        diffvar => scalar_xi
        DA_scalar_xi(l,:) = diffadv(l,-1)

!       Derivatives of pi.

        diffvar => scalar_pi
        DA_scalar_pi(l,:) = diffadv(l,+1)

     end if

  end if


! ******************************
! ***   GHOST SCALAR FIELD   ***
! ******************************

  if (contains(mattertype,"ghost")) then

!    Derivatives of phi.

     diffvar => ghost_phi
     D1_ghost_phi(l,:) = diff1(l,+1)
     D2_ghost_phi(l,:) = diff2(l,+1)

!    Derivatives of xi.

     diffvar => ghost_xi
     D1_ghost_xi(l,:) = diff1(l,-1)

!    Derivatives of pi.

     diffvar => ghost_pi
     D1_ghost_pi(l,:) = diff1(l,+1)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of phi.

        diffvar => ghost_phi
        DA_ghost_phi(l,:) = diffadv(l,+1)

!       Derivartives of xi.

        diffvar => ghost_xi
        DA_ghost_xi(l,:) = diffadv(l,-1)

!       Derivatives of pi.

        diffvar => ghost_pi
        DA_ghost_pi(l,:) = diffadv(l,+1)

     end if

  end if


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

  if (contains(mattertype,"complex")) then

!    The symmetry depends on the value of complex_l.

     sym = int(1.d0-2.d0*mod(complex_l,2))

!    Derivatives of phiR.

     diffvar => complex_phiR
     D1_complex_phiR(l,:) = diff1(l,+sym)
     D2_complex_phiR(l,:) = diff2(l,+sym)

!    Derivatives of phiI.

     diffvar => complex_phiI
     D1_complex_phiI(l,:) = diff1(l,+sym)
     D2_complex_phiI(l,:) = diff2(l,+sym)

!    Derivatives of xiR.

     diffvar => complex_xiR
     D1_complex_xiR(l,:) = diff1(l,-sym)

!    Derivatives of xiI

     diffvar => complex_xiI
     D1_complex_xiI(l,:) = diff1(l,-sym)

!    Derivatives of piR.

     diffvar => complex_piR
     D1_complex_piR(l,:) = diff1(l,+sym)

!    Derivatives of piI.

     diffvar => complex_piI
     D1_complex_piI(l,:) = diff1(l,+sym)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of phi.

        diffvar => complex_phiR
        DA_complex_phiR(l,:) = diffadv(l,+sym)

        diffvar => complex_phiI
        DA_complex_phiI(l,:) = diffadv(l,+sym)

!       Derivatives of xi.

        diffvar => complex_xiR
        DA_complex_xiR(l,:) = diffadv(l,-sym)

        diffvar => complex_xiI
        DA_complex_xiI(l,:) = diffadv(l,-sym)

!       Derivatives of pi.

        diffvar => complex_piR
        DA_complex_piR(l,:) = diffadv(l,+sym)

        diffvar => complex_piI
        DA_complex_piI(l,:) = diffadv(l,+sym)

     end if

  end if


! **************************************
! ***   COMPLEX GHOST SCALAR FIELD   ***
! **************************************

  if (contains(mattertype,"complexghost")) then

!    The symmetry depends on the value of complex_l.

     sym = int(1.d0-2.d0*mod(complexghost_l,2))

!    Derivatives of phiR.

     diffvar => complexghost_phiR
     D1_complexghost_phiR(l,:) = diff1(l,+sym)
     D2_complexghost_phiR(l,:) = diff2(l,+sym)

!    Derivatives of phiI.

     diffvar => complexghost_phiI
     D1_complexghost_phiI(l,:) = diff1(l,+sym)
     D2_complexghost_phiI(l,:) = diff2(l,+sym)

!    Derivatives of xiR.

     diffvar => complexghost_xiR
     D1_complexghost_xiR(l,:) = diff1(l,-sym)

!    Derivatives of xiI

     diffvar => complexghost_xiI
     D1_complexghost_xiI(l,:) = diff1(l,-sym)

!    Derivatives of piR.

     diffvar => complexghost_piR
     D1_complexghost_piR(l,:) = diff1(l,+sym)

!    Derivatives of piI.

     diffvar => complexghost_piI
     D1_complexghost_piI(l,:) = diff1(l,+sym)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of phi.

        diffvar => complexghost_phiR
        DA_complexghost_phiR(l,:) = diffadv(l,+sym)

        diffvar => complexghost_phiI
        DA_complexghost_phiI(l,:) = diffadv(l,+sym)

!       Derivatives of xi.

        diffvar => complexghost_xiR
        DA_complexghost_xiR(l,:) = diffadv(l,-sym)

        diffvar => complexghost_xiI
        DA_complexghost_xiI(l,:) = diffadv(l,-sym)

!       Derivatives of pi.

        diffvar => complexghost_piR
        DA_complexghost_piR(l,:) = diffadv(l,+sym)

        diffvar => complexghost_piI
        DA_complexghost_piI(l,:) = diffadv(l,+sym)

     end if

  end if


! ************************
! ***   NONMIN FIELD   ***
! ************************

  if (contains(mattertype,"nonmin")) then

!    Derivatives of phi.

     diffvar => nonmin_phi
     D1_nonmin_phi(l,:) = diff1(l,+1)
     D2_nonmin_phi(l,:) = diff2(l,+1)

!    Derivatives of xi.

     diffvar => nonmin_xi
     D1_nonmin_xi(l,:) = diff1(l,-1)

!    Derivatives of pi.

     diffvar => nonmin_pi
     D1_nonmin_pi(l,:) = diff1(l,+1)

!    Regularization term.

     auxarray(l,:) = nonmin_xi(l,:)/r(l,:)
     diffvar => auxarray
     DD_nonmin_xir(l,:) = diff1(l,+1)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of phi.

        diffvar => nonmin_phi
        DA_nonmin_phi(l,:) = diffadv(l,+1)

!       Derivartives of xi.

        diffvar => nonmin_xi
        DA_nonmin_xi(l,:) = diffadv(l,-1)

!       Derivatives of pi.

        diffvar => nonmin_pi
        DA_nonmin_pi(l,:) = diffadv(l,+1)

     end if

  end if


! *********************************
! ***   ELECTROMAGNETIC FIELD   ***
! *********************************

  if (contains(mattertype,"electric")) then

!    Derivatives of electric field.

     diffvar => electric
     D1_electric(l,:) = diff1(l,-1)

!    Derivatives of ePhi.

     diffvar => ePhi
     D1_ePhi(l,:) = diff1(l,+1)

!    Derivarives of eAr.

     diffvar => eAr
     D1_eAr(l,:) = diff1(l,-1)

!    Rescaled scalar and vector potential.

     eF(l,:) = alpha(l,:)*ePhi(l,:)
     eH(l,:) = alpha(l,:)*eAr(l,:)

     diffvar => eF
     D1_eF(l,:) = diff1(l,+1)

     diffvar => eH
     D1_eH(l,:) = diff1(l,-1)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of electric field.

        diffvar => electric
        DA_electric(l,:) = diffadv(l,-1)

!       Derivatives of scalar potential ePhi.

        diffvar => ePhi
        DA_ePhi(l,:) = diffadv(l,+1)

!       Derivarives of vector potential eAr.

        diffvar => eAr
        DA_eAr(l,:) = diffadv(l,-1)

     end if

  end if


! ***********************
! ***   PROCA FIELD   ***
! ***********************

  if (contains(mattertype,"proca")) then

!    Derivatives of Proca electric field.

     diffvar => procaE
     D1_procaE(l,:) = diff1(l,-1)

!    Derivatives of Proca scalar potential.

     diffvar => procaPhi
     D1_procaPhi(l,:) = diff1(l,+1)

!    Derivatives of Proca vector potential.

     diffvar => procaA
     D1_procaA(l,:) = diff1(l,-1)

!    Rescaled scalar and vector potentials.

     procaF(l,:) = alpha(l,:)*procaPhi(l,:)
     procaH(l,:) = alpha(l,:)*procaA(l,:)

     diffvar => procaF
     D1_procaF(l,:) = diff1(l,+1)

     diffvar => procaH
     D1_procaH(l,:) = diff1(l,+1)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of Proca electric field.

        diffvar => procaE
        DA_procaE(l,:) = diffadv(l,-1)

!       Derivatives of Proca scalar potential.

        diffvar => procaPhi
        DA_procaPhi(l,:) = diffadv(l,+1)

!       Derivatives of Proca vector potential.

        diffvar => procaA
        DA_procaA(l,:) = diffadv(l,-1)

     end if

  end if


! *******************************
! ***   COMPLEX PROCA FIELD   ***
! *******************************

  if (contains(mattertype,"complexproca")) then

!    Derivatives of Proca electric field.

     diffvar => cprocaE_R
     D1_cprocaE_R(l,:) = diff1(l,-1)

     diffvar => cprocaE_I
     D1_cprocaE_I(l,:) = diff1(l,-1)

!    Derivatives of Proca scalar potential.

     diffvar => cprocaPhi_R
     D1_cprocaPhi_R(l,:) = diff1(l,+1)

     diffvar => cprocaPhi_I
     D1_cprocaPhi_I(l,:) = diff1(l,+1)

!    Derivatives of Proca vector potential.

     diffvar => cprocaA_R
     D1_cprocaA_R(l,:) = diff1(l,-1)

     diffvar => cprocaA_I
     D1_cprocaA_I(l,:) = diff1(l,-1)

!    Rescaled scalar and vector potentials.

     cprocaF_R(l,:) = alpha(l,:)*cprocaPhi_R(l,:)
     cprocaF_I(l,:) = alpha(l,:)*cprocaPhi_I(l,:)

     cprocaH_R(l,:) = alpha(l,:)*cprocaA_R(l,:)
     cprocaH_I(l,:) = alpha(l,:)*cprocaA_I(l,:)

     diffvar => cprocaF_R
     D1_cprocaF_R(l,:) = diff1(l,+1)
     diffvar => cprocaF_I
     D1_cprocaF_I(l,:) = diff1(l,+1)

     diffvar => cprocaH_R
     D1_cprocaH_R(l,:) = diff1(l,+1)
     diffvar => cprocaH_I
     D1_cprocaH_I(l,:) = diff1(l,+1)

!    Advective derivatives.

     if (shift/="none") then

!       Derivatives of Proca electric field.

        diffvar => cprocaE_R
        DA_cprocaE_R(l,:) = diffadv(l,-1)

        diffvar => cprocaE_I
        DA_cprocaE_I(l,:) = diffadv(l,-1)

!       Derivatives of Proca scalar potential.

        diffvar => cprocaPhi_R
        DA_cprocaPhi_R(l,:) = diffadv(l,+1)

        diffvar => cprocaPhi_I
        DA_cprocaPhi_I(l,:) = diffadv(l,+1)

!       Derivatives of Proca vector potential.

        diffvar => cprocaA_R
        DA_cprocaA_R(l,:) = diffadv(l,-1)

        diffvar => cprocaA_I
        DA_cprocaA_I(l,:) = diffadv(l,-1)

     end if

  end if


! ***********************
! ***   DIRAC FIELD   ***
! ***********************

  if (contains(mattertype,"dirac")) then

!    Radial derivatives of Dirac fields (F,G).

     diffvar => dirac_FR
     D1_dirac_FR(l,:) = diff1(l,+1)

     diffvar => dirac_FI
     D1_dirac_FI(l,:) = diff1(l,+1)

     diffvar => dirac_GR
     D1_dirac_GR(l,:) = diff1(l,-1)

     diffvar => dirac_GI
     D1_dirac_GI(l,:) = diff1(l,-1)

!    Auxiliary quantity H = G/r.

     dirac_HR = dirac_GR/r
     dirac_HI = dirac_GI/r

     diffvar => dirac_HR
     D1_dirac_HR(l,:) = diff1(l,+1)

     diffvar => dirac_HI
     D1_dirac_HI(l,:) = diff1(l,+1)

!    Advective derivatives.

     if (shift/="none") then

        diffvar => dirac_FR
        DA_dirac_FR(l,:) = diffadv(l,+1)

        diffvar => dirac_FI
        DA_dirac_FI(l,:) = diffadv(l,+1)

        diffvar => dirac_GR
        D1_dirac_GR(l,:) = diffadv(l,-1)

        diffvar => dirac_GI
        D1_dirac_GI(l,:) = diffadv(l,-1)

     end if

  end if


! ****************
! ***   DUST   ***
! ****************

  if (contains(mattertype,"dust")) then

     diffvar => dust_cD
     D1_dust_cD(l,:) = diff1(l,+1)

     diffvar => dust_cE
     D1_dust_cE(l,:) = diff1(l,+1)

     diffvar => dust_cS
     D1_dust_cS(l,:) = diff1(l,-1)

!    Advective derivatives.

     if (shift/="none") then

        diffvar => dust_cD
        DA_dust_cD(l,:) = diffadv(l,+1)

        diffvar => dust_cE
        DA_dust_cE(l,:) = diffadv(l,+1)

        diffvar => dust_cS
        DA_dust_cS(l,:) = diffadv(l,-1)

     end if

!    Recover primitive variables:
!
!    Fluid speed:        v    =  S / [ (E + D) A psi^4 ]
!
!    Lorentz factor:     W    =  1 / [ 1 - A psi^4 v^2 ]^(1/2)
!
!    Rest-mass density:  rho  =  D / W
!
!    Notice that in order to prevent possible divisions by
!    very small quantities when calculating v, if the density
!    is too small we just set it to a small number and set
!    E, S and v to zero.
 
     do i=1-ghost,Nrtotal
        if (dust_cD(l,i)<=dust_atmos) then
           dust_cD(l,i) = dust_atmos/10.d0
           dust_cE(l,i) = 0.d0
           dust_cS(l,i) = 0.d0
           dust_v(l,i)  = 0.d0
        else
           dust_v(l,:) = dust_cS(l,:)/(dust_cE(l,:) + dust_cD(l,:))/(A(l,:)*psi4(l,:))
        end if
     end do

     dust_W(l,:)   = 1.d0/sqrt(abs(1.d0 - A(l,:)*psi4(l,:)*dust_v(l,:)**2))
     dust_rho(l,:) = dust_cD(l,:)/dust_W(l,:)

  end if


! *************************
! ***   PERFECT FLUID   ***
! *************************

  if (contains(mattertype,"fluid")) then

     diffvar => fluid_cD
     D1_fluid_cD(l,:) = diff1(l,+1)

     diffvar => fluid_cE
     D1_fluid_cE(l,:) = diff1(l,+1)

     diffvar => fluid_cS
     D1_fluid_cS(l,:) = diff1(l,-1)

     diffvar => fluid_p
     D1_fluid_p(l,:) = diff1(l,+1)

!    Advective derivatives.

     if (shift/="none") then

        diffvar => fluid_cD
        DA_fluid_cD(l,:) = diffadv(l,+1)

        diffvar => fluid_cE
        DA_fluid_cE(l,:) = diffadv(l,+1)

        diffvar => fluid_cS
        DA_fluid_cS(l,:) = diffadv(l,-1)

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine auxiliary_matter

