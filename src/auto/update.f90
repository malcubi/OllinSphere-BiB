! Automatically generated file.  Do not edit!

  subroutine update(l,dtw)

  use param
  use arrays

  implicit none

  logical contains

  integer l

  real(8) dtw 

  alpha(l,:) = alpha_p(l,:) + dtw*salpha(l,:)

  tau_origin(l) = tau_origin_p(l) + dtw*stau_origin(l)

  if (shift/="none") then
     beta(l,:) = beta_p(l,:) + dtw*sbeta(l,:)
  end if

  if (shift/="none") then
     dtbeta(l,:) = dtbeta_p(l,:) + dtw*sdtbeta(l,:)
  end if

  if (shift/="none") then
     fdriver(l,:) = fdriver_p(l,:) + dtw*sfdriver(l,:)
  end if

  phi(l,:) = phi_p(l,:) + dtw*sphi(l,:)

  chi(l,:) = chi_p(l,:) + dtw*schi(l,:)

  A(l,:) = A_p(l,:) + dtw*sA(l,:)

  B(l,:) = B_p(l,:) + dtw*sB(l,:)

  trK(l,:) = trK_p(l,:) + dtw*strK(l,:)

  KTA(l,:) = KTA_p(l,:) + dtw*sKTA(l,:)

  Deltar(l,:) = Deltar_p(l,:) + dtw*sDeltar(l,:)

  lambda(l,:) = lambda_p(l,:) + dtw*slambda(l,:)

  Klambda(l,:) = Klambda_p(l,:) + dtw*sKlambda(l,:)

  lambda2(l,:) = lambda2_p(l,:) + dtw*slambda2(l,:)

  Klambda2(l,:) = Klambda2_p(l,:) + dtw*sKlambda2(l,:)

  z4theta(l,:) = z4theta_p(l,:) + dtw*sz4theta(l,:)

  if (mattertype /= "vacuum") then
     P_Kodama(l,:) = P_Kodama_p(l,:) + dtw*sP_Kodama(l,:)
  end if

  if (contains(mattertype,"scalar")) then
     scalar_phi(l,:) = scalar_phi_p(l,:) + dtw*sscalar_phi(l,:)
  end if

  if (contains(mattertype,"scalar")) then
     scalar_xi(l,:) = scalar_xi_p(l,:) + dtw*sscalar_xi(l,:)
  end if

  if (contains(mattertype,"scalar")) then
     scalar_pi(l,:) = scalar_pi_p(l,:) + dtw*sscalar_pi(l,:)
  end if

  if (contains(mattertype,"ghost")) then
     ghost_phi(l,:) = ghost_phi_p(l,:) + dtw*sghost_phi(l,:)
  end if

  if (contains(mattertype,"ghost")) then
     ghost_xi(l,:) = ghost_xi_p(l,:) + dtw*sghost_xi(l,:)
  end if

  if (contains(mattertype,"ghost")) then
     ghost_pi(l,:) = ghost_pi_p(l,:) + dtw*sghost_pi(l,:)
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_phiR(l,:) = complexghost_phiR_p(l,:) + dtw*scomplexghost_phiR(l,:)
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_phiI(l,:) = complexghost_phiI_p(l,:) + dtw*scomplexghost_phiI(l,:)
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_xiR(l,:) = complexghost_xiR_p(l,:) + dtw*scomplexghost_xiR(l,:)
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_xiI(l,:) = complexghost_xiI_p(l,:) + dtw*scomplexghost_xiI(l,:)
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_piR(l,:) = complexghost_piR_p(l,:) + dtw*scomplexghost_piR(l,:)
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_piI(l,:) = complexghost_piI_p(l,:) + dtw*scomplexghost_piI(l,:)
  end if

  if (contains(mattertype,"complex")) then
     complex_phiR(l,:) = complex_phiR_p(l,:) + dtw*scomplex_phiR(l,:)
  end if

  if (contains(mattertype,"complex")) then
     complex_phiI(l,:) = complex_phiI_p(l,:) + dtw*scomplex_phiI(l,:)
  end if

  if (contains(mattertype,"complex")) then
     complex_xiR(l,:) = complex_xiR_p(l,:) + dtw*scomplex_xiR(l,:)
  end if

  if (contains(mattertype,"complex")) then
     complex_xiI(l,:) = complex_xiI_p(l,:) + dtw*scomplex_xiI(l,:)
  end if

  if (contains(mattertype,"complex")) then
     complex_piR(l,:) = complex_piR_p(l,:) + dtw*scomplex_piR(l,:)
  end if

  if (contains(mattertype,"complex")) then
     complex_piI(l,:) = complex_piI_p(l,:) + dtw*scomplex_piI(l,:)
  end if

  if (contains(mattertype,"nonmin")) then
     nonmin_phi(l,:) = nonmin_phi_p(l,:) + dtw*snonmin_phi(l,:)
  end if

  if (contains(mattertype,"nonmin")) then
     nonmin_xi(l,:) = nonmin_xi_p(l,:) + dtw*snonmin_xi(l,:)
  end if

  if (contains(mattertype,"nonmin")) then
     nonmin_pi(l,:) = nonmin_pi_p(l,:) + dtw*snonmin_pi(l,:)
  end if

  if (contains(mattertype,"electric")) then
     electric(l,:) = electric_p(l,:) + dtw*selectric(l,:)
  end if

  if (contains(mattertype,"electric")) then
     ePhi(l,:) = ePhi_p(l,:) + dtw*sePhi(l,:)
  end if

  if (contains(mattertype,"electric")) then
     eAr(l,:) = eAr_p(l,:) + dtw*seAr(l,:)
  end if

  if (contains(mattertype,"proca")) then
     procaE(l,:) = procaE_p(l,:) + dtw*sprocaE(l,:)
  end if

  if (contains(mattertype,"proca")) then
     procaPhi(l,:) = procaPhi_p(l,:) + dtw*sprocaPhi(l,:)
  end if

  if (contains(mattertype,"proca")) then
     procaA(l,:) = procaA_p(l,:) + dtw*sprocaA(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaE_R(l,:) = cprocaE_R_p(l,:) + dtw*scprocaE_R(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaE_I(l,:) = cprocaE_I_p(l,:) + dtw*scprocaE_I(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaXi_R(l,:) = cprocaXi_R_p(l,:) + dtw*scprocaXi_R(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaXi_I(l,:) = cprocaXi_I_p(l,:) + dtw*scprocaXi_I(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaPhi_R(l,:) = cprocaPhi_R_p(l,:) + dtw*scprocaPhi_R(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaPhi_I(l,:) = cprocaPhi_I_p(l,:) + dtw*scprocaPhi_I(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaA_R(l,:) = cprocaA_R_p(l,:) + dtw*scprocaA_R(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaA_I(l,:) = cprocaA_I_p(l,:) + dtw*scprocaA_I(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaB_R(l,:) = cprocaB_R_p(l,:) + dtw*scprocaB_R(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaB_I(l,:) = cprocaB_I_p(l,:) + dtw*scprocaB_I(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaG_R(l,:) = cprocaG_R_p(l,:) + dtw*scprocaG_R(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaG_I(l,:) = cprocaG_I_p(l,:) + dtw*scprocaG_I(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaL_R(l,:) = cprocaL_R_p(l,:) + dtw*scprocaL_R(l,:)
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaL_I(l,:) = cprocaL_I_p(l,:) + dtw*scprocaL_I(l,:)
  end if

  if (contains(mattertype,"dirac")) then
     dirac_FR(l,:) = dirac_FR_p(l,:) + dtw*sdirac_FR(l,:)
  end if

  if (contains(mattertype,"dirac")) then
     dirac_FI(l,:) = dirac_FI_p(l,:) + dtw*sdirac_FI(l,:)
  end if

  if (contains(mattertype,"dirac")) then
     dirac_GR(l,:) = dirac_GR_p(l,:) + dtw*sdirac_GR(l,:)
  end if

  if (contains(mattertype,"dirac")) then
     dirac_GI(l,:) = dirac_GI_p(l,:) + dtw*sdirac_GI(l,:)
  end if

  if (contains(mattertype,"dust")) then
     dust_cD(l,:) = dust_cD_p(l,:) + dtw*sdust_cD(l,:)
  end if

  if (contains(mattertype,"dust")) then
     dust_cE(l,:) = dust_cE_p(l,:) + dtw*sdust_cE(l,:)
  end if

  if (contains(mattertype,"dust")) then
     dust_cS(l,:) = dust_cS_p(l,:) + dtw*sdust_cS(l,:)
  end if

  if (contains(mattertype,"fluid")) then
     fluid_cD(l,:) = fluid_cD_p(l,:) + dtw*sfluid_cD(l,:)
  end if

  if (contains(mattertype,"fluid")) then
     fluid_cE(l,:) = fluid_cE_p(l,:) + dtw*sfluid_cE(l,:)
  end if

  if (contains(mattertype,"fluid")) then
     fluid_cS(l,:) = fluid_cS_p(l,:) + dtw*sfluid_cS(l,:)
  end if

  if (cosmic_run) then
     cosmobg_tau(l) = cosmobg_tau_p(l) + dtw*scosmobg_tau(l)
  end if

  if (cosmic_run) then
     cosmobg_H(l) = cosmobg_H_p(l) + dtw*scosmobg_H(l)
  end if

  if (cosmic_run) then
     cosmobg_a(l) = cosmobg_a_p(l) + dtw*scosmobg_a(l)
  end if

  if (cosmic_run) then
     cosmobg_phi(l) = cosmobg_phi_p(l) + dtw*scosmobg_phi(l)
  end if

  if (cosmic_run) then
     cosmobg_trK(l) = cosmobg_trK_p(l) + dtw*scosmobg_trK(l)
  end if

  if (cosmic_run) then
     cosmobg_alpha(l) = cosmobg_alpha_p(l) + dtw*scosmobg_alpha(l)
  end if

  if (cosmic_run) then
     cosmobg_scalar_phi(l) = cosmobg_scalar_phi_p(l) + dtw*scosmobg_scalar_phi(l)
  end if

  if (cosmic_run) then
     cosmobg_scalar_pi(l) = cosmobg_scalar_pi_p(l) + dtw*scosmobg_scalar_pi(l)
  end if

  if (cosmic_run) then
     cosmobg_complex_phiR(l) = cosmobg_complex_phiR_p(l) + dtw*scosmobg_complex_phiR(l)
  end if

  if (cosmic_run) then
     cosmobg_complex_phiI(l) = cosmobg_complex_phiI_p(l) + dtw*scosmobg_complex_phiI(l)
  end if

  if (cosmic_run) then
     cosmobg_complex_piR(l) = cosmobg_complex_piR_p(l) + dtw*scosmobg_complex_piR(l)
  end if

  if (cosmic_run) then
     cosmobg_complex_piI(l) = cosmobg_complex_piI_p(l) + dtw*scosmobg_complex_piI(l)
  end if

  end subroutine update

