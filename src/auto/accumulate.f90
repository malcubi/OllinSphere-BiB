! Automatically generated file.  Do not edit!

  subroutine accumulate(l,k,niter,w)

  use param
  use arrays

  implicit none

  logical contains

  integer l,k,niter

  real(8) w

  if (k==1) then
     alpha_a(l,:) = w*salpha(l,:)
  else if (k<niter) then
     alpha_a(l,:) = alpha_a(l,:) + w*salpha(l,:)
  else
     salpha(l,:)  = alpha_a(l,:) + w*salpha(l,:)
  end if

  if (k==1) then
     tau_origin_a(l) = w*stau_origin(l)
  else if (k<niter) then
     tau_origin_a(l) = tau_origin_a(l) + w*stau_origin(l)
  else
     stau_origin(l)  = tau_origin_a(l) + w*stau_origin(l)
  end if

  if (shift/="none") then
     if (k==1) then
        beta_a(l,:) = w*sbeta(l,:)
     else if (k<niter) then
        beta_a(l,:) = beta_a(l,:) + w*sbeta(l,:)
     else
        sbeta(l,:)  = beta_a(l,:) + w*sbeta(l,:)
     end if
  end if

  if (shift/="none") then
     if (k==1) then
        dtbeta_a(l,:) = w*sdtbeta(l,:)
     else if (k<niter) then
        dtbeta_a(l,:) = dtbeta_a(l,:) + w*sdtbeta(l,:)
     else
        sdtbeta(l,:)  = dtbeta_a(l,:) + w*sdtbeta(l,:)
     end if
  end if

  if (shift/="none") then
     if (k==1) then
        fdriver_a(l,:) = w*sfdriver(l,:)
     else if (k<niter) then
        fdriver_a(l,:) = fdriver_a(l,:) + w*sfdriver(l,:)
     else
        sfdriver(l,:)  = fdriver_a(l,:) + w*sfdriver(l,:)
     end if
  end if

  if (k==1) then
     phi_a(l,:) = w*sphi(l,:)
  else if (k<niter) then
     phi_a(l,:) = phi_a(l,:) + w*sphi(l,:)
  else
     sphi(l,:)  = phi_a(l,:) + w*sphi(l,:)
  end if

  if (k==1) then
     chi_a(l,:) = w*schi(l,:)
  else if (k<niter) then
     chi_a(l,:) = chi_a(l,:) + w*schi(l,:)
  else
     schi(l,:)  = chi_a(l,:) + w*schi(l,:)
  end if

  if (k==1) then
     A_a(l,:) = w*sA(l,:)
  else if (k<niter) then
     A_a(l,:) = A_a(l,:) + w*sA(l,:)
  else
     sA(l,:)  = A_a(l,:) + w*sA(l,:)
  end if

  if (k==1) then
     B_a(l,:) = w*sB(l,:)
  else if (k<niter) then
     B_a(l,:) = B_a(l,:) + w*sB(l,:)
  else
     sB(l,:)  = B_a(l,:) + w*sB(l,:)
  end if

  if (k==1) then
     trK_a(l,:) = w*strK(l,:)
  else if (k<niter) then
     trK_a(l,:) = trK_a(l,:) + w*strK(l,:)
  else
     strK(l,:)  = trK_a(l,:) + w*strK(l,:)
  end if

  if (k==1) then
     KTA_a(l,:) = w*sKTA(l,:)
  else if (k<niter) then
     KTA_a(l,:) = KTA_a(l,:) + w*sKTA(l,:)
  else
     sKTA(l,:)  = KTA_a(l,:) + w*sKTA(l,:)
  end if

  if (k==1) then
     Deltar_a(l,:) = w*sDeltar(l,:)
  else if (k<niter) then
     Deltar_a(l,:) = Deltar_a(l,:) + w*sDeltar(l,:)
  else
     sDeltar(l,:)  = Deltar_a(l,:) + w*sDeltar(l,:)
  end if

  if (k==1) then
     lambda_a(l,:) = w*slambda(l,:)
  else if (k<niter) then
     lambda_a(l,:) = lambda_a(l,:) + w*slambda(l,:)
  else
     slambda(l,:)  = lambda_a(l,:) + w*slambda(l,:)
  end if

  if (k==1) then
     Klambda_a(l,:) = w*sKlambda(l,:)
  else if (k<niter) then
     Klambda_a(l,:) = Klambda_a(l,:) + w*sKlambda(l,:)
  else
     sKlambda(l,:)  = Klambda_a(l,:) + w*sKlambda(l,:)
  end if

  if (k==1) then
     lambda2_a(l,:) = w*slambda2(l,:)
  else if (k<niter) then
     lambda2_a(l,:) = lambda2_a(l,:) + w*slambda2(l,:)
  else
     slambda2(l,:)  = lambda2_a(l,:) + w*slambda2(l,:)
  end if

  if (k==1) then
     Klambda2_a(l,:) = w*sKlambda2(l,:)
  else if (k<niter) then
     Klambda2_a(l,:) = Klambda2_a(l,:) + w*sKlambda2(l,:)
  else
     sKlambda2(l,:)  = Klambda2_a(l,:) + w*sKlambda2(l,:)
  end if

  if (k==1) then
     z4theta_a(l,:) = w*sz4theta(l,:)
  else if (k<niter) then
     z4theta_a(l,:) = z4theta_a(l,:) + w*sz4theta(l,:)
  else
     sz4theta(l,:)  = z4theta_a(l,:) + w*sz4theta(l,:)
  end if

  if (mattertype /= "vacuum") then
     if (k==1) then
        P_Kodama_a(l,:) = w*sP_Kodama(l,:)
     else if (k<niter) then
        P_Kodama_a(l,:) = P_Kodama_a(l,:) + w*sP_Kodama(l,:)
     else
        sP_Kodama(l,:)  = P_Kodama_a(l,:) + w*sP_Kodama(l,:)
     end if
  end if

  if (contains(mattertype,"scalar")) then
     if (k==1) then
        scalar_phi_a(l,:) = w*sscalar_phi(l,:)
     else if (k<niter) then
        scalar_phi_a(l,:) = scalar_phi_a(l,:) + w*sscalar_phi(l,:)
     else
        sscalar_phi(l,:)  = scalar_phi_a(l,:) + w*sscalar_phi(l,:)
     end if
  end if

  if (contains(mattertype,"scalar")) then
     if (k==1) then
        scalar_xi_a(l,:) = w*sscalar_xi(l,:)
     else if (k<niter) then
        scalar_xi_a(l,:) = scalar_xi_a(l,:) + w*sscalar_xi(l,:)
     else
        sscalar_xi(l,:)  = scalar_xi_a(l,:) + w*sscalar_xi(l,:)
     end if
  end if

  if (contains(mattertype,"scalar")) then
     if (k==1) then
        scalar_pi_a(l,:) = w*sscalar_pi(l,:)
     else if (k<niter) then
        scalar_pi_a(l,:) = scalar_pi_a(l,:) + w*sscalar_pi(l,:)
     else
        sscalar_pi(l,:)  = scalar_pi_a(l,:) + w*sscalar_pi(l,:)
     end if
  end if

  if (contains(mattertype,"ghost")) then
     if (k==1) then
        ghost_phi_a(l,:) = w*sghost_phi(l,:)
     else if (k<niter) then
        ghost_phi_a(l,:) = ghost_phi_a(l,:) + w*sghost_phi(l,:)
     else
        sghost_phi(l,:)  = ghost_phi_a(l,:) + w*sghost_phi(l,:)
     end if
  end if

  if (contains(mattertype,"ghost")) then
     if (k==1) then
        ghost_xi_a(l,:) = w*sghost_xi(l,:)
     else if (k<niter) then
        ghost_xi_a(l,:) = ghost_xi_a(l,:) + w*sghost_xi(l,:)
     else
        sghost_xi(l,:)  = ghost_xi_a(l,:) + w*sghost_xi(l,:)
     end if
  end if

  if (contains(mattertype,"ghost")) then
     if (k==1) then
        ghost_pi_a(l,:) = w*sghost_pi(l,:)
     else if (k<niter) then
        ghost_pi_a(l,:) = ghost_pi_a(l,:) + w*sghost_pi(l,:)
     else
        sghost_pi(l,:)  = ghost_pi_a(l,:) + w*sghost_pi(l,:)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     if (k==1) then
        complexghost_phiR_a(l,:) = w*scomplexghost_phiR(l,:)
     else if (k<niter) then
        complexghost_phiR_a(l,:) = complexghost_phiR_a(l,:) + w*scomplexghost_phiR(l,:)
     else
        scomplexghost_phiR(l,:)  = complexghost_phiR_a(l,:) + w*scomplexghost_phiR(l,:)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     if (k==1) then
        complexghost_phiI_a(l,:) = w*scomplexghost_phiI(l,:)
     else if (k<niter) then
        complexghost_phiI_a(l,:) = complexghost_phiI_a(l,:) + w*scomplexghost_phiI(l,:)
     else
        scomplexghost_phiI(l,:)  = complexghost_phiI_a(l,:) + w*scomplexghost_phiI(l,:)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     if (k==1) then
        complexghost_xiR_a(l,:) = w*scomplexghost_xiR(l,:)
     else if (k<niter) then
        complexghost_xiR_a(l,:) = complexghost_xiR_a(l,:) + w*scomplexghost_xiR(l,:)
     else
        scomplexghost_xiR(l,:)  = complexghost_xiR_a(l,:) + w*scomplexghost_xiR(l,:)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     if (k==1) then
        complexghost_xiI_a(l,:) = w*scomplexghost_xiI(l,:)
     else if (k<niter) then
        complexghost_xiI_a(l,:) = complexghost_xiI_a(l,:) + w*scomplexghost_xiI(l,:)
     else
        scomplexghost_xiI(l,:)  = complexghost_xiI_a(l,:) + w*scomplexghost_xiI(l,:)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     if (k==1) then
        complexghost_piR_a(l,:) = w*scomplexghost_piR(l,:)
     else if (k<niter) then
        complexghost_piR_a(l,:) = complexghost_piR_a(l,:) + w*scomplexghost_piR(l,:)
     else
        scomplexghost_piR(l,:)  = complexghost_piR_a(l,:) + w*scomplexghost_piR(l,:)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     if (k==1) then
        complexghost_piI_a(l,:) = w*scomplexghost_piI(l,:)
     else if (k<niter) then
        complexghost_piI_a(l,:) = complexghost_piI_a(l,:) + w*scomplexghost_piI(l,:)
     else
        scomplexghost_piI(l,:)  = complexghost_piI_a(l,:) + w*scomplexghost_piI(l,:)
     end if
  end if

  if (contains(mattertype,"complex")) then
     if (k==1) then
        complex_phiR_a(l,:) = w*scomplex_phiR(l,:)
     else if (k<niter) then
        complex_phiR_a(l,:) = complex_phiR_a(l,:) + w*scomplex_phiR(l,:)
     else
        scomplex_phiR(l,:)  = complex_phiR_a(l,:) + w*scomplex_phiR(l,:)
     end if
  end if

  if (contains(mattertype,"complex")) then
     if (k==1) then
        complex_phiI_a(l,:) = w*scomplex_phiI(l,:)
     else if (k<niter) then
        complex_phiI_a(l,:) = complex_phiI_a(l,:) + w*scomplex_phiI(l,:)
     else
        scomplex_phiI(l,:)  = complex_phiI_a(l,:) + w*scomplex_phiI(l,:)
     end if
  end if

  if (contains(mattertype,"complex")) then
     if (k==1) then
        complex_xiR_a(l,:) = w*scomplex_xiR(l,:)
     else if (k<niter) then
        complex_xiR_a(l,:) = complex_xiR_a(l,:) + w*scomplex_xiR(l,:)
     else
        scomplex_xiR(l,:)  = complex_xiR_a(l,:) + w*scomplex_xiR(l,:)
     end if
  end if

  if (contains(mattertype,"complex")) then
     if (k==1) then
        complex_xiI_a(l,:) = w*scomplex_xiI(l,:)
     else if (k<niter) then
        complex_xiI_a(l,:) = complex_xiI_a(l,:) + w*scomplex_xiI(l,:)
     else
        scomplex_xiI(l,:)  = complex_xiI_a(l,:) + w*scomplex_xiI(l,:)
     end if
  end if

  if (contains(mattertype,"complex")) then
     if (k==1) then
        complex_piR_a(l,:) = w*scomplex_piR(l,:)
     else if (k<niter) then
        complex_piR_a(l,:) = complex_piR_a(l,:) + w*scomplex_piR(l,:)
     else
        scomplex_piR(l,:)  = complex_piR_a(l,:) + w*scomplex_piR(l,:)
     end if
  end if

  if (contains(mattertype,"complex")) then
     if (k==1) then
        complex_piI_a(l,:) = w*scomplex_piI(l,:)
     else if (k<niter) then
        complex_piI_a(l,:) = complex_piI_a(l,:) + w*scomplex_piI(l,:)
     else
        scomplex_piI(l,:)  = complex_piI_a(l,:) + w*scomplex_piI(l,:)
     end if
  end if

  if (contains(mattertype,"nonmin")) then
     if (k==1) then
        nonmin_phi_a(l,:) = w*snonmin_phi(l,:)
     else if (k<niter) then
        nonmin_phi_a(l,:) = nonmin_phi_a(l,:) + w*snonmin_phi(l,:)
     else
        snonmin_phi(l,:)  = nonmin_phi_a(l,:) + w*snonmin_phi(l,:)
     end if
  end if

  if (contains(mattertype,"nonmin")) then
     if (k==1) then
        nonmin_xi_a(l,:) = w*snonmin_xi(l,:)
     else if (k<niter) then
        nonmin_xi_a(l,:) = nonmin_xi_a(l,:) + w*snonmin_xi(l,:)
     else
        snonmin_xi(l,:)  = nonmin_xi_a(l,:) + w*snonmin_xi(l,:)
     end if
  end if

  if (contains(mattertype,"nonmin")) then
     if (k==1) then
        nonmin_pi_a(l,:) = w*snonmin_pi(l,:)
     else if (k<niter) then
        nonmin_pi_a(l,:) = nonmin_pi_a(l,:) + w*snonmin_pi(l,:)
     else
        snonmin_pi(l,:)  = nonmin_pi_a(l,:) + w*snonmin_pi(l,:)
     end if
  end if

  if (contains(mattertype,"electric")) then
     if (k==1) then
        electric_a(l,:) = w*selectric(l,:)
     else if (k<niter) then
        electric_a(l,:) = electric_a(l,:) + w*selectric(l,:)
     else
        selectric(l,:)  = electric_a(l,:) + w*selectric(l,:)
     end if
  end if

  if (contains(mattertype,"electric")) then
     if (k==1) then
        ePhi_a(l,:) = w*sePhi(l,:)
     else if (k<niter) then
        ePhi_a(l,:) = ePhi_a(l,:) + w*sePhi(l,:)
     else
        sePhi(l,:)  = ePhi_a(l,:) + w*sePhi(l,:)
     end if
  end if

  if (contains(mattertype,"electric")) then
     if (k==1) then
        eAr_a(l,:) = w*seAr(l,:)
     else if (k<niter) then
        eAr_a(l,:) = eAr_a(l,:) + w*seAr(l,:)
     else
        seAr(l,:)  = eAr_a(l,:) + w*seAr(l,:)
     end if
  end if

  if (contains(mattertype,"proca")) then
     if (k==1) then
        procaE_a(l,:) = w*sprocaE(l,:)
     else if (k<niter) then
        procaE_a(l,:) = procaE_a(l,:) + w*sprocaE(l,:)
     else
        sprocaE(l,:)  = procaE_a(l,:) + w*sprocaE(l,:)
     end if
  end if

  if (contains(mattertype,"proca")) then
     if (k==1) then
        procaPhi_a(l,:) = w*sprocaPhi(l,:)
     else if (k<niter) then
        procaPhi_a(l,:) = procaPhi_a(l,:) + w*sprocaPhi(l,:)
     else
        sprocaPhi(l,:)  = procaPhi_a(l,:) + w*sprocaPhi(l,:)
     end if
  end if

  if (contains(mattertype,"proca")) then
     if (k==1) then
        procaA_a(l,:) = w*sprocaA(l,:)
     else if (k<niter) then
        procaA_a(l,:) = procaA_a(l,:) + w*sprocaA(l,:)
     else
        sprocaA(l,:)  = procaA_a(l,:) + w*sprocaA(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaE_R_a(l,:) = w*scprocaE_R(l,:)
     else if (k<niter) then
        cprocaE_R_a(l,:) = cprocaE_R_a(l,:) + w*scprocaE_R(l,:)
     else
        scprocaE_R(l,:)  = cprocaE_R_a(l,:) + w*scprocaE_R(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaE_I_a(l,:) = w*scprocaE_I(l,:)
     else if (k<niter) then
        cprocaE_I_a(l,:) = cprocaE_I_a(l,:) + w*scprocaE_I(l,:)
     else
        scprocaE_I(l,:)  = cprocaE_I_a(l,:) + w*scprocaE_I(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaXi_R_a(l,:) = w*scprocaXi_R(l,:)
     else if (k<niter) then
        cprocaXi_R_a(l,:) = cprocaXi_R_a(l,:) + w*scprocaXi_R(l,:)
     else
        scprocaXi_R(l,:)  = cprocaXi_R_a(l,:) + w*scprocaXi_R(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaXi_I_a(l,:) = w*scprocaXi_I(l,:)
     else if (k<niter) then
        cprocaXi_I_a(l,:) = cprocaXi_I_a(l,:) + w*scprocaXi_I(l,:)
     else
        scprocaXi_I(l,:)  = cprocaXi_I_a(l,:) + w*scprocaXi_I(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaPhi_R_a(l,:) = w*scprocaPhi_R(l,:)
     else if (k<niter) then
        cprocaPhi_R_a(l,:) = cprocaPhi_R_a(l,:) + w*scprocaPhi_R(l,:)
     else
        scprocaPhi_R(l,:)  = cprocaPhi_R_a(l,:) + w*scprocaPhi_R(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaPhi_I_a(l,:) = w*scprocaPhi_I(l,:)
     else if (k<niter) then
        cprocaPhi_I_a(l,:) = cprocaPhi_I_a(l,:) + w*scprocaPhi_I(l,:)
     else
        scprocaPhi_I(l,:)  = cprocaPhi_I_a(l,:) + w*scprocaPhi_I(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaA_R_a(l,:) = w*scprocaA_R(l,:)
     else if (k<niter) then
        cprocaA_R_a(l,:) = cprocaA_R_a(l,:) + w*scprocaA_R(l,:)
     else
        scprocaA_R(l,:)  = cprocaA_R_a(l,:) + w*scprocaA_R(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaA_I_a(l,:) = w*scprocaA_I(l,:)
     else if (k<niter) then
        cprocaA_I_a(l,:) = cprocaA_I_a(l,:) + w*scprocaA_I(l,:)
     else
        scprocaA_I(l,:)  = cprocaA_I_a(l,:) + w*scprocaA_I(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaB_R_a(l,:) = w*scprocaB_R(l,:)
     else if (k<niter) then
        cprocaB_R_a(l,:) = cprocaB_R_a(l,:) + w*scprocaB_R(l,:)
     else
        scprocaB_R(l,:)  = cprocaB_R_a(l,:) + w*scprocaB_R(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaB_I_a(l,:) = w*scprocaB_I(l,:)
     else if (k<niter) then
        cprocaB_I_a(l,:) = cprocaB_I_a(l,:) + w*scprocaB_I(l,:)
     else
        scprocaB_I(l,:)  = cprocaB_I_a(l,:) + w*scprocaB_I(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaG_R_a(l,:) = w*scprocaG_R(l,:)
     else if (k<niter) then
        cprocaG_R_a(l,:) = cprocaG_R_a(l,:) + w*scprocaG_R(l,:)
     else
        scprocaG_R(l,:)  = cprocaG_R_a(l,:) + w*scprocaG_R(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaG_I_a(l,:) = w*scprocaG_I(l,:)
     else if (k<niter) then
        cprocaG_I_a(l,:) = cprocaG_I_a(l,:) + w*scprocaG_I(l,:)
     else
        scprocaG_I(l,:)  = cprocaG_I_a(l,:) + w*scprocaG_I(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaL_R_a(l,:) = w*scprocaL_R(l,:)
     else if (k<niter) then
        cprocaL_R_a(l,:) = cprocaL_R_a(l,:) + w*scprocaL_R(l,:)
     else
        scprocaL_R(l,:)  = cprocaL_R_a(l,:) + w*scprocaL_R(l,:)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (k==1) then
        cprocaL_I_a(l,:) = w*scprocaL_I(l,:)
     else if (k<niter) then
        cprocaL_I_a(l,:) = cprocaL_I_a(l,:) + w*scprocaL_I(l,:)
     else
        scprocaL_I(l,:)  = cprocaL_I_a(l,:) + w*scprocaL_I(l,:)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     if (k==1) then
        dirac_FR_a(l,:) = w*sdirac_FR(l,:)
     else if (k<niter) then
        dirac_FR_a(l,:) = dirac_FR_a(l,:) + w*sdirac_FR(l,:)
     else
        sdirac_FR(l,:)  = dirac_FR_a(l,:) + w*sdirac_FR(l,:)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     if (k==1) then
        dirac_FI_a(l,:) = w*sdirac_FI(l,:)
     else if (k<niter) then
        dirac_FI_a(l,:) = dirac_FI_a(l,:) + w*sdirac_FI(l,:)
     else
        sdirac_FI(l,:)  = dirac_FI_a(l,:) + w*sdirac_FI(l,:)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     if (k==1) then
        dirac_GR_a(l,:) = w*sdirac_GR(l,:)
     else if (k<niter) then
        dirac_GR_a(l,:) = dirac_GR_a(l,:) + w*sdirac_GR(l,:)
     else
        sdirac_GR(l,:)  = dirac_GR_a(l,:) + w*sdirac_GR(l,:)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     if (k==1) then
        dirac_GI_a(l,:) = w*sdirac_GI(l,:)
     else if (k<niter) then
        dirac_GI_a(l,:) = dirac_GI_a(l,:) + w*sdirac_GI(l,:)
     else
        sdirac_GI(l,:)  = dirac_GI_a(l,:) + w*sdirac_GI(l,:)
     end if
  end if

  if (contains(mattertype,"dust")) then
     if (k==1) then
        dust_cD_a(l,:) = w*sdust_cD(l,:)
     else if (k<niter) then
        dust_cD_a(l,:) = dust_cD_a(l,:) + w*sdust_cD(l,:)
     else
        sdust_cD(l,:)  = dust_cD_a(l,:) + w*sdust_cD(l,:)
     end if
  end if

  if (contains(mattertype,"dust")) then
     if (k==1) then
        dust_cE_a(l,:) = w*sdust_cE(l,:)
     else if (k<niter) then
        dust_cE_a(l,:) = dust_cE_a(l,:) + w*sdust_cE(l,:)
     else
        sdust_cE(l,:)  = dust_cE_a(l,:) + w*sdust_cE(l,:)
     end if
  end if

  if (contains(mattertype,"dust")) then
     if (k==1) then
        dust_cS_a(l,:) = w*sdust_cS(l,:)
     else if (k<niter) then
        dust_cS_a(l,:) = dust_cS_a(l,:) + w*sdust_cS(l,:)
     else
        sdust_cS(l,:)  = dust_cS_a(l,:) + w*sdust_cS(l,:)
     end if
  end if

  if (contains(mattertype,"fluid")) then
     if (k==1) then
        fluid_cD_a(l,:) = w*sfluid_cD(l,:)
     else if (k<niter) then
        fluid_cD_a(l,:) = fluid_cD_a(l,:) + w*sfluid_cD(l,:)
     else
        sfluid_cD(l,:)  = fluid_cD_a(l,:) + w*sfluid_cD(l,:)
     end if
  end if

  if (contains(mattertype,"fluid")) then
     if (k==1) then
        fluid_cE_a(l,:) = w*sfluid_cE(l,:)
     else if (k<niter) then
        fluid_cE_a(l,:) = fluid_cE_a(l,:) + w*sfluid_cE(l,:)
     else
        sfluid_cE(l,:)  = fluid_cE_a(l,:) + w*sfluid_cE(l,:)
     end if
  end if

  if (contains(mattertype,"fluid")) then
     if (k==1) then
        fluid_cS_a(l,:) = w*sfluid_cS(l,:)
     else if (k<niter) then
        fluid_cS_a(l,:) = fluid_cS_a(l,:) + w*sfluid_cS(l,:)
     else
        sfluid_cS(l,:)  = fluid_cS_a(l,:) + w*sfluid_cS(l,:)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_tau_a(l) = w*scosmobg_tau(l)
     else if (k<niter) then
        cosmobg_tau_a(l) = cosmobg_tau_a(l) + w*scosmobg_tau(l)
     else
        scosmobg_tau(l)  = cosmobg_tau_a(l) + w*scosmobg_tau(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_H_a(l) = w*scosmobg_H(l)
     else if (k<niter) then
        cosmobg_H_a(l) = cosmobg_H_a(l) + w*scosmobg_H(l)
     else
        scosmobg_H(l)  = cosmobg_H_a(l) + w*scosmobg_H(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_a_a(l) = w*scosmobg_a(l)
     else if (k<niter) then
        cosmobg_a_a(l) = cosmobg_a_a(l) + w*scosmobg_a(l)
     else
        scosmobg_a(l)  = cosmobg_a_a(l) + w*scosmobg_a(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_phi_a(l) = w*scosmobg_phi(l)
     else if (k<niter) then
        cosmobg_phi_a(l) = cosmobg_phi_a(l) + w*scosmobg_phi(l)
     else
        scosmobg_phi(l)  = cosmobg_phi_a(l) + w*scosmobg_phi(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_trK_a(l) = w*scosmobg_trK(l)
     else if (k<niter) then
        cosmobg_trK_a(l) = cosmobg_trK_a(l) + w*scosmobg_trK(l)
     else
        scosmobg_trK(l)  = cosmobg_trK_a(l) + w*scosmobg_trK(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_alpha_a(l) = w*scosmobg_alpha(l)
     else if (k<niter) then
        cosmobg_alpha_a(l) = cosmobg_alpha_a(l) + w*scosmobg_alpha(l)
     else
        scosmobg_alpha(l)  = cosmobg_alpha_a(l) + w*scosmobg_alpha(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_scalar_phi_a(l) = w*scosmobg_scalar_phi(l)
     else if (k<niter) then
        cosmobg_scalar_phi_a(l) = cosmobg_scalar_phi_a(l) + w*scosmobg_scalar_phi(l)
     else
        scosmobg_scalar_phi(l)  = cosmobg_scalar_phi_a(l) + w*scosmobg_scalar_phi(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_scalar_pi_a(l) = w*scosmobg_scalar_pi(l)
     else if (k<niter) then
        cosmobg_scalar_pi_a(l) = cosmobg_scalar_pi_a(l) + w*scosmobg_scalar_pi(l)
     else
        scosmobg_scalar_pi(l)  = cosmobg_scalar_pi_a(l) + w*scosmobg_scalar_pi(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_complex_phiR_a(l) = w*scosmobg_complex_phiR(l)
     else if (k<niter) then
        cosmobg_complex_phiR_a(l) = cosmobg_complex_phiR_a(l) + w*scosmobg_complex_phiR(l)
     else
        scosmobg_complex_phiR(l)  = cosmobg_complex_phiR_a(l) + w*scosmobg_complex_phiR(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_complex_phiI_a(l) = w*scosmobg_complex_phiI(l)
     else if (k<niter) then
        cosmobg_complex_phiI_a(l) = cosmobg_complex_phiI_a(l) + w*scosmobg_complex_phiI(l)
     else
        scosmobg_complex_phiI(l)  = cosmobg_complex_phiI_a(l) + w*scosmobg_complex_phiI(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_complex_piR_a(l) = w*scosmobg_complex_piR(l)
     else if (k<niter) then
        cosmobg_complex_piR_a(l) = cosmobg_complex_piR_a(l) + w*scosmobg_complex_piR(l)
     else
        scosmobg_complex_piR(l)  = cosmobg_complex_piR_a(l) + w*scosmobg_complex_piR(l)
     end if
  end if

  if (cosmic_run) then
     if (k==1) then
        cosmobg_complex_piI_a(l) = w*scosmobg_complex_piI(l)
     else if (k<niter) then
        cosmobg_complex_piI_a(l) = cosmobg_complex_piI_a(l) + w*scosmobg_complex_piI(l)
     else
        scosmobg_complex_piI(l)  = cosmobg_complex_piI_a(l) + w*scosmobg_complex_piI(l)
     end if
  end if

  end subroutine accumulate

