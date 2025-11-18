! Automatically generated file.  Do not edit!

  subroutine saveold(l)

  use param
  use arrays

  implicit none

  logical contains

  integer l,i

  alpha_p(l,:) = alpha(l,:)
  do i=0,ghost-1
     alpha_bound(l,i,3) = alpha_bound(l,i,2)
     alpha_bound(l,i,2) = alpha_bound(l,i,1)
     alpha_bound(l,i,1) = alpha(l,Nr-i)
  end do

  tau_origin_p(l) = tau_origin(l)
  if (shift/="none") then
     beta_p(l,:) = beta(l,:)
     do i=0,ghost-1
        beta_bound(l,i,3) = beta_bound(l,i,2)
        beta_bound(l,i,2) = beta_bound(l,i,1)
        beta_bound(l,i,1) = beta(l,Nr-i)
     end do
  end if

  if (shift/="none") then
     dtbeta_p(l,:) = dtbeta(l,:)
     do i=0,ghost-1
        dtbeta_bound(l,i,3) = dtbeta_bound(l,i,2)
        dtbeta_bound(l,i,2) = dtbeta_bound(l,i,1)
        dtbeta_bound(l,i,1) = dtbeta(l,Nr-i)
     end do
  end if

  if (shift/="none") then
     fdriver_p(l,:) = fdriver(l,:)
     do i=0,ghost-1
        fdriver_bound(l,i,3) = fdriver_bound(l,i,2)
        fdriver_bound(l,i,2) = fdriver_bound(l,i,1)
        fdriver_bound(l,i,1) = fdriver(l,Nr-i)
     end do
  end if

  phi_p(l,:) = phi(l,:)
  do i=0,ghost-1
     phi_bound(l,i,3) = phi_bound(l,i,2)
     phi_bound(l,i,2) = phi_bound(l,i,1)
     phi_bound(l,i,1) = phi(l,Nr-i)
  end do

  chi_p(l,:) = chi(l,:)
  do i=0,ghost-1
     chi_bound(l,i,3) = chi_bound(l,i,2)
     chi_bound(l,i,2) = chi_bound(l,i,1)
     chi_bound(l,i,1) = chi(l,Nr-i)
  end do

  A_p(l,:) = A(l,:)
  do i=0,ghost-1
     A_bound(l,i,3) = A_bound(l,i,2)
     A_bound(l,i,2) = A_bound(l,i,1)
     A_bound(l,i,1) = A(l,Nr-i)
  end do

  B_p(l,:) = B(l,:)
  do i=0,ghost-1
     B_bound(l,i,3) = B_bound(l,i,2)
     B_bound(l,i,2) = B_bound(l,i,1)
     B_bound(l,i,1) = B(l,Nr-i)
  end do

  trK_p(l,:) = trK(l,:)
  do i=0,ghost-1
     trK_bound(l,i,3) = trK_bound(l,i,2)
     trK_bound(l,i,2) = trK_bound(l,i,1)
     trK_bound(l,i,1) = trK(l,Nr-i)
  end do

  KTA_p(l,:) = KTA(l,:)
  do i=0,ghost-1
     KTA_bound(l,i,3) = KTA_bound(l,i,2)
     KTA_bound(l,i,2) = KTA_bound(l,i,1)
     KTA_bound(l,i,1) = KTA(l,Nr-i)
  end do

  Deltar_p(l,:) = Deltar(l,:)
  do i=0,ghost-1
     Deltar_bound(l,i,3) = Deltar_bound(l,i,2)
     Deltar_bound(l,i,2) = Deltar_bound(l,i,1)
     Deltar_bound(l,i,1) = Deltar(l,Nr-i)
  end do

  lambda_p(l,:) = lambda(l,:)
  do i=0,ghost-1
     lambda_bound(l,i,3) = lambda_bound(l,i,2)
     lambda_bound(l,i,2) = lambda_bound(l,i,1)
     lambda_bound(l,i,1) = lambda(l,Nr-i)
  end do

  Klambda_p(l,:) = Klambda(l,:)
  do i=0,ghost-1
     Klambda_bound(l,i,3) = Klambda_bound(l,i,2)
     Klambda_bound(l,i,2) = Klambda_bound(l,i,1)
     Klambda_bound(l,i,1) = Klambda(l,Nr-i)
  end do

  lambda2_p(l,:) = lambda2(l,:)
  do i=0,ghost-1
     lambda2_bound(l,i,3) = lambda2_bound(l,i,2)
     lambda2_bound(l,i,2) = lambda2_bound(l,i,1)
     lambda2_bound(l,i,1) = lambda2(l,Nr-i)
  end do

  Klambda2_p(l,:) = Klambda2(l,:)
  do i=0,ghost-1
     Klambda2_bound(l,i,3) = Klambda2_bound(l,i,2)
     Klambda2_bound(l,i,2) = Klambda2_bound(l,i,1)
     Klambda2_bound(l,i,1) = Klambda2(l,Nr-i)
  end do

  z4theta_p(l,:) = z4theta(l,:)
  do i=0,ghost-1
     z4theta_bound(l,i,3) = z4theta_bound(l,i,2)
     z4theta_bound(l,i,2) = z4theta_bound(l,i,1)
     z4theta_bound(l,i,1) = z4theta(l,Nr-i)
  end do

  if (mattertype /= "vacuum") then
     P_Kodama_p(l,:) = P_Kodama(l,:)
     do i=0,ghost-1
        P_Kodama_bound(l,i,3) = P_Kodama_bound(l,i,2)
        P_Kodama_bound(l,i,2) = P_Kodama_bound(l,i,1)
        P_Kodama_bound(l,i,1) = P_Kodama(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"scalar")) then
     scalar_phi_p(l,:) = scalar_phi(l,:)
     do i=0,ghost-1
        scalar_phi_bound(l,i,3) = scalar_phi_bound(l,i,2)
        scalar_phi_bound(l,i,2) = scalar_phi_bound(l,i,1)
        scalar_phi_bound(l,i,1) = scalar_phi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"scalar")) then
     scalar_xi_p(l,:) = scalar_xi(l,:)
     do i=0,ghost-1
        scalar_xi_bound(l,i,3) = scalar_xi_bound(l,i,2)
        scalar_xi_bound(l,i,2) = scalar_xi_bound(l,i,1)
        scalar_xi_bound(l,i,1) = scalar_xi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"scalar")) then
     scalar_pi_p(l,:) = scalar_pi(l,:)
     do i=0,ghost-1
        scalar_pi_bound(l,i,3) = scalar_pi_bound(l,i,2)
        scalar_pi_bound(l,i,2) = scalar_pi_bound(l,i,1)
        scalar_pi_bound(l,i,1) = scalar_pi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"ghost")) then
     ghost_phi_p(l,:) = ghost_phi(l,:)
     do i=0,ghost-1
        ghost_phi_bound(l,i,3) = ghost_phi_bound(l,i,2)
        ghost_phi_bound(l,i,2) = ghost_phi_bound(l,i,1)
        ghost_phi_bound(l,i,1) = ghost_phi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"ghost")) then
     ghost_xi_p(l,:) = ghost_xi(l,:)
     do i=0,ghost-1
        ghost_xi_bound(l,i,3) = ghost_xi_bound(l,i,2)
        ghost_xi_bound(l,i,2) = ghost_xi_bound(l,i,1)
        ghost_xi_bound(l,i,1) = ghost_xi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"ghost")) then
     ghost_pi_p(l,:) = ghost_pi(l,:)
     do i=0,ghost-1
        ghost_pi_bound(l,i,3) = ghost_pi_bound(l,i,2)
        ghost_pi_bound(l,i,2) = ghost_pi_bound(l,i,1)
        ghost_pi_bound(l,i,1) = ghost_pi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_phiR_p(l,:) = complexghost_phiR(l,:)
     do i=0,ghost-1
        complexghost_phiR_bound(l,i,3) = complexghost_phiR_bound(l,i,2)
        complexghost_phiR_bound(l,i,2) = complexghost_phiR_bound(l,i,1)
        complexghost_phiR_bound(l,i,1) = complexghost_phiR(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_phiI_p(l,:) = complexghost_phiI(l,:)
     do i=0,ghost-1
        complexghost_phiI_bound(l,i,3) = complexghost_phiI_bound(l,i,2)
        complexghost_phiI_bound(l,i,2) = complexghost_phiI_bound(l,i,1)
        complexghost_phiI_bound(l,i,1) = complexghost_phiI(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_xiR_p(l,:) = complexghost_xiR(l,:)
     do i=0,ghost-1
        complexghost_xiR_bound(l,i,3) = complexghost_xiR_bound(l,i,2)
        complexghost_xiR_bound(l,i,2) = complexghost_xiR_bound(l,i,1)
        complexghost_xiR_bound(l,i,1) = complexghost_xiR(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_xiI_p(l,:) = complexghost_xiI(l,:)
     do i=0,ghost-1
        complexghost_xiI_bound(l,i,3) = complexghost_xiI_bound(l,i,2)
        complexghost_xiI_bound(l,i,2) = complexghost_xiI_bound(l,i,1)
        complexghost_xiI_bound(l,i,1) = complexghost_xiI(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_piR_p(l,:) = complexghost_piR(l,:)
     do i=0,ghost-1
        complexghost_piR_bound(l,i,3) = complexghost_piR_bound(l,i,2)
        complexghost_piR_bound(l,i,2) = complexghost_piR_bound(l,i,1)
        complexghost_piR_bound(l,i,1) = complexghost_piR(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     complexghost_piI_p(l,:) = complexghost_piI(l,:)
     do i=0,ghost-1
        complexghost_piI_bound(l,i,3) = complexghost_piI_bound(l,i,2)
        complexghost_piI_bound(l,i,2) = complexghost_piI_bound(l,i,1)
        complexghost_piI_bound(l,i,1) = complexghost_piI(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complex")) then
     complex_phiR_p(l,:) = complex_phiR(l,:)
     do i=0,ghost-1
        complex_phiR_bound(l,i,3) = complex_phiR_bound(l,i,2)
        complex_phiR_bound(l,i,2) = complex_phiR_bound(l,i,1)
        complex_phiR_bound(l,i,1) = complex_phiR(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complex")) then
     complex_phiI_p(l,:) = complex_phiI(l,:)
     do i=0,ghost-1
        complex_phiI_bound(l,i,3) = complex_phiI_bound(l,i,2)
        complex_phiI_bound(l,i,2) = complex_phiI_bound(l,i,1)
        complex_phiI_bound(l,i,1) = complex_phiI(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complex")) then
     complex_xiR_p(l,:) = complex_xiR(l,:)
     do i=0,ghost-1
        complex_xiR_bound(l,i,3) = complex_xiR_bound(l,i,2)
        complex_xiR_bound(l,i,2) = complex_xiR_bound(l,i,1)
        complex_xiR_bound(l,i,1) = complex_xiR(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complex")) then
     complex_xiI_p(l,:) = complex_xiI(l,:)
     do i=0,ghost-1
        complex_xiI_bound(l,i,3) = complex_xiI_bound(l,i,2)
        complex_xiI_bound(l,i,2) = complex_xiI_bound(l,i,1)
        complex_xiI_bound(l,i,1) = complex_xiI(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complex")) then
     complex_piR_p(l,:) = complex_piR(l,:)
     do i=0,ghost-1
        complex_piR_bound(l,i,3) = complex_piR_bound(l,i,2)
        complex_piR_bound(l,i,2) = complex_piR_bound(l,i,1)
        complex_piR_bound(l,i,1) = complex_piR(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complex")) then
     complex_piI_p(l,:) = complex_piI(l,:)
     do i=0,ghost-1
        complex_piI_bound(l,i,3) = complex_piI_bound(l,i,2)
        complex_piI_bound(l,i,2) = complex_piI_bound(l,i,1)
        complex_piI_bound(l,i,1) = complex_piI(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"nonmin")) then
     nonmin_phi_p(l,:) = nonmin_phi(l,:)
     do i=0,ghost-1
        nonmin_phi_bound(l,i,3) = nonmin_phi_bound(l,i,2)
        nonmin_phi_bound(l,i,2) = nonmin_phi_bound(l,i,1)
        nonmin_phi_bound(l,i,1) = nonmin_phi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"nonmin")) then
     nonmin_xi_p(l,:) = nonmin_xi(l,:)
     do i=0,ghost-1
        nonmin_xi_bound(l,i,3) = nonmin_xi_bound(l,i,2)
        nonmin_xi_bound(l,i,2) = nonmin_xi_bound(l,i,1)
        nonmin_xi_bound(l,i,1) = nonmin_xi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"nonmin")) then
     nonmin_pi_p(l,:) = nonmin_pi(l,:)
     do i=0,ghost-1
        nonmin_pi_bound(l,i,3) = nonmin_pi_bound(l,i,2)
        nonmin_pi_bound(l,i,2) = nonmin_pi_bound(l,i,1)
        nonmin_pi_bound(l,i,1) = nonmin_pi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"electric")) then
     electric_p(l,:) = electric(l,:)
     do i=0,ghost-1
        electric_bound(l,i,3) = electric_bound(l,i,2)
        electric_bound(l,i,2) = electric_bound(l,i,1)
        electric_bound(l,i,1) = electric(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"electric")) then
     ePhi_p(l,:) = ePhi(l,:)
     do i=0,ghost-1
        ePhi_bound(l,i,3) = ePhi_bound(l,i,2)
        ePhi_bound(l,i,2) = ePhi_bound(l,i,1)
        ePhi_bound(l,i,1) = ePhi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"electric")) then
     eAr_p(l,:) = eAr(l,:)
     do i=0,ghost-1
        eAr_bound(l,i,3) = eAr_bound(l,i,2)
        eAr_bound(l,i,2) = eAr_bound(l,i,1)
        eAr_bound(l,i,1) = eAr(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"proca")) then
     procaE_p(l,:) = procaE(l,:)
     do i=0,ghost-1
        procaE_bound(l,i,3) = procaE_bound(l,i,2)
        procaE_bound(l,i,2) = procaE_bound(l,i,1)
        procaE_bound(l,i,1) = procaE(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"proca")) then
     procaPhi_p(l,:) = procaPhi(l,:)
     do i=0,ghost-1
        procaPhi_bound(l,i,3) = procaPhi_bound(l,i,2)
        procaPhi_bound(l,i,2) = procaPhi_bound(l,i,1)
        procaPhi_bound(l,i,1) = procaPhi(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"proca")) then
     procaA_p(l,:) = procaA(l,:)
     do i=0,ghost-1
        procaA_bound(l,i,3) = procaA_bound(l,i,2)
        procaA_bound(l,i,2) = procaA_bound(l,i,1)
        procaA_bound(l,i,1) = procaA(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaE_R_p(l,:) = cprocaE_R(l,:)
     do i=0,ghost-1
        cprocaE_R_bound(l,i,3) = cprocaE_R_bound(l,i,2)
        cprocaE_R_bound(l,i,2) = cprocaE_R_bound(l,i,1)
        cprocaE_R_bound(l,i,1) = cprocaE_R(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaE_I_p(l,:) = cprocaE_I(l,:)
     do i=0,ghost-1
        cprocaE_I_bound(l,i,3) = cprocaE_I_bound(l,i,2)
        cprocaE_I_bound(l,i,2) = cprocaE_I_bound(l,i,1)
        cprocaE_I_bound(l,i,1) = cprocaE_I(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaXi_R_p(l,:) = cprocaXi_R(l,:)
     do i=0,ghost-1
        cprocaXi_R_bound(l,i,3) = cprocaXi_R_bound(l,i,2)
        cprocaXi_R_bound(l,i,2) = cprocaXi_R_bound(l,i,1)
        cprocaXi_R_bound(l,i,1) = cprocaXi_R(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaXi_I_p(l,:) = cprocaXi_I(l,:)
     do i=0,ghost-1
        cprocaXi_I_bound(l,i,3) = cprocaXi_I_bound(l,i,2)
        cprocaXi_I_bound(l,i,2) = cprocaXi_I_bound(l,i,1)
        cprocaXi_I_bound(l,i,1) = cprocaXi_I(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaPhi_R_p(l,:) = cprocaPhi_R(l,:)
     do i=0,ghost-1
        cprocaPhi_R_bound(l,i,3) = cprocaPhi_R_bound(l,i,2)
        cprocaPhi_R_bound(l,i,2) = cprocaPhi_R_bound(l,i,1)
        cprocaPhi_R_bound(l,i,1) = cprocaPhi_R(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaPhi_I_p(l,:) = cprocaPhi_I(l,:)
     do i=0,ghost-1
        cprocaPhi_I_bound(l,i,3) = cprocaPhi_I_bound(l,i,2)
        cprocaPhi_I_bound(l,i,2) = cprocaPhi_I_bound(l,i,1)
        cprocaPhi_I_bound(l,i,1) = cprocaPhi_I(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaA_R_p(l,:) = cprocaA_R(l,:)
     do i=0,ghost-1
        cprocaA_R_bound(l,i,3) = cprocaA_R_bound(l,i,2)
        cprocaA_R_bound(l,i,2) = cprocaA_R_bound(l,i,1)
        cprocaA_R_bound(l,i,1) = cprocaA_R(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaA_I_p(l,:) = cprocaA_I(l,:)
     do i=0,ghost-1
        cprocaA_I_bound(l,i,3) = cprocaA_I_bound(l,i,2)
        cprocaA_I_bound(l,i,2) = cprocaA_I_bound(l,i,1)
        cprocaA_I_bound(l,i,1) = cprocaA_I(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaB_R_p(l,:) = cprocaB_R(l,:)
     do i=0,ghost-1
        cprocaB_R_bound(l,i,3) = cprocaB_R_bound(l,i,2)
        cprocaB_R_bound(l,i,2) = cprocaB_R_bound(l,i,1)
        cprocaB_R_bound(l,i,1) = cprocaB_R(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaB_I_p(l,:) = cprocaB_I(l,:)
     do i=0,ghost-1
        cprocaB_I_bound(l,i,3) = cprocaB_I_bound(l,i,2)
        cprocaB_I_bound(l,i,2) = cprocaB_I_bound(l,i,1)
        cprocaB_I_bound(l,i,1) = cprocaB_I(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaG_R_p(l,:) = cprocaG_R(l,:)
     do i=0,ghost-1
        cprocaG_R_bound(l,i,3) = cprocaG_R_bound(l,i,2)
        cprocaG_R_bound(l,i,2) = cprocaG_R_bound(l,i,1)
        cprocaG_R_bound(l,i,1) = cprocaG_R(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaG_I_p(l,:) = cprocaG_I(l,:)
     do i=0,ghost-1
        cprocaG_I_bound(l,i,3) = cprocaG_I_bound(l,i,2)
        cprocaG_I_bound(l,i,2) = cprocaG_I_bound(l,i,1)
        cprocaG_I_bound(l,i,1) = cprocaG_I(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaL_R_p(l,:) = cprocaL_R(l,:)
     do i=0,ghost-1
        cprocaL_R_bound(l,i,3) = cprocaL_R_bound(l,i,2)
        cprocaL_R_bound(l,i,2) = cprocaL_R_bound(l,i,1)
        cprocaL_R_bound(l,i,1) = cprocaL_R(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     cprocaL_I_p(l,:) = cprocaL_I(l,:)
     do i=0,ghost-1
        cprocaL_I_bound(l,i,3) = cprocaL_I_bound(l,i,2)
        cprocaL_I_bound(l,i,2) = cprocaL_I_bound(l,i,1)
        cprocaL_I_bound(l,i,1) = cprocaL_I(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     dirac_FR_p(l,:) = dirac_FR(l,:)
     do i=0,ghost-1
        dirac_FR_bound(l,i,3) = dirac_FR_bound(l,i,2)
        dirac_FR_bound(l,i,2) = dirac_FR_bound(l,i,1)
        dirac_FR_bound(l,i,1) = dirac_FR(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     dirac_FI_p(l,:) = dirac_FI(l,:)
     do i=0,ghost-1
        dirac_FI_bound(l,i,3) = dirac_FI_bound(l,i,2)
        dirac_FI_bound(l,i,2) = dirac_FI_bound(l,i,1)
        dirac_FI_bound(l,i,1) = dirac_FI(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     dirac_GR_p(l,:) = dirac_GR(l,:)
     do i=0,ghost-1
        dirac_GR_bound(l,i,3) = dirac_GR_bound(l,i,2)
        dirac_GR_bound(l,i,2) = dirac_GR_bound(l,i,1)
        dirac_GR_bound(l,i,1) = dirac_GR(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     dirac_GI_p(l,:) = dirac_GI(l,:)
     do i=0,ghost-1
        dirac_GI_bound(l,i,3) = dirac_GI_bound(l,i,2)
        dirac_GI_bound(l,i,2) = dirac_GI_bound(l,i,1)
        dirac_GI_bound(l,i,1) = dirac_GI(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"dust")) then
     dust_cD_p(l,:) = dust_cD(l,:)
     do i=0,ghost-1
        dust_cD_bound(l,i,3) = dust_cD_bound(l,i,2)
        dust_cD_bound(l,i,2) = dust_cD_bound(l,i,1)
        dust_cD_bound(l,i,1) = dust_cD(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"dust")) then
     dust_cE_p(l,:) = dust_cE(l,:)
     do i=0,ghost-1
        dust_cE_bound(l,i,3) = dust_cE_bound(l,i,2)
        dust_cE_bound(l,i,2) = dust_cE_bound(l,i,1)
        dust_cE_bound(l,i,1) = dust_cE(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"dust")) then
     dust_cS_p(l,:) = dust_cS(l,:)
     do i=0,ghost-1
        dust_cS_bound(l,i,3) = dust_cS_bound(l,i,2)
        dust_cS_bound(l,i,2) = dust_cS_bound(l,i,1)
        dust_cS_bound(l,i,1) = dust_cS(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"fluid")) then
     fluid_cD_p(l,:) = fluid_cD(l,:)
     do i=0,ghost-1
        fluid_cD_bound(l,i,3) = fluid_cD_bound(l,i,2)
        fluid_cD_bound(l,i,2) = fluid_cD_bound(l,i,1)
        fluid_cD_bound(l,i,1) = fluid_cD(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"fluid")) then
     fluid_cE_p(l,:) = fluid_cE(l,:)
     do i=0,ghost-1
        fluid_cE_bound(l,i,3) = fluid_cE_bound(l,i,2)
        fluid_cE_bound(l,i,2) = fluid_cE_bound(l,i,1)
        fluid_cE_bound(l,i,1) = fluid_cE(l,Nr-i)
     end do
  end if

  if (contains(mattertype,"fluid")) then
     fluid_cS_p(l,:) = fluid_cS(l,:)
     do i=0,ghost-1
        fluid_cS_bound(l,i,3) = fluid_cS_bound(l,i,2)
        fluid_cS_bound(l,i,2) = fluid_cS_bound(l,i,1)
        fluid_cS_bound(l,i,1) = fluid_cS(l,Nr-i)
     end do
  end if

  if (cosmic_run) then
     cosmobg_tau_p(l) = cosmobg_tau(l)
  end if

  if (cosmic_run) then
     cosmobg_H_p(l) = cosmobg_H(l)
  end if

  if (cosmic_run) then
     cosmobg_a_p(l) = cosmobg_a(l)
  end if

  if (cosmic_run) then
     cosmobg_phi_p(l) = cosmobg_phi(l)
  end if

  if (cosmic_run) then
     cosmobg_trK_p(l) = cosmobg_trK(l)
  end if

  if (cosmic_run) then
     cosmobg_alpha_p(l) = cosmobg_alpha(l)
  end if

  if (cosmic_run) then
     cosmobg_scalar_phi_p(l) = cosmobg_scalar_phi(l)
  end if

  if (cosmic_run) then
     cosmobg_scalar_pi_p(l) = cosmobg_scalar_pi(l)
  end if

  if (cosmic_run) then
     cosmobg_complex_phiR_p(l) = cosmobg_complex_phiR(l)
  end if

  if (cosmic_run) then
     cosmobg_complex_phiI_p(l) = cosmobg_complex_phiI(l)
  end if

  if (cosmic_run) then
     cosmobg_complex_piR_p(l) = cosmobg_complex_piR(l)
  end if

  if (cosmic_run) then
     cosmobg_complex_piI_p(l) = cosmobg_complex_piI(l)
  end if

  end subroutine saveold

