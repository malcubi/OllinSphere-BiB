! Automatically generated file.  Do not edit!

  subroutine grabarray(varname)

  use param
  use arrays
  use procinfo

  implicit none
  logical exists
  character(len=*) varname

  exists = .false.

  if (varname=='auxarray') then
     exists = .true.
     savevar => auxarray
  end if

  if (varname=='auxarray2') then
     exists = .true.
     savevar => auxarray2
  end if

  if (varname=='r') then
     exists = .true.
     savevar => r
  end if

  if (varname=='r_area') then
     exists = .true.
     savevar => r_area
  end if

  if (varname=='D1_r_area') then
     exists = .true.
     savevar => D1_r_area
  end if

  if (varname=='alpha') then
     exists = .true.
     savevar => alpha
     grabvar_bound => alpha_bound
  end if

  if (varname=='salpha') then
     exists = .true.
     savevar => salpha
  end if

  if (varname=='D1_alpha') then
     exists = .true.
     savevar => D1_alpha
  end if

  if (varname=='D2_alpha') then
     exists = .true.
     savevar => D2_alpha
  end if

  if (varname=='DA_alpha') then
     exists = .true.
     savevar => DA_alpha
  end if

  if (varname=='Dcov2_alpha') then
     exists = .true.
     savevar => Dcov2_alpha
  end if

  if (varname=='Lapla_alpha') then
     exists = .true.
     savevar => Lapla_alpha
  end if

  if (varname=='DD_alphar') then
     exists = .true.
     savevar => DD_alphar
  end if

  if (varname=='falpha') then
     exists = .true.
     savevar => falpha
  end if

  if (varname=='tau_origin') then
     exists = .true.
     savevar0D => tau_origin
     nullify(savevar)
  end if

  if (varname=='beta') then
     exists = .true.
     savevar => beta
     grabvar_bound => beta_bound
  end if

  if (varname=='sbeta') then
     exists = .true.
     savevar => sbeta
  end if

  if (varname=='D1_beta') then
     exists = .true.
     savevar => D1_beta
  end if

  if (varname=='D2_beta') then
     exists = .true.
     savevar => D2_beta
  end if

  if (varname=='DA_beta') then
     exists = .true.
     savevar => DA_beta
  end if

  if (varname=='dtbeta') then
     exists = .true.
     savevar => dtbeta
     grabvar_bound => dtbeta_bound
  end if

  if (varname=='sdtbeta') then
     exists = .true.
     savevar => sdtbeta
  end if

  if (varname=='D1_dtbeta') then
     exists = .true.
     savevar => D1_dtbeta
  end if

  if (varname=='DA_dtbeta') then
     exists = .true.
     savevar => DA_dtbeta
  end if

  if (varname=='DIV_beta') then
     exists = .true.
     savevar => DIV_beta
  end if

  if (varname=='D1_DIV_beta') then
     exists = .true.
     savevar => D1_DIV_beta
  end if

  if (varname=='DD_beta') then
     exists = .true.
     savevar => DD_beta
  end if

  if (varname=='fdriver') then
     exists = .true.
     savevar => fdriver
     grabvar_bound => fdriver_bound
  end if

  if (varname=='sfdriver') then
     exists = .true.
     savevar => sfdriver
  end if

  if (varname=='D1_fdriver') then
     exists = .true.
     savevar => D1_fdriver
  end if

  if (varname=='DA_fdriver') then
     exists = .true.
     savevar => DA_fdriver
  end if

  if (varname=='phi') then
     exists = .true.
     savevar => phi
     grabvar_bound => phi_bound
  end if

  if (varname=='sphi') then
     exists = .true.
     savevar => sphi
  end if

  if (varname=='D1_phi') then
     exists = .true.
     savevar => D1_phi
  end if

  if (varname=='D2_phi') then
     exists = .true.
     savevar => D2_phi
  end if

  if (varname=='DA_phi') then
     exists = .true.
     savevar => DA_phi
  end if

  if (varname=='psi') then
     exists = .true.
     savevar => psi
  end if

  if (varname=='psi2') then
     exists = .true.
     savevar => psi2
  end if

  if (varname=='psi4') then
     exists = .true.
     savevar => psi4
  end if

  if (varname=='D1_psi') then
     exists = .true.
     savevar => D1_psi
  end if

  if (varname=='D2_psi') then
     exists = .true.
     savevar => D2_psi
  end if

  if (varname=='chi') then
     exists = .true.
     savevar => chi
     grabvar_bound => chi_bound
  end if

  if (varname=='schi') then
     exists = .true.
     savevar => schi
  end if

  if (varname=='D1_chi') then
     exists = .true.
     savevar => D1_chi
  end if

  if (varname=='D2_chi') then
     exists = .true.
     savevar => D2_chi
  end if

  if (varname=='DA_chi') then
     exists = .true.
     savevar => DA_chi
  end if

  if (varname=='DD_phir') then
     exists = .true.
     savevar => DD_phir
  end if

  if (varname=='A') then
     exists = .true.
     savevar => A
     grabvar_bound => A_bound
  end if

  if (varname=='sA') then
     exists = .true.
     savevar => sA
  end if

  if (varname=='D1_A') then
     exists = .true.
     savevar => D1_A
  end if

  if (varname=='D2_A') then
     exists = .true.
     savevar => D2_A
  end if

  if (varname=='B') then
     exists = .true.
     savevar => B
     grabvar_bound => B_bound
  end if

  if (varname=='sB') then
     exists = .true.
     savevar => sB
  end if

  if (varname=='D1_B') then
     exists = .true.
     savevar => D1_B
  end if

  if (varname=='D2_B') then
     exists = .true.
     savevar => D2_B
  end if

  if (varname=='DA_A') then
     exists = .true.
     savevar => DA_A
  end if

  if (varname=='DA_B') then
     exists = .true.
     savevar => DA_B
  end if

  if (varname=='DD_Ar') then
     exists = .true.
     savevar => DD_Ar
  end if

  if (varname=='DD_Br') then
     exists = .true.
     savevar => DD_Br
  end if

  if (varname=='APHYS') then
     exists = .true.
     savevar => APHYS
  end if

  if (varname=='BPHYS') then
     exists = .true.
     savevar => BPHYS
  end if

  if (varname=='AB2') then
     exists = .true.
     savevar => AB2
  end if

  if (varname=='D1_AB2') then
     exists = .true.
     savevar => D1_AB2
  end if

  if (varname=='D2_AB2') then
     exists = .true.
     savevar => D2_AB2
  end if

  if (varname=='GRR') then
     exists = .true.
     savevar => GRR
  end if

  if (varname=='hembed') then
     exists = .true.
     savevar => hembed
  end if

  if (varname=='trK') then
     exists = .true.
     savevar => trK
     grabvar_bound => trK_bound
  end if

  if (varname=='strK') then
     exists = .true.
     savevar => strK
  end if

  if (varname=='D1_trK') then
     exists = .true.
     savevar => D1_trK
  end if

  if (varname=='DA_trK') then
     exists = .true.
     savevar => DA_trK
  end if

  if (varname=='KTA') then
     exists = .true.
     savevar => KTA
     grabvar_bound => KTA_bound
  end if

  if (varname=='sKTA') then
     exists = .true.
     savevar => sKTA
  end if

  if (varname=='D1_KTA') then
     exists = .true.
     savevar => D1_KTA
  end if

  if (varname=='DA_KTA') then
     exists = .true.
     savevar => DA_KTA
  end if

  if (varname=='KTB') then
     exists = .true.
     savevar => KTB
  end if

  if (varname=='D1_KTB') then
     exists = .true.
     savevar => D1_KTB
  end if

  if (varname=='DA_KTB') then
     exists = .true.
     savevar => DA_KTB
  end if

  if (varname=='KAPHYS') then
     exists = .true.
     savevar => KAPHYS
  end if

  if (varname=='KBPHYS') then
     exists = .true.
     savevar => KBPHYS
  end if

  if (varname=='K2') then
     exists = .true.
     savevar => K2
  end if

  if (varname=='Deltar') then
     exists = .true.
     savevar => Deltar
     grabvar_bound => Deltar_bound
  end if

  if (varname=='sDeltar') then
     exists = .true.
     savevar => sDeltar
  end if

  if (varname=='D1_Deltar') then
     exists = .true.
     savevar => D1_Deltar
  end if

  if (varname=='DA_Deltar') then
     exists = .true.
     savevar => DA_Deltar
  end if

  if (varname=='Deltar0') then
     exists = .true.
     savevar => Deltar0
  end if

  if (varname=='DD_Deltar') then
     exists = .true.
     savevar => DD_Deltar
  end if

  if (varname=='DeltaAB') then
     exists = .true.
     savevar => DeltaAB
  end if

  if (varname=='D1_DeltaAB') then
     exists = .true.
     savevar => D1_DeltaAB
  end if

  if (varname=='lambda') then
     exists = .true.
     savevar => lambda
     grabvar_bound => lambda_bound
  end if

  if (varname=='slambda') then
     exists = .true.
     savevar => slambda
  end if

  if (varname=='D1_lambda') then
     exists = .true.
     savevar => D1_lambda
  end if

  if (varname=='D2_lambda') then
     exists = .true.
     savevar => D2_lambda
  end if

  if (varname=='DA_lambda') then
     exists = .true.
     savevar => DA_lambda
  end if

  if (varname=='Klambda') then
     exists = .true.
     savevar => Klambda
     grabvar_bound => Klambda_bound
  end if

  if (varname=='sKlambda') then
     exists = .true.
     savevar => sKlambda
  end if

  if (varname=='D1_Klambda') then
     exists = .true.
     savevar => D1_Klambda
  end if

  if (varname=='DA_Klambda') then
     exists = .true.
     savevar => DA_Klambda
  end if

  if (varname=='lambda2') then
     exists = .true.
     savevar => lambda2
     grabvar_bound => lambda2_bound
  end if

  if (varname=='slambda2') then
     exists = .true.
     savevar => slambda2
  end if

  if (varname=='D1_lambda2') then
     exists = .true.
     savevar => D1_lambda2
  end if

  if (varname=='D2_lambda2') then
     exists = .true.
     savevar => D2_lambda2
  end if

  if (varname=='DA_lambda2') then
     exists = .true.
     savevar => DA_lambda2
  end if

  if (varname=='Klambda2') then
     exists = .true.
     savevar => Klambda2
     grabvar_bound => Klambda2_bound
  end if

  if (varname=='sKlambda2') then
     exists = .true.
     savevar => sKlambda2
  end if

  if (varname=='D1_Klambda2') then
     exists = .true.
     savevar => D1_Klambda2
  end if

  if (varname=='DA_Klambda2') then
     exists = .true.
     savevar => DA_Klambda2
  end if

  if (varname=='RICA') then
     exists = .true.
     savevar => RICA
  end if

  if (varname=='RSCAL') then
     exists = .true.
     savevar => RSCAL
  end if

  if (varname=='z4theta') then
     exists = .true.
     savevar => z4theta
     grabvar_bound => z4theta_bound
  end if

  if (varname=='sz4theta') then
     exists = .true.
     savevar => sz4theta
  end if

  if (varname=='D1_z4theta') then
     exists = .true.
     savevar => D1_z4theta
  end if

  if (varname=='DA_z4theta') then
     exists = .true.
     savevar => DA_z4theta
  end if

  if (varname=='EWEYLA') then
     exists = .true.
     savevar => EWEYLA
  end if

  if (varname=='EWEYLB') then
     exists = .true.
     savevar => EWEYLB
  end if

  if (varname=='invariantI') then
     exists = .true.
     savevar => invariantI
  end if

  if (varname=='invariantJ') then
     exists = .true.
     savevar => invariantJ
  end if

  if (varname=='Ricci4D') then
     exists = .true.
     savevar => Ricci4D
  end if

  if (varname=='Ricci4D_E') then
     exists = .true.
     savevar => Ricci4D_E
  end if

  if (varname=='ham') then
     exists = .true.
     savevar => ham
  end if

  if (varname=='hamabs') then
     exists = .true.
     savevar => hamabs
  end if

  if (varname=='mom') then
     exists = .true.
     savevar => mom
  end if

  if (varname=='momabs') then
     exists = .true.
     savevar => momabs
  end if

  if (varname=='CDeltar') then
     exists = .true.
     savevar => CDeltar
  end if

  if (varname=='Clambda') then
     exists = .true.
     savevar => Clambda
  end if

  if (varname=='CKlambda') then
     exists = .true.
     savevar => CKlambda
  end if

  if (varname=='virial1') then
     exists = .true.
     savevar => virial1
  end if

  if (varname=='virial2') then
     exists = .true.
     savevar => virial2
  end if

  if (varname=='expansion') then
     exists = .true.
     savevar => expansion
  end if

  if (varname=='expansion_in') then
     exists = .true.
     savevar => expansion_in
  end if

  if (varname=='vlight') then
     exists = .true.
     savevar => vlight
  end if

  if (varname=='vgauge') then
     exists = .true.
     savevar => vgauge
  end if

  if (varname=='wp_gauge') then
     exists = .true.
     savevar => wp_gauge
  end if

  if (varname=='wm_gauge') then
     exists = .true.
     savevar => wm_gauge
  end if

  if (varname=='wp_eta') then
     exists = .true.
     savevar => wp_eta
  end if

  if (varname=='wm_eta') then
     exists = .true.
     savevar => wm_eta
  end if

  if (varname=='wDelta') then
     exists = .true.
     savevar => wDelta
  end if

  if (varname=='rho') then
     exists = .true.
     savevar => rho
  end if

  if (varname=='JA') then
     exists = .true.
     savevar => JA
  end if

  if (varname=='SAA') then
     exists = .true.
     savevar => SAA
  end if

  if (varname=='SBB') then
     exists = .true.
     savevar => SBB
  end if

  if (varname=='trS') then
     exists = .true.
     savevar => trS
  end if

  if (varname=='SLL') then
     exists = .true.
     savevar => SLL
  end if

  if (varname=='rho_r2') then
     exists = .true.
     savevar => rho_r2
  end if

  if (varname=='mass_sch') then
     exists = .true.
     savevar => mass_sch
  end if

  if (varname=='mass_int') then
     exists = .true.
     savevar => mass_int
  end if

  if (varname=='mass_rn') then
     exists = .true.
     savevar => mass_rn
  end if

  if (varname=='Kodama_mass') then
     exists = .true.
     savevar => Kodama_mass
  end if

  if (varname=='P_Kodama') then
     exists = .true.
     savevar => P_Kodama
     grabvar_bound => P_Kodama_bound
  end if

  if (varname=='sP_Kodama') then
     exists = .true.
     savevar => sP_Kodama
  end if

  if (varname=='rho_int') then
     exists = .true.
     savevar => rho_int
  end if

  if (varname=='compactness') then
     exists = .true.
     savevar => compactness
  end if

  if (varname=='trumpet_u') then
     exists = .true.
     savevar => trumpet_u
  end if

  if (varname=='trumpet_v') then
     exists = .true.
     savevar => trumpet_v
  end if

  if (varname=='scalar_phi') then
     exists = .true.
     savevar => scalar_phi
     grabvar_bound => scalar_phi_bound
  end if

  if (varname=='sscalar_phi') then
     exists = .true.
     savevar => sscalar_phi
  end if

  if (varname=='D1_scalar_phi') then
     exists = .true.
     savevar => D1_scalar_phi
  end if

  if (varname=='D2_scalar_phi') then
     exists = .true.
     savevar => D2_scalar_phi
  end if

  if (varname=='DA_scalar_phi') then
     exists = .true.
     savevar => DA_scalar_phi
  end if

  if (varname=='scalar_xi') then
     exists = .true.
     savevar => scalar_xi
     grabvar_bound => scalar_xi_bound
  end if

  if (varname=='sscalar_xi') then
     exists = .true.
     savevar => sscalar_xi
  end if

  if (varname=='scalar_pi') then
     exists = .true.
     savevar => scalar_pi
     grabvar_bound => scalar_pi_bound
  end if

  if (varname=='sscalar_pi') then
     exists = .true.
     savevar => sscalar_pi
  end if

  if (varname=='D1_scalar_xi') then
     exists = .true.
     savevar => D1_scalar_xi
  end if

  if (varname=='D1_scalar_pi') then
     exists = .true.
     savevar => D1_scalar_pi
  end if

  if (varname=='DA_scalar_xi') then
     exists = .true.
     savevar => DA_scalar_xi
  end if

  if (varname=='DA_scalar_pi') then
     exists = .true.
     savevar => DA_scalar_pi
  end if

  if (varname=='scalar_V') then
     exists = .true.
     savevar => scalar_V
  end if

  if (varname=='scalar_VP') then
     exists = .true.
     savevar => scalar_VP
  end if

  if (varname=='wp_scalar') then
     exists = .true.
     savevar => wp_scalar
  end if

  if (varname=='wm_scalar') then
     exists = .true.
     savevar => wm_scalar
  end if

  if (varname=='ghost_phi') then
     exists = .true.
     savevar => ghost_phi
     grabvar_bound => ghost_phi_bound
  end if

  if (varname=='sghost_phi') then
     exists = .true.
     savevar => sghost_phi
  end if

  if (varname=='D1_ghost_phi') then
     exists = .true.
     savevar => D1_ghost_phi
  end if

  if (varname=='D2_ghost_phi') then
     exists = .true.
     savevar => D2_ghost_phi
  end if

  if (varname=='DA_ghost_phi') then
     exists = .true.
     savevar => DA_ghost_phi
  end if

  if (varname=='ghost_xi') then
     exists = .true.
     savevar => ghost_xi
     grabvar_bound => ghost_xi_bound
  end if

  if (varname=='sghost_xi') then
     exists = .true.
     savevar => sghost_xi
  end if

  if (varname=='ghost_pi') then
     exists = .true.
     savevar => ghost_pi
     grabvar_bound => ghost_pi_bound
  end if

  if (varname=='sghost_pi') then
     exists = .true.
     savevar => sghost_pi
  end if

  if (varname=='D1_ghost_xi') then
     exists = .true.
     savevar => D1_ghost_xi
  end if

  if (varname=='D1_ghost_pi') then
     exists = .true.
     savevar => D1_ghost_pi
  end if

  if (varname=='DA_ghost_xi') then
     exists = .true.
     savevar => DA_ghost_xi
  end if

  if (varname=='DA_ghost_pi') then
     exists = .true.
     savevar => DA_ghost_pi
  end if

  if (varname=='ghost_V') then
     exists = .true.
     savevar => ghost_V
  end if

  if (varname=='ghost_VP') then
     exists = .true.
     savevar => ghost_VP
  end if

  if (varname=='complexghost_phiR') then
     exists = .true.
     savevar => complexghost_phiR
     grabvar_bound => complexghost_phiR_bound
  end if

  if (varname=='scomplexghost_phiR') then
     exists = .true.
     savevar => scomplexghost_phiR
  end if

  if (varname=='complexghost_phiI') then
     exists = .true.
     savevar => complexghost_phiI
     grabvar_bound => complexghost_phiI_bound
  end if

  if (varname=='scomplexghost_phiI') then
     exists = .true.
     savevar => scomplexghost_phiI
  end if

  if (varname=='D1_complexghost_phiR') then
     exists = .true.
     savevar => D1_complexghost_phiR
  end if

  if (varname=='D2_complexghost_phiR') then
     exists = .true.
     savevar => D2_complexghost_phiR
  end if

  if (varname=='DA_complexghost_phiR') then
     exists = .true.
     savevar => DA_complexghost_phiR
  end if

  if (varname=='D1_complexghost_phiI') then
     exists = .true.
     savevar => D1_complexghost_phiI
  end if

  if (varname=='D2_complexghost_phiI') then
     exists = .true.
     savevar => D2_complexghost_phiI
  end if

  if (varname=='DA_complexghost_phiI') then
     exists = .true.
     savevar => DA_complexghost_phiI
  end if

  if (varname=='complexghost_xiR') then
     exists = .true.
     savevar => complexghost_xiR
     grabvar_bound => complexghost_xiR_bound
  end if

  if (varname=='scomplexghost_xiR') then
     exists = .true.
     savevar => scomplexghost_xiR
  end if

  if (varname=='complexghost_xiI') then
     exists = .true.
     savevar => complexghost_xiI
     grabvar_bound => complexghost_xiI_bound
  end if

  if (varname=='scomplexghost_xiI') then
     exists = .true.
     savevar => scomplexghost_xiI
  end if

  if (varname=='complexghost_piR') then
     exists = .true.
     savevar => complexghost_piR
     grabvar_bound => complexghost_piR_bound
  end if

  if (varname=='scomplexghost_piR') then
     exists = .true.
     savevar => scomplexghost_piR
  end if

  if (varname=='complexghost_piI') then
     exists = .true.
     savevar => complexghost_piI
     grabvar_bound => complexghost_piI_bound
  end if

  if (varname=='scomplexghost_piI') then
     exists = .true.
     savevar => scomplexghost_piI
  end if

  if (varname=='D1_complexghost_xiR') then
     exists = .true.
     savevar => D1_complexghost_xiR
  end if

  if (varname=='D1_complexghost_xiI') then
     exists = .true.
     savevar => D1_complexghost_xiI
  end if

  if (varname=='DA_complexghost_xiR') then
     exists = .true.
     savevar => DA_complexghost_xiR
  end if

  if (varname=='DA_complexghost_xiI') then
     exists = .true.
     savevar => DA_complexghost_xiI
  end if

  if (varname=='D1_complexghost_piR') then
     exists = .true.
     savevar => D1_complexghost_piR
  end if

  if (varname=='D1_complexghost_piI') then
     exists = .true.
     savevar => D1_complexghost_piI
  end if

  if (varname=='DA_complexghost_piR') then
     exists = .true.
     savevar => DA_complexghost_piR
  end if

  if (varname=='DA_complexghost_piI') then
     exists = .true.
     savevar => DA_complexghost_piI
  end if

  if (varname=='complexghost_V') then
     exists = .true.
     savevar => complexghost_V
  end if

  if (varname=='complexghost_VPR') then
     exists = .true.
     savevar => complexghost_VPR
  end if

  if (varname=='complexghost_VPI') then
     exists = .true.
     savevar => complexghost_VPI
  end if

  if (varname=='complexghost_phi_norm') then
     exists = .true.
     savevar => complexghost_phi_norm
  end if

  if (varname=='complex_rhotot') then
     exists = .true.
     savevar => complex_rhotot
  end if

  if (varname=='complex_phiR') then
     exists = .true.
     savevar => complex_phiR
     grabvar_bound => complex_phiR_bound
  end if

  if (varname=='scomplex_phiR') then
     exists = .true.
     savevar => scomplex_phiR
  end if

  if (varname=='complex_phiI') then
     exists = .true.
     savevar => complex_phiI
     grabvar_bound => complex_phiI_bound
  end if

  if (varname=='scomplex_phiI') then
     exists = .true.
     savevar => scomplex_phiI
  end if

  if (varname=='D1_complex_phiR') then
     exists = .true.
     savevar => D1_complex_phiR
  end if

  if (varname=='D2_complex_phiR') then
     exists = .true.
     savevar => D2_complex_phiR
  end if

  if (varname=='DA_complex_phiR') then
     exists = .true.
     savevar => DA_complex_phiR
  end if

  if (varname=='D1_complex_phiI') then
     exists = .true.
     savevar => D1_complex_phiI
  end if

  if (varname=='D2_complex_phiI') then
     exists = .true.
     savevar => D2_complex_phiI
  end if

  if (varname=='DA_complex_phiI') then
     exists = .true.
     savevar => DA_complex_phiI
  end if

  if (varname=='complex_xiR') then
     exists = .true.
     savevar => complex_xiR
     grabvar_bound => complex_xiR_bound
  end if

  if (varname=='scomplex_xiR') then
     exists = .true.
     savevar => scomplex_xiR
  end if

  if (varname=='complex_xiI') then
     exists = .true.
     savevar => complex_xiI
     grabvar_bound => complex_xiI_bound
  end if

  if (varname=='scomplex_xiI') then
     exists = .true.
     savevar => scomplex_xiI
  end if

  if (varname=='complex_piR') then
     exists = .true.
     savevar => complex_piR
     grabvar_bound => complex_piR_bound
  end if

  if (varname=='scomplex_piR') then
     exists = .true.
     savevar => scomplex_piR
  end if

  if (varname=='complex_piI') then
     exists = .true.
     savevar => complex_piI
     grabvar_bound => complex_piI_bound
  end if

  if (varname=='scomplex_piI') then
     exists = .true.
     savevar => scomplex_piI
  end if

  if (varname=='complex_gxiR') then
     exists = .true.
     savevar => complex_gxiR
  end if

  if (varname=='complex_gxiI') then
     exists = .true.
     savevar => complex_gxiI
  end if

  if (varname=='complex_gpiR') then
     exists = .true.
     savevar => complex_gpiR
  end if

  if (varname=='complex_gpiI') then
     exists = .true.
     savevar => complex_gpiI
  end if

  if (varname=='D1_complex_xiR') then
     exists = .true.
     savevar => D1_complex_xiR
  end if

  if (varname=='D1_complex_xiI') then
     exists = .true.
     savevar => D1_complex_xiI
  end if

  if (varname=='DA_complex_xiR') then
     exists = .true.
     savevar => DA_complex_xiR
  end if

  if (varname=='DA_complex_xiI') then
     exists = .true.
     savevar => DA_complex_xiI
  end if

  if (varname=='D1_complex_piR') then
     exists = .true.
     savevar => D1_complex_piR
  end if

  if (varname=='D1_complex_piI') then
     exists = .true.
     savevar => D1_complex_piI
  end if

  if (varname=='DA_complex_piR') then
     exists = .true.
     savevar => DA_complex_piR
  end if

  if (varname=='DA_complex_piI') then
     exists = .true.
     savevar => DA_complex_piI
  end if

  if (varname=='complex_V') then
     exists = .true.
     savevar => complex_V
  end if

  if (varname=='complex_VPR') then
     exists = .true.
     savevar => complex_VPR
  end if

  if (varname=='complex_VPI') then
     exists = .true.
     savevar => complex_VPI
  end if

  if (varname=='complex_phi_norm') then
     exists = .true.
     savevar => complex_phi_norm
  end if

  if (varname=='complex_Bdens') then
     exists = .true.
     savevar => complex_Bdens
  end if

  if (varname=='complex_Bdens_r2') then
     exists = .true.
     savevar => complex_Bdens_r2
  end if

  if (varname=='complex_Bflux') then
     exists = .true.
     savevar => complex_Bflux
  end if

  if (varname=='complex_NB') then
     exists = .true.
     savevar => complex_NB
  end if

  if (varname=='complex_phiaux') then
     exists = .true.
     savevar => complex_phiaux
  end if

  if (varname=='complex_phiR_l0') then
     exists = .true.
     savevar => complex_phiR_l0
  end if

  if (varname=='complex_phiR_l1') then
     exists = .true.
     savevar => complex_phiR_l1
  end if

  if (varname=='complex_phiR_l2') then
     exists = .true.
     savevar => complex_phiR_l2
  end if

  if (varname=='complex_phiR_l3') then
     exists = .true.
     savevar => complex_phiR_l3
  end if

  if (varname=='nonmin_phi') then
     exists = .true.
     savevar => nonmin_phi
     grabvar_bound => nonmin_phi_bound
  end if

  if (varname=='snonmin_phi') then
     exists = .true.
     savevar => snonmin_phi
  end if

  if (varname=='D1_nonmin_phi') then
     exists = .true.
     savevar => D1_nonmin_phi
  end if

  if (varname=='D2_nonmin_phi') then
     exists = .true.
     savevar => D2_nonmin_phi
  end if

  if (varname=='DA_nonmin_phi') then
     exists = .true.
     savevar => DA_nonmin_phi
  end if

  if (varname=='nonmin_xi') then
     exists = .true.
     savevar => nonmin_xi
     grabvar_bound => nonmin_xi_bound
  end if

  if (varname=='snonmin_xi') then
     exists = .true.
     savevar => snonmin_xi
  end if

  if (varname=='nonmin_pi') then
     exists = .true.
     savevar => nonmin_pi
     grabvar_bound => nonmin_pi_bound
  end if

  if (varname=='snonmin_pi') then
     exists = .true.
     savevar => snonmin_pi
  end if

  if (varname=='D1_nonmin_xi') then
     exists = .true.
     savevar => D1_nonmin_xi
  end if

  if (varname=='D1_nonmin_pi') then
     exists = .true.
     savevar => D1_nonmin_pi
  end if

  if (varname=='DA_nonmin_xi') then
     exists = .true.
     savevar => DA_nonmin_xi
  end if

  if (varname=='DA_nonmin_pi') then
     exists = .true.
     savevar => DA_nonmin_pi
  end if

  if (varname=='nonmin_V') then
     exists = .true.
     savevar => nonmin_V
  end if

  if (varname=='nonmin_VP') then
     exists = .true.
     savevar => nonmin_VP
  end if

  if (varname=='nonmin_Q') then
     exists = .true.
     savevar => nonmin_Q
  end if

  if (varname=='nonmin_f') then
     exists = .true.
     savevar => nonmin_f
  end if

  if (varname=='nonmin_fp') then
     exists = .true.
     savevar => nonmin_fp
  end if

  if (varname=='nonmin_fpp') then
     exists = .true.
     savevar => nonmin_fpp
  end if

  if (varname=='nonmin_rhs') then
     exists = .true.
     savevar => nonmin_rhs
  end if

  if (varname=='DD_nonmin_xir') then
     exists = .true.
     savevar => DD_nonmin_xir
  end if

  if (varname=='electric') then
     exists = .true.
     savevar => electric
     grabvar_bound => electric_bound
  end if

  if (varname=='selectric') then
     exists = .true.
     savevar => selectric
  end if

  if (varname=='D1_electric') then
     exists = .true.
     savevar => D1_electric
  end if

  if (varname=='DA_electric') then
     exists = .true.
     savevar => DA_electric
  end if

  if (varname=='ePhi') then
     exists = .true.
     savevar => ePhi
     grabvar_bound => ePhi_bound
  end if

  if (varname=='sePhi') then
     exists = .true.
     savevar => sePhi
  end if

  if (varname=='D1_ePhi') then
     exists = .true.
     savevar => D1_ePhi
  end if

  if (varname=='DA_ePhi') then
     exists = .true.
     savevar => DA_ePhi
  end if

  if (varname=='eAr') then
     exists = .true.
     savevar => eAr
     grabvar_bound => eAr_bound
  end if

  if (varname=='seAr') then
     exists = .true.
     savevar => seAr
  end if

  if (varname=='D1_eAr') then
     exists = .true.
     savevar => D1_eAr
  end if

  if (varname=='DA_eAr') then
     exists = .true.
     savevar => DA_eAr
  end if

  if (varname=='eF') then
     exists = .true.
     savevar => eF
  end if

  if (varname=='D1_eF') then
     exists = .true.
     savevar => D1_eF
  end if

  if (varname=='eH') then
     exists = .true.
     savevar => eH
  end if

  if (varname=='D1_eH') then
     exists = .true.
     savevar => D1_eH
  end if

  if (varname=='echarge') then
     exists = .true.
     savevar => echarge
  end if

  if (varname=='ecurrent') then
     exists = .true.
     savevar => ecurrent
  end if

  if (varname=='Celectric') then
     exists = .true.
     savevar => Celectric
  end if

  if (varname=='eQ_int') then
     exists = .true.
     savevar => eQ_int
  end if

  if (varname=='eQ_surf') then
     exists = .true.
     savevar => eQ_surf
  end if

  if (varname=='wp_electric') then
     exists = .true.
     savevar => wp_electric
  end if

  if (varname=='wm_electric') then
     exists = .true.
     savevar => wm_electric
  end if

  if (varname=='procaE') then
     exists = .true.
     savevar => procaE
     grabvar_bound => procaE_bound
  end if

  if (varname=='sprocaE') then
     exists = .true.
     savevar => sprocaE
  end if

  if (varname=='D1_procaE') then
     exists = .true.
     savevar => D1_procaE
  end if

  if (varname=='DA_procaE') then
     exists = .true.
     savevar => DA_procaE
  end if

  if (varname=='procaPhi') then
     exists = .true.
     savevar => procaPhi
     grabvar_bound => procaPhi_bound
  end if

  if (varname=='sprocaPhi') then
     exists = .true.
     savevar => sprocaPhi
  end if

  if (varname=='D1_procaPhi') then
     exists = .true.
     savevar => D1_procaPhi
  end if

  if (varname=='DA_procaPhi') then
     exists = .true.
     savevar => DA_procaPhi
  end if

  if (varname=='procaA') then
     exists = .true.
     savevar => procaA
     grabvar_bound => procaA_bound
  end if

  if (varname=='sprocaA') then
     exists = .true.
     savevar => sprocaA
  end if

  if (varname=='D1_procaA') then
     exists = .true.
     savevar => D1_procaA
  end if

  if (varname=='DA_procaA') then
     exists = .true.
     savevar => DA_procaA
  end if

  if (varname=='procaF') then
     exists = .true.
     savevar => procaF
  end if

  if (varname=='D1_procaF') then
     exists = .true.
     savevar => D1_procaF
  end if

  if (varname=='procaH') then
     exists = .true.
     savevar => procaH
  end if

  if (varname=='D1_procaH') then
     exists = .true.
     savevar => D1_procaH
  end if

  if (varname=='procaV') then
     exists = .true.
     savevar => procaV
  end if

  if (varname=='Cproca') then
     exists = .true.
     savevar => Cproca
  end if

  if (varname=='wp_proca') then
     exists = .true.
     savevar => wp_proca
  end if

  if (varname=='wm_proca') then
     exists = .true.
     savevar => wm_proca
  end if

  if (varname=='cprocaE_R') then
     exists = .true.
     savevar => cprocaE_R
     grabvar_bound => cprocaE_R_bound
  end if

  if (varname=='scprocaE_R') then
     exists = .true.
     savevar => scprocaE_R
  end if

  if (varname=='D1_cprocaE_R') then
     exists = .true.
     savevar => D1_cprocaE_R
  end if

  if (varname=='DA_cprocaE_R') then
     exists = .true.
     savevar => DA_cprocaE_R
  end if

  if (varname=='cprocaE_I') then
     exists = .true.
     savevar => cprocaE_I
     grabvar_bound => cprocaE_I_bound
  end if

  if (varname=='scprocaE_I') then
     exists = .true.
     savevar => scprocaE_I
  end if

  if (varname=='D1_cprocaE_I') then
     exists = .true.
     savevar => D1_cprocaE_I
  end if

  if (varname=='DA_cprocaE_I') then
     exists = .true.
     savevar => DA_cprocaE_I
  end if

  if (varname=='cprocaXi_R') then
     exists = .true.
     savevar => cprocaXi_R
     grabvar_bound => cprocaXi_R_bound
  end if

  if (varname=='scprocaXi_R') then
     exists = .true.
     savevar => scprocaXi_R
  end if

  if (varname=='D1_cprocaXi_R') then
     exists = .true.
     savevar => D1_cprocaXi_R
  end if

  if (varname=='DA_cprocaXi_R') then
     exists = .true.
     savevar => DA_cprocaXi_R
  end if

  if (varname=='cprocaXi_I') then
     exists = .true.
     savevar => cprocaXi_I
     grabvar_bound => cprocaXi_I_bound
  end if

  if (varname=='scprocaXi_I') then
     exists = .true.
     savevar => scprocaXi_I
  end if

  if (varname=='D1_cprocaXi_I') then
     exists = .true.
     savevar => D1_cprocaXi_I
  end if

  if (varname=='DA_cprocaXi_I') then
     exists = .true.
     savevar => DA_cprocaXi_I
  end if

  if (varname=='cprocaPhi_R') then
     exists = .true.
     savevar => cprocaPhi_R
     grabvar_bound => cprocaPhi_R_bound
  end if

  if (varname=='scprocaPhi_R') then
     exists = .true.
     savevar => scprocaPhi_R
  end if

  if (varname=='D1_cprocaPhi_R') then
     exists = .true.
     savevar => D1_cprocaPhi_R
  end if

  if (varname=='DA_cprocaPhi_R') then
     exists = .true.
     savevar => DA_cprocaPhi_R
  end if

  if (varname=='cprocaPhi_I') then
     exists = .true.
     savevar => cprocaPhi_I
     grabvar_bound => cprocaPhi_I_bound
  end if

  if (varname=='scprocaPhi_I') then
     exists = .true.
     savevar => scprocaPhi_I
  end if

  if (varname=='D1_cprocaPhi_I') then
     exists = .true.
     savevar => D1_cprocaPhi_I
  end if

  if (varname=='DA_cprocaPhi_I') then
     exists = .true.
     savevar => DA_cprocaPhi_I
  end if

  if (varname=='cprocaA_R') then
     exists = .true.
     savevar => cprocaA_R
     grabvar_bound => cprocaA_R_bound
  end if

  if (varname=='scprocaA_R') then
     exists = .true.
     savevar => scprocaA_R
  end if

  if (varname=='D1_cprocaA_R') then
     exists = .true.
     savevar => D1_cprocaA_R
  end if

  if (varname=='DA_cprocaA_R') then
     exists = .true.
     savevar => DA_cprocaA_R
  end if

  if (varname=='cprocaA_I') then
     exists = .true.
     savevar => cprocaA_I
     grabvar_bound => cprocaA_I_bound
  end if

  if (varname=='scprocaA_I') then
     exists = .true.
     savevar => scprocaA_I
  end if

  if (varname=='D1_cprocaA_I') then
     exists = .true.
     savevar => D1_cprocaA_I
  end if

  if (varname=='DA_cprocaA_I') then
     exists = .true.
     savevar => DA_cprocaA_I
  end if

  if (varname=='cprocaF_R') then
     exists = .true.
     savevar => cprocaF_R
  end if

  if (varname=='D1_cprocaF_R') then
     exists = .true.
     savevar => D1_cprocaF_R
  end if

  if (varname=='cprocaF_I') then
     exists = .true.
     savevar => cprocaF_I
  end if

  if (varname=='D1_cprocaF_I') then
     exists = .true.
     savevar => D1_cprocaF_I
  end if

  if (varname=='cprocaH_R') then
     exists = .true.
     savevar => cprocaH_R
  end if

  if (varname=='D1_cprocaH_R') then
     exists = .true.
     savevar => D1_cprocaH_R
  end if

  if (varname=='cprocaH_I') then
     exists = .true.
     savevar => cprocaH_I
  end if

  if (varname=='D1_cprocaH_I') then
     exists = .true.
     savevar => D1_cprocaH_I
  end if

  if (varname=='cprocaB_R') then
     exists = .true.
     savevar => cprocaB_R
     grabvar_bound => cprocaB_R_bound
  end if

  if (varname=='scprocaB_R') then
     exists = .true.
     savevar => scprocaB_R
  end if

  if (varname=='D1_cprocaB_R') then
     exists = .true.
     savevar => D1_cprocaB_R
  end if

  if (varname=='D2_cprocaB_R') then
     exists = .true.
     savevar => D2_cprocaB_R
  end if

  if (varname=='DA_cprocaB_R') then
     exists = .true.
     savevar => DA_cprocaB_R
  end if

  if (varname=='cprocaB_I') then
     exists = .true.
     savevar => cprocaB_I
     grabvar_bound => cprocaB_I_bound
  end if

  if (varname=='scprocaB_I') then
     exists = .true.
     savevar => scprocaB_I
  end if

  if (varname=='D1_cprocaB_I') then
     exists = .true.
     savevar => D1_cprocaB_I
  end if

  if (varname=='D2_cprocaB_I') then
     exists = .true.
     savevar => D2_cprocaB_I
  end if

  if (varname=='DA_cprocaB_I') then
     exists = .true.
     savevar => DA_cprocaB_I
  end if

  if (varname=='cprocaG_R') then
     exists = .true.
     savevar => cprocaG_R
     grabvar_bound => cprocaG_R_bound
  end if

  if (varname=='scprocaG_R') then
     exists = .true.
     savevar => scprocaG_R
  end if

  if (varname=='D1_cprocaG_R') then
     exists = .true.
     savevar => D1_cprocaG_R
  end if

  if (varname=='DA_cprocaG_R') then
     exists = .true.
     savevar => DA_cprocaG_R
  end if

  if (varname=='cprocaG_I') then
     exists = .true.
     savevar => cprocaG_I
     grabvar_bound => cprocaG_I_bound
  end if

  if (varname=='scprocaG_I') then
     exists = .true.
     savevar => scprocaG_I
  end if

  if (varname=='D1_cprocaG_I') then
     exists = .true.
     savevar => D1_cprocaG_I
  end if

  if (varname=='DA_cprocaG_I') then
     exists = .true.
     savevar => DA_cprocaG_I
  end if

  if (varname=='cprocaL_R') then
     exists = .true.
     savevar => cprocaL_R
     grabvar_bound => cprocaL_R_bound
  end if

  if (varname=='scprocaL_R') then
     exists = .true.
     savevar => scprocaL_R
  end if

  if (varname=='cprocaL_I') then
     exists = .true.
     savevar => cprocaL_I
     grabvar_bound => cprocaL_I_bound
  end if

  if (varname=='scprocaL_I') then
     exists = .true.
     savevar => scprocaL_I
  end if

  if (varname=='Ccomplexproca_R') then
     exists = .true.
     savevar => Ccomplexproca_R
  end if

  if (varname=='Ccomplexproca_I') then
     exists = .true.
     savevar => Ccomplexproca_I
  end if

  if (varname=='cprocaPhi_norm') then
     exists = .true.
     savevar => cprocaPhi_norm
  end if

  if (varname=='cprocaA_norm') then
     exists = .true.
     savevar => cprocaA_norm
  end if

  if (varname=='cprocaB_norm') then
     exists = .true.
     savevar => cprocaB_norm
  end if

  if (varname=='cprocaE_norm') then
     exists = .true.
     savevar => cprocaE_norm
  end if

  if (varname=='cprocaXi_norm') then
     exists = .true.
     savevar => cprocaXi_norm
  end if

  if (varname=='Ccomplexproca_norm') then
     exists = .true.
     savevar => Ccomplexproca_norm
  end if

  if (varname=='cproca_Qdens') then
     exists = .true.
     savevar => cproca_Qdens
  end if

  if (varname=='cproca_Qdens_r2') then
     exists = .true.
     savevar => cproca_Qdens_r2
  end if

  if (varname=='cproca_Qflux') then
     exists = .true.
     savevar => cproca_Qflux
  end if

  if (varname=='cproca_Qint') then
     exists = .true.
     savevar => cproca_Qint
  end if

  if (varname=='dirac_FR') then
     exists = .true.
     savevar => dirac_FR
     grabvar_bound => dirac_FR_bound
  end if

  if (varname=='sdirac_FR') then
     exists = .true.
     savevar => sdirac_FR
  end if

  if (varname=='dirac_FI') then
     exists = .true.
     savevar => dirac_FI
     grabvar_bound => dirac_FI_bound
  end if

  if (varname=='sdirac_FI') then
     exists = .true.
     savevar => sdirac_FI
  end if

  if (varname=='dirac_GR') then
     exists = .true.
     savevar => dirac_GR
     grabvar_bound => dirac_GR_bound
  end if

  if (varname=='sdirac_GR') then
     exists = .true.
     savevar => sdirac_GR
  end if

  if (varname=='dirac_GI') then
     exists = .true.
     savevar => dirac_GI
     grabvar_bound => dirac_GI_bound
  end if

  if (varname=='sdirac_GI') then
     exists = .true.
     savevar => sdirac_GI
  end if

  if (varname=='D1_dirac_FR') then
     exists = .true.
     savevar => D1_dirac_FR
  end if

  if (varname=='D1_dirac_FI') then
     exists = .true.
     savevar => D1_dirac_FI
  end if

  if (varname=='DA_dirac_FR') then
     exists = .true.
     savevar => DA_dirac_FR
  end if

  if (varname=='DA_dirac_FI') then
     exists = .true.
     savevar => DA_dirac_FI
  end if

  if (varname=='D1_dirac_GR') then
     exists = .true.
     savevar => D1_dirac_GR
  end if

  if (varname=='D1_dirac_GI') then
     exists = .true.
     savevar => D1_dirac_GI
  end if

  if (varname=='DA_dirac_GR') then
     exists = .true.
     savevar => DA_dirac_GR
  end if

  if (varname=='DA_dirac_GI') then
     exists = .true.
     savevar => DA_dirac_GI
  end if

  if (varname=='dirac_HR') then
     exists = .true.
     savevar => dirac_HR
  end if

  if (varname=='dirac_HI') then
     exists = .true.
     savevar => dirac_HI
  end if

  if (varname=='D1_dirac_HR') then
     exists = .true.
     savevar => D1_dirac_HR
  end if

  if (varname=='D1_dirac_HI') then
     exists = .true.
     savevar => D1_dirac_HI
  end if

  if (varname=='dirac_PiFR') then
     exists = .true.
     savevar => dirac_PiFR
  end if

  if (varname=='dirac_PiFI') then
     exists = .true.
     savevar => dirac_PiFI
  end if

  if (varname=='dirac_PiGR') then
     exists = .true.
     savevar => dirac_PiGR
  end if

  if (varname=='dirac_PiGI') then
     exists = .true.
     savevar => dirac_PiGI
  end if

  if (varname=='dirac_F_norm') then
     exists = .true.
     savevar => dirac_F_norm
  end if

  if (varname=='dirac_G_norm') then
     exists = .true.
     savevar => dirac_G_norm
  end if

  if (varname=='dirac_dens') then
     exists = .true.
     savevar => dirac_dens
  end if

  if (varname=='dirac_flux') then
     exists = .true.
     savevar => dirac_flux
  end if

  if (varname=='dirac_Nint') then
     exists = .true.
     savevar => dirac_Nint
  end if

  if (varname=='wpR_dirac') then
     exists = .true.
     savevar => wpR_dirac
  end if

  if (varname=='wmR_dirac') then
     exists = .true.
     savevar => wmR_dirac
  end if

  if (varname=='wpI_dirac') then
     exists = .true.
     savevar => wpI_dirac
  end if

  if (varname=='wmI_dirac') then
     exists = .true.
     savevar => wmI_dirac
  end if

  if (varname=='dust_rho') then
     exists = .true.
     savevar => dust_rho
  end if

  if (varname=='dust_v') then
     exists = .true.
     savevar => dust_v
  end if

  if (varname=='dust_W') then
     exists = .true.
     savevar => dust_W
  end if

  if (varname=='dust_cD') then
     exists = .true.
     savevar => dust_cD
     grabvar_bound => dust_cD_bound
  end if

  if (varname=='sdust_cD') then
     exists = .true.
     savevar => sdust_cD
  end if

  if (varname=='D1_dust_cD') then
     exists = .true.
     savevar => D1_dust_cD
  end if

  if (varname=='DA_dust_cD') then
     exists = .true.
     savevar => DA_dust_cD
  end if

  if (varname=='dust_cE') then
     exists = .true.
     savevar => dust_cE
     grabvar_bound => dust_cE_bound
  end if

  if (varname=='sdust_cE') then
     exists = .true.
     savevar => sdust_cE
  end if

  if (varname=='D1_dust_cE') then
     exists = .true.
     savevar => D1_dust_cE
  end if

  if (varname=='DA_dust_cE') then
     exists = .true.
     savevar => DA_dust_cE
  end if

  if (varname=='dust_cS') then
     exists = .true.
     savevar => dust_cS
     grabvar_bound => dust_cS_bound
  end if

  if (varname=='sdust_cS') then
     exists = .true.
     savevar => sdust_cS
  end if

  if (varname=='D1_dust_cS') then
     exists = .true.
     savevar => D1_dust_cS
  end if

  if (varname=='DA_dust_cS') then
     exists = .true.
     savevar => DA_dust_cS
  end if

  if (varname=='dust_restmass') then
     exists = .true.
     savevar => dust_restmass
  end if

  if (varname=='fluid_rho') then
     exists = .true.
     savevar => fluid_rho
  end if

  if (varname=='fluid_rhotot') then
     exists = .true.
     savevar => fluid_rhotot
  end if

  if (varname=='fluid_e') then
     exists = .true.
     savevar => fluid_e
  end if

  if (varname=='fluid_p') then
     exists = .true.
     savevar => fluid_p
  end if

  if (varname=='D1_fluid_p') then
     exists = .true.
     savevar => D1_fluid_p
  end if

  if (varname=='fluid_h') then
     exists = .true.
     savevar => fluid_h
  end if

  if (varname=='fluid_v') then
     exists = .true.
     savevar => fluid_v
  end if

  if (varname=='fluid_W') then
     exists = .true.
     savevar => fluid_W
  end if

  if (varname=='fluid_vs') then
     exists = .true.
     savevar => fluid_vs
  end if

  if (varname=='fluid_Mach') then
     exists = .true.
     savevar => fluid_Mach
  end if

  if (varname=='fluid_vcp') then
     exists = .true.
     savevar => fluid_vcp
  end if

  if (varname=='fluid_vcm') then
     exists = .true.
     savevar => fluid_vcm
  end if

  if (varname=='fluid_wb_rest') then
     exists = .true.
     savevar => fluid_wb_rest
  end if

  if (varname=='fluid_wb_tot') then
     exists = .true.
     savevar => fluid_wb_tot
  end if

  if (varname=='fluid_cD') then
     exists = .true.
     savevar => fluid_cD
     grabvar_bound => fluid_cD_bound
  end if

  if (varname=='sfluid_cD') then
     exists = .true.
     savevar => sfluid_cD
  end if

  if (varname=='D1_fluid_cD') then
     exists = .true.
     savevar => D1_fluid_cD
  end if

  if (varname=='DA_fluid_cD') then
     exists = .true.
     savevar => DA_fluid_cD
  end if

  if (varname=='fluid_cE') then
     exists = .true.
     savevar => fluid_cE
     grabvar_bound => fluid_cE_bound
  end if

  if (varname=='sfluid_cE') then
     exists = .true.
     savevar => sfluid_cE
  end if

  if (varname=='D1_fluid_cE') then
     exists = .true.
     savevar => D1_fluid_cE
  end if

  if (varname=='DA_fluid_cE') then
     exists = .true.
     savevar => DA_fluid_cE
  end if

  if (varname=='fluid_cS') then
     exists = .true.
     savevar => fluid_cS
     grabvar_bound => fluid_cS_bound
  end if

  if (varname=='sfluid_cS') then
     exists = .true.
     savevar => sfluid_cS
  end if

  if (varname=='D1_fluid_cS') then
     exists = .true.
     savevar => D1_fluid_cS
  end if

  if (varname=='DA_fluid_cS') then
     exists = .true.
     savevar => DA_fluid_cS
  end if

  if (varname=='fluid_q') then
     exists = .true.
     savevar => fluid_q
  end if

  if (varname=='fluid_restmass') then
     exists = .true.
     savevar => fluid_restmass
  end if

  if (varname=='fluid_totalmass') then
     exists = .true.
     savevar => fluid_totalmass
  end if

  if (varname=='cosmobg_tau') then
     exists = .true.
     savevar0D => cosmobg_tau
     nullify(savevar)
  end if

  if (varname=='cosmobg_H') then
     exists = .true.
     savevar0D => cosmobg_H
     nullify(savevar)
  end if

  if (varname=='cosmobg_a') then
     exists = .true.
     savevar0D => cosmobg_a
     nullify(savevar)
  end if

  if (varname=='cosmobg_phi') then
     exists = .true.
     savevar0D => cosmobg_phi
     nullify(savevar)
  end if

  if (varname=='cosmobg_trK') then
     exists = .true.
     savevar0D => cosmobg_trK
     nullify(savevar)
  end if

  if (varname=='cosmobg_alpha') then
     exists = .true.
     savevar0D => cosmobg_alpha
     nullify(savevar)
  end if

  if (varname=='cosmobg_falpha') then
     exists = .true.
     savevar0D => cosmobg_falpha
     nullify(savevar)
  end if

  if (varname=='cosmobg_vlight') then
     exists = .true.
     savevar0D => cosmobg_vlight
     nullify(savevar)
  end if

  if (varname=='cosmobg_rho') then
     exists = .true.
     savevar0D => cosmobg_rho
     nullify(savevar)
  end if

  if (varname=='cosmobg_p') then
     exists = .true.
     savevar0D => cosmobg_p
     nullify(savevar)
  end if

  if (varname=='alpha_pert') then
     exists = .true.
     savevar => alpha_pert
  end if

  if (varname=='phi_pert') then
     exists = .true.
     savevar => phi_pert
  end if

  if (varname=='trK_pert') then
     exists = .true.
     savevar => trK_pert
  end if

  if (varname=='rho_pert') then
     exists = .true.
     savevar => rho_pert
  end if

  if (varname=='rho_contrast') then
     exists = .true.
     savevar => rho_contrast
  end if

  if (varname=='rho_pert_int') then
     exists = .true.
     savevar => rho_pert_int
  end if

  if (varname=='rho_contrast_int') then
     exists = .true.
     savevar => rho_contrast_int
  end if

  if (varname=='cosmo_compactness') then
     exists = .true.
     savevar => cosmo_compactness
  end if

  if (varname=='cosmobg_scalar_phi') then
     exists = .true.
     savevar0D => cosmobg_scalar_phi
     nullify(savevar)
  end if

  if (varname=='cosmobg_scalar_pi') then
     exists = .true.
     savevar0D => cosmobg_scalar_pi
     nullify(savevar)
  end if

  if (varname=='cosmobg_scalar_V') then
     exists = .true.
     savevar0D => cosmobg_scalar_V
     nullify(savevar)
  end if

  if (varname=='cosmobg_scalar_VP') then
     exists = .true.
     savevar0D => cosmobg_scalar_VP
     nullify(savevar)
  end if

  if (varname=='cosmobg_complex_phiR') then
     exists = .true.
     savevar0D => cosmobg_complex_phiR
     nullify(savevar)
  end if

  if (varname=='cosmobg_complex_phiI') then
     exists = .true.
     savevar0D => cosmobg_complex_phiI
     nullify(savevar)
  end if

  if (varname=='cosmobg_complex_piR') then
     exists = .true.
     savevar0D => cosmobg_complex_piR
     nullify(savevar)
  end if

  if (varname=='cosmobg_complex_piI') then
     exists = .true.
     savevar0D => cosmobg_complex_piI
     nullify(savevar)
  end if

  if (varname=='cosmobg_complex_V') then
     exists = .true.
     savevar0D => cosmobg_complex_V
     nullify(savevar)
  end if

  if (varname=='cosmobg_complex_VPR') then
     exists = .true.
     savevar0D => cosmobg_complex_VPR
     nullify(savevar)
  end if

  if (varname=='cosmobg_complex_VPI') then
     exists = .true.
     savevar0D => cosmobg_complex_VPI
     nullify(savevar)
  end if

  if (varname=='r_trans') then
     exists = .true.
     savevar => r_trans
  end if

  if (varname=='R_MINK') then
     exists = .true.
     savevar => R_MINK
  end if

  if (varname=='R_MINK_P') then
     exists = .true.
     savevar => R_MINK_P
  end if

  if (varname=='R_MINK_S') then
     exists = .true.
     savevar => R_MINK_S
  end if

  if (varname=='T_MINK') then
     exists = .true.
     savevar => T_MINK
  end if

  if (varname=='T_MINK_P') then
     exists = .true.
     savevar => T_MINK_P
  end if

  if (varname=='T_MINK_S') then
     exists = .true.
     savevar => T_MINK_S
  end if

  if (varname=='CSI_SCHWARZ') then
     exists = .true.
     savevar => CSI_SCHWARZ
  end if

  if (varname=='CSI_SCHWARZ_P') then
     exists = .true.
     savevar => CSI_SCHWARZ_P
  end if

  if (varname=='ETA_SCHWARZ') then
     exists = .true.
     savevar => ETA_SCHWARZ
  end if

  if (varname=='ETA_SCHWARZ_P') then
     exists = .true.
     savevar => ETA_SCHWARZ_P
  end if

  if (.not.exists) then
     if (rank==0) then
        print *
        print *, 'Error in parfile, non-existent array: ',varname
        print *, 'Aborting! (subroutine grabarray.f90)'
        print *
     end if
     call die
  end if

  end subroutine grabarray

