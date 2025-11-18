! Automatically generated file.  Do not edit!

  subroutine simpleboundary(l)

  use param
  use arrays

  implicit none

  logical contains

  integer l

  if (shift/="none") then
     if (boundtype=='static') then
        sdtbeta(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sdtbeta(l,Nr) = sdtbeta(l,Nr-1)
     end if
  end if

  if (boundtype=='static') then
     strK(l,Nr) = 0.d0
  else if (boundtype=='flat') then
     strK(l,Nr) = strK(l,Nr-1)
  end if

  if (boundtype=='static') then
     sKTA(l,Nr) = 0.d0
  else if (boundtype=='flat') then
     sKTA(l,Nr) = sKTA(l,Nr-1)
  end if

  if (boundtype=='static') then
     sDeltar(l,Nr) = 0.d0
  else if (boundtype=='flat') then
     sDeltar(l,Nr) = sDeltar(l,Nr-1)
  end if

  if (boundtype=='static') then
     sKlambda(l,Nr) = 0.d0
  else if (boundtype=='flat') then
     sKlambda(l,Nr) = sKlambda(l,Nr-1)
  end if

  if (boundtype=='static') then
     sKlambda2(l,Nr) = 0.d0
  else if (boundtype=='flat') then
     sKlambda2(l,Nr) = sKlambda2(l,Nr-1)
  end if

  if (boundtype=='static') then
     sz4theta(l,Nr) = 0.d0
  else if (boundtype=='flat') then
     sz4theta(l,Nr) = sz4theta(l,Nr-1)
  end if

  if (mattertype /= "vacuum") then
     if (boundtype=='static') then
        sP_Kodama(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sP_Kodama(l,Nr) = sP_Kodama(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"scalar")) then
     if (boundtype=='static') then
        sscalar_pi(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sscalar_pi(l,Nr) = sscalar_pi(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"ghost")) then
     if (boundtype=='static') then
        sghost_pi(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sghost_pi(l,Nr) = sghost_pi(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     if (boundtype=='static') then
        scomplexghost_xiR(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scomplexghost_xiR(l,Nr) = scomplexghost_xiR(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     if (boundtype=='static') then
        scomplexghost_xiI(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scomplexghost_xiI(l,Nr) = scomplexghost_xiI(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     if (boundtype=='static') then
        scomplexghost_piR(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scomplexghost_piR(l,Nr) = scomplexghost_piR(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     if (boundtype=='static') then
        scomplexghost_piI(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scomplexghost_piI(l,Nr) = scomplexghost_piI(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complex")) then
     if (boundtype=='static') then
        scomplex_piR(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scomplex_piR(l,Nr) = scomplex_piR(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complex")) then
     if (boundtype=='static') then
        scomplex_piI(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scomplex_piI(l,Nr) = scomplex_piI(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"nonmin")) then
     if (boundtype=='static') then
        snonmin_pi(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        snonmin_pi(l,Nr) = snonmin_pi(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"electric")) then
     if (boundtype=='static') then
        sePhi(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sePhi(l,Nr) = sePhi(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"electric")) then
     if (boundtype=='static') then
        seAr(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        seAr(l,Nr) = seAr(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"proca")) then
     if (boundtype=='static') then
        sprocaPhi(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sprocaPhi(l,Nr) = sprocaPhi(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"proca")) then
     if (boundtype=='static') then
        sprocaA(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sprocaA(l,Nr) = sprocaA(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (boundtype=='static') then
        scprocaPhi_R(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scprocaPhi_R(l,Nr) = scprocaPhi_R(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (boundtype=='static') then
        scprocaPhi_I(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scprocaPhi_I(l,Nr) = scprocaPhi_I(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (boundtype=='static') then
        scprocaA_R(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scprocaA_R(l,Nr) = scprocaA_R(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (boundtype=='static') then
        scprocaA_I(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scprocaA_I(l,Nr) = scprocaA_I(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (boundtype=='static') then
        scprocaB_R(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scprocaB_R(l,Nr) = scprocaB_R(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (boundtype=='static') then
        scprocaB_I(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scprocaB_I(l,Nr) = scprocaB_I(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (boundtype=='static') then
        scprocaG_R(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scprocaG_R(l,Nr) = scprocaG_R(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (boundtype=='static') then
        scprocaG_I(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scprocaG_I(l,Nr) = scprocaG_I(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (boundtype=='static') then
        scprocaL_R(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scprocaL_R(l,Nr) = scprocaL_R(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (boundtype=='static') then
        scprocaL_I(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        scprocaL_I(l,Nr) = scprocaL_I(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     if (boundtype=='static') then
        sdirac_FR(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sdirac_FR(l,Nr) = sdirac_FR(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     if (boundtype=='static') then
        sdirac_FI(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sdirac_FI(l,Nr) = sdirac_FI(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     if (boundtype=='static') then
        sdirac_GR(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sdirac_GR(l,Nr) = sdirac_GR(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     if (boundtype=='static') then
        sdirac_GI(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sdirac_GI(l,Nr) = sdirac_GI(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"dust")) then
     if (boundtype=='static') then
        sdust_cD(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sdust_cD(l,Nr) = sdust_cD(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"dust")) then
     if (boundtype=='static') then
        sdust_cE(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sdust_cE(l,Nr) = sdust_cE(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"dust")) then
     if (boundtype=='static') then
        sdust_cS(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sdust_cS(l,Nr) = sdust_cS(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"fluid")) then
     if (boundtype=='static') then
        sfluid_cD(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sfluid_cD(l,Nr) = sfluid_cD(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"fluid")) then
     if (boundtype=='static') then
        sfluid_cE(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sfluid_cE(l,Nr) = sfluid_cE(l,Nr-1)
     end if
  end if

  if (contains(mattertype,"fluid")) then
     if (boundtype=='static') then
        sfluid_cS(l,Nr) = 0.d0
     else if (boundtype=='flat') then
        sfluid_cS(l,Nr) = sfluid_cS(l,Nr-1)
     end if
  end if

  end subroutine simpleboundary

