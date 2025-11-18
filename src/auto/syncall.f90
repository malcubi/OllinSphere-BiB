! Automatically generated file.  Do not edit!

  subroutine syncall(l)

  use param
  use arrays

  implicit none

  logical contains

  integer l

  syncvar => alpha(l,:)
  call sync

  if (shift/="none") then
     syncvar => beta(l,:)
     call sync
  end if

  if (shift/="none") then
     syncvar => dtbeta(l,:)
     call sync
  end if

  if (shift/="none") then
     syncvar => fdriver(l,:)
     call sync
  end if

  syncvar => phi(l,:)
  call sync

  syncvar => chi(l,:)
  call sync

  syncvar => A(l,:)
  call sync

  syncvar => B(l,:)
  call sync

  syncvar => trK(l,:)
  call sync

  syncvar => KTA(l,:)
  call sync

  syncvar => Deltar(l,:)
  call sync

  syncvar => lambda(l,:)
  call sync

  syncvar => Klambda(l,:)
  call sync

  syncvar => lambda2(l,:)
  call sync

  syncvar => Klambda2(l,:)
  call sync

  syncvar => z4theta(l,:)
  call sync

  if (mattertype /= "vacuum") then
     syncvar => P_Kodama(l,:)
     call sync
  end if

  if (contains(mattertype,"scalar")) then
     syncvar => scalar_phi(l,:)
     call sync
  end if

  if (contains(mattertype,"scalar")) then
     syncvar => scalar_xi(l,:)
     call sync
  end if

  if (contains(mattertype,"scalar")) then
     syncvar => scalar_pi(l,:)
     call sync
  end if

  if (contains(mattertype,"ghost")) then
     syncvar => ghost_phi(l,:)
     call sync
  end if

  if (contains(mattertype,"ghost")) then
     syncvar => ghost_xi(l,:)
     call sync
  end if

  if (contains(mattertype,"ghost")) then
     syncvar => ghost_pi(l,:)
     call sync
  end if

  if (contains(mattertype,"complexghost")) then
     syncvar => complexghost_phiR(l,:)
     call sync
  end if

  if (contains(mattertype,"complexghost")) then
     syncvar => complexghost_phiI(l,:)
     call sync
  end if

  if (contains(mattertype,"complexghost")) then
     syncvar => complexghost_xiR(l,:)
     call sync
  end if

  if (contains(mattertype,"complexghost")) then
     syncvar => complexghost_xiI(l,:)
     call sync
  end if

  if (contains(mattertype,"complexghost")) then
     syncvar => complexghost_piR(l,:)
     call sync
  end if

  if (contains(mattertype,"complexghost")) then
     syncvar => complexghost_piI(l,:)
     call sync
  end if

  if (contains(mattertype,"complex")) then
     syncvar => complex_phiR(l,:)
     call sync
  end if

  if (contains(mattertype,"complex")) then
     syncvar => complex_phiI(l,:)
     call sync
  end if

  if (contains(mattertype,"complex")) then
     syncvar => complex_xiR(l,:)
     call sync
  end if

  if (contains(mattertype,"complex")) then
     syncvar => complex_xiI(l,:)
     call sync
  end if

  if (contains(mattertype,"complex")) then
     syncvar => complex_piR(l,:)
     call sync
  end if

  if (contains(mattertype,"complex")) then
     syncvar => complex_piI(l,:)
     call sync
  end if

  if (contains(mattertype,"nonmin")) then
     syncvar => nonmin_phi(l,:)
     call sync
  end if

  if (contains(mattertype,"nonmin")) then
     syncvar => nonmin_xi(l,:)
     call sync
  end if

  if (contains(mattertype,"nonmin")) then
     syncvar => nonmin_pi(l,:)
     call sync
  end if

  if (contains(mattertype,"electric")) then
     syncvar => electric(l,:)
     call sync
  end if

  if (contains(mattertype,"electric")) then
     syncvar => ePhi(l,:)
     call sync
  end if

  if (contains(mattertype,"electric")) then
     syncvar => eAr(l,:)
     call sync
  end if

  if (contains(mattertype,"proca")) then
     syncvar => procaE(l,:)
     call sync
  end if

  if (contains(mattertype,"proca")) then
     syncvar => procaPhi(l,:)
     call sync
  end if

  if (contains(mattertype,"proca")) then
     syncvar => procaA(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaE_R(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaE_I(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaXi_R(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaXi_I(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaPhi_R(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaPhi_I(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaA_R(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaA_I(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaB_R(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaB_I(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaG_R(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaG_I(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaL_R(l,:)
     call sync
  end if

  if (contains(mattertype,"complexproca")) then
     syncvar => cprocaL_I(l,:)
     call sync
  end if

  if (contains(mattertype,"dirac")) then
     syncvar => dirac_FR(l,:)
     call sync
  end if

  if (contains(mattertype,"dirac")) then
     syncvar => dirac_FI(l,:)
     call sync
  end if

  if (contains(mattertype,"dirac")) then
     syncvar => dirac_GR(l,:)
     call sync
  end if

  if (contains(mattertype,"dirac")) then
     syncvar => dirac_GI(l,:)
     call sync
  end if

  if (contains(mattertype,"dust")) then
     syncvar => dust_cD(l,:)
     call sync
  end if

  if (contains(mattertype,"dust")) then
     syncvar => dust_cE(l,:)
     call sync
  end if

  if (contains(mattertype,"dust")) then
     syncvar => dust_cS(l,:)
     call sync
  end if

  if (contains(mattertype,"fluid")) then
     syncvar => fluid_cD(l,:)
     call sync
  end if

  if (contains(mattertype,"fluid")) then
     syncvar => fluid_cE(l,:)
     call sync
  end if

  if (contains(mattertype,"fluid")) then
     syncvar => fluid_cS(l,:)
     call sync
  end if

  end subroutine syncall

