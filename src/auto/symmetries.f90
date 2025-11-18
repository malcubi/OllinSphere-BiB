! Automatically generated file.  Do not edit!

  subroutine symmetries(l)

  use param
  use arrays

  implicit none

  logical contains

  integer i,l

  do i=1,ghost

     alpha(l,1-i) = + alpha(l,i)

     if (shift/="none") then
        beta(l,1-i) = - beta(l,i)
     end if

     if (shift/="none") then
        dtbeta(l,1-i) = - dtbeta(l,i)
     end if

     if (shift/="none") then
        fdriver(l,1-i) = + fdriver(l,i)
     end if

     phi(l,1-i) = + phi(l,i)

     chi(l,1-i) = + chi(l,i)

     A(l,1-i) = + A(l,i)

     B(l,1-i) = + B(l,i)

     trK(l,1-i) = + trK(l,i)

     KTA(l,1-i) = + KTA(l,i)

     Deltar(l,1-i) = - Deltar(l,i)

     lambda(l,1-i) = + lambda(l,i)

     Klambda(l,1-i) = + Klambda(l,i)

     lambda2(l,1-i) = + lambda2(l,i)

     Klambda2(l,1-i) = + Klambda2(l,i)

     z4theta(l,1-i) = + z4theta(l,i)

     if (mattertype /= "vacuum") then
        P_Kodama(l,1-i) = + P_Kodama(l,i)
     end if

     if (contains(mattertype,"scalar")) then
        scalar_phi(l,1-i) = + scalar_phi(l,i)
     end if

     if (contains(mattertype,"scalar")) then
        scalar_xi(l,1-i) = - scalar_xi(l,i)
     end if

     if (contains(mattertype,"scalar")) then
        scalar_pi(l,1-i) = + scalar_pi(l,i)
     end if

     if (contains(mattertype,"ghost")) then
        ghost_phi(l,1-i) = + ghost_phi(l,i)
     end if

     if (contains(mattertype,"ghost")) then
        ghost_xi(l,1-i) = - ghost_xi(l,i)
     end if

     if (contains(mattertype,"ghost")) then
        ghost_pi(l,1-i) = + ghost_pi(l,i)
     end if

     if (contains(mattertype,"complexghost")) then
        complexghost_phiR(l,1-i) = +(1.d0-2.d0*mod(complexghost_l,2))*complexghost_phiR(l,i)
     end if

     if (contains(mattertype,"complexghost")) then
        complexghost_phiI(l,1-i) = +(1.d0-2.d0*mod(complexghost_l,2))*complexghost_phiI(l,i)
     end if

     if (contains(mattertype,"complexghost")) then
        complexghost_xiR(l,1-i) = -(1.d0-2.d0*mod(complexghost_l,2))*complexghost_xiR(l,i)
     end if

     if (contains(mattertype,"complexghost")) then
        complexghost_xiI(l,1-i) = -(1.d0-2.d0*mod(complexghost_l,2))*complexghost_xiI(l,i)
     end if

     if (contains(mattertype,"complexghost")) then
        complexghost_piR(l,1-i) = +(1.d0-2.d0*mod(complexghost_l,2))*complexghost_piR(l,i)
     end if

     if (contains(mattertype,"complexghost")) then
        complexghost_piI(l,1-i) = +(1.d0-2.d0*mod(complexghost_l,2))*complexghost_piI(l,i)
     end if

     if (contains(mattertype,"complex")) then
        complex_phiR(l,1-i) = +(1.d0-2.d0*mod(complex_l,2))*complex_phiR(l,i)
     end if

     if (contains(mattertype,"complex")) then
        complex_phiI(l,1-i) = +(1.d0-2.d0*mod(complex_l,2))*complex_phiI(l,i)
     end if

     if (contains(mattertype,"complex")) then
        complex_xiR(l,1-i) = -(1.d0-2.d0*mod(complex_l,2))*complex_xiR(l,i)
     end if

     if (contains(mattertype,"complex")) then
        complex_xiI(l,1-i) = -(1.d0-2.d0*mod(complex_l,2))*complex_xiI(l,i)
     end if

     if (contains(mattertype,"complex")) then
        complex_piR(l,1-i) = +(1.d0-2.d0*mod(complex_l,2))*complex_piR(l,i)
     end if

     if (contains(mattertype,"complex")) then
        complex_piI(l,1-i) = +(1.d0-2.d0*mod(complex_l,2))*complex_piI(l,i)
     end if

     if (contains(mattertype,"nonmin")) then
        nonmin_phi(l,1-i) = + nonmin_phi(l,i)
     end if

     if (contains(mattertype,"nonmin")) then
        nonmin_xi(l,1-i) = - nonmin_xi(l,i)
     end if

     if (contains(mattertype,"nonmin")) then
        nonmin_pi(l,1-i) = + nonmin_pi(l,i)
     end if

     if (contains(mattertype,"electric")) then
        electric(l,1-i) = - electric(l,i)
     end if

     if (contains(mattertype,"electric")) then
        ePhi(l,1-i) = + ePhi(l,i)
     end if

     if (contains(mattertype,"electric")) then
        eAr(l,1-i) = - eAr(l,i)
     end if

     if (contains(mattertype,"proca")) then
        procaE(l,1-i) = - procaE(l,i)
     end if

     if (contains(mattertype,"proca")) then
        procaPhi(l,1-i) = + procaPhi(l,i)
     end if

     if (contains(mattertype,"proca")) then
        procaA(l,1-i) = - procaA(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaE_R(l,1-i) = +(1.d0-2.d0*mod(cproca_l+1,2))*cprocaE_R(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaE_I(l,1-i) = +(1.d0-2.d0*mod(cproca_l+1,2))*cprocaE_I(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaXi_R(l,1-i) = +(1.d0-2.d0*mod(cproca_l,2))*cprocaXi_R(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaXi_I(l,1-i) = +(1.d0-2.d0*mod(cproca_l,2))*cprocaXi_I(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaPhi_R(l,1-i) = +(1.d0-2.d0*mod(cproca_l,2))*cprocaPhi_R(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaPhi_I(l,1-i) = +(1.d0-2.d0*mod(cproca_l,2))*cprocaPhi_I(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaA_R(l,1-i) = +(1.d0-2.d0*mod(cproca_l+1,2))*cprocaA_R(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaA_I(l,1-i) = +(1.d0-2.d0*mod(cproca_l+1,2))*cprocaA_I(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaB_R(l,1-i) = +(1.d0-2.d0*mod(cproca_l+2,2))*cprocaB_R(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaB_I(l,1-i) = +(1.d0-2.d0*mod(cproca_l+2,2))*cprocaB_I(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaG_R(l,1-i) = -(1.d0-2.d0*mod(cproca_l+2,2))*cprocaG_R(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaG_I(l,1-i) = -(1.d0-2.d0*mod(cproca_l+2,2))*cprocaG_I(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaL_R(l,1-i) = +(1.d0-2.d0*mod(cproca_l+1,2))*cprocaL_R(l,i)
     end if

     if (contains(mattertype,"complexproca")) then
        cprocaL_I(l,1-i) = +(1.d0-2.d0*mod(cproca_l+1,2))*cprocaL_I(l,i)
     end if

     if (contains(mattertype,"dirac")) then
        dirac_FR(l,1-i) = + dirac_FR(l,i)
     end if

     if (contains(mattertype,"dirac")) then
        dirac_FI(l,1-i) = + dirac_FI(l,i)
     end if

     if (contains(mattertype,"dirac")) then
        dirac_GR(l,1-i) = - dirac_GR(l,i)
     end if

     if (contains(mattertype,"dirac")) then
        dirac_GI(l,1-i) = - dirac_GI(l,i)
     end if

     if (contains(mattertype,"dust")) then
        dust_cD(l,1-i) = + dust_cD(l,i)
     end if

     if (contains(mattertype,"dust")) then
        dust_cE(l,1-i) = + dust_cE(l,i)
     end if

     if (contains(mattertype,"dust")) then
        dust_cS(l,1-i) = - dust_cS(l,i)
     end if

     if (contains(mattertype,"fluid")) then
        fluid_cD(l,1-i) = + fluid_cD(l,i)
     end if

     if (contains(mattertype,"fluid")) then
        fluid_cE(l,1-i) = + fluid_cE(l,i)
     end if

     if (contains(mattertype,"fluid")) then
        fluid_cS(l,1-i) = - fluid_cS(l,i)
     end if

  end do

  end subroutine symmetries

