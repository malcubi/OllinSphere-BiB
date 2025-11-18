! Automatically generated file.  Do not edit!

  if (slicing/='maximal') then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => alpha
        alpha(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (shift/="none") then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => beta
        beta(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (shift/="none") then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dtbeta
        dtbeta(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (shift/="none") then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => fdriver
        fdriver(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => phi
     phi(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => chi
     chi(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => A
     A(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => B
     B(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => trK
     trK(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => KTA
     KTA(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => Deltar
     Deltar(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => lambda
     lambda(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => Klambda
     Klambda(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => lambda2
     lambda2(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => Klambda2
     Klambda2(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  do i=1,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => z4theta
     z4theta(l-1,i/2+1) = interp(l,r0,.true.)
  end do

  if (mattertype /= "vacuum") then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => P_Kodama
        P_Kodama(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"scalar")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => scalar_phi
        scalar_phi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"scalar")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => scalar_xi
        scalar_xi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"scalar")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => scalar_pi
        scalar_pi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"ghost")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => ghost_phi
        ghost_phi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"ghost")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => ghost_xi
        ghost_xi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"ghost")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => ghost_pi
        ghost_pi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_phiR
        complexghost_phiR(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_phiI
        complexghost_phiI(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_xiR
        complexghost_xiR(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_xiI
        complexghost_xiI(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_piR
        complexghost_piR(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_piI
        complexghost_piI(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complex")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_phiR
        complex_phiR(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complex")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_phiI
        complex_phiI(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complex")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_xiR
        complex_xiR(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complex")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_xiI
        complex_xiI(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complex")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_piR
        complex_piR(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complex")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_piI
        complex_piI(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"nonmin")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => nonmin_phi
        nonmin_phi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"nonmin")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => nonmin_xi
        nonmin_xi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"nonmin")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => nonmin_pi
        nonmin_pi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"electric")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => electric
        electric(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"electric")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => ePhi
        ePhi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"electric")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => eAr
        eAr(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"proca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => procaE
        procaE(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"proca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => procaPhi
        procaPhi(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"proca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => procaA
        procaA(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaE_R
        cprocaE_R(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaE_I
        cprocaE_I(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaXi_R
        cprocaXi_R(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaXi_I
        cprocaXi_I(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaPhi_R
        cprocaPhi_R(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaPhi_I
        cprocaPhi_I(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaA_R
        cprocaA_R(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaA_I
        cprocaA_I(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaB_R
        cprocaB_R(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaB_I
        cprocaB_I(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaG_R
        cprocaG_R(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaG_I
        cprocaG_I(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaL_R
        cprocaL_R(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaL_I
        cprocaL_I(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dirac_FR
        dirac_FR(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dirac_FI
        dirac_FI(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dirac_GR
        dirac_GR(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dirac_GI
        dirac_GI(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"dust")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dust_cD
        dust_cD(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"dust")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dust_cE
        dust_cE(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"dust")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dust_cS
        dust_cS(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"fluid")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => fluid_cD
        fluid_cD(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"fluid")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => fluid_cE
        fluid_cE(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

  if (contains(mattertype,"fluid")) then
     do i=1,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => fluid_cS
        fluid_cS(l-1,i/2+1) = interp(l,r0,.true.)
     end do
  end if

