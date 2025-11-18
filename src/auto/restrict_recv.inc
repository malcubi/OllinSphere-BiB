! Automatically generated file.  Do not edit!

     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        alpha(l-1,i/2+k) = w(i/2)
     end do

  if (shift/="none") then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        beta(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (shift/="none") then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        dtbeta(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (shift/="none") then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        fdriver(l-1,i/2+k) = w(i/2)
     end do
  end if

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     phi(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     chi(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     A(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     B(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     trK(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     KTA(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     Deltar(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     lambda(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     Klambda(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     lambda2(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     Klambda2(l-1,i/2+k) = w(i/2)
  end do

  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
  do i=imin,Nrl(p)-(ghost+1),2
     z4theta(l-1,i/2+k) = w(i/2)
  end do

  if (mattertype /= "vacuum") then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        P_Kodama(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"scalar")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        scalar_phi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"scalar")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        scalar_xi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"scalar")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        scalar_pi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"ghost")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        ghost_phi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"ghost")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        ghost_xi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"ghost")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        ghost_pi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complexghost_phiR(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complexghost_phiI(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complexghost_xiR(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complexghost_xiI(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complexghost_piR(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexghost")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complexghost_piI(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complex")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complex_phiR(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complex")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complex_phiI(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complex")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complex_xiR(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complex")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complex_xiI(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complex")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complex_piR(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complex")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        complex_piI(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"nonmin")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        nonmin_phi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"nonmin")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        nonmin_xi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"nonmin")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        nonmin_pi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"electric")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        electric(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"electric")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        ePhi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"electric")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        eAr(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"proca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        procaE(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"proca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        procaPhi(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"proca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        procaA(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaE_R(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaE_I(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaXi_R(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaXi_I(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaPhi_R(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaPhi_I(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaA_R(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaA_I(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaB_R(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaB_I(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaG_R(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaG_I(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaL_R(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"complexproca")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        cprocaL_I(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        dirac_FR(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        dirac_FI(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        dirac_GR(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"dirac")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        dirac_GI(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"dust")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        dust_cD(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"dust")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        dust_cE(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"dust")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        dust_cS(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"fluid")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        fluid_cD(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"fluid")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        fluid_cE(l-1,i/2+k) = w(i/2)
     end do
  end if

  if (contains(mattertype,"fluid")) then
     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
     do i=imin,Nrl(p)-(ghost+1),2
        fluid_cS(l-1,i/2+k) = w(i/2)
     end do
  end if

