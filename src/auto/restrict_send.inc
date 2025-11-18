! Automatically generated file.  Do not edit!

  if (slicing/='maximal') then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => alpha
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (slicing/='maximal') then
  end if

  if (shift/="none") then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => beta
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (shift/="none") then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dtbeta
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (shift/="none") then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => fdriver
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => phi
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => chi
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => A
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => B
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => trK
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => KTA
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => Deltar
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => lambda
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => Klambda
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => lambda2
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => Klambda2
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  do i=imin,Nr-(ghost+1),2
     r0 = r(l,i) + delta
     interpvar => z4theta
     w(i/2) = interp(l,r0,.true.)
  end do
  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

  if (mattertype /= "vacuum") then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => P_Kodama
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"scalar")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => scalar_phi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"scalar")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => scalar_xi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"scalar")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => scalar_pi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"ghost")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => ghost_phi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"ghost")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => ghost_xi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"ghost")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => ghost_pi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexghost")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_phiR
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexghost")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_phiI
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexghost")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_xiR
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexghost")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_xiI
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexghost")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_piR
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexghost")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complexghost_piI
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complex")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_phiR
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complex")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_phiI
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complex")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_xiR
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complex")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_xiI
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complex")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_piR
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complex")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => complex_piI
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"nonmin")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => nonmin_phi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"nonmin")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => nonmin_xi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"nonmin")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => nonmin_pi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"electric")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => electric
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"electric")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => ePhi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"electric")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => eAr
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"proca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => procaE
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"proca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => procaPhi
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"proca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => procaA
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaE_R
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaE_I
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaXi_R
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaXi_I
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaPhi_R
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaPhi_I
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaA_R
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaA_I
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaB_R
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaB_I
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaG_R
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaG_I
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaL_R
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"complexproca")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => cprocaL_I
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"dirac")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dirac_FR
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"dirac")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dirac_FI
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"dirac")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dirac_GR
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"dirac")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dirac_GI
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"dust")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dust_cD
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"dust")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dust_cE
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"dust")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => dust_cS
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"fluid")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => fluid_cD
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"fluid")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => fluid_cE
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

  if (contains(mattertype,"fluid")) then
     do i=imin,Nr-(ghost+1),2
        r0 = r(l,i) + delta
        interpvar => fluid_cS
        w(i/2) = interp(l,r0,.true.)
     end do
     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
  end if

