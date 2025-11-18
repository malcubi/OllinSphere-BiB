! Automatically generated file.  Do not edit!

  if (slicing/='maximal') then
     interpvar => alpha
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        alpha(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*alpha_p(l,Nr-i)
     else
        alpha(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*alpha_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*alpha_bound(l,i,2)
     end if
  end if

  if (shift/="none") then
     interpvar => beta
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        beta(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*beta_p(l,Nr-i)
     else
        beta(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*beta_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*beta_bound(l,i,2)
     end if
  end if

  if (shift/="none") then
     interpvar => dtbeta
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        dtbeta(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*dtbeta_p(l,Nr-i)
     else
        dtbeta(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*dtbeta_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*dtbeta_bound(l,i,2)
     end if
  end if

  if (shift/="none") then
     interpvar => fdriver
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        fdriver(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*fdriver_p(l,Nr-i)
     else
        fdriver(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*fdriver_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*fdriver_bound(l,i,2)
     end if
  end if

  interpvar => phi
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     phi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*phi_p(l,Nr-i)
  else
     phi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*phi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*phi_bound(l,i,2)
  end if

  interpvar => chi
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     chi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*chi_p(l,Nr-i)
  else
     chi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*chi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*chi_bound(l,i,2)
  end if

  interpvar => A
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     A(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*A_p(l,Nr-i)
  else
     A(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*A_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*A_bound(l,i,2)
  end if

  interpvar => B
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     B(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*B_p(l,Nr-i)
  else
     B(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*B_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*B_bound(l,i,2)
  end if

  interpvar => trK
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     trK(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*trK_p(l,Nr-i)
  else
     trK(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*trK_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*trK_bound(l,i,2)
  end if

  interpvar => KTA
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     KTA(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*KTA_p(l,Nr-i)
  else
     KTA(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*KTA_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*KTA_bound(l,i,2)
  end if

  interpvar => Deltar
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     Deltar(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*Deltar_p(l,Nr-i)
  else
     Deltar(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*Deltar_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*Deltar_bound(l,i,2)
  end if

  interpvar => lambda
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     lambda(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*lambda_p(l,Nr-i)
  else
     lambda(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*lambda_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*lambda_bound(l,i,2)
  end if

  interpvar => Klambda
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     Klambda(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*Klambda_p(l,Nr-i)
  else
     Klambda(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*Klambda_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*Klambda_bound(l,i,2)
  end if

  interpvar => lambda2
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     lambda2(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*lambda2_p(l,Nr-i)
  else
     lambda2(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*lambda2_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*lambda2_bound(l,i,2)
  end if

  interpvar => Klambda2
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     Klambda2(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*Klambda2_p(l,Nr-i)
  else
     Klambda2(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*Klambda2_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*Klambda2_bound(l,i,2)
  end if

  interpvar => z4theta
  aux1 = interp(l-1,r0,.false.)
  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (border==1) then
     z4theta(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*z4theta_p(l,Nr-i)
  else
     z4theta(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*z4theta_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*z4theta_bound(l,i,2)
  end if

  if (mattertype /= "vacuum") then
     interpvar => P_Kodama
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        P_Kodama(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*P_Kodama_p(l,Nr-i)
     else
        P_Kodama(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*P_Kodama_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*P_Kodama_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"scalar")) then
     interpvar => scalar_phi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        scalar_phi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*scalar_phi_p(l,Nr-i)
     else
        scalar_phi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*scalar_phi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*scalar_phi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"scalar")) then
     interpvar => scalar_xi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        scalar_xi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*scalar_xi_p(l,Nr-i)
     else
        scalar_xi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*scalar_xi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*scalar_xi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"scalar")) then
     interpvar => scalar_pi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        scalar_pi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*scalar_pi_p(l,Nr-i)
     else
        scalar_pi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*scalar_pi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*scalar_pi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"ghost")) then
     interpvar => ghost_phi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        ghost_phi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*ghost_phi_p(l,Nr-i)
     else
        ghost_phi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*ghost_phi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*ghost_phi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"ghost")) then
     interpvar => ghost_xi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        ghost_xi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*ghost_xi_p(l,Nr-i)
     else
        ghost_xi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*ghost_xi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*ghost_xi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"ghost")) then
     interpvar => ghost_pi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        ghost_pi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*ghost_pi_p(l,Nr-i)
     else
        ghost_pi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*ghost_pi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*ghost_pi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     interpvar => complexghost_phiR
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complexghost_phiR(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complexghost_phiR_p(l,Nr-i)
     else
        complexghost_phiR(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complexghost_phiR_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complexghost_phiR_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     interpvar => complexghost_phiI
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complexghost_phiI(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complexghost_phiI_p(l,Nr-i)
     else
        complexghost_phiI(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complexghost_phiI_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complexghost_phiI_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     interpvar => complexghost_xiR
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complexghost_xiR(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complexghost_xiR_p(l,Nr-i)
     else
        complexghost_xiR(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complexghost_xiR_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complexghost_xiR_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     interpvar => complexghost_xiI
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complexghost_xiI(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complexghost_xiI_p(l,Nr-i)
     else
        complexghost_xiI(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complexghost_xiI_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complexghost_xiI_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     interpvar => complexghost_piR
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complexghost_piR(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complexghost_piR_p(l,Nr-i)
     else
        complexghost_piR(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complexghost_piR_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complexghost_piR_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexghost")) then
     interpvar => complexghost_piI
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complexghost_piI(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complexghost_piI_p(l,Nr-i)
     else
        complexghost_piI(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complexghost_piI_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complexghost_piI_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complex")) then
     interpvar => complex_phiR
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complex_phiR(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complex_phiR_p(l,Nr-i)
     else
        complex_phiR(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complex_phiR_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complex_phiR_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complex")) then
     interpvar => complex_phiI
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complex_phiI(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complex_phiI_p(l,Nr-i)
     else
        complex_phiI(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complex_phiI_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complex_phiI_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complex")) then
     interpvar => complex_xiR
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complex_xiR(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complex_xiR_p(l,Nr-i)
     else
        complex_xiR(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complex_xiR_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complex_xiR_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complex")) then
     interpvar => complex_xiI
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complex_xiI(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complex_xiI_p(l,Nr-i)
     else
        complex_xiI(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complex_xiI_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complex_xiI_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complex")) then
     interpvar => complex_piR
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complex_piR(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complex_piR_p(l,Nr-i)
     else
        complex_piR(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complex_piR_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complex_piR_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complex")) then
     interpvar => complex_piI
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        complex_piI(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*complex_piI_p(l,Nr-i)
     else
        complex_piI(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*complex_piI_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*complex_piI_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"nonmin")) then
     interpvar => nonmin_phi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        nonmin_phi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*nonmin_phi_p(l,Nr-i)
     else
        nonmin_phi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*nonmin_phi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*nonmin_phi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"nonmin")) then
     interpvar => nonmin_xi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        nonmin_xi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*nonmin_xi_p(l,Nr-i)
     else
        nonmin_xi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*nonmin_xi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*nonmin_xi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"nonmin")) then
     interpvar => nonmin_pi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        nonmin_pi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*nonmin_pi_p(l,Nr-i)
     else
        nonmin_pi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*nonmin_pi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*nonmin_pi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"electric")) then
     interpvar => electric
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        electric(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*electric_p(l,Nr-i)
     else
        electric(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*electric_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*electric_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"electric")) then
     interpvar => ePhi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        ePhi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*ePhi_p(l,Nr-i)
     else
        ePhi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*ePhi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*ePhi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"electric")) then
     interpvar => eAr
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        eAr(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*eAr_p(l,Nr-i)
     else
        eAr(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*eAr_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*eAr_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"proca")) then
     interpvar => procaE
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        procaE(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*procaE_p(l,Nr-i)
     else
        procaE(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*procaE_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*procaE_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"proca")) then
     interpvar => procaPhi
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        procaPhi(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*procaPhi_p(l,Nr-i)
     else
        procaPhi(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*procaPhi_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*procaPhi_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"proca")) then
     interpvar => procaA
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        procaA(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*procaA_p(l,Nr-i)
     else
        procaA(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*procaA_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*procaA_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaE_R
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaE_R(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaE_R_p(l,Nr-i)
     else
        cprocaE_R(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaE_R_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaE_R_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaE_I
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaE_I(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaE_I_p(l,Nr-i)
     else
        cprocaE_I(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaE_I_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaE_I_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaXi_R
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaXi_R(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaXi_R_p(l,Nr-i)
     else
        cprocaXi_R(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaXi_R_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaXi_R_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaXi_I
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaXi_I(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaXi_I_p(l,Nr-i)
     else
        cprocaXi_I(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaXi_I_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaXi_I_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaPhi_R
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaPhi_R(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaPhi_R_p(l,Nr-i)
     else
        cprocaPhi_R(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaPhi_R_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaPhi_R_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaPhi_I
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaPhi_I(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaPhi_I_p(l,Nr-i)
     else
        cprocaPhi_I(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaPhi_I_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaPhi_I_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaA_R
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaA_R(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaA_R_p(l,Nr-i)
     else
        cprocaA_R(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaA_R_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaA_R_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaA_I
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaA_I(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaA_I_p(l,Nr-i)
     else
        cprocaA_I(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaA_I_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaA_I_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaB_R
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaB_R(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaB_R_p(l,Nr-i)
     else
        cprocaB_R(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaB_R_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaB_R_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaB_I
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaB_I(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaB_I_p(l,Nr-i)
     else
        cprocaB_I(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaB_I_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaB_I_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaG_R
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaG_R(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaG_R_p(l,Nr-i)
     else
        cprocaG_R(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaG_R_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaG_R_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaG_I
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaG_I(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaG_I_p(l,Nr-i)
     else
        cprocaG_I(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaG_I_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaG_I_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaL_R
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaL_R(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaL_R_p(l,Nr-i)
     else
        cprocaL_R(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaL_R_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaL_R_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     interpvar => cprocaL_I
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        cprocaL_I(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*cprocaL_I_p(l,Nr-i)
     else
        cprocaL_I(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*cprocaL_I_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*cprocaL_I_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     interpvar => dirac_FR
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        dirac_FR(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*dirac_FR_p(l,Nr-i)
     else
        dirac_FR(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*dirac_FR_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*dirac_FR_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     interpvar => dirac_FI
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        dirac_FI(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*dirac_FI_p(l,Nr-i)
     else
        dirac_FI(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*dirac_FI_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*dirac_FI_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     interpvar => dirac_GR
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        dirac_GR(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*dirac_GR_p(l,Nr-i)
     else
        dirac_GR(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*dirac_GR_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*dirac_GR_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"dirac")) then
     interpvar => dirac_GI
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        dirac_GI(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*dirac_GI_p(l,Nr-i)
     else
        dirac_GI(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*dirac_GI_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*dirac_GI_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"dust")) then
     interpvar => dust_cD
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        dust_cD(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*dust_cD_p(l,Nr-i)
     else
        dust_cD(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*dust_cD_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*dust_cD_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"dust")) then
     interpvar => dust_cE
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        dust_cE(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*dust_cE_p(l,Nr-i)
     else
        dust_cE(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*dust_cE_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*dust_cE_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"dust")) then
     interpvar => dust_cS
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        dust_cS(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*dust_cS_p(l,Nr-i)
     else
        dust_cS(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*dust_cS_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*dust_cS_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"fluid")) then
     interpvar => fluid_cD
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        fluid_cD(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*fluid_cD_p(l,Nr-i)
     else
        fluid_cD(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*fluid_cD_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*fluid_cD_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"fluid")) then
     interpvar => fluid_cE
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        fluid_cE(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*fluid_cE_p(l,Nr-i)
     else
        fluid_cE(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*fluid_cE_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*fluid_cE_bound(l,i,2)
     end if
  end if

  if (contains(mattertype,"fluid")) then
     interpvar => fluid_cS
     aux1 = interp(l-1,r0,.false.)
     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (border==1) then
        fluid_cS(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*fluid_cS_p(l,Nr-i)
     else
        fluid_cS(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*fluid_cS_bound(l,i,1) &
              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*fluid_cS_bound(l,i,2)
     end if
  end if

