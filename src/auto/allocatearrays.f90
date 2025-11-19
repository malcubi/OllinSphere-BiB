! Automatically generated file.  Do not edit!

  subroutine allocatearrays(status)

  use param
  use arrays
  use procinfo

  implicit none
  logical contains
  character(len=*) status

  if (trim(status)=='on') then
     allocate(s(0:Nl-1),t(0:Nl-1),t1(0:Nl-1),t2(0:Nl-1),dr(0:Nl-1),dt(0:Nl-1))
     s = 0; t = 0.d0; dr=0.d0; dt=0d0
  else
     deallocate(t,dr,dt)
  end if

  if (trim(status)=='on') then
     allocate(auxarray(0:Nl-1,1-ghost:Nrmax))
     auxarray = 0.d0
  else
     deallocate(auxarray)
  end if

  if (trim(status)=='on') then
     allocate(auxarray2(0:Nl-1,1-ghost:Nrmax))
     auxarray2 = 0.d0
  else
     deallocate(auxarray2)
  end if

  if (trim(status)=='on') then
     allocate(r(0:Nl-1,1-ghost:Nrmax))
     r = 0.d0
  else
     deallocate(r)
  end if

  if (trim(status)=='on') then
     allocate(r_area(0:Nl-1,1-ghost:Nrmax))
     r_area = 0.d0
  else
     deallocate(r_area)
  end if

  if (trim(status)=='on') then
     allocate(D1_r_area(0:Nl-1,1-ghost:Nrmax))
     D1_r_area = 0.d0
  else
     deallocate(D1_r_area)
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',alpha'
     allocate(alpha(0:Nl-1,1-ghost:Nrmax))
     allocate(salpha(0:Nl-1,1-ghost:Nrmax))
     allocate(alpha_p(0:Nl-1,1-ghost:Nrmax))
     allocate(alpha_a(0:Nl-1,1-ghost:Nrmax))
     allocate(alpha_bound(0:Nl-1,0:ghost-1,0:3))
     alpha   = 0.d0
     salpha  = 0.d0
     alpha_p = 0.d0
     alpha_a = 0.d0
     alpha_bound = 0.d0
  else
     deallocate(alpha,salpha,alpha_p,alpha_a,alpha_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_alpha(0:Nl-1,1-ghost:Nrmax))
     D1_alpha = 0.d0
  else
     deallocate(D1_alpha)
  end if

  if (trim(status)=='on') then
     allocate(D2_alpha(0:Nl-1,1-ghost:Nrmax))
     D2_alpha = 0.d0
  else
     deallocate(D2_alpha)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_alpha(0:Nl-1,1-ghost:Nrmax))
        DA_alpha = 0.d0
     else
        deallocate(DA_alpha)
     end if
  else if (contains(outvars0D,"DA_alpha").or. &
           contains(outvars1D,"DA_alpha")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_alpha has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     allocate(Dcov2_alpha(0:Nl-1,1-ghost:Nrmax))
     Dcov2_alpha = 0.d0
  else
     deallocate(Dcov2_alpha)
  end if

  if (trim(status)=='on') then
     allocate(Lapla_alpha(0:Nl-1,1-ghost:Nrmax))
     Lapla_alpha = 0.d0
  else
     deallocate(Lapla_alpha)
  end if

  if (trim(status)=='on') then
     allocate(DD_alphar(0:Nl-1,1-ghost:Nrmax))
     DD_alphar = 0.d0
  else
     deallocate(DD_alphar)
  end if

  if (trim(status)=='on') then
     allocate(falpha(0:Nl-1,1-ghost:Nrmax))
     falpha = 0.d0
  else
     deallocate(falpha)
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',tau_origin'
     allocate(tau_origin(0:Nl-1))
     allocate(stau_origin(0:Nl-1))
     allocate(tau_origin_p(0:Nl-1))
     allocate(tau_origin_a(0:Nl-1))
     tau_origin   = 0.d0
     stau_origin  = 0.d0
     tau_origin_p = 0.d0
     tau_origin_a = 0.d0
  else
     deallocate(tau_origin,stau_origin,tau_origin_p,tau_origin_a)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',beta'
        allocate(beta(0:Nl-1,1-ghost:Nrmax))
        allocate(sbeta(0:Nl-1,1-ghost:Nrmax))
        allocate(beta_p(0:Nl-1,1-ghost:Nrmax))
        allocate(beta_a(0:Nl-1,1-ghost:Nrmax))
        allocate(beta_bound(0:Nl-1,0:ghost-1,0:3))
        beta   = 0.d0
        sbeta  = 0.d0
        beta_p = 0.d0
        beta_a = 0.d0
        beta_bound = 0.d0
     else
        deallocate(beta,sbeta,beta_p,beta_a,beta_bound)
     end if
  else if (contains(outvars0D,"beta").or. &
           contains(outvars1D,"beta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array beta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(D1_beta(0:Nl-1,1-ghost:Nrmax))
        D1_beta = 0.d0
     else
        deallocate(D1_beta)
     end if
  else if (contains(outvars0D,"D1_beta").or. &
           contains(outvars1D,"D1_beta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_beta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(D2_beta(0:Nl-1,1-ghost:Nrmax))
        D2_beta = 0.d0
     else
        deallocate(D2_beta)
     end if
  else if (contains(outvars0D,"D2_beta").or. &
           contains(outvars1D,"D2_beta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D2_beta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_beta(0:Nl-1,1-ghost:Nrmax))
        DA_beta = 0.d0
     else
        deallocate(DA_beta)
     end if
  else if (contains(outvars0D,"DA_beta").or. &
           contains(outvars1D,"DA_beta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_beta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',dtbeta'
        allocate(dtbeta(0:Nl-1,1-ghost:Nrmax))
        allocate(sdtbeta(0:Nl-1,1-ghost:Nrmax))
        allocate(dtbeta_p(0:Nl-1,1-ghost:Nrmax))
        allocate(dtbeta_a(0:Nl-1,1-ghost:Nrmax))
        allocate(dtbeta_bound(0:Nl-1,0:ghost-1,0:3))
        dtbeta   = 0.d0
        sdtbeta  = 0.d0
        dtbeta_p = 0.d0
        dtbeta_a = 0.d0
        dtbeta_bound = 0.d0
     else
        deallocate(dtbeta,sdtbeta,dtbeta_p,dtbeta_a,dtbeta_bound)
     end if
  else if (contains(outvars0D,"dtbeta").or. &
           contains(outvars1D,"dtbeta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dtbeta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(D1_dtbeta(0:Nl-1,1-ghost:Nrmax))
        D1_dtbeta = 0.d0
     else
        deallocate(D1_dtbeta)
     end if
  else if (contains(outvars0D,"D1_dtbeta").or. &
           contains(outvars1D,"D1_dtbeta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_dtbeta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_dtbeta(0:Nl-1,1-ghost:Nrmax))
        DA_dtbeta = 0.d0
     else
        deallocate(DA_dtbeta)
     end if
  else if (contains(outvars0D,"DA_dtbeta").or. &
           contains(outvars1D,"DA_dtbeta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_dtbeta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DIV_beta(0:Nl-1,1-ghost:Nrmax))
        DIV_beta = 0.d0
     else
        deallocate(DIV_beta)
     end if
  else if (contains(outvars0D,"DIV_beta").or. &
           contains(outvars1D,"DIV_beta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DIV_beta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(D1_DIV_beta(0:Nl-1,1-ghost:Nrmax))
        D1_DIV_beta = 0.d0
     else
        deallocate(D1_DIV_beta)
     end if
  else if (contains(outvars0D,"D1_DIV_beta").or. &
           contains(outvars1D,"D1_DIV_beta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_DIV_beta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DD_beta(0:Nl-1,1-ghost:Nrmax))
        DD_beta = 0.d0
     else
        deallocate(DD_beta)
     end if
  else if (contains(outvars0D,"DD_beta").or. &
           contains(outvars1D,"DD_beta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DD_beta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',fdriver'
        allocate(fdriver(0:Nl-1,1-ghost:Nrmax))
        allocate(sfdriver(0:Nl-1,1-ghost:Nrmax))
        allocate(fdriver_p(0:Nl-1,1-ghost:Nrmax))
        allocate(fdriver_a(0:Nl-1,1-ghost:Nrmax))
        allocate(fdriver_bound(0:Nl-1,0:ghost-1,0:3))
        fdriver   = 0.d0
        sfdriver  = 0.d0
        fdriver_p = 0.d0
        fdriver_a = 0.d0
        fdriver_bound = 0.d0
     else
        deallocate(fdriver,sfdriver,fdriver_p,fdriver_a,fdriver_bound)
     end if
  else if (contains(outvars0D,"fdriver").or. &
           contains(outvars1D,"fdriver")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fdriver has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(D1_fdriver(0:Nl-1,1-ghost:Nrmax))
        D1_fdriver = 0.d0
     else
        deallocate(D1_fdriver)
     end if
  else if (contains(outvars0D,"D1_fdriver").or. &
           contains(outvars1D,"D1_fdriver")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_fdriver has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_fdriver(0:Nl-1,1-ghost:Nrmax))
        DA_fdriver = 0.d0
     else
        deallocate(DA_fdriver)
     end if
  else if (contains(outvars0D,"DA_fdriver").or. &
           contains(outvars1D,"DA_fdriver")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_fdriver has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',phi'
     allocate(phi(0:Nl-1,1-ghost:Nrmax))
     allocate(sphi(0:Nl-1,1-ghost:Nrmax))
     allocate(phi_p(0:Nl-1,1-ghost:Nrmax))
     allocate(phi_a(0:Nl-1,1-ghost:Nrmax))
     allocate(phi_bound(0:Nl-1,0:ghost-1,0:3))
     phi   = 0.d0
     sphi  = 0.d0
     phi_p = 0.d0
     phi_a = 0.d0
     phi_bound = 0.d0
  else
     deallocate(phi,sphi,phi_p,phi_a,phi_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_phi(0:Nl-1,1-ghost:Nrmax))
     D1_phi = 0.d0
  else
     deallocate(D1_phi)
  end if

  if (trim(status)=='on') then
     allocate(D2_phi(0:Nl-1,1-ghost:Nrmax))
     D2_phi = 0.d0
  else
     deallocate(D2_phi)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_phi(0:Nl-1,1-ghost:Nrmax))
        DA_phi = 0.d0
     else
        deallocate(DA_phi)
     end if
  else if (contains(outvars0D,"DA_phi").or. &
           contains(outvars1D,"DA_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     allocate(psi(0:Nl-1,1-ghost:Nrmax))
     psi = 0.d0
  else
     deallocate(psi)
  end if

  if (trim(status)=='on') then
     allocate(psi2(0:Nl-1,1-ghost:Nrmax))
     psi2 = 0.d0
  else
     deallocate(psi2)
  end if

  if (trim(status)=='on') then
     allocate(psi4(0:Nl-1,1-ghost:Nrmax))
     psi4 = 0.d0
  else
     deallocate(psi4)
  end if

  if (trim(status)=='on') then
     allocate(D1_psi(0:Nl-1,1-ghost:Nrmax))
     D1_psi = 0.d0
  else
     deallocate(D1_psi)
  end if

  if (trim(status)=='on') then
     allocate(D2_psi(0:Nl-1,1-ghost:Nrmax))
     D2_psi = 0.d0
  else
     deallocate(D2_psi)
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',chi'
     allocate(chi(0:Nl-1,1-ghost:Nrmax))
     allocate(schi(0:Nl-1,1-ghost:Nrmax))
     allocate(chi_p(0:Nl-1,1-ghost:Nrmax))
     allocate(chi_a(0:Nl-1,1-ghost:Nrmax))
     allocate(chi_bound(0:Nl-1,0:ghost-1,0:3))
     chi   = 0.d0
     schi  = 0.d0
     chi_p = 0.d0
     chi_a = 0.d0
     chi_bound = 0.d0
  else
     deallocate(chi,schi,chi_p,chi_a,chi_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_chi(0:Nl-1,1-ghost:Nrmax))
     D1_chi = 0.d0
  else
     deallocate(D1_chi)
  end if

  if (trim(status)=='on') then
     allocate(D2_chi(0:Nl-1,1-ghost:Nrmax))
     D2_chi = 0.d0
  else
     deallocate(D2_chi)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_chi(0:Nl-1,1-ghost:Nrmax))
        DA_chi = 0.d0
     else
        deallocate(DA_chi)
     end if
  else if (contains(outvars0D,"DA_chi").or. &
           contains(outvars1D,"DA_chi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_chi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     allocate(DD_phir(0:Nl-1,1-ghost:Nrmax))
     DD_phir = 0.d0
  else
     deallocate(DD_phir)
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',A'
     allocate(A(0:Nl-1,1-ghost:Nrmax))
     allocate(sA(0:Nl-1,1-ghost:Nrmax))
     allocate(A_p(0:Nl-1,1-ghost:Nrmax))
     allocate(A_a(0:Nl-1,1-ghost:Nrmax))
     allocate(A_bound(0:Nl-1,0:ghost-1,0:3))
     A   = 0.d0
     sA  = 0.d0
     A_p = 0.d0
     A_a = 0.d0
     A_bound = 0.d0
  else
     deallocate(A,sA,A_p,A_a,A_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_A(0:Nl-1,1-ghost:Nrmax))
     D1_A = 0.d0
  else
     deallocate(D1_A)
  end if

  if (trim(status)=='on') then
     allocate(D2_A(0:Nl-1,1-ghost:Nrmax))
     D2_A = 0.d0
  else
     deallocate(D2_A)
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',B'
     allocate(B(0:Nl-1,1-ghost:Nrmax))
     allocate(sB(0:Nl-1,1-ghost:Nrmax))
     allocate(B_p(0:Nl-1,1-ghost:Nrmax))
     allocate(B_a(0:Nl-1,1-ghost:Nrmax))
     allocate(B_bound(0:Nl-1,0:ghost-1,0:3))
     B   = 0.d0
     sB  = 0.d0
     B_p = 0.d0
     B_a = 0.d0
     B_bound = 0.d0
  else
     deallocate(B,sB,B_p,B_a,B_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_B(0:Nl-1,1-ghost:Nrmax))
     D1_B = 0.d0
  else
     deallocate(D1_B)
  end if

  if (trim(status)=='on') then
     allocate(D2_B(0:Nl-1,1-ghost:Nrmax))
     D2_B = 0.d0
  else
     deallocate(D2_B)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_A(0:Nl-1,1-ghost:Nrmax))
        DA_A = 0.d0
     else
        deallocate(DA_A)
     end if
  else if (contains(outvars0D,"DA_A").or. &
           contains(outvars1D,"DA_A")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_A has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_B(0:Nl-1,1-ghost:Nrmax))
        DA_B = 0.d0
     else
        deallocate(DA_B)
     end if
  else if (contains(outvars0D,"DA_B").or. &
           contains(outvars1D,"DA_B")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_B has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     allocate(DD_Ar(0:Nl-1,1-ghost:Nrmax))
     DD_Ar = 0.d0
  else
     deallocate(DD_Ar)
  end if

  if (trim(status)=='on') then
     allocate(DD_Br(0:Nl-1,1-ghost:Nrmax))
     DD_Br = 0.d0
  else
     deallocate(DD_Br)
  end if

  if (contains(outvars0D,"APHYS").or. &
      contains(outvars1D,"APHYS")) then
     if (trim(status)=='on') then
        allocate(APHYS(0:Nl-1,1-ghost:Nrmax))
        APHYS = 0.d0
     else
        deallocate(APHYS)
     end if
  end if

  if (contains(outvars0D,"BPHYS").or. &
      contains(outvars1D,"BPHYS")) then
     if (trim(status)=='on') then
        allocate(BPHYS(0:Nl-1,1-ghost:Nrmax))
        BPHYS = 0.d0
     else
        deallocate(BPHYS)
     end if
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',AB2'
     allocate(AB2(0:Nl-1,1-ghost:Nrmax))
     AB2 = 0.d0
  else
     deallocate(AB2)
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',D1_AB2'
     allocate(D1_AB2(0:Nl-1,1-ghost:Nrmax))
     D1_AB2 = 0.d0
  else
     deallocate(D1_AB2)
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',D2_AB2'
     allocate(D2_AB2(0:Nl-1,1-ghost:Nrmax))
     D2_AB2 = 0.d0
  else
     deallocate(D2_AB2)
  end if

  if (contains(outvars0D,"GRR").or. &
      contains(outvars1D,"GRR")) then
     if (trim(status)=='on') then
        allocate(GRR(0:Nl-1,1-ghost:Nrmax))
        GRR = 0.d0
     else
        deallocate(GRR)
     end if
  end if

  if (contains(outvars0D,"hembed").or. &
      contains(outvars1D,"hembed")) then
     if (trim(status)=='on') then
        allocate(hembed(0:Nl-1,1-ghost:Nrmax))
        hembed = 0.d0
     else
        deallocate(hembed)
     end if
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',trK'
     allocate(trK(0:Nl-1,1-ghost:Nrmax))
     allocate(strK(0:Nl-1,1-ghost:Nrmax))
     allocate(trK_p(0:Nl-1,1-ghost:Nrmax))
     allocate(trK_a(0:Nl-1,1-ghost:Nrmax))
     allocate(trK_bound(0:Nl-1,0:ghost-1,0:3))
     trK   = 0.d0
     strK  = 0.d0
     trK_p = 0.d0
     trK_a = 0.d0
     trK_bound = 0.d0
  else
     deallocate(trK,strK,trK_p,trK_a,trK_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_trK(0:Nl-1,1-ghost:Nrmax))
     D1_trK = 0.d0
  else
     deallocate(D1_trK)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_trK(0:Nl-1,1-ghost:Nrmax))
        DA_trK = 0.d0
     else
        deallocate(DA_trK)
     end if
  else if (contains(outvars0D,"DA_trK").or. &
           contains(outvars1D,"DA_trK")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_trK has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',KTA'
     allocate(KTA(0:Nl-1,1-ghost:Nrmax))
     allocate(sKTA(0:Nl-1,1-ghost:Nrmax))
     allocate(KTA_p(0:Nl-1,1-ghost:Nrmax))
     allocate(KTA_a(0:Nl-1,1-ghost:Nrmax))
     allocate(KTA_bound(0:Nl-1,0:ghost-1,0:3))
     KTA   = 0.d0
     sKTA  = 0.d0
     KTA_p = 0.d0
     KTA_a = 0.d0
     KTA_bound = 0.d0
  else
     deallocate(KTA,sKTA,KTA_p,KTA_a,KTA_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_KTA(0:Nl-1,1-ghost:Nrmax))
     D1_KTA = 0.d0
  else
     deallocate(D1_KTA)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_KTA(0:Nl-1,1-ghost:Nrmax))
        DA_KTA = 0.d0
     else
        deallocate(DA_KTA)
     end if
  else if (contains(outvars0D,"DA_KTA").or. &
           contains(outvars1D,"DA_KTA")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_KTA has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     allocate(KTB(0:Nl-1,1-ghost:Nrmax))
     KTB = 0.d0
  else
     deallocate(KTB)
  end if

  if (trim(status)=='on') then
     allocate(D1_KTB(0:Nl-1,1-ghost:Nrmax))
     D1_KTB = 0.d0
  else
     deallocate(D1_KTB)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_KTB(0:Nl-1,1-ghost:Nrmax))
        DA_KTB = 0.d0
     else
        deallocate(DA_KTB)
     end if
  else if (contains(outvars0D,"DA_KTB").or. &
           contains(outvars1D,"DA_KTB")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_KTB has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     allocate(KAPHYS(0:Nl-1,1-ghost:Nrmax))
     KAPHYS = 0.d0
  else
     deallocate(KAPHYS)
  end if

  if (trim(status)=='on') then
     allocate(KBPHYS(0:Nl-1,1-ghost:Nrmax))
     KBPHYS = 0.d0
  else
     deallocate(KBPHYS)
  end if

  if (trim(status)=='on') then
     allocate(K2(0:Nl-1,1-ghost:Nrmax))
     K2 = 0.d0
  else
     deallocate(K2)
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',Deltar'
     allocate(Deltar(0:Nl-1,1-ghost:Nrmax))
     allocate(sDeltar(0:Nl-1,1-ghost:Nrmax))
     allocate(Deltar_p(0:Nl-1,1-ghost:Nrmax))
     allocate(Deltar_a(0:Nl-1,1-ghost:Nrmax))
     allocate(Deltar_bound(0:Nl-1,0:ghost-1,0:3))
     Deltar   = 0.d0
     sDeltar  = 0.d0
     Deltar_p = 0.d0
     Deltar_a = 0.d0
     Deltar_bound = 0.d0
  else
     deallocate(Deltar,sDeltar,Deltar_p,Deltar_a,Deltar_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_Deltar(0:Nl-1,1-ghost:Nrmax))
     D1_Deltar = 0.d0
  else
     deallocate(D1_Deltar)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_Deltar(0:Nl-1,1-ghost:Nrmax))
        DA_Deltar = 0.d0
     else
        deallocate(DA_Deltar)
     end if
  else if (contains(outvars0D,"DA_Deltar").or. &
           contains(outvars1D,"DA_Deltar")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_Deltar has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',Deltar0'
     allocate(Deltar0(0:Nl-1,1-ghost:Nrmax))
     Deltar0 = 0.d0
  else
     deallocate(Deltar0)
  end if

  if (trim(status)=='on') then
     allocate(DD_Deltar(0:Nl-1,1-ghost:Nrmax))
     DD_Deltar = 0.d0
  else
     deallocate(DD_Deltar)
  end if

  if (trim(status)=='on') then
     allocate(DeltaAB(0:Nl-1,1-ghost:Nrmax))
     DeltaAB = 0.d0
  else
     deallocate(DeltaAB)
  end if

  if (trim(status)=='on') then
     allocate(D1_DeltaAB(0:Nl-1,1-ghost:Nrmax))
     D1_DeltaAB = 0.d0
  else
     deallocate(D1_DeltaAB)
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',lambda'
     allocate(lambda(0:Nl-1,1-ghost:Nrmax))
     allocate(slambda(0:Nl-1,1-ghost:Nrmax))
     allocate(lambda_p(0:Nl-1,1-ghost:Nrmax))
     allocate(lambda_a(0:Nl-1,1-ghost:Nrmax))
     allocate(lambda_bound(0:Nl-1,0:ghost-1,0:3))
     lambda   = 0.d0
     slambda  = 0.d0
     lambda_p = 0.d0
     lambda_a = 0.d0
     lambda_bound = 0.d0
  else
     deallocate(lambda,slambda,lambda_p,lambda_a,lambda_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_lambda(0:Nl-1,1-ghost:Nrmax))
     D1_lambda = 0.d0
  else
     deallocate(D1_lambda)
  end if

  if (trim(status)=='on') then
     allocate(D2_lambda(0:Nl-1,1-ghost:Nrmax))
     D2_lambda = 0.d0
  else
     deallocate(D2_lambda)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_lambda(0:Nl-1,1-ghost:Nrmax))
        DA_lambda = 0.d0
     else
        deallocate(DA_lambda)
     end if
  else if (contains(outvars0D,"DA_lambda").or. &
           contains(outvars1D,"DA_lambda")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_lambda has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',Klambda'
     allocate(Klambda(0:Nl-1,1-ghost:Nrmax))
     allocate(sKlambda(0:Nl-1,1-ghost:Nrmax))
     allocate(Klambda_p(0:Nl-1,1-ghost:Nrmax))
     allocate(Klambda_a(0:Nl-1,1-ghost:Nrmax))
     allocate(Klambda_bound(0:Nl-1,0:ghost-1,0:3))
     Klambda   = 0.d0
     sKlambda  = 0.d0
     Klambda_p = 0.d0
     Klambda_a = 0.d0
     Klambda_bound = 0.d0
  else
     deallocate(Klambda,sKlambda,Klambda_p,Klambda_a,Klambda_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_Klambda(0:Nl-1,1-ghost:Nrmax))
     D1_Klambda = 0.d0
  else
     deallocate(D1_Klambda)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_Klambda(0:Nl-1,1-ghost:Nrmax))
        DA_Klambda = 0.d0
     else
        deallocate(DA_Klambda)
     end if
  else if (contains(outvars0D,"DA_Klambda").or. &
           contains(outvars1D,"DA_Klambda")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_Klambda has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',lambda2'
     allocate(lambda2(0:Nl-1,1-ghost:Nrmax))
     allocate(slambda2(0:Nl-1,1-ghost:Nrmax))
     allocate(lambda2_p(0:Nl-1,1-ghost:Nrmax))
     allocate(lambda2_a(0:Nl-1,1-ghost:Nrmax))
     allocate(lambda2_bound(0:Nl-1,0:ghost-1,0:3))
     lambda2   = 0.d0
     slambda2  = 0.d0
     lambda2_p = 0.d0
     lambda2_a = 0.d0
     lambda2_bound = 0.d0
  else
     deallocate(lambda2,slambda2,lambda2_p,lambda2_a,lambda2_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_lambda2(0:Nl-1,1-ghost:Nrmax))
     D1_lambda2 = 0.d0
  else
     deallocate(D1_lambda2)
  end if

  if (trim(status)=='on') then
     allocate(D2_lambda2(0:Nl-1,1-ghost:Nrmax))
     D2_lambda2 = 0.d0
  else
     deallocate(D2_lambda2)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_lambda2(0:Nl-1,1-ghost:Nrmax))
        DA_lambda2 = 0.d0
     else
        deallocate(DA_lambda2)
     end if
  else if (contains(outvars0D,"DA_lambda2").or. &
           contains(outvars1D,"DA_lambda2")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_lambda2 has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',Klambda2'
     allocate(Klambda2(0:Nl-1,1-ghost:Nrmax))
     allocate(sKlambda2(0:Nl-1,1-ghost:Nrmax))
     allocate(Klambda2_p(0:Nl-1,1-ghost:Nrmax))
     allocate(Klambda2_a(0:Nl-1,1-ghost:Nrmax))
     allocate(Klambda2_bound(0:Nl-1,0:ghost-1,0:3))
     Klambda2   = 0.d0
     sKlambda2  = 0.d0
     Klambda2_p = 0.d0
     Klambda2_a = 0.d0
     Klambda2_bound = 0.d0
  else
     deallocate(Klambda2,sKlambda2,Klambda2_p,Klambda2_a,Klambda2_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_Klambda2(0:Nl-1,1-ghost:Nrmax))
     D1_Klambda2 = 0.d0
  else
     deallocate(D1_Klambda2)
  end if

  if (shift/="none") then
     if (trim(status)=='on') then
        allocate(DA_Klambda2(0:Nl-1,1-ghost:Nrmax))
        DA_Klambda2 = 0.d0
     else
        deallocate(DA_Klambda2)
     end if
  else if (contains(outvars0D,"DA_Klambda2").or. &
           contains(outvars1D,"DA_Klambda2")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_Klambda2 has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     allocate(RICA(0:Nl-1,1-ghost:Nrmax))
     RICA = 0.d0
  else
     deallocate(RICA)
  end if

  if (trim(status)=='on') then
     allocate(RSCAL(0:Nl-1,1-ghost:Nrmax))
     RSCAL = 0.d0
  else
     deallocate(RSCAL)
  end if

  if (trim(status)=='on') then
     checkvars = trim(checkvars) // ',z4theta'
     allocate(z4theta(0:Nl-1,1-ghost:Nrmax))
     allocate(sz4theta(0:Nl-1,1-ghost:Nrmax))
     allocate(z4theta_p(0:Nl-1,1-ghost:Nrmax))
     allocate(z4theta_a(0:Nl-1,1-ghost:Nrmax))
     allocate(z4theta_bound(0:Nl-1,0:ghost-1,0:3))
     z4theta   = 0.d0
     sz4theta  = 0.d0
     z4theta_p = 0.d0
     z4theta_a = 0.d0
     z4theta_bound = 0.d0
  else
     deallocate(z4theta,sz4theta,z4theta_p,z4theta_a,z4theta_bound)
  end if

  if (trim(status)=='on') then
     allocate(D1_z4theta(0:Nl-1,1-ghost:Nrmax))
     D1_z4theta = 0.d0
  else
     deallocate(D1_z4theta)
  end if

  if ((formulation=="z4c").and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_z4theta(0:Nl-1,1-ghost:Nrmax))
        DA_z4theta = 0.d0
     else
        deallocate(DA_z4theta)
     end if
  else if (contains(outvars0D,"DA_z4theta").or. &
           contains(outvars1D,"DA_z4theta")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_z4theta has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     allocate(EWEYLA(0:Nl-1,1-ghost:Nrmax))
     EWEYLA = 0.d0
  else
     deallocate(EWEYLA)
  end if

  if (trim(status)=='on') then
     allocate(EWEYLB(0:Nl-1,1-ghost:Nrmax))
     EWEYLB = 0.d0
  else
     deallocate(EWEYLB)
  end if

  if (contains(outvars0D,"invariantI").or. &
      contains(outvars1D,"invariantI")) then
     if (trim(status)=='on') then
        allocate(invariantI(0:Nl-1,1-ghost:Nrmax))
        invariantI = 0.d0
     else
        deallocate(invariantI)
     end if
  end if

  if (contains(outvars0D,"invariantJ").or. &
      contains(outvars1D,"invariantJ")) then
     if (trim(status)=='on') then
        allocate(invariantJ(0:Nl-1,1-ghost:Nrmax))
        invariantJ = 0.d0
     else
        deallocate(invariantJ)
     end if
  end if

  if (contains(outvars0D,"Ricci4D").or. &
      contains(outvars1D,"Ricci4D")) then
     if (trim(status)=='on') then
        allocate(Ricci4D(0:Nl-1,1-ghost:Nrmax))
        Ricci4D = 0.d0
     else
        deallocate(Ricci4D)
     end if
  end if

  if (contains(outvars0D,"Ricci4D_E").or. &
      contains(outvars1D,"Ricci4D_E")) then
     if (trim(status)=='on') then
        allocate(Ricci4D_E(0:Nl-1,1-ghost:Nrmax))
        Ricci4D_E = 0.d0
     else
        deallocate(Ricci4D_E)
     end if
  end if

  if (contains(outvars0D,"ham").or. &
      contains(outvars1D,"ham")) then
     if (trim(status)=='on') then
        allocate(ham(0:Nl-1,1-ghost:Nrmax))
        ham = 0.d0
     else
        deallocate(ham)
     end if
  end if

  if (contains(outvars0D,"hamabs").or. &
      contains(outvars1D,"hamabs")) then
     if (trim(status)=='on') then
        allocate(hamabs(0:Nl-1,1-ghost:Nrmax))
        hamabs = 0.d0
     else
        deallocate(hamabs)
     end if
  end if

  if (contains(outvars0D,"mom").or. &
      contains(outvars1D,"mom")) then
     if (trim(status)=='on') then
        allocate(mom(0:Nl-1,1-ghost:Nrmax))
        mom = 0.d0
     else
        deallocate(mom)
     end if
  end if

  if (contains(outvars0D,"momabs").or. &
      contains(outvars1D,"momabs")) then
     if (trim(status)=='on') then
        allocate(momabs(0:Nl-1,1-ghost:Nrmax))
        momabs = 0.d0
     else
        deallocate(momabs)
     end if
  end if

  if (contains(outvars0D,"CDeltar").or. &
      contains(outvars1D,"CDeltar")) then
     if (trim(status)=='on') then
        allocate(CDeltar(0:Nl-1,1-ghost:Nrmax))
        CDeltar = 0.d0
     else
        deallocate(CDeltar)
     end if
  end if

  if (contains(outvars0D,"Clambda").or. &
      contains(outvars1D,"Clambda")) then
     if (trim(status)=='on') then
        allocate(Clambda(0:Nl-1,1-ghost:Nrmax))
        Clambda = 0.d0
     else
        deallocate(Clambda)
     end if
  end if

  if (contains(outvars0D,"CKlambda").or. &
      contains(outvars1D,"CKlambda")) then
     if (trim(status)=='on') then
        allocate(CKlambda(0:Nl-1,1-ghost:Nrmax))
        CKlambda = 0.d0
     else
        deallocate(CKlambda)
     end if
  end if

  if (contains(outvars0D,"virial1").or. &
      contains(outvars1D,"virial1")) then
     if (trim(status)=='on') then
        allocate(virial1(0:Nl-1,1-ghost:Nrmax))
        virial1 = 0.d0
     else
        deallocate(virial1)
     end if
  end if

  if (contains(outvars0D,"virial2").or. &
      contains(outvars1D,"virial2")) then
     if (trim(status)=='on') then
        allocate(virial2(0:Nl-1,1-ghost:Nrmax))
        virial2 = 0.d0
     else
        deallocate(virial2)
     end if
  end if

  if (contains(outvars0D,"expansion").or. &
      contains(outvars1D,"expansion")) then
     if (trim(status)=='on') then
        allocate(expansion(0:Nl-1,1-ghost:Nrmax))
        expansion = 0.d0
     else
        deallocate(expansion)
     end if
  end if

  if (contains(outvars0D,"expansion_in").or. &
      contains(outvars1D,"expansion_in")) then
     if (trim(status)=='on') then
        allocate(expansion_in(0:Nl-1,1-ghost:Nrmax))
        expansion_in = 0.d0
     else
        deallocate(expansion_in)
     end if
  end if

  if (contains(outvars0D,"vlight").or. &
      contains(outvars1D,"vlight")) then
     if (trim(status)=='on') then
        allocate(vlight(0:Nl-1,1-ghost:Nrmax))
        vlight = 0.d0
     else
        deallocate(vlight)
     end if
  end if

  if (contains(outvars0D,"vgauge").or. &
      contains(outvars1D,"vgauge")) then
     if (trim(status)=='on') then
        allocate(vgauge(0:Nl-1,1-ghost:Nrmax))
        vgauge = 0.d0
     else
        deallocate(vgauge)
     end if
  end if

  if (contains(outvars0D,"wp_gauge").or. &
      contains(outvars1D,"wp_gauge")) then
     if (trim(status)=='on') then
        allocate(wp_gauge(0:Nl-1,1-ghost:Nrmax))
        wp_gauge = 0.d0
     else
        deallocate(wp_gauge)
     end if
  end if

  if (contains(outvars0D,"wm_gauge").or. &
      contains(outvars1D,"wm_gauge")) then
     if (trim(status)=='on') then
        allocate(wm_gauge(0:Nl-1,1-ghost:Nrmax))
        wm_gauge = 0.d0
     else
        deallocate(wm_gauge)
     end if
  end if

  if (contains(outvars0D,"wp_eta").or. &
      contains(outvars1D,"wp_eta")) then
     if (trim(status)=='on') then
        allocate(wp_eta(0:Nl-1,1-ghost:Nrmax))
        wp_eta = 0.d0
     else
        deallocate(wp_eta)
     end if
  end if

  if (contains(outvars0D,"wm_eta").or. &
      contains(outvars1D,"wm_eta")) then
     if (trim(status)=='on') then
        allocate(wm_eta(0:Nl-1,1-ghost:Nrmax))
        wm_eta = 0.d0
     else
        deallocate(wm_eta)
     end if
  end if

  if (contains(outvars0D,"wDelta").or. &
      contains(outvars1D,"wDelta")) then
     if (trim(status)=='on') then
        allocate(wDelta(0:Nl-1,1-ghost:Nrmax))
        wDelta = 0.d0
     else
        deallocate(wDelta)
     end if
  end if

  if (mattertype/="vacuum") then
     if (trim(status)=='on') then
        allocate(rho(0:Nl-1,1-ghost:Nrmax))
        rho = 0.d0
     else
        deallocate(rho)
     end if
  else if (contains(outvars0D,"rho").or. &
           contains(outvars1D,"rho")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array rho has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (mattertype/="vacuum") then
     if (trim(status)=='on') then
        allocate(JA(0:Nl-1,1-ghost:Nrmax))
        JA = 0.d0
     else
        deallocate(JA)
     end if
  else if (contains(outvars0D,"JA").or. &
           contains(outvars1D,"JA")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array JA has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (mattertype/="vacuum") then
     if (trim(status)=='on') then
        allocate(SAA(0:Nl-1,1-ghost:Nrmax))
        SAA = 0.d0
     else
        deallocate(SAA)
     end if
  else if (contains(outvars0D,"SAA").or. &
           contains(outvars1D,"SAA")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array SAA has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (mattertype/="vacuum") then
     if (trim(status)=='on') then
        allocate(SBB(0:Nl-1,1-ghost:Nrmax))
        SBB = 0.d0
     else
        deallocate(SBB)
     end if
  else if (contains(outvars0D,"SBB").or. &
           contains(outvars1D,"SBB")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array SBB has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (mattertype/="vacuum") then
     if (trim(status)=='on') then
        allocate(trS(0:Nl-1,1-ghost:Nrmax))
        trS = 0.d0
     else
        deallocate(trS)
     end if
  else if (contains(outvars0D,"trS").or. &
           contains(outvars1D,"trS")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array trS has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (mattertype/="vacuum") then
     if (trim(status)=='on') then
        allocate(SLL(0:Nl-1,1-ghost:Nrmax))
        SLL = 0.d0
     else
        deallocate(SLL)
     end if
  else if (contains(outvars0D,"SLL").or. &
           contains(outvars1D,"SLL")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array SLL has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (mattertype/="vacuum") then
     if (contains(outvars0D,"rho_r2").or. &
         contains(outvars1D,"rho_r2")) then
        if (trim(status)=='on') then
           allocate(rho_r2(0:Nl-1,1-ghost:Nrmax))
           rho_r2 = 0.d0
        else
           deallocate(rho_r2)
        end if
     end if
  else if (contains(outvars0D,"rho_r2").or. &
           contains(outvars1D,"rho_r2")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array rho_r2 has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(outvars0D,"mass_sch").or. &
      contains(outvars1D,"mass_sch")) then
     if (trim(status)=='on') then
        allocate(mass_sch(0:Nl-1,1-ghost:Nrmax))
        mass_sch = 0.d0
     else
        deallocate(mass_sch)
     end if
  end if

  if (mattertype /= "vacuum") then
     if (trim(status)=='on') then
        allocate(mass_int(0:Nl-1,1-ghost:Nrmax))
        mass_int = 0.d0
     else
        deallocate(mass_int)
     end if
  else if (contains(outvars0D,"mass_int").or. &
           contains(outvars1D,"mass_int")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array mass_int has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(mass_rn(0:Nl-1,1-ghost:Nrmax))
        mass_rn = 0.d0
     else
        deallocate(mass_rn)
     end if
  else if (contains(outvars0D,"mass_rn").or. &
           contains(outvars1D,"mass_rn")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array mass_rn has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (mattertype /= "vacuum") then
     if (trim(status)=='on') then
        allocate(Kodama_mass(0:Nl-1,1-ghost:Nrmax))
        Kodama_mass = 0.d0
     else
        deallocate(Kodama_mass)
     end if
  else if (contains(outvars0D,"Kodama_mass").or. &
           contains(outvars1D,"Kodama_mass")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array Kodama_mass has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (mattertype /= "vacuum") then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',P_Kodama'
        allocate(P_Kodama(0:Nl-1,1-ghost:Nrmax))
        allocate(sP_Kodama(0:Nl-1,1-ghost:Nrmax))
        allocate(P_Kodama_p(0:Nl-1,1-ghost:Nrmax))
        allocate(P_Kodama_a(0:Nl-1,1-ghost:Nrmax))
        allocate(P_Kodama_bound(0:Nl-1,0:ghost-1,0:3))
        P_Kodama   = 0.d0
        sP_Kodama  = 0.d0
        P_Kodama_p = 0.d0
        P_Kodama_a = 0.d0
        P_Kodama_bound = 0.d0
     else
        deallocate(P_Kodama,sP_Kodama,P_Kodama_p,P_Kodama_a,P_Kodama_bound)
     end if
  else if (contains(outvars0D,"P_Kodama").or. &
           contains(outvars1D,"P_Kodama")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array P_Kodama has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (mattertype /= "vacuum") then
     if (trim(status)=='on') then
        allocate(rho_int(0:Nl-1,1-ghost:Nrmax))
        rho_int = 0.d0
     else
        deallocate(rho_int)
     end if
  else if (contains(outvars0D,"rho_int").or. &
           contains(outvars1D,"rho_int")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array rho_int has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (mattertype /="vacuum") then
     if (trim(status)=='on') then
        allocate(compactness(0:Nl-1,1-ghost:Nrmax))
        compactness = 0.d0
     else
        deallocate(compactness)
     end if
  else if (contains(outvars0D,"compactness").or. &
           contains(outvars1D,"compactness")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array compactness has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (idata=="trumpetBH") then
     if (contains(outvars0D,"trumpet_u").or. &
         contains(outvars1D,"trumpet_u")) then
        if (trim(status)=='on') then
           allocate(trumpet_u(0:Nl-1,1-ghost:Nrmax))
           trumpet_u = 0.d0
        else
           deallocate(trumpet_u)
        end if
     end if
  else if (contains(outvars0D,"trumpet_u").or. &
           contains(outvars1D,"trumpet_u")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array trumpet_u has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (idata=="trumpetBH") then
     if (contains(outvars0D,"trumpet_v").or. &
         contains(outvars1D,"trumpet_v")) then
        if (trim(status)=='on') then
           allocate(trumpet_v(0:Nl-1,1-ghost:Nrmax))
           trumpet_v = 0.d0
        else
           deallocate(trumpet_v)
        end if
     end if
  else if (contains(outvars0D,"trumpet_v").or. &
           contains(outvars1D,"trumpet_v")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array trumpet_v has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"scalar")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',scalar_phi'
        allocate(scalar_phi(0:Nl-1,1-ghost:Nrmax))
        allocate(sscalar_phi(0:Nl-1,1-ghost:Nrmax))
        allocate(scalar_phi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(scalar_phi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(scalar_phi_bound(0:Nl-1,0:ghost-1,0:3))
        scalar_phi   = 0.d0
        sscalar_phi  = 0.d0
        scalar_phi_p = 0.d0
        scalar_phi_a = 0.d0
        scalar_phi_bound = 0.d0
     else
        deallocate(scalar_phi,sscalar_phi,scalar_phi_p,scalar_phi_a,scalar_phi_bound)
     end if
  else if (contains(outvars0D,"scalar_phi").or. &
           contains(outvars1D,"scalar_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array scalar_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"scalar")) then
     if (trim(status)=='on') then
        allocate(D1_scalar_phi(0:Nl-1,1-ghost:Nrmax))
        D1_scalar_phi = 0.d0
     else
        deallocate(D1_scalar_phi)
     end if
  else if (contains(outvars0D,"D1_scalar_phi").or. &
           contains(outvars1D,"D1_scalar_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_scalar_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"scalar")) then
     if (trim(status)=='on') then
        allocate(D2_scalar_phi(0:Nl-1,1-ghost:Nrmax))
        D2_scalar_phi = 0.d0
     else
        deallocate(D2_scalar_phi)
     end if
  else if (contains(outvars0D,"D2_scalar_phi").or. &
           contains(outvars1D,"D2_scalar_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D2_scalar_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"scalar")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_scalar_phi(0:Nl-1,1-ghost:Nrmax))
        DA_scalar_phi = 0.d0
     else
        deallocate(DA_scalar_phi)
     end if
  else if (contains(outvars0D,"DA_scalar_phi").or. &
           contains(outvars1D,"DA_scalar_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_scalar_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"scalar")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',scalar_xi'
        allocate(scalar_xi(0:Nl-1,1-ghost:Nrmax))
        allocate(sscalar_xi(0:Nl-1,1-ghost:Nrmax))
        allocate(scalar_xi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(scalar_xi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(scalar_xi_bound(0:Nl-1,0:ghost-1,0:3))
        scalar_xi   = 0.d0
        sscalar_xi  = 0.d0
        scalar_xi_p = 0.d0
        scalar_xi_a = 0.d0
        scalar_xi_bound = 0.d0
     else
        deallocate(scalar_xi,sscalar_xi,scalar_xi_p,scalar_xi_a,scalar_xi_bound)
     end if
  else if (contains(outvars0D,"scalar_xi").or. &
           contains(outvars1D,"scalar_xi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array scalar_xi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"scalar")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',scalar_pi'
        allocate(scalar_pi(0:Nl-1,1-ghost:Nrmax))
        allocate(sscalar_pi(0:Nl-1,1-ghost:Nrmax))
        allocate(scalar_pi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(scalar_pi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(scalar_pi_bound(0:Nl-1,0:ghost-1,0:3))
        scalar_pi   = 0.d0
        sscalar_pi  = 0.d0
        scalar_pi_p = 0.d0
        scalar_pi_a = 0.d0
        scalar_pi_bound = 0.d0
     else
        deallocate(scalar_pi,sscalar_pi,scalar_pi_p,scalar_pi_a,scalar_pi_bound)
     end if
  else if (contains(outvars0D,"scalar_pi").or. &
           contains(outvars1D,"scalar_pi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array scalar_pi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"scalar")) then
     if (trim(status)=='on') then
        allocate(D1_scalar_xi(0:Nl-1,1-ghost:Nrmax))
        D1_scalar_xi = 0.d0
     else
        deallocate(D1_scalar_xi)
     end if
  else if (contains(outvars0D,"D1_scalar_xi").or. &
           contains(outvars1D,"D1_scalar_xi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_scalar_xi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"scalar")) then
     if (trim(status)=='on') then
        allocate(D1_scalar_pi(0:Nl-1,1-ghost:Nrmax))
        D1_scalar_pi = 0.d0
     else
        deallocate(D1_scalar_pi)
     end if
  else if (contains(outvars0D,"D1_scalar_pi").or. &
           contains(outvars1D,"D1_scalar_pi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_scalar_pi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"scalar")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_scalar_xi(0:Nl-1,1-ghost:Nrmax))
        DA_scalar_xi = 0.d0
     else
        deallocate(DA_scalar_xi)
     end if
  else if (contains(outvars0D,"DA_scalar_xi").or. &
           contains(outvars1D,"DA_scalar_xi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_scalar_xi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"scalar")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_scalar_pi(0:Nl-1,1-ghost:Nrmax))
        DA_scalar_pi = 0.d0
     else
        deallocate(DA_scalar_pi)
     end if
  else if (contains(outvars0D,"DA_scalar_pi").or. &
           contains(outvars1D,"DA_scalar_pi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_scalar_pi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"scalar")) then
     if (trim(status)=='on') then
        allocate(scalar_V(0:Nl-1,1-ghost:Nrmax))
        scalar_V = 0.d0
     else
        deallocate(scalar_V)
     end if
  else if (contains(outvars0D,"scalar_V").or. &
           contains(outvars1D,"scalar_V")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array scalar_V has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"scalar")) then
     if (trim(status)=='on') then
        allocate(scalar_VP(0:Nl-1,1-ghost:Nrmax))
        scalar_VP = 0.d0
     else
        deallocate(scalar_VP)
     end if
  else if (contains(outvars0D,"scalar_VP").or. &
           contains(outvars1D,"scalar_VP")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array scalar_VP has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(outvars0D,"wp_scalar").or. &
      contains(outvars1D,"wp_scalar")) then
     if (trim(status)=='on') then
        allocate(wp_scalar(0:Nl-1,1-ghost:Nrmax))
        wp_scalar = 0.d0
     else
        deallocate(wp_scalar)
     end if
  end if

  if (contains(outvars0D,"wm_scalar").or. &
      contains(outvars1D,"wm_scalar")) then
     if (trim(status)=='on') then
        allocate(wm_scalar(0:Nl-1,1-ghost:Nrmax))
        wm_scalar = 0.d0
     else
        deallocate(wm_scalar)
     end if
  end if

  if (contains(mattertype,"ghost")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',ghost_phi'
        allocate(ghost_phi(0:Nl-1,1-ghost:Nrmax))
        allocate(sghost_phi(0:Nl-1,1-ghost:Nrmax))
        allocate(ghost_phi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(ghost_phi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(ghost_phi_bound(0:Nl-1,0:ghost-1,0:3))
        ghost_phi   = 0.d0
        sghost_phi  = 0.d0
        ghost_phi_p = 0.d0
        ghost_phi_a = 0.d0
        ghost_phi_bound = 0.d0
     else
        deallocate(ghost_phi,sghost_phi,ghost_phi_p,ghost_phi_a,ghost_phi_bound)
     end if
  else if (contains(outvars0D,"ghost_phi").or. &
           contains(outvars1D,"ghost_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array ghost_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"ghost")) then
     if (trim(status)=='on') then
        allocate(D1_ghost_phi(0:Nl-1,1-ghost:Nrmax))
        D1_ghost_phi = 0.d0
     else
        deallocate(D1_ghost_phi)
     end if
  else if (contains(outvars0D,"D1_ghost_phi").or. &
           contains(outvars1D,"D1_ghost_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_ghost_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"ghost")) then
     if (trim(status)=='on') then
        allocate(D2_ghost_phi(0:Nl-1,1-ghost:Nrmax))
        D2_ghost_phi = 0.d0
     else
        deallocate(D2_ghost_phi)
     end if
  else if (contains(outvars0D,"D2_ghost_phi").or. &
           contains(outvars1D,"D2_ghost_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D2_ghost_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"ghost")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_ghost_phi(0:Nl-1,1-ghost:Nrmax))
        DA_ghost_phi = 0.d0
     else
        deallocate(DA_ghost_phi)
     end if
  else if (contains(outvars0D,"DA_ghost_phi").or. &
           contains(outvars1D,"DA_ghost_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_ghost_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"ghost")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',ghost_xi'
        allocate(ghost_xi(0:Nl-1,1-ghost:Nrmax))
        allocate(sghost_xi(0:Nl-1,1-ghost:Nrmax))
        allocate(ghost_xi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(ghost_xi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(ghost_xi_bound(0:Nl-1,0:ghost-1,0:3))
        ghost_xi   = 0.d0
        sghost_xi  = 0.d0
        ghost_xi_p = 0.d0
        ghost_xi_a = 0.d0
        ghost_xi_bound = 0.d0
     else
        deallocate(ghost_xi,sghost_xi,ghost_xi_p,ghost_xi_a,ghost_xi_bound)
     end if
  else if (contains(outvars0D,"ghost_xi").or. &
           contains(outvars1D,"ghost_xi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array ghost_xi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"ghost")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',ghost_pi'
        allocate(ghost_pi(0:Nl-1,1-ghost:Nrmax))
        allocate(sghost_pi(0:Nl-1,1-ghost:Nrmax))
        allocate(ghost_pi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(ghost_pi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(ghost_pi_bound(0:Nl-1,0:ghost-1,0:3))
        ghost_pi   = 0.d0
        sghost_pi  = 0.d0
        ghost_pi_p = 0.d0
        ghost_pi_a = 0.d0
        ghost_pi_bound = 0.d0
     else
        deallocate(ghost_pi,sghost_pi,ghost_pi_p,ghost_pi_a,ghost_pi_bound)
     end if
  else if (contains(outvars0D,"ghost_pi").or. &
           contains(outvars1D,"ghost_pi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array ghost_pi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"ghost")) then
     if (trim(status)=='on') then
        allocate(D1_ghost_xi(0:Nl-1,1-ghost:Nrmax))
        D1_ghost_xi = 0.d0
     else
        deallocate(D1_ghost_xi)
     end if
  else if (contains(outvars0D,"D1_ghost_xi").or. &
           contains(outvars1D,"D1_ghost_xi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_ghost_xi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"ghost")) then
     if (trim(status)=='on') then
        allocate(D1_ghost_pi(0:Nl-1,1-ghost:Nrmax))
        D1_ghost_pi = 0.d0
     else
        deallocate(D1_ghost_pi)
     end if
  else if (contains(outvars0D,"D1_ghost_pi").or. &
           contains(outvars1D,"D1_ghost_pi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_ghost_pi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"ghost")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_ghost_xi(0:Nl-1,1-ghost:Nrmax))
        DA_ghost_xi = 0.d0
     else
        deallocate(DA_ghost_xi)
     end if
  else if (contains(outvars0D,"DA_ghost_xi").or. &
           contains(outvars1D,"DA_ghost_xi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_ghost_xi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"ghost")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_ghost_pi(0:Nl-1,1-ghost:Nrmax))
        DA_ghost_pi = 0.d0
     else
        deallocate(DA_ghost_pi)
     end if
  else if (contains(outvars0D,"DA_ghost_pi").or. &
           contains(outvars1D,"DA_ghost_pi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_ghost_pi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"ghost")) then
     if (trim(status)=='on') then
        allocate(ghost_V(0:Nl-1,1-ghost:Nrmax))
        ghost_V = 0.d0
     else
        deallocate(ghost_V)
     end if
  else if (contains(outvars0D,"ghost_V").or. &
           contains(outvars1D,"ghost_V")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array ghost_V has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"ghost")) then
     if (trim(status)=='on') then
        allocate(ghost_VP(0:Nl-1,1-ghost:Nrmax))
        ghost_VP = 0.d0
     else
        deallocate(ghost_VP)
     end if
  else if (contains(outvars0D,"ghost_VP").or. &
           contains(outvars1D,"ghost_VP")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array ghost_VP has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complexghost_phiR'
        allocate(complexghost_phiR(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplexghost_phiR(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_phiR_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_phiR_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_phiR_bound(0:Nl-1,0:ghost-1,0:3))
        complexghost_phiR   = 0.d0
        scomplexghost_phiR  = 0.d0
        complexghost_phiR_p = 0.d0
        complexghost_phiR_a = 0.d0
        complexghost_phiR_bound = 0.d0
     else
        deallocate(complexghost_phiR,scomplexghost_phiR,complexghost_phiR_p,complexghost_phiR_a,complexghost_phiR_bound)
     end if
  else if (contains(outvars0D,"complexghost_phiR").or. &
           contains(outvars1D,"complexghost_phiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complexghost_phiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complexghost_phiI'
        allocate(complexghost_phiI(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplexghost_phiI(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_phiI_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_phiI_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_phiI_bound(0:Nl-1,0:ghost-1,0:3))
        complexghost_phiI   = 0.d0
        scomplexghost_phiI  = 0.d0
        complexghost_phiI_p = 0.d0
        complexghost_phiI_a = 0.d0
        complexghost_phiI_bound = 0.d0
     else
        deallocate(complexghost_phiI,scomplexghost_phiI,complexghost_phiI_p,complexghost_phiI_a,complexghost_phiI_bound)
     end if
  else if (contains(outvars0D,"complexghost_phiI").or. &
           contains(outvars1D,"complexghost_phiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complexghost_phiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(D1_complexghost_phiR(0:Nl-1,1-ghost:Nrmax))
        D1_complexghost_phiR = 0.d0
     else
        deallocate(D1_complexghost_phiR)
     end if
  else if (contains(outvars0D,"D1_complexghost_phiR").or. &
           contains(outvars1D,"D1_complexghost_phiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complexghost_phiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(D2_complexghost_phiR(0:Nl-1,1-ghost:Nrmax))
        D2_complexghost_phiR = 0.d0
     else
        deallocate(D2_complexghost_phiR)
     end if
  else if (contains(outvars0D,"D2_complexghost_phiR").or. &
           contains(outvars1D,"D2_complexghost_phiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D2_complexghost_phiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost").and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complexghost_phiR(0:Nl-1,1-ghost:Nrmax))
        DA_complexghost_phiR = 0.d0
     else
        deallocate(DA_complexghost_phiR)
     end if
  else if (contains(outvars0D,"DA_complexghost_phiR").or. &
           contains(outvars1D,"DA_complexghost_phiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complexghost_phiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(D1_complexghost_phiI(0:Nl-1,1-ghost:Nrmax))
        D1_complexghost_phiI = 0.d0
     else
        deallocate(D1_complexghost_phiI)
     end if
  else if (contains(outvars0D,"D1_complexghost_phiI").or. &
           contains(outvars1D,"D1_complexghost_phiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complexghost_phiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(D2_complexghost_phiI(0:Nl-1,1-ghost:Nrmax))
        D2_complexghost_phiI = 0.d0
     else
        deallocate(D2_complexghost_phiI)
     end if
  else if (contains(outvars0D,"D2_complexghost_phiI").or. &
           contains(outvars1D,"D2_complexghost_phiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D2_complexghost_phiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost").and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complexghost_phiI(0:Nl-1,1-ghost:Nrmax))
        DA_complexghost_phiI = 0.d0
     else
        deallocate(DA_complexghost_phiI)
     end if
  else if (contains(outvars0D,"DA_complexghost_phiI").or. &
           contains(outvars1D,"DA_complexghost_phiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complexghost_phiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complexghost_xiR'
        allocate(complexghost_xiR(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplexghost_xiR(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_xiR_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_xiR_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_xiR_bound(0:Nl-1,0:ghost-1,0:3))
        complexghost_xiR   = 0.d0
        scomplexghost_xiR  = 0.d0
        complexghost_xiR_p = 0.d0
        complexghost_xiR_a = 0.d0
        complexghost_xiR_bound = 0.d0
     else
        deallocate(complexghost_xiR,scomplexghost_xiR,complexghost_xiR_p,complexghost_xiR_a,complexghost_xiR_bound)
     end if
  else if (contains(outvars0D,"complexghost_xiR").or. &
           contains(outvars1D,"complexghost_xiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complexghost_xiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complexghost_xiI'
        allocate(complexghost_xiI(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplexghost_xiI(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_xiI_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_xiI_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_xiI_bound(0:Nl-1,0:ghost-1,0:3))
        complexghost_xiI   = 0.d0
        scomplexghost_xiI  = 0.d0
        complexghost_xiI_p = 0.d0
        complexghost_xiI_a = 0.d0
        complexghost_xiI_bound = 0.d0
     else
        deallocate(complexghost_xiI,scomplexghost_xiI,complexghost_xiI_p,complexghost_xiI_a,complexghost_xiI_bound)
     end if
  else if (contains(outvars0D,"complexghost_xiI").or. &
           contains(outvars1D,"complexghost_xiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complexghost_xiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complexghost_piR'
        allocate(complexghost_piR(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplexghost_piR(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_piR_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_piR_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_piR_bound(0:Nl-1,0:ghost-1,0:3))
        complexghost_piR   = 0.d0
        scomplexghost_piR  = 0.d0
        complexghost_piR_p = 0.d0
        complexghost_piR_a = 0.d0
        complexghost_piR_bound = 0.d0
     else
        deallocate(complexghost_piR,scomplexghost_piR,complexghost_piR_p,complexghost_piR_a,complexghost_piR_bound)
     end if
  else if (contains(outvars0D,"complexghost_piR").or. &
           contains(outvars1D,"complexghost_piR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complexghost_piR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complexghost_piI'
        allocate(complexghost_piI(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplexghost_piI(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_piI_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_piI_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complexghost_piI_bound(0:Nl-1,0:ghost-1,0:3))
        complexghost_piI   = 0.d0
        scomplexghost_piI  = 0.d0
        complexghost_piI_p = 0.d0
        complexghost_piI_a = 0.d0
        complexghost_piI_bound = 0.d0
     else
        deallocate(complexghost_piI,scomplexghost_piI,complexghost_piI_p,complexghost_piI_a,complexghost_piI_bound)
     end if
  else if (contains(outvars0D,"complexghost_piI").or. &
           contains(outvars1D,"complexghost_piI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complexghost_piI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(D1_complexghost_xiR(0:Nl-1,1-ghost:Nrmax))
        D1_complexghost_xiR = 0.d0
     else
        deallocate(D1_complexghost_xiR)
     end if
  else if (contains(outvars0D,"D1_complexghost_xiR").or. &
           contains(outvars1D,"D1_complexghost_xiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complexghost_xiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(D1_complexghost_xiI(0:Nl-1,1-ghost:Nrmax))
        D1_complexghost_xiI = 0.d0
     else
        deallocate(D1_complexghost_xiI)
     end if
  else if (contains(outvars0D,"D1_complexghost_xiI").or. &
           contains(outvars1D,"D1_complexghost_xiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complexghost_xiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost").and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complexghost_xiR(0:Nl-1,1-ghost:Nrmax))
        DA_complexghost_xiR = 0.d0
     else
        deallocate(DA_complexghost_xiR)
     end if
  else if (contains(outvars0D,"DA_complexghost_xiR").or. &
           contains(outvars1D,"DA_complexghost_xiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complexghost_xiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost").and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complexghost_xiI(0:Nl-1,1-ghost:Nrmax))
        DA_complexghost_xiI = 0.d0
     else
        deallocate(DA_complexghost_xiI)
     end if
  else if (contains(outvars0D,"DA_complexghost_xiI").or. &
           contains(outvars1D,"DA_complexghost_xiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complexghost_xiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(D1_complexghost_piR(0:Nl-1,1-ghost:Nrmax))
        D1_complexghost_piR = 0.d0
     else
        deallocate(D1_complexghost_piR)
     end if
  else if (contains(outvars0D,"D1_complexghost_piR").or. &
           contains(outvars1D,"D1_complexghost_piR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complexghost_piR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(D1_complexghost_piI(0:Nl-1,1-ghost:Nrmax))
        D1_complexghost_piI = 0.d0
     else
        deallocate(D1_complexghost_piI)
     end if
  else if (contains(outvars0D,"D1_complexghost_piI").or. &
           contains(outvars1D,"D1_complexghost_piI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complexghost_piI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost").and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complexghost_piR(0:Nl-1,1-ghost:Nrmax))
        DA_complexghost_piR = 0.d0
     else
        deallocate(DA_complexghost_piR)
     end if
  else if (contains(outvars0D,"DA_complexghost_piR").or. &
           contains(outvars1D,"DA_complexghost_piR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complexghost_piR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost").and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complexghost_piI(0:Nl-1,1-ghost:Nrmax))
        DA_complexghost_piI = 0.d0
     else
        deallocate(DA_complexghost_piI)
     end if
  else if (contains(outvars0D,"DA_complexghost_piI").or. &
           contains(outvars1D,"DA_complexghost_piI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complexghost_piI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(complexghost_V(0:Nl-1,1-ghost:Nrmax))
        complexghost_V = 0.d0
     else
        deallocate(complexghost_V)
     end if
  else if (contains(outvars0D,"complexghost_V").or. &
           contains(outvars1D,"complexghost_V")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complexghost_V has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(complexghost_VPR(0:Nl-1,1-ghost:Nrmax))
        complexghost_VPR = 0.d0
     else
        deallocate(complexghost_VPR)
     end if
  else if (contains(outvars0D,"complexghost_VPR").or. &
           contains(outvars1D,"complexghost_VPR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complexghost_VPR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (trim(status)=='on') then
        allocate(complexghost_VPI(0:Nl-1,1-ghost:Nrmax))
        complexghost_VPI = 0.d0
     else
        deallocate(complexghost_VPI)
     end if
  else if (contains(outvars0D,"complexghost_VPI").or. &
           contains(outvars1D,"complexghost_VPI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complexghost_VPI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexghost")) then
     if (contains(outvars0D,"complexghost_phi_norm").or. &
         contains(outvars1D,"complexghost_phi_norm")) then
        if (trim(status)=='on') then
           allocate(complexghost_phi_norm(0:Nl-1,1-ghost:Nrmax))
           complexghost_phi_norm = 0.d0
        else
           deallocate(complexghost_phi_norm)
        end if
     end if
  else if (contains(outvars0D,"complexghost_phi_norm").or. &
           contains(outvars1D,"complexghost_phi_norm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complexghost_phi_norm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_rhotot(0:Nl-1,1-ghost:Nrmax))
        complex_rhotot = 0.d0
     else
        deallocate(complex_rhotot)
     end if
  else if (contains(outvars0D,"complex_rhotot").or. &
           contains(outvars1D,"complex_rhotot")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_rhotot has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complex_phiR'
        allocate(complex_phiR(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplex_phiR(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_phiR_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_phiR_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_phiR_bound(0:Nl-1,0:ghost-1,0:3))
        complex_phiR   = 0.d0
        scomplex_phiR  = 0.d0
        complex_phiR_p = 0.d0
        complex_phiR_a = 0.d0
        complex_phiR_bound = 0.d0
     else
        deallocate(complex_phiR,scomplex_phiR,complex_phiR_p,complex_phiR_a,complex_phiR_bound)
     end if
  else if (contains(outvars0D,"complex_phiR").or. &
           contains(outvars1D,"complex_phiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_phiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complex_phiI'
        allocate(complex_phiI(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplex_phiI(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_phiI_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_phiI_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_phiI_bound(0:Nl-1,0:ghost-1,0:3))
        complex_phiI   = 0.d0
        scomplex_phiI  = 0.d0
        complex_phiI_p = 0.d0
        complex_phiI_a = 0.d0
        complex_phiI_bound = 0.d0
     else
        deallocate(complex_phiI,scomplex_phiI,complex_phiI_p,complex_phiI_a,complex_phiI_bound)
     end if
  else if (contains(outvars0D,"complex_phiI").or. &
           contains(outvars1D,"complex_phiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_phiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(D1_complex_phiR(0:Nl-1,1-ghost:Nrmax))
        D1_complex_phiR = 0.d0
     else
        deallocate(D1_complex_phiR)
     end if
  else if (contains(outvars0D,"D1_complex_phiR").or. &
           contains(outvars1D,"D1_complex_phiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complex_phiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(D2_complex_phiR(0:Nl-1,1-ghost:Nrmax))
        D2_complex_phiR = 0.d0
     else
        deallocate(D2_complex_phiR)
     end if
  else if (contains(outvars0D,"D2_complex_phiR").or. &
           contains(outvars1D,"D2_complex_phiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D2_complex_phiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complex")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complex_phiR(0:Nl-1,1-ghost:Nrmax))
        DA_complex_phiR = 0.d0
     else
        deallocate(DA_complex_phiR)
     end if
  else if (contains(outvars0D,"DA_complex_phiR").or. &
           contains(outvars1D,"DA_complex_phiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complex_phiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(D1_complex_phiI(0:Nl-1,1-ghost:Nrmax))
        D1_complex_phiI = 0.d0
     else
        deallocate(D1_complex_phiI)
     end if
  else if (contains(outvars0D,"D1_complex_phiI").or. &
           contains(outvars1D,"D1_complex_phiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complex_phiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(D2_complex_phiI(0:Nl-1,1-ghost:Nrmax))
        D2_complex_phiI = 0.d0
     else
        deallocate(D2_complex_phiI)
     end if
  else if (contains(outvars0D,"D2_complex_phiI").or. &
           contains(outvars1D,"D2_complex_phiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D2_complex_phiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complex")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complex_phiI(0:Nl-1,1-ghost:Nrmax))
        DA_complex_phiI = 0.d0
     else
        deallocate(DA_complex_phiI)
     end if
  else if (contains(outvars0D,"DA_complex_phiI").or. &
           contains(outvars1D,"DA_complex_phiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complex_phiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complex_xiR'
        allocate(complex_xiR(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplex_xiR(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_xiR_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_xiR_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_xiR_bound(0:Nl-1,0:ghost-1,0:3))
        complex_xiR   = 0.d0
        scomplex_xiR  = 0.d0
        complex_xiR_p = 0.d0
        complex_xiR_a = 0.d0
        complex_xiR_bound = 0.d0
     else
        deallocate(complex_xiR,scomplex_xiR,complex_xiR_p,complex_xiR_a,complex_xiR_bound)
     end if
  else if (contains(outvars0D,"complex_xiR").or. &
           contains(outvars1D,"complex_xiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_xiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complex_xiI'
        allocate(complex_xiI(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplex_xiI(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_xiI_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_xiI_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_xiI_bound(0:Nl-1,0:ghost-1,0:3))
        complex_xiI   = 0.d0
        scomplex_xiI  = 0.d0
        complex_xiI_p = 0.d0
        complex_xiI_a = 0.d0
        complex_xiI_bound = 0.d0
     else
        deallocate(complex_xiI,scomplex_xiI,complex_xiI_p,complex_xiI_a,complex_xiI_bound)
     end if
  else if (contains(outvars0D,"complex_xiI").or. &
           contains(outvars1D,"complex_xiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_xiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complex_piR'
        allocate(complex_piR(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplex_piR(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_piR_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_piR_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_piR_bound(0:Nl-1,0:ghost-1,0:3))
        complex_piR   = 0.d0
        scomplex_piR  = 0.d0
        complex_piR_p = 0.d0
        complex_piR_a = 0.d0
        complex_piR_bound = 0.d0
     else
        deallocate(complex_piR,scomplex_piR,complex_piR_p,complex_piR_a,complex_piR_bound)
     end if
  else if (contains(outvars0D,"complex_piR").or. &
           contains(outvars1D,"complex_piR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_piR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',complex_piI'
        allocate(complex_piI(0:Nl-1,1-ghost:Nrmax))
        allocate(scomplex_piI(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_piI_p(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_piI_a(0:Nl-1,1-ghost:Nrmax))
        allocate(complex_piI_bound(0:Nl-1,0:ghost-1,0:3))
        complex_piI   = 0.d0
        scomplex_piI  = 0.d0
        complex_piI_p = 0.d0
        complex_piI_a = 0.d0
        complex_piI_bound = 0.d0
     else
        deallocate(complex_piI,scomplex_piI,complex_piI_p,complex_piI_a,complex_piI_bound)
     end if
  else if (contains(outvars0D,"complex_piI").or. &
           contains(outvars1D,"complex_piI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_piI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_gxiR(0:Nl-1,1-ghost:Nrmax))
        complex_gxiR = 0.d0
     else
        deallocate(complex_gxiR)
     end if
  else if (contains(outvars0D,"complex_gxiR").or. &
           contains(outvars1D,"complex_gxiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_gxiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_gxiI(0:Nl-1,1-ghost:Nrmax))
        complex_gxiI = 0.d0
     else
        deallocate(complex_gxiI)
     end if
  else if (contains(outvars0D,"complex_gxiI").or. &
           contains(outvars1D,"complex_gxiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_gxiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_gpiR(0:Nl-1,1-ghost:Nrmax))
        complex_gpiR = 0.d0
     else
        deallocate(complex_gpiR)
     end if
  else if (contains(outvars0D,"complex_gpiR").or. &
           contains(outvars1D,"complex_gpiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_gpiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_gpiI(0:Nl-1,1-ghost:Nrmax))
        complex_gpiI = 0.d0
     else
        deallocate(complex_gpiI)
     end if
  else if (contains(outvars0D,"complex_gpiI").or. &
           contains(outvars1D,"complex_gpiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_gpiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(D1_complex_xiR(0:Nl-1,1-ghost:Nrmax))
        D1_complex_xiR = 0.d0
     else
        deallocate(D1_complex_xiR)
     end if
  else if (contains(outvars0D,"D1_complex_xiR").or. &
           contains(outvars1D,"D1_complex_xiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complex_xiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(D1_complex_xiI(0:Nl-1,1-ghost:Nrmax))
        D1_complex_xiI = 0.d0
     else
        deallocate(D1_complex_xiI)
     end if
  else if (contains(outvars0D,"D1_complex_xiI").or. &
           contains(outvars1D,"D1_complex_xiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complex_xiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complex")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complex_xiR(0:Nl-1,1-ghost:Nrmax))
        DA_complex_xiR = 0.d0
     else
        deallocate(DA_complex_xiR)
     end if
  else if (contains(outvars0D,"DA_complex_xiR").or. &
           contains(outvars1D,"DA_complex_xiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complex_xiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complex")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complex_xiI(0:Nl-1,1-ghost:Nrmax))
        DA_complex_xiI = 0.d0
     else
        deallocate(DA_complex_xiI)
     end if
  else if (contains(outvars0D,"DA_complex_xiI").or. &
           contains(outvars1D,"DA_complex_xiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complex_xiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(D1_complex_piR(0:Nl-1,1-ghost:Nrmax))
        D1_complex_piR = 0.d0
     else
        deallocate(D1_complex_piR)
     end if
  else if (contains(outvars0D,"D1_complex_piR").or. &
           contains(outvars1D,"D1_complex_piR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complex_piR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(D1_complex_piI(0:Nl-1,1-ghost:Nrmax))
        D1_complex_piI = 0.d0
     else
        deallocate(D1_complex_piI)
     end if
  else if (contains(outvars0D,"D1_complex_piI").or. &
           contains(outvars1D,"D1_complex_piI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_complex_piI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complex")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complex_piR(0:Nl-1,1-ghost:Nrmax))
        DA_complex_piR = 0.d0
     else
        deallocate(DA_complex_piR)
     end if
  else if (contains(outvars0D,"DA_complex_piR").or. &
           contains(outvars1D,"DA_complex_piR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complex_piR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complex")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_complex_piI(0:Nl-1,1-ghost:Nrmax))
        DA_complex_piI = 0.d0
     else
        deallocate(DA_complex_piI)
     end if
  else if (contains(outvars0D,"DA_complex_piI").or. &
           contains(outvars1D,"DA_complex_piI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_complex_piI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_V(0:Nl-1,1-ghost:Nrmax))
        complex_V = 0.d0
     else
        deallocate(complex_V)
     end if
  else if (contains(outvars0D,"complex_V").or. &
           contains(outvars1D,"complex_V")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_V has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_VPR(0:Nl-1,1-ghost:Nrmax))
        complex_VPR = 0.d0
     else
        deallocate(complex_VPR)
     end if
  else if (contains(outvars0D,"complex_VPR").or. &
           contains(outvars1D,"complex_VPR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_VPR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_VPI(0:Nl-1,1-ghost:Nrmax))
        complex_VPI = 0.d0
     else
        deallocate(complex_VPI)
     end if
  else if (contains(outvars0D,"complex_VPI").or. &
           contains(outvars1D,"complex_VPI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_VPI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (contains(outvars0D,"complex_phi_norm").or. &
         contains(outvars1D,"complex_phi_norm")) then
        if (trim(status)=='on') then
           allocate(complex_phi_norm(0:Nl-1,1-ghost:Nrmax))
           complex_phi_norm = 0.d0
        else
           deallocate(complex_phi_norm)
        end if
     end if
  else if (contains(outvars0D,"complex_phi_norm").or. &
           contains(outvars1D,"complex_phi_norm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_phi_norm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_Bdens(0:Nl-1,1-ghost:Nrmax))
        complex_Bdens = 0.d0
     else
        deallocate(complex_Bdens)
     end if
  else if (contains(outvars0D,"complex_Bdens").or. &
           contains(outvars1D,"complex_Bdens")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_Bdens has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (contains(outvars0D,"complex_Bdens_r2").or. &
         contains(outvars1D,"complex_Bdens_r2")) then
        if (trim(status)=='on') then
           allocate(complex_Bdens_r2(0:Nl-1,1-ghost:Nrmax))
           complex_Bdens_r2 = 0.d0
        else
           deallocate(complex_Bdens_r2)
        end if
     end if
  else if (contains(outvars0D,"complex_Bdens_r2").or. &
           contains(outvars1D,"complex_Bdens_r2")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_Bdens_r2 has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_Bflux(0:Nl-1,1-ghost:Nrmax))
        complex_Bflux = 0.d0
     else
        deallocate(complex_Bflux)
     end if
  else if (contains(outvars0D,"complex_Bflux").or. &
           contains(outvars1D,"complex_Bflux")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_Bflux has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (contains(outvars0D,"complex_NB").or. &
         contains(outvars1D,"complex_NB")) then
        if (trim(status)=='on') then
           allocate(complex_NB(0:Nl-1,1-ghost:Nrmax))
           complex_NB = 0.d0
        else
           deallocate(complex_NB)
        end if
     end if
  else if (contains(outvars0D,"complex_NB").or. &
           contains(outvars1D,"complex_NB")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_NB has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (contains(outvars0D,"complex_phiaux").or. &
         contains(outvars1D,"complex_phiaux")) then
        if (trim(status)=='on') then
           allocate(complex_phiaux(0:Nl-1,1-ghost:Nrmax))
           complex_phiaux = 0.d0
        else
           deallocate(complex_phiaux)
        end if
     end if
  else if (contains(outvars0D,"complex_phiaux").or. &
           contains(outvars1D,"complex_phiaux")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_phiaux has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_phiR_l0(0:Nl-1,1-ghost:Nrmax))
        complex_phiR_l0 = 0.d0
     else
        deallocate(complex_phiR_l0)
     end if
  else if (contains(outvars0D,"complex_phiR_l0").or. &
           contains(outvars1D,"complex_phiR_l0")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_phiR_l0 has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_phiR_l1(0:Nl-1,1-ghost:Nrmax))
        complex_phiR_l1 = 0.d0
     else
        deallocate(complex_phiR_l1)
     end if
  else if (contains(outvars0D,"complex_phiR_l1").or. &
           contains(outvars1D,"complex_phiR_l1")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_phiR_l1 has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_phiR_l2(0:Nl-1,1-ghost:Nrmax))
        complex_phiR_l2 = 0.d0
     else
        deallocate(complex_phiR_l2)
     end if
  else if (contains(outvars0D,"complex_phiR_l2").or. &
           contains(outvars1D,"complex_phiR_l2")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_phiR_l2 has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complex")) then
     if (trim(status)=='on') then
        allocate(complex_phiR_l3(0:Nl-1,1-ghost:Nrmax))
        complex_phiR_l3 = 0.d0
     else
        deallocate(complex_phiR_l3)
     end if
  else if (contains(outvars0D,"complex_phiR_l3").or. &
           contains(outvars1D,"complex_phiR_l3")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array complex_phiR_l3 has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',nonmin_phi'
        allocate(nonmin_phi(0:Nl-1,1-ghost:Nrmax))
        allocate(snonmin_phi(0:Nl-1,1-ghost:Nrmax))
        allocate(nonmin_phi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(nonmin_phi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(nonmin_phi_bound(0:Nl-1,0:ghost-1,0:3))
        nonmin_phi   = 0.d0
        snonmin_phi  = 0.d0
        nonmin_phi_p = 0.d0
        nonmin_phi_a = 0.d0
        nonmin_phi_bound = 0.d0
     else
        deallocate(nonmin_phi,snonmin_phi,nonmin_phi_p,nonmin_phi_a,nonmin_phi_bound)
     end if
  else if (contains(outvars0D,"nonmin_phi").or. &
           contains(outvars1D,"nonmin_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array nonmin_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(D1_nonmin_phi(0:Nl-1,1-ghost:Nrmax))
        D1_nonmin_phi = 0.d0
     else
        deallocate(D1_nonmin_phi)
     end if
  else if (contains(outvars0D,"D1_nonmin_phi").or. &
           contains(outvars1D,"D1_nonmin_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_nonmin_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(D2_nonmin_phi(0:Nl-1,1-ghost:Nrmax))
        D2_nonmin_phi = 0.d0
     else
        deallocate(D2_nonmin_phi)
     end if
  else if (contains(outvars0D,"D2_nonmin_phi").or. &
           contains(outvars1D,"D2_nonmin_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D2_nonmin_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"nonmin")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_nonmin_phi(0:Nl-1,1-ghost:Nrmax))
        DA_nonmin_phi = 0.d0
     else
        deallocate(DA_nonmin_phi)
     end if
  else if (contains(outvars0D,"DA_nonmin_phi").or. &
           contains(outvars1D,"DA_nonmin_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_nonmin_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',nonmin_xi'
        allocate(nonmin_xi(0:Nl-1,1-ghost:Nrmax))
        allocate(snonmin_xi(0:Nl-1,1-ghost:Nrmax))
        allocate(nonmin_xi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(nonmin_xi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(nonmin_xi_bound(0:Nl-1,0:ghost-1,0:3))
        nonmin_xi   = 0.d0
        snonmin_xi  = 0.d0
        nonmin_xi_p = 0.d0
        nonmin_xi_a = 0.d0
        nonmin_xi_bound = 0.d0
     else
        deallocate(nonmin_xi,snonmin_xi,nonmin_xi_p,nonmin_xi_a,nonmin_xi_bound)
     end if
  else if (contains(outvars0D,"nonmin_xi").or. &
           contains(outvars1D,"nonmin_xi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array nonmin_xi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',nonmin_pi'
        allocate(nonmin_pi(0:Nl-1,1-ghost:Nrmax))
        allocate(snonmin_pi(0:Nl-1,1-ghost:Nrmax))
        allocate(nonmin_pi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(nonmin_pi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(nonmin_pi_bound(0:Nl-1,0:ghost-1,0:3))
        nonmin_pi   = 0.d0
        snonmin_pi  = 0.d0
        nonmin_pi_p = 0.d0
        nonmin_pi_a = 0.d0
        nonmin_pi_bound = 0.d0
     else
        deallocate(nonmin_pi,snonmin_pi,nonmin_pi_p,nonmin_pi_a,nonmin_pi_bound)
     end if
  else if (contains(outvars0D,"nonmin_pi").or. &
           contains(outvars1D,"nonmin_pi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array nonmin_pi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(D1_nonmin_xi(0:Nl-1,1-ghost:Nrmax))
        D1_nonmin_xi = 0.d0
     else
        deallocate(D1_nonmin_xi)
     end if
  else if (contains(outvars0D,"D1_nonmin_xi").or. &
           contains(outvars1D,"D1_nonmin_xi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_nonmin_xi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(D1_nonmin_pi(0:Nl-1,1-ghost:Nrmax))
        D1_nonmin_pi = 0.d0
     else
        deallocate(D1_nonmin_pi)
     end if
  else if (contains(outvars0D,"D1_nonmin_pi").or. &
           contains(outvars1D,"D1_nonmin_pi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_nonmin_pi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"nonmin")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_nonmin_xi(0:Nl-1,1-ghost:Nrmax))
        DA_nonmin_xi = 0.d0
     else
        deallocate(DA_nonmin_xi)
     end if
  else if (contains(outvars0D,"DA_nonmin_xi").or. &
           contains(outvars1D,"DA_nonmin_xi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_nonmin_xi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"nonmin")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_nonmin_pi(0:Nl-1,1-ghost:Nrmax))
        DA_nonmin_pi = 0.d0
     else
        deallocate(DA_nonmin_pi)
     end if
  else if (contains(outvars0D,"DA_nonmin_pi").or. &
           contains(outvars1D,"DA_nonmin_pi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_nonmin_pi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(nonmin_V(0:Nl-1,1-ghost:Nrmax))
        nonmin_V = 0.d0
     else
        deallocate(nonmin_V)
     end if
  else if (contains(outvars0D,"nonmin_V").or. &
           contains(outvars1D,"nonmin_V")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array nonmin_V has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(nonmin_VP(0:Nl-1,1-ghost:Nrmax))
        nonmin_VP = 0.d0
     else
        deallocate(nonmin_VP)
     end if
  else if (contains(outvars0D,"nonmin_VP").or. &
           contains(outvars1D,"nonmin_VP")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array nonmin_VP has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(nonmin_Q(0:Nl-1,1-ghost:Nrmax))
        nonmin_Q = 0.d0
     else
        deallocate(nonmin_Q)
     end if
  else if (contains(outvars0D,"nonmin_Q").or. &
           contains(outvars1D,"nonmin_Q")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array nonmin_Q has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(nonmin_f(0:Nl-1,1-ghost:Nrmax))
        nonmin_f = 0.d0
     else
        deallocate(nonmin_f)
     end if
  else if (contains(outvars0D,"nonmin_f").or. &
           contains(outvars1D,"nonmin_f")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array nonmin_f has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(nonmin_fp(0:Nl-1,1-ghost:Nrmax))
        nonmin_fp = 0.d0
     else
        deallocate(nonmin_fp)
     end if
  else if (contains(outvars0D,"nonmin_fp").or. &
           contains(outvars1D,"nonmin_fp")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array nonmin_fp has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(nonmin_fpp(0:Nl-1,1-ghost:Nrmax))
        nonmin_fpp = 0.d0
     else
        deallocate(nonmin_fpp)
     end if
  else if (contains(outvars0D,"nonmin_fpp").or. &
           contains(outvars1D,"nonmin_fpp")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array nonmin_fpp has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"nonmin")) then
     if (trim(status)=='on') then
        allocate(nonmin_rhs(0:Nl-1,1-ghost:Nrmax))
        nonmin_rhs = 0.d0
     else
        deallocate(nonmin_rhs)
     end if
  else if (contains(outvars0D,"nonmin_rhs").or. &
           contains(outvars1D,"nonmin_rhs")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array nonmin_rhs has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trim(status)=='on') then
     allocate(DD_nonmin_xir(0:Nl-1,1-ghost:Nrmax))
     DD_nonmin_xir = 0.d0
  else
     deallocate(DD_nonmin_xir)
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',electric'
        allocate(electric(0:Nl-1,1-ghost:Nrmax))
        allocate(selectric(0:Nl-1,1-ghost:Nrmax))
        allocate(electric_p(0:Nl-1,1-ghost:Nrmax))
        allocate(electric_a(0:Nl-1,1-ghost:Nrmax))
        allocate(electric_bound(0:Nl-1,0:ghost-1,0:3))
        electric   = 0.d0
        selectric  = 0.d0
        electric_p = 0.d0
        electric_a = 0.d0
        electric_bound = 0.d0
     else
        deallocate(electric,selectric,electric_p,electric_a,electric_bound)
     end if
  else if (contains(outvars0D,"electric").or. &
           contains(outvars1D,"electric")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array electric has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(D1_electric(0:Nl-1,1-ghost:Nrmax))
        D1_electric = 0.d0
     else
        deallocate(D1_electric)
     end if
  else if (contains(outvars0D,"D1_electric").or. &
           contains(outvars1D,"D1_electric")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_electric has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"electric")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_electric(0:Nl-1,1-ghost:Nrmax))
        DA_electric = 0.d0
     else
        deallocate(DA_electric)
     end if
  else if (contains(outvars0D,"DA_electric").or. &
           contains(outvars1D,"DA_electric")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_electric has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',ePhi'
        allocate(ePhi(0:Nl-1,1-ghost:Nrmax))
        allocate(sePhi(0:Nl-1,1-ghost:Nrmax))
        allocate(ePhi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(ePhi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(ePhi_bound(0:Nl-1,0:ghost-1,0:3))
        ePhi   = 0.d0
        sePhi  = 0.d0
        ePhi_p = 0.d0
        ePhi_a = 0.d0
        ePhi_bound = 0.d0
     else
        deallocate(ePhi,sePhi,ePhi_p,ePhi_a,ePhi_bound)
     end if
  else if (contains(outvars0D,"ePhi").or. &
           contains(outvars1D,"ePhi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array ePhi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(D1_ePhi(0:Nl-1,1-ghost:Nrmax))
        D1_ePhi = 0.d0
     else
        deallocate(D1_ePhi)
     end if
  else if (contains(outvars0D,"D1_ePhi").or. &
           contains(outvars1D,"D1_ePhi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_ePhi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"electric")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_ePhi(0:Nl-1,1-ghost:Nrmax))
        DA_ePhi = 0.d0
     else
        deallocate(DA_ePhi)
     end if
  else if (contains(outvars0D,"DA_ePhi").or. &
           contains(outvars1D,"DA_ePhi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_ePhi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',eAr'
        allocate(eAr(0:Nl-1,1-ghost:Nrmax))
        allocate(seAr(0:Nl-1,1-ghost:Nrmax))
        allocate(eAr_p(0:Nl-1,1-ghost:Nrmax))
        allocate(eAr_a(0:Nl-1,1-ghost:Nrmax))
        allocate(eAr_bound(0:Nl-1,0:ghost-1,0:3))
        eAr   = 0.d0
        seAr  = 0.d0
        eAr_p = 0.d0
        eAr_a = 0.d0
        eAr_bound = 0.d0
     else
        deallocate(eAr,seAr,eAr_p,eAr_a,eAr_bound)
     end if
  else if (contains(outvars0D,"eAr").or. &
           contains(outvars1D,"eAr")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array eAr has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(D1_eAr(0:Nl-1,1-ghost:Nrmax))
        D1_eAr = 0.d0
     else
        deallocate(D1_eAr)
     end if
  else if (contains(outvars0D,"D1_eAr").or. &
           contains(outvars1D,"D1_eAr")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_eAr has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"electric")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_eAr(0:Nl-1,1-ghost:Nrmax))
        DA_eAr = 0.d0
     else
        deallocate(DA_eAr)
     end if
  else if (contains(outvars0D,"DA_eAr").or. &
           contains(outvars1D,"DA_eAr")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_eAr has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(eF(0:Nl-1,1-ghost:Nrmax))
        eF = 0.d0
     else
        deallocate(eF)
     end if
  else if (contains(outvars0D,"eF").or. &
           contains(outvars1D,"eF")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array eF has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(D1_eF(0:Nl-1,1-ghost:Nrmax))
        D1_eF = 0.d0
     else
        deallocate(D1_eF)
     end if
  else if (contains(outvars0D,"D1_eF").or. &
           contains(outvars1D,"D1_eF")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_eF has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(eH(0:Nl-1,1-ghost:Nrmax))
        eH = 0.d0
     else
        deallocate(eH)
     end if
  else if (contains(outvars0D,"eH").or. &
           contains(outvars1D,"eH")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array eH has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(D1_eH(0:Nl-1,1-ghost:Nrmax))
        D1_eH = 0.d0
     else
        deallocate(D1_eH)
     end if
  else if (contains(outvars0D,"D1_eH").or. &
           contains(outvars1D,"D1_eH")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_eH has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(echarge(0:Nl-1,1-ghost:Nrmax))
        echarge = 0.d0
     else
        deallocate(echarge)
     end if
  else if (contains(outvars0D,"echarge").or. &
           contains(outvars1D,"echarge")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array echarge has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(ecurrent(0:Nl-1,1-ghost:Nrmax))
        ecurrent = 0.d0
     else
        deallocate(ecurrent)
     end if
  else if (contains(outvars0D,"ecurrent").or. &
           contains(outvars1D,"ecurrent")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array ecurrent has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(Celectric(0:Nl-1,1-ghost:Nrmax))
        Celectric = 0.d0
     else
        deallocate(Celectric)
     end if
  else if (contains(outvars0D,"Celectric").or. &
           contains(outvars1D,"Celectric")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array Celectric has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(eQ_int(0:Nl-1,1-ghost:Nrmax))
        eQ_int = 0.d0
     else
        deallocate(eQ_int)
     end if
  else if (contains(outvars0D,"eQ_int").or. &
           contains(outvars1D,"eQ_int")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array eQ_int has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"electric")) then
     if (trim(status)=='on') then
        allocate(eQ_surf(0:Nl-1,1-ghost:Nrmax))
        eQ_surf = 0.d0
     else
        deallocate(eQ_surf)
     end if
  else if (contains(outvars0D,"eQ_surf").or. &
           contains(outvars1D,"eQ_surf")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array eQ_surf has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(outvars0D,"wp_electric").or. &
      contains(outvars1D,"wp_electric")) then
     if (trim(status)=='on') then
        allocate(wp_electric(0:Nl-1,1-ghost:Nrmax))
        wp_electric = 0.d0
     else
        deallocate(wp_electric)
     end if
  end if

  if (contains(outvars0D,"wm_electric").or. &
      contains(outvars1D,"wm_electric")) then
     if (trim(status)=='on') then
        allocate(wm_electric(0:Nl-1,1-ghost:Nrmax))
        wm_electric = 0.d0
     else
        deallocate(wm_electric)
     end if
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',procaE'
        allocate(procaE(0:Nl-1,1-ghost:Nrmax))
        allocate(sprocaE(0:Nl-1,1-ghost:Nrmax))
        allocate(procaE_p(0:Nl-1,1-ghost:Nrmax))
        allocate(procaE_a(0:Nl-1,1-ghost:Nrmax))
        allocate(procaE_bound(0:Nl-1,0:ghost-1,0:3))
        procaE   = 0.d0
        sprocaE  = 0.d0
        procaE_p = 0.d0
        procaE_a = 0.d0
        procaE_bound = 0.d0
     else
        deallocate(procaE,sprocaE,procaE_p,procaE_a,procaE_bound)
     end if
  else if (contains(outvars0D,"procaE").or. &
           contains(outvars1D,"procaE")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array procaE has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        allocate(D1_procaE(0:Nl-1,1-ghost:Nrmax))
        D1_procaE = 0.d0
     else
        deallocate(D1_procaE)
     end if
  else if (contains(outvars0D,"D1_procaE").or. &
           contains(outvars1D,"D1_procaE")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_procaE has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"proca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_procaE(0:Nl-1,1-ghost:Nrmax))
        DA_procaE = 0.d0
     else
        deallocate(DA_procaE)
     end if
  else if (contains(outvars0D,"DA_procaE").or. &
           contains(outvars1D,"DA_procaE")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_procaE has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',procaPhi'
        allocate(procaPhi(0:Nl-1,1-ghost:Nrmax))
        allocate(sprocaPhi(0:Nl-1,1-ghost:Nrmax))
        allocate(procaPhi_p(0:Nl-1,1-ghost:Nrmax))
        allocate(procaPhi_a(0:Nl-1,1-ghost:Nrmax))
        allocate(procaPhi_bound(0:Nl-1,0:ghost-1,0:3))
        procaPhi   = 0.d0
        sprocaPhi  = 0.d0
        procaPhi_p = 0.d0
        procaPhi_a = 0.d0
        procaPhi_bound = 0.d0
     else
        deallocate(procaPhi,sprocaPhi,procaPhi_p,procaPhi_a,procaPhi_bound)
     end if
  else if (contains(outvars0D,"procaPhi").or. &
           contains(outvars1D,"procaPhi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array procaPhi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        allocate(D1_procaPhi(0:Nl-1,1-ghost:Nrmax))
        D1_procaPhi = 0.d0
     else
        deallocate(D1_procaPhi)
     end if
  else if (contains(outvars0D,"D1_procaPhi").or. &
           contains(outvars1D,"D1_procaPhi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_procaPhi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"proca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_procaPhi(0:Nl-1,1-ghost:Nrmax))
        DA_procaPhi = 0.d0
     else
        deallocate(DA_procaPhi)
     end if
  else if (contains(outvars0D,"DA_procaPhi").or. &
           contains(outvars1D,"DA_procaPhi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_procaPhi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',procaA'
        allocate(procaA(0:Nl-1,1-ghost:Nrmax))
        allocate(sprocaA(0:Nl-1,1-ghost:Nrmax))
        allocate(procaA_p(0:Nl-1,1-ghost:Nrmax))
        allocate(procaA_a(0:Nl-1,1-ghost:Nrmax))
        allocate(procaA_bound(0:Nl-1,0:ghost-1,0:3))
        procaA   = 0.d0
        sprocaA  = 0.d0
        procaA_p = 0.d0
        procaA_a = 0.d0
        procaA_bound = 0.d0
     else
        deallocate(procaA,sprocaA,procaA_p,procaA_a,procaA_bound)
     end if
  else if (contains(outvars0D,"procaA").or. &
           contains(outvars1D,"procaA")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array procaA has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        allocate(D1_procaA(0:Nl-1,1-ghost:Nrmax))
        D1_procaA = 0.d0
     else
        deallocate(D1_procaA)
     end if
  else if (contains(outvars0D,"D1_procaA").or. &
           contains(outvars1D,"D1_procaA")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_procaA has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"proca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_procaA(0:Nl-1,1-ghost:Nrmax))
        DA_procaA = 0.d0
     else
        deallocate(DA_procaA)
     end if
  else if (contains(outvars0D,"DA_procaA").or. &
           contains(outvars1D,"DA_procaA")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_procaA has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        allocate(procaF(0:Nl-1,1-ghost:Nrmax))
        procaF = 0.d0
     else
        deallocate(procaF)
     end if
  else if (contains(outvars0D,"procaF").or. &
           contains(outvars1D,"procaF")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array procaF has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        allocate(D1_procaF(0:Nl-1,1-ghost:Nrmax))
        D1_procaF = 0.d0
     else
        deallocate(D1_procaF)
     end if
  else if (contains(outvars0D,"D1_procaF").or. &
           contains(outvars1D,"D1_procaF")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_procaF has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        allocate(procaH(0:Nl-1,1-ghost:Nrmax))
        procaH = 0.d0
     else
        deallocate(procaH)
     end if
  else if (contains(outvars0D,"procaH").or. &
           contains(outvars1D,"procaH")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array procaH has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        allocate(D1_procaH(0:Nl-1,1-ghost:Nrmax))
        D1_procaH = 0.d0
     else
        deallocate(D1_procaH)
     end if
  else if (contains(outvars0D,"D1_procaH").or. &
           contains(outvars1D,"D1_procaH")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_procaH has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (trim(status)=='on') then
        allocate(procaV(0:Nl-1,1-ghost:Nrmax))
        procaV = 0.d0
     else
        deallocate(procaV)
     end if
  else if (contains(outvars0D,"procaV").or. &
           contains(outvars1D,"procaV")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array procaV has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"proca")) then
     if (contains(outvars0D,"Cproca").or. &
         contains(outvars1D,"Cproca")) then
        if (trim(status)=='on') then
           allocate(Cproca(0:Nl-1,1-ghost:Nrmax))
           Cproca = 0.d0
        else
           deallocate(Cproca)
        end if
     end if
  else if (contains(outvars0D,"Cproca").or. &
           contains(outvars1D,"Cproca")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array Cproca has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(outvars0D,"wp_proca").or. &
      contains(outvars1D,"wp_proca")) then
     if (trim(status)=='on') then
        allocate(wp_proca(0:Nl-1,1-ghost:Nrmax))
        wp_proca = 0.d0
     else
        deallocate(wp_proca)
     end if
  end if

  if (contains(outvars0D,"wm_proca").or. &
      contains(outvars1D,"wm_proca")) then
     if (trim(status)=='on') then
        allocate(wm_proca(0:Nl-1,1-ghost:Nrmax))
        wm_proca = 0.d0
     else
        deallocate(wm_proca)
     end if
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaE_R'
        allocate(cprocaE_R(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaE_R(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaE_R_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaE_R_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaE_R_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaE_R   = 0.d0
        scprocaE_R  = 0.d0
        cprocaE_R_p = 0.d0
        cprocaE_R_a = 0.d0
        cprocaE_R_bound = 0.d0
     else
        deallocate(cprocaE_R,scprocaE_R,cprocaE_R_p,cprocaE_R_a,cprocaE_R_bound)
     end if
  else if (contains(outvars0D,"cprocaE_R").or. &
           contains(outvars1D,"cprocaE_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaE_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaE_R(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaE_R = 0.d0
     else
        deallocate(D1_cprocaE_R)
     end if
  else if (contains(outvars0D,"D1_cprocaE_R").or. &
           contains(outvars1D,"D1_cprocaE_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaE_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaE_R(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaE_R = 0.d0
     else
        deallocate(DA_cprocaE_R)
     end if
  else if (contains(outvars0D,"DA_cprocaE_R").or. &
           contains(outvars1D,"DA_cprocaE_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaE_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaE_I'
        allocate(cprocaE_I(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaE_I(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaE_I_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaE_I_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaE_I_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaE_I   = 0.d0
        scprocaE_I  = 0.d0
        cprocaE_I_p = 0.d0
        cprocaE_I_a = 0.d0
        cprocaE_I_bound = 0.d0
     else
        deallocate(cprocaE_I,scprocaE_I,cprocaE_I_p,cprocaE_I_a,cprocaE_I_bound)
     end if
  else if (contains(outvars0D,"cprocaE_I").or. &
           contains(outvars1D,"cprocaE_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaE_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaE_I(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaE_I = 0.d0
     else
        deallocate(D1_cprocaE_I)
     end if
  else if (contains(outvars0D,"D1_cprocaE_I").or. &
           contains(outvars1D,"D1_cprocaE_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaE_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaE_I(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaE_I = 0.d0
     else
        deallocate(DA_cprocaE_I)
     end if
  else if (contains(outvars0D,"DA_cprocaE_I").or. &
           contains(outvars1D,"DA_cprocaE_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaE_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaXi_R'
        allocate(cprocaXi_R(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaXi_R(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaXi_R_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaXi_R_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaXi_R_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaXi_R   = 0.d0
        scprocaXi_R  = 0.d0
        cprocaXi_R_p = 0.d0
        cprocaXi_R_a = 0.d0
        cprocaXi_R_bound = 0.d0
     else
        deallocate(cprocaXi_R,scprocaXi_R,cprocaXi_R_p,cprocaXi_R_a,cprocaXi_R_bound)
     end if
  else if (contains(outvars0D,"cprocaXi_R").or. &
           contains(outvars1D,"cprocaXi_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaXi_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaXi_R(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaXi_R = 0.d0
     else
        deallocate(D1_cprocaXi_R)
     end if
  else if (contains(outvars0D,"D1_cprocaXi_R").or. &
           contains(outvars1D,"D1_cprocaXi_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaXi_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaXi_R(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaXi_R = 0.d0
     else
        deallocate(DA_cprocaXi_R)
     end if
  else if (contains(outvars0D,"DA_cprocaXi_R").or. &
           contains(outvars1D,"DA_cprocaXi_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaXi_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaXi_I'
        allocate(cprocaXi_I(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaXi_I(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaXi_I_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaXi_I_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaXi_I_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaXi_I   = 0.d0
        scprocaXi_I  = 0.d0
        cprocaXi_I_p = 0.d0
        cprocaXi_I_a = 0.d0
        cprocaXi_I_bound = 0.d0
     else
        deallocate(cprocaXi_I,scprocaXi_I,cprocaXi_I_p,cprocaXi_I_a,cprocaXi_I_bound)
     end if
  else if (contains(outvars0D,"cprocaXi_I").or. &
           contains(outvars1D,"cprocaXi_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaXi_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaXi_I(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaXi_I = 0.d0
     else
        deallocate(D1_cprocaXi_I)
     end if
  else if (contains(outvars0D,"D1_cprocaXi_I").or. &
           contains(outvars1D,"D1_cprocaXi_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaXi_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaXi_I(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaXi_I = 0.d0
     else
        deallocate(DA_cprocaXi_I)
     end if
  else if (contains(outvars0D,"DA_cprocaXi_I").or. &
           contains(outvars1D,"DA_cprocaXi_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaXi_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaPhi_R'
        allocate(cprocaPhi_R(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaPhi_R(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaPhi_R_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaPhi_R_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaPhi_R_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaPhi_R   = 0.d0
        scprocaPhi_R  = 0.d0
        cprocaPhi_R_p = 0.d0
        cprocaPhi_R_a = 0.d0
        cprocaPhi_R_bound = 0.d0
     else
        deallocate(cprocaPhi_R,scprocaPhi_R,cprocaPhi_R_p,cprocaPhi_R_a,cprocaPhi_R_bound)
     end if
  else if (contains(outvars0D,"cprocaPhi_R").or. &
           contains(outvars1D,"cprocaPhi_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaPhi_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaPhi_R(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaPhi_R = 0.d0
     else
        deallocate(D1_cprocaPhi_R)
     end if
  else if (contains(outvars0D,"D1_cprocaPhi_R").or. &
           contains(outvars1D,"D1_cprocaPhi_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaPhi_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaPhi_R(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaPhi_R = 0.d0
     else
        deallocate(DA_cprocaPhi_R)
     end if
  else if (contains(outvars0D,"DA_cprocaPhi_R").or. &
           contains(outvars1D,"DA_cprocaPhi_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaPhi_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaPhi_I'
        allocate(cprocaPhi_I(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaPhi_I(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaPhi_I_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaPhi_I_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaPhi_I_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaPhi_I   = 0.d0
        scprocaPhi_I  = 0.d0
        cprocaPhi_I_p = 0.d0
        cprocaPhi_I_a = 0.d0
        cprocaPhi_I_bound = 0.d0
     else
        deallocate(cprocaPhi_I,scprocaPhi_I,cprocaPhi_I_p,cprocaPhi_I_a,cprocaPhi_I_bound)
     end if
  else if (contains(outvars0D,"cprocaPhi_I").or. &
           contains(outvars1D,"cprocaPhi_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaPhi_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaPhi_I(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaPhi_I = 0.d0
     else
        deallocate(D1_cprocaPhi_I)
     end if
  else if (contains(outvars0D,"D1_cprocaPhi_I").or. &
           contains(outvars1D,"D1_cprocaPhi_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaPhi_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaPhi_I(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaPhi_I = 0.d0
     else
        deallocate(DA_cprocaPhi_I)
     end if
  else if (contains(outvars0D,"DA_cprocaPhi_I").or. &
           contains(outvars1D,"DA_cprocaPhi_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaPhi_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaA_R'
        allocate(cprocaA_R(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaA_R(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaA_R_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaA_R_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaA_R_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaA_R   = 0.d0
        scprocaA_R  = 0.d0
        cprocaA_R_p = 0.d0
        cprocaA_R_a = 0.d0
        cprocaA_R_bound = 0.d0
     else
        deallocate(cprocaA_R,scprocaA_R,cprocaA_R_p,cprocaA_R_a,cprocaA_R_bound)
     end if
  else if (contains(outvars0D,"cprocaA_R").or. &
           contains(outvars1D,"cprocaA_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaA_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaA_R(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaA_R = 0.d0
     else
        deallocate(D1_cprocaA_R)
     end if
  else if (contains(outvars0D,"D1_cprocaA_R").or. &
           contains(outvars1D,"D1_cprocaA_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaA_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaA_R(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaA_R = 0.d0
     else
        deallocate(DA_cprocaA_R)
     end if
  else if (contains(outvars0D,"DA_cprocaA_R").or. &
           contains(outvars1D,"DA_cprocaA_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaA_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaA_I'
        allocate(cprocaA_I(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaA_I(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaA_I_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaA_I_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaA_I_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaA_I   = 0.d0
        scprocaA_I  = 0.d0
        cprocaA_I_p = 0.d0
        cprocaA_I_a = 0.d0
        cprocaA_I_bound = 0.d0
     else
        deallocate(cprocaA_I,scprocaA_I,cprocaA_I_p,cprocaA_I_a,cprocaA_I_bound)
     end if
  else if (contains(outvars0D,"cprocaA_I").or. &
           contains(outvars1D,"cprocaA_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaA_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaA_I(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaA_I = 0.d0
     else
        deallocate(D1_cprocaA_I)
     end if
  else if (contains(outvars0D,"D1_cprocaA_I").or. &
           contains(outvars1D,"D1_cprocaA_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaA_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaA_I(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaA_I = 0.d0
     else
        deallocate(DA_cprocaA_I)
     end if
  else if (contains(outvars0D,"DA_cprocaA_I").or. &
           contains(outvars1D,"DA_cprocaA_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaA_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(cprocaF_R(0:Nl-1,1-ghost:Nrmax))
        cprocaF_R = 0.d0
     else
        deallocate(cprocaF_R)
     end if
  else if (contains(outvars0D,"cprocaF_R").or. &
           contains(outvars1D,"cprocaF_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaF_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaF_R(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaF_R = 0.d0
     else
        deallocate(D1_cprocaF_R)
     end if
  else if (contains(outvars0D,"D1_cprocaF_R").or. &
           contains(outvars1D,"D1_cprocaF_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaF_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(cprocaF_I(0:Nl-1,1-ghost:Nrmax))
        cprocaF_I = 0.d0
     else
        deallocate(cprocaF_I)
     end if
  else if (contains(outvars0D,"cprocaF_I").or. &
           contains(outvars1D,"cprocaF_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaF_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaF_I(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaF_I = 0.d0
     else
        deallocate(D1_cprocaF_I)
     end if
  else if (contains(outvars0D,"D1_cprocaF_I").or. &
           contains(outvars1D,"D1_cprocaF_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaF_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(cprocaH_R(0:Nl-1,1-ghost:Nrmax))
        cprocaH_R = 0.d0
     else
        deallocate(cprocaH_R)
     end if
  else if (contains(outvars0D,"cprocaH_R").or. &
           contains(outvars1D,"cprocaH_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaH_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaH_R(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaH_R = 0.d0
     else
        deallocate(D1_cprocaH_R)
     end if
  else if (contains(outvars0D,"D1_cprocaH_R").or. &
           contains(outvars1D,"D1_cprocaH_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaH_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(cprocaH_I(0:Nl-1,1-ghost:Nrmax))
        cprocaH_I = 0.d0
     else
        deallocate(cprocaH_I)
     end if
  else if (contains(outvars0D,"cprocaH_I").or. &
           contains(outvars1D,"cprocaH_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaH_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaH_I(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaH_I = 0.d0
     else
        deallocate(D1_cprocaH_I)
     end if
  else if (contains(outvars0D,"D1_cprocaH_I").or. &
           contains(outvars1D,"D1_cprocaH_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaH_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaB_R'
        allocate(cprocaB_R(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaB_R(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaB_R_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaB_R_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaB_R_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaB_R   = 0.d0
        scprocaB_R  = 0.d0
        cprocaB_R_p = 0.d0
        cprocaB_R_a = 0.d0
        cprocaB_R_bound = 0.d0
     else
        deallocate(cprocaB_R,scprocaB_R,cprocaB_R_p,cprocaB_R_a,cprocaB_R_bound)
     end if
  else if (contains(outvars0D,"cprocaB_R").or. &
           contains(outvars1D,"cprocaB_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaB_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaB_R(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaB_R = 0.d0
     else
        deallocate(D1_cprocaB_R)
     end if
  else if (contains(outvars0D,"D1_cprocaB_R").or. &
           contains(outvars1D,"D1_cprocaB_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaB_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D2_cprocaB_R(0:Nl-1,1-ghost:Nrmax))
        D2_cprocaB_R = 0.d0
     else
        deallocate(D2_cprocaB_R)
     end if
  else if (contains(outvars0D,"D2_cprocaB_R").or. &
           contains(outvars1D,"D2_cprocaB_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D2_cprocaB_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaB_R(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaB_R = 0.d0
     else
        deallocate(DA_cprocaB_R)
     end if
  else if (contains(outvars0D,"DA_cprocaB_R").or. &
           contains(outvars1D,"DA_cprocaB_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaB_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaB_I'
        allocate(cprocaB_I(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaB_I(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaB_I_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaB_I_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaB_I_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaB_I   = 0.d0
        scprocaB_I  = 0.d0
        cprocaB_I_p = 0.d0
        cprocaB_I_a = 0.d0
        cprocaB_I_bound = 0.d0
     else
        deallocate(cprocaB_I,scprocaB_I,cprocaB_I_p,cprocaB_I_a,cprocaB_I_bound)
     end if
  else if (contains(outvars0D,"cprocaB_I").or. &
           contains(outvars1D,"cprocaB_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaB_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaB_I(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaB_I = 0.d0
     else
        deallocate(D1_cprocaB_I)
     end if
  else if (contains(outvars0D,"D1_cprocaB_I").or. &
           contains(outvars1D,"D1_cprocaB_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaB_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D2_cprocaB_I(0:Nl-1,1-ghost:Nrmax))
        D2_cprocaB_I = 0.d0
     else
        deallocate(D2_cprocaB_I)
     end if
  else if (contains(outvars0D,"D2_cprocaB_I").or. &
           contains(outvars1D,"D2_cprocaB_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D2_cprocaB_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaB_I(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaB_I = 0.d0
     else
        deallocate(DA_cprocaB_I)
     end if
  else if (contains(outvars0D,"DA_cprocaB_I").or. &
           contains(outvars1D,"DA_cprocaB_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaB_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaG_R'
        allocate(cprocaG_R(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaG_R(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaG_R_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaG_R_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaG_R_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaG_R   = 0.d0
        scprocaG_R  = 0.d0
        cprocaG_R_p = 0.d0
        cprocaG_R_a = 0.d0
        cprocaG_R_bound = 0.d0
     else
        deallocate(cprocaG_R,scprocaG_R,cprocaG_R_p,cprocaG_R_a,cprocaG_R_bound)
     end if
  else if (contains(outvars0D,"cprocaG_R").or. &
           contains(outvars1D,"cprocaG_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaG_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaG_R(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaG_R = 0.d0
     else
        deallocate(D1_cprocaG_R)
     end if
  else if (contains(outvars0D,"D1_cprocaG_R").or. &
           contains(outvars1D,"D1_cprocaG_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaG_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaG_R(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaG_R = 0.d0
     else
        deallocate(DA_cprocaG_R)
     end if
  else if (contains(outvars0D,"DA_cprocaG_R").or. &
           contains(outvars1D,"DA_cprocaG_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaG_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaG_I'
        allocate(cprocaG_I(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaG_I(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaG_I_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaG_I_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaG_I_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaG_I   = 0.d0
        scprocaG_I  = 0.d0
        cprocaG_I_p = 0.d0
        cprocaG_I_a = 0.d0
        cprocaG_I_bound = 0.d0
     else
        deallocate(cprocaG_I,scprocaG_I,cprocaG_I_p,cprocaG_I_a,cprocaG_I_bound)
     end if
  else if (contains(outvars0D,"cprocaG_I").or. &
           contains(outvars1D,"cprocaG_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaG_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(D1_cprocaG_I(0:Nl-1,1-ghost:Nrmax))
        D1_cprocaG_I = 0.d0
     else
        deallocate(D1_cprocaG_I)
     end if
  else if (contains(outvars0D,"D1_cprocaG_I").or. &
           contains(outvars1D,"D1_cprocaG_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_cprocaG_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if ((contains(mattertype,"complexproca")).and.(shift/="none")) then
     if (trim(status)=='on') then
        allocate(DA_cprocaG_I(0:Nl-1,1-ghost:Nrmax))
        DA_cprocaG_I = 0.d0
     else
        deallocate(DA_cprocaG_I)
     end if
  else if (contains(outvars0D,"DA_cprocaG_I").or. &
           contains(outvars1D,"DA_cprocaG_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_cprocaG_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaL_R'
        allocate(cprocaL_R(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaL_R(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaL_R_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaL_R_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaL_R_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaL_R   = 0.d0
        scprocaL_R  = 0.d0
        cprocaL_R_p = 0.d0
        cprocaL_R_a = 0.d0
        cprocaL_R_bound = 0.d0
     else
        deallocate(cprocaL_R,scprocaL_R,cprocaL_R_p,cprocaL_R_a,cprocaL_R_bound)
     end if
  else if (contains(outvars0D,"cprocaL_R").or. &
           contains(outvars1D,"cprocaL_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaL_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cprocaL_I'
        allocate(cprocaL_I(0:Nl-1,1-ghost:Nrmax))
        allocate(scprocaL_I(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaL_I_p(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaL_I_a(0:Nl-1,1-ghost:Nrmax))
        allocate(cprocaL_I_bound(0:Nl-1,0:ghost-1,0:3))
        cprocaL_I   = 0.d0
        scprocaL_I  = 0.d0
        cprocaL_I_p = 0.d0
        cprocaL_I_a = 0.d0
        cprocaL_I_bound = 0.d0
     else
        deallocate(cprocaL_I,scprocaL_I,cprocaL_I_p,cprocaL_I_a,cprocaL_I_bound)
     end if
  else if (contains(outvars0D,"cprocaL_I").or. &
           contains(outvars1D,"cprocaL_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaL_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (contains(outvars0D,"Ccomplexproca_R").or. &
         contains(outvars1D,"Ccomplexproca_R")) then
        if (trim(status)=='on') then
           allocate(Ccomplexproca_R(0:Nl-1,1-ghost:Nrmax))
           Ccomplexproca_R = 0.d0
        else
           deallocate(Ccomplexproca_R)
        end if
     end if
  else if (contains(outvars0D,"Ccomplexproca_R").or. &
           contains(outvars1D,"Ccomplexproca_R")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array Ccomplexproca_R has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (contains(outvars0D,"Ccomplexproca_I").or. &
         contains(outvars1D,"Ccomplexproca_I")) then
        if (trim(status)=='on') then
           allocate(Ccomplexproca_I(0:Nl-1,1-ghost:Nrmax))
           Ccomplexproca_I = 0.d0
        else
           deallocate(Ccomplexproca_I)
        end if
     end if
  else if (contains(outvars0D,"Ccomplexproca_I").or. &
           contains(outvars1D,"Ccomplexproca_I")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array Ccomplexproca_I has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (contains(outvars0D,"cprocaPhi_norm").or. &
         contains(outvars1D,"cprocaPhi_norm")) then
        if (trim(status)=='on') then
           allocate(cprocaPhi_norm(0:Nl-1,1-ghost:Nrmax))
           cprocaPhi_norm = 0.d0
        else
           deallocate(cprocaPhi_norm)
        end if
     end if
  else if (contains(outvars0D,"cprocaPhi_norm").or. &
           contains(outvars1D,"cprocaPhi_norm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaPhi_norm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (contains(outvars0D,"cprocaA_norm").or. &
         contains(outvars1D,"cprocaA_norm")) then
        if (trim(status)=='on') then
           allocate(cprocaA_norm(0:Nl-1,1-ghost:Nrmax))
           cprocaA_norm = 0.d0
        else
           deallocate(cprocaA_norm)
        end if
     end if
  else if (contains(outvars0D,"cprocaA_norm").or. &
           contains(outvars1D,"cprocaA_norm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaA_norm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (contains(outvars0D,"cprocaB_norm").or. &
         contains(outvars1D,"cprocaB_norm")) then
        if (trim(status)=='on') then
           allocate(cprocaB_norm(0:Nl-1,1-ghost:Nrmax))
           cprocaB_norm = 0.d0
        else
           deallocate(cprocaB_norm)
        end if
     end if
  else if (contains(outvars0D,"cprocaB_norm").or. &
           contains(outvars1D,"cprocaB_norm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaB_norm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (contains(outvars0D,"cprocaE_norm").or. &
         contains(outvars1D,"cprocaE_norm")) then
        if (trim(status)=='on') then
           allocate(cprocaE_norm(0:Nl-1,1-ghost:Nrmax))
           cprocaE_norm = 0.d0
        else
           deallocate(cprocaE_norm)
        end if
     end if
  else if (contains(outvars0D,"cprocaE_norm").or. &
           contains(outvars1D,"cprocaE_norm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaE_norm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (contains(outvars0D,"cprocaXi_norm").or. &
         contains(outvars1D,"cprocaXi_norm")) then
        if (trim(status)=='on') then
           allocate(cprocaXi_norm(0:Nl-1,1-ghost:Nrmax))
           cprocaXi_norm = 0.d0
        else
           deallocate(cprocaXi_norm)
        end if
     end if
  else if (contains(outvars0D,"cprocaXi_norm").or. &
           contains(outvars1D,"cprocaXi_norm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cprocaXi_norm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (contains(outvars0D,"Ccomplexproca_norm").or. &
         contains(outvars1D,"Ccomplexproca_norm")) then
        if (trim(status)=='on') then
           allocate(Ccomplexproca_norm(0:Nl-1,1-ghost:Nrmax))
           Ccomplexproca_norm = 0.d0
        else
           deallocate(Ccomplexproca_norm)
        end if
     end if
  else if (contains(outvars0D,"Ccomplexproca_norm").or. &
           contains(outvars1D,"Ccomplexproca_norm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array Ccomplexproca_norm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(cproca_Qdens(0:Nl-1,1-ghost:Nrmax))
        cproca_Qdens = 0.d0
     else
        deallocate(cproca_Qdens)
     end if
  else if (contains(outvars0D,"cproca_Qdens").or. &
           contains(outvars1D,"cproca_Qdens")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cproca_Qdens has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (contains(outvars0D,"cproca_Qdens_r2").or. &
         contains(outvars1D,"cproca_Qdens_r2")) then
        if (trim(status)=='on') then
           allocate(cproca_Qdens_r2(0:Nl-1,1-ghost:Nrmax))
           cproca_Qdens_r2 = 0.d0
        else
           deallocate(cproca_Qdens_r2)
        end if
     end if
  else if (contains(outvars0D,"cproca_Qdens_r2").or. &
           contains(outvars1D,"cproca_Qdens_r2")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cproca_Qdens_r2 has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (trim(status)=='on') then
        allocate(cproca_Qflux(0:Nl-1,1-ghost:Nrmax))
        cproca_Qflux = 0.d0
     else
        deallocate(cproca_Qflux)
     end if
  else if (contains(outvars0D,"cproca_Qflux").or. &
           contains(outvars1D,"cproca_Qflux")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cproca_Qflux has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"complexproca")) then
     if (contains(outvars0D,"cproca_Qint").or. &
         contains(outvars1D,"cproca_Qint")) then
        if (trim(status)=='on') then
           allocate(cproca_Qint(0:Nl-1,1-ghost:Nrmax))
           cproca_Qint = 0.d0
        else
           deallocate(cproca_Qint)
        end if
     end if
  else if (contains(outvars0D,"cproca_Qint").or. &
           contains(outvars1D,"cproca_Qint")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cproca_Qint has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',dirac_FR'
        allocate(dirac_FR(0:Nl-1,1-ghost:Nrmax))
        allocate(sdirac_FR(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_FR_p(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_FR_a(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_FR_bound(0:Nl-1,0:ghost-1,0:3))
        dirac_FR   = 0.d0
        sdirac_FR  = 0.d0
        dirac_FR_p = 0.d0
        dirac_FR_a = 0.d0
        dirac_FR_bound = 0.d0
     else
        deallocate(dirac_FR,sdirac_FR,dirac_FR_p,dirac_FR_a,dirac_FR_bound)
     end if
  else if (contains(outvars0D,"dirac_FR").or. &
           contains(outvars1D,"dirac_FR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_FR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',dirac_FI'
        allocate(dirac_FI(0:Nl-1,1-ghost:Nrmax))
        allocate(sdirac_FI(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_FI_p(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_FI_a(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_FI_bound(0:Nl-1,0:ghost-1,0:3))
        dirac_FI   = 0.d0
        sdirac_FI  = 0.d0
        dirac_FI_p = 0.d0
        dirac_FI_a = 0.d0
        dirac_FI_bound = 0.d0
     else
        deallocate(dirac_FI,sdirac_FI,dirac_FI_p,dirac_FI_a,dirac_FI_bound)
     end if
  else if (contains(outvars0D,"dirac_FI").or. &
           contains(outvars1D,"dirac_FI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_FI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',dirac_GR'
        allocate(dirac_GR(0:Nl-1,1-ghost:Nrmax))
        allocate(sdirac_GR(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_GR_p(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_GR_a(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_GR_bound(0:Nl-1,0:ghost-1,0:3))
        dirac_GR   = 0.d0
        sdirac_GR  = 0.d0
        dirac_GR_p = 0.d0
        dirac_GR_a = 0.d0
        dirac_GR_bound = 0.d0
     else
        deallocate(dirac_GR,sdirac_GR,dirac_GR_p,dirac_GR_a,dirac_GR_bound)
     end if
  else if (contains(outvars0D,"dirac_GR").or. &
           contains(outvars1D,"dirac_GR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_GR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',dirac_GI'
        allocate(dirac_GI(0:Nl-1,1-ghost:Nrmax))
        allocate(sdirac_GI(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_GI_p(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_GI_a(0:Nl-1,1-ghost:Nrmax))
        allocate(dirac_GI_bound(0:Nl-1,0:ghost-1,0:3))
        dirac_GI   = 0.d0
        sdirac_GI  = 0.d0
        dirac_GI_p = 0.d0
        dirac_GI_a = 0.d0
        dirac_GI_bound = 0.d0
     else
        deallocate(dirac_GI,sdirac_GI,dirac_GI_p,dirac_GI_a,dirac_GI_bound)
     end if
  else if (contains(outvars0D,"dirac_GI").or. &
           contains(outvars1D,"dirac_GI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_GI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(D1_dirac_FR(0:Nl-1,1-ghost:Nrmax))
        D1_dirac_FR = 0.d0
     else
        deallocate(D1_dirac_FR)
     end if
  else if (contains(outvars0D,"D1_dirac_FR").or. &
           contains(outvars1D,"D1_dirac_FR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_dirac_FR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(D1_dirac_FI(0:Nl-1,1-ghost:Nrmax))
        D1_dirac_FI = 0.d0
     else
        deallocate(D1_dirac_FI)
     end if
  else if (contains(outvars0D,"D1_dirac_FI").or. &
           contains(outvars1D,"D1_dirac_FI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_dirac_FI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(DA_dirac_FR(0:Nl-1,1-ghost:Nrmax))
        DA_dirac_FR = 0.d0
     else
        deallocate(DA_dirac_FR)
     end if
  else if (contains(outvars0D,"DA_dirac_FR").or. &
           contains(outvars1D,"DA_dirac_FR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_dirac_FR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(DA_dirac_FI(0:Nl-1,1-ghost:Nrmax))
        DA_dirac_FI = 0.d0
     else
        deallocate(DA_dirac_FI)
     end if
  else if (contains(outvars0D,"DA_dirac_FI").or. &
           contains(outvars1D,"DA_dirac_FI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_dirac_FI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(D1_dirac_GR(0:Nl-1,1-ghost:Nrmax))
        D1_dirac_GR = 0.d0
     else
        deallocate(D1_dirac_GR)
     end if
  else if (contains(outvars0D,"D1_dirac_GR").or. &
           contains(outvars1D,"D1_dirac_GR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_dirac_GR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(D1_dirac_GI(0:Nl-1,1-ghost:Nrmax))
        D1_dirac_GI = 0.d0
     else
        deallocate(D1_dirac_GI)
     end if
  else if (contains(outvars0D,"D1_dirac_GI").or. &
           contains(outvars1D,"D1_dirac_GI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_dirac_GI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(DA_dirac_GR(0:Nl-1,1-ghost:Nrmax))
        DA_dirac_GR = 0.d0
     else
        deallocate(DA_dirac_GR)
     end if
  else if (contains(outvars0D,"DA_dirac_GR").or. &
           contains(outvars1D,"DA_dirac_GR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_dirac_GR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(DA_dirac_GI(0:Nl-1,1-ghost:Nrmax))
        DA_dirac_GI = 0.d0
     else
        deallocate(DA_dirac_GI)
     end if
  else if (contains(outvars0D,"DA_dirac_GI").or. &
           contains(outvars1D,"DA_dirac_GI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_dirac_GI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(dirac_HR(0:Nl-1,1-ghost:Nrmax))
        dirac_HR = 0.d0
     else
        deallocate(dirac_HR)
     end if
  else if (contains(outvars0D,"dirac_HR").or. &
           contains(outvars1D,"dirac_HR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_HR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(dirac_HI(0:Nl-1,1-ghost:Nrmax))
        dirac_HI = 0.d0
     else
        deallocate(dirac_HI)
     end if
  else if (contains(outvars0D,"dirac_HI").or. &
           contains(outvars1D,"dirac_HI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_HI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(D1_dirac_HR(0:Nl-1,1-ghost:Nrmax))
        D1_dirac_HR = 0.d0
     else
        deallocate(D1_dirac_HR)
     end if
  else if (contains(outvars0D,"D1_dirac_HR").or. &
           contains(outvars1D,"D1_dirac_HR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_dirac_HR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(D1_dirac_HI(0:Nl-1,1-ghost:Nrmax))
        D1_dirac_HI = 0.d0
     else
        deallocate(D1_dirac_HI)
     end if
  else if (contains(outvars0D,"D1_dirac_HI").or. &
           contains(outvars1D,"D1_dirac_HI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_dirac_HI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(dirac_PiFR(0:Nl-1,1-ghost:Nrmax))
        dirac_PiFR = 0.d0
     else
        deallocate(dirac_PiFR)
     end if
  else if (contains(outvars0D,"dirac_PiFR").or. &
           contains(outvars1D,"dirac_PiFR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_PiFR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(dirac_PiFI(0:Nl-1,1-ghost:Nrmax))
        dirac_PiFI = 0.d0
     else
        deallocate(dirac_PiFI)
     end if
  else if (contains(outvars0D,"dirac_PiFI").or. &
           contains(outvars1D,"dirac_PiFI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_PiFI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(dirac_PiGR(0:Nl-1,1-ghost:Nrmax))
        dirac_PiGR = 0.d0
     else
        deallocate(dirac_PiGR)
     end if
  else if (contains(outvars0D,"dirac_PiGR").or. &
           contains(outvars1D,"dirac_PiGR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_PiGR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(dirac_PiGI(0:Nl-1,1-ghost:Nrmax))
        dirac_PiGI = 0.d0
     else
        deallocate(dirac_PiGI)
     end if
  else if (contains(outvars0D,"dirac_PiGI").or. &
           contains(outvars1D,"dirac_PiGI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_PiGI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (contains(outvars0D,"dirac_F_norm").or. &
         contains(outvars1D,"dirac_F_norm")) then
        if (trim(status)=='on') then
           allocate(dirac_F_norm(0:Nl-1,1-ghost:Nrmax))
           dirac_F_norm = 0.d0
        else
           deallocate(dirac_F_norm)
        end if
     end if
  else if (contains(outvars0D,"dirac_F_norm").or. &
           contains(outvars1D,"dirac_F_norm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_F_norm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (contains(outvars0D,"dirac_G_norm").or. &
         contains(outvars1D,"dirac_G_norm")) then
        if (trim(status)=='on') then
           allocate(dirac_G_norm(0:Nl-1,1-ghost:Nrmax))
           dirac_G_norm = 0.d0
        else
           deallocate(dirac_G_norm)
        end if
     end if
  else if (contains(outvars0D,"dirac_G_norm").or. &
           contains(outvars1D,"dirac_G_norm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_G_norm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(dirac_dens(0:Nl-1,1-ghost:Nrmax))
        dirac_dens = 0.d0
     else
        deallocate(dirac_dens)
     end if
  else if (contains(outvars0D,"dirac_dens").or. &
           contains(outvars1D,"dirac_dens")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_dens has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (trim(status)=='on') then
        allocate(dirac_flux(0:Nl-1,1-ghost:Nrmax))
        dirac_flux = 0.d0
     else
        deallocate(dirac_flux)
     end if
  else if (contains(outvars0D,"dirac_flux").or. &
           contains(outvars1D,"dirac_flux")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_flux has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dirac")) then
     if (contains(outvars0D,"dirac_Nint").or. &
         contains(outvars1D,"dirac_Nint")) then
        if (trim(status)=='on') then
           allocate(dirac_Nint(0:Nl-1,1-ghost:Nrmax))
           dirac_Nint = 0.d0
        else
           deallocate(dirac_Nint)
        end if
     end if
  else if (contains(outvars0D,"dirac_Nint").or. &
           contains(outvars1D,"dirac_Nint")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dirac_Nint has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(outvars0D,"wpR_dirac").or. &
      contains(outvars1D,"wpR_dirac")) then
     if (trim(status)=='on') then
        allocate(wpR_dirac(0:Nl-1,1-ghost:Nrmax))
        wpR_dirac = 0.d0
     else
        deallocate(wpR_dirac)
     end if
  end if

  if (contains(outvars0D,"wmR_dirac").or. &
      contains(outvars1D,"wmR_dirac")) then
     if (trim(status)=='on') then
        allocate(wmR_dirac(0:Nl-1,1-ghost:Nrmax))
        wmR_dirac = 0.d0
     else
        deallocate(wmR_dirac)
     end if
  end if

  if (contains(outvars0D,"wpI_dirac").or. &
      contains(outvars1D,"wpI_dirac")) then
     if (trim(status)=='on') then
        allocate(wpI_dirac(0:Nl-1,1-ghost:Nrmax))
        wpI_dirac = 0.d0
     else
        deallocate(wpI_dirac)
     end if
  end if

  if (contains(outvars0D,"wmI_dirac").or. &
      contains(outvars1D,"wmI_dirac")) then
     if (trim(status)=='on') then
        allocate(wmI_dirac(0:Nl-1,1-ghost:Nrmax))
        wmI_dirac = 0.d0
     else
        deallocate(wmI_dirac)
     end if
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        allocate(dust_rho(0:Nl-1,1-ghost:Nrmax))
        dust_rho = 0.d0
     else
        deallocate(dust_rho)
     end if
  else if (contains(outvars0D,"dust_rho").or. &
           contains(outvars1D,"dust_rho")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dust_rho has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        allocate(dust_v(0:Nl-1,1-ghost:Nrmax))
        dust_v = 0.d0
     else
        deallocate(dust_v)
     end if
  else if (contains(outvars0D,"dust_v").or. &
           contains(outvars1D,"dust_v")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dust_v has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        allocate(dust_W(0:Nl-1,1-ghost:Nrmax))
        dust_W = 0.d0
     else
        deallocate(dust_W)
     end if
  else if (contains(outvars0D,"dust_W").or. &
           contains(outvars1D,"dust_W")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dust_W has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',dust_cD'
        allocate(dust_cD(0:Nl-1,1-ghost:Nrmax))
        allocate(sdust_cD(0:Nl-1,1-ghost:Nrmax))
        allocate(dust_cD_p(0:Nl-1,1-ghost:Nrmax))
        allocate(dust_cD_a(0:Nl-1,1-ghost:Nrmax))
        allocate(dust_cD_bound(0:Nl-1,0:ghost-1,0:3))
        dust_cD   = 0.d0
        sdust_cD  = 0.d0
        dust_cD_p = 0.d0
        dust_cD_a = 0.d0
        dust_cD_bound = 0.d0
     else
        deallocate(dust_cD,sdust_cD,dust_cD_p,dust_cD_a,dust_cD_bound)
     end if
  else if (contains(outvars0D,"dust_cD").or. &
           contains(outvars1D,"dust_cD")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dust_cD has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        allocate(D1_dust_cD(0:Nl-1,1-ghost:Nrmax))
        D1_dust_cD = 0.d0
     else
        deallocate(D1_dust_cD)
     end if
  else if (contains(outvars0D,"D1_dust_cD").or. &
           contains(outvars1D,"D1_dust_cD")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_dust_cD has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        allocate(DA_dust_cD(0:Nl-1,1-ghost:Nrmax))
        DA_dust_cD = 0.d0
     else
        deallocate(DA_dust_cD)
     end if
  else if (contains(outvars0D,"DA_dust_cD").or. &
           contains(outvars1D,"DA_dust_cD")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_dust_cD has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',dust_cE'
        allocate(dust_cE(0:Nl-1,1-ghost:Nrmax))
        allocate(sdust_cE(0:Nl-1,1-ghost:Nrmax))
        allocate(dust_cE_p(0:Nl-1,1-ghost:Nrmax))
        allocate(dust_cE_a(0:Nl-1,1-ghost:Nrmax))
        allocate(dust_cE_bound(0:Nl-1,0:ghost-1,0:3))
        dust_cE   = 0.d0
        sdust_cE  = 0.d0
        dust_cE_p = 0.d0
        dust_cE_a = 0.d0
        dust_cE_bound = 0.d0
     else
        deallocate(dust_cE,sdust_cE,dust_cE_p,dust_cE_a,dust_cE_bound)
     end if
  else if (contains(outvars0D,"dust_cE").or. &
           contains(outvars1D,"dust_cE")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dust_cE has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        allocate(D1_dust_cE(0:Nl-1,1-ghost:Nrmax))
        D1_dust_cE = 0.d0
     else
        deallocate(D1_dust_cE)
     end if
  else if (contains(outvars0D,"D1_dust_cE").or. &
           contains(outvars1D,"D1_dust_cE")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_dust_cE has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        allocate(DA_dust_cE(0:Nl-1,1-ghost:Nrmax))
        DA_dust_cE = 0.d0
     else
        deallocate(DA_dust_cE)
     end if
  else if (contains(outvars0D,"DA_dust_cE").or. &
           contains(outvars1D,"DA_dust_cE")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_dust_cE has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',dust_cS'
        allocate(dust_cS(0:Nl-1,1-ghost:Nrmax))
        allocate(sdust_cS(0:Nl-1,1-ghost:Nrmax))
        allocate(dust_cS_p(0:Nl-1,1-ghost:Nrmax))
        allocate(dust_cS_a(0:Nl-1,1-ghost:Nrmax))
        allocate(dust_cS_bound(0:Nl-1,0:ghost-1,0:3))
        dust_cS   = 0.d0
        sdust_cS  = 0.d0
        dust_cS_p = 0.d0
        dust_cS_a = 0.d0
        dust_cS_bound = 0.d0
     else
        deallocate(dust_cS,sdust_cS,dust_cS_p,dust_cS_a,dust_cS_bound)
     end if
  else if (contains(outvars0D,"dust_cS").or. &
           contains(outvars1D,"dust_cS")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dust_cS has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        allocate(D1_dust_cS(0:Nl-1,1-ghost:Nrmax))
        D1_dust_cS = 0.d0
     else
        deallocate(D1_dust_cS)
     end if
  else if (contains(outvars0D,"D1_dust_cS").or. &
           contains(outvars1D,"D1_dust_cS")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_dust_cS has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        allocate(DA_dust_cS(0:Nl-1,1-ghost:Nrmax))
        DA_dust_cS = 0.d0
     else
        deallocate(DA_dust_cS)
     end if
  else if (contains(outvars0D,"DA_dust_cS").or. &
           contains(outvars1D,"DA_dust_cS")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_dust_cS has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"dust")) then
     if (trim(status)=='on') then
        allocate(dust_restmass(0:Nl-1,1-ghost:Nrmax))
        dust_restmass = 0.d0
     else
        deallocate(dust_restmass)
     end if
  else if (contains(outvars0D,"dust_restmass").or. &
           contains(outvars1D,"dust_restmass")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array dust_restmass has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',fluid_rho'
        allocate(fluid_rho(0:Nl-1,1-ghost:Nrmax))
        fluid_rho = 0.d0
     else
        deallocate(fluid_rho)
     end if
  else if (contains(outvars0D,"fluid_rho").or. &
           contains(outvars1D,"fluid_rho")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_rho has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(fluid_rhotot(0:Nl-1,1-ghost:Nrmax))
        fluid_rhotot = 0.d0
     else
        deallocate(fluid_rhotot)
     end if
  else if (contains(outvars0D,"fluid_rhotot").or. &
           contains(outvars1D,"fluid_rhotot")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_rhotot has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',fluid_e'
        allocate(fluid_e(0:Nl-1,1-ghost:Nrmax))
        fluid_e = 0.d0
     else
        deallocate(fluid_e)
     end if
  else if (contains(outvars0D,"fluid_e").or. &
           contains(outvars1D,"fluid_e")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_e has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',fluid_p'
        allocate(fluid_p(0:Nl-1,1-ghost:Nrmax))
        fluid_p = 0.d0
     else
        deallocate(fluid_p)
     end if
  else if (contains(outvars0D,"fluid_p").or. &
           contains(outvars1D,"fluid_p")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_p has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(D1_fluid_p(0:Nl-1,1-ghost:Nrmax))
        D1_fluid_p = 0.d0
     else
        deallocate(D1_fluid_p)
     end if
  else if (contains(outvars0D,"D1_fluid_p").or. &
           contains(outvars1D,"D1_fluid_p")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_fluid_p has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',fluid_h'
        allocate(fluid_h(0:Nl-1,1-ghost:Nrmax))
        fluid_h = 0.d0
     else
        deallocate(fluid_h)
     end if
  else if (contains(outvars0D,"fluid_h").or. &
           contains(outvars1D,"fluid_h")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_h has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',fluid_v'
        allocate(fluid_v(0:Nl-1,1-ghost:Nrmax))
        fluid_v = 0.d0
     else
        deallocate(fluid_v)
     end if
  else if (contains(outvars0D,"fluid_v").or. &
           contains(outvars1D,"fluid_v")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_v has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',fluid_W'
        allocate(fluid_W(0:Nl-1,1-ghost:Nrmax))
        fluid_W = 0.d0
     else
        deallocate(fluid_W)
     end if
  else if (contains(outvars0D,"fluid_W").or. &
           contains(outvars1D,"fluid_W")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_W has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(fluid_vs(0:Nl-1,1-ghost:Nrmax))
        fluid_vs = 0.d0
     else
        deallocate(fluid_vs)
     end if
  else if (contains(outvars0D,"fluid_vs").or. &
           contains(outvars1D,"fluid_vs")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_vs has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(fluid_Mach(0:Nl-1,1-ghost:Nrmax))
        fluid_Mach = 0.d0
     else
        deallocate(fluid_Mach)
     end if
  else if (contains(outvars0D,"fluid_Mach").or. &
           contains(outvars1D,"fluid_Mach")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_Mach has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(fluid_vcp(0:Nl-1,1-ghost:Nrmax))
        fluid_vcp = 0.d0
     else
        deallocate(fluid_vcp)
     end if
  else if (contains(outvars0D,"fluid_vcp").or. &
           contains(outvars1D,"fluid_vcp")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_vcp has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(fluid_vcm(0:Nl-1,1-ghost:Nrmax))
        fluid_vcm = 0.d0
     else
        deallocate(fluid_vcm)
     end if
  else if (contains(outvars0D,"fluid_vcm").or. &
           contains(outvars1D,"fluid_vcm")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_vcm has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (contains(outvars0D,"fluid_wb_rest").or. &
         contains(outvars1D,"fluid_wb_rest")) then
        if (trim(status)=='on') then
           allocate(fluid_wb_rest(0:Nl-1,1-ghost:Nrmax))
           fluid_wb_rest = 0.d0
        else
           deallocate(fluid_wb_rest)
        end if
     end if
  else if (contains(outvars0D,"fluid_wb_rest").or. &
           contains(outvars1D,"fluid_wb_rest")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_wb_rest has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (contains(outvars0D,"fluid_wb_tot").or. &
         contains(outvars1D,"fluid_wb_tot")) then
        if (trim(status)=='on') then
           allocate(fluid_wb_tot(0:Nl-1,1-ghost:Nrmax))
           fluid_wb_tot = 0.d0
        else
           deallocate(fluid_wb_tot)
        end if
     end if
  else if (contains(outvars0D,"fluid_wb_tot").or. &
           contains(outvars1D,"fluid_wb_tot")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_wb_tot has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',fluid_cD'
        allocate(fluid_cD(0:Nl-1,1-ghost:Nrmax))
        allocate(sfluid_cD(0:Nl-1,1-ghost:Nrmax))
        allocate(fluid_cD_p(0:Nl-1,1-ghost:Nrmax))
        allocate(fluid_cD_a(0:Nl-1,1-ghost:Nrmax))
        allocate(fluid_cD_bound(0:Nl-1,0:ghost-1,0:3))
        fluid_cD   = 0.d0
        sfluid_cD  = 0.d0
        fluid_cD_p = 0.d0
        fluid_cD_a = 0.d0
        fluid_cD_bound = 0.d0
     else
        deallocate(fluid_cD,sfluid_cD,fluid_cD_p,fluid_cD_a,fluid_cD_bound)
     end if
  else if (contains(outvars0D,"fluid_cD").or. &
           contains(outvars1D,"fluid_cD")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_cD has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(D1_fluid_cD(0:Nl-1,1-ghost:Nrmax))
        D1_fluid_cD = 0.d0
     else
        deallocate(D1_fluid_cD)
     end if
  else if (contains(outvars0D,"D1_fluid_cD").or. &
           contains(outvars1D,"D1_fluid_cD")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_fluid_cD has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(DA_fluid_cD(0:Nl-1,1-ghost:Nrmax))
        DA_fluid_cD = 0.d0
     else
        deallocate(DA_fluid_cD)
     end if
  else if (contains(outvars0D,"DA_fluid_cD").or. &
           contains(outvars1D,"DA_fluid_cD")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_fluid_cD has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',fluid_cE'
        allocate(fluid_cE(0:Nl-1,1-ghost:Nrmax))
        allocate(sfluid_cE(0:Nl-1,1-ghost:Nrmax))
        allocate(fluid_cE_p(0:Nl-1,1-ghost:Nrmax))
        allocate(fluid_cE_a(0:Nl-1,1-ghost:Nrmax))
        allocate(fluid_cE_bound(0:Nl-1,0:ghost-1,0:3))
        fluid_cE   = 0.d0
        sfluid_cE  = 0.d0
        fluid_cE_p = 0.d0
        fluid_cE_a = 0.d0
        fluid_cE_bound = 0.d0
     else
        deallocate(fluid_cE,sfluid_cE,fluid_cE_p,fluid_cE_a,fluid_cE_bound)
     end if
  else if (contains(outvars0D,"fluid_cE").or. &
           contains(outvars1D,"fluid_cE")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_cE has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(D1_fluid_cE(0:Nl-1,1-ghost:Nrmax))
        D1_fluid_cE = 0.d0
     else
        deallocate(D1_fluid_cE)
     end if
  else if (contains(outvars0D,"D1_fluid_cE").or. &
           contains(outvars1D,"D1_fluid_cE")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_fluid_cE has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(DA_fluid_cE(0:Nl-1,1-ghost:Nrmax))
        DA_fluid_cE = 0.d0
     else
        deallocate(DA_fluid_cE)
     end if
  else if (contains(outvars0D,"DA_fluid_cE").or. &
           contains(outvars1D,"DA_fluid_cE")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_fluid_cE has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',fluid_cS'
        allocate(fluid_cS(0:Nl-1,1-ghost:Nrmax))
        allocate(sfluid_cS(0:Nl-1,1-ghost:Nrmax))
        allocate(fluid_cS_p(0:Nl-1,1-ghost:Nrmax))
        allocate(fluid_cS_a(0:Nl-1,1-ghost:Nrmax))
        allocate(fluid_cS_bound(0:Nl-1,0:ghost-1,0:3))
        fluid_cS   = 0.d0
        sfluid_cS  = 0.d0
        fluid_cS_p = 0.d0
        fluid_cS_a = 0.d0
        fluid_cS_bound = 0.d0
     else
        deallocate(fluid_cS,sfluid_cS,fluid_cS_p,fluid_cS_a,fluid_cS_bound)
     end if
  else if (contains(outvars0D,"fluid_cS").or. &
           contains(outvars1D,"fluid_cS")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_cS has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(D1_fluid_cS(0:Nl-1,1-ghost:Nrmax))
        D1_fluid_cS = 0.d0
     else
        deallocate(D1_fluid_cS)
     end if
  else if (contains(outvars0D,"D1_fluid_cS").or. &
           contains(outvars1D,"D1_fluid_cS")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array D1_fluid_cS has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(DA_fluid_cS(0:Nl-1,1-ghost:Nrmax))
        DA_fluid_cS = 0.d0
     else
        deallocate(DA_fluid_cS)
     end if
  else if (contains(outvars0D,"DA_fluid_cS").or. &
           contains(outvars1D,"DA_fluid_cS")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array DA_fluid_cS has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(fluid_q(0:Nl-1,1-ghost:Nrmax))
        fluid_q = 0.d0
     else
        deallocate(fluid_q)
     end if
  else if (contains(outvars0D,"fluid_q").or. &
           contains(outvars1D,"fluid_q")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_q has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(fluid_restmass(0:Nl-1,1-ghost:Nrmax))
        fluid_restmass = 0.d0
     else
        deallocate(fluid_restmass)
     end if
  else if (contains(outvars0D,"fluid_restmass").or. &
           contains(outvars1D,"fluid_restmass")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_restmass has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (contains(mattertype,"fluid")) then
     if (trim(status)=='on') then
        allocate(fluid_totalmass(0:Nl-1,1-ghost:Nrmax))
        fluid_totalmass = 0.d0
     else
        deallocate(fluid_totalmass)
     end if
  else if (contains(outvars0D,"fluid_totalmass").or. &
           contains(outvars1D,"fluid_totalmass")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array fluid_totalmass has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_tau'
        allocate(cosmobg_tau(0:Nl-1))
        allocate(scosmobg_tau(0:Nl-1))
        allocate(cosmobg_tau_p(0:Nl-1))
        allocate(cosmobg_tau_a(0:Nl-1))
        cosmobg_tau   = 0.d0
        scosmobg_tau  = 0.d0
        cosmobg_tau_p = 0.d0
        cosmobg_tau_a = 0.d0
     else
        deallocate(cosmobg_tau,scosmobg_tau,cosmobg_tau_p,cosmobg_tau_a)
     end if
  else if (contains(outvars0D,"cosmobg_tau").or. &
           contains(outvars1D,"cosmobg_tau")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_tau has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_H'
        allocate(cosmobg_H(0:Nl-1))
        allocate(scosmobg_H(0:Nl-1))
        allocate(cosmobg_H_p(0:Nl-1))
        allocate(cosmobg_H_a(0:Nl-1))
        cosmobg_H   = 0.d0
        scosmobg_H  = 0.d0
        cosmobg_H_p = 0.d0
        cosmobg_H_a = 0.d0
     else
        deallocate(cosmobg_H,scosmobg_H,cosmobg_H_p,cosmobg_H_a)
     end if
  else if (contains(outvars0D,"cosmobg_H").or. &
           contains(outvars1D,"cosmobg_H")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_H has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_a'
        allocate(cosmobg_a(0:Nl-1))
        allocate(scosmobg_a(0:Nl-1))
        allocate(cosmobg_a_p(0:Nl-1))
        allocate(cosmobg_a_a(0:Nl-1))
        cosmobg_a   = 0.d0
        scosmobg_a  = 0.d0
        cosmobg_a_p = 0.d0
        cosmobg_a_a = 0.d0
     else
        deallocate(cosmobg_a,scosmobg_a,cosmobg_a_p,cosmobg_a_a)
     end if
  else if (contains(outvars0D,"cosmobg_a").or. &
           contains(outvars1D,"cosmobg_a")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_a has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_phi'
        allocate(cosmobg_phi(0:Nl-1))
        allocate(scosmobg_phi(0:Nl-1))
        allocate(cosmobg_phi_p(0:Nl-1))
        allocate(cosmobg_phi_a(0:Nl-1))
        cosmobg_phi   = 0.d0
        scosmobg_phi  = 0.d0
        cosmobg_phi_p = 0.d0
        cosmobg_phi_a = 0.d0
     else
        deallocate(cosmobg_phi,scosmobg_phi,cosmobg_phi_p,cosmobg_phi_a)
     end if
  else if (contains(outvars0D,"cosmobg_phi").or. &
           contains(outvars1D,"cosmobg_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_trK'
        allocate(cosmobg_trK(0:Nl-1))
        allocate(scosmobg_trK(0:Nl-1))
        allocate(cosmobg_trK_p(0:Nl-1))
        allocate(cosmobg_trK_a(0:Nl-1))
        cosmobg_trK   = 0.d0
        scosmobg_trK  = 0.d0
        cosmobg_trK_p = 0.d0
        cosmobg_trK_a = 0.d0
     else
        deallocate(cosmobg_trK,scosmobg_trK,cosmobg_trK_p,cosmobg_trK_a)
     end if
  else if (contains(outvars0D,"cosmobg_trK").or. &
           contains(outvars1D,"cosmobg_trK")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_trK has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_alpha'
        allocate(cosmobg_alpha(0:Nl-1))
        allocate(scosmobg_alpha(0:Nl-1))
        allocate(cosmobg_alpha_p(0:Nl-1))
        allocate(cosmobg_alpha_a(0:Nl-1))
        cosmobg_alpha   = 0.d0
        scosmobg_alpha  = 0.d0
        cosmobg_alpha_p = 0.d0
        cosmobg_alpha_a = 0.d0
     else
        deallocate(cosmobg_alpha,scosmobg_alpha,cosmobg_alpha_p,cosmobg_alpha_a)
     end if
  else if (contains(outvars0D,"cosmobg_alpha").or. &
           contains(outvars1D,"cosmobg_alpha")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_alpha has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(cosmobg_falpha(0:Nl-1))
        cosmobg_falpha = 0.d0
     else
        deallocate(cosmobg_falpha)
     end if
  else if (contains(outvars0D,"cosmobg_falpha").or. &
           contains(outvars1D,"cosmobg_falpha")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_falpha has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(cosmobg_vlight(0:Nl-1))
        cosmobg_vlight = 0.d0
     else
        deallocate(cosmobg_vlight)
     end if
  else if (contains(outvars0D,"cosmobg_vlight").or. &
           contains(outvars1D,"cosmobg_vlight")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_vlight has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(cosmobg_rho(0:Nl-1))
        cosmobg_rho = 0.d0
     else
        deallocate(cosmobg_rho)
     end if
  else if (contains(outvars0D,"cosmobg_rho").or. &
           contains(outvars1D,"cosmobg_rho")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_rho has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(cosmobg_p(0:Nl-1))
        cosmobg_p = 0.d0
     else
        deallocate(cosmobg_p)
     end if
  else if (contains(outvars0D,"cosmobg_p").or. &
           contains(outvars1D,"cosmobg_p")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_p has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(alpha_pert(0:Nl-1,1-ghost:Nrmax))
        alpha_pert = 0.d0
     else
        deallocate(alpha_pert)
     end if
  else if (contains(outvars0D,"alpha_pert").or. &
           contains(outvars1D,"alpha_pert")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array alpha_pert has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(phi_pert(0:Nl-1,1-ghost:Nrmax))
        phi_pert = 0.d0
     else
        deallocate(phi_pert)
     end if
  else if (contains(outvars0D,"phi_pert").or. &
           contains(outvars1D,"phi_pert")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array phi_pert has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(trK_pert(0:Nl-1,1-ghost:Nrmax))
        trK_pert = 0.d0
     else
        deallocate(trK_pert)
     end if
  else if (contains(outvars0D,"trK_pert").or. &
           contains(outvars1D,"trK_pert")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array trK_pert has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(rho_pert(0:Nl-1,1-ghost:Nrmax))
        rho_pert = 0.d0
     else
        deallocate(rho_pert)
     end if
  else if (contains(outvars0D,"rho_pert").or. &
           contains(outvars1D,"rho_pert")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array rho_pert has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(rho_contrast(0:Nl-1,1-ghost:Nrmax))
        rho_contrast = 0.d0
     else
        deallocate(rho_contrast)
     end if
  else if (contains(outvars0D,"rho_contrast").or. &
           contains(outvars1D,"rho_contrast")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array rho_contrast has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(rho_pert_int(0:Nl-1,1-ghost:Nrmax))
        rho_pert_int = 0.d0
     else
        deallocate(rho_pert_int)
     end if
  else if (contains(outvars0D,"rho_pert_int").or. &
           contains(outvars1D,"rho_pert_int")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array rho_pert_int has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(rho_contrast_int(0:Nl-1,1-ghost:Nrmax))
        rho_contrast_int = 0.d0
     else
        deallocate(rho_contrast_int)
     end if
  else if (contains(outvars0D,"rho_contrast_int").or. &
           contains(outvars1D,"rho_contrast_int")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array rho_contrast_int has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(cosmo_compactness(0:Nl-1,1-ghost:Nrmax))
        cosmo_compactness = 0.d0
     else
        deallocate(cosmo_compactness)
     end if
  else if (contains(outvars0D,"cosmo_compactness").or. &
           contains(outvars1D,"cosmo_compactness")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmo_compactness has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_scalar_phi'
        allocate(cosmobg_scalar_phi(0:Nl-1))
        allocate(scosmobg_scalar_phi(0:Nl-1))
        allocate(cosmobg_scalar_phi_p(0:Nl-1))
        allocate(cosmobg_scalar_phi_a(0:Nl-1))
        cosmobg_scalar_phi   = 0.d0
        scosmobg_scalar_phi  = 0.d0
        cosmobg_scalar_phi_p = 0.d0
        cosmobg_scalar_phi_a = 0.d0
     else
        deallocate(cosmobg_scalar_phi,scosmobg_scalar_phi,cosmobg_scalar_phi_p,cosmobg_scalar_phi_a)
     end if
  else if (contains(outvars0D,"cosmobg_scalar_phi").or. &
           contains(outvars1D,"cosmobg_scalar_phi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_scalar_phi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_scalar_pi'
        allocate(cosmobg_scalar_pi(0:Nl-1))
        allocate(scosmobg_scalar_pi(0:Nl-1))
        allocate(cosmobg_scalar_pi_p(0:Nl-1))
        allocate(cosmobg_scalar_pi_a(0:Nl-1))
        cosmobg_scalar_pi   = 0.d0
        scosmobg_scalar_pi  = 0.d0
        cosmobg_scalar_pi_p = 0.d0
        cosmobg_scalar_pi_a = 0.d0
     else
        deallocate(cosmobg_scalar_pi,scosmobg_scalar_pi,cosmobg_scalar_pi_p,cosmobg_scalar_pi_a)
     end if
  else if (contains(outvars0D,"cosmobg_scalar_pi").or. &
           contains(outvars1D,"cosmobg_scalar_pi")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_scalar_pi has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(cosmobg_scalar_V(0:Nl-1))
        cosmobg_scalar_V = 0.d0
     else
        deallocate(cosmobg_scalar_V)
     end if
  else if (contains(outvars0D,"cosmobg_scalar_V").or. &
           contains(outvars1D,"cosmobg_scalar_V")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_scalar_V has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(cosmobg_scalar_VP(0:Nl-1))
        cosmobg_scalar_VP = 0.d0
     else
        deallocate(cosmobg_scalar_VP)
     end if
  else if (contains(outvars0D,"cosmobg_scalar_VP").or. &
           contains(outvars1D,"cosmobg_scalar_VP")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_scalar_VP has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_complex_phiR'
        allocate(cosmobg_complex_phiR(0:Nl-1))
        allocate(scosmobg_complex_phiR(0:Nl-1))
        allocate(cosmobg_complex_phiR_p(0:Nl-1))
        allocate(cosmobg_complex_phiR_a(0:Nl-1))
        cosmobg_complex_phiR   = 0.d0
        scosmobg_complex_phiR  = 0.d0
        cosmobg_complex_phiR_p = 0.d0
        cosmobg_complex_phiR_a = 0.d0
     else
        deallocate(cosmobg_complex_phiR,scosmobg_complex_phiR,cosmobg_complex_phiR_p,cosmobg_complex_phiR_a)
     end if
  else if (contains(outvars0D,"cosmobg_complex_phiR").or. &
           contains(outvars1D,"cosmobg_complex_phiR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_complex_phiR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_complex_phiI'
        allocate(cosmobg_complex_phiI(0:Nl-1))
        allocate(scosmobg_complex_phiI(0:Nl-1))
        allocate(cosmobg_complex_phiI_p(0:Nl-1))
        allocate(cosmobg_complex_phiI_a(0:Nl-1))
        cosmobg_complex_phiI   = 0.d0
        scosmobg_complex_phiI  = 0.d0
        cosmobg_complex_phiI_p = 0.d0
        cosmobg_complex_phiI_a = 0.d0
     else
        deallocate(cosmobg_complex_phiI,scosmobg_complex_phiI,cosmobg_complex_phiI_p,cosmobg_complex_phiI_a)
     end if
  else if (contains(outvars0D,"cosmobg_complex_phiI").or. &
           contains(outvars1D,"cosmobg_complex_phiI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_complex_phiI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_complex_piR'
        allocate(cosmobg_complex_piR(0:Nl-1))
        allocate(scosmobg_complex_piR(0:Nl-1))
        allocate(cosmobg_complex_piR_p(0:Nl-1))
        allocate(cosmobg_complex_piR_a(0:Nl-1))
        cosmobg_complex_piR   = 0.d0
        scosmobg_complex_piR  = 0.d0
        cosmobg_complex_piR_p = 0.d0
        cosmobg_complex_piR_a = 0.d0
     else
        deallocate(cosmobg_complex_piR,scosmobg_complex_piR,cosmobg_complex_piR_p,cosmobg_complex_piR_a)
     end if
  else if (contains(outvars0D,"cosmobg_complex_piR").or. &
           contains(outvars1D,"cosmobg_complex_piR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_complex_piR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        checkvars = trim(checkvars) // ',cosmobg_complex_piI'
        allocate(cosmobg_complex_piI(0:Nl-1))
        allocate(scosmobg_complex_piI(0:Nl-1))
        allocate(cosmobg_complex_piI_p(0:Nl-1))
        allocate(cosmobg_complex_piI_a(0:Nl-1))
        cosmobg_complex_piI   = 0.d0
        scosmobg_complex_piI  = 0.d0
        cosmobg_complex_piI_p = 0.d0
        cosmobg_complex_piI_a = 0.d0
     else
        deallocate(cosmobg_complex_piI,scosmobg_complex_piI,cosmobg_complex_piI_p,cosmobg_complex_piI_a)
     end if
  else if (contains(outvars0D,"cosmobg_complex_piI").or. &
           contains(outvars1D,"cosmobg_complex_piI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_complex_piI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(cosmobg_complex_V(0:Nl-1))
        cosmobg_complex_V = 0.d0
     else
        deallocate(cosmobg_complex_V)
     end if
  else if (contains(outvars0D,"cosmobg_complex_V").or. &
           contains(outvars1D,"cosmobg_complex_V")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_complex_V has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(cosmobg_complex_VPR(0:Nl-1))
        cosmobg_complex_VPR = 0.d0
     else
        deallocate(cosmobg_complex_VPR)
     end if
  else if (contains(outvars0D,"cosmobg_complex_VPR").or. &
           contains(outvars1D,"cosmobg_complex_VPR")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_complex_VPR has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (cosmic_run) then
     if (trim(status)=='on') then
        allocate(cosmobg_complex_VPI(0:Nl-1))
        cosmobg_complex_VPI = 0.d0
     else
        deallocate(cosmobg_complex_VPI)
     end if
  else if (contains(outvars0D,"cosmobg_complex_VPI").or. &
           contains(outvars1D,"cosmobg_complex_VPI")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array cosmobg_complex_VPI has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (newr) then
     if (trim(status)=='on') then
        allocate(r_trans(0:Nl-1,1-ghost:Nrmax))
        r_trans = 0.d0
     else
        deallocate(r_trans)
     end if
  else if (contains(outvars0D,"r_trans").or. &
           contains(outvars1D,"r_trans")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array r_trans has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trackminkowski) then
     if (trim(status)=='on') then
        allocate(R_MINK(0:Nl-1,1-ghost:Nrmax))
        R_MINK = 0.d0
     else
        deallocate(R_MINK)
     end if
  else if (contains(outvars0D,"R_MINK").or. &
           contains(outvars1D,"R_MINK")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array R_MINK has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trackminkowski) then
     if (trim(status)=='on') then
        allocate(R_MINK_P(0:Nl-1,1-ghost:Nrmax))
        R_MINK_P = 0.d0
     else
        deallocate(R_MINK_P)
     end if
  else if (contains(outvars0D,"R_MINK_P").or. &
           contains(outvars1D,"R_MINK_P")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array R_MINK_P has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trackminkowski) then
     if (trim(status)=='on') then
        allocate(R_MINK_S(0:Nl-1,1-ghost:Nrmax))
        R_MINK_S = 0.d0
     else
        deallocate(R_MINK_S)
     end if
  else if (contains(outvars0D,"R_MINK_S").or. &
           contains(outvars1D,"R_MINK_S")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array R_MINK_S has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trackminkowski) then
     if (trim(status)=='on') then
        allocate(T_MINK(0:Nl-1,1-ghost:Nrmax))
        T_MINK = 0.d0
     else
        deallocate(T_MINK)
     end if
  else if (contains(outvars0D,"T_MINK").or. &
           contains(outvars1D,"T_MINK")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array T_MINK has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trackminkowski) then
     if (trim(status)=='on') then
        allocate(T_MINK_P(0:Nl-1,1-ghost:Nrmax))
        T_MINK_P = 0.d0
     else
        deallocate(T_MINK_P)
     end if
  else if (contains(outvars0D,"T_MINK_P").or. &
           contains(outvars1D,"T_MINK_P")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array T_MINK_P has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trackminkowski) then
     if (trim(status)=='on') then
        allocate(T_MINK_S(0:Nl-1,1-ghost:Nrmax))
        T_MINK_S = 0.d0
     else
        deallocate(T_MINK_S)
     end if
  else if (contains(outvars0D,"T_MINK_S").or. &
           contains(outvars1D,"T_MINK_S")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array T_MINK_S has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trackschwarzschild) then
     if (trim(status)=='on') then
        allocate(CSI_SCHWARZ(0:Nl-1,1-ghost:Nrmax))
        CSI_SCHWARZ = 0.d0
     else
        deallocate(CSI_SCHWARZ)
     end if
  else if (contains(outvars0D,"CSI_SCHWARZ").or. &
           contains(outvars1D,"CSI_SCHWARZ")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array CSI_SCHWARZ has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trackschwarzschild) then
     if (trim(status)=='on') then
        allocate(CSI_SCHWARZ_P(0:Nl-1,1-ghost:Nrmax))
        CSI_SCHWARZ_P = 0.d0
     else
        deallocate(CSI_SCHWARZ_P)
     end if
  else if (contains(outvars0D,"CSI_SCHWARZ_P").or. &
           contains(outvars1D,"CSI_SCHWARZ_P")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array CSI_SCHWARZ_P has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trackschwarzschild) then
     if (trim(status)=='on') then
        allocate(ETA_SCHWARZ(0:Nl-1,1-ghost:Nrmax))
        ETA_SCHWARZ = 0.d0
     else
        deallocate(ETA_SCHWARZ)
     end if
  else if (contains(outvars0D,"ETA_SCHWARZ").or. &
           contains(outvars1D,"ETA_SCHWARZ")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array ETA_SCHWARZ has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  if (trackschwarzschild) then
     if (trim(status)=='on') then
        allocate(ETA_SCHWARZ_P(0:Nl-1,1-ghost:Nrmax))
        ETA_SCHWARZ_P = 0.d0
     else
        deallocate(ETA_SCHWARZ_P)
     end if
  else if (contains(outvars0D,"ETA_SCHWARZ_P").or. &
           contains(outvars1D,"ETA_SCHWARZ_P")) then
     if (rank==0) then
        print *
        print *, 'Error in parfile: array ETA_SCHWARZ_P has no storage,'
        print *, 'so no output for it is possible.'
        print *, 'Aborting! (subroutine allocatearrays.f90)'
        print *
     end if
     call die
  end if

  end subroutine allocatearrays

