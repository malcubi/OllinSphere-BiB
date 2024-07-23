!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/analysis_matter.f90,v 1.53 2024/07/23 20:11:51 malcubi Exp $

  subroutine analysis_matter

! ********************************************************
! ***   CALCULATION OF ANALYSIS VARIABLES FOR OUTPUT   ***
! ********************************************************

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical contains

  integer i,l

  real(8) half,one,two,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  one = 1.d0
  two = 2.d0

  smallpi = acos(-one)


! **************************************************
! ***   INTEGRATED MASS AND DENSITY TIMES r**2   ***
! **************************************************

! Density times r**2.

  if (allocated(rho_r2)) then
     rho_r2 = rho*r**2
  end if

! The integrated mass is defined as:
!
!      /                          2        2           2       2
! m  = | [ 4 pi rho  + (1/4) ( KTA  + 2 KTB   - 2/3 trK ) ] r_a  dr_a
!      /
!
! with r_a the areal radius.
!
! Notice that we don't call the mass integral routine
! for cosmological runs, as in that case the expression
! above reduces to zero.

  if (.not.cosmic_run) then
     call massintegral
  end if


! ******************
! ***   VIRIAL   ***
! ******************

! Virial integrals. Look at the routine virial.f90
! to see how they are defined.

  if (allocated(virial1)) then
     call virialintegral1
  end if

  if (allocated(virial2)) then
     call virialintegral2
  end if


! ***********************
! ***   KODAMA MASS   ***
! ***********************

! Kodama mass.

  if (allocated(Kodama_mass)) then
     call kodamamass
  end if


! ************************
! ***   SCALAR FIELD   ***
! ************************

  if (contains(mattertype,"scalar")) then

!    Scalar eigenfields.

     if (allocated(wp_scalar)) then
        wp_scalar = scalar_pi - scalar_xi/(sqrt(A)*psi2)
     end if

     if (allocated(wm_scalar)) then
        wm_scalar = scalar_pi + scalar_xi/(sqrt(A)*psi2)
     end if

  end if


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

  if (contains(mattertype,"complex")) then

!    Norm of complex scalar field: sqrt(phiR**2 + phiI**2).

     if (allocated(complex_phi_norm)) then
        complex_phi_norm = sqrt(complex_phiR**2 + complex_phiI**2)
     end if

!    Boson density times r**2.

     if (allocated(complex_Bdens_r2)) then
        complex_Bdens_r2 = complex_Bdens*r**2
     end if

!    The total boson number (integrated boson charge) is then:
!
!                        /r                2
!    complex_NB  =  4 pi |  complex_Bdens R  dR
!                        /0
!
!    with R the Schwarzschild radius.

     if (allocated(complex_NB)) then

        call bosonintegral

!       Output total binding energy, but only at t=0 and
!       for certain type of initial data.

        if ((t(0)==0.d0).and.(rank==0)) then
           if (idata=="bosonstar") then
              write(*,'(A,E19.12)') ' Binding energy (M-m*Q) = ',mass_int(0,Nr) - complex_mass*complex_NB(0,Nr)
           end if
        end if

     end if

  end if


! **************************************
! ***   COMPLEX GHOST SCALAR FIELD   ***
! **************************************

  if (contains(mattertype,"complexghost")) then

!    Norm of complex ghost scalar field: sqrt(phiR**2 + phiI**2).

     if (allocated(complexghost_phi_norm)) then
        complexghost_phi_norm = sqrt(complexghost_phiR**2 + complexghost_phiI**2)
     end if

  end if


! ***********************
! ***   PROCA FIELD   ***
! ***********************

  if (contains(mattertype,"proca")) then

!    Proca Gauss constraint:  Div(E) + m^2 Phi = 0.

     if (allocated(Cproca)) then
        Cproca = D1_procaE + procaE*(half*D1_AB2/AB2 &
               + 6.d0*D1_phi + two/r) + proca_mass**2*procaPhi
     end if

!    Proca eigenfields.

     if (allocated(wp_proca)) then
        wp_proca = procaPhi + procaA/(sqrt(A)*psi2)
     end if

     if (allocated(wm_proca)) then
        wm_proca = procaPhi - procaA/(sqrt(A)*psi2)
     end if

  end if


! *******************************
! ***   COMPLEX PROCA FIELD   ***
! *******************************

  if (contains(mattertype,"complexproca")) then

!    Complex Proca Gauss constraint:  Div(E) + m^2 Phi = 0.
!    Real part.

     if (allocated(Ccomplexproca_R)) then

         Ccomplexproca_R = D1_cprocaE_R + cprocaE_R*(half*D1_AB2/AB2 &
               + 6.d0*D1_phi + two/r) + cproca_mass**2*cprocaPhi_R

!        Aditional term added for non-zero angular momentum of constituent fields.

         if (cproca_l/=0) then
             Ccomplexproca_R = Ccomplexproca_R - dble(cproca_l*(cproca_l+1))*cprocaXi_R/(B*psi4)/r**2
         end if

!        Charge terms if needed.  For a charged Proca field we need
!        to add a term of the form:
!
!        Ccomplexproca_R = Ccomplexproca_R - q eAr cprocaE_I

         if (cproca_q/=0.d0) then
            Ccomplexproca_R = Ccomplexproca_R - cproca_q*eAr*cprocaE_I
         end if

     end if

!    Imaginary part.

     if (allocated(Ccomplexproca_I)) then

         Ccomplexproca_I = D1_cprocaE_I + cprocaE_I*(half*D1_AB2/AB2 &
               + 6.d0*D1_phi + two/r) + cproca_mass**2*cprocaPhi_I
               
!        Aditional term added for non-zero angular momentum of constituent fields.

         if (cproca_l/=0) then
             Ccomplexproca_I = Ccomplexproca_I - dble(cproca_l*(cproca_l+1))*cprocaXi_I/(B*psi4)/r**2
         end if

!        Charge terms if needed.  For a charged Proca field we need
!        to add a term of the form:
!
!        Ccomplexproca_I = Ccomplexproca_I + q eAr cprocaE_R

         if (cproca_q/=0.d0) then
             Ccomplexproca_I = Ccomplexproca_I + cproca_q*eAr*cprocaE_R
         end if

     end if

!    Norms.

     if (allocated(cprocaPhi_norm)) then
        cprocaPhi_norm = sqrt(cprocaPhi_R**2 + cprocaPhi_I**2)
     end if

     if (allocated(cprocaA_norm)) then
        cprocaA_norm = sqrt(cprocaA_R**2 + cprocaA_I**2)
     end if

     if (allocated(cprocaB_norm)) then
        cprocaA_norm = sqrt(cprocaB_R**2 + cprocaB_I**2)
     end if

     if (allocated(cprocaE_norm)) then
        cprocaE_norm = sqrt(cprocaE_R**2 + cprocaE_I**2)
     end if

     if (allocated(cprocaXi_norm)) then
        cprocaE_norm = sqrt(cprocaXi_R**2 + cprocaXi_I**2)
     end if

!    Proca charge density times r**2.

     if (allocated(cproca_Qdens_r2)) then
        cproca_Qdens_r2 = cproca_Qdens*r**2
     end if

!    Call subroutine for integrated Proca charge.

     if (allocated(cproca_Qint)) then

        call procaintegral

!       Output total binding energy, but only at t=0 and
!       for certain type of initial data.

        if ((t(0)==0.d0).and.(rank==0)) then
           if (index(adjustl(idata),"procastar")/=0) then
              write(*,'(A,E19.12)') ' Binding energy (M-m*Q) = ',mass_int(0,Nr) - cproca_mass*cproca_Qint(0,Nr)
           end if
        end if

     end if

  end if


! ***********************
! ***   DIRAC FIELD   ***
! ***********************

  if (contains(mattertype,"dirac")) then

!    Scalar eigenfields.

     if (allocated(wpR_dirac)) then
        wpR_dirac = dirac_FR + dirac_GR
     end if

     if (allocated(wpI_dirac)) then
        wpI_dirac = dirac_FI + dirac_GI
     end if

     if (allocated(wmR_dirac)) then
        wmR_dirac = dirac_FR - dirac_GR
     end if

     if (allocated(wmI_dirac)) then
        wmI_dirac = dirac_FI - dirac_GI
     end if

!    Total conserved charge.

!    The total boson number (integrated boson charge) is then:
!
!                        /r             2
!    dirac_Nint  =  4 pi |  dirac_dens R  dR
!                        /0
!
!    with R the Schwarzschild radius.

     if (allocated(dirac_Nint)) then

        call diracintegral

!       Output total binding energy, but only at t=0 and
!       for certain type of initial data.

        if ((t(0)==0.d0).and.(rank==0)) then
           if (idata=="diracstar") then
              write(*,'(A,E19.12)') ' Binding energy (M-m*Q) = ',mass_int(0,Nr) - dirac_mass*dirac_Nint(0,Nr)
           end if
        end if

     end if

  end if


! **************************
! ***   ELECTRIC FIELD   ***
! **************************

! Calculation of quantities related to the electric field.
!
! Note:  This section must always be AFTER any type of field
! (scalar, complex, Proca, Dirac, etc.) that might have an
! electric charge.

  if (contains(mattertype,"electric")) then

!    Call subroutine for integrated electric charge.

     call chargeintegral

!    Calculate electric charge from surface integral (divergence theorem).

     eQ_surf = r**2*sqrt(A)*B*psi**6*electric

     if (rank==0) then
        do l=0,Nl-1
           do i=1,ghost
              eQ_surf(l,1-i) = eQ_surf(l,i)
           end do
        end do
     end if

!    Electric constraint (Gauss equation).

     if (allocated(Celectric)) then
        Celectric = D1_electric + electric*(half*D1_AB2/AB2 &
                  + 6.d0*D1_phi + two/r) - 4.d0*smallpi*echarge
     end if

!    Electric eigenfields.

     if (allocated(wp_electric)) then
        wp_electric = ePhi + eAr/(sqrt(A)*psi2)
     end if

     if (allocated(wm_electric)) then
        wm_electric = ePhi - eAr/(sqrt(A)*psi2)
     end if


!    ***********************************
!    ***   REISSNER-NORDSTROM MASS   ***
!    ***********************************

!    This is the Reissner-Nordstrom mass.  It is obtained by noticing
!    that for the Reissner-Nordstrom metric we have:
!
!    g   =  1 / (1 - 2M/R + (Q/R)^2)
!     RR
!
!    with R=r_area the Schwarzschild or areal radius (see above).
!    This implies:
!
!                                                         2       4
!    M  =  R/2 ( 1 - 1/g  + Q^2/R^2)  =  R/2 ( 1 - (dR/dr) / A psi  + (Q/R)^2)
!                       RR
!
!    with:
!
!               1/2    2
!    dR/dr  =  B    psi  ( 1 + r d B / (2 B)  +  2 r d psi / psi )
!                                 r                   r
!
!    So that finally:
!
!               2  1/2                                                  2
!    M  =  r psi  B   / 2  [ 1 - (B/A) ( 1 + r d B / (2 B) + 2 r d phi )  + (Q/R)**2]
!                                               r                 r
!
!    Note:  This must ALWAYS be after the section for the elecgric field.

     if (allocated(mass_rn)) then


!       Check if eQ_int is allocated.

        if (allocated(eQ_int)) then

           mass_rn = half*r_area*(one - B/A*(one + r*(half*D1_B/B + two*D1_phi))**2 + (eQ_int/r_area)**2)

!       If eQ_int is not allocated, check if the initial data
!       is Reissner-Nordstrom for which Q is constant.

        else if (idata=="reissnernordstrom") then

           mass_rn = half*r_area*(one - B/A*(one + r*(half*D1_B/B + two*D1_phi))**2 + (BHcharge/r_area)**2)

!       Otherwise stop.

        else

           print *, 'In order to calculate the Reissner-Nodstrom mass "mass_rn"'
           print *, 'you also need to output the integrated charge "eQ_int".'
           print *, 'Aborting!  (subroutine analysis_matter)'
           print *
           call die

        end if

!       Fix boundary.

        mass_rn(:,Nr) = 2.d0*mass_rn(:,Nr-1) - mass_rn(:,Nr-2)

!       Output total RN mass, but only at t=0 and
!       for certain type of initial data.

        if ((t(0)==0).and.(rank==0)) then
           if ((idata=="reissnernordstrom").or.(idata=="chargedboson").or.(idata=="chargedproca"))  then
              write(*,'(A,E14.8)') ' Total Reissner-Nordstrom mass = ',mass_rn(0,Nr)
              print *
           end if
        end if

!       Output total binding energy, but only at t=0 and
!       for certain type of initial data.

        if ((t(0)==0.d0).and.(rank==0)) then
           if (idata=="chargedboson") then
              write(*,'(A,E19.12)') ' Binding energy (M-m*Q) = ',mass_rn(0,Nr) - complex_mass*complex_NB(0,Nr)
           else if (idata=="chargedproca") then
              write(*,'(A,E19.12)') ' Binding energy (M-m*Q) = ',mass_rn(0,Nr) - cproca_mass*cproca_Qint(0,Nr)
           end if
        end if

!       Restrict.

        restrictvar => mass_rn

        do l=Nl-1,1,-1
           call restrict(l,.false.)
        end do

     end if

  end if


! *********************
! ***   COSMOLOGY   ***
! *********************

! For cosmological runs, we define some quantities
! subtracting the background for easier analysis.

  if (cosmic_run) then

     do l=0,Nl-1

!       Perturbed density.

        rho_pert(l,:) = rho(l,:) - cosmobg_rho(l)

!       Density constrast.

        rho_contrast(l,:) = (rho(l,:) - cosmobg_rho(l))/cosmobg_rho(l)

     end do

!    Integrated mass and density constrast.

     call cosmo_massintegral

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine analysis_matter


