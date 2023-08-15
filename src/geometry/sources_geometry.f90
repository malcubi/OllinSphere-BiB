!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/sources_geometry.f90,v 1.72 2023/01/27 16:21:27 malcubi Exp $

  subroutine sources_geometry(l)

! *******************************************
! ***   SOURCES FOR EVOLUTION EQUATIONS   ***
! *******************************************

! This routine calculates the sources for the evolution
! equations of the different geometrical variables.

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  logical contains

  integer i,l

  real(8) sigma                 ! BSSN flavor: 0 for eulerian, 1 for lagrangian.
  real(8) zero,sixth,third
  real(8) half,one,two,smallpi
  real(8) aux,r0,interp


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  half  = 0.5d0
  third = 1.d0/3.d0
  sixth = 1.d0/6.d0

  one   = 1.d0
  two   = 2.d0

  smallpi = acos(-one)


! ***********************
! ***   BSSN FLAVOR   ***
! ***********************

! Depending on the choice for the evolution equation
! of the conformal volume elements, we have two
! different BSSN schemes, Eulerian and Lagrangian,
! that change the form of the shift terms in the
! sources of the geometric variables.

  if (bssnflavor=='eulerian') then
     sigma = 0.0
  else if (bssnflavor=='lagrangian') then
     sigma = 1.0
  end if


! *********************************
! ***   PROPER TIME AT ORIGIN   ***
! *********************************

! The source for the proper time at the origin
! is just the lapse.

  r0 = zero
  interpvar => alpha
  aux = interp(l,r0,.false.)

  stau_origin(l) = aux


! ****************************************
! ***   SOURCES FOR CONFORMAL FACTOR   ***
! ****************************************

! Sources for conformal factor.  Notice that if
! chimethod is true we evolve chi instead of phi.
!
! Notice also that sphi must always be defined, even
! if we are using chimethod=true.  This is because
! it might later used to construct other sources.

! Sources for phi.

  sphi(l,:) = - sixth*alpha(l,:)*trK(l,:)

  if (formulation=="z4c") then
     sphi(l,:) = sphi(l,:) - two*sixth*alpha(l,:)*z4theta(l,:)
  end if

  if (shift/="none") then
     sphi(l,:) = sphi(l,:) + beta(l,:)*DA_phi(l,:) + sigma*sixth*DIV_beta(l,:)
  end if

  !if (geodiss/=0.d0) then
  !   dissipvar => phi
  !   sourcevar => sphi
  !   call dissipation(l,+1,geodiss)
  !end if

! Sources for chi = 1/psi**n = exp(-n*phi).

  if (chimethod) then

     schi(l,:) = dble(chipower)*sixth*alpha(l,:)*chi(l,:)*trK(l,:)

     if (formulation=="z4c") then
        schi(l,:) = schi(l,:) + two*dble(chipower)*sixth*alpha(l,:)*chi(l,:)*z4theta(l,:)
     end if

     if (shift/="none") then
        schi(l,:) = schi(l,:) + beta(l,:)*DA_chi(l,:) &
                  - sigma*dble(chipower)*sixth*chi(l,:)*DIV_beta(l,:)
     end if

     !if (geodiss/=0.d0) then
     !   dissipvar => chi
     !   sourcevar => schi
     !   call dissipation(l,+1,geodiss)
     !end if

  end if


! ****************************************
! ***   SOURCES FOR CONFORMAL METRIC   ***
! ****************************************

! Sources for conformal metric.

  sA(l,:) = - two*alpha(l,:)*A(l,:)*KTA(l,:)
  sB(l,:) = - two*alpha(l,:)*B(l,:)*KTB(l,:)

! Shift terms.

  if (shift/="none") then
     sA(l,:) = sA(l,:) + beta(l,:)*DA_A(l,:) + two*A(l,:)*D1_beta(l,:) &
             - two*third*sigma*A(l,:)*DIV_beta(l,:)
     sB(l,:) = sB(l,:) + beta(l,:)*DA_B(l,:) + two*beta(l,:)*B(l,:)/r(l,:) &
             - two*third*sigma*B(l,:)*DIV_beta(l,:)
  end if

! Dissipation.

  !if (geodiss/=0.d0) then
  !   dissipvar => A
  !   sourcevar => sA
  !   call dissipation(l,+1,geodiss)
  !   dissipvar => B
  !   sourcevar => sB
  !   call dissipation(l,+1,geodiss)
  !end if


! ****************************************************
! ***   SOURCES FOR TRACE OF EXTRINSIC CURVATURE   ***
! ****************************************************

  if (slicing=="maximal") then

!    For maximal slicing the source for trK is set to zero.

     strK(l,:) = zero

  else

!    For any other slicing we calculate the full source term:
!    which is given by:
!
!                 2                 ij
!    strK  =  -  D alpha  +  alpha K  K   +  4 pi alpha ( rho + S )
!                                      ij
!
!           2
!    where D  is the physical Laplacian.
!
!    Notice also that the square of Kij is in fact calculated as an
!    auxiliary variable named "K2" in auxiliary_geometry.f90.

!    Terms coming from Laplacian of lapse.

     strK(l,:) = - Lapla_alpha(l,:)

!    Quadratic terms.

     strK(l,:) = strK(l,:) + alpha(l,:)*K2(l,:)

!    Damping term for Z4c.

     if (formulation=="z4c") then
        strK(l,:) = strK(l,:) + alpha(l,:)*kappa1*(one - kappa2)*z4theta(l,:)
     end if

!    Shift terms.

     if (shift/="none") then
        strK(l,:) = strK(l,:) + beta(l,:)*DA_trK(l,:)
     end if

!    Matter terms.

     if (mattertype/="vacuum") then
        if (contains(mattertype,"nonmin")) then
           strK(l,:) = strK(l,:) + 0.5d0/nonmin_f(l,:)*alpha(l,:)*(rho(l,:) + trS(l,:))
        else
           strK(l,:) = strK(l,:) + 4.d0*smallpi*alpha(l,:)*(rho(l,:) + trS(l,:))
        end if
     end if

!    Dissipation.

     if (geodiss/=0.d0) then
        dissipvar => trK
        sourcevar => strK
        call dissipation(l,+1,geodiss)
     end if

  end if


! *****************************************************
! ***   SOURCES FOR TRACELESS EXTRINSIC CURVATURE   ***
! *****************************************************

! The full source term is:
!
!                   r                 2
!    sKTA  =  -  ( D D alpha - (1/3) D alpha )
!                     r
!
!                     r
!          + alpha ( R  - (1/3) R )  +  alpha trK KTA
!                     r
!
!
!          - (16/3) pi alpha  ( SAA - SBB )
!
!
! where D denotes physical covariant derivatives, R^r_r
! is the mixed radial component of the physical Ricci 
! tensor and R is the physical scalar curvature.
!
! Do notice that there is a typo in the paper with Martha Mendez
! (M. Alcubierre and M. Mendez, Gen. Rel. Grav. 43, p. 2769 (2011))
! in equation (102). The last term of that equation in the paper
! has 16*pi instead of the correct coefficient (16/3)*pi.

  if (.true.) then
  !if (.not.noKTA) then

!    Terms coming from second derivatives of the lapse.
!    The expression used below is equivalent to the one commented out,
!    but the numerical error is considerable smaller since we have
!    compacted the terms using: DD_alphar = d_r (d_r alpha / r psi4)

!    sKTA(l,:) = - Dcov2_alpha(l,:) + third*Lapla_alpha(l,:)
     sKTA(l,:) = - 2.d0*third/A(l,:)*(r(l,:)*DD_alphar(l,:) &
               - half*D1_alpha(l,:)/psi4(l,:)*(D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:)))

!    Terms coming from the (physical) Ricci tensor.
!    As before, the expression used below is equivalent to the one commented
!    out, but the numerical error is considerable smaller since we have
!    compacted the terms using: DD_phir = d_r (d_r phi / r psi2)

!    sKTA(l,:) = sKTA(l,:) + alpha(l,:)*(RICA(l,:) - third*RSCAL(l,:))

     !$OMP PARALLEL DO SCHEDULE(GUIDED)
     do i=1-ghost,Nrmax
        sKTA(l,i) = sKTA(l,i) - alpha(l,i)/(A(l,i)*psi4(l,i)) &
             *(third*(D2_A(l,i)/A(l,i) - D2_B(l,i)/B(l,i) - 2.d0*A(l,i)*D1_Deltar(l,i) &
             - 2.d0*(D1_A(l,i)/A(l,i))**2 + (D1_B(l,i)/B(l,i))**2) &
             + half*D1_A(l,i)*D1_B(l,i)/A(l,i)/B(l,i) &
             + (D1_A(l,i)/A(l,i) - 4.d0*third*A(l,i)*D1_B(l,i)/B(l,i)**2)/r(l,i) &
             + two*third*lambda(l,i) + 4.d0*third*(r(l,i)*psi2(l,i)*DD_phir(l,i) &
             - half*D1_phi(l,i)*(D1_A(l,i)/A(l,i) + D1_B(l,i)/B(l,i))))
     end do
     !$OMP END PARALLEL DO

!    Quadratic terms.

     sKTA(l,:) = sKTA(l,:) + alpha(l,:)*trK(l,:)*KTA(l,:)

!    Shift terms.

     if (shift/="none") then
        sKTA(l,:) = sKTA(l,:) + beta(l,:)*DA_KTA(l,:)
     end if

!    Matter terms.

     if (mattertype/="vacuum") then
        if (contains(mattertype,"nonmin")) then
           sKTA(l,:) = sKTA(l,:) - 2.d0*third/nonmin_f(l,:)*alpha(l,:)*(SAA(l,:) - SBB(l,:))
        else
           sKTA(l,:) = sKTA(l,:) - 16.d0*third*smallpi*alpha(l,:)*(SAA(l,:) - SBB(l,:))
        end if
     end if

!    Dissipation.

     if (geodiss/=0.d0) then
        dissipvar => KTA
        sourcevar => sKTA
        call dissipation(l,+1,geodiss)
     end if

  end if


! **********************************
! ***   SOURCES FOR BSSN DELTA   ***
! **********************************

  if (.not.noDeltar) then

!    From the definition of Deltar one obtains the following
!    source term:
!
!    sDeltar  =  alpha ( 2 KTA Deltar - 4 r Klambda / B )
!
!             -  (2/A) ( KTA d alpha + alpha d KTA )
!                             r               r
!
!    To this one can add a multiple of the momentum constraint:
!                                
!    sDeltar  =  sDeltar + eta alpha /A  ( d KTA - (2/3) d trK
!                                           r             r
!                                2
!             + 6 KTA d phi  +  r ( d B / B + 2/r ) Klambda  - 8 pi JA)
!                      r             r

!    Terms coming from definition of Deltar.

     sDeltar(l,:) = alpha(l,:)*(two*KTA(l,:)*Deltar(l,:) - 4.0d0*r(l,:)*Klambda(l,:)/B(l,:)) &
             - two/A(l,:)*(KTA(l,:)*D1_alpha(l,:) + alpha(l,:)*D1_KTA(l,:))

!    Terms coming from adding the momentum constraint.

     sDeltar(l,:) = sDeltar(l,:) + eta*alpha(l,:)/A(l,:)*(D1_KTA(l,:) - two*third*D1_trK(l,:) &
             + 6.d0*KTA(l,:)*D1_phi(l,:) + r(l,:)*(r(l,:)*D1_B(l,:)/B(l,:) + two)*Klambda(l,:))

!    Terms for Z4c formulation. This term seems to cause an instability,
!    so I commented it out. Not sure why this is, or what "theoretical"
!    consequences this might have.

     !if (formulation=="z4c") then
        !sDeltar(l,:) = sDeltar(l,:) - third*eta*alpha(l,:)*D1_z4theta(l,:) &
        !             - two*alpha(l,:)*kappa1*(Deltar(l,:) - DeltaAB(l,:))
     !end if

!    Shift terms. Here we add the shift terms:
!
!    sDeltar  =  sDeltar  +  beta d Deltar  -  Deltar d beta
!                                  r                   r
!                  2
!             +  (d beta) / A  +  (2/B) d (beta/r)
!                  r                     r
!
!             +  (sigma/3) [ d (div(beta)) / A  +  2 Deltar div(beta) ]
!                             r

     if (shift/="none") then
        !$OMP PARALLEL DO SCHEDULE(GUIDED)
        do i=1-ghost,Nrmax
           sDeltar(l,i) = sDeltar(l,i) + beta(l,i)*DA_Deltar(l,i) - Deltar(l,i)*D1_beta(l,i) &
                + D2_beta(l,i)/A(l,i) + two/B(l,i)*DD_beta(l,i) &
                + sigma*third*(D1_DIV_beta(l,i)/A(l,i) + two*Deltar(l,i)*DIV_beta(l,i))
        end do
        !$OMP END PARALLEL DO
     end if

!    Matter terms.  The matter contributions only appear
!    if we have added a multiple of the momentum constraint.

     if (mattertype/="vacuum") then
        if (contains(mattertype,"nonmin")) then
           sDeltar(l,:) = sDeltar(l,:) - eta/nonmin_f(l,:)*alpha(l,:)/A(l,:)*JA(l,:)
        else
           sDeltar(l,:) = sDeltar(l,:) - 8.d0*smallpi*eta*alpha(l,:)/A(l,:)*JA(l,:)
        end if
     end if

!    Dissipation.

     if (geodiss/=0.d0) then
        dissipvar => Deltar
        sourcevar => sDeltar
        call dissipation(l,-1,geodiss)
     end if

  end if


! ******************************
! ***   SOURCES FOR LAMBDA   ***
! ******************************

! For regularization we introduce lambda and lambda2
! defined as:
!
! lambda  := (1 - A/B)/r**2
! lambda2 := lambda/psi**N

  if (.not.nolambda) then

     if (.not.regular2) then

!       Sources for lambda:
!
!       slambda  =  2 alpha (A/B)  Klambda
!
!
!                +  beta d lambda  +  2/r ( beta lambda - (A/B) d (beta/r) )
!                         r                                      r

        slambda(l,:) = two*alpha(l,:)*A(l,:)/B(l,:)*Klambda(l,:)

        if (shift/="none") then
           slambda(l,:) = slambda(l,:) + beta(l,:)*DA_lambda(l,:) &
                        + two/r(l,:)*(beta(l,:)*lambda(l,:) &
                        - (A(l,:)/B(l,:))*DD_beta(l,:))
        end if

     else

!       Sources for lambda2:
!
!       slambda  =  2 alpha (A/B)  Klambda2  +  (N/6) alpha lambda2 trK
!
!
!                +  beta d lambda2  +  2/r ( beta lambda2 - (A/B/psi**N) d (beta/r) )
!                         r                                               r
!
!                -  N/6 sigma lambda2 DIV_beta

        slambda2(l,:) = two*alpha(l,:)*A(l,:)/B(l,:)*Klambda2(l,:) &
                      + (lambdapower/6.d0)*alpha(l,:)*lambda2(l,:)*trK(l,:)

        if (shift/="none") then
           slambda2(l,:) = slambda2(l,:) + beta(l,:)*DA_lambda2(l,:) &
                         + two/r(l,:)*(beta(l,:)*lambda2(l,:) &
                         - (A(l,:)/B(l,:)/psi(l,:)**lambdapower)*DD_beta(l,:)) &
                         - lambdapower/6.d0*sigma*lambda2(l,:)*DIV_beta(l,:)
        end if

     end if

  end if


! *******************************
! ***   SOURCES FOR KLAMBDA   ***
! *******************************

! For regularization we introduce Klambda and Klambda2
! defined as:
!
! Klambda  := (KTA - KTB)/r**2 = (3/2) KTA / r**2
! Klambda2 := Klambda/psi**N

  if (.not.nolambda) then

     if (.not.regular2) then

!       Sources for Klambda.  They are separated by type of contribution.
!
!       I) Terms coming from derivatives of lapse and conformal factor:
!
!                                                     4                   2                    2
!       sKlambda  =  - 1 / (r A) [ d ( d alpha / r psi ) + ( 2 alpha / psi  ) d ( d phi / r psi )
!                                   r   r                                      r   r
!                             4
!                 - (1 / r psi ) ( d alpha / 2 + alpha d phi ) ( d A / A + d B / B ) ]
!                                   r                   r         r         r

        sKlambda(l,:) = - one/r(l,:)/A(l,:)*(DD_alphar(l,:) + (two*alpha(l,:)/psi2(l,:))*DD_phir(l,:) &
              - (half*D1_alpha(l,:) + alpha(l,:)*D1_phi(l,:))*(D1_A(l,:)/A(l,:) &
              + D1_B(l,:)/B(l,:))/r(l,:)/psi4(l,:))

!       II) Terms coming from the conformal Ricci tensor. These terms are
!       quite complicated since one needs to regularize them and rewrite
!       them in terms of lambda.
!
!                                               4            2
!       sKlambda  =  sKlambda  +  alpha / (A psi ) [ (B/2A) d lambda
!                                                            r
!                                                          
!                 + (A/r) d ( Deltar / r ) + d lambda / r ( 1 + 2 B/A - r B Deltar / 2 )
!                          r                  r
!                             2
!                 +  d A / ( r A )  ( (3/4) d A / A  -  d B / B )
!                     r                      r           r
!                                                                         2
!                 -  lambda / r  ( B Deltar  +  2 d B / B ) + (B/A) lambda  ]
!                                                  r

        sKlambda(l,:) = sKlambda(l,:) + alpha(l,:)/A(l,:)/psi4(l,:) &
              *(0.5d0*B(l,:)/A(l,:)*D2_lambda(l,:) + A(l,:)*DD_Deltar(l,:)/r(l,:) &
              + D1_lambda(l,:)/r(l,:)*(1.d0 + two*B(l,:)/A(l,:) - 0.5d0*r(l,:)*B(l,:)*Deltar(l,:)) &
              + D1_A(l,:)/(A(l,:)*r(l,:)**2)*(0.75d0*D1_A(l,:)/A(l,:) - D1_B(l,:)/B(l,:)) &
              - lambda(l,:)/r(l,:)*(B(l,:)*Deltar(l,:) + two*D1_B(l,:)/B(l,:)) + B(l,:)*lambda(l,:)**2/A(l,:))

!       III) Quadratic terms:  sKlambda = sKlambda + alpha trK Klambda

        sKlambda(l,:) = sKlambda(l,:) + alpha(l,:)*trK(l,:)*Klambda(l,:)

!       IV) Shift terms:
!
!       sKlambda  =  sKlambda  +  beta ( d Klambda  +  2 Klambda/r )
!                                         r
!                                               2
!                 =  sKlambda  +  (3/2) beta / r  ( d KTA )
!                                                    r

        if (shift/="none") then
!          sKlambda(l,:) = sKlambda(l,:) + beta(l,:)*(DA_Klambda(l,:) + two*Klambda(l,:)/r(l,:))
           sKlambda(l,:) = sKlambda(l,:) + 1.5d0*beta(l,:)*DA_KTA(l,:)/r(l,:)**2
        end if

!       V) Matter terms:  sKlambda = sKlambda - 8 pi alpha SLL

        if (mattertype/="vacuum") then
           if (contains(mattertype,"nonmin")) then
              sKlambda(l,:) = sKlambda(l,:) - alpha(l,:)/nonmin_f(l,:)*SLL(l,:)
           else
              sKlambda(l,:) = sKlambda(l,:) - 8.d0*smallpi*alpha(l,:)*SLL(l,:)
           end if
        end if

!       Dissipation.

        if (geodiss/=0.d0) then
           dissipvar => Klambda
           sourcevar => sKlambda
           call dissipation(l,+1,geodiss)
        end if

     else

!       Sources for Klambda2.  They are separated by type of contribution.

        !$OMP PARALLEL DO SCHEDULE(GUIDED)
        do i=1-ghost,Nrmax

!          I) Terms coming from derivatives of lapse and conformal factor:

           sKlambda2(l,i) = - one/r(l,i)/A(l,i)/psi(l,i)**lambdapower*(DD_alphar(l,i) &
                       + (two*alpha(l,i)/psi2(l,i))*DD_phir(l,i) &
                       - (half*D1_alpha(l,i) + alpha(l,i)*D1_phi(l,i)) &
                       *(D1_A(l,i)/A(l,i) + D1_B(l,i)/B(l,i))/r(l,i)/psi4(l,i))

!          II) Terms coming from the conformal Ricci tensor:

           sKlambda2(l,i) = sKlambda2(l,i) + alpha(l,i)/r(l,i)/exp((4.d0+lambdapower)*phi(l,i))*DD_Deltar(l,i) &
                       + alpha(l,i)/A(l,i)**2/psi4(l,i)*(B(l,i) &
                       *(half*lambdapower*lambda2(l,i)*D2_phi(l,i) + 0.5d0*D2_lambda2(l,i) &
                       + lambdapower*D1_lambda2(l,i)*D1_phi(l,i) + half*lambdapower**2*lambda2(l,i)*D1_phi(l,i)**2) &
                       + (A(l,i) + 2.d0*B(l,i))/r(l,i)*(lambdapower*lambda2(l,i)*D1_phi(l,i) + D1_lambda2(l,i)) &
                       + D1_A(l,i)/r(l,i)**2/psi(l,i)**lambdapower*(0.75d0*D1_A(l,i)/A(l,i) - D1_B(l,i)/B(l,i)) &
                       - 2.d0*A(l,i)/B(l,i)/r(l,i)*lambda2(l,i)*D1_B(l,i)) &
                       - alpha(l,i)*B(l,i)/A(l,i)/psi4(l,i)*Deltar(l,i) &
                       *(half*lambdapower*lambda2(l,i)*D1_phi(l,i) + half*D1_lambda2(l,i) + lambda2(l,i)/r(l,i)) &
                       + exp((lambdapower-4.d0)*phi(l,i))*alpha(l,i)*B(l,i)*(lambda2(l,i)/A(l,i))**2

!          II) Quadratic terms:

           sKlambda2(l,i) = sKlambda2(l,i) + (1.d0 + lambdapower/6.d0)*alpha(l,i)*trK(l,i)*Klambda2(l,i)

        end do
        !$OMP END PARALLEL DO

!       IV) Shift terms.

        if (shift/="none") then
           sKlambda2(l,:) = sKlambda2(l,:) + beta(l,:)*(DA_Klambda2(l,:) &
                          + 2.d0*Klambda2(l,:)/r(l,:)) &
                          - lambdapower/6.d0*sigma*Klambda2(l,:)*DIV_beta(l,:)
        end if

!       V) Matter terms:

        if (mattertype/="vacuum") then
           if (contains(mattertype,"nonmin")) then
              sKlambda2(l,:) = sKlambda2(l,:) - alpha(l,:)/nonmin_f(l,:)/psi(l,:)**lambdapower*SLL(l,:)
           else
              sKlambda2(l,:) = sKlambda2(l,:) - 8.d0*smallpi*alpha(l,:)/psi(l,:)**lambdapower*SLL(l,:)
           end if
        end if

!       Dissipation.

        if (geodiss/=0.d0) then
           dissipvar => Klambda2
           sourcevar => sKlambda2
           call dissipation(l,+1,geodiss)
        end if

     end if

  end if


! ********************************
! ***   SOURCES FOR Z4 THETA   ***
! ********************************

! Sources for z4theta. The source is essentially
! the Hamiltonian constraint.

  if (formulation=="z4c") then

!    Hamiltonian constraint term.

     sz4theta(l,:) = half*alpha(l,:)*(RSCAL(l,:) - (KTA(l,:)**2 + two*KTB(l,:)**2) &
                   + two*third*trK(l,:)**2)

     if (mattertype/="vacuum") then
        if (contains(mattertype,"nonmin")) then
           sz4theta(l,:) = sz4theta(l,:) - alpha(l,:)/nonmin_f(l,:)*rho(l,:)
        else
           sz4theta(l,:) = sz4theta(l,:) - 8.d0*smallpi*alpha(l,:)*rho(l,:)
        end if
     end if

!    Shift terms.

     if (shift/="none") then
        sz4theta(l,:) = sz4theta(l,:) + beta(l,:)*DA_z4theta(l,:)
     end if

!    Damping terms.

     sz4theta(l,:) = sz4theta(l,:) - alpha(l,:)*kappa1*(two + kappa2)*z4theta(l,:)

!    Dissipation.

     if (geodiss/=0.d0) then
        dissipvar => z4theta
        sourcevar => sz4theta
        call dissipation(l,+1,geodiss)
     end if

! For BSSN the source is set to 0.

  else

     sz4theta(l,:) = zero

  end if


! *********************
! ***   COSMOLOGY   ***
! *********************

! When needed, sources for cosmological background quantities.

  if (cosmic_run) then

!    Source for background proper time, dtau = alpha dt.

     scosmobg_tau(l) = cosmobg_alpha(l)

!    Source for background scale factor.  This is just the definition of H:
!
!    H  :=  d[ln(a)]/dt  =>  da/dt  =  a H
!
!    But notice that this derivative is with respect to cosmological comoving
!    time (i.e. proper time of the comoving observers).  This means that
!    with respect to the code time it actually becomes:
!
!    da/dt  =  alpha a H

     scosmobg_a(l) = cosmobg_alpha(l)*cosmobg_a(l)*cosmobg_H(l)

!    Source for Hubble parameter.  This is just the Friedmann equation, but
!    again with a term proportional to the lapse to convert from comoving
!    time to code time:
!                         2
!    dH/dt  =  - alpha [ H  + 4 pi ( rho/3 + p ) ]

     scosmobg_H(l) = - cosmobg_alpha(l)*(cosmobg_H(l)**2 &
                   + 4.d0*smallpi*(third*cosmobg_rho(l) + cosmobg_p(l)))

!    Source for background conformal factor.

     scosmobg_phi(l) = - sixth*cosmobg_alpha(l)*cosmobg_trK(l)

!    Source for background extrinsic curvature. Since K = -3H, this is just
!    essentally the same as the above equation except for a factor.

     scosmobg_trK(l) = - 3.d0*scosmobg_H(l)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_geometry

