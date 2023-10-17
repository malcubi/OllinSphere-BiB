
  subroutine sources_dirac(l)

! ************************************
! ***   SOURCES FOR SCALAR FIELD   ***
! ************************************

! This routine calculates the sources for
! a Dirac field.  The dynamical variables
! are: dirac_FR, dirac_FI, dirac_GR, dirac_GI.
!
! The evolution equations for these quantities come
! directly from the Dirac equation:
!
! dFR/dt  =  -  alpha / (sqrt(A) psi**2) { d GR + GR [ d alpha / 2 alpha + d B / 2B + 2 d phi
!                                           r           r                   r            r
!
!            + ( 1 + sqrt(A/B) ) / r ] }  +  alpha ( trK FR / 2  +  m FI )
!
!
! dFI/dt  =  -  alpha / (sqrt(A) psi**2) { d GI + GI [ d alpha / 2 alpha + d B / 2B + 2 d phi
!                                           r           r                   r            r
!
!            + ( 1 + sqrt(A/B) ) / r ] }  +  alpha ( trK FI / 2  -  m FR )
!
!
! dGR/dt  =  -  alpha / (sqrt(A) psi**2) { d FR + FR [ d alpha / 2 alpha + d B / 2B + 2 d phi
!                                           r           r                   r            r
!
!            + ( 1 - sqrt(A/B) ) / r ] }  +  alpha ( trK GR / 2  -  m GI )
!
!
! dGI/dt  =  -  alpha / (sqrt(A) psi**2) { d FI + FI [ d alpha / 2 alpha + d B / 2B + 2 d phi
!                                           r           r                   r            r
!
!            + ( 1 - sqrt(A/B) ) / r ] }  +  alpha ( trK GI / 2  +  m GR )
!
!
! Notice that since (FR,FI) are even, the terms with  F*(1-sqrt(A/B))/r  in the
! evolution equations for (GR,GI) would seem to be singular.  However, they
! can be written in a regular way in terms of lambda as:
!
! (1 - sqrt(A/B))/r  =  (1 - A/B) / (1 + sqrt(A/B)) / r  = r lambda / (1 + sqrt(A/B))
!
! This is what we do below.

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer l

  real(8) DFR,DFI


! *******************
! ***   SOURCES   ***
! *******************

! Sources for (FR,FI).

  sdirac_FR(l,:) = - alpha(l,:)/(sqrt(A(l,:))*psi2(l,:))*(D1_dirac_GR(l,:) &
                 + dirac_GR(l,:)*(0.5d0*D1_alpha(l,:)/alpha(l,:) + 0.5d0*D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:) &
                 + (1.d0 + sqrt(A(l,:)/B(l,:)))/r(l,:))) &
                 + alpha(l,:)*(0.5d0*trK(l,:)*dirac_FR(l,:) + dirac_mass*dirac_FI(l,:))

  sdirac_FI(l,:) = - alpha(l,:)/(sqrt(A(l,:))*psi2(l,:))*(D1_dirac_GI(l,:) &
                 + dirac_GI(l,:)*(0.5d0*D1_alpha(l,:)/alpha(l,:) + 0.5d0*D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:) &
                 + (1.d0 + sqrt(A(l,:)/B(l,:)))/r(l,:))) &
                 + alpha(l,:)*(0.5d0*trK(l,:)*dirac_FI(l,:) - dirac_mass*dirac_FR(l,:))

! Sources for (GR,GI).

  sdirac_GR(l,:) = - alpha(l,:)/(sqrt(A(l,:))*psi2(l,:))*(D1_dirac_FR(l,:) &
                 + dirac_FR(l,:)*(0.5d0*D1_alpha(l,:)/alpha(l,:) + 0.5d0*D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:) &
                 + r(l,:)*lambda(l,:)/(1.d0 + sqrt(A(l,:)/B(l,:))))) &
                 + alpha(l,:)*(0.5d0*trK(l,:)*dirac_GR(l,:) - dirac_mass*dirac_GI(l,:))

  sdirac_GI(l,:) = - alpha(l,:)/(sqrt(A(l,:))*psi2(l,:))*(D1_dirac_FI(l,:) &
                 + dirac_FI(l,:)*(0.5d0*D1_alpha(l,:)/alpha(l,:) + 0.5d0*D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:) &
                 + r(l,:)*lambda(l,:)/(1.d0 + sqrt(A(l,:)/B(l,:))))) &
                 + alpha(l,:)*(0.5d0*trK(l,:)*dirac_GI(l,:) + dirac_mass*dirac_GR(l,:))

! Shift terms.

  if (shift/="none") then
     sdirac_FR(l,:) =  sdirac_FR(l,:) + beta(l,:)*DA_dirac_FR(l,:)
     sdirac_FI(l,:) =  sdirac_FI(l,:) + beta(l,:)*DA_dirac_FI(l,:)
     sdirac_GR(l,:) =  sdirac_GR(l,:) + beta(l,:)*DA_dirac_GR(l,:)
     sdirac_GI(l,:) =  sdirac_GI(l,:) + beta(l,:)*DA_dirac_GI(l,:)
  end if

! Dissipation.

  if (diracdiss/=0.d0) then

     dissipvar => dirac_FR
     sourcevar => sdirac_FR
     call dissipation(l,+1,diracdiss)

     dissipvar => dirac_FI
     sourcevar => sdirac_FI
     call dissipation(l,+1,diracdiss)

     dissipvar => dirac_GR
     sourcevar => sdirac_GR
     call dissipation(l,-1,diracdiss)

     dissipvar => dirac_GI
     sourcevar => sdirac_GI
     call dissipation(l,-1,diracdiss)

  end if


! **********************
! ***   BOUNDARIES   ***
! **********************

! At the moment I use a simple outgoing wave boundary condition
! for (FR,FI) assuming that far away they behave as f(r-t)/r.
! This implies a condition of the form:
!
! dF/dt  ~  - (dF/dr + F/r)
!
! Notice that we do not need a boundary condition for (GR,GI),
! in fact trying to impose a boundary condition for them is
! inconsistent.

  if ((boundtype/="static").and.(boundtype/="flat")) then
     sdirac_FR(l,Nr) = - D1_dirac_FR(l,Nr) - dirac_FR(l,Nr)/r(l,Nr)
     sdirac_FI(l,Nr) = - D1_dirac_FI(l,Nr) - dirac_FI(l,Nr)/r(l,Nr)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_dirac
