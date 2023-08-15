!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/complexscalar.f90,v 1.28 2023/02/20 19:48:12 malcubi Exp $

  subroutine sources_complex(l)

! ********************************************
! ***   SOURCES FOR COMPLEX SCALAR FIELD   ***
! ********************************************

! This routine calculates the sources for a complex
! scalar field.  The Klein-Gordon equation with
! an arbitrary potential V(phi) has the form:
!
! Box (phi)  -  (dV/d|phi|**2) phi  =  0
!
! Remember that we have defined the space and time
! derivative arrays as:
!        
! pi  =  ( d phi  -  beta d phi ) / alpha
!           t              r
!
! xi  =  d phi
!         r
!
! The evolution equations for {phi,xi,pi} in first
! order form are:
!
!
! d phi  =  beta d phi  +  alpha pi
!  t              r
!
!
! d xi   =  beta d xi  +  xi d beta  +  alpha d pi  +  pi d alpha
!  t              r           r                r           r
!
!                                        4   /   
! d pi   =  beta d pi  +  ( alpha / A psi )  |  d xi 
!  t              r                          \   r 
!                                                               \
!        +  xi ( 2/r -  d A / 2A   +  d B / B  +  2 d ln(psi) ) |
!                        r             r             r          /
!                      4
!        +  ( 1 / A psi ) xi d alpha  +  alpha trK pi  -  alpha dV/dphi
!                             r
!
! If we want to evolve in second order form we just replace d_r(xi)
! in the last equation for the second derivative of phi.
!
! ELECTRO-MAGNETIC FIELD:
!
! Including the EM interaction terms only modifies the source of pi:
!
!                                       2    mu             mu               mu
! Box (phi)  =  (dV/d|phi|**2) phi  +  q A  A  phi  -  2iq A  d  phi  -  iq A      phi
!                                         mu                   mu              ;mu
!
! Note that choosing the Lorenz gauge cancels the last term.
!
! In 3+1 language we substitute the source for pi with:
!
!                          2      2        i                            i
!  (dV/d|phi|**2) phi  +  q (-ePhi  +  A  A ) phi  -  2iq (ePhi pi  +  A  xi )
!                                       i                                   i
!
! when decomposing the field in real and imaginary part the last term couples 
! the equations:
!
!                                         2      2        i                             i
!  Box (phiR)  = (dV/d|phi|**2) phiR  +  q (-ePhi  +  A  A ) phiR  +  2q (ePhi piI  +  A xiI )
!                                                      i                                    i
!
!                                         2      2        i                             i
!  Box (phiI)  = (dV/d|phi|**2) phiI  +  q (-ePhi  +  A  A ) phiI  -  2q (ePhi piR  +  A xiR )
!                                                      i                                    i


! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  logical contains

  integer i,l,sym


! *******************
! ***   SOURCES   ***
! *******************

! Source for scalar field.

  scomplex_phiR(l,:) = alpha(l,:)*complex_piR(l,:)
  scomplex_phiI(l,:) = alpha(l,:)*complex_piI(l,:)

! Source for radial derivative.

  scomplex_xiR(l,:) = alpha(l,:)*D1_complex_piR(l,:) + complex_piR(l,:)*D1_alpha(l,:)
  scomplex_xiI(l,:) = alpha(l,:)*D1_complex_piI(l,:) + complex_piI(l,:)*D1_alpha(l,:)

! Source for time derivative. Notice that there are
! two ways of doing this: using only derivatives of
! phi, or using the xi.

  if (complexmethod == "first") then

!    Term with derivative of xi.

     scomplex_piR(l,:) = alpha(l,:)/(A(l,:)*psi4(l,:))*D1_complex_xiR(l,:)
     scomplex_piI(l,:) = alpha(l,:)/(A(l,:)*psi4(l,:))*D1_complex_xiI(l,:)

!    Terms coming from Christoffel symbols.

     scomplex_piR(l,:) = scomplex_piR(l,:) + (alpha(l,:)*complex_xiR(l,:)*(2.d0/r(l,:) &
                       - 0.5d0*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:))  &
                       + complex_xiR(l,:)*D1_alpha(l,:))/(A(l,:)*psi4(l,:))
     scomplex_piI(l,:) = scomplex_piI(l,:) + (alpha(l,:)*complex_xiI(l,:)*(2.d0/r(l,:) &
                       - 0.5d0*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:))  &
                       + complex_xiI(l,:)*D1_alpha(l,:))/(A(l,:)*psi4(l,:))

  else

!    Term with second derivative of phi.

     !$OMP PARALLEL DO SCHEDULE(GUIDED)
     do i=1-ghost,Nrmax
        scomplex_piR(l,i) = alpha(l,i)/(A(l,i)*psi4(l,i))*D2_complex_phiR(l,i)
        scomplex_piI(l,i) = alpha(l,i)/(A(l,i)*psi4(l,i))*D2_complex_phiI(l,i)
     end do
     !$OMP END PARALLEL DO

!    Terms coming from Christoffel symbols.

     !$OMP PARALLEL DO SCHEDULE(GUIDED)
     do i=1-ghost,Nrmax
        scomplex_piR(l,i) = scomplex_piR(l,i) + (alpha(l,i)*D1_complex_phiR(l,i)*(2.d0/r(l,i) &
                          - 0.5d0*D1_A(l,i)/A(l,i) + D1_B(l,i)/B(l,i) + 2.d0*D1_phi(l,i))  &
                          + D1_complex_phiR(l,i)*D1_alpha(l,i))/(A(l,i)*psi4(l,i))
        scomplex_piI(l,i) = scomplex_piI(l,i) + (alpha(l,i)*D1_complex_phiI(l,i)*(2.d0/r(l,i) &
                          - 0.5d0*D1_A(l,i)/A(l,i) + D1_B(l,i)/B(l,i) + 2.d0*D1_phi(l,i))  &
                          + D1_complex_phiI(l,i)*D1_alpha(l,i))/(A(l,i)*psi4(l,i))
     end do
    !$OMP END PARALLEL DO

  end if

! Term proportional to trK.

  scomplex_piR(l,:) = scomplex_piR(l,:) + alpha(l,:)*trK(l,:)*complex_piR(l,:)
  scomplex_piI(l,:) = scomplex_piI(l,:) + alpha(l,:)*trK(l,:)*complex_piI(l,:)

! Potential term.

  if (complexpotential/="none") then
     scomplex_piR(l,:) = scomplex_piR(l,:) - alpha(l,:)*complex_VPR(l,:)
     scomplex_piI(l,:) = scomplex_piI(l,:) - alpha(l,:)*complex_VPI(l,:)
  end if

! Shift terms.

  if (shift/="none") then

!    Shift terms for phi.

     scomplex_phiR(l,:) = scomplex_phiR(l,:) + beta(l,:)*DA_complex_phiR(l,:)
     scomplex_phiI(l,:) = scomplex_phiI(l,:) + beta(l,:)*DA_complex_phiI(l,:)

!    Shift terms for xi.

     scomplex_xiR(l,:) = scomplex_xiR(l,:) + beta(l,:)*DA_complex_xiR(l,:) &
                       + complex_xiR(l,:)*D1_beta(l,:)
     scomplex_xiI(l,:) = scomplex_xiI(l,:) + beta(l,:)*DA_complex_xiI(l,:) &
                       + complex_xiI(l,:)*D1_beta(l,:)

!    Shift terms for pi.

     scomplex_piR(l,:) = scomplex_piR(l,:) + beta(l,:)*DA_complex_piR(l,:)
     scomplex_piI(l,:) = scomplex_piI(l,:) + beta(l,:)*DA_complex_piI(l,:)

  end if

! Angular momentum. For non-zero angular momentum
! we add to the Klein-Gordon equation the term:
!
!  - alpha l(l+1) phi / (r^2 B psi^4)

  if (complex_l>0) then

     scomplex_piR(l,:) = scomplex_piR(l,:) - alpha(l,:)*dble(complex_l*(complex_l+1)) &
                       * complex_phiR(l,:)/(B(l,:)*psi4(l,:))/r(l,:)**2

     scomplex_piI(l,:) = scomplex_piI(l,:) - alpha(l,:)*dble(complex_l*(complex_l+1)) &
                       * complex_phiI(l,:)/(B(l,:)*psi4(l,:))/r(l,:)**2

  end if

! Electromagnetic sources.

  if (contains(mattertype,"electric")) then

     scomplex_piR(l,:) = scomplex_piR(l,:) - alpha(l,:) &
          *(complex_q**2*(-ePhi(l,:)**2 + eAr(l,:)**2/(psi4(l,:)*A(l,:)))*complex_phiR(l,:)  &
          + 2.d0*complex_q*(ePhi(l,:)*complex_piI(l,:) + eAr(l,:)*complex_xiI(l,:)/(psi4(l,:)*A(l,:))))

     scomplex_piI(l,:) = scomplex_piI(l,:) - alpha(l,:) &
          *(complex_q**2*(-ePhi(l,:)**2 + eAr(l,:)**2/(psi4(l,:)*A(l,:)) )*complex_phiI(l,:) & 
          - 2.d0*complex_q*(ePhi(l,:)*complex_piR(l,:) + eAr(l,:)*complex_xiR(l,:)/(psi4(l,:)*A(l,:))))

  end if

! Dissipation.

  if (scalardiss/=0.d0) then

!    The symmetry depends on the value of complex_l.

     sym = int(1.d0-2.d0*mod(complex_l,2))

     dissipvar => complex_phiR
     sourcevar => scomplex_phiR
     call dissipation(l,sym,scalardiss)

     dissipvar => complex_phiI
     sourcevar => scomplex_phiI
     call dissipation(l,sym,scalardiss)

     dissipvar => complex_piR
     sourcevar => scomplex_piR
     call dissipation(l,sym,scalardiss)

     dissipvar => complex_piI
     sourcevar => scomplex_piI
     call dissipation(l,sym,scalardiss)

     dissipvar => complex_xiR
     sourcevar => scomplex_xiR
     call dissipation(l,-sym,scalardiss)

     dissipvar => complex_xiI
     sourcevar => scomplex_xiI
     call dissipation(l,-sym,scalardiss)

  end if

! Cosmological background quantities.

  if (cosmic_run) then

!    Sources for (phiR,phiI).

     scosmobg_complex_phiR(l) = cosmobg_alpha(l)*cosmobg_complex_piR(l)
     scosmobg_complex_phiI(l) = cosmobg_alpha(l)*cosmobg_complex_piI(l)

!    Sources for (piR,piI).

     scosmobg_complex_piR(l) = cosmobg_alpha(l)*cosmobg_trK(l)*cosmobg_complex_piR(l)
     scosmobg_complex_piI(l) = cosmobg_alpha(l)*cosmobg_trK(l)*cosmobg_complex_piI(l)

     if (complexpotential/="none") then
        scosmobg_complex_piR(l) = scosmobg_complex_piR(l) - cosmobg_alpha(l)*cosmobg_complex_VPR(l)
        scosmobg_complex_piI(l) = scosmobg_complex_piI(l) - cosmobg_alpha(l)*cosmobg_complex_VPI(l)
     end if

  end if


! **********************
! ***   BOUNDARIES   ***
! **********************

! For the boundary condition see the comment at the end
! of the routine "scalar.f90".  The condition here is
! the same.

  if ((boundtype/="static").and.(boundtype/="flat")) then

     if (.not.cosmic_run) then

        scomplex_piR(l,Nr) = - (D1_complex_piR(l,Nr) + complex_piR(l,Nr)/r(l,Nr))
        scomplex_piI(l,Nr) = - (D1_complex_piI(l,Nr) + complex_piI(l,Nr)/r(l,Nr))

!       Potential term.

        if (complexpotential/="none") then
           scomplex_piR(l,Nr) = scomplex_piR(l,Nr) - 0.5d0*complex_VPR(l,Nr)
           scomplex_piI(l,Nr) = scomplex_piI(l,Nr) - 0.5d0*complex_VPI(l,Nr)
        end if

!    Cosmological runs need to include the lapse and metric coefficients
!    that are ignored in the asymptotically flat case.

     else

        scomplex_piR(l,Nr) = scosmobg_complex_piR(l) - cosmobg_alpha(l)/cosmobg_a(l) &
                           *(D1_complex_piR(l,Nr) + (complex_piR(l,Nr)-cosmobg_complex_piR(l))/r(l,Nr))
        scomplex_piI(l,Nr) = scosmobg_complex_piI(l) - cosmobg_alpha(l)/cosmobg_a(l) &
                           *(D1_complex_piI(l,Nr) + (complex_piI(l,Nr)-cosmobg_complex_piI(l))/r(l,Nr))

!       Potential term.

        if (complexpotential/="none") then

        end if

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_complex

