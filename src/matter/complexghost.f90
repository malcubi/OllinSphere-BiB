!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/complexghost.f90,v 1.4 2023/03/02 19:49:55 malcubi Exp $

  subroutine sources_complexghost(l)

! *******************************************
! ***   SOURCES FOR COMPLEX GHOST FIELD   ***
! *******************************************

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
!        +  ( 1 / A psi ) xi d alpha  +  alpha trK pi  +  alpha dV/dphi
!                             r
!
!
! Notice that the convention of the code is such that the potential
! term in the Klein-Gordon has THE OPPOSITE sign as for a standard scalar
! field, but the potential is in fact negative, so we end up with
! exactly the same Klein-Gordon equation.
!
! If we want to evolve in second order form we just replace d_r(xi)
! in the last equation for the second derivative of phi.

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer i,l,sym


! *******************
! ***   SOURCES   ***
! *******************

! Source for scalar field.

  scomplexghost_phiR(l,:) = alpha(l,:)*complexghost_piR(l,:)
  scomplexghost_phiI(l,:) = alpha(l,:)*complexghost_piI(l,:)

! Source for radial derivative.

  scomplexghost_xiR(l,:) = alpha(l,:)*D1_complexghost_piR(l,:) + complexghost_piR(l,:)*D1_alpha(l,:)
  scomplexghost_xiI(l,:) = alpha(l,:)*D1_complexghost_piI(l,:) + complexghost_piI(l,:)*D1_alpha(l,:)

! Source for time derivative. Notice that there are
! two ways of doing this: using only derivatives of
! phi, or using the xi.

  if (complexghostmethod == "first") then

!    Term with derivative of xi.

     scomplexghost_piR(l,:) = alpha(l,:)/(A(l,:)*psi4(l,:))*D1_complexghost_xiR(l,:)
     scomplexghost_piI(l,:) = alpha(l,:)/(A(l,:)*psi4(l,:))*D1_complexghost_xiI(l,:)

!    Terms coming from Christoffel symbols.

     scomplexghost_piR(l,:) = scomplexghost_piR(l,:) + (alpha(l,:)*complexghost_xiR(l,:)*(2.d0/r(l,:) &
                       - 0.5d0*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:))  &
                       + complexghost_xiR(l,:)*D1_alpha(l,:))/(A(l,:)*psi4(l,:))
     scomplexghost_piI(l,:) = scomplexghost_piI(l,:) + (alpha(l,:)*complexghost_xiI(l,:)*(2.d0/r(l,:) &
                       - 0.5d0*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:))  &
                       + complexghost_xiI(l,:)*D1_alpha(l,:))/(A(l,:)*psi4(l,:))

  else

!    Term with second derivative of phi.

     !$OMP PARALLEL DO SCHEDULE(GUIDED)
     do i=1-ghost,Nrmax
        scomplexghost_piR(l,i) = alpha(l,i)/(A(l,i)*psi4(l,i))*D2_complexghost_phiR(l,i)
        scomplexghost_piI(l,i) = alpha(l,i)/(A(l,i)*psi4(l,i))*D2_complexghost_phiI(l,i)
     end do
     !$OMP END PARALLEL DO

!    Terms coming from Christoffel symbols.

     !$OMP PARALLEL DO SCHEDULE(GUIDED)
     do i=1-ghost,Nrmax
        scomplexghost_piR(l,i) = scomplexghost_piR(l,i) + (alpha(l,i)*D1_complexghost_phiR(l,i)*(2.d0/r(l,i) &
                          - 0.5d0*D1_A(l,i)/A(l,i) + D1_B(l,i)/B(l,i) + 2.d0*D1_phi(l,i))  &
                          + D1_complexghost_phiR(l,i)*D1_alpha(l,i))/(A(l,i)*psi4(l,i))
        scomplexghost_piI(l,i) = scomplexghost_piI(l,i) + (alpha(l,i)*D1_complexghost_phiI(l,i)*(2.d0/r(l,i) &
                          - 0.5d0*D1_A(l,i)/A(l,i) + D1_B(l,i)/B(l,i) + 2.d0*D1_phi(l,i))  &
                          + D1_complexghost_phiI(l,i)*D1_alpha(l,i))/(A(l,i)*psi4(l,i))
     end do
    !$OMP END PARALLEL DO

  end if

! Term proportional to trK.

  scomplexghost_piR(l,:) = scomplexghost_piR(l,:) + alpha(l,:)*trK(l,:)*complexghost_piR(l,:)
  scomplexghost_piI(l,:) = scomplexghost_piI(l,:) + alpha(l,:)*trK(l,:)*complexghost_piI(l,:)

! Potential term.

  if (complexghostpotential/="none") then
     scomplexghost_piR(l,:) = scomplexghost_piR(l,:) + alpha(l,:)*complexghost_VPR(l,:)
     scomplexghost_piI(l,:) = scomplexghost_piI(l,:) + alpha(l,:)*complexghost_VPI(l,:)
  end if

! Shift terms.

  if (shift/="none") then

!    Shift terms for phi.

     scomplexghost_phiR(l,:) = scomplexghost_phiR(l,:) + beta(l,:)*DA_complexghost_phiR(l,:)
     scomplexghost_phiI(l,:) = scomplexghost_phiI(l,:) + beta(l,:)*DA_complexghost_phiI(l,:)

!    Shift terms for xi.

     scomplexghost_xiR(l,:) = scomplexghost_xiR(l,:) + beta(l,:)*DA_complexghost_xiR(l,:) &
                            + complexghost_xiR(l,:)*D1_beta(l,:)
     scomplexghost_xiI(l,:) = scomplexghost_xiI(l,:) + beta(l,:)*DA_complexghost_xiI(l,:) &
                            + complexghost_xiI(l,:)*D1_beta(l,:)

!    Shift terms for pi.

     scomplexghost_piR(l,:) = scomplexghost_piR(l,:) + beta(l,:)*DA_complexghost_piR(l,:)
     scomplexghost_piI(l,:) = scomplexghost_piI(l,:) + beta(l,:)*DA_complexghost_piI(l,:)

  end if

! Dissipation.

  if (scalardiss/=0.d0) then

!    The symmetry depends on the value of complexghost_l.

     sym = int(1.d0-2.d0*mod(complexghost_l,2))

     dissipvar => complexghost_phiR
     sourcevar => scomplexghost_phiR
     call dissipation(l,sym,scalardiss)

     dissipvar => complexghost_phiI
     sourcevar => scomplexghost_phiI
     call dissipation(l,sym,scalardiss)

     dissipvar => complexghost_piR
     sourcevar => scomplexghost_piR
     call dissipation(l,sym,scalardiss)

     dissipvar => complexghost_piI
     sourcevar => scomplexghost_piI
     call dissipation(l,sym,scalardiss)

     dissipvar => complexghost_xiR
     sourcevar => scomplexghost_xiR
     call dissipation(l,-sym,scalardiss)

     dissipvar => complexghost_xiI
     sourcevar => scomplexghost_xiI
     call dissipation(l,-sym,scalardiss)

  end if

! Cosmological background quantities.

  if (cosmic_run) then

  end if


! **********************
! ***   BOUNDARIES   ***
! **********************

! For the boundary condition see the comment at the end
! of the routine "scalar.f90".  The condition here is
! the same.

  if ((boundtype/="static").and.(boundtype/="flat")) then

     scomplexghost_piR(l,Nr) = - (D1_complexghost_piR(l,Nr) + complexghost_piR(l,Nr)/r(l,Nr))
     scomplexghost_piI(l,Nr) = - (D1_complexghost_piI(l,Nr) + complexghost_piI(l,Nr)/r(l,Nr))

!    Potential term.

     if (complexghostpotential/="none") then
        scomplexghost_piR(l,Nr) = scomplexghost_piR(l,Nr) - 0.5d0*complexghost_VPR(l,Nr)
        scomplexghost_piI(l,Nr) = scomplexghost_piI(l,Nr) - 0.5d0*complexghost_VPI(l,Nr)
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_complexghost

