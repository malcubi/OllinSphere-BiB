!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/nonmin.f90,v 1.4 2022/01/26 19:01:56 erik Exp $

  subroutine sources_nonmin(l)

! ************************************
! ***   SOURCES FOR NONMIN FIELD   ***
! ************************************

! This routine calculates the sources for a non-minimally
! coupled scalar field.  Notice that this field obeys an
! evolution equation of the form:
!
! Box(phi)  =  RHS
!
! where RHS is a function of the nonmininal coupling f(phi)
! that has already been calculated in "stressenergy.f90" and
! has the form:
!
!         /             2        \ -1
! RHS  =  | f ( 1 + 3 fp / 2 f ) |
!         \                      /
!
!         /                                        2        4      2
!         | f VP - 2 fp V - (fp/2) (1 + 3 fpp) ( xi / (A psi ) - pi )
!         \
!
!                        \
!         + (fp/2) T     |
!                   matt /
!
! with T_matt the trace of the stress energy tensor of all
! matter fields NOT including the non-minimally coupled field
! itself (see the paper by M. Salgado: Class. Quant. Grav. 23,
! 4719-4741, 2006).
!
! The nonmin field is evolved in first order form.
! Remember that we have defined the space and time
! derivative arrays as:
!        
! pi  =  ( d phi  -  beta d phi ) / alpha
!           t              r
!
! xi  =  d phi
!         r
!
! The evolution equations for {phi,xi,pi} are:
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
!                                                             \
!        +  xi ( 2/r -  d A / 2A   +  d B / B  +  2 d ln(psi) |
!                        r             r             r        /
!                      4
!        +  ( 1 / A psi ) xi d alpha  +  alpha trK pi  -  alpha RHS
!                             r

! Include modules.
  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer i,l


! *******************
! ***   SOURCES   ***
! *******************

! Source for nonmin field.

  snonmin_phi(l,:) = alpha(l,:)*nonmin_pi(l,:)

! Source for radial derivative.

  snonmin_xi(l,:) = alpha(l,:)*D1_nonmin_pi(l,:) + nonmin_pi(l,:)*D1_alpha(l,:)

! Source for time derivative. Notice that there are
! two ways of doing this: using only derivatives of
! phi, or using the xi.

  if (nonminmethod=="first") then

!    Term with derivative of xi.

     snonmin_pi(l,:) = alpha(l,:)*D1_nonmin_xi(l,:)/(A(l,:)*psi4(l,:))

!    Terms coming from Christoffel symbols.

     snonmin_pi(l,:) = snonmin_pi(l,:) + (alpha(l,:)*nonmin_xi(l,:)*(2.d0/r(l,:) &
             - 0.5d0*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:))  &
             + nonmin_xi(l,:)*D1_alpha(l,:))/(A(l,:)*psi4(l,:))

  else
  !$OMP PARALLEL DO SCHEDULE(GUIDED)
  do i=1-ghost,Nrmax
!    Term with second derivative of phi.

     snonmin_pi(l,i) = alpha(l,i)*D2_nonmin_phi(l,i)/(A(l,i)*psi4(l,i))

!    Terms coming from Christoffel symbols.

     snonmin_pi(l,i) = snonmin_pi(l,i) + (alpha(l,i)*D1_nonmin_phi(l,i)*(2.d0/r(l,i) &
             - 0.5d0*D1_A(l,i)/A(l,i) + D1_B(l,i)/B(l,i) + 2.d0*D1_phi(l,i))  &
             + D1_nonmin_phi(l,i)*D1_alpha(l,i))/(A(l,i)*psi4(l,i))
  end do
  !$OMP END PARALLEL DO
  end if

! Terms proportional to trK and nonmin_rhs.

  snonmin_pi(l,:) = snonmin_pi(l,:) + alpha(l,:)*(trK(l,:)*nonmin_pi(l,:) - nonmin_rhs(l,:))

! Potential term.

  if (nonminpotential/="none") then
     snonmin_pi(l,:) = snonmin_pi(l,:) - alpha(l,:)*nonmin_VP(l,:)
  end if

! Shift terms.

  if (shift/="none") then

!    Shift terms for phi.

     snonmin_phi(l,:) = snonmin_phi(l,:) + beta(l,:)*DA_nonmin_phi(l,:)

!    Shift terms for xi.

     snonmin_xi(l,:) = snonmin_xi(l,:) + beta(l,:)*DA_nonmin_xi(l,:) + nonmin_xi(l,:)*D1_beta(l,:)

!    Shift terms for pi.

     snonmin_pi(l,:) = snonmin_pi(l,:) + beta(l,:)*DA_nonmin_pi(l,:)

  end if

! Dissipation.

  if (nonmindiss/=0.d0) then

     dissipvar => nonmin_phi
     sourcevar => snonmin_phi
     call dissipation(l,+1,nonmindiss)

     dissipvar => nonmin_pi
     sourcevar => snonmin_pi
     call dissipation(l,+1,nonmindiss)

     dissipvar => nonmin_xi
     sourcevar => snonmin_xi
     call dissipation(l,-1,nonmindiss)

  end if


! **********************
! ***   BOUNDARIES   ***
! **********************

! The incoming eigenfield is:
!
! wm  =  nonmin_pi + nonmin_xi / (sqrt(A) psi**2)
!
! In order to impose the boundary condition we assume
! that at the boundary the scalar field behaves as:
!
! nonmin_phi = f(r-vt)/r
!
! with v = alpha sqrt(1/A/psi**2) the speed of light.
! This can be shown to imply:
!
! wm = - nonmin_phi / r / (sqrt(A) psi**2)
!
! Solving now for nonmin_pi we find:
!
! nonmin_pi  =  - ( nonmin_phi/r + nonmin_xi ) / (sqrt(A) psi**2)
!
! Taking the time derivative and ignoring terms that are
! quadratic in small quantites we find:
!
! snonmin_pi ~ - ( snonmin_phi/r + snonmin_xi )

  snonmin_pi(l,Nr) = - (snonmin_phi(l,Nr)/r(l,Nr) + snonmin_xi(l,Nr))


! ***************
! ***   END   ***
! ***************

  end subroutine sources_nonmin

