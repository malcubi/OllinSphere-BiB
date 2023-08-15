!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/scalar.f90,v 1.20 2023/02/20 19:48:12 malcubi Exp $

  subroutine sources_scalar(l)

! ************************************
! ***   SOURCES FOR SCALAR FIELD   ***
! ************************************

! This routine calculates the sources for a real
! scalar field.  The Klein-Gordon equation with
! an arbitrary potential V(phi) has the form:
!
! Box (phi)  -  (dV/dphi)  =  0
!
! Remember that we have defined the space and time
! derivatives of the scalar field as:
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
!                                                                \
!        +  xi ( 2/r  -  d A / 2A   +  d B / B  +  2 d ln(psi) ) |
!                         r             r             r          /
!                      4
!        +  ( 1 / A psi ) xi d alpha  +  alpha trK pi  -  alpha dV/dphi
!                             r
!
! If we want to evolve in second order form we just replace d_r(xi)
! in the last equation for the second derivative of phi.

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

! Source for scalar field.

  sscalar_phi(l,:) = alpha(l,:)*scalar_pi(l,:)

! Source for radial derivative.

  sscalar_xi(l,:) = alpha(l,:)*D1_scalar_pi(l,:) + scalar_pi(l,:)*D1_alpha(l,:)

! Source for time derivative. Notice that there are
! two ways of doing this: using only derivatives of
! phi, or using the xi.

  if (scalarmethod=="first") then

!    Term with derivative of xi.

     sscalar_pi(l,:) = alpha(l,:)/(A(l,:)*psi4(l,:))*D1_scalar_xi(l,:)

!    Terms coming from Christoffel symbols.

     sscalar_pi(l,:) = sscalar_pi(l,:) + (alpha(l,:)*scalar_xi(l,:)*(2.d0/r(l,:) &
                  - 0.5d0*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:)) &
                  + scalar_xi(l,:)*D1_alpha(l,:))/(A(l,:)*psi4(l,:))

  else

!    Term with second derivative of phi.

     !$OMP PARALLEL DO SCHEDULE(GUIDED)
     do i=1-ghost,Nrmax
        sscalar_pi(l,i) = alpha(l,i)/(A(l,i)*psi4(l,i))*D2_scalar_phi(l,i)
     end do
     !$OMP END PARALLEL DO

!    Terms coming from Christoffel symbols.

     !$OMP PARALLEL DO SCHEDULE(GUIDED)
     do i=1-ghost,Nrmax
        sscalar_pi(l,i) = sscalar_pi(l,i) + (alpha(l,i)*D1_scalar_phi(l,i)*(2.d0/r(l,i) &
                        - 0.5d0*D1_A(l,i)/A(l,i) + D1_B(l,i)/B(l,i) + 2.d0*D1_phi(l,i)) &
                        + D1_scalar_phi(l,i)*D1_alpha(l,i))/(A(l,i)*psi4(l,i))
     end do
     !$OMP END PARALLEL DO

  end if
                  
! Term proportional to trK.

  sscalar_pi(l,:) = sscalar_pi(l,:) + alpha(l,:)*trK(l,:)*scalar_pi(l,:)

! Potential term.

  if (scalarpotential/="none") then
     sscalar_pi(l,:) = sscalar_pi(l,:) - alpha(l,:)*scalar_VP(l,:)
  end if

! Shift terms.

  if (shift/="none") then

!    Shift terms for phi.

     sscalar_phi(l,:) = sscalar_phi(l,:) + beta(l,:)*DA_scalar_phi(l,:)

!    Shift terms for xi.

     sscalar_xi(l,:) = sscalar_xi(l,:) + beta(l,:)*DA_scalar_xi(l,:) &
                     + scalar_xi(l,:)*D1_beta(l,:)

!    Shift terms for pi.

     sscalar_pi(l,:) = sscalar_pi(l,:) + beta(l,:)*DA_scalar_pi(l,:)

  end if

! Dissipation.

  if (scalardiss/=0.d0) then

     dissipvar => scalar_phi
     sourcevar => sscalar_phi
     call dissipation(l,+1,scalardiss)

     dissipvar => scalar_pi
     sourcevar => sscalar_pi
     call dissipation(l,+1,scalardiss)

     dissipvar => scalar_xi
     sourcevar => sscalar_xi
     call dissipation(l,-1,scalardiss)

  end if

! Cosmological background quantities.

  if (cosmic_run) then

!    Source for phi.

     scosmobg_scalar_phi(l) = cosmobg_alpha(l)*cosmobg_scalar_pi(l)

!    Source for pi.

     scosmobg_scalar_pi(l) = cosmobg_alpha(l)*cosmobg_trK(l)*cosmobg_scalar_pi(l)

     if (scalarpotential/="none") then
        scosmobg_scalar_pi(l) = scosmobg_scalar_pi(l) - cosmobg_alpha(l)*cosmobg_scalar_VP(l)
     end if

  end if


! **********************
! ***   BOUNDARIES   ***
! **********************

! Radiative boundary for "Pi".
!
! In order to impose the boundary condition we assume that
! far away "Pi" behaves as an outgoing spherical wave:
!
! Pi  ~  g(r-t)/r
!
! This can be shown to imply:
!
! dPi/dt  ~  - (dPi/dr + Pi/r)                  (1)
!
! This condition is esentially equivalent to assuming that
! far away the "incoming" w^(-) field is much smaller that the
! "outgoing field" w^(+), where these fields are given by:
!
! w^(+)  =~  Pi - Xi
! w^(-)  =~  Pi + Xi
!
! which implies:
!
! Pi  =~  (w^(+) + w^(-))/2
! Xi  =~  (w^(-) - w^(+))/2
!
! Notice that the "incoming" field is not actually zero far
! away because of the term 2*Xi/r and the potential terms
! in the evolution equation for Pi, which couple the incoming
! and outgoing fields.
!
! In practice I have found that when there is no potential the
! incoming field is indeed much smaller than the outgoing one
! far away, so the boundary condition (1) works very well.
!
! However, when there is a non-zero potential we need to modify
! the boundary condition. It is not difficult to show that in
! the asymptotically flat region the incoming and outgoing fields
! evolve as:
!
! dw^(+-)/dt  +-  dw^(+-)/dr  =  [w^(-) - w^(+)]/r  -  VP
!
! with VP=dV/dphi.  If we again assume that the incoming field is
! much smaler than the outgoing field the above equations imply:
!
! dPi/dt  ~  - (dPi/dr + Pi/r + VP/2)           (2) 
!
! In practice I have found that, even though it is indeed true that
! the incoming field is smaller than the outgoing field far away,
! is is not "much smaller" unless the boundaries are really far.
! Condition (2) then works reasonably well, but it is not perfect.

  if ((boundtype/="static").and.(boundtype/="flat")) then

     if (.not.cosmic_run) then

        sscalar_pi(l,Nr) = - (D1_scalar_pi(l,Nr) + scalar_pi(l,Nr)/r(l,Nr))

!       Potential term.

        if (scalarpotential/="none") then
           sscalar_pi(l,Nr) = sscalar_pi(l,Nr) - 0.5d0*scalar_VP(l,Nr)
        end if

!    Cosmological runs need to include the lapse and metric coefficients
!    that are ignored in the asymptotically flat case.

     else

        sscalar_pi(l,Nr) = scosmobg_scalar_pi(l) - cosmobg_alpha(l)/cosmobg_a(l) &
                         *(D1_scalar_pi(l,Nr) + (scalar_pi(l,Nr)-cosmobg_scalar_pi(l))/r(l,Nr))

!       Potential term.

        if (scalarpotential/="none") then
           sscalar_pi(l,Nr) = sscalar_pi(l,Nr) - 0.5d0*cosmobg_alpha(l)*(scalar_VP(l,Nr)-cosmobg_scalar_VP(l))
        end if

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_scalar


