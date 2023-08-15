!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/ghost.f90,v 1.5 2023/02/20 19:48:12 malcubi Exp $

  subroutine sources_ghost(l)

! ***********************************
! ***   SOURCES FOR GHOST FIELD   ***
! ***********************************

! This routine calculates the sources for a real
! ghost field, which are identical to those of
! the standard scalar field.
!
! The Klein-Gordon equation with an arbitrary
! potential V(phi) has the form:
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
!        +  ( 1 / A psi ) xi d alpha  +  alpha trK pi  +  alpha dV/dphi
!                             r
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

  integer i,l


! *******************
! ***   SOURCES   ***
! *******************

! Source for ghost field.

  sghost_phi(l,:) = alpha(l,:)*ghost_pi(l,:)

! Source for radial derivative.

  sghost_xi(l,:) = alpha(l,:)*D1_ghost_pi(l,:) + ghost_pi(l,:)*D1_alpha(l,:)

! Source for time derivative. Notice that there are
! two ways of doing this: using only derivatives of
! phi, or using the xi.

  if (ghostmethod=="first") then

!    Term with derivative of xi.

     sghost_pi(l,:) = alpha(l,:)/(A(l,:)*psi4(l,:))*D1_ghost_xi(l,:)

!    Terms coming from Christoffel symbols.

     sghost_pi(l,:) = sghost_pi(l,:) + (alpha(l,:)*ghost_xi(l,:)*(2.d0/r(l,:) &
                  - 0.5d0*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 2.d0*D1_phi(l,:)) &
                  + ghost_xi(l,:)*D1_alpha(l,:))/(A(l,:)*psi4(l,:))
  else

!    Term with second derivative of phi.

     !$OMP PARALLEL DO SCHEDULE(GUIDED)
     do i=1-ghost,Nrmax
        sghost_pi(l,i) = alpha(l,i)/(A(l,i)*psi4(l,i))*D2_ghost_phi(l,i)
     end do
     !$OMP END PARALLEL DO

!    Terms coming from Christoffel symbols.

     !$OMP PARALLEL DO SCHEDULE(GUIDED)
     do i=1-ghost,Nrmax
        sghost_pi(l,i) = sghost_pi(l,i) + (alpha(l,i)*D1_ghost_phi(l,i)*(2.d0/r(l,i) &
                       - 0.5d0*D1_A(l,i)/A(l,i) + D1_B(l,i)/B(l,i) + 2.d0*D1_phi(l,i)) &
                       + D1_ghost_phi(l,i)*D1_alpha(l,i))/(A(l,i)*psi4(l,i))
     end do
     !$OMP END PARALLEL DO

  end if
                  
! Term proportional to trK.

  sghost_pi(l,:) = sghost_pi(l,:) + alpha(l,:)*trK(l,:)*ghost_pi(l,:)

! Potential term.

  if (ghostpotential/="none") then
     sghost_pi(l,:) = sghost_pi(l,:) + alpha(l,:)*ghost_VP(l,:)
  end if

! Shift terms.

  if (shift/="none") then

!    Shift terms for phi.

     sghost_phi(l,:) = sghost_phi(l,:) + beta(l,:)*DA_ghost_phi(l,:)

!    Shift terms for xi.

     sghost_xi(l,:) = sghost_xi(l,:) + beta(l,:)*DA_ghost_xi(l,:) &
                    + ghost_xi(l,:)*D1_beta(l,:)

!    Shift terms for pi.

     sghost_pi(l,:) = sghost_pi(l,:) + beta(l,:)*DA_ghost_pi(l,:)

  end if

! Dissipation.

  if (scalardiss/=0.d0) then

     dissipvar => ghost_phi
     sourcevar => sghost_phi
     call dissipation(l,+1,scalardiss)

     dissipvar => ghost_pi
     sourcevar => sghost_pi
     call dissipation(l,+1,scalardiss)

     dissipvar => ghost_xi
     sourcevar => sghost_xi
     call dissipation(l,-1,scalardiss)

  end if

! Cosmological background quantities.

  if (cosmic_run) then
     if (rank==0) then
        print *, 'Cosmic runs not yet implemented for ghost field.'
        print *, 'Aborting! (subroutine ghost.f90)'
        print *
     end if
     call die
  end if


! **********************
! ***   BOUNDARIES   ***
! **********************

! For the boundary condition see the comment at the end
! of the routine "scalar.f90".  The condition here is
! the same (but do remember that for a ghost field the
! sign of the potential is reversed).

  if ((boundtype/="static").and.(boundtype/="flat")) then

     if (.not.cosmic_run) then

        sghost_pi(l,Nr) = - (D1_ghost_pi(l,Nr) + ghost_pi(l,Nr)/r(l,Nr))

!       Potential term.

        if (ghostpotential/="none") then
           sghost_pi(l,Nr) = sghost_pi(l,Nr) + 0.5d0*ghost_VP(l,Nr)
        end if

!    Cosmological runs need to include the lapse and metric coefficients
!    that are ignored in the asymptotically flat case.

     else

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_ghost

