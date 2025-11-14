
  subroutine sources_electric(l)

! **************************************
! ***   SOURCES FOR ELECTRIC FIELD   ***
! **************************************

! This routine calculates the sources for the electric
! field and the potentials ePhi and eAr.

! Include modules.

  use param
  use arrays
  use derivatives

! Extra variables.

  implicit none

  integer l,i,imax

  real(8) half,one,smallpi
  real(8) aux


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0
  one  = 1.d0

  smallpi = acos(-1.d0)


! **************************************
! ***   SOURCES FOR ELECTRIC FIELD   ***
! **************************************

! The evolution equation for the electric field is:
!
!
! d electric  =  beta d electric  -  electric d beta
!  t                   r                       r
!
!             +  alpha trK electric  -  4 pi alpha ecurrent
!
!
! This evolution equation is quite simple because
! in spherical symmetry we have no magnetic field.

! Term with trK.

  selectric(l,:) = alpha(l,:)*trK(l,:)*electric(l,:)

! Electric current term.

  selectric(l,:) = selectric(l,:) - 4.d0*smallpi*alpha(l,:)*ecurrent(l,:)

! Shift terms (Lie derivative).

  if (shift/="none") then
     selectric(l,:) = selectric(l,:) + beta(l,:)*DA_electric(l,:) - electric(l,:)*D1_beta(l,:)
  end if

! Dissipation.

  if (elecdiss/=0.d0) then
     dissipvar => electric
     sourcevar => selectric
     call dissipation(l,-1,elecdiss)
  end if


! **********************************
! ***   SOURCES FOR POTENTIALS   ***
! **********************************

! Notice that the evolution of the potentials does not feed back
! into the evolution for the electric field.  Still, the potentials
! are needed for coupling to different matter models (e.g. a complex
! scalar field).
!
! In the Lorenz gauge the evolution equations for the potentials are:
!
!                                                 4
! d ePhi  =  beta d ePhi  -  d (alpha eAr) / A psi
!  t               r          r
!                                              
!                             4  /                                    \
!         +  alpha eAr / A psi   | d A / 2A - d B / B - 2 d phi - 2/r |
!                                \  r          r           r          /
!
!
!         +  alpha ePhi trK
!
!
!
! d eAr   =  beta d eAr  +  eAr d beta
!  t               r             r
!
!                       4
!         -  alpha A psi  electric  -  d ( alpha ePhi )
!                                       r
!
! Notice that eAr has the index down in the above expressions, while the
! electric field has index up, so we need to be careful with this
! (that is why we have the A*psi**4 in the term with alpha*electric).

! Source for ePhi (remember that eF = alpha*ePhi and eH = alpha*eAr).

  sePhi(l,:) = one/(A(l,:)*psi4(l,:))*(-D1_eH(l,:) &
             + eH(l,:)*(half*D1_A(l,:)/A(l,:) - D1_B(l,:)/B(l,:) &
             - 2.d0*D1_phi(l,:) - 2.d0/r(l,:))) + trK(l,:)*eF(l,:)

  if (shift/="none") then
     sePhi(l,:) = sePhi(l,:) + beta(l,:)*DA_ePhi(l,:)
  end if

! Source of eAr.

  seAr(l,:) = - alpha(l,:)*A(l,:)*psi4(l,:)*electric(l,:) - D1_eF(l,:)

  if (shift/="none") then
     seAr(l,:) = seAr(l,:) + beta(l,:)*DA_eAr(l,:) + eAr(l,:)*D1_beta(l,:)
  end if

! Dissipation.

  if (elecdiss/=0.d0) then

     dissipvar => ePhi
     sourcevar => sePhi
     call dissipation(l,+1,elecdiss)

     dissipvar => eAr
     sourcevar => seAr
     call dissipation(l,-1,elecdiss)

  end if


! **********************
! ***   BOUNDARIES   ***
! **********************

! The electric field evolves only through sources (since there
! is no magnetic field), so it needs no boundary condition.
! For the potentials we do need a boundary condition.

! Boundary conditions for ePhi:  We assume a simple outgoing
! spherical wave. Notice that for a localized charge
! distribution that does not extend all the way to the boundary
! the right hand side below should vanish since we expect
! to have sPhi ~ 1/r (of course this is only exactly true
! for Minkowski).

  sePhi(l,Nr) = - (D1_ePhi(l,Nr) + ePhi(l,Nr)/r(l,Nr))

! The variable eAr in fact does not need a boundary condition,
! but ... it turns out that for this type of first order system
! the higher order one-sided derivatives cause instabilities,
! so here I recalculate the sources close to the boundary to
! fourth order.

  if ((order/="two").or.(order/="four")) then

     do i=Nr-ghost+2,Nr-1

        aux = (3.d0*eH(l,i+1) + 10.d0*eH(l,i) - 18.d0*eH(l,i-1) + 6.d0*eH(l,i-2) - eH(l,i-3))/12.d0/dr(l)

        sePhi(l,i) = one/(A(l,i)*psi4(l,i))*(eH(l,i)*(half*D1_A(l,i)/A(l,i) - D1_B(l,i)/B(l,i) &
                   - 2.d0*D1_phi(l,i) - 2.d0/r(l,i)) - aux) + trK(l,i)*eF(l,i)

        aux = (3.d0*eF(l,i+1) + 10.d0*eF(l,i) - 18.d0*eF(l,i-1) + 6.d0*eF(l,i-2) - eF(l,i-3))/12.d0/dr(l)

        seAr(l,i) = - alpha(l,i)*A(l,i)*psi4(l,i)*electric(l,i) - aux

        if (shift/="none") then
           sePhi(l,i) = sePhi(l,i) + beta(l,i)*DA_ePhi(l,i)
           seAr(l,i)  = seAr(l,i) + beta(l,i)*DA_eAr(l,i) + eAr(l,i)*D1_beta(l,i)
        end if

     end do

!    Source for eAr at last point.

     seAr(l,Nr) = - alpha(l,Nr)*A(l,Nr)*psi4(l,Nr)*electric(l,Nr) - D1_eF(l,Nr)

     if (shift/="none") then
        seAr(l,Nr) = seAr(l,Nr) + beta(l,Nr)*DA_eAr(l,Nr) + eAr(l,Nr)*D1_beta(l,Nr)
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_electric

