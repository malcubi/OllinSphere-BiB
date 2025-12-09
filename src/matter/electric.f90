
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

! In order to impose the boundary condition we assume that far
! away ePhi and eAr behave as outgoing spherical waves:
!
! u(r,t)  ~  g(r-vt)/r
!
! This can be shown to imply:
!
! du/dt  ~  - v (du/dr + u/r)
!
! where v is the coordinate speed of light:  v = alpha / (sqrt(A)*psi**2)
!
! I have experimented with radiative boundary conditions for both
! ePhi and eAr.  It turns out that the best well behaved case
! is to apply the radiative condition only to procaA, and leave procaPhi
! to evolve freely all the way to the boundary. This seems to be stable
! and the Gauss constraint converges to zero.
!
! The main reason why applying the condition to A works better in this
! case is that for stationary charged configurations Phi is typically
! non-trivial (it is the electrostatic potential), while A is zero.

  if ((boundtype/="static").and.(boundtype/="flat")) then

     aux = alpha(l,Nr)/sqrt(A(l,Nr))/psi2(l,Nr)

!    Boundary conditions for ePhi.

     !sePhi(l,Nr) = - aux*(D1_ePhi(l,Nr) + ePhi(l,Nr)/r(l,Nr))

!    Boundary condition for eAr.

     seAr(l,Nr) = - aux*(D1_eAr(l,Nr) + eAr(l,Nr)/r(l,Nr))

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_electric

