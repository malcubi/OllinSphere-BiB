
  subroutine sources_complexproca(l)

! *******************************************
! ***   SOURCES FOR COMPLEX PROCA FIELD   ***
! *******************************************

! This routine calculates the sources for the
! different quantities associated to the
! complex Proca field.

! Include modules.

  use param
  use arrays
  use derivatives

! Extra variables.

  implicit none

  integer l,i

  real(8) half,one
  real(8) aux


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0
  one  = 1.d0


! *******************************
! ***   SOURCES FOR PROCA E   ***
! *******************************

! The evolution equation for the Proca electric field E (index up) is:
!
!
! d procaE  =  beta d procaE  -  procaE d beta
!  t                 r                   r
!
!                                      2              4
!           +  alpha ( trK procaE  +  m procaA / A psi  )
!
! This equation comes from:
!
! __   mu nu     2  nu
! \/  W      =  m  X
!  mu
!
! These are essentially Maxwell's equations in vacuum but with a mass term.

  scprocaE_R(l,:) = alpha(l,:)*(trK(l,:)*cprocaE_R(l,:) &
                  + cproca_mass**2*cprocaA_R(l,:)/(A(l,:)*psi4(l,:)))
  scprocaE_I(l,:) = alpha(l,:)*(trK(l,:)*cprocaE_I(l,:) &
                  + cproca_mass**2*cprocaA_I(l,:)/(A(l,:)*psi4(l,:)))

! Shift terms (Lie derivative).

  if (shift/="none") then
     scprocaE_R(l,:) = scprocaE_R(l,:) + beta(l,:)*DA_cprocaE_R(l,:) - cprocaE_R(l,:)*D1_beta(l,:)
     scprocaE_I(l,:) = scprocaE_I(l,:) + beta(l,:)*DA_cprocaE_I(l,:) - cprocaE_I(l,:)*D1_beta(l,:)
  end if

! Charge terms if needed. In the case of a charged Proca field we need
! to add the following term to the source:
!
! scprocaE_R  =  scprocaE_R - q alpha ePhi cprocaE_I
!
! scprocaE_I  =  scprocaE_R + q alpha ePhi cprocaE_R

  if (cproca_q/=0.d0) then
     scprocaE_R(l,:) = scprocaE_R(l,:) - cproca_q*alpha(l,:)*ePhi(l,:)*cprocaE_I(l,:)
     scprocaE_I(l,:) = scprocaE_I(l,:) + cproca_q*alpha(l,:)*ePhi(l,:)*cprocaE_R(l,:)
  end if

! Dissipation.

  if (procadiss/=0.d0) then

     dissipvar => cprocaE_R
     sourcevar => scprocaE_R
     call dissipation(l,-1,procadiss)

     dissipvar => cprocaE_I
     sourcevar => scprocaE_I
     call dissipation(l,-1,procadiss)

  end if


! ****************************************
! ***   SOURCES FOR SCALAR POTENTIAL   ***
! ****************************************

! The evolution equation for the Proca scalar potential Phi
! is given by the Lorenz gauge equation, which in the
! case of a Proca field is not a choice but must always
! be satisfied. It has the form:
!
!                                                           4
! d procaPhi  =  beta d procaPhi  -  d (alpha procaA) / A psi
!  t                   r              r
!
!                               4  /                                    \
!        +  alpha procaA / A psi   | d A / 2A - d B / B - 2 d phi - 2/r |
!                                  \  r          r           r          /
!
!        +  alpha procaPhi trK
!
!
! (remember that procaF = alpha*procaPhi and procaH = alpha*procaA).

  scprocaPhi_R(l,:) = one/(A(l,:)*psi4(l,:))*(-D1_cprocaH_R(l,:) &
                    + cprocaH_R(l,:)*(half*D1_A(l,:)/A(l,:) - D1_B(l,:)/B(l,:) &
                    - 2.d0*D1_phi(l,:) - 2.d0/r(l,:))) + trK(l,:)*cprocaF_R(l,:)

  scprocaPhi_I(l,:) = one/(A(l,:)*psi4(l,:))*(-D1_cprocaH_I(l,:) &
                    + cprocaH_I(l,:)*(half*D1_A(l,:)/A(l,:) - D1_B(l,:)/B(l,:) &
                    - 2.d0*D1_phi(l,:) - 2.d0/r(l,:))) + trK(l,:)*cprocaF_I(l,:)

! Shift terms (Lie derivative).

  if (shift/="none") then
     scprocaPhi_R(l,:) = scprocaPhi_R(l,:) + beta(l,:)*DA_cprocaPhi_R(l,:)
     scprocaPhi_I(l,:) = scprocaPhi_I(l,:) + beta(l,:)*DA_cprocaPhi_I(l,:)
  end if

! Charge terms if needed.  In the case of a charged Proca field we need
! to add the following term to the source:
!
! scprocaPhi_R  =  scprocaPhi_R - q alpha [ ePhi cprocaPhi_I - eAr cprocaA_I /(A psi4)
!               + (A psi4 / cproca_mass**2) electric cprocaE_I ]
!
! scprocaPhi_I  =  scprocaPhi_R + q alpha [ ePhi cprocaPhi_R - eAr cprocaA_R /(A psi4)
!               + (A psi4 / cproca_mass**2) electric cprocaE_R ]

  if (cproca_q/=0.d0) then
     scprocaPhi_R(l,:) = scprocaPhi_R(l,:) - cproca_q*alpha(l,:) &
        *(ePhi(l,:)*cprocaPhi_I(l,:) - eAr(l,:)*cprocaA_I(l,:)/(A(l,:)*psi4(l,:)) &
        + A(l,:)*psi4(l,:)/cproca_mass**2*electric(l,:)*cprocaE_I(l,:))
     scprocaPhi_I(l,:) = scprocaPhi_I(l,:) + cproca_q*alpha(l,:) &
        *(ePhi(l,:)*cprocaPhi_R(l,:) - eAr(l,:)*cprocaA_R(l,:)/(A(l,:)*psi4(l,:)) &
        + A(l,:)*psi4(l,:)/cproca_mass**2*electric(l,:)*cprocaE_R(l,:))
  end if

! Dissipation.

  if (procadiss/=0.d0) then

     dissipvar => cprocaPhi_R
     sourcevar => scprocaPhi_R
     call dissipation(l,+1,procadiss)

     dissipvar => cprocaPhi_I
     sourcevar => scprocaPhi_I
     call dissipation(l,+1,procadiss)

  end if


! ****************************************
! ***   SOURCES FOR VECTOR POTENTIAL   ***
! ****************************************

! The evolution equation for the Proca vector potential A (index down)
! comes directly from the definition of the electric field, and takes
! the form:
!
! d procaA  =   beta d procaA  +  procaA d beta
!  t                  r                   r
!
!                         4
!           -  alpha A psi procaE  -  d ( alpha procaPhi )
!                                      r
!
! (remember that procaF = alpha*procaPhi)

  scprocaA_R(l,:) = - alpha(l,:)*A(l,:)*psi4(l,:)*cprocaE_R(l,:) - D1_cprocaF_R(l,:)
  scprocaA_I(l,:) = - alpha(l,:)*A(l,:)*psi4(l,:)*cprocaE_I(l,:) - D1_cprocAF_I(l,:)

! Shift terms (Lie derivative).

  if (shift/="none") then
     scprocaA_R(l,:) = scprocaA_R(l,:) + beta(l,:)*DA_cprocaA_R(l,:) + cprocaA_R(l,:)*D1_beta(l,:)
     scprocaA_I(l,:) = scprocaA_I(l,:) + beta(l,:)*DA_cprocaA_I(l,:) + cprocaA_I(l,:)*D1_beta(l,:)
  end if

! Charge terms if needed. In the case of a charged Proca field we need
! to add the following term to the source:
!
! scprocaA_R  =  scprocaA_R - q alpha ( ePhi cprocaA_I - eAr cprocaPhi_I )
!
! scprocaA_I  =  scprocaA_I + q alpha ( ePhi cprocaA_R - eAr cprocaPhi_R )

  if (cproca_q/=0.d0) then
     scprocaA_R(l,:) = scprocaA_R(l,:) - cproca_q*alpha(l,:) &
                     *(ePhi(l,:)*cprocaA_I(l,:) - eAr(l,:)*cprocaPhi_I(l,:))
     scprocaA_I(l,:) = scprocaA_I(l,:) + cproca_q*alpha(l,:) &
                     *(ePhi(l,:)*cprocaA_R(l,:) - eAr(l,:)*cprocaPhi_R(l,:))
  end if

! Dissipation.

  if (procadiss/=0.d0) then

     dissipvar => cprocaA_R
     sourcevar => scprocaA_R
     call dissipation(l,-1,procadiss)

     dissipvar => cprocaA_I
     sourcevar => scprocaA_I
     call dissipation(l,-1,procadiss)

  end if


! **********************
! ***   BOUNDARIES   ***
! **********************

! The electric field evolves only through sources (since there
! is no magnetic field), so it needs no boundary condition.
! For the potentials we do need a boundary condition.

! In order to impose the boundary condition we assume that far
! away procaPhi and procaA behave as outgoing spherical waves:
!
! u(r,t)  ~  g(r-vt)/r
!
! This can be shown to imply:
!
! du/dt   ~  - v (du/dr + u/r)
!
! where v is the coordinate speed of light:  v = alpha / (sqrt(A)*psi**2)
!
! I have experimented with radiative boundary conditions for both
! procaPhi and ProcaA. It turns out that the best well behaved case
! is to aply the radiative condition only to procaPhi, and leave procaA
! to evolve freely all the way to the boundary. This seems to be stable
! and the Gauss constraint converges to zero.
!
! This is the opposite from what happened with the Maxwell field,
! where the best behaved case was applying the radiative condition
! tio the vector potential.  This might be because the Maxwell field
! is long range while the Proca field is not (it is massive).
!
! The extra second term when applying the radiative condition
! to A improves it, and makes when used makes the radiative
! condition for A almost identical to the radiative condition
! for Phi.  But it seems that it is easier to just apply the
! condition to Phi.  I leave it there just in case it works
! better in some cases.

  if ((boundtype/="static").and.(boundtype/="flat")) then

     aux = alpha(l,Nr)/sqrt(A(l,Nr))/psi2(l,Nr)

!    Boundary condition for procaPhi.

     if (procabound=="radPhi") then

        scprocaPhi_R(l,Nr) = - aux*(D1_cprocaPhi_R(l,Nr) + cprocaPhi_R(l,Nr)/r(l,Nr))
        scprocaPhi_I(l,Nr) = - aux*(D1_cprocaPhi_I(l,Nr) + cprocaPhi_I(l,Nr)/r(l,Nr))

!    Boundary condition for procaA.

     else if (procabound=="radA") then

        scprocaA_R(l,Nr) = - aux*(D1_cprocaA_R(l,Nr) + cprocaA_R(l,Nr)/r(l,Nr))
        scprocaA_I(l,Nr) = - aux*(D1_cprocaA_I(l,Nr) + cprocaA_I(l,Nr)/r(l,Nr))

        scprocaA_R(l,Nr) = scprocaA_R(l,Nr) - alpha(l,Nr)*A(l,Nr)*psi4(l,Nr)*cprocaE_R(l,Nr) &
                         + (alpha(l,Nr)*cprocaPhi_R(l,Nr) - aux*cprocaA_R(l,Nr))/r(l,Nr)
        scprocaA_I(l,Nr) = scprocaA_I(l,Nr) - alpha(l,Nr)*A(l,Nr)*psi4(l,Nr)*cprocaE_I(l,Nr) &
                         + (alpha(l,Nr)*cprocaPhi_I(l,Nr) - aux*cprocaA_I(l,Nr))/r(l,Nr)

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_complexproca



