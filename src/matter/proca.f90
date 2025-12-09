
  subroutine sources_proca(l)

! ***********************************
! ***   SOURCES FOR PROCA FIELD   ***
! ***********************************

! This routine calculates the sources for the
! different Proca field quantities.

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

  sprocaE(l,:) = alpha(l,:)*(trK(l,:)*procaE(l,:) &
               + proca_mass**2*procaA(l,:)/(A(l,:)*psi4(l,:)))

! Shift terms (Lie derivative).

  if (shift/="none") then
     sprocaE(l,:) = sprocaE(l,:) + beta(l,:)*DA_procaE(l,:) - procaE(l,:)*D1_beta(l,:)
  end if

! Dissipation.

  if (procadiss/=0.d0) then
     dissipvar => procaE
     sourcevar => sprocaE
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
!                                                            4
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

  sprocaPhi(l,:) = one/(A(l,:)*psi4(l,:))*(-D1_procaH(l,:)  &
                 + procaH(l,:)*(half*D1_A(l,:)/A(l,:) - D1_B(l,:)/B(l,:) &
                 - 2.d0*D1_phi(l,:) - 2.d0/r(l,:))) + trK(l,:)*procaF(l,:)

! Shift terms (Lie derivative).

  if (shift/="none") then
     sprocaPhi(l,:) = sprocaPhi(l,:) + beta(l,:)*DA_procaPhi(l,:)
  end if

! Dissipation.

  if (procadiss/=0.d0) then
     dissipvar => procaPhi
     sourcevar => sprocaPhi
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

  sprocaA(l,:) = - alpha(l,:)*A(l,:)*psi4(l,:)*procaE(l,:) - D1_procaF(l,:)

! Shift terms (Lie derivative).

  if (shift/="none") then
     sprocaA(l,:) = sprocaA(l,:) + beta(l,:)*DA_procaA(l,:) + procaA(l,:)*D1_beta(l,:)
  end if

! Dissipation.

  if (procadiss/=0.d0) then
     dissipvar => procaA
     sourcevar => sprocaA
     call dissipation(l,-1,procadiss)
  end if


! **********************
! ***   BOUNDARIES   ***
! **********************

! The proca electric field evolves only through sources (since
! there is no magnetic field), so it needs no boundary condition.
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
! procaPhi and ProcaA.  Applying the radiative condition to Phi
! works reasonably well, but applying it to A directly introduces
! larger errors.  Applying it to both at the same time is bad,
! it introduces constraint violations (which makes sense since
! the system becomes overdetermined).
!
! One can improve the condition for A by adding the source term
! coming from the electric field. When doing this the condition
! for A seems to work better that the condition for Phi.

  if ((boundtype/="static").and.(boundtype/="flat")) then

     aux = alpha(l,Nr)/sqrt(A(l,Nr))/psi2(l,Nr)

!    Boundary condition for procaPhi.

     if (procabound=="radPhi") then

        sprocaPhi(l,Nr) = - aux*(D1_procaPhi(l,Nr) + procaPhi(l,Nr)/r(l,Nr))

!    Boundary condition for procaA.

     else if (procabound=="radA") then

        sprocaA(l,Nr) = - aux*(D1_procaA(l,Nr) + procaA(l,Nr)/r(l,Nr)) &
                      - alpha(l,Nr)*A(l,Nr)*psi4(l,Nr)*procaE(l,Nr)

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_proca

