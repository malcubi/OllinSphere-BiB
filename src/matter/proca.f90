!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/proca.f90,v 1.21 2023/02/22 21:29:19 malcubi Exp $

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

! The electric field evolves only through sources (since there
! is no magnetic field), so it needs no boundary condition.
! For the potentials we do need a boundary condition.

! Boundary condition for procaPhi: We assume a simple outgoing
! spherical wave.  We also add a contribution to the source (-E/2).
! This term comes from expressing ePhi in terms of incoming and
! outgoing eigenfields, and then assuming the incoming eigenfield
! is very small.

  sprocaPhi(l,Nr) = - (D1_procaPhi(l,Nr) + procaPhi(l,Nr)/r(l,Nr)) - 0.5d0*procaE(l,Nr)

! The variable procaA in fact does not need a boundary condition,
! but ... it turns out that for this type of first order system
! the higher order one-sided derivatives cause instabilities,
! so here I recalculate the sources close to the boundary to
! fourth order.

  do i=Nr-ghost+2,Nr-1

     aux = (3.d0*procaH(l,i+1) + 10.d0*procaH(l,i) - 18.d0*procaH(l,i-1) &
         + 6.d0*procaH(l,i-2) - procaH(l,i-3))/12.d0/dr(l)

     sprocaPhi(l,i) = one/(A(l,i)*psi4(l,i))*(procaH(l,i)*(half*D1_A(l,i)/A(l,i) - D1_B(l,i)/B(l,i) &
                    - 2.d0*D1_phi(l,i) - 2.d0/r(l,i)) - aux) + trK(l,i)*procaF(l,i)

     aux = (3.d0*procaF(l,i+1) + 10.d0*procaF(l,i) - 18.d0*procaF(l,i-1) &
         + 6.d0*procaF(l,i-2) - procaF(l,i-3))/12.d0/dr(l)

     sprocaA(l,i) = - alpha(l,i)*A(l,i)*psi4(l,i)*procaE(l,i) - aux

     if (shift/="none") then
        sprocaPhi(l,i) = sprocaPhi(l,i) + beta(l,i)*DA_procaPhi(l,i)
        sprocaA(l,i) = sprocaA(l,i) + beta(l,i)*DA_procaA(l,i) + procaA(l,i)*D1_beta(l,i)
     end if

  end do

! Source for procaA at last point.

  sprocaA(l,Nr) = - alpha(l,Nr)*A(l,Nr)*psi4(l,Nr)*procaE(l,Nr) - D1_procaF(l,Nr)

  if (shift/="none") then
     sprocaA(l,Nr) = sprocaA(l,Nr) + beta(l,Nr)*DA_procaA(l,Nr) + procaA(l,Nr)*D1_beta(l,Nr)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_proca

