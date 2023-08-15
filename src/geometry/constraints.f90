!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/constraints.f90,v 1.13 2022/08/03 19:39:17 malcubi Exp $

  subroutine constraints

! *************************************
! ***   EVALUATION OF CONSTRAINTS   ***
! *************************************

! This routine evaluates the hamiltonian, momentum and
! other constraints.
!
! Remember to check if the arrays have been allocated,
! since the constraints are declared as AUXILIARY and will
! only have memory allocated if one wants output for them.

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  logical contains

  real(8) third,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  third = 1.d0/3.d0

  smallpi = acos(-1.d0)


! ***********************
! ***   HAMILTONIAN   ***
! ***********************

! The hamiltonian constraint in this case takes the form:
!
!                    2        2               2
! ham  =  R  -  ( KTA  + 2 KTB  )  +  2/3  trK  -  16 pi rho
!
! Where R is the physical Ricci scalar.

  if (allocated(ham)) then

     ham = RSCAL - (KTA**2 + 2.d0*KTB**2) + 2.d0*third*trK**2

     if (mattertype /= "vacuum") then
        if (contains(mattertype,"nonmin")) then
           ham = ham - 2.d0/nonmin_f*rho
        else
           ham = ham - 16.d0*smallpi*rho
        end if
     end if

     !ham = psi**5*ham

!    Boundary. The calculation above is done all the way to
!    the boundary with one-sided differences, but this introduces
!    large artificial errors, so I fix it here.
  
     if (rank==size-1) then
        ham(:,Nr  ) = ham(:,Nr-3)
        ham(:,Nr-1) = ham(:,Nr-3)
        ham(:,Nr-2) = ham(:,Nr-3)
     end if

  end if

! Absolute value of Hamiltonian constraint (for log plots).

  if (allocated(hamabs)) then
     if (allocated(ham)) then
        hamabs = abs(ham)
     else
        print *
        print *, 'hamabs can not be output without ham!'
        print *
        call die
     end if
  end if


! ********************
! ***   MOMENTUM   ***
! ********************

! The momentum constraint in this case takes the form:
!
!
! mom  =  d KTA  -  (2/3) d trK  +  6 KTA d phi
!          r               r               r
!
!          2
!      +  r ( d B / B + 2/r ) Klambda  -  8 pi JA
!              r
!
! Notice that this expression in fact corresponds to mom_r,
! with a covariant (down) index.

  if (allocated(mom)) then

     mom = D1_KTA - 2.d0*third*D1_trK + 6.d0*KTA*D1_phi  &
         + r*(r*D1_B/B + 2.d0)*Klambda

     if (mattertype /= "vacuum") then
        if (contains(mattertype,"nonmin")) then
           mom = mom - JA/nonmin_f
        else
           mom = mom - 8.d0*smallpi*JA
        end if
     end if

!    Boundary. The calculation above is done all the way to
!    the boundary with one-sided differences, but this introduces
!    large artificial errors, so I fix it here.
  
     if (rank==size-1) then
        mom(:,Nr  ) = mom(:,Nr-3)
        mom(:,Nr-1) = mom(:,Nr-3)
        mom(:,Nr-2) = mom(:,Nr-3)
     end if

  end if

! Absolute value of momentum constraint (for log plots).

  if (allocated(momabs)) then
     if (allocated(mom)) then
        momabs = abs(mom)
     else
        print *
        print *, 'momabs can not be output without mom!'
        print *
        call die
     end if
  end if


! ****************************
! ***   DELTA CONSTRAINT   ***
! ****************************

! CDeltar:   A Deltar - [ (dA/dr)/(2A) - (dB/dr)/B - 2 r lambda ]

  if (allocated(CDeltar)) then
     CDeltar = A*Deltar - (0.5d0*D1_A/A - D1_B/B - 2.d0*r*lambda)
  end if


! **************************************
! ***   REGULARIZATION CONSTRAINTS   ***
! **************************************

! Clambda:   r**2*lambda - (1 - A/B)

  if (allocated(Clambda)) then
     Clambda = r**2*lambda - (1.d0 - A/B)
  end if

! CKlambda:  r**2*Klambda - (KTA - KTB)

  if (allocated(CKlambda)) then
     CKlambda = r**2*Klambda - (KTA - KTB)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine constraints
