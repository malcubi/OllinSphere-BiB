!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/constraintbound.f90,v 1.32 2023/10/04 19:42:33 malcubi Exp $

  subroutine constraintbound(l)

! ********************************************
! ***   IMPOSE CONSTRAINTS AT BOUNDARIES   ***
! ********************************************

! This routine imposes the constraints at the boundaries.
!
! Notice that the boundary conditions here are applied at the level
! of the source terms, and not directly to the dynamical quantities
! themselves.
!
! Also, this routine only applies a boundary condition to KTA (and Klambda),
! as the boundary conditions for trK and Deltar have already been applied
! when we get here.
!
! NOTE: Even though the grid level l is passed as an argument,
! the routine is really only ever called for l=0.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer l

  real(8) aux
  real(8) zero,one,two,half,third,smallpi


! **********************************
! ***   DO WE WANT TO BE HERE?   ***
! **********************************

  if (spacetime/="dynamic") return


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  one   = 1.d0
  two   = 2.d0
  half  = 0.5d0
  third = 1.d0/3.d0

  smallpi = acos(-1.d0)


! ************************
! ***   SANITY CHECK   ***
! ************************

  if (rank==0) then

!    Only allow standard BSSN.

     if (eta/=2.d0) then
        print *
        print *, 'boundtype=constraint only works for eta=2.'
        print *, 'Aborting! (subroutine constraintbound)'
        print *
        call die
     end if

  end if


! **********************************
! ***   ONLY PHYSICAL BOUNDARY   ***
! **********************************

! Only the processor at the physical boundary applies boundary conditions.

  if ((rank==size-1).and.(l==0)) then


!    ********************************************
!    ***   APPLY BOUNDARY CONDITION FOR KTA   ***
!    ********************************************

!    Here is where we impose the constraint preserving boundary
!    conditions. At the moment this always assumes eta=2!
!
!    The outgoing and incoming constraint eigenfields at the boundary are:
!
!    wp = sqrt(A) psi**2 ( KTA - 2/3 trK ) + ( D1_A/A - D1_B/B - 2 A Deltar + 4 D1_phi ) / 3
!    wm = sqrt(A) psi**2 ( KTA - 2/3 trK ) - ( D1_A/A - D1_B/B - 2 A Deltar + 4 D1_phi ) / 3
!
!    From this we easily find:
!
!    KTA = 2/3 trK + (whp + whm) / (2 sqrt(A) psi**2)
!
!    We can then calculate the sources for whp and whm at the boundary in terms
!    of the constraints, assume that the constraints vanish, and keep the remaining
!    terms to find the source of KTA at the boundary. I do this with Maple, and
!    the result is what is coded below.
!
!    The last element is to add again a multiple of the momentum constraint to model
!    an outgoing wave moving at the speed of light at the boundary. This reduces
!    the late-time drift considerably.
!
!    Notice that by the time we get here we already have the source for trK at the
!    boundary which is calculated in the routine "radiative_geometry.f90".

     aux = one/A(l,Nr)/psi4(l,Nr)

     sKTA(l,Nr) = two*third*strK(l,Nr) &
                + aux*alpha(l,Nr)*((4.d0*D1_phi(l,Nr) - 0.5d0*A(l,Nr)*Deltar(l,Nr) &
                + 0.25d0*D1_AB2(l,Nr)/AB2(l,Nr))/r(l,Nr) &
                + (two*D1_phi(l,Nr) + half*D1_B(l,Nr)/B(l,Nr))**2) &
                + alpha(l,Nr)*(KTA(l,Nr)*trK(l,Nr) - 0.75d0*KTA(l,Nr)**2 - third*trK(l,Nr)**2) &
                + aux*D1_alpha(l,Nr)*(two/r(l,Nr) + 4.d0*D1_phi(l,Nr) + D1_B(l,Nr)/B(l,Nr))

!    Add shift terms.

     if (shift/="none") then
        sKTA(l,Nr) = sKTA(l,Nr) - 3.d0*beta(l,Nr)*KTA(l,Nr) &
                   *(2.d0*D1_phi(l,Nr) + half*D1_B(l,Nr)/B(l,Nr) + 1.d0/r(l,Nr))
     end if

!    Add matter terms.

     if (mattertype/="vacuum") then

        sKTA(l,Nr) = sKTA(l,Nr) - 8.d0*smallpi*alpha(l,Nr)*SAA(l,Nr)

!       Term with matter and shift. This term is only here because
!       we eliminated the momentum constraint.

        if (shift/="none") then
           sKTA(l,Nr) = sKTA(l,Nr) + 8.d0*smallpi*beta(l,Nr)*JA(l,Nr)
        end if

     end if

!    Now add again a multiple of the momentum constraint in order to have
!    an outgoing wave moving at the speed of light at the boundary.
!
!    For ordinary runs we take the speed of light to be vl=1 at the boundary,
!    and for cosmological runs we take vl=(alpha_bg/a_bg), where alpha_bg and
!    a_bg are the background lapse and scale factor.

     if (.not.cosmic_run) then
        aux = 1.d0
     else
        aux = cosmobg_alpha(l)/cosmobg_a(l)
     end if

     sKTA(l,Nr) = sKTA(l,Nr) - aux*(D1_KTA(l,Nr) &
                - 2.d0*third*D1_trK(l,Nr) + 6.d0*KTA(l,Nr)*D1_phi(l,Nr) &
                + 1.5d0*(D1_B(l,Nr)/B(l,Nr) + 2.d0/r(l,Nr))*KTA(l,Nr))

     if (mattertype/="vacuum") then
        sKTA(l,Nr) = sKTA(l,Nr) + 8.d0*smallpi*JA(l,Nr)*aux
     end if

!    Add multiple of Delta constraint.  This also seems to reduce late time drift.
!    I found this by trial and error, not sure why it works.

     sKTA(l,Nr) = sKTA(l,Nr) - (Deltar(l,Nr) - DeltaAB(l,Nr))


!    ******************************************
!    ***   BOUNDARY CONDITION FOR Klambda   ***
!    ******************************************

!    Having found the source for KTA at the boundary,
!    we now recover the source for Klambda.

     if (.not.nolambda) then
        sKlambda(l,Nr) = 1.5d0*sKTA(l,Nr)/r(l,Nr)**2
        if (regular2) then
           sKlambda2(l,Nr) = (sKlambda(l,Nr) - 4.d0*Klambda(l,Nr)*sphi(l,Nr))/psi4(l,Nr)
        end if
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine constraintbound

