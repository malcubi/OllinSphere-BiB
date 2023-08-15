!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/invariants.f90,v 1.1 2019/08/27 17:15:46 malcubi Exp $

  subroutine invariants

! ********************************
! ***   CURVATURE INVARIANTS   ***
! ********************************

! This routine calculates the curvature invariants I and J.
! These invariantes are calculated in terms of the electric
! part of Weyl (notice that the magnetic part of Weyl vanishes
! in spherical symmetry).

! Include modules.

  use arrays
  use param

! Extra variables.

  implicit none


! **************************
! ***   REAL PART OF I   ***
! **************************

! Real part of I:
!
!                    i   j
!    Re(I)  =  1/2  E   E
!                     j   i

  if (allocated(invariantI)) then
     invariantI = 0.5d0*(EWEYLA**2 + 2.d0*EWEYLB**2)
  end if


! **************************
! ***   REAL PART OF J   ***
! **************************

! Real part of J:
!                       j   i   k
!    Re(J)  =  - 1/6  E   E   E
!                      i   k   j

  if (allocated(invariantJ)) then
     invariantJ = - 1.d0/6.d0*(EWEYLA**3 + 2.d0*EWEYLB**3)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine invariants
