!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/auxiliary.f90,v 1.7 2024/09/23 18:49:10 malcubi Exp $

  subroutine auxiliary(l)

! *******************************************************
! ***   AUXILARY QUANTITIES FOR GEOMETRY AND MATTER   ***
! *******************************************************

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  logical contains

  integer l


! ********************
! ***   GEOMETRY   ***
! ********************

! Update auxiliary geometric variables.  For initial
! data we always have to do it, even if we are evolving
! on a background spacetime, since we need to calculate
! things like "psi" and spatial dervatives of geometric
! quantities, etc.

  if ((spacetime=="dynamic").or.(t(0)==0.d0)) then
     call auxiliary_geometry(l)
  end if


! ******************
! ***   MATTER   ***
! ******************

! Update auxiliary matter variables.

  if (mattertype/="vacuum") then
     call auxiliary_matter(l)
  end if


! ********************************
! ***   STRESS-ENERGY TENSOR   ***
! ********************************

  if (mattertype/="vacuum") then

!    Potential for scalar fields.

     if ((index(mattertype,"scalar") /=0).or.(index(mattertype,"ghost") /=0).or. &
         (index(mattertype,"complex")/=0).or.(index(mattertype,"nonmin")/=0)) then
        call potential(l)
     end if

!    Recover fluid primitive variables.

     if (contains(mattertype,"fluid")) then
        call fluidprimitive(l)
     end if

!    Non-minimal coupling (before stress-energy).

     if (index(mattertype,"nonmin")/=0) then
        call nonmincoupling
     end if

!    Calculate stress-energy variables.

     call stressenergy(l)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine auxiliary
