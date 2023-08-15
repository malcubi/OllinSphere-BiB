!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/virial.f90,v 1.3 2022/05/24 16:55:52 malcubi Exp $

  subroutine virialintegral1

! *************************************************
! ***   CALCULATION OF VIRIAL INTEGRAL TYPE 1   ***
! *************************************************

! Here we calculate a virial-like integral that must vanish for
! regular static spacetimes with a metric in the areal gauge:
!
!   2        2     2       2
! dl  =  A dr  +  r d Omega
!
! That is, B=psi=1 (we can probably modify the expression so that
! it works in the more general case, I will check this later.)
!
! The virial-like integral in this case is:
!
!        infinity
!       /
! I  =  |  [ (d alpha / alpha + d A / 2A) d alpha / (alpha A)  -  4 pi tr(S) ] dV
!       /      r                 r         r
!        0
!
! with tr(S) = SAA + 2*SBB the trace of the stress tensor.
!
! Notice that the integration above must be done using the
! FLAT volume element.  We ignore the integral over the
! angles that just gives 4*pi, so our volume element is:
!
!         2
! dV  =  r  dr
!
! Altough the integral above is not directly the virial theorem
! (which is much better expressed in the conformally flat case),
! it is a related integral condition.
!
! The viriao integral is useful to check initial data for solutions
! that are supposed to be static.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use integrals

! Extra variables.

  implicit none

  integer l

  real(8) one,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  one = 1.d0

  smallpi = acos(-one)


! *********************
! ***   INTEGRATE   ***
! *********************

! Find integrand.

  auxarray = r**2*((D1_alpha/alpha + 0.5d0*D1_A/A)*D1_alpha/(alpha*A) &
           - 4.d0*smallpi*(SAA + 2.d0*SBB))

! Integrate.

  intvar => auxarray

  do l=Nl-1,0,-1
     virial1(l,:) = integral(l)
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine virialintegral1









  subroutine virialintegral2

! *************************************************
! ***   CALCULATION OF VIRIAL INTEGRAL TYPE 2   ***
! *************************************************

! Here we calculate a virial integral that must vanish for
! conformally-flat regular static spacetimes with spatial
! metric of the form:
!
!   2       4     2     2       2
! dl  =  psi  ( dr  +  r d Omega )
!
! That is, A=B=1 (we can probably modify the expression so that
! it works in the more general case, I will check this later.)
!
! The virial integral in this case has the form:
!
!        infinity
!       /             r             2             r         2
! I  =  |  [ d alpha d alpha / alpha  -  2 d psi d psi / psi  -  4 pi tr(S) ] dV
!       /     r                             r
!       0
!
! with tr(S) = SAA + 2*SBB the trace of the stress tensor.
!
! Notice that the integration abovr must be done using the 
! PHYSICAL volume element.  We ignore the integral over the
! angles that just gives 4*pi, so our volume element is:
!
!         2    6
! dV  =  r  psi  dr
!
! The integral above is in fact directly the virial theorem
! for the GR case.  Do notice that there are extra factors
! of 1/psi**4 in the first two terms due to the fact that
! we need to raise the index in some of the derivatives.
!
! The virial integral is useful to check initial data for
! solutions that are supposed to be static.

! Include modules.

  use param
  use arrays
  use integrals

! Extra variables.

  implicit none

  integer l

  real(8) one,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  one = 1.d0

  smallpi = acos(-one)


! *********************
! ***   INTEGRATE   ***
! *********************

  auxarray = r**2*psi**6*(((D1_alpha/alpha)**2 - 2.d0*(D1_psi/psi)**2)/psi**4 &
           - 4.d0*smallpi*(SAA + 2.d0*SBB))

! Integrate.

  intvar => auxarray

  do l=Nl-1,0,-1
     virial2(l,:) = integral(l)
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine virialintegral2
