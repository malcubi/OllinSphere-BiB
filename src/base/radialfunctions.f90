!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/radialfunctions.f90,v 1.6 2021/09/03 17:17:19 malcubi Exp $

  module radialfunctions

  contains

! I use a Fortran module for array-valued radial functions. 
!
! It turns out that defining array-valued functions is not
! trivial, and putting them inside a module seems to be
! a way to solve the problem.



! ********************
! ***   GAUSSIAN   ***
! ********************

  function gaussian(ga_a0,ga_r0,ga_s0)

! This is a simple gaussian profile:
!
! F  =  a0 exp[-(r-r0)**2/s0**2]
!
! Here r0 is the position of the center of the gaussian,
! a0 is the height, and s0 the width.

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  real(8) ga_a0,ga_r0,ga_s0

  real(8), dimension (0:Nl-1,1-ghost:Nrmax), target :: gaussian

! Gaussian function.

  if (newr) then
     gaussian = ga_a0*exp(-(r_trans-ga_r0)**2/ga_s0**2)
  else
     gaussian = ga_a0*exp(-(r-ga_r0)**2/ga_s0**2)
  end if

! End.

  end function gaussian





! *********************
! ***   GAUSSIAN 2  ***
! *********************

  function r2gaussian(ga_a,ga_r0,ga_s)

! This is a gaussian profile multiplied with r**2.
!
! F  =  a r**2 exp[-(r-r0)**2/s**2]
!
! Here r0 is the position of the center of the gaussian,
! a is the height, and s the width.

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  real(8) ga_a,ga_r0,ga_s

  real(8), dimension (0:Nl-1,1-ghost:Nrmax), target :: r2gaussian

! Gaussian function.

  if (newr) then
     r2gaussian = ga_a*r_trans**2*exp(-(r_trans-ga_r0)**2/ga_s**2)
  else
     r2gaussian = ga_a*r**2*exp(-(r-ga_r0)**2/ga_s**2)
  end if

! End.

  end function r2gaussian





! ****************************
! ***   TOP HAT FUNCTION   ***
! ****************************

  function tophat(th_a0,th_r0,th_s0,th_t0)

! This is a simple smooth top-hat profile built using
! hyperbolic tangents and defined as:
!
! F  =  (a0/2) [ tanh((r-r0+s0)/t0) - tanh((r-r0-s0)/t0) ] / tanh(s0/t0)
!
! Here a is the height of the function at r=r0, t is the
! transition radius, and s the width of the transition.

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  real(8) th_a0,th_r0,th_s0,th_t0

  real(8), dimension (0:Nl-1,1-ghost:Nrmax), target :: tophat

! Sanity check.

  if (th_t0==0.d0) then
    if (rank==0) then
        print *, 'For a tophat profile t0 marks the transition radius and can not be zero'
        print *, 'Aborting!  (subroutine radialfunctions.f90)'
        print *
     end if
     call die
  end if

! Tophat function.

  if (newr) then
     tophat = 0.5d0*th_a0*(tanh((r_trans-th_r0+th_s0)/th_t0) - tanh((r_trans-th_r0-th_s0)/th_t0))/tanh(th_s0/th_t0)
  else
     tophat = 0.5d0*th_a0*(tanh((r-th_r0+th_s0)/th_t0) - tanh((r-th_r0-th_s0)/th_t0))/tanh(th_s0/th_t0)
  end if

! End.

  end function tophat






! **********************
! ***   END MODULE   ***
! **********************

  end module radialfunctions
