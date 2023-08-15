!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/integrals.f90,v 1.10 2022/07/12 23:58:01 malcubi Exp $

  module integrals

  contains

! I use a Fortran module for the array-valued function that
! calculates the integral of arrays.
!
! It turns out that defining array-valued functions is not
! trivial, and putting them inside a module seems to be
! a way to solve the problem.


  function integral(l)

! **************************
! ***   INTEGRATE IN r   ***
! **************************

! This subroutine does an integral in r of an array.
! The integrated array is always assumed to be associated
! with the pointer "intvar".

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer l,p
  integer i,i0

  real(8) f0,f1,f2
  real(8) c0,c1,c2
  real(8) half

  real(8), dimension (1-ghost:Nrmax), target :: integral


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0


! **********************
! ***   INITIALIZE   ***
! **********************

  integral = 0.d0


! *********************************
! ***   SECOND ORDER INTEGRAL   ***
! *********************************

! For second order we integrate using the trapezoidal rule.
! Notice that even though I have coded this, by default
! the integral is always donde to 4th order. If you want
! to test the second order case, you can uncomment
! the if statement below.

! if (order=="two") then
  if (.false.) then

     !print *, 'Second order integration'

!    Integral from 0 to dr/2:  We assume the integrand is a
!    linear function of the form:
!
!    f  =  c0 + c1 x
!
!    We obtain the coefficients c1,c2 from the first two grid
!    points and integrate exactly:
!
!    int(f) =  c0 x + c1 x^2/2
!
!    This should be exact for functions that are linear close
!    to the origin, but will have a second order error for
!    quadratic functions.

     if (rank==0) then
        c0 = half*(3.d0*intvar(l,1) - intvar(l,2))
        c1 = (intvar(l,2) - intvar(l,1))/dr(l)
        integral(1) = c0*dr(l)/2.d0 + c1*dr(l)**2/8.d0
     end if

!    Integrate over interior points.

     if (rank==0) then
        i0 = 2
     else
        i0 = 1
     end if

     do i=i0,Nrl(rank)
        f0 = intvar(l,i-1)
        f1 = intvar(l,i  )
        integral(i) = integral(i-1) + half*dr(l)*(f0 + f1)
     end do


! *********************************
! ***   FOURTH ORDER INTEGRAL   ***
! *********************************

! For fourth order we integrate using Simpson's rule,
! but I do it by interpolating the integrand at the mid-point
! to fourth order.  If I try it with the standard 2*dr interval
! it is not very accurate and I get red-black noise.
!
! Notice that by default the integral is always done to fourth
! order.

! else if (order=="four") then
  else if (.true.) then

     !print *, 'Fourth order integration'

!    At the origin we must take care, since we should not assume
!    that the integrand makes sense for the ghost points with
!    negative r.  Here we use cubic extrapolation to find the
!    value of the integrand at the origin r=0 and at r=dr/4 using
!    the first four grid points to the right of the origin.

     if (rank==0) then

!       Integral from 0 to r=dr/2 (i=1).

        f0 = (35.d0*(intvar(l,1) - intvar(l,2)) &
           + 21.d0*intvar(l,3) - 5.d0*intvar(l,4))/16.d0    ! Extrapolation to r=0.
        f2 = intvar(l,1)

        f1 = (195.d0*intvar(l,1) - 117.d0*intvar(l,2) &
           + 65.d0*intvar(l,3) - 15.d0*intvar(l,4))/128.d0  ! Extrapolation to r=dr/4.

        integral(1) = dr(l)*(f0 + 4.d0*f1 + f2)/12.d0

!       Integral from dr/2 to 3dr/2 (i=1 to i=2).

        f0 = intvar(l,1)
        f2 = intvar(l,2)

        f1 = (5.d0*intvar(l,1) + 15.d0*intvar(l,2) &
            - 5.d0*intvar(l,3) + intvar(l,4))/16.d0         ! One-sided interpolation.

        integral(2) = integral(1) + dr(l)*(f0 + 4.d0*f1 + f2)/6.d0

     end if

!    Integrate over interior points using Simpson's rule.

     if (rank==0) then
        i0 = 3
     else
        i0 = 1
     end if

     do i=i0,Nrl(rank)

        f0 = intvar(l,i-1)
        f2 = intvar(l,i  )

        if (i<Nrl(rank)) then
           f1 = (9.d0*(intvar(l,i-1)+intvar(l,i)) &
                    - (intvar(l,i-2)+intvar(l,i+1)))/16.d0  ! Centered cubic interpolation.
        else
           f1 = (5.d0*intvar(l,i) + 15.d0*intvar(l,i-1) &
               - 5.d0*intvar(l,i-2) + intvar(l,i-3))/16.d0  ! One-sided interpolation.
        end if

        integral(i) = integral(i-1) + dr(l)*(f0 + 4.d0*f1 + f2)/6.d0

     end do


! ********************************
! ***   SIXTH ORDER INTEGRAL   ***
! ********************************

  else if (order=="six") then

     print *
     print *, 'Sixth order integral not yet implemented'
     print *
     call die


! ********************************
! ***   EIGHT ORDER INTEGRAL   ***
! ********************************

  else if (order=="eight") then

     print *
     print *, 'Eighth order integral not yet implemented'
     print *
     call die

  end if


! ***********************
! ***   GHOST ZONES   ***
! ***********************

! Ghost zones using symmetries (only for processor zero).
! We assume the integral is even, but it really shouldn't
! matter, and can be fixed later if we need to.

  if (rank==0) then
     do i=1,ghost
        integral(1-i) = integral(i)
     end do
  end if


! ***********************************************
! ***   CORRECTION FOR MULTI-PROCESSOR RUNS   ***
! ***********************************************

! When we have more than one processor, each processor
! has now calculated definite integrals within its
! boundaries.  We now need to add the constant contributions
! from the processors to the left.

  if (size>1) then

     syncvar => integral

!    Loop over processors with rank>0.

     do p=1,size

!       Synchronize.

        call sync

!       If rank=p, it is time to add the accumulated constants
!       that we have picked up from the previous syncs.

        if (rank==p) then

           do i=1,Nrl(rank)
              integral(i) = integral(i) + integral(0)
           end do

        end if

     end do

!    One final sync to be sure all is OK.

     call sync

  end if


! ***************
! ***   END   ***
! ***************

  end function integral




! **********************
! ***   END MODULE   ***
! **********************

  end module integrals

