!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/derivadvect.f90,v 1.17 2023/02/14 18:50:38 malcubi Exp $

  module derivadvect

  contains

! I use a Fortran module for the array-valued function that
! calculates the advective derivative of arrays.
!
! It turns out that defining array-valued functions is not
! trivial, and putting them inside a module seems to be
! a way to solve the problem.

  function diffadv(l,sym,origin)

! ***************************************************************
! ***   CALCULATE ADVECTIVE FIRST DERIVATIVE TO GIVEN ORDER   ***
! ***************************************************************

! This routine calculates the advective (one-sided) first derivative
! of the array "var" based on the sign of the array "beta".

  use procinfo
  use param
  use arrays

  implicit none

  logical oflag
  logical,optional :: origin

  integer i,l,sym

  real(8) idr,hidr

  real(8), dimension (1-ghost:Nrmax), target :: diffadv


! *******************
! ***   NUMBERS   ***
! *******************

  idr  = 1.d0/dr(l)
  hidr = 0.5d0*idr


! ***************************
! ***   ORIGIN BEHAVIOR   ***
! ***************************

  if (present(origin)) then
     oflag = origin
  else
     oflag = .false.
  end if


! **********************
! ***   INITIALIZE   ***
! **********************

  diffadv = 0.d0


! ************************
! ***   SECOND ORDER   ***
! ************************

  if (order=="two") then

!    Interior points:  Second order fully one-sided derivatives
!    depending on sign of the shift.

     do i=1,Nr-2

        if (beta(l,i)<=0.d0) then
           diffadv(i) = + hidr*(3.d0*diffvar(l,i) - 4.d0*diffvar(l,i-1) + diffvar(l,i-2))
        else
           diffadv(i) = - hidr*(3.d0*diffvar(l,i) - 4.d0*diffvar(l,i+1) + diffvar(l,i+2))
        end if

     end do

!    Point Nr-1:  Second order fully one-sided for negative shift,
!    second order centered for positive shift.

     i = Nr-1

     if (beta(l,i)<=0.d0) then
        diffadv(i) = hidr*(3.d0*diffvar(l,i) - 4.d0*diffvar(l,i-1) + diffvar(l,i-2))
     else
        diffadv(i) = hidr*(diffvar(l,i+1) - diffvar(l,i-1))
     end if

!    Point Nr:  Second order fully one-sided always.

     i = Nr
     diffadv(i) = hidr*(3.d0*diffvar(l,i) - 4.d0*diffvar(l,i-1) + diffvar(l,i-2))

!    Check if we want one-sided derivatives at origin.

     if (oflag) then
        i = 1
        diffadv(i) = - hidr*(3.d0*diffvar(l,i) - 4.d0*diffvar(l,i+1) + diffvar(l,i+2))
     end if


! ************************
! ***   FOURTH ORDER   ***
! ************************

  else if (order=="four") then

!    Interior points:  Fourth order semi one-sided derivatives (1 point upwind,
!    3 downwind) depending on sign of the shift.  This works well and is more
!    accurate than using fully one-sided derivatives.

     do i=1,Nr-3

        if (beta(l,i)<=0.d0) then
           diffadv(i) = + idr*(3.d0*diffvar(l,i+1) + 10.d0*diffvar(l,i) &
                      - 18.d0*diffvar(l,i-1) + 6.d0*diffvar(l,i-2) - diffvar(l,i-3))/12.d0
        else
           diffadv(i) = - idr*(3.d0*diffvar(l,i-1) + 10.d0*diffvar(l,i) &
                      - 18.d0*diffvar(l,i+1) + 6.d0*diffvar(l,i+2) - diffvar(l,i+3))/12.d0
        end if

     end do

!    Point Nr-2:  Fourth order semi one-sided derivatives for negative shift,
!    fourth order centered for positive shift.

     i = Nr-2

     if (beta(l,i)<=0.d0) then
        diffadv(i) = idr*(3.d0*diffvar(l,i+1) + 10.d0*diffvar(l,i) &
                   - 18.d0*diffvar(l,i-1) + 6.d0*diffvar(l,i-2) - diffvar(l,i-3))/12.d0
     else
        diffadv(i) = idr*(8.d0*(diffvar(l,i+1) - diffvar(l,i-1)) &
                   - (diffvar(l,i+2) - diffvar(l,i-2)))/12.d0
     end if

!    Point Nr-1: Sixth order semi one-sided always (1 point to the right, 5 to the left).
!    I use sixth order to improve boundary behaviour.

     i = Nr-1

     diffadv(i) = idr*(1.d0/6.d0*diffvar(l,i+1) + 77.d0/60.d0*diffvar(l,i) - 5.d0/2.d0*diffvar(l,i-1) &
                + 5.d0/3.d0*diffvar(l,i-2) - 5.d0/6.d0*diffvar(l,i-3) + 0.25*diffvar(l,i-4) &
                - 1.d0/30.d0*diffvar(l,i-5)) ! Sixth order semi one-sided

!    Point Nr:  Sixth order fully one-sided always. 
!    I use sixth order to improve boundary behaviour.

     i = Nr

     diffadv(i) = idr*(49.d0/20.d0*diffvar(l,i) - 6.d0*diffvar(l,i-1) + 7.5d0*diffvar(l,i-2) &
                - 20.d0/3.d0*diffvar(l,i-3) + 15.d0/4.d0*diffvar(l,i-4) - 6.d0/5.d0*diffvar(l,i-5) &
                + 1.d0/6.d0*diffvar(l,i-6)) ! Sixth order fully one-sided

!    Check if we want one-sided derivatives at origin.

     if (oflag) then

        i = 1
        diffadv(i) = - idr*(49.d0/20.d0*diffvar(l,i) - 6.d0*diffvar(l,i+1) + 7.5d0*diffvar(l,i+2) &
                   - 20.d0/3.d0*diffvar(l,i+3) + 15.d0/4.d0*diffvar(l,i+4) - 6.d0/5.d0*diffvar(l,i+5) &
                   + 1.d0/6.d0*diffvar(l,i+6)) ! Sixth order fully one-sided

        i = 2
        diffadv(i) = - idr*(1.d0/6.d0*diffvar(l,i-1) + 77.d0/60.d0*diffvar(l,i) - 5.d0/2.d0*diffvar(l,i+1) &
                   + 5.d0/3.d0*diffvar(l,i+2) - 5.d0/6.d0*diffvar(l,i+3) + 0.25*diffvar(l,i+4) &
                   - 1.d0/30.d0*diffvar(l,i+5)) ! Sixth order semi one-sided

     end if


! ***********************
! ***   SIXTH ORDER   ***
! ***********************

  else if (order=="six") then

!    Interior points:  Sixth order semi-semi one-sided derivatives (2 points upwind,
!    4 downwind) depending on sign of the shift.  This works well and is more
!    accurate than using fully one-sided derivatives.

     do i=1,Nr-4
        if (beta(l,i)<=0.d0) then
           diffadv(i) = + idr*(-1.d0/30.d0*diffvar(l,i+2) + 2.d0/5.d0*diffvar(l,i+1) + 7.d0/12.d0*diffvar(l,i) &
                      - 4.d0/3.d0*diffvar(l,i-1) + 0.5d0*diffvar(l,i-2) - 2.d0/15.d0*diffvar(l,i-3) &
                      + 1.d0/60.d0*diffvar(l,i-4))
        else
           diffadv(i) = - idr*(-1.d0/30.d0*diffvar(l,i-2) + 2.d0/5.d0*diffvar(l,i-1) + 7.d0/12.d0*diffvar(l,i) &
                      - 4.d0/3.d0*diffvar(l,i+1) + 0.5d0*diffvar(l,i+2) - 2.d0/15.d0*diffvar(l,i+3) &
                      + 1.d0/60.d0*diffvar(l,i+4))
        end if
     end do

!    Point Nr-3:  Sixth order semi-semi one-sided derivatives for negative shift,
!    sixth order centered for positive shift.

     i = Nr-3

     if (beta(l,i)<=0.d0) then
        diffadv(i) = idr*(-1.d0/30.d0*diffvar(l,i+2) + 2.d0/5.d0*diffvar(l,i+1) + 7.d0/12.d0*diffvar(l,i) &
                   - 4.d0/3.d0*diffvar(l,i-1) + 0.5d0*diffvar(l,i-2) - 2.d0/15.d0*diffvar(l,i-3) &
                   + 1.d0/60.d0*diffvar(l,i-4))
     else
        diffadv(i) = idr*(45.d0*(diffvar(l,i+1) - diffvar(l,i-1)) &
                   -       9.d0*(diffvar(l,i+2) - diffvar(l,i-2)) &
                   +            (diffvar(l,i+3) - diffvar(l,i-3)))/60.d0
     end if

!    Point Nr-2:  Sixth order semi-semi one-sided always (2 points to the right, 4 to the left).

     i = Nr-2

     diffadv(i) = idr*(-1.d0/30.d0*diffvar(l,i+2) + 2.d0/5.d0*diffvar(l,i+1) + 7.d0/12.d0*diffvar(l,i) &
                - 4.d0/3.d0*diffvar(l,i-1) + 0.5d0*diffvar(l,i-2) - 2.d0/15.d0*diffvar(l,i-3) &
                + 1.d0/60.d0*diffvar(l,i-4))

!    Point Nr-1:  Sixth order semi one-sided always (1 point to the right, 5 to the left).

     i = Nr-1

     diffadv(i) = idr*(1.d0/6.d0*diffvar(l,i+1) + 77.d0/60.d0*diffvar(l,i) - 5.d0/2.d0*diffvar(l,i-1) &
                + 5.d0/3.d0*diffvar(l,i-2) - 5.d0/6.d0*diffvar(l,i-3) + 0.25*diffvar(l,i-4) &
                - 1.d0/30.d0*diffvar(l,i-5))

!    Point Nr:  Sixth order fully one-sided always.

     i = Nr

     diffadv(i) = idr*(49.d0/20.d0*diffvar(l,i) - 6.d0*diffvar(l,i-1) + 7.5d0*diffvar(l,i-2) &
                - 20.d0/3.d0*diffvar(l,i-3) + 15.d0/4.d0*diffvar(l,i-4) - 6.d0/5.d0*diffvar(l,i-5) &
                + 1.d0/6.d0*diffvar(l,i-6))

!    Check if we want one-sided derivatives at origin.

     if (oflag) then

        i = 1
        diffadv(i) = - idr*(49.d0/20.d0*diffvar(l,i) - 6.d0*diffvar(l,i+1) + 7.5d0*diffvar(l,i+2) &
                   - 20.d0/3.d0*diffvar(l,i+3) + 15.d0/4.d0*diffvar(l,i+4) - 6.d0/5.d0*diffvar(l,i+5) &
                   + 1.d0/6.d0*diffvar(l,i+6))

        i = 2
        diffadv(i) = - idr*(1.d0/6.d0*diffvar(l,i-1) + 77.d0/60.d0*diffvar(l,i) - 5.d0/2.d0*diffvar(l,i+1) &
                   + 5.d0/3.d0*diffvar(l,i+2) - 5.d0/6.d0*diffvar(l,i+3) + 0.25*diffvar(l,i+4) &
                   - 1.d0/30.d0*diffvar(l,i+5))

        i = 3
        diffadv(i) = - idr*(-1.d0/30.d0*diffvar(l,i-2) + 2.d0/5.d0*diffvar(l,i-1) + 7.d0/12.d0*diffvar(l,i) &
                   - 4.d0/3.d0*diffvar(l,i+1) + 0.5d0*diffvar(l,i+2) - 2.d0/15.d0*diffvar(l,i+3) &
                   + 1.d0/60.d0*diffvar(l,i+4))

     end if


! ************************
! ***   EIGHTH ORDER   ***
! ************************

  else if (order=="eight") then

!    Interior points:  Eighth order semi-semi-semi left-sided (3 points upwind,
!    5 downwind) depending on sign of the shift.  This works well and is more
!    accurate than using fully one-sided derivatives.

     do i=1,Nr-5
        if (beta(l,i)<=0.d0) then
           diffadv(i) = + idr*(1.d0/168.d0*diffvar(l,i+3) - 1.d0/14.d0*diffvar(l,i+2) + 0.5d0*diffvar(l,i+1) &
                      + 9.d0/20.d0*diffvar(l,i) - 5.d0/4.d0*diffvar(l,i-1) + 0.5d0*diffvar(l,i-2) &
                      - 1.d0/6.d0*diffvar(l,i-3) + 1.d0/28.d0*diffvar(l,i-4) - 1.d0/280.d0*diffvar(l,i-5))
        else
           diffadv(i) = - idr*(1.d0/168.d0*diffvar(l,i-3) - 1.d0/14.d0*diffvar(l,i-2) + 0.5d0*diffvar(l,i-1) &
                      + 9.d0/20.d0*diffvar(l,i) - 5.d0/4.d0*diffvar(l,i+1) + 0.5d0*diffvar(l,i+2) &
                      - 1.d0/6.d0*diffvar(l,i+3) + 1.d0/28.d0*diffvar(l,i+4) - 1.d0/280.d0*diffvar(l,i+5))
        end if
     end do

!    Point Nr-4:  Eighth order semi-semi-semi one-sided derivatives for negative shift,
!    eight order centered for positive shift.

     i = Nr-4

     if (beta(l,i)<=0.d0) then
        diffadv(i) = + idr*(1.d0/168.d0*diffvar(l,i+3) - 1.d0/14.d0*diffvar(l,i+2) + 0.5d0*diffvar(l,i+1) &
                   + 9.d0/20.d0*diffvar(l,i) - 5.d0/4.d0*diffvar(l,i-1) + 0.5d0*diffvar(l,i-2) &
                   - 1.d0/6.d0*diffvar(l,i-3) + 1.d0/28.d0*diffvar(l,i-4) - 1.d0/280.d0*diffvar(l,i-5))
     else
        diffadv(i) = idr*( 0.8d0*(diffvar(l,i+1) - diffvar(l,i-1)) &
                   -       0.2d0*(diffvar(l,i+2) - diffvar(l,i-2)) &
                   + 4.d0/105.d0*(diffvar(l,i+3) - diffvar(l,i-3)) &
                   - 1.d0/280.d0*(diffvar(l,i+4) - diffvar(l,i-4)))
     end if

!    Point Nr-3:  Eighth order semi-semi-semi one-sided always (3 points to the right, 5 to the left).

     i = Nr-3

     diffadv(i) = idr*(1.d0/168.d0*diffvar(l,i+3) - 1.d0/14.d0*diffvar(l,i+2) + 0.5d0*diffvar(l,i+1) &
                + 9.d0/20.d0*diffvar(l,i) - 5.d0/4.d0*diffvar(l,i-1) + 0.5d0*diffvar(l,i-2) &
                - 1.d0/6.d0*diffvar(l,i-3) + 1.d0/28.d0*diffvar(l,i-4) - 1.d0/280.d0*diffvar(l,i-5))

!    Point Nr-2:  Eighth order semi-semi one-sided always (2 points to the right, 6 to the left).

     i = Nr-2

     diffadv(i) = idr*(-1.d0/56.d0*diffvar(l,i+2) + 2.d0/7.d0*diffvar(l,i+1) + 19.d0/20.d0*diffvar(l,i) &
                - 2.d0*diffvar(l,i-1) + 1.25d0*diffvar(l,i-2) - 2.d0/3.d0*diffvar(l,i-3) &
                + 0.25d0*diffvar(l,i-4) - 2.d0/35.d0*diffvar(l,i-5) + 1.d0/168.d0*diffvar(l,i-6))

!    Point Nr-1:  Eighth order semi one-sided always (1 point to the right, 7 to the left).

     i = Nr-1

     diffadv(i) = idr*(1.d0/8.d0*diffvar(l,i+1) + 223.d0/140.d0*diffvar(l,i) - 3.5d0*diffvar(l,i-1) &
                + 3.5d0*diffvar(l,i-2) - 35.d0/12.d0*diffvar(l,i-3) + 1.75d0*diffvar(l,i-4) &
                - 0.7d0*diffvar(l,i-5) + 1.d0/6.d0*diffvar(l,i-6) - 1.d0/56.d0*diffvar(l,i-7))

!    Point Nr:  Sixth order fully one-sided always.  The eighth order fully one-sided difference
!    (commented out) causes the runs to go unstable.

     i = Nr

     !diffadv(i) = idr*(761.d0/280.d0*diffvar(l,i) - 8.d0*diffvar(l,i-1) + 14.d0*diffvar(l,i-2) &
     !           - 56.d0/3.d0*diffvar(l,i-3) + 35.d0/2.d0*diffvar(l,i-4) - 56.d0/5.d0*diffvar(l,i-5) &
     !           + 14.d0/3.d0*diffvar(l,i-6) - 8.d0/7.d0*diffvar(l,i-7) + 1.d0/8.d0*diffvar(l,i-8)) ! Eighth order

     diffadv(i) = idr*(49.d0/20.d0*diffvar(l,i) - 6.d0*diffvar(l,i-1) + 7.5d0*diffvar(l,i-2) &
                - 20.d0/3.d0*diffvar(l,i-3) + 15.d0/4.d0*diffvar(l,i-4) - 6.d0/5.d0*diffvar(l,i-5) &
                + 1.d0/6.d0*diffvar(l,i-6)) ! Sixth order fully one-sided

!    Check if we want one-sided derivatives at origin.

     if (oflag) then

        i = 1
        diffadv(i) = - idr*(49.d0/20.d0*diffvar(l,i) - 6.d0*diffvar(l,i+1) + 7.5d0*diffvar(l,i+2) &
                   - 20.d0/3.d0*diffvar(l,i+3) + 15.d0/4.d0*diffvar(l,i+4) - 6.d0/5.d0*diffvar(l,i+5) &
                   + 1.d0/6.d0*diffvar(l,i+6)) ! Sixth order fully one-sided

        i = 2
        diffadv(i) = - idr*(1.d0/8.d0*diffvar(l,i-1) + 223.d0/140.d0*diffvar(l,i) - 3.5d0*diffvar(l,i+1) &
                   + 3.5d0*diffvar(l,i+2) - 35.d0/12.d0*diffvar(l,i+3) + 1.75d0*diffvar(l,i+4) &
                   - 0.7d0*diffvar(l,i+5) + 1.d0/6.d0*diffvar(l,i+6) - 1.d0/56.d0*diffvar(l,i+7))

        i = 3
        diffadv(i) = - idr*(-1.d0/56.d0*diffvar(l,i-2) + 2.d0/7.d0*diffvar(l,i-1) + 19.d0/20.d0*diffvar(l,i) &
                   - 2.d0*diffvar(l,i+1) + 1.25d0*diffvar(l,i+2) - 2.d0/3.d0*diffvar(l,i+3) &
                   + 0.25d0*diffvar(l,i+4) - 2.d0/35.d0*diffvar(l,i+5) + 1.d0/168.d0*diffvar(l,i+6))

        i = 4
        diffadv(i) = - idr*(1.d0/168.d0*diffvar(l,i-3) - 1.d0/14.d0*diffvar(l,i-2) + 0.5d0*diffvar(l,i-1) &
                   + 9.d0/20.d0*diffvar(l,i) - 5.d0/4.d0*diffvar(l,i+1) + 0.5d0*diffvar(l,i+2) &
                   - 1.d0/6.d0*diffvar(l,i+3) + 1.d0/28.d0*diffvar(l,i+4) - 1.d0/280.d0*diffvar(l,i+5))

     end if

  end if


! ********************************
! ***   SYMMETRIES AT ORIGIN   ***
! ********************************

 if (rank==0) then
     do i=1,ghost
        diffadv(1-i) = - sym*diffadv(i)
     end do
  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

  if (size>1) then
     syncvar => diffadv
     call sync
  end if


! ***************
! ***   END   ***
! ***************

  end function diffadv




! **********************
! ***   END MODULE   ***
! **********************

  end module derivadvect
