
  module derivatives

  contains

! I use a Fortran module for the array-valued functions that
! calculate the first and second derivatives of arrays.
!
! It turns out that defining array-valued functions is not
! trivial, and putting them inside a module seems to be
! a way to solve the problem.
!
! Notice that the array for which the derivatives are
! calculated is always the one that corresponds to the
! pointer "diffvar".
!
! Also, on entry "l" is the grid level and "sym" the
! symmetry of the array at the origin (+1,-1).

  function diff1(l,sym,origin)

! *****************************************************
! ***   CALCULATE FIRST DERIVATIVE TO GIVEN ORDER   ***
! *****************************************************

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  logical oflag
  logical,optional :: origin

  integer i,l,sym

  real(8) idr,hidr

  real(8), dimension (1-ghost:Nrmax), target :: diff1


! *******************
! ***   NUMBERS   ***
! *******************

  idr  = 1.d0/dr(l)
  hidr = 0.5d0*idr


! **********************
! ***   INITIALIZE   ***
! **********************

  diff1 = 0.d0


! ***************************
! ***   ORIGIN BEHAVIOR   ***
! ***************************

  if (present(origin)) then
     oflag = origin
  else
     oflag = .false.
  end if


! *********************
! ***   2ND ORDER   ***
! *********************

  if (order=="two") then

!    Interior points:  2nd order centered first derivative.

     do i=1,Nr-1
        diff1(i) = hidr*(diffvar(l,i+1) - diffvar(l,i-1))
     end do

!    Point Nr:  2nd order fully one-sided.

     i = Nr

     diff1(i) = hidr*(3.d0*diffvar(l,i) - 4.d0*diffvar(l,i-1) + diffvar(l,i-2))

!    Check if we want one-sided derivatives at origin.

     if (oflag) then
        diff1(1) = - hidr*(3.d0*diffvar(l,1) - 4.d0*diffvar(l,2) + diffvar(l,3))
     end if


! *********************
! ***   4TH ORDER   ***
! *********************

  else if (order=="four") then

!    Interior points:  4th order centered first derivative.

     do i=1,Nr-2
        diff1(i) = idr*(8.d0*(diffvar(l,i+1) - diffvar(l,i-1)) &
                 -           (diffvar(l,i+2) - diffvar(l,i-2)))/12.d0
     end do

!    Point Nr-1:  4th order semi one-sided (1 point to the right, 5¡3 to the left).

     i = Nr-1

     diff1(i) = idr*(3.d0*diffvar(l,i+1) + 10.d0*diffvar(l,i) - 18.d0*diffvar(l,i-1) &
              + 6.d0*diffvar(l,i-2) - diffvar(l,i-3))/12.d0

!    Point Nr:  4th order fully one-sided.

     i = Nr

     diff1(i) = idr*(25.d0*diffvar(l,i) - 48.d0*diffvar(l,i-1) &
              + 36.d0*diffvar(l,i-2) - 16.d0*diffvar(l,i-3) + 3.d0*diffvar(l,i-4))/12.d0

!    Check if we want one-sided derivatives at origin.

     if (oflag) then

        i = 1
        diff1(i) = - idr*(25.d0*diffvar(l,i) - 48.d0*diffvar(l,i+1) &
                + 36.d0*diffvar(l,i+2) - 16.d0*diffvar(l,i+3) + 3.d0*diffvar(l,i+4))/12.d0

        i = 2
        diff1(i) = - idr*(3.d0*diffvar(l,i+1) + 10.d0*diffvar(l,i) - 18.d0*diffvar(l,i-1) &
              + 6.d0*diffvar(l,i-2) - diffvar(l,i-3))/12.d0

     end if


! *********************
! ***   6TH ORDER   ***
! *********************

  else if (order=="six") then

!    Interior points:  6th order centered first derivative.

     do i=1,Nr-3
        diff1(i) = idr*(45.d0*(diffvar(l,i+1) - diffvar(l,i-1)) &
                 -       9.d0*(diffvar(l,i+2) - diffvar(l,i-2)) &
                 +            (diffvar(l,i+3) - diffvar(l,i-3)))/60.d0
     end do

!    Point Nr-2:  6th order semi-semi left-sided differences (2 points to the right, 4 to the left).

     i = Nr-2

     diff1(i) = idr*(-1.d0/30.d0*diffvar(l,i+2) + 2.d0/5.d0*diffvar(l,i+1) + 7.d0/12.d0*diffvar(l,i) &
              - 4.d0/3.d0*diffvar(l,i-1) + 0.5d0*diffvar(l,i-2) - 2.d0/15.d0*diffvar(l,i-3) &
              + 1.d0/60.d0*diffvar(l,i-4))

!    Point Nr-1:  6th order semi one-sided (1 point to the right, 5 to the left).

     i = Nr-1

     diff1(i) = idr*(1.d0/6.d0*diffvar(l,i+1) + 77.d0/60.d0*diffvar(l,i) - 5.d0/2.d0*diffvar(l,i-1) &
              + 5.d0/3.d0*diffvar(l,i-2) - 5.d0/6.d0*diffvar(l,i-3) + 0.25d0*diffvar(l,i-4) &
              - 1.d0/30.d0*diffvar(l,i-5))

!    Point Nr:  6th order fully one-sided.

     i = Nr

     diff1(i) = idr*(49.d0/20.d0*diffvar(l,i) - 6.d0*diffvar(l,i-1) + 7.5d0*diffvar(l,i-2) &
              - 20.d0/3.d0*diffvar(l,i-3) + 15.d0/4.d0*diffvar(l,i-4) - 6.d0/5.d0*diffvar(l,i-5) &
              + 1.d0/6.d0*diffvar(l,i-6))

!    Check if we want one-sided derivatives at origin.

     if (oflag) then

        i = 1
        diff1(i) = - idr*(49.d0/20.d0*diffvar(l,i) - 6.d0*diffvar(l,i+1) + 7.5d0*diffvar(l,i+2) &
                 - 20.d0/3.d0*diffvar(l,i+3) + 15.d0/4.d0*diffvar(l,i+4) - 6.d0/5.d0*diffvar(l,i+5) &
                 + 1.d0/6.d0*diffvar(l,i+6))

        i = 2
        diff1(i) = - idr*(1.d0/6.d0*diffvar(l,i-1) + 77.d0/60.d0*diffvar(l,i) - 5.d0/2.d0*diffvar(l,i+1) &
                 + 5.d0/3.d0*diffvar(l,i+2) - 5.d0/6.d0*diffvar(l,i+3) + 0.25d0*diffvar(l,i+4) &
                 - 1.d0/30.d0*diffvar(l,i+5))

        i = 3
        diff1(i) = - idr*(-1.d0/30.d0*diffvar(l,i-2) + 2.d0/5.d0*diffvar(l,i-1) + 7.d0/12.d0*diffvar(l,i) &
                 - 4.d0/3.d0*diffvar(l,i+1) + 0.5d0*diffvar(l,i+2) - 2.d0/15.d0*diffvar(l,i+3) &
                 + 1.d0/60.d0*diffvar(l,i+4))

     end if


! *********************
! ***   8TH ORDER   ***
! *********************

  else if (order=="eight") then

!    Interior points:  8th order centered first derivative.

     do i=1,Nr-4
        diff1(i) = idr*( 0.8d0*(diffvar(l,i+1) - diffvar(l,i-1)) &
                 -       0.2d0*(diffvar(l,i+2) - diffvar(l,i-2)) &
                 + 4.d0/105.d0*(diffvar(l,i+3) - diffvar(l,i-3)) &
                 - 1.d0/280.d0*(diffvar(l,i+4) - diffvar(l,i-4)))
     end do

!    Point Nr-3:  8th order semi-semi-semi left-sided (3 points to the right, 5 to the left).

     i = Nr-3

     diff1(i) = idr*(1.d0/168.d0*diffvar(l,i+3) - 1.d0/14.d0*diffvar(l,i+2) + 0.5d0*diffvar(l,i+1) &
               + 9.d0/20.d0*diffvar(l,i) - 5.d0/4.d0*diffvar(l,i-1) + 0.5d0*diffvar(l,i-2) &
               - 1.d0/6.d0*diffvar(l,i-3) + 1.d0/28.d0*diffvar(l,i-4) - 1.d0/280.d0*diffvar(l,i-5))

!    Point Nr-2:  8th order semi-semi left-sided differences (2 points to the right, 6 to the left).

     i = Nr-2

     diff1(i) = idr*(-1.d0/56.d0*diffvar(l,i+2) + 2.d0/7.d0*diffvar(l,i+1) + 19.d0/20.d0*diffvar(l,i) &
              - 2.d0*diffvar(l,i-1) + 1.25d0*diffvar(l,i-2) - 2.d0/3.d0*diffvar(l,i-3) &
              + 0.25d0*diffvar(l,i-4) - 2.d0/35.d0*diffvar(l,i-5) + 1.d0/168.d0*diffvar(l,i-6))

!    Point Nr-1:  8th order semi one-sided (1 point to the right, 7 to the left).

     i = Nr-1

     diff1(i) = idr*(1.d0/8.d0*diffvar(l,i+1) + 223.d0/140.d0*diffvar(l,i) - 3.5d0*diffvar(l,i-1) &
              + 3.5d0*diffvar(l,i-2) - 35.d0/12.d0*diffvar(l,i-3) + 1.75d0*diffvar(l,i-4) &
              - 0.7d0*diffvar(l,i-5) + 1.d0/6.d0*diffvar(l,i-6) - 1.d0/56.d0*diffvar(l,i-7))

!    Point Nr:  6th order fully one-sided.  The 8th order fully one-sided difference
!    (commented out) causes the runs to go unstable.

     i = Nr

     diff1(i) = idr*(49.d0/20.d0*diffvar(l,i) - 6.d0*diffvar(l,i-1) + 7.5d0*diffvar(l,i-2) &
              - 20.d0/3.d0*diffvar(l,i-3) + 15.d0/4.d0*diffvar(l,i-4) - 6.d0/5.d0*diffvar(l,i-5) &
              + 1.d0/6.d0*diffvar(l,i-6)) ! 6th order fully one-sided.

     !diff1(i) = idr*(761.d0/280.d0*diffvar(l,i) - 8.d0*diffvar(l,i-1) + 14.d0*diffvar(l,i-2) &
     !         - 56.d0/3.d0*diffvar(l,i-3) + 35.d0/2.d0*diffvar(l,i-4) - 56.d0/5.d0*diffvar(l,i-5) &
     !         + 14.d0/3.d0*diffvar(l,i-6) - 8.d0/7.d0*diffvar(l,i-7) + 1.d0/8.d0*diffvar(l,i-8)) ! 8th order

!    Check if we want one-sided derivatives at origin.

     if (oflag) then

        i = 1
        diff1(i) = - idr*(49.d0/20.d0*diffvar(l,i) - 6.d0*diffvar(l,i+1) + 7.5d0*diffvar(l,i+2) &
                 - 20.d0/3.d0*diffvar(l,i+3) + 15.d0/4.d0*diffvar(l,i+4) - 6.d0/5.d0*diffvar(l,i+5) &
                 + 1.d0/6.d0*diffvar(l,i+6)) ! 6th order fully one-sided.

        i = 2
        diff1(i) = - idr*(1.d0/8.d0*diffvar(l,i-1) + 223.d0/140.d0*diffvar(l,i) - 3.5d0*diffvar(l,i+1) &
                 + 3.5d0*diffvar(l,i+2) - 35.d0/12.d0*diffvar(l,i+3) + 1.75d0*diffvar(l,i+4) &
                 - 0.7d0*diffvar(l,i+5) + 1.d0/6.d0*diffvar(l,i+6) - 1.d0/56.d0*diffvar(l,i+7))

        i = 3
        diff1(i) = - idr*(-1.d0/56.d0*diffvar(l,i-2) + 2.d0/7.d0*diffvar(l,i-1) + 19.d0/20.d0*diffvar(l,i) &
                 - 2.d0*diffvar(l,i+1) + 1.25d0*diffvar(l,i+2) - 2.d0/3.d0*diffvar(l,i+3) &
                 + 0.25d0*diffvar(l,i+4) - 2.d0/35.d0*diffvar(l,i+5) + 1.d0/168.d0*diffvar(l,i+6))

        i = 4
        diff1(i) = - idr*(1.d0/168.d0*diffvar(l,i-3) - 1.d0/14.d0*diffvar(l,i-2) + 0.5d0*diffvar(l,i-1) &
                 + 9.d0/20.d0*diffvar(l,i) - 5.d0/4.d0*diffvar(l,i+1) + 0.5d0*diffvar(l,i+2) &
                 - 1.d0/6.d0*diffvar(l,i+3) + 1.d0/28.d0*diffvar(l,i+4) - 1.d0/280.d0*diffvar(l,i+5))

     end if

  end if


! ********************************
! ***   SYMMETRIES AT ORIGIN   ***
! ********************************

  if (rank==0) then
     do i=1,ghost
        diff1(1-i) = - dble(sym)*diff1(i)
     end do
  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

  if (size>1) then
     syncvar => diff1
     call sync
  end if


! ******************************
! ***   END FUNCTION DIFF1   ***
! ******************************

  end function diff1







  function diff2(l,sym,origin)

! ******************************************************
! ***   CALCULATE SECOND DERIVATIVE TO GIVEN ORDER   ***
! ******************************************************

  use procinfo
  use param
  use arrays

  implicit none

  logical oflag
  logical,optional :: origin

  integer i,l,sym

  real(8) idr,idr2

  real(8), dimension (1-ghost:Nrmax), target :: diff2


! *******************
! ***   NUMBERS   ***
! *******************

  idr  = 1.d0/dr(l)
  idr2 = idr**2


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

  diff2 = 0.d0


! *********************
! ***   2ND ORDER   ***
! *********************

  if (order=="two") then

!    Interior points:  2nd order centered second derivative.

     do i=1,Nr-1
        diff2(i) = idr2*(diffvar(l,i+1) - 2.d0*diffvar(l,i) + diffvar(l,i-1))
     end do

!    Point Nr:  2nd order fully one-sided.

     i = Nr

     diff2(i) = idr2*(2.d0*diffvar(l,i) - 5.d0*diffvar(l,i-1) &
               + 4.d0*diffvar(l,i-2) - diffvar(l,i-3))

!    Check if we want one-sided derivatives at origin.

     if (oflag) then
        diff2(1) = idr2*(2.d0*diffvar(l,1) - 5.d0*diffvar(l,2) &
                 + 4.d0*diffvar(l,3) - diffvar(l,4))
     end if


! *********************
! ***   4TH ORDER   ***
! *********************

  else if (order=="four") then

!    Interior points:  4th order centered second derivative.

     do i=1,Nr-2
        diff2(i) = - idr2*(30.d0*diffvar(l,i) &
                 - 16.d0*(diffvar(l,i+1) + diffvar(l,i-1)) &
                 +        diffvar(l,i+2) + diffvar(l,i-2))/12.d0
     end do

!    Point Nr-1:  4th order semi one-sided (1 point to the right, 4 to the left).

     i = Nr-1

     diff2(i) = idr2*(10.d0*diffvar(l,i+1) - 15.d0*diffvar(l,i) - 4.d0*diffvar(l,i-1) &
              + 14.d0*diffvar(l,i-2) - 6.d0*diffvar(l,i-3) + diffvar(l,i-4))/12.d0

!    Point Nr:  4th order fully one-sided.

     i = Nr

     diff2(i) = idr2*(45.d0*diffvar(l,i) - 154.d0*diffvar(l,i-1) + 214.d0*diffvar(l,i-2) &
              - 156.d0*diffvar(l,i-3) + 61.d0*diffvar(l,i-4) - 10.d0*diffvar(l,i-5))/12.d0

!    Check if we want one-sided derivatives at origin.

     if (oflag) then

        i = 1
        diff2(i) = idr2*(45.d0*diffvar(l,i) - 154.d0*diffvar(l,i+1) + 214.d0*diffvar(l,i+2) &
                 - 156.d0*diffvar(l,i+3) + 61.d0*diffvar(l,i+4) - 10.d0*diffvar(l,i+5))/12.d0

        i = 2
        diff2(i) = idr2*(10.d0*diffvar(l,i-1) - 15.d0*diffvar(l,i) - 4.d0*diffvar(l,i+1) &
                 + 14.d0*diffvar(l,i+2) - 6.d0*diffvar(l,i+3) + diffvar(l,i+4))/12.d0

     end if


! ***********************
! ***   6TH ORDER   ***
! ***********************

  else if (order=="six") then

!    Interior points:  6th order centered second derivative.

     do i=1,Nr-3
        diff2(i) = - idr2*(49.d0/18.d0*diffvar(l,i) &
                 -      1.5d0*(diffvar(l,i+1) + diffvar(l,i-1)) &
                 +     0.15d0*(diffvar(l,i+2) + diffvar(l,i-2)) &
                 - 1.d0/90.d0*(diffvar(l,i+3) + diffvar(l,i-3)))
     end do

!    Point Nr-2:  6th order semi-semi one-sided (2 points to the right, 5 to the left).

     i = Nr-2

     diff2(i) = idr2*(-11.d0/180.d0*diffvar(l,i+2) + 107.d0/90.d0*diffvar(l,i+1) - 2.1d0*diffvar(l,i) &
              + 13.d0/18.d0*diffvar(l,i-1) + 17.d0/36.d0*diffvar(l,i-2) - 0.3d0*diffvar(l,i-3) &
              + 4.d0/45.d0*diffvar(l,i-4) - 1.d0/90.d0*diffvar(l,i-5))

!    Point Nr-1:  4th order semi one-sided difference, the 6th order semi one-sided difference
!    (commented out) causes the runs to go unstable.  I also tried 5th order, but it too goes
!    ustable (thought not as badly).

     i = Nr-1

     diff2(i) = idr2*(10.d0*diffvar(l,i+1) - 15.d0*diffvar(l,i) - 4.d0*diffvar(l,i-1) &
              + 14.d0*diffvar(l,i-2) - 6.d0*diffvar(l,i-3) + diffvar(l,i-4))/12.d0 ! 4th order

     !diff2(i) = idr2*(0.7d0*diffvar(l,i+1) - 7.d0/18.d0*diffvar(l,i) - 2.7d0*diffvar(l,i-1) &
     !         + 19.d0/4.d0*diffvar(l,i-2) - 67.d0/18.d0*diffvar(l,i-3) + 9.d0/5.d0*diffvar(l,i-4) &
     !         - 0.5d0*diffvar(l,i-5) + 11.d0/180.d0*diffvar(l,i-6)) ! 6th order

!    Point Nr:  6th order fully one-sided.

     i = Nr

     diff2(i) = idr2*(469.d0/90.d0*diffvar(l,i) - 22.3d0*diffvar(l,i-1) + 879.d0/20.d0*diffvar(l,i-2) &
              - 949.d0/18.d0*diffvar(l,i-3) + 41.d0*diffvar(l,i-4) - 20.1d0*diffvar(l,i-5) &
              + 1019.d0/180.d0*diffvar(l,i-6) - 0.7d0*diffvar(l,i-7))

!    Check if we want one-sided derivatives at origin.

     if (oflag) then

        i = 1
        diff2(i) = idr2*(469.d0/90.d0*diffvar(l,i) - 22.3d0*diffvar(l,i+1) + 879.d0/20.d0*diffvar(l,i+2) &
                 - 949.d0/18.d0*diffvar(l,i+3) + 41.d0*diffvar(l,i+4) - 20.1d0*diffvar(l,i+5) &
                 + 1019.d0/180.d0*diffvar(l,i+6) - 0.7d0*diffvar(l,i+7))

        i = 2
        diff2(i) = idr2*(10.d0*diffvar(l,i-1) - 15.d0*diffvar(l,i) - 4.d0*diffvar(l,i+1) &
                 + 14.d0*diffvar(l,i+2) - 6.d0*diffvar(l,i+3) + diffvar(l,i+4))/12.d0 ! 4th order

        i = 3
        diff2(i) = idr2*(-11.d0/180.d0*diffvar(l,i-2) + 107.d0/90.d0*diffvar(l,i-1) - 2.1d0*diffvar(l,i) &
                 + 13.d0/18.d0*diffvar(l,i+1) + 17.d0/36.d0*diffvar(l,i+2) - 0.3d0*diffvar(l,i+3) &
                 + 4.d0/45.d0*diffvar(l,i+4) - 1.d0/90.d0*diffvar(l,i+5))

     end if


! *********************
! ***   8TH ORDER   ***
! *********************

  else if (order=="eight") then

!    Interior points:  8th order centered second derivative.

     do i=1,Nr-4
        diff2(i) = - idr2*(205.d0/72.d0*diffvar(l,i) &
                 -   8.d0/5.d0*(diffvar(l,i+1) + diffvar(l,i-1)) &
                 +       0.2d0*(diffvar(l,i+2) + diffvar(l,i-2)) &
                 - 8.d0/315.d0*(diffvar(l,i+3) + diffvar(l,i-3)) &
                 + 1.d0/560.d0*(diffvar(l,i+4) + diffvar(l,i-4)))
     end do

!    Point Nr-3:  8th order semi-semi-semi left-sided (3 points to the right, 6 to the left).

     i = Nr-3

     diff2(i) = idr2*(19.d0/2520.d0*diffvar(l,i+3) - 67.d0/560.d0*diffvar(l,i+2) + 97.d0/70.d0*diffvar(l,i+1) &
              - 89.d0/36.d0*diffvar(l,i) + 23.d0/20.d0*diffvar(l,i-1) + 7.d0/40.d0*diffvar(l,i-2) &
              - 17.d0/90.d0*diffvar(l,i-3) + 11.d0/140.d0*diffvar(l,i-4) - 1.d0/56.d0*diffvar(l,i-5) &
              + 1.d0/560.d0*diffvar(l,i-6))

!    Point Nr-2:  8th order semi-semi one-sided (2 points to the right, 7 to the left).

     i = Nr-2

     diff2(i) = idr2*(-223.d0/5040.d0*diffvar(l,i+2) + 293.d0/280.d0*diffvar(l,i+1) - 395.d0/252.d0*diffvar(l,i) &
              - 13.d0/30.d0*diffvar(l,i-1) + 83.d0/40.d0*diffvar(l,i-2) - 319.d0/180.d0*diffvar(l,i-3) &
              + 59.d0/60.d0*diffvar(l,i-4) - 5.d0/14.d0*diffvar(l,i-5) + 389.d0/5040.d0*diffvar(l,i-6) &
              - 19.d0/2520.d0*diffvar(l,i-7))

!    Point Nr-1:  4th order semi one-sided difference, the 6th and 8th order semi one-sided
!    differences (commented out) cause the runs to go unstable.

     i = Nr-1

     diff2(i) = idr2*(10.d0*diffvar(l,i+1) - 15.d0*diffvar(l,i) - 4.d0*diffvar(l,i-1) &
              + 14.d0*diffvar(l,i-2) - 6.d0*diffvar(l,i-3) + diffvar(l,i-4))/12.d0 ! 4th order

     !diff2(i) = idr2*(0.7d0*diffvar(l,i+1) - 7.d0/18.d0*diffvar(l,i) - 2.7d0*diffvar(l,i-1) &
     !         + 19.d0/4.d0*diffvar(l,i-2) - 67.d0/18.d0*diffvar(l,i-3) + 9.d0/5.d0*diffvar(l,i-4) &
     !         - 0.5d0*diffvar(l,i-5) + 11.d0/180.d0*diffvar(l,i-6)) ! 6th order

     !diff2(i) = idr2*(761.d0/1260.d0*diffvar(l,i+1) + 61.d0/144.d0*diffvar(l,i) - 201.d0/35.d0*diffvar(l,i-1) &
     !         + 341.d0/30.d0*diffvar(l,i-2) - 1163.d0/90.d0*diffvar(l,i-3) + 411.d0/40.d0*diffvar(l,i-4) &
     !         - 17.d0/3.d0*diffvar(l,i-5) + 1303/630.d0*diffvar(l,i-6) - 9.d0/20.d0*diffvar(l,i-7) &
     !         + 223.d0/5040.d0*diffvar(l,i-8)) ! 8th order

!    Point Nr:  8th order fully one-sided.

     i = Nr

     diff2(i) = idr2*(6515.d0/1008.d0*diffvar(l,i) - 4609.d0/140.d0*diffvar(l,i-1) + 5869.d0/70.d0*diffvar(l,i-2) &
              - 6289.d0/45.d0*diffvar(l,i-3) + 6499.d0/40.d0*diffvar(l,i-4) - 265.d0/2.d0*diffvar(l,i-5) &
              + 6709.d0/90.d0*diffvar(l,i-6) - 967.d0/35.d0*diffvar(l,i-7) + 3407.d0/560.d0*diffvar(l,i-8) &
              - 761.d0/1260.d0*diffvar(l,i-9))

!    Check if we want one-sided derivatives at origin.

     if (oflag) then

        i = 1
        diff2(i) = idr2*(6515.d0/1008.d0*diffvar(l,i) - 4609.d0/140.d0*diffvar(l,i+1) + 5869.d0/70.d0*diffvar(l,i+2) &
                 - 6289.d0/45.d0*diffvar(l,i+3) + 6499.d0/40.d0*diffvar(l,i+4) - 265.d0/2.d0*diffvar(l,i+5) &
                 + 6709.d0/90.d0*diffvar(l,i+6) - 967.d0/35.d0*diffvar(l,i+7) + 3407.d0/560.d0*diffvar(l,i+8) &
                 - 761.d0/1260.d0*diffvar(l,i+9))

        i = 2
        diff2(i) = idr2*(10.d0*diffvar(l,i-1) - 15.d0*diffvar(l,i) - 4.d0*diffvar(l,i+1) &
                 + 14.d0*diffvar(l,i+2) - 6.d0*diffvar(l,i+3) + diffvar(l,i+4))/12.d0 ! 4th order

        i = 3
        diff2(i) = idr2*(-223.d0/5040.d0*diffvar(l,i-2) + 293.d0/280.d0*diffvar(l,i-1) - 395.d0/252.d0*diffvar(l,i) &
                 - 13.d0/30.d0*diffvar(l,i+1) + 83.d0/40.d0*diffvar(l,i+2) - 319.d0/180.d0*diffvar(l,i+3) &
                 + 59.d0/60.d0*diffvar(l,i+4) - 5.d0/14.d0*diffvar(l,i+5) + 389.d0/5040.d0*diffvar(l,i+6) &
                 - 19.d0/2520.d0*diffvar(l,i+7))

        i = 4
        diff2(i) = idr2*(19.d0/2520.d0*diffvar(l,i-3) - 67.d0/560.d0*diffvar(l,i-2) + 97.d0/70.d0*diffvar(l,i-1) &
                 - 89.d0/36.d0*diffvar(l,i) + 23.d0/20.d0*diffvar(l,i+1) + 7.d0/40.d0*diffvar(l,i+2) &
                 - 17.d0/90.d0*diffvar(l,i+3) + 11.d0/140.d0*diffvar(l,i+4) - 1.d0/56.d0*diffvar(l,i+5) &
                 + 1.d0/560.d0*diffvar(l,i+6))

     end if

  end if


! ********************************
! ***   SYMMETRIES AT ORIGIN   ***
! ********************************

  if (rank==0) then
     do i=1,ghost
        diff2(1-i) = sym*diff2(i)
     end do
  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

  if (size>1) then
     syncvar => diff2
     call sync
  end if


! ******************************
! ***   END FUNCTION DIFF2   ***
! ******************************

  end function diff2







! **********************
! ***   END MODULE   ***
! **********************

  end module derivatives

