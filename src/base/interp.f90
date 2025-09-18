!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/interp.f90,v 1.7 2019/08/02 16:02:48 malcubi Exp $

  real(8) function interp(l,r0,errorflag)

! *************************
! ***   INTERPOLATION   ***
! *************************

! This function does linear or cubic interpolation.
! Notice that the array that is interpolated is always
! the pointer array "interpvar", so one needs to make sure
! that this points to the variable one really wants to
! interpolate before calling this function.
!
! On entry r0 is the point where we want to interpolate,
! and l is the grid level at which we want to interpolate.
!
! The third argument is a logical flag that, if true,
! causes the code to abort if the point r0 is outside
! the local grid (for multiprocessor runs), and if false 
! allows the code to continue and returns 0 as the
! "interpolated" value. 
!
! On output, 'interp' contains the interpolated value at r0.

! Include modules.

  use param
  use arrays
  use procinfo

! Declare variables
  
  implicit none

  logical errorflag       ! Flag to control the behaviour on error.

  integer l               ! Grid level to use for interpolation.
  integer i0              ! Grid point immediatly to the left.

  real(8) r0              ! Point at which we are interpolating.
  real(8) ddr             ! Current grid spacing.
  real(8) delta           ! Distance to grid point to the left.
  real(8) rmin,rmax       ! Minimum and maximum radii.
  real(8) aux

  real(8) :: small = 1.d-10

  real(8) fa(-1:2)        ! Small local array with values at nearest neighbors.


! ************************
! ***   SANITY CHECK   ***
! ************************

! First check that the interpolating point is inside the grid,
! and that there are enough points to either side for cubic
! interpolation.

  rmin = (dble(1-ghost)-0.5d0)*dr(l)
  rmax = (dble(Nrtotal)-0.5d0)*dr(l)

  if ((r0<rmin).or.(r0>rmax)) then
     if (rank==0) then
        print *
        print *, 'interp.f90: Interpolating point is outside the grid'
        write(*,'(A,ES12.4)') ' Interpolating position =',r0
        write(*,'(A,ES12.4)') ' Left boundary          =',rmin
        write(*,'(A,ES12.4)') ' Right boundary         =',rmax
        write(*,'(A,I0)') ' Grid level = ',l
        write(*,'(A,I0)') ' Processor  = ',rank
        print *
        print *, 'Aborting! (subroutine interp)'
        print *
     end if
     call die
  else if (r0<rmin+dr(l)) then
     if (rank==0) then
        print *
        print *, 'interp.f90: Not enough points to left side for cubic interpolation'
        write(*,'(A,ES12.4)') ' Interpolating position =',r0
        write(*,'(A,ES12.4)') ' Left boundary          =',rmin
        write(*,'(A,ES12.4)') ' Grid spacing           =',dr(l)
        write(*,'(A,I0)') ' Grid level = ',l
        print *
        print *, 'Aborting! (subroutine interp)'
        print *
     end if
     call die
  else if (r0>rmax-dr(l)) then
     if (rank==0) then
        print *
        print *, 'interp.f90: Not enough points to right side for cubic interpolation'
        write(*,'(A,ES12.4)') ' Interpolating position =',r0
        write(*,'(A,ES12.4)') ' Right boundary         =',rmax
        write(*,'(A,ES12.4)') ' Grid spacing           =',dr(l)
        write(*,'(A,I0)')     ' Grid level = ',l
        print *
        print *, 'Aborting! (subroutine interp)'
        print *
     end if
     call die
  end if


! **************************************
! ***   CHECK IF YOU OWN THE POINT   ***
! **************************************

! For multiple processor runs check that you
! own the point to be interpolated and that
! there are enough points for cubic interpolation.
!
! Notice that we need to make sure that only one
! processor owns the point, which might be a problem
! close to inter-processor boundaries where there is
! an overlap.  When that happens here I choose the
! processor to the left as owner of the point (but
! carefull with processor zero that has no neighbour
! to the left).

  if (size>1) then

!    Find out minimum and maximum radius owned by
!    a given processor.  Be careful for ghost>2.

     if (ghost<=2) then
        aux = 2.d0
     else
        aux = dble(ghost)
     end if

     if (rank==0) then
        rmin = rleft(rank,l) + dr(l) + small
     else
        rmin = rleft(rank,l) + aux*dr(l) + small
     end if

     rmax = rright(rank,l) - (aux-1.d0)*dr(l) + small

!    If you don't own the grid point then either die,
!    or return 0 depending on the value of 'errorflag'.

     if ((r0<=rmin).or.(r0>rmax)) then

        if (errorflag) then
            write(*,'(A,3ES12.4)') ' I do not own the point (r0,rmin,rmax):',r0,rmin,rmax
            write(*,'(A,I0,A,I0)') ' Processor: ',rank,'   Grid level: ',l
            print *, 'Aborting! (subroutine interp)'
            print *
            call die
        else
           interp = 0.d0
           return
        end if

     end if

  end if


! *********************************************
! ***   FIND GRID POINT TO THE LEFT OF r0   ***
! *********************************************

! Grid spacing.

  ddr = dr(l)

! Find grid point to the left of r0.

  i0 = int((r0-r(l,0))/ddr)
  delta = r0 - r(l,i0)

! Check that the point to be interpolated
! is indeed between grid points (i0,i0+1).

  if ((delta<-small).or.(delta>dr(l)+small)) then
     print *
     print *, 'This should never happen!'
     print *, r0,r(l,i0),r(l,i0+1),l,rank,delta,dr(l),i0,(r0-r(l,0))/ddr
     print *, 'Aborting! (subroutine interp)'
     print *
     call die
  end if


! ********************************************************
! ***   COPY THE VALUES AT NEAREST GRID POINTS TO fa   ***
! ********************************************************

! Copy the values of the array to be interpolated on the
! nearest grid points to either side of r0 to the small
! local array 'fa':
!
!    o        o        o        o
!   i0-1     i0       i0+1     i0+2
!  fa(-1)   fa(0 )   fa(+1)   fa(+2)

  fa(-1) = interpvar(l,i0-1)
  fa( 0) = interpvar(l,i0  )
  fa(+1) = interpvar(l,i0+1)
  fa(+2) = interpvar(l,i0+2)


! *******************************
! ***   LINEAR INTERPOLATON   ***
! *******************************

  if (intorder==1) then

     aux = delta/ddr

     interp = fa(0) + aux*(fa(1) - fa(0))


! ******************************
! ***   CUBIC INTERPOLATON   ***
! ******************************

  else if (intorder==3) then

     aux = delta/ddr

     interp = fa(0) - aux*(fa(2) - 6.d0*fa(1) + 3.d0*fa(0) + 2.d0*fa(-1))/6.d0 &
            + aux**2*(3.d0*fa(1) - 6.d0*fa(0) + 3.d0*fa(-1))/6.d0 &
            + aux**3*(fa(2) - 3.d0*fa(1) + 3.d0*fa(0) - fa(-1))/6.d0


! **************************
! ***    UNKNOWN ORDER   ***
! **************************

  else

     print *
     print *, 'The order of interpolation should be 1 or 3.'
     print *, 'Aborting! (subroutine interp)'
     print *
     call die

  end if


! ***************
! ***   END   ***
! ***************

  end function interp
