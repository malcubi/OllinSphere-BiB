!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/evolve.f90,v 1.54 2025/10/13 18:44:37 malcubi Exp $

  subroutine evolve

! **********************************
! ***   MAIN ITERATION ROUTINE   ***
! **********************************

! Include modules.

  use param
  use arrays
  use procinfo
  use derivatives

! Extra variables.

  implicit none

  logical contains

  integer i,h,l,p     ! Counters.

  real(8) zero,half   ! Numbers.


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  half  = 0.5d0


! *****************************************************
! ***   GRID SPACING AND TIME STEP FOR EACH LEVEL   ***
! *****************************************************

! Find grid spacing and time step for each level.

  dt0 = dtfac*dr0

  do l=0,Nl-1
     dr(l) = dr0/2**l
     dt(l) = dt0/2**l
     ! print *,l,dr(l),dt(l)
  end do


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

! Notice that we stagger the origin and take a number
! of ghost zones on the negative side given by the
! parameter "ghost".

! Find grid point positions.

  do l=0,Nl-1
     do i=1-ghost,Nrmax
        r(l,i) = (dble(Nmin(rank) + i) - half)*dr(l)
     end do
  end do

! Make sure that all processors know the position
! of the first and last grid point of every other
! processor at all grid levels.

  do p=0,size-1
     do l=0,Nl-1
        rleft(p,l)  = (dble(Nmin(p) + 1 - ghost) - half)*dr(l)
        rright(p,l) = (dble(Nmin(p) + Nrl(p)   ) - half)*dr(l)
     end do
  end do

! Find position of external boundary.

  rbound = rright(size-1,0)


! *****************************************************
! ***   SET TIME AND NUMBER OF TIME STEPS TO ZERO   ***
! *****************************************************

! Step counter and current time.

  s = 0
  t = 0.d0

! Old times.

  t1 = 0.d0
  t2 = 0.d0


! *****************************
! ***   FIND INITIAL DATA   ***
! *****************************

! Call initial data routine.

  call initial

! Make sure we have correct symmetries at the origin.
! Symmetries at the origin are only needed for processor 0
! who always owns the origin.

  if (rank==0) then
     do l=0,Nl-1
        call symmetries(l)
     end do
  end if

! If we have more than one processor synchronize ghost zones.

  if (size>1) then
     do l=0,Nl-1
        call syncall(l)
     end do
  end if


! ******************************************
! ***   CALCULATE AUXILIARY QUANTITIES   ***
! ******************************************

! Auxiliary quantities for matter and geometry.

  do l=0,Nl-1
     call auxiliary(l)
  end do


! ********************************************
! ***   MAXIMAL SLICING FOR INITIAL DATA   ***
! ********************************************

! For maximal slicing we have to solve for the lapse
! at t=0 also (must be done after call to stressenergy).

  if ((slicing=="maximal").or.(ilapse=="maximal")) then
     call alphamaximal(0,Nl-1,maximalbound,1.d0)
  end if


! ********************
! ***   ANALYSIS   ***
! ********************

! Calculate constraints.

  call constraints

! Analysis.

  call analysis_geometry

  if (mattertype/="vacuum") then
     call analysis_matter
  end if

! Find apparent horizon.

  if (ahfind.and.(t(0)>=ahafter)) then
     call horizon_finder
  end if


! *************************************
! ***   TRACK ANALYTIC SPACETIMES   ***
! *************************************

! Minkowski.

  if (TrackMinkowski) then

     if (idata=="minkowski") then

        do l=0,Nl-1
           call trackmink(l)
        end do

     else

        print *, 'Tracking Minkowski requires idata=minkowski ...'
        print *, 'Aborting! (subroutine base/evolve.f90)'
        print *
        call die

     end if

  end if


! ************************
! ***   FIND SOURCES   ***
! ************************

! Here we calculate the sources for all evolution
! equations once.  This is not really necessary,
! it is only done so that we can output them at
! t=0 if needed.

  if (spacetime=="dynamic") then

     do l=0,Nl-1

!       Sources for geometry.

        call sources_geometry(l)

!       Sources for lapse.

        call sources_lapse(l,dt(l))

!       Sources for shift.

        if (shift/="none") then
           call sources_shift(l)
        end if

     end do

  end if

! Sources for matter.

  if (mattertype/="vacuum") then
     do l=0,Nl-1
        call sources_matter(l)
     end do
  end if


! *********************************
! ***   SAVE THE INITIAL DATA   ***
! *********************************

  call save0Ddata
  call save1Ddata

! Checkpoint initial data, except when we are
! restarting from a checkpoint file (since we
! already have it).

  if (checkpointinitial.and.(idata/="checkpoint")) then
     call checkpointsave
  end if


! ****************************
! ***   OUTPUT TO SCREEN   ***
! ****************************

  if (rank==0) then
     print *
     print *,'------------------------------'
     print *,'|  Time step  |     Time     |'
     print *,'------------------------------'
     write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',0,'   | ',t(0),'  | '
  end if


! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

  do h=1,Nt

!    ***************************************
!    ***   DO WE ADJUST THE TIME STEP?   ***
!    ***************************************

!    If needed, adjust the time step according to
!    the values of the propagation speeds on the
!    base grid.

     if (adjuststep) then
        call cfl(h)
     end if


!    ****************************************
!    ***   ADVANCE ONE GLOBAL TIME STEP   ***
!    ****************************************

!    Advance one time step.  We only call the
!    time stepping routine for the coarse
!    grid, it is then called recursively from
!    inside for the other levels.

     call onestep(0)


!    ***************************
!    ***   MAXIMAL SLICING   ***
!    ***************************

!    If we are using maximal slicing, we must now
!    solve for the lapse. We also need to calculate
!    the derivatives of the lapse.

     if (slicing=="maximal") then
        call alphamaximal(0,Nl-1,maximalbound,1.d0)
     end if


!    ********************
!    ***   ANALYSIS   ***
!    ********************

!    Calculate constraints.

     call constraints

!    Analysis.

     call analysis_geometry

     if (mattertype /= "vacuum") then
        call analysis_matter
     end if


!    ********************
!    ***   HORIZONS   ***
!    ********************

!    We only look for apparent horizons every "ahfind_every"
!    time steps.  Notice that even if ahfind_every=1, we
!    only look for horizons at the coarse time steps.

     if (ahfind.and.(mod(h,ahfind_every).eq.0)) then
        call horizon_finder
     end if


!    ******************************
!    ***   SAVE DATA TO FILES   ***
!    ******************************

!    Save 0D data (only every Noutput0D time steps).

     if (Noutput0D>0) then
        if (mod(h,Noutput0D)==0) then
           call save0Ddata
        end if
     end if

!    Save 1D data (only every Noutput1D time steps).

     if (Noutput1D>0) then
        if (mod(h,Noutput1D)==0) then
           call save1Ddata
        end if
     end if

!    Save checkpoint data.

     if (checkpoint.and.(mod(h,Ncheckpoint)==0)) then
        call checkpointsave
     end if


!    ***********************************
!    ***   END MAIN EVOLUTION LOOP   ***
!    ***********************************

!    Time step information to screen.

     if (rank==0) then
        if (mod(h,Ninfo).eq.0) then
           write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',h,'   | ',t(0),'  | '
        end if
     end if

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine evolve
