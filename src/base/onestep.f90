
  recursive subroutine onestep(l)

! *********************************
! ***   ADVANCE ONE TIME STEP   ***
! *********************************

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use derivatives

! Extra variables.

  implicit none

  logical contains

  integer i,imax           ! Grid point counter.
  integer l                ! Refinement level counter.
  integer k                ! Counter for internal iterations.
  integer niter            ! Number of internal iterations.
  integer border           ! Order of interpolation at fine boundary.

  real(8) dtw              ! Internal time step.
  real(8) tp,t0,tm1,tm2,tl ! Local time at different time steps (for interpolation).
  real(8) weight           ! Weight for rk4.
  real(8) r0               ! Interpolation point.
  real(8) interp           ! Interpolation function.
  real(8) aux1,aux2,aux3


! ******************************
! ***   SAVE OLD TIME STEP   ***
! ******************************

  call saveold(l)


! **********************************************
! ***   FIND NUMBER OF INTERNAL ITERATIONS   ***
! **********************************************

  if (integrator=="icn") then
     niter = icniter
  else if (integrator=="rk4") then
     niter = 4
  end if


! **************************************
! ***   ADVANCE ONE FULL TIME STEP   ***
! **************************************

  do k=1,niter

!    ************************
!    ***   FIND WEIGHTS   ***
!    ************************

!    Find out weights for each iteration for the
!    different time integration schemes.

     if (integrator=="icn") then

!       In normal ICN, all iterations except the last one
!       jump only half a time step.

        if (k<niter) then
           dtw = 0.5d0*dt(l)
        else
           dtw = dt(l)
        end if

!    Fourth order Runge-Kutta.

     else if (integrator=="rk4") then

!       In fourth order Runge-Kutta the first two iterations
!       jump half a time step and the last two a full time step.
!       Here we also set the weights with which intermediate
!       results contribute to final answer: 1/6 for first and
!       last intermediate results and 1/3 for the two middle ones.

        if (k==1) then
           dtw = 0.5d0*dt(l)
           weight = 1.d0/6.d0
        else if (k==2) then
           dtw = 0.5d0*dt(l)
           weight = 1.d0/3.d0
        else if (k==3) then
           dtw = dt(l)
           weight = 1.d0/3.d0
        else
           dtw = dt(l)
           weight = 1.d0/6.d0
        end if

     end if


!    *********************************
!    ***   SOURCES FOR SPACETIME   ***
!    *********************************

     if (spacetime=="dynamic") then

!       Sources for geometry.

        call sources_geometry(l)

!       Sources for lapse. The sources for the lapse should
!       be calculated after the sources for the geometry,
!       since some lapse conditions might require strK.

        call sources_lapse(l,dtw)

!       Sources for shift. The sources for the shift should
!       be calculated after the sources for the geometry
!       and the sources for the lapse, since the shift
!       conditions use sDeltar and some use salpha.

        if (shift/="none") then
           call sources_shift(l)
        end if

     end if


!    ******************************
!    ***   SOURCES FOR MATTER   ***
!    ******************************

     if (mattertype/="vacuum") then
        call sources_matter(l)
     end if


!    *******************************
!    ***   BOUNDARY CONDITIONS   ***
!    *******************************

!    Outer boundaries are only needed for level 0
!    and the last processor.  They are applied directly
!    to the sources before updating the variables.

     if ((l==0).and.(rank==size-1)) then

!       Static and flat boundaries.  The routine "simpleboundary" applies
!       either static or flat boundaries to the sources of all those arrays
!       that are not declared with the attribute NOBOUND.
!
!       These type of boundary conditions are really just for testing,
!       but "flat" boundaries can in fact be useful for cosmological
!       spacetimes, as long they are are sufficiently far away.

        if ((boundtype=="static").or.(boundtype=="flat")) then

           call simpleboundary(l)

!       Radiative boundaries for geometry. The routine "radiative_geometry"
!       applies outgoing wave boundary conditions to the geometric
!       variables (trK,KTA,Klambda,Deltar).  These boundary conditions take
!       into account the characteristic structure, but not the constraints.
!       Notice that matter variables should define their own radiative
!       boundaries in their source routines.

        else if (boundtype=="radiative") then

           call radiative_geometry(l)

!       Constraint preserving boundary conditions. In this case we first
!       apply the standard radiative boundaries for (trK,Deltar), and
!       after that we apply the constraint preserving boundary conditions
!       for (KTA,Klambda).

        else if (boundtype=="constraint") then

           call radiative_geometry(l)
           call constraintbound(l)

        end if

     end if


!    *****************************************************
!    ***   FOR RUNGE-KUTTA ADD TO ACCUMULATOR ARRAYS   ***
!    *****************************************************

!    The accumulator arrays add the contributions from the
!    different Runge-Kutta iterations with their corresponding
!    weights.  The routine "accumulate" (generated by perl)
!    works in the following way:
!
!    1) On the first step (k=1) it just copies the source array times
!       its weight to the accumulator array.
!    2) On all subsequent steps except the last (1<k<niter), the routine
!       adds the source times its weight to the accumulator array.
!    3) On the last call (k=niter), it adds the final source times its
!       weight, but stores the result back into the source arrays so
!       that the last update will work correctly.

     if (integrator=="rk4") then
        call accumulate(l,k,niter,weight)
     end if


!    ****************************
!    ***   UPDATE VARIABLES   ***
!    ****************************

     call update(l,dtw)


!    *************************************************
!    ***   FOR FINE GRIDS INTERPOLATE BOUNDARIES   ***
!    *************************************************

!    For fine grids we need to interpolate from the new
!    time level of the coarse grid to get boundary data.
!
!    Remember that the coarse grid has already advanced
!    to the next time level.

     if (l>0) then

!       Figure out times for interpolation:
!
!       tp  = Local time on coarse grid (which has already advanced).
!       t0  = Local time on fine grid.
!       tm1 = Time 1 step  to the past for fine grid.
!       tm2 = Time 2 steps to the past for fine grid.
!       tl  = t0 + dtw (time to which we want to interpolate).
!
!       Notice that I am not assuming here that the time steps are uniform.

        tp  = t(l-1)
        t0  = t(l  )

        tm1 = t1(l)
        tm2 = t2(l)

        tl = t0 + dtw

!       Order of interpolation at boundary.

        if (t1(l)<=dt(l)) then
           border = 1
        else
           border = 2
        end if

!       Apply boundary condition. Since here I don't know how many
!       variables are evolving, I generate the necessary code
!       automatically at compile time and include it here.
!       The basic block has the form:
!
!       interpvar => array
!       aux1 = interp(l-1,r0,.false.)
!       call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
!       if (border==1) then
!          array(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*array_p(l,Nr-i)
!       else
!          array(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &
!                   + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*array_bound(l,i,1) &
!                   + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*array_bound(l,i,2)
!       end if
!
!       The interpolation above is done to second order, except
!       for the first two timesteps for which we don't have enough
!       points to the past.
!
!       Notice that for multi-processor runs all processors will
!       interpolate to the real boundary of the fine grid (aux1)
!       (processors that don't own that point will just return 0).
!       The MPI_ALLREDUCE call above makes all processors end up
!       with the same boundary value (aux2), but this doesn't matter
!       as it will be fixed when we synchronize below.

        if (tp==tl) then
           imax = ghost+2
        else
           imax = ghost-1
        end if

        do i=0,imax
           r0 = (dble(Nrtotal-i)-0.5d0)*dr(l)
           include '../auto/boundinterp.inc'
        end do

     end if


!    **********************
!    ***   SYMMETRIES   ***
!    **********************

!    Apply symmetries at origin.

     if (rank==0) then
        call symmetries(l)
     end if


!    ***********************
!    ***   SYNCHRONIZE   ***
!    ***********************

!    If we have more than one processor we must now
!    synchronize ghost zones.

     if (size>1) then
        call syncall(l)
     end if


!    *************************************************************
!    ***   UPDATE AUXILIARY VARIABLES FOR MATTER AND GEOMETRY  ***
!    *************************************************************

!    Auxiliary quantities for matter and geometry.

     call auxiliary(l)


!    ***********************************
!    ***   END INTERNAL ITERATIONS   ***
!    ***********************************

  end do


! ****************************************************
! ***   ADVANCE LOCAL TIME AND TIME STEP COUNTER   ***
! ****************************************************

! Save old local times.

  t2(l) = t1(l)
  t1(l) = t(l)

! Advance time step counter and local time.

  s(l) = s(l) + 1
  t(l) = t(l) + dt(l)


! *************************************
! ***   TRACK ANALYTIC SPACETIMES   ***
! *************************************

! Minkowski.

  if (TrackMinkowski) then
     call trackmink(l)
  end if

! Schwarzschild.

  if (TrackSchwarzschild) then
     call trackschwarz(l)
  end if


! **********************************
! ***   ARE THERE FINER GRIDS?   ***
! **********************************

! If there is a finer grid we need to advance it twice
! to catch up.  Notice that here I am calling the
! current subroutine "onestep" recursively.

  if (l<Nl-1) then
     call onestep(l+1)
     call onestep(l+1)
  end if


! ****************************************************
! ***   RESTRICT FINE GRID DATA INTO COARSE GRID   ***
! ****************************************************

! Restrict the data from the fine to the coarse grid
! when both levels coincide in time.
! This restriction does not change data in the current
! grid level, but rather in the coarser level.
!
! Remember to fix the symmetries at the origin for level (l-1)
! and synchronize again for multi-processor runs!

  if ((l>0).and.(mod(s(l),2)==0)) then

!    Restrict.

     call restrict(l,.true.)

!    Symmetries.

     if (rank==0) then
        call symmetries(l-1)
     end if

!    Sync.

     if (size>1) then
        call syncall(l-1)
     end if

!    Restrict background cosmological variables.

     if (cosmic_run) then
        include '../auto/restrict_cosmo.inc'
     end if

!    Auxiliary quantities for matter and geometry on level (l-1).

     call auxiliary(l-1)

  end if


! **************************
! ***   SPECIAL OUTPUT   ***
! **************************

! Here we do output of ALL time steps for the current level,
! and not only at the coarse level. Notice that the output times
! will then not coincide, as we will have more output for fine
! levels than for coarse levels.
!
! This is really intended only for testing, and is done when the
! corresponding parameter Noutput is equal to 0.

  if (Noutput0D==0) then
     do i=1,nvars0D
        call grabarray(trim(outvars0Darray(i)))
        call save0Dvariable(directory,trim(outvars0Darray(i)),i,s(l),t(l),'old',l)
     end do
  end if

  if (Noutput1D==0) then
     do i=1,nvars1D
        call grabarray(trim(outvars1Darray(i)))
        call save1Dvariable(directory,trim(outvars1Darray(i)),i,s(l),'old',l)
     end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine onestep

