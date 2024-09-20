!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/TrackMink.f90,v 1.1 2024/09/03 17:12:30 malcubi Exp $

  subroutine trackmink(l)

! ****************************************
! ***   TRACKING MINKOWSKI SPACETIME   ***
! ****************************************

! This subroutine tracks the position of the coordinate lines
! for Minkowski spacetime in the standard Minkowski coordinates.
!
! We assume that the initial data is Minkowski and the initial
! slice is a usual flat Minkowski slice.  But we won't stay
! in Minkowski coordinates because of the non-trivial slicing
! and/or shift conditions.
!
! Notice that since we interpolate the metric functions only
! on two time slices, the resulting solutions are only second
! order accurate even if the rest of the code is 4th order.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  logical firstcall(0:Nl-1)        ! Is this the first call?
  logical error                    ! Error in finding solution.

  integer i,l,p                    ! Counters
  integer imax,Naux                ! Auxiliary.
  integer status(MPI_STATUS_SIZE)

  real(8) R_guess,R_left,R_right   ! R positions.
  real(8) T_guess,T_left,T_right   ! T positions.
  real(8) alpha_avg,beta_avg       ! Average lapse and shift.
  real(8) phi_avg,A_avg            ! Average metric components.
  real(8) gammarr                  ! Physical radial metric: grr = exp(4*phi)*A
  real(8) dleft,dright             ! Left and right physical intervals.

  real(8), dimension (1-ghost:Nrmax) :: var1,var2,var3  !  Arrays for output.

  character(5)  filen
  character(20) filestatus


! **********************
! ***   INITIALIZE   ***
! **********************

! Initialize error flag.

  error = .false.

! Initialize output files flags.

  firstcall(l) = .false.

! At t=0 initialize R_MINK and T_MINK to (r,0).  Remember that
! we assume we start with a standard flat Minkowski slice.

  if (t(l)==0.d0) then
     firstcall(l) = .true.
     R_MINK(l,:) = r(l,:)
     T_MINK(l,:) = 0.d0
     goto 100
  end if


! ***********************************************
! ***   UPDATE POSITION OF COORDINATE LINES   ***
! ***********************************************

! To update the position of the coordinate lines I first calculate
! the interval of the point we are considering to the points to
! the left and right onn the previous time step (forming a triangle):
!
!                     ...   o  o  x  o  o  ...
!                     ...   o  x  o  x  o  ...

! The idea is then to find that point in Minkowski spacetime that is
! at precisely those distances from the previous points whose positions
! are known.  Notice that one would expect two possible solutions,
! one in the past and one on the future, but I lock onto the correct
! one by choosing a initial guess in the future.

! Copy old positions.

  R_MINK_P = R_MINK
  T_MINK_P = T_MINK

! Loop over grid points.

  do i=1,Nr

!    As initial guess I assume that the observer did not move in space
!    and advanced one coarse time step in time.
 
     R_guess = R_MINK(l,i)
     T_guess = T_MINK(l,i) + dt(0)

!    Calculate interval to point on the left on previous time step.
!    First we average the lapse, shift and metric components, and then
!    we find the interval.  Notice the sign of the mixed term, it is
!    negative since seen from the current point, the point in the
!    last time level is in the past (-dt) and to the left (-dr).
!    Also remember that beta is the contravariant shift, and 
!    the interval requires the covariant shift, so we need to
!    multiply with the radial metric.

     alpha_avg = 0.5d0*(alpha(l,i) + alpha_p(l,i-1))
     phi_avg   = 0.5d0*(phi(l,i) + phi_p(l,i-1))
     A_avg     = 0.5d0*(A(l,i) + A_p(l,i-1))

     gammarr = exp(4.d0*phi_avg)*A_avg

     dleft = - (alpha_avg*dt(l))**2 + gammarr*dr(l)**2

     if (shift/="none") then
        beta_avg = 0.5d0*(beta(l,i) + beta_p(l,i-1))
        dleft = dleft + gammarr*(beta_avg*dt(l))**2 + 2.d0*gammarr*beta_avg*dr(l)*dt(l)
     end if

!    Calculate interval to point on the right on previous time step.
!    First we average the lapse, shift and metric components, and then
!    we find the interval.  Notice the sign of the mixed term, it is
!    negative since seen from the current point, the point in the
!    last time level is in the past (-dt) and to the right (+dr).
!    Also remember that beta is the contravariant shift, and 
!    the interval requires the covariant shift, so we need to
!    multiply with the radial metric.

     if (i<Nr) then

        alpha_avg = 0.5d0*(alpha(l,i) + alpha_p(l,i+1))
        phi_avg   = 0.5d0*(phi(l,i) + phi_p(l,i+1))
        A_avg     = 0.5d0*(A(l,i) + A_p(l,i+1))

        gammarr = exp(4.d0*phi_avg)*A_avg

        dright = - (alpha_avg*dt(l))**2 + gammarr*dr(l)**2

        if (shift/="none") then
           beta_avg = 0.5d0*(beta(l,i) + beta_p(l,i+1))
           dright = dright + gammarr*(beta_avg*dt(l))**2 - 2.d0*gammarr*beta_avg*dr(l)*dt(l)
        end if

     else

        alpha_avg = 0.5d0*(alpha(l,Nr) + alpha_p(l,Nr))
        phi_avg   = 0.5d0*(phi(l,Nr) + phi_p(l,Nr))
        A_avg     = 0.5d0*(A(l,Nr) + A_p(l,Nr))

        dright = - (alpha_avg*dt(l))**2

       if (shift/="none") then
           beta_avg = 0.5d0*(beta(l,Nr) + beta_p(l,Nr))
           dright = dright + gammarr*(beta_avg*dt(l))**2
        end if

     end if

!    Now solve for the position of the point that has precisely those
!    intervals in Minkowski coordinates.

     R_left  = R_MINK_P(l,i-1)
     T_left  = T_MINK_P(l,i-1)

     if (i<Nr) then
        R_right = R_MINK_P(l,i+1)
        T_right = T_MINK_P(l,i+1)
     else
        R_right = R_MINK_P(l,Nr)
        T_right = T_MINK_P(l,Nr)
     end if

     call solvepointmink(R_guess,T_guess,R_left,T_left,R_right,T_right,dleft,dright,error)

!    Error message.

     if (error) then
          print *, 'Newton-Raphson iterations did not converge at r =',r(l,i), ' level l =',l
          print *, 'Aborting! (subroutine geometry/trackmink.f90)'
          print *
          call die
     end if

!    Save new coordinates.

     R_MINK(l,i) = R_guess
     T_MINK(l,i) = T_guess

  end do

! Ghost zones (only processor 0 owns the origin).

  if (rank==0) then
     do i=1,ghost
        R_MINK(l,1-i) = - R_MINK(l,i)
        T_MINK(l,1-i) = + T_MINK(l,i)
     end do
  end if

! For parallel runs synchronize across processors.

  if (size>1) then
     syncvar => R_MINK(l,:)
     call sync
     syncvar => T_MINK(l,:)
     call sync
  end if

! Restrict fine data to coarse grid.

  if ((l>0).and.(mod(s(l),2)==0)) then

!    Restrict.

     restrictvar => R_MINK
     call restrict(l,.false.)

     restrictvar => T_MINK
     call restrict(l,.false.)

!    Fix ghost zones

     if (rank==0) then
        do i=1,ghost
           R_MINK(l,1-i) = - R_MINK(l,i)
           T_MINK(l,1-i) = + T_MINK(l,i)
        end do
     end if

!    Sync.

     if (size>1) then
        syncvar => R_MINK(l,:)
        call sync
        syncvar => T_MINK(l,:)
        call sync
     end if

  end if

! Now subtract r and t from their respective arrays.
! This is in order to see more clearly the details.

  R_MINK_S(l,:) = R_MINK(l,:) - r(l,:)
  T_MINK_S(l,:) = T_MINK(l,:) - t(l)


! ***********************
! ***   OUTPUT DATA   ***
! ***********************

  100 continue

! Do output only every Noutput1D time steps.
! Notice that for refinement levels we need
! to multipy the Noutput with 2**l to make
! sure all level do output at the same times.

  if (mod(s(l),Noutput1D*(2**l))/=0) return

! Find total number of points.

  Naux = Nrmax + ghost

! On first call, replace file. Otherwise just append to it.

  if (firstcall(l)) then
     firstcall(l) = .false.
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Processor 0 does output.

  if (rank==0) then

     if (l<10) then
        write(filen,'(i1)') l
     else
        write(filen,'(i2)') l
     end if

!    Open files.

     if (filestatus == 'replace') then
        open(1,file=trim(directory)//'/MinkSliceT'//trim(filen)//'.rl',form='formatted', &
             status=filestatus)
        open(2,file=trim(directory)//'/MinkSliceS'//trim(filen)//'.rl',form='formatted', &
             status=filestatus)
     else
        open(1,file=trim(directory)//'/MinkSliceT'//trim(filen)//'.rl',form='formatted', &
                status=filestatus,position='append')
        open(2,file=trim(directory)//'/MinkSliceS'//trim(filen)//'.rl',form='formatted', &
                status=filestatus,position='append')
     end if

!    Write current time.

     if (commenttype=='xgraph') then
        write(1,"(A8,ES18.10)") '"Time = ',t(l)
        write(2,"(A8,ES18.10)") '"Time = ',t(l)
     else
        write(1,"(A8,ES18.10)") '#Time = ',t(l)
        write(2,"(A8,ES18.10)") '#Time = ',t(l)
     end if

!    Save data from processor 0.

     if (size==1) then
        imax = Nrl(0)
     else
        imax = Nrl(0) - ghost
     end if

     do i=1-ghost,imax
        write(1,"(2ES16.8)") R_MINK(l,i),T_MINK(l,i)
        write(2,"(2ES16.8)") R_MINK(l,i),T_MINK_S(l,i)
     end do

!    Iterate over other processors.

     do p=1,size-1

!       Receive and save data from other processors.

        call MPI_RECV(var1,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(var2,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(var3,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

        if (p==size-1) then
           imax = Nrl(p)
        else
           imax = Nrl(p) - ghost
        end if

        do i=1,imax
           write(1,"(2ES16.8)") var1(i),var2(i)
           write(2,"(2ES16.8)") var1(i),var3(i)
        end do

     end do

!    Leave two blank spaces before next time.
!    The reason to leave two spaces is that 'gnuplot' asks
!    for two spaces to distinguish different records.

     write(1,*)
     write(1,*)

     write(2,*)
     write(2,*)

!    Close files.

     close(1)
     close(2)

! Other processors send data to processor 0.

  else

     call MPI_SEND(R_MINK(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
     call MPI_SEND(T_MINK(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
     call MPI_SEND(T_MINK_S(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine trackmink









  subroutine solvepointmink(r,t,r1,t1,r2,t2,d1,d2,error)

! ***************************
! ***   SOLVE FOR POINT   ***
! ***************************

! This subroutine solves for the position in Minkowski
! coordinates of a point (r,t) that has an interval d1
! with respect to the point (r1,t1) and an interval d2
! with respect to the point (r2,t2).
!
! The routine does that by using Newton-Raphson.  Notice
! that on entry (r,t) contains the initial guess and on
! exit the solution.

  implicit none

  logical error

  integer k,Nmax

  real(8) r,r1,r2
  real(8) t,t1,t2
  real(8) d1,d2
  real(8) det
  real(8) epsilon

  real(8) res(1:2),delta(1:2),jac(1:2,1:2)

! The variables introduced above are:
!
! k          Counter.
! Nmax       Maximum number of iterations.
!
! (r,t)      Solution (and initial guess).
!
! (r1,t1)    Coordinates of point 1.
! (r2,t2)    Coordinates of point 2. 
!
! d1         Distance to point 1.
! d2         Distance to point 2.
!
! res        Vector of residuals.
! delta      Vector of increments.
! jac        Jacobian matrix.
!
! det        Determinant of Jacobian.
!
! epsilon    Toletance.


! **************************
! ***   NEWTON-RAPHSON   ***
! **************************

! Start iterations.

  Nmax = 1000

  do k=1,Nmax

!    Find residuals using Minkowski metric in standard coordinates.
!    Notice that we do not need the angular coordinates.

     res(1) = (r - r1)**2 - (t - t1)**2 - d1
     res(2) = (r - r2)**2 - (t - t2)**2 - d2

!    Find Jacobian matrix and its determinant.

     jac(1,1) = + 2.0d0*(r - r1)
     jac(1,2) = - 2.0d0*(t - t1)

     jac(2,1) = + 2.0d0*(r - r2)
     jac(2,2) = - 2.0d0*(t - t2)

     det = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)

!    Solve for increments (J delta = - res).

     delta(1) = (res(2)*jac(1,2) - res(1)*jac(2,2))/det
     delta(2) = (res(1)*jac(2,1) - res(2)*jac(1,1))/det

!    Correct solution.

     r = r + delta(1)
     t = t + delta(2)

!    Check if we achieved desired tolerance.

     epsilon = dabs(delta(1)) + dabs(delta(2))

     if (epsilon < 1.0d-10) return

  end do

! If we get here we the algorithm did not converge, so we
! set the error flag to .true.

  error = .true.


! ***************
! ***   END   ***
! ***************

  end subroutine solvepointmink
