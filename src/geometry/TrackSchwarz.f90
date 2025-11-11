
  subroutine trackschwarz(l)

! ********************************************
! ***   TRACKING SCHWARZSCHILD SPACETIME   ***
! ********************************************

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

  logical firstcall(0:Nl-1)              ! Is this the first call?
  logical error                          ! Error in finding solution.

  integer i,l,p                          ! Counters
  integer imin,imax,Naux                 ! Auxiliary.
  integer status(MPI_STATUS_SIZE)

  real(8) CSI_guess,CSI_left,CSI_right   ! R positions.
  real(8) ETA_guess,ETA_left,ETA_right   ! T positions.
  real(8) alpha_avg,beta_avg             ! Average lapse and shift.
  real(8) phi_avg,A_avg,B_avg            ! Average metric components.
  real(8) ra_left,ra_right               ! Average areal radius to left and right.
  real(8) gammarr                        ! Physical radial metric: grr = exp(4*phi)*A
  real(8) dleft,dright                   ! Left and right physical intervals.

  real(8), dimension (1-ghost:Nrmax) :: var1,var2  !  Arrays for output.

  character(5)  filen
  character(20) filestatus


! **********************
! ***   INITIALIZE   ***
! **********************

! Initialize output files flags.

  firstcall(l) = .false.

! At t=0 initialize CSI_SCHWARZ and ETA_MINK  Remember that we
! assume that we start with Schwarzschild in isotrtopic coordinates.

  if (t(l)==0.d0) then

!    Set firstcall flag to true.

     firstcall(l) = .true.

!    Find Schwarzschild areal radius.

     r_area(l,:) = r(l,:)*psi(l,:)**2*sqrt(B(l,:))

!    Set Kruskal ETA to zero.

     ETA_SCHWARZ(l,:) = 0.d0

!    Find Kruskal CSI. At t=0 the transformation of
!    coordinates is, for r_area > 2M:
!
!                           1/2
!    CSI  =  (r_area/2M - 1)    exp(r_area/4M)
!
!    where r_area is the usual Schwarzschild areal
!    radius calculated above.  The expression above
!    is valid outside the horizon, which in the isotropic
!    coordinates of our initial data corresponds to r>M/2.
!
!    Remember also that r>M/2 is not the interior of
!    the black hole, but rather the other side of the
!    Einstein-Rosen bridge, so in that case we just
!    change the sign of CSI.
!
!    Notice that we ignore the ghost point to the left
!    of the origin since they make no sense here.
!    Also, for very small r we fiund very large CSI
!    because as r goes to zero we are in fact approaching
!   infinity on the other side of the wormhole.

     CSI_SCHWARZ(l,:) = 0.d0

     do i=1,Nr
        if (r(l,i)>=0.5d0*BHmass) then
           CSI_SCHWARZ(l,i) = + dsqrt(r_area(l,i)/2.d0/BHmass - 1.d0) &
                                *dexp(r_area(l,i)/4.d0/BHmass)
        else
           CSI_SCHWARZ(l,i) = - dsqrt(r_area(l,i)/2.d0/BHmass - 1.d0) &
                                *dexp(r_area(l,i)/4.d0/BHmass)
        end if
     end do

!    Jump to output.

     goto 100

  end if


! ***********************************************
! ***   UPDATE POSITION OF COORDINATE LINES   ***
! ***********************************************

! To update the position of the coordinate lines I first calculate
! the interval of the point we are considering to the points to
! the left and right on the previous time step (forming a triangle):
!
!                     ...   o  o  x  o  o  ...
!                     ...   o  x  o  x  o  ...

! The idea is then to find that point in Kruskal coordinates that is
! at precisely those distances from the previous points whose positions
! are known.  Notice that one would expect two possible solutions,
! one in the past and one on the future, but I lock onto the correct
! one by choosing a initial guess in the future.
!
! At the last grid point we modify the molecule above slighlty to:
!
!                     ...   o  o  x
!                     ...   o  x  x

! Copy old positions.

  CSI_SCHWARZ_p(l,:) = CSI_SCHWARZ(l,:)
  ETA_SCHWARZ_p(l,:) = ETA_SCHWARZ(l,:)

! Loop over grid points.

  imin = 30

  do i=imin,Nr

!    As initial guess I assume that the coordinate lines did not move
!    in space and advanced one time step in time.
 
     CSI_guess = CSI_SCHWARZ_P(l,i)
     ETA_guess = ETA_SCHWARZ_P(l,i) + dt(l)

!    Calculate interval to point on the left on previous time step.
!    First we average the lapse, shift and metric components, and then
!    we find the interval. Notice the sign of the mixed term, it is
!    negative since seen from the current point, the point in the
!    last time level is in the past (-dt) and to the left (-dr).
!    Also remember that beta is the contravariant shift, and 
!    the interval requires the covariant shift, so we need to
!    multiply with the radial metric.

     if (i==imin) then

        alpha_avg = 0.5d0*(alpha(l,imin) + alpha_p(l,imin))
        phi_avg   = 0.5d0*(  phi(l,imin) +   phi_p(l,imin))
        A_avg     = 0.5d0*(    A(l,imin) +     A_p(l,imin))
        B_avg     = 0.5d0*(    B(l,imin) +     B_p(l,imin))

        ra_left = r(l,imin)*dexp(2.d0*phi_avg)*dsqrt(B_avg)

        gammarr = exp(4.d0*phi_avg)*A_avg

        dleft = - (alpha_avg*dt(l))**2

        if (shift/="none") then
           beta_avg = 0.5d0*(beta(l,i) + beta_p(l,imin))
           dleft = dleft + gammarr*(beta_avg*dt(l))**2 + 2.d0*gammarr*beta_avg*dr(l)*dt(l)
        end if

        CSI_left = CSI_SCHWARZ_P(l,imin)
        ETA_left = ETA_SCHWARZ_P(l,imin)

     else

        alpha_avg = 0.5d0*(alpha(l,i) + alpha_p(l,i-1))
        phi_avg   = 0.5d0*(  phi(l,i) +   phi_p(l,i-1))
        A_avg     = 0.5d0*(    A(l,i) +     A_p(l,i-1))
        B_avg     = 0.5d0*(    B(l,i) +     B_p(l,i-1))

        ra_left = 0.5d0*(r(l,i)+r(l,i-1))*dexp(2.d0*phi_avg)*dsqrt(B_avg)

        gammarr = exp(4.d0*phi_avg)*A_avg

        dleft = - (alpha_avg*dt(l))**2 + gammarr*dr(l)**2

        CSI_left = CSI_SCHWARZ_P(l,i-1)
        ETA_left = ETA_SCHWARZ_P(l,i-1)

     end if

!    Calculate interval to point on the right on previous time step.
!    First we average the lapse, shift and metric components, and then
!    we find the interval.  Notice the sign of the mixed term, it is
!    negative since seen from the current point, the point in the
!    last time level is in the past (-dt) and to the right (+dr).
!    Also remember that beta is the contravariant shift, and 
!    the interval requires the covariant shift, so we need to
!    multiply with the radial metric.

     if (i==Nr) then

        alpha_avg = 0.5d0*(alpha(l,Nr) + alpha_p(l,Nr))
        phi_avg   = 0.5d0*(  phi(l,Nr) +   phi_p(l,Nr))
        A_avg     = 0.5d0*(    A(l,Nr) +     A_p(l,Nr))
        B_avg     = 0.5d0*(    B(l,Nr) +     B_p(l,Nr))

        ra_right = r(l,Nr)*dexp(2.d0*phi_avg)*dsqrt(B_avg)

        gammarr = exp(4.d0*phi_avg)*A_avg

        dright = - (alpha_avg*dt(l))**2

        if (shift/="none") then
           beta_avg = 0.5d0*(beta(l,Nr) + beta_p(l,Nr))
           dright = dright + gammarr*(beta_avg*dt(l))**2
        end if

        CSI_right = CSI_SCHWARZ_P(l,Nr)
        ETA_right = ETA_SCHWARZ_P(l,Nr)

     else

        alpha_avg = 0.5d0*(alpha(l,i) + alpha_p(l,i+1))
        phi_avg   = 0.5d0*(  phi(l,i) +   phi_p(l,i+1))
        A_avg     = 0.5d0*(    A(l,i) +     A_p(l,i+1))
        B_avg     = 0.5d0*(    B(l,i) +     B_p(l,i+1))

        ra_right = 0.5d0*(r(l,i)+r(l,i+1))*dexp(2.d0*phi_avg)*dsqrt(B_avg)

        gammarr = exp(4.d0*phi_avg)*A_avg

        dright = - (alpha_avg*dt(l))**2 + gammarr*dr(l)**2

        if (shift/="none") then
           beta_avg = 0.5d0*(beta(l,i) + beta_p(l,i+1))
           dright = dright + gammarr*(beta_avg*dt(l))**2 - 2.d0*gammarr*beta_avg*dr(l)*dt(l)
        end if

        CSI_right = CSI_SCHWARZ_P(l,i+1)
        ETA_right = ETA_SCHWARZ_P(l,i+1)

     end if

!    Now solve for the position of the point that has precisely those
!    intervals in Kruskal-Szekeres coordinates.

     call solvepointkruskal(BHmass,CSI_guess,ETA_guess,CSI_left,ETA_left,CSI_right,ETA_right,dleft,dright,ra_left,ra_right,error)

!    Error message.

     if (error) then
          print *, 'Newton-Raphson iterations did not converge at point',i,'  , r =',r(l,i), ' , level l =',l
          print *, 'Aborting! (subroutine geometry/trackschwarz.f90)'
          print *
          !call die
     end if

!    Save new coordinates.

     CSI_SCHWARZ(l,i) = CSI_guess
     ETA_SCHWARZ(l,i) = ETA_guess

  end do

  goto 10

! For parallel runs synchronize across processors.

  if (size>1) then
     syncvar => CSI_SCHWARZ(l,:)
     call sync
     syncvar => ETA_SCHWARZ(l,:)
     call sync
  end if

! Restrict fine data to coarse grid.

  if ((l>0).and.(mod(s(l),2)==0)) then

!    Restrict.

     restrictvar => CSI_SCHWARZ
     call restrict(l,.false.)

     restrictvar => ETA_SCHWARZ
     call restrict(l,.false.)

!    Sync.

     if (size>1) then
        syncvar => CSI_SCHWARZ(l,:)
        call sync
        syncvar => ETA_SCHWARZ(l,:)
        call sync
     end if

  end if

  10 continue


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
        open(1,file=trim(directory)//'/SchwarzSlice'//trim(filen)//'.rl',form='formatted', &
             status=filestatus)
     else
        open(1,file=trim(directory)//'/SchwarzSlice'//trim(filen)//'.rl',form='formatted', &
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

!    Save data from processor 0. Notice that we ignore
!    the ghost points left of the origin for this.

     if (size==1) then
        imax = Nrl(0)
     else
        imax = Nrl(0) - ghost
     end if

     do i=imin,imax
        write(1,"(2ES16.8)") CSI_SCHWARZ(l,i),ETA_SCHWARZ(l,i)
     end do

!    Iterate over other processors.

     do p=1,size-1

!       Receive and save data from other processors.

        call MPI_RECV(var1,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(var2,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

        if (p==size-1) then
           imax = Nrl(p)
        else
           imax = Nrl(p) - ghost
        end if

        do i=1,imax
           write(1,"(2ES16.8)") var1(i),var2(i)
        end do

     end do

!    Leave two blank spaces before next time.
!    The reason to leave two spaces is that 'gnuplot' asks
!    for two spaces to distinguish different records.

     write(1,*)
     write(1,*)

!    Close files.

     close(1)

! Other processors send data to processor 0.

  else

     call MPI_SEND(CSI_SCHWARZ(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
     call MPI_SEND(ETA_SCHWARZ(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine trackschwarz







  subroutine solvepointkruskal(M,csi,eta,csi1,eta1,csi2,eta2,d1,d2,ra1,ra2,error)

! ***************************
! ***   SOLVE FOR POINT   ***
! ***************************

! This subroutine solves for the position in Kruskal
! coordinates of a point (csi,eta) that has an interval d1
! with respect to the point (csi1,eta1) and an interval d2
! with respect to the point (csi2,eta2).
!
! The routine does that by using Newton-Raphson.  Notice
! that on entry (r,t) contains the initial guess and on
! exit the solution.

  implicit none

  logical error

  integer k,Nmax

  real(8) M
  real(8) csi,csi1,csi2
  real(8) eta,eta1,eta2
  real(8) d1,d2,ra1,ra2
  real(8) det,norm
  real(8) epsilon
  real(8) aux1,aux2

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
! epsilon    Tolerance.


! **************************
! ***   NEWTON-RAPHSON   ***
! **************************

  epsilon = 1.d-10

! Start iterations.

  Nmax = 5000

  do k=1,Nmax

     aux1 = (32.d0*M**3/ra1)*dexp(-0.5d0*ra1/M)
     aux2 = (32.d0*M**3/ra2)*dexp(-0.5d0*ra2/M)

!    Find residuals using Kruskal metric in standard coordinates.
!    Notice that we do not need the angular coordinates.

     res(1) = aux1*((csi - csi1)**2 - (eta - eta1)**2) - d1
     res(2) = aux2*((csi - csi2)**2 - (eta - eta2)**2) - d2

!    Find Jacobian matrix and its determinant.

     jac(1,1) = + 2.0d0*aux1*(csi - csi1)
     jac(1,2) = - 2.0d0*aux1*(eta - eta1)

     jac(2,1) = + 2.0d0*aux2*(csi - csi2)
     jac(2,2) = - 2.0d0*aux2*(eta - eta2)

     det = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)

!    Solve for increments (J delta = - res).

     delta(1) = (res(2)*jac(1,2) - res(1)*jac(2,2))/det
     delta(2) = (res(1)*jac(2,1) - res(2)*jac(1,1))/det

!    Correct solution.

     csi = csi + delta(1)
     eta = eta + delta(2)

!    Check if we achieved desired tolerance.

     norm = dabs(delta(1)) + dabs(delta(2))

     if (norm<epsilon) then
        !d1 = aux1*((csi - csi1)**2 - (eta - eta1)**2)
        !d2 = aux2*((csi - csi2)**2 - (eta - eta2)**2)
        return
     end if

  end do

! If we get here we the algorithm did not converge, so we
! set the error flag to .true.

  error = .true.


! ***************
! ***   END   ***
! ***************

  end subroutine solvepointkruskal


