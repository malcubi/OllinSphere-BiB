!$Header: $

  subroutine idata_chargedprocastar

! ***************************************************************
! ***   CHARGED PROCA STAR INITIAL DATA IN POLAR AREAL GAUGE  ***
! ***************************************************************

! This subroutine calculates initial data for a charged
! Proca star using a shooting method in the polar areal gauge.

! MAIN PROBLEM: The main problem with this routine is that
! it requires a very fine tuned initial guess since we are
! fighting a growing exponential.  This gets worse as
! we move the boundary further away. So the way to
! proceed is to put the boundary close in, find omega,
! and then move it slowly further out adjusting omega
! to high precision as we move out.

! Charged Proca stars are self-gravitating solutions of a complex
! massive vector field (Proca field) such that the spacetime
! is static and the Proca field has a harmonic time dependence.

! NOTE FOR PARALLEL RUNS:  The initial data is not really
! solved in parallel.  It is in fact solved only on processor
! zero on a full size array, and then it is distributed
! among the other processors.  This is slow, but works.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays
  use derivatives

! Extra variables.

  implicit none

  logical :: left  = .false.             ! This flag is true if the left  value of omega is well behaved.
  logical :: right = .false.             ! This flag is true if the right value of omega is well behaved.

  integer i,l,iter                       ! Counters.
  integer imin                           ! Leftmost grid point.
  integer iaux                           ! Auxiliary quantity.
  integer :: maxiter = 200               ! Maximum number of iterations.

  real(8) r0,delta                       ! Local radius and grid spacing.
  real(8) A0,alpha0                      ! Initial values of (A,alpha).
  real(8) procaF0,procaA0                ! Initial values of Proca variables.
  real(8) maxwellF0,maxwellE0            ! Initial values of Maxwell variables. 

  real(8) k11,k12,k13,k14                ! Runge-Kutta sources for A.
  real(8) k21,k22,k23,k24                ! Runge-Kutta sources for alpha.
  real(8) k31,k32,k33,k34                ! Runge-Kutta sources for procaF.
  real(8) k41,k42,k43,k44                ! Runge-Kutta sources for procaA.
  real(8) k51,k52,k53,k54                ! Runge-Kutta sources for maxwellE.
  real(8) k61,k62,k63,k64                ! Runge-Kutta sources for maxwellF.

  real(8) A_rk,alpha_rk                  ! Runge-Kutta values of (A,alpha).
  real(8) procaF_rk,procaA_rk,procaE_rk  ! Runge-Kutta values of Proca variables.
  real(8) maxwellF_rk,maxwellE_rk        ! Runge-Kutta values of Maxwell variables.
  real(8) cps_omega_rk                   ! Modified frequency value 
  real(8) J1_CPS,J2_CPS,J3_CPS           ! Functions for sources of differential equations.
  real(8) J4_CPS,J5_CPS,J6_CPS           ! Functions for sources of differential equations.

  real(8) res,res_old                    ! Residual.
  real(8) omega_new,omega_old,domega     ! Trial frequency and frequency interval. 
  real(8) half,smallpi                   ! Numbers.
  real(8) rm,alphafac,Ffac,aux           ! Auxiliary quantities.
  real(8) :: epsilon = 1.d-8             ! Tolerance.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                        ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g,alpha_g               ! Radial metric and lapse global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaF_g,procaPhi_g       ! Global arrays for Proca fields (F,Phi).
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaA_g,procaE_g         ! Global arrays for Proca fields (X,E).
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: maxwellF_g,maxwellPhi_g   ! Global arrays for Maxwell fields (F,Phi)
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: maxwellE_g                ! Global array for Maxwell field (E)


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  smallpi = acos(-1.d0)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for a charged Proca star in polar-areal gauge using shooting method'
     print *
  end if

! Sanity check.

  if (spacetime=="minkowski") then
     print *, 'Charged Proca star initial data is not compatible with a Minkowski background ...'
     print *, 'Aborting! (subroutine idata_ChargedProcastar)'
     print *
     call die
  end if


! *************************
! ***   NORMALIZATION   ***
! *************************

! Physical normalization.

  if (proca_factor=="physical") then

     if (rank==0) then
        write(*,'(A)') ' Using physical normalization at origin: procaPhi(r=0) = proca_phi0'
     end if

! Alternative normalization.

  else

     proca_phi0  = proca_phi0 /sqrt(4.d0*smallpi)

     if (rank==0) then
        write(*,'(A,E16.10)') ' Using harmonic normalization at origin: procaPhi(r=0) = proca_phi0 /sqrt(4pi) = ', proca_phi0
     end if

  end if

! Rescale the tolerance with proca_phi0, so the final tolerance
! is measured with respect to this value.

  epsilon = epsilon*proca_phi0


! ********************************
! ***   CHARGE NORMALIZATION   ***
! ********************************


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

  do l=0,Nl-1
     do i=1-ghost,Nrtotal
        rr(l,i) = (dble(Nmin(rank) + i) - half)*dr(l)
     end do
  end do


! ***********************
! ***   INTEGRATION   ***
! ***********************

! Here we use a shooting method to get the correct
! asymptotic behaviour.  On each iteration, the
! integration is done with fourth order Runge-Kutta.

! For the moment only processor 0 solves for the
! initial data.  The solution is later distributed
! among processors.

! Initialize arrays.

  A_g      = 1.d0
  alpha_g  = 1.d0

  procaF_g = 0.d0
  procaA_g = 0.d0
  procaE_g = 0.d0

  maxwellF_g = 0.d0
  maxwellE_g = 0.d0


! *********************************************
! ***   ONLY PROCESSOR 0 SOLVES THE ODE's   ***
! *********************************************

  if (rank==0) then


!    **********************
!    ***   INITIALIZE   ***
!    **********************

!    Find initial domega.

     domega = (omega_right - omega_left)/20.0d0

!    Sanity check.

     if (domega<0.d0) then
        print *
        print *, 'For Proca star initial data we must have omega_right>omega_left.'
        print *, 'Aborting! (subroutine idata_ChargedProcastar)'
        print *
        call die
     end if

!    If omega_right=omega_left it means we want a fixed frequency with NO iterations.

     if (domega==0.d0) then
        maxiter = 1
     end if

!    Initialize residual.

     res = 1.d0


!    ****************************
!    ***   BEGIN ITERATIONS   ***
!    ****************************

     proca_omega = omega_left

     iter = 0

     do while ((abs(res).gt.epsilon).and.(iter.lt.maxiter))

        iter = iter + 1


!       *********************************
!       ***   LOOP OVER GRID LEVELS   ***
!       *********************************

!       We solve from fine to coarse grid.

        do l=Nl-1,0,-1


!          *******************************
!          ***   FIND INITIAL VALUES   ***
!          *******************************

!          Find initial point. Only the finest grid
!          integrates from the origin.

           if (l==Nl-1) then
              imin = 1
           else
              imin = Nrtotal/2
           end if

!          For coarse grids we interpolate the initial point
!          (which is half way through the grid) from the last
!          points of the higher resolution grid.  Here we use
!          cubic interpolation.

           if (l<Nl-1) then

              A_g(l,imin-1)         = (9.d0*(A_g(l+1,Nrtotal-2)+A_g(l+1,Nrtotal-3)) &
                                      - (A_g(l+1,Nrtotal-4)+A_g(l+1,Nrtotal-1)))/16.d0
              alpha_g(l,imin-1)     = (9.d0*(alpha_g(l+1,Nrtotal-2)+alpha_g(l+1,Nrtotal-3)) &
                                      - (alpha_g(l+1,Nrtotal-4)+alpha_g(l+1,Nrtotal-1)))/16.d0

              procaF_g(l,imin-1)    = (9.d0*(procaF_g(l+1,Nrtotal-2)+procaF_g(l+1,Nrtotal-3)) &
                                      - (procaF_g(l+1,Nrtotal-4)+procaF_g(l+1,Nrtotal-1)))/16.d0
              procaA_g(l,imin-1)    = (9.d0*(procaA_g(l+1,Nrtotal-2)+procaA_g(l+1,Nrtotal-3)) &
                                      - (procaA_g(l+1,Nrtotal-4)+procaA_g(l+1,Nrtotal-1)))/16.d0

              maxwellE_g(l,imin-1) = (9.d0*(maxwellE_g(l+1,Nrtotal-2)+maxwellE_g(l+1,Nrtotal-3)) &
                                      - (maxwellE_g(l+1,Nrtotal-4)+maxwellE_g(l+1,Nrtotal-1)))/16.d0
              maxwellF_g(l,imin-1) = (9.d0*(maxwellF_g(l+1,Nrtotal-2)+maxwellF_g(l+1,Nrtotal-3)) &
                                      - (maxwellF_g(l+1,Nrtotal-4)+maxwellF_g(l+1,Nrtotal-1)))/16.d0
 
           end if


!          ************************************
!          ***   FOURTH ORDER RUNGE-KUTTA   ***
!          ************************************

           do i=imin,Nrtotal

!             Grid spacing and values at first point
!             if we start from the origin (finer grid).

              if (i==1) then

!                For the first point we use dr/2.

                 delta = half*dr(l)
                 r0    = 0.d0

!                Values of (alpha,A) at origin.

                 A0     = 1.d0
                 alpha0 = 1.d0

!                Values of (procaF,procaA) at origin.

                 procaF0 = proca_phi0 
                 procaA0 = 0.d0

!                Values of (maxwellF,maxwellE) at origin.

                 maxwellF0 = 0.d0
                 maxwellE0 = 0.d0

!             Grid spacing and values at previous grid point.

              else

                 delta = dr(l)
                 r0    = rr(l,i-1)

                 A0     = A_g(l,i-1)
                 alpha0 = alpha_g(l,i-1)

                 procaF0 = procaF_g(l,i-1)
                 procaA0 = procaA_g(l,i-1)
                 
                 maxwellE0 = maxwellE_g(l,i-1)
                 maxwellF0 = maxwellF_g(l,i-1)                

              end if

!             I) First Runge-Kutta step.

!             Sources at first grid point if we start
!             from the origin (for finer grid).

              if (i==1) then

!                At the origin we have:  A' = alpha' = procaF' = 0.

                 k11 = 0.d0
                 k21 = 0.d0
                 k31 = 0.d0

!                For procaA we have:  procaA' = omega*procaF0/3

                 k41 = delta*(proca_omega*procaF0/3.d0)
                 
!                For maxwellE and maxwellF we have we have:  maxwellE' = maxwellF' = 0

                 k51 = 0.d0
                 k61 = 0.d0

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk     = A0
                 alpha_rk = alpha0

                 procaF_rk = procaF0
                 procaA_rk = procaA0

                 maxwellE_rk = maxwellE0
                 maxwellF_rk = maxwellF0

!                Sources.

                 cps_omega_rk = proca_omega + cproca_q*maxwellF_rk
                 procaE_rk = - cproca_mass**2*alpha_rk*procaA_rk/(cps_omega_rk*A_rk)

                 k11 = delta*J1_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
                 k21 = delta*J2_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
                 k31 = delta*J3_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
                 k41 = delta*J4_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
                 k51 = delta*J5_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
                 k61 = delta*J6_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk     = A0     + half*k11
              alpha_rk = alpha0 + half*k21

              procaF_rk = procaF0 + half*k31
              procaA_rk = procaA0 + half*k41

              maxwellE_rk = maxwellE0 + half*k51
              maxwellF_rk = maxwellF0 + half*k61

              cps_omega_rk = proca_omega + cproca_q*maxwellF_rk
              procaE_rk = - cproca_mass**2*alpha_rk*procaA_rk/(cps_omega_rk*A_rk)

              k12 = delta*J1_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k22 = delta*J2_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k32 = delta*J3_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k42 = delta*J4_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k52 = delta*J5_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k62 = delta*J6_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk     = A0     + half*k12
              alpha_rk = alpha0 + half*k22

              procaF_rk = procaF0 + half*k32
              procaA_rk = procaA0 + half*k42

              maxwellE_rk = maxwellE0 + half*k52
              maxwellF_rk = maxwellF0 + half*k62

!             Sources.

              cps_omega_rk = proca_omega + cproca_q*maxwellF_rk
              procaE_rk = - cproca_mass**2*alpha_rk*procaA_rk/(cps_omega_rk*A_rk)

              k13 = delta*J1_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k23 = delta*J2_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k33 = delta*J3_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k43 = delta*J4_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k53 = delta*J5_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k63 = delta*J6_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk     = A0     + k13
              alpha_rk = alpha0 + k23

              procaF_rk = procaF0 + k33
              procaA_rk = procaA0 + k43
              
              maxwellE_rk = maxwellE0 + k53
              maxwellF_rk = maxwellF0 + k63

!             Sources.

              cps_omega_rk = proca_omega + cproca_q*maxwellF_rk
              procaE_rk = - cproca_mass**2*alpha_rk*procaA_rk/(cps_omega_rk*A_rk)

              k14 = delta*J1_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k24 = delta*J2_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k34 = delta*J3_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k44 = delta*J4_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k54 = delta*J5_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k64 = delta*J6_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)

!             Advance variables to next grid point.

              A_g(l,i)     = A0     + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i) = alpha0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

              procaF_g(l,i) = procaF0 + (k31 + 2.d0*(k32 + k33) + k34)/6.d0
              procaA_g(l,i) = procaA0 + (k41 + 2.d0*(k42 + k43) + k44)/6.d0
              
              maxwellE_g(l,i) = maxwellE0 + (k51 + 2.d0*(k52 + k53) + k54)/6.d0
              maxwellF_g(l,i) = maxwellF0 + (k61 + 2.d0*(k62 + k63) + k64)/6.d0

!             Check if solution is blowing up. This helps to reduce
!             the need for a very fine tuned initial guess.

              if (abs(procaF_g(l,i))>2.d0*proca_phi0 ) then
                 if (.not.left) then
                    omega_left = omega_left + domega
                    boson_omega = omega_left
                    goto 100
                 else if (.not.right) then
                    omega_right = omega_right - domega
                    boson_omega = omega_right
                    goto 100
                 end if
              end if

!             Check if solution is already very small.

              if (abs(procaF_g(l,i))+abs(procaF_g(l,i-1))<epsilon) then
                 procaF_g(l,i) = 0.d0
              end if

           end do


!          ***********************
!          ***   GHOST ZONES   ***
!          ***********************

           do i=1,ghost

              A_g(l,1-i)      = + A_g(l,i)
              alpha_g(l,1-i)  = + alpha_g(l,i)

              procaF_g(l,1-i) = + procaF_g(l,i)
              procaA_g(l,1-i) = - procaA_g(l,i)

              maxwellF_g(l,1-i) = + maxwellF_g(l,i)      
              maxwellE_g(l,1-i) = - maxwellE_g(l,i)

           end do

        end do


!       *******************************************
!       ***   FIND RESIDUAL AND CORRECT OMEGA   ***
!       *******************************************

!       Find the difference between the solution at the coarsest level
!       and the expected asymptotic behaviour. The solution should
!       decay exponetially far away as:
!
!       F ~ exp(-k*r)
!
!       so that:
!
!       F' + k F = 0
!
!       Notice that far away we must have:  k = sqrt(m**2 - (omega/alpha)**2)

        res_old = res

        if (abs(procaF_g(0,Nrtotal))+abs(procaF_g(0,Nrtotal-1))<epsilon) then
           res = epsilon/2.d0
           goto 100
        else
           aux = abs(cproca_mass**2 - (proca_omega/alpha_g(0,Nrtotal))**2)
           res = (procaF_g(0,Nrtotal) - procaF_g(0,Nrtotal-1))/dr(0) + dsqrt(aux)*procaF_g(0,Nrtotal)
        end if

!       Secant method:  Having found the difference for the two values of omega
!       we can use it to extrapolate linearly to the next best guess.

        if (.not.left) then

           left = .true.

           omega_old = proca_omega
           proca_omega = omega_right

        else

           if (.not.right) right = .true.

           aux = res - res_old

           if ((aux/=0.d0).and.(abs(aux)<1.d-1)) then
              omega_new = proca_omega - res*(proca_omega - omega_old)/aux
           else
              omega_new = 0.5d0*(omega_old + proca_omega)
           end if

           omega_old = proca_omega
           proca_omega = omega_new

        end if

100     continue


!       **************************
!       ***   END ITERATIONS   ***
!       **************************

!       Output data to screen.

        write(*,"(A,I4,A,ES22.16,A,ES9.2)") ' Iteration: ',iter,'    Frequency: ',proca_omega,'    Residual: ',res

     end do

!    Message in case we reached the maximum number of iterations.

     if (iter>=maxiter) then
        print *
        print *, 'Maximum iteration number reached, initial data solver did not converge.'
        print *
     else
        print *
        print *, 'Done!'
        print *
     end if


!    ************************************
!    ***   RESTRICT TO COARSE GRIDS   ***
!    ************************************

!    Restrict solution from fine to coarse grid.
!    We don't call the subroutine "restrict"
!    since here we are running only on processor 0.
!    We use cubic interpolation.

     do l=Nl-1,1,-1

        do i=1,Nrtotal-ghost,2

           iaux = i/2 + 1
           rm = rr(l-1,iaux)

           A_g(l-1,iaux)      = (9.d0*(A_g(l,i)+A_g(l,i+1))         - (A_g(l,i-1)+A_g(l,i+2)))/16.d0
           alpha_g(l-1,iaux)  = (9.d0*(alpha_g(l,i)+alpha_g(l,i+1)) - (alpha_g(l,i-1)+alpha_g(l,i+2)))/16.d0

           procaF_g(l-1,iaux) = (9.d0*(procaF_g(l,i)+procaF_g(l,i+1)) - (procaF_g(l,i-1)+procaF_g(l,i+2)))/16.d0
           procaA_g(l-1,iaux) = (9.d0*(procaA_g(l,i)+procaA_g(l,i+1)) - (procaA_g(l,i-1)+procaA_g(l,i+2)))/16.d0

           maxwellE_g(l-1,iaux) = (9.d0*(maxwellE_g(l,i)+procaF_g(l,i+1)) - (maxwellE_g(l,i-1)+maxwellE_g(l,i+2)))/16.d0
           maxwellF_g(l-1,iaux) = (9.d0*(maxwellF_g(l,i)+maxwellF_g(l,i+1)) - (maxwellF_g(l,i-1)+maxwellF_g(l,i+2)))/16.d0

        end do

!       Fix ghost zones.

        do i=1,ghost

           A_g(l-1,1-i)      = + A_g(l-1,i)
           alpha_g(l-1,1-i)  = + alpha_g(l-1,i)

           procaF_g(l-1,1-i) = + procaF_g(l-1,i)
           procaA_g(l-1,1-i) = - procaA_g(l-1,i)

           maxwellF_g(l-1,1-i) = + maxwellF_g(l-1,i)
           maxwellE_g(l-1,1-i) = - maxwellE_g(l-1,i)

        end do

     end do


!    ***********************************
!    ***   RESCALE (ALPHA,OMEGA,F)   ***
!    ***********************************

!    Normalize lapse and omega.  Since we integrated out
!    assuming alpha(r=0)=1, we now need to rescale the lapse.
!    This can be done easily since its equation is linear.
!    But notice that this also rescales the final frequency omega.
!
!    In practice we just divide both alpha and omega by
!    the asymptotic value of the lapse.  In order to obtain
!    this value we assume that far away the lapse behaves
!    as alpha = alpha0 + const/r.  This implies that:
!
!    alpha0  =  alpha + r dalpha/dr
!
!    We calculate dalpha/dr at the boundary with one sided
!    differences and solve for alpha0.

     if (order=="two") then
        alphafac = alpha_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
            *(3.d0*alpha_g(0,Nrtotal) - 4.d0*alpha_g(0,Nrtotal-1) + alpha_g(0,Nrtotal-2))
     else
        alphafac = alpha_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
            *(25.d0*alpha_g(0,Nrtotal) - 48.d0*alpha_g(0,Nrtotal-1) &
            + 36.d0*alpha_g(0,Nrtotal-2) - 16.d0*alpha_g(0,Nrtotal-3) + 3.d0*alpha_g(0,Nrtotal-4))/3.d0
     end if

     alpha_g = alpha_g/alphafac

!    Write value of omega to screen and rescale it.

     omega_new = proca_omega/alphafac

     if (rank==0) then
        write(*,'(A,E22.16)') ' Omega (not-rescaled) = ', proca_omega
        write(*,'(A,E22.16)') ' Omega (rescaled)     = ', omega_new
        print *
     end if

!    Rescale procaF and maxwellF.

     procaF_g   = procaF_g/alphafac
     maxwellF_g = maxwellF_g/alphafac


!    *******************************************
!    ***   RECONSTRUCT procaPhi AND procaE   ***
!    *******************************************

!    Having found the solution, the scalar potential
!    is recovered as: phi = F/alpha.

     procaPhi_g   = procaF_g/alpha_g
     maxwellPhi_g = maxwellF_g/alpha_g

!    The Proca electric field is given by:
! 
!    E = - m**2 alpha procaA / [A (omega + q*maxwellF)]

     procaE_g = - cproca_mass**2*alpha_g*procaA_g/(A_g*(omega_new + cproca_q*maxwellF_g))


!    ***********************************
!    ***   PROCA STAR PERTURBATION   ***
!    ***********************************


! *************************************
! ***   FINISHED FINDING SOLUTION   ***
! *************************************

! When we get here the solution has been found,
! so we close the "if" statement for processor 0.

  end if


! ************************************************
! ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
! ************************************************

! For single processor runs just copy the arrays.
! Remember that at t=0 procaA is purely imaginary.
 
  if (size==1) then

     A     = A_g
     alpha = alpha_g

     cprocaPhi_R = procaPhi_g
     cprocaA_I   = procaA_g
     cprocaE_R   = procaE_g

     ePhi     = maxwellPhi_g
     electric = maxwellE_g

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! arrays with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  else

     call distribute(0,Nl-1,A,A_g)
     call distribute(0,Nl-1,alpha,alpha_g)

     call distribute(0,Nl-1,cprocaPhi_R,procaPhi_g)
     call distribute(0,Nl-1,cprocaA_I,procaA_g)
     call distribute(0,Nl-1,cprocaE_R,procaE_g)

     call distribute(0,Nl-1,ePhi,maxwellPhi_g)
     call distribute(0,Nl-1,electric,maxwellE_g)

  end if


! ************************************************************************
! ***   IMAGINARY PARTS OF (procaPhi,procaE) AND REAL PART OF procaA   ***
! ************************************************************************

! Set imaginary parts of procaPhi and procaE to zero.

  cprocaPhi_I = 0.d0
  cprocaE_I   = 0.d0

! Set real part of procaA to zero.

  cprocaA_R = 0.d0

! Also set Maxwell vector potential to zero.

  eAr = 0.d0


! ***************
! ***   END   ***
! ***************

  end subroutine idata_ChargedProcastar







! **********************************
! ***   RADIAL DERIVATIVE OF A   ***
! **********************************

! The radial derivative of A comes from the
! Hamiltonian constraint.

  function J1_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J1_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm
  real(8) rho

!                              2           2        2                   2          2  
! rho = + 1/(8 pi) { A ( procaE  + maxwellE  )  +  m  [ (procaF / alpha)  +  procaA / A ] }
!
! Notice that we don't divide by 8*pi since it cancels.

  rho = A*(procaE**2 + maxwellE**2) + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A)

! dA/dr  =  A [ (1 - A)/r + 8 pi r A rho ]
! 
! Notice that we don't multiply the last term with 8*pi since it cancels.

  J1_CPS = A*((1.d0-A)/rm + rm*A*rho)

  end function J1_CPS







! **************************************
! ***   RADIAL DERIVATIVE OF ALPHA   ***
! **************************************

! The radial derivative of alpha comes from
! the polar slicing condition:  dK_{theta,theta}/dt=0.

  function J2_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J2_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm
  real(8) SA

!                              2           2        2                   2          2
! SA  = - 1/(8 pi) { A ( procaE  + maxwellE  )  -  m  [ (procaF / alpha)  +  procaA / A ] }
!
! Notice that we don't divide by 8*pi since it cancels. 

  SA = - A*(procaE**2 + maxwellE**2) + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A)

! dalpha/dr  =  alpha [ (A - 1)/(2r) + 4 pi r A SA ]
! 
! Notice that we don't multiply the last term with 4*pi*A since
! it cancels (we are left only with 1/2).

  J2_CPS = alpha*(0.5d0*(A-1.d0)/rm + 0.5d0*rm*A*SA)

  end function J2_CPS







! ***************************************
! ***   RADIAL DERIVATIVE OF procaF   ***
! ***************************************

! The radial derivative of procaF from the
! Proca evolution equations.

  function J3_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J3_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm
  real(8) W

! dprocaF/dr = procaA W { [(m alpha/ W)**2 ] - 1 }
!
! where:  W  :=  omega + q maxwellF

  W = proca_omega + cproca_q*maxwellF

  J3_CPS = procaA*W*((cproca_mass*alpha/W)**2 - 1.d0)

  end function J3_CPS







! ***************************************
! ***   RADIAL DERIVATIVE OF procaA   ***
! ***************************************

! The radial derivative of procaA from
! the Procan evolution equations.

  function J4_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J4_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm
  real(8) aux,W,correc
  
!                         2            2
! SA - rho  = -  (A procaE + A maxwellE ) / 4 pi
! 
! Notice that we don't divide by 4*pi since it cancels.
 
  aux = - A*(procaE**2 + maxwellE**2)
              
! dprocaA/dr =  W F A /alpha**2  -  procaA [ (A+1)/r + 4 pi r A (SA - rho) ]
!
!            - (q alpha A maxwellE procaA) / W
!
! where:  W  :=  omega + q maxwellF

  W = proca_omega + cproca_q*maxwellF
  correc = cproca_q*alpha*A*maxwellE*procaA/W

  J4_CPS = W*procaF*A/alpha**2 - procaA*((A + 1.d0)/rm + rm*A*aux) - correc

  end function J4_CPS










! *****************************************
! ***   RADIAL DERIVATIVE OF maxwellE   ***
! *****************************************

! The radial derivative of maxwellF comes from
! Maxwell's equations.

  function J5_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J5_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm
  real(8) W,rho

!                              2           2        2                   2          2  
! rho = + 1/(8 pi) { A ( procaE  + maxwellE  )  +  m  [ (procaF / alpha)  +  procaA / A ] }
!
! Notice that we don't divide by 8*pi since it cancels.

  rho = A*(procaE**2 + maxwellE**2) + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A)

! dmaxwellE/dr =  -  2 q alpha m**2 procaA**2 / (A W)  -  maxwellE [ (5-A) / 2r + 4 pi A rho ]
!
! where:  W  :=  omega + q maxwellF

  W = proca_omega + cproca_q*procaF

! Notice that we don't multiply the last term with 4*pi*A since
! it cancels.

  J5_CPS = - 2.d0*cproca_q*alpha*cproca_mass**2*procaA**2/(A*W) &
         - maxwellE*((5.d0-A)/(2.d0*rm) + 0.5d0*rm*A*rho)

  end function J5_CPS







! *****************************************
! ***   RADIAL DERIVATIVE OF maxwellF   ***
! *****************************************

! The radial derivative of procaA from the Proca
! evolution equations.

  function J6_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J6_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm

! dmaxwellF/dr  =  - alpha A maxwellE

  J6_CPS = - alpha*A*maxwellE

  end function J6_CPS



