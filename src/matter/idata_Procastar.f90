
  subroutine idata_Procastar

! *************************************************************
! ***   PROCA STAR INITIAL DATA IN CONFORMALLY FLAT GAUGE   ***
! *************************************************************

! This subroutine calculates initial data for a Proca star
! using a shooting method in the conformally flat gauge.

! MAIN PROBLEM: The main problem with this routine is that
! it requires a very fine tuned initial guess since we are
! fighting a growing exponential.  This gets worse as
! we move the boundary further away. So the way to
! proceed is to put the boundary close in, find omega,
! and then move it slowly further out adjusting omega
! to high precision as we move out.

! Proca stars are self-gravitating solutions of a complex
! massive vector field (Proca field) such that the spacetime
! is static and the Proca field has a harmonic time dependence.
!
! To obtain the initial data we assume that spacetime is
! static (K_ij=0), and also that the complex scalar and
! vector potentials (Phi,A) have the form:
!
! procaPhi(t,r) = phi(r) exp(-i omega t)
! procaA(t,r)   = i a(r) exp(-i omega t)
!
! Notice that with this ansatz the stress-energy tensor
! is time-independent so the metric can be static. Also,
! at t=0 procaPhi is purely real while procaA is purely
! imaginary.
!
! The standard ansatz for the metric for boson stars is:
!
!   2          2   2        2     2      2
! ds  = - alpha  dt  +  A dr  +  r dOmega
!
! In other words, we are assuming psi=B=1, beta=0, and
! the radial coordinate is the area radius.
!
! It is in fact much better to work with the quantity F := alpha*phi
! instead of phi itself, and then reconstruct phi at the end.
! With the above ansatz, and using F, the Proca evolution
! equations reduce to:
!
!                                       2
! dF/dr  =  a omega[ ( m alpha / omega )  -  1 ]
!
!                            2
! da/dr  =  omega F A / alpha  -  a [ (A+1)/r + 4 pi r A (SA - rho) ]
!
!
! with rho the energy density given by:
!
!                           2     2           2    2
! rho  =  + 1 / (8 pi) { A E  +  m [ (F/alpha)  + a / A ] ]
!
!
! and SA the radial stress given by:
!
!                           2     2           2    2
! SA   =  - 1 / (8 pi) { A E  -  m [ (F/alpha)  + a / A ] ]
!
!
! with:
!
!          2
! E  =  - m  alpha a / omega A
!
!
! In particular:
!
!                   2
! SA - rho  = -  A E / 4 pi
!
!
! On the other hand, the Hamiltonian constraint takes the form:
!
!                      
! dA/dr  =  A [ (1 - A)/r + 8 pi r A rho ]
!
!
! For the lapse we use the polar slicing condition
! K_{theta,theta}=0, which when combined with the
! hamiltonian constraint above takes the form:
!
!
! dalpha/dr  =  alpha [ (A - 1)/(2r) + 4 pi r A SA ]
!
!
! This routine takes as input parameter the value of the scalar
! potential at the origin "proca_phi0", and uses a shooting method
! (with Runge-Kutta) to determine the value of the frequency
! "proca_omega" that is compatible with the asymptotically
! decaying solution.
!
! For the boundary conditions at the origin we take:
!
! A(r=0)     = 1,     d A(r=0)     = 0
!                      r
!
! alpha(r=0) = 1,     d alpha(r=0) = 0
!                      r
!
! phi(r=0)   = phi0,  d phi(r=0)   = 0
!                      r     
!
! a(r=0)     = 0,     d a(r=0)     = constant
!                      r     
!
! Notice that "a", as the radial component of the vector
! potential, must go as  a ~ (k_a) r  close to the origin, with
! k_a constant. The value of this constant can be determined
! by evaluating the equations above for small r, and can be
! shown to be:  k_a = omega phi0 / 3.
!
! Also, in order to avoid division by zero we need to expand
! the radial metric close to the origin as  A ~ 1  + k_A r^2.
! The value of k_A can also be found by evaluating the equations
! close to the origin ang is given by:  k_A = (8pi/3) rho0.
!
! Notice that in the final solution we don't want alpha(r=0)=1,
! but rather alpha=1 at infinity.  But this is no problem as the
! slicing condition above is linear in alpha so we can always
! just rescale the lapse at the end, but in order not to affect
! the solution for (A,phi,a) we must also rescale the final value
! of the frequency omega by the same factor.

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

  logical :: left  = .false.            ! This flag is true if the left  value of omega is well behaved.
  logical :: right = .false.            ! This flag is true if the right value of omega is well behaved.

  integer i,l,iter                      ! Counters.
  integer imin                          ! Leftmost grid point.
  integer iaux                          ! Auxiliary quantity.
  integer :: maxiter = 200              ! Maximum number of iterations.

  real(8) r0,delta                      ! Local radius and grid spacing.
  real(8) A0,alpha0                     ! Initial values of variables.
  real(8) procaF0,procaA0,procaE0       ! Initial values of variables.
  real(8) k11,k12,k13,k14               ! Runge-Kutta sources for A.
  real(8) k21,k22,k23,k24               ! Runge-Kutta sources for alpha.
  real(8) k31,k32,k33,k34               ! Runge-Kutta sources for procaF.
  real(8) k41,k42,k43,k44               ! Runge-Kutta sources for procaA.
  real(8) k51,k52,k53,k54               ! Runge-Kutta sources for procaE.
  real(8) A_rk,alpha_rk                 ! Runge-Kutta values of variables.
  real(8) procaF_rk,procaA_rk,procaE_rk ! Runge-Kutta values of variables.
  real(8) J1_proca,J2_proca             ! Functions for sources of differential equations.
  real(8) J3_proca,J4_proca,J5_proca    ! Functions for sources of differential equations.
  real(8) res,res_old                   ! Residual.
  real(8) omega_new,omega_old,domega    ! Trial frequency and frequency interval. 
  real(8) half,smallpi                  ! Numbers.
  real(8) rm,alphafac,aux               ! Auxiliary quantities.
  real(8) :: epsilon = 1.d-8            ! Tolerance.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                    ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g,alpha_g           ! Global arrays for radial metric and lapse.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaPhi_g,procaA_g   ! Global arrays for potentials.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaF_g              ! Global array for F = alpha*phi.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaE_g              ! Global array for electric field.


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  smallpi = acos(-1.d0)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for a Proca star in polar-areal gauge using shooting method ...'
     print *
  end if

! Sanity check.

  if (spacetime=="minkowski") then
     print *, 'Proca star initial data is not compatible with a Minkowski background ...'
     print *, 'Aborting! (subroutine idata_Procastar)'
     print *
     call die
  end if


! *************************
! ***   NORMALIZATION   ***
! *************************

! Physical normalization.

  if (proca_factor=="physical") then

     if (rank==0) then
        write(*,'(A)') ' Using physical normalization at origin: phi(r=0) = phi0'
     end if

! Alternative normalization.

  else

     proca_phi0 = proca_phi0/sqrt(4.d0*smallpi)

     if (rank==0) then
        write(*,'(A,E16.10)') ' Using harmonic normalization at origin: phi(r=0) = phi0/sqrt(4pi) = ', proca_phi0
     end if

  end if

! Rescale the tolerance with boson_phi0, so the final tolerance
! is measured with respect to this value.

  epsilon = epsilon*proca_phi0


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
        print *, 'Aborting! (subroutine idata_Procastar)'
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

              A_g(l,imin-1)     = (9.d0*(A_g(l+1,Nrtotal-2)+A_g(l+1,Nrtotal-3)) &
                                - (A_g(l+1,Nrtotal-4)+A_g(l+1,Nrtotal-1)))/16.d0
              alpha_g(l,imin-1) = (9.d0*(alpha_g(l+1,Nrtotal-2)+alpha_g(l+1,Nrtotal-3)) &
                                - (alpha_g(l+1,Nrtotal-4)+alpha_g(l+1,Nrtotal-1)))/16.d0

              procaF_g(l,imin-1) = (9.d0*(procaF_g(l+1,Nrtotal-2)+procaF_g(l+1,Nrtotal-3)) &
                                 - (procaF_g(l+1,Nrtotal-4)+procaF_g(l+1,Nrtotal-1)))/16.d0
              procaA_g(l,imin-1) = (9.d0*(procaA_g(l+1,Nrtotal-2)+procaA_g(l+1,Nrtotal-3)) &
                                 - (procaA_g(l+1,Nrtotal-4)+procaA_g(l+1,Nrtotal-1)))/16.d0

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

!             Grid spacing and values at previous grid point.

              else

                 delta = dr(l)
                 r0    = rr(l,i-1)

                 A0     = A_g(l,i-1)
                 alpha0 = alpha_g(l,i-1)

                 procaF0 = procaF_g(l,i-1)
                 procaA0 = procaA_g(l,i-1)

              end if

!             I) First Runge-Kutta step.

!             Sources at first grid point if we start
!             from the origin (for finer grid).

              if (i==1) then

!                At the origin we have:  A' = alpha' = procaF' = 0.

                 k11 = 0.d0
                 k21 = 0.d0
                 k31 = 0.d0

!                For procaA we have:  procaA' = omega*F0/3.

                 k41 = delta*(proca_omega*procaF0/3.d0)

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk     = A0
                 alpha_rk = alpha0

                 procaF_rk = procaF0
                 procaA_rk = procaA0

!                Sources.

                 procaE_rk = - cproca_mass**2*alpha_rk*procaA_rk/(proca_omega*A_rk)

                 k11 = delta*J1_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
                 k21 = delta*J2_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
                 k31 = delta*J3_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
                 k41 = delta*J4_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk     = A0     + half*k11
              alpha_rk = alpha0 + half*k21

              procaF_rk = procaF0 + half*k31
              procaA_rk = procaA0 + half*k41

!             Sources.

              procaE_rk = - cproca_mass**2*alpha_rk*procaA_rk/(proca_omega*A_rk)

              k12 = delta*J1_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k22 = delta*J2_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k32 = delta*J3_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k42 = delta*J4_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk     = A0     + half*k12
              alpha_rk = alpha0 + half*k22

              procaF_rk = procaF0 + half*k32
              procaA_rk = procaA0 + half*k42

!             Sources.

              procaE_rk = - cproca_mass**2*alpha_rk*procaA_rk/(proca_omega*A_rk)

              k13 = delta*J1_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k23 = delta*J2_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k33 = delta*J3_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k43 = delta*J4_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk     = A0     + k13
              alpha_rk = alpha0 + k23

              procaF_rk = procaF0 + k33
              procaA_rk = procaA0 + k43

!             Sources.

              procaE_rk = - cproca_mass**2*alpha_rk*procaA_rk/(proca_omega*A_rk)

              k14 = delta*J1_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k24 = delta*J2_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k34 = delta*J3_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k44 = delta*J4_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)

!             Advance variables to next grid point.

              A_g(l,i)     = A0     + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i) = alpha0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

              procaF_g(l,i) = procaF0 + (k31 + 2.d0*(k32 + k33) + k34)/6.d0
              procaA_g(l,i) = procaA0 + (k41 + 2.d0*(k42 + k43) + k44)/6.d0

!             Check if solution is blowing up. This helps to reduce
!             the need for a very fine tuned initial guess.

              if (abs(procaF_g(l,i))>2.d0*proca_phi0) then
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

        if (abs(procaF_g(0,Nrtotal)) + abs(procaF_g(0,Nrtotal-1))<epsilon) then
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

        end do

!       Fix ghost zones.

        do i=1,ghost
           A_g(l-1,1-i)      = + A_g(l-1,i)
           alpha_g(l-1,1-i)  = + alpha_g(l-1,i)
           procaF_g(l-1,1-i) = + procaF_g(l-1,i)
           procaA_g(l-1,1-i) = - procaA_g(l-1,i)
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

!    Rescale procaF.

     procaF_g = procaF_g/alphafac


!    *******************************************
!    ***   RECONSTRUCT procaPhi AND procaE   ***
!    *******************************************

!    Having found the solution, the scalar potential
!    is recovered as: phi = F/alpha.

     procaPhi_g = procaF_g/alpha_g

!    The electric field is given by:
!
!    E = - (omega*procaA + dF/dr)/(alpha*A)
!
!      = - m**2 alpha procaA / A omega

     procaE_g = - cproca_mass**2*alpha_g*procaA_g/(A_g*omega_new)


!    ***********************************
!    ***   PROCA STAR PERTURBATION   ***
!    ***********************************

!    Here we add a gaussian perturbation to the solution for the
!    Proca scalar potential procaPhi leaving procaA unchanged.
!    We then solve again the hamiltonian constraint for the radial
!    metric A, and the Gauss constraint for the electric field procaE.
!    We also solve again the polar slicing condition for the lapse.
!    
!    The perturbation should be small, and since it must
!    be even we take the sum of two symmetric gaussians.

     if ((procagauss).and.(proca_phiR_a0/=0.d0)) then

!       Message to screen.

        print *, 'Adding gaussian perturbation to Proca star ...'
        print *

!       Rescale the amplitude of the gaussian with max(phi).

        aux = maxval(abs(procaPhi_g))
        proca_phiR_a0 = aux*proca_phiR_a0

!       Initialize again A, alpha, E.

        A_g      = 1.d0
        alpha_g  = 1.d0
        procaE_g = 0.d0

!       Rescale back procaF (otherwise it won't be consistent any more).

        procaF_g = procaF_g*alphafac

!       Add perturbation to proca_Phi.

        if (proca_phiR_r0==0.d0) then
           procaF_g = procaF_g + proca_phiR_a0*exp(-rr**2/proca_phiR_s0**2)
        else
           procaF_g = procaF_g + proca_phiR_a0 &
                    *(exp(-(rr-proca_phiR_r0)**2/proca_phiR_s0**2) &
                    + exp(-(rr+proca_phiR_r0)**2/proca_phiR_s0**2))
        end if

!       Loop over grid levels. We solve from fine to coarse grid.

        do l=Nl-1,0,-1

!          Find initial point. Only the finest grid
!          integrates from the origin.

           if (l==Nl-1) then
              imin = 1
           else
              imin = Nrtotal/2
           end if

!          For coarse grids we interpolate the initial point.

           if (l<Nl-1) then

              A_g(l,imin-1) = (9.d0*(A_g(l+1,Nrtotal-2)+A_g(l+1,Nrtotal-3)) &
                            - (A_g(l+1,Nrtotal-4)+A_g(l+1,Nrtotal-1)))/16.d0
              alpha_g(l,imin-1) = (9.d0*(alpha_g(l+1,Nrtotal-2)+alpha_g(l+1,Nrtotal-3)) &
                            - (alpha_g(l+1,Nrtotal-4)+alpha_g(l+1,Nrtotal-1)))/16.d0

              procaE_g(l,imin-1) = (9.d0*(procaE_g(l+1,Nrtotal-2)+procaE_g(l+1,Nrtotal-3)) &
                                 - (procaE_g(l+1,Nrtotal-4)+procaE_g(l+1,Nrtotal-1)))/16.d0

           end if

!          Fourth order Runge-Kutta.

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

!                Values of (procaF,procaA,procaE) at origin.

                 procaF0 = (9.d0*(procaF_g(l,0)+procaF_g(l,1)) - (procaF_g(l,-1)+procaF_g(l,2)))/16.d0
                 procaA0 = 0.d0
                 procaE0 = 0.d0

!             Grid spacing and values at previous grid point.

              else

                 delta = dr(l)
                 r0    = rr(l,i-1)

                 A0     = A_g(l,i-1)
                 alpha0 = alpha_g(l,i-1)

                 procaF0 = procaF_g(l,i-1)
                 procaA0 = procaA_g(l,i-1)
                 procaE0 = procaE_g(l,i-1)

              end if

!             I) First Runge-Kutta step.

!             Sources at first grid point if we start
!             from the origin (for finer grid).

              if (i==1) then

!                At the origin we have:  A' = alpha' = 0.

                 k11 = 0.d0
                 k21 = 0.d0

!                If we take E ~ k r with k constant close to
!                the origin we find from the Gauss constraint:
!
!                k = - m^2 Phi(r=0) / 3.

                 k51 = - delta*cproca_mass**2*procaF0/3.d0

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk     = A0
                 alpha_rk = alpha0

                 procaF_rk = procaF0
                 procaA_rk = procaA0
                 procaE_rk = procaE0

!                Sources.

                 k11 = delta*J1_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
                 k21 = delta*J2_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
                 k51 = delta*J5_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk      = A0      + half*k11
              alpha_rk  = alpha0  + half*k21
              procaE_rk = procaE0 + half*k51

              if (i==1) then  ! Linear interpolation for first point.
                 procaF_rk = 0.5d0*(procaF0 + procaF_g(l,1))
                 procaA_rk = 0.5d0*(procaA0 + procaA_g(l,1))
              else            ! Cubic interpolation for the rest.
                 procaF_rk = (9.d0*(procaF_g(l,i)+procaF_g(l,i-1)) - (procaF_g(l,i-2)+procaF_g(l,i+1)))/16.d0
                 procaA_rk = (9.d0*(procaA_g(l,i)+procaA_g(l,i-1)) - (procaA_g(l,i-2)+procaA_g(l,i+1)))/16.d0
              end if

!             Sources.

              k12 = delta*J1_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k22 = delta*J2_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k52 = delta*J5_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk      = A0      + half*k12
              alpha_rk  = alpha0  + half*k22
              procaE_rk = procaE0 + half*k52

!             Sources.

              k13 = delta*J1_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k23 = delta*J2_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k53 = delta*J5_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk      = A0      + k13
              alpha_rk  = alpha0  + k23
              procaE_rk = procaE0 + k53

              procaF_rk = procaF_g(l,i)
              procaA_rk = procaA_g(l,i)

!             Sources.

              k14 = delta*J1_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k24 = delta*J2_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)
              k54 = delta*J5_proca(A_rk,alpha_rk,procaF_rk,procaA_rk,procaE_rk,rm)

!             Advance variables to next grid point.

              A_g(l,i)      = A0      + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i)  = alpha0  + (k21 + 2.d0*(k22 + k23) + k24)/6.d0
              procaE_g(l,i) = procaE0 + (k51 + 2.d0*(k52 + k53) + k54)/6.d0

           end do

!          Fix ghost zones.

           do i=1,ghost
              A_g(l,1-i)      = + A_g(l,i)
              alpha_g(l,1-i)  = + alpha_g(l,i)
              procaE_g(l,1-i) = - procaE_g(l,i)
           end do

        end do

!       Restrict solution from fine to coarse grid.

        do l=Nl-1,1,-1

           do i=1,Nrtotal-ghost,2

              iaux = i/2 + 1
              rm = rr(l-1,iaux)

              A_g(l-1,iaux)      = (9.d0*(A_g(l,i)+A_g(l,i+1)) - (A_g(l,i-1)+A_g(l,i+2)))/16.d0
              alpha_g(l-1,iaux)  = (9.d0*(alpha_g(l,i)+alpha_g(l,i+1)) - (alpha_g(l,i-1)+alpha_g(l,i+2)))/16.d0
              procaE_g(l-1,iaux) = (9.d0*(procaE_g(l,i)+procaE_g(l,i+1)) - (procaE_g(l,i-1)+procaE_g(l,i+2)))/16.d0

           end do

!          Fix ghost zones.

           do i=1,ghost
              A_g(l-1,1-i)      = + A_g(l-1,i)
              alpha_g(l-1,1-i)  = + alpha_g(l-1,i)
              procaE_g(l-1,1-i) = - procaE_g(l-1,i)
           end do

        end do

!       Rescale lapse and F again.

        if (order=="two") then
           alphafac = alpha_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
               *(3.d0*alpha_g(0,Nrtotal) - 4.d0*alpha_g(0,Nrtotal-1) + alpha_g(0,Nrtotal-2))
        else
           alphafac = alpha_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
               *(25.d0*alpha_g(0,Nrtotal) - 48.d0*alpha_g(0,Nrtotal-1) &
               + 36.d0*alpha_g(0,Nrtotal-2) - 16.d0*alpha_g(0,Nrtotal-3) + 3.d0*alpha_g(0,Nrtotal-4))/3.d0
        end if

        alpha_g  = alpha_g/alphafac
        procaF_g = procaF_g/alphafac

!       Find perturbed Phi = F/alpha.

        procaPhi_g = procaF_g/alpha_g

     end if


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

  end if


! ************************************************************************
! ***   IMAGINARY PARTS OF (procaPhi,procaE) AND REAL PART OF procaA   ***
! ************************************************************************

! Set imaginary parts of procaPhi and procaE to zero.

  cprocaPhi_I = 0.d0
  cprocaE_I   = 0.d0

! Set real part of procaA to zero.

  cprocaA_R = 0.d0


! ***************
! ***   END   ***
! ***************

  end subroutine idata_Procastar








! **********************************
! ***   RADIAL DERIVATIVE OF A   ***
! **********************************

! The radial derivative of A comes from the
! Hamiltonian constraint.

  function J1_proca(A,alpha,procaF,procaA,procaE,rm)

  use param

  implicit none

  real(8) J1_proca
  real(8) A,alpha,procaF,procaA,procaE,rm
  real(8) rho

!                           2     2           2    2
! rho  =  + 1 / (8 pi) { A E  +  m [ (F/alpha)  + a / A ] ]
!
! Notice that we don't divide by 8*pi since it cancels.

  rho = A*procaE**2 + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A)

! dA/dr  =  A [ (1 - A)/r + 8 pi r A rho ]
!
! Notice that we don't multiply the last term with 8*pi since it cancels.

  J1_proca = A*((1.d0-A)/rm + rm*A*rho)

  end function J1_proca








! **************************************
! ***   RADIAL DERIVATIVE OF ALPHA   ***
! **************************************

! The radial derivative of alpha comes from
! the polar slicing condition:  dK_{theta,theta}/dt=0.

  function J2_proca(A,alpha,procaF,procaA,procaE,rm)

  use param

  implicit none

  real(8) J2_proca
  real(8) A,alpha,procaF,procaA,procaE,rm
  real(8) SA

!                           2     2           2    2
! SA   =  - 1 / (8 pi) { A E  -  m [ (F/alpha)  + a / A ] ]
!
!
! Notice that we don't divide by 8*pi since it cancels. 

  SA = - A*procaE**2 + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A)

! dalpha/dr  =  alpha [ (A - 1)/(2r) + 4 pi r A SA ]
!
! Notice that we don't multiply the last term with 4*pi*A since
! it cancels (we are left only with 1/2).

  J2_proca = alpha*(0.5d0*(A-1.d0)/rm + 0.5d0*rm*A*SA)

  end function J2_proca








! ***************************************
! ***   RADIAL DERIVATIVE OF procaF   ***
! ***************************************

! The radial derivative of procaF from the
! Proca evolution equations.

  function J3_proca(A,alpha,procaF,procaA,procaE,rm)

  use param

  implicit none

  real(8) J3_proca
  real(8) A,alpha,procaF,procaA,procaE,rm

!                                       2
! dF/dr  =  a omega [ ( m alpha / omega )  -  1 ]

  J3_proca = procaA*proca_omega*((cproca_mass*alpha/proca_omega)**2 - 1.d0)

  end function J3_proca








! ***************************************
! ***   RADIAL DERIVATIVE OF procaA   ***
! ***************************************

! The radial derivative of procaA from the Proca
! evolution equations.

  function J4_proca(A,alpha,procaF,procaA,procaE,rm)

  use param

  implicit none

  real(8) J4_proca
  real(8) A,alpha,procaF,procaA,procaE,rm
  real(8) aux

!                   2
! SA - rho  = -  A E / 4 pi
!
! Notice that we don't divide by 4*pi since it cancels.

  aux = - A*procaE**2

!                            2
! da/dr  =  omega F A / alpha  -  a [ (A+1)/r + 4 pi r A (SA - rho) ]
!
! Notice that we don't multiply the last term with 4*pi*A since
! it cancels.

  J4_proca = proca_omega*procaF*A/alpha**2 - procaA*((A + 1.d0)/rm + rm*A*aux)

  end function J4_proca







! ***************************************
! ***   RADIAL DERIVATIVE OF procaE   ***
! ***************************************

! The radial derivative of procaE. From the
! Gauss constraint.

  function J5_proca(A,alpha,procaF,procaA,procaE,rm)

  use param

  implicit none

  real(8) J5_proca
  real(8) A,alpha,procaF,procaA,procaE,rm
  real(8) rho,aux

! We use the fact that:
!
! dA/dr  =  A [ (1 - A)/r + 8 pi r A rho ]
!
! with rho given by:
!                           2     2           2    2
! rho  =  + 1 / (8 pi) { A E  +  m [ (F/alpha)  + a / A ] ]

  rho = A*procaE**2 + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A)
  aux = (1.d0-A)/rm + rm*A*rho

! Gauss constraint:
!
!                                         2
! dE/dr  =  - E [ 2/r + (1/2A) dA/dr ] - m  Phi

  J5_proca = - procaE*(2.d0/rm + 0.5d0*aux) - cproca_mass**2*procaF/alpha

  end function J5_proca

