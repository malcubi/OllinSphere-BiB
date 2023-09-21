!$Header: $

  subroutine idata_diracstar

! *******************************************************
! ***   DIRAC STAR INITIAL DATA IN POLAR-AREAL GAUGE  ***
! *******************************************************

! This subroutine calculates initial data for a Dirac star
! using a shooting method in the polar-areal gauge (PA).
! Dirac stars are solutions such that the spacetime is static
! and the Dirac fields have a harmonic dependence on time.

! MAIN PROBLEM: The main problem with this routine is that
! it requires a very fine tuned initial guess since we are
! fighting a growing exponential.  This gets worse as
! we move the boundary further away. So the way to
! proceed is to put the boundary close in, find omega,
! and then move it slowly further out adjusting omega
! to high precision as we move out.

! To obtain the initial data we assume that spacetime is
! static (K_ij=0), and also that the Dirac fields have
! the form:
!
! F(t,r) = f(r) exp(-i omega t)
! G(t,r) = i g(r) exp(-i omega t)
!
! Notice that with this ansatz the stress-energy tensor
! is time-independent so the metric can be static.
!
! The standard ansatz for the metric for boson stars is:
!
!   2          2   2        2     2      2
! ds  = - alpha  dt  +  A dr  +  r dOmega
!
! In other words, we are assuming psi=B=1, beta=0, and
! the radial coordinate is the area radius.
!
! We substitute the above metric and this ansatz into the
! Hamiltonian constraint to find:
!
!
! dA/dr  =  A [ (1 - A)/r + 8 pi r A rho ]
!
!
! with rho the energy density of the Dirac field which
! now becomes:
!
! rho  =  omega / (2 pi alpha)  [ f^2 + g^2 ]
!
! For the lapse we use the polar slicing condition
! K_{theta,theta}=0, which implies:
!
!
! dalpha/dr  =  alpha [ (A - 1)/2r + 4 pi r A SA ]
!
!
! with SA given by:
!
! SA  =  rho - 1/pi [ f g / r + m/2 (f^2 - g^2) ]
!
! Notice that for a static solution this should be equivalent
! to maximal slicing but is easier to solve.
!
! On the other hand, substituting our ansatz in the Dirac
! equations for (F,G) we find:
!
!                                              1/2             1/2
! dF/dr  =  - F [ d alpha / 2 alpha  +  ( 1 - A   ) / r ]  -  A   G ( m + omega / alpha )
!                  r
!
!                                              1/2             1/2
! dG/dr  =  - G [ d alpha / 2 alpha  +  ( 1 + A   ) / r ]  -  A   F ( m - omega / alpha )
!                  r
!
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
  real(8) A0,alpha0,F0,G0               ! Initial values of variables.
  real(8) A_rk,alpha_rk,F_rk,G_rk       ! Runge-Kutta values of variables.
  real(8) k11,k12,k13,k14               ! Runge-Kutta sources for A.
  real(8) k21,k22,k23,k24               ! Runge-Kutta sources for alpha.
  real(8) k31,k32,k33,k34               ! Runge-Kutta sources for F.
  real(8) k41,k42,k43,k44               ! Runge-Kutta sources for G.
  real(8) J1_DIR,J2_DIR,J3_DIR,J4_DIR   ! Functions for sources of differential equations.
  real(8) res,res_old                   ! Residual.
  real(8) omega_new,omega_old,domega    ! Trial frequency and frequency interval.
  real(8) half,smallpi                  ! Numbers.
  real(8) rm,alphafac,aux               ! Auxiliary quantities.
  real(8) :: epsilon = 1.d-8            ! Tolerance.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                 ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g,alpha_g        ! Radial metric and lapse global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: F_g,G_g            ! Global arrays for Dirac functions (f,g).


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  smallpi = acos(-1.d0)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for a Dirac star in polar-areal gauge using shooting method ...'
     print *
  end if


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

  A_g     = 1.d0
  alpha_g = 1.d0

  F_g = 0.d0
  G_g = 0.d0


! *********************************************
! ***   ONLY PROCESSOR 0 SOLVES THE ODE's   ***
! *********************************************

! Rescale the tolerance with dirac_f0, so the final tolerance
! is measured with respect to this value.

  epsilon = epsilon*dirac_f0

  if (rank==0) then


!    **********************
!    ***   INITIALIZE   ***
!    **********************

!    Find initial domega.

     domega = (omega_right - omega_left)/20.d0

!    Sanity check.

     if (domega<0.d0) then
        print *
        print *, 'For Dirac star initial data we must have omega_right>omega_left.'
        print *, 'Aborting! (subroutine idata_diracstar)'
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

     dirac_omega = omega_left

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

              F_g(l,imin-1)     = (9.d0*(F_g(l+1,Nrtotal-2)+F_g(l+1,Nrtotal-3)) &
                                - (F_g(l+1,Nrtotal-4)+F_g(l+1,Nrtotal-1)))/16.d0

              G_g(l,imin-1)     = (9.d0*(G_g(l+1,Nrtotal-2)+G_g(l+1,Nrtotal-3)) &
                                - (G_g(l+1,Nrtotal-4)+G_g(l+1,Nrtotal-1)))/16.d0

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

!                Values of (F,G) at origin.

                 F0 = dirac_f0
                 G0 = 0.d0
 
!             Grid spacing and values at previous grid point.

              else

                 delta  = dr(l)
                 r0     = rr(l,i-1)

                 A0     = A_g(l,i-1)
                 alpha0 = alpha_g(l,i-1)

                 F0     = F_g(l,i-1)
                 G0     = G_g(l,i-1)

              end if

!             I) First Runge-Kutta step.

!             Sources at first grid point if we start
!             from the origin (for finer grid).

              if (i==1) then

!                At the origin we have:  A' = alpha' = F' = 0.

                 k11 = 0.d0
                 k21 = 0.d0
                 k31 = 0.d0

!                For G at the origin we have:  G' = f0 (omega - m)/3.
!                (Remember that k41 must be multiplied with delta.)

                 k41 = delta*dirac_f0*(dirac_omega - dirac_mass)/3.d0

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk     = A0
                 alpha_rk = alpha0

                 F_rk = F0
                 G_rk = G0

!                Sources.

                 k11 = delta*J1_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
                 k21 = delta*J2_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
                 k31 = delta*J3_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
                 k41 = delta*J4_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk     = A0     + half*k11
              alpha_rk = alpha0 + half*k21
              F_rk     = F0     + half*k31
              G_rk     = G0     + half*k41

!             Sources.

              k12 = delta*J1_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
              k22 = delta*J2_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
              k32 = delta*J3_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
              k42 = delta*J4_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk     = A0     + half*k12
              alpha_rk = alpha0 + half*k22
              F_rk     = F0     + half*k32
              G_rk     = G0     + half*k42

!             Sources.

              k13 = delta*J1_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
              k23 = delta*J2_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
              k33 = delta*J3_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
              k43 = delta*J4_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk     = A0     + k13
              alpha_rk = alpha0 + k23
              F_rk     = F0     + k33
              G_rk     = G0     + k43

!             Sources.

              k14 = delta*J1_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
              k24 = delta*J2_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
              k34 = delta*J3_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)
              k44 = delta*J4_DIR(A_rk,alpha_rk,F_rk,G_rk,rm)

!             Advance variables to next grid point.

              A_g(l,i)     = A0     + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i) = alpha0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0
              F_g(l,i)     = F0     + (k31 + 2.d0*(k32 + k33) + k34)/6.d0
              G_g(l,i)     = G0     + (k41 + 2.d0*(k42 + k43) + k44)/6.d0

!             Check if solution is blowing up. This helps to reduce
!             the need for a very fine tuned initial guess.

              if (abs(F_g(l,i))>2.d0*dirac_f0) then
                 if (.not.left) then
                    omega_left = omega_left + domega
                    dirac_omega = omega_left
                    goto 100
                 else if (.not.right) then
                    omega_right = omega_right - domega
                    dirac_omega = omega_right
                    goto 100
                 end if
              end if

!             Check if solution is already very small.

              if (abs(F_g(l,i))+abs(F_g(l,i-1))<epsilon) then
                 F_g(l,i) = 0.d0
                 G_g(l,i) = 0.d0
              end if

           end do


!          ***********************
!          ***   GHOST ZONES   ***
!          ***********************

           do i=1,ghost

              A_g(l,1-i)     = + A_g(l,i)
              alpha_g(l,1-i) = + alpha_g(l,i)

              F_g(l,1-i)     = + F_g(l,i)
              G_g(l,1-i)     = - G_g(l,i)

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
!       But in principle we don't know the value of k.  However, we can obtain
!       it as: k = - log|F|/r.  So the final condition is:
!
!       F' - ( log|F| / r ) phi = 0
!
!       Also, if the solution at the last two points is already smaller than
!       the tolerance we just assume we have found the solution and jump out.

        res_old = res

        if (abs(F_g(0,Nrtotal))+abs(F_g(0,Nrtotal-1))<epsilon) then
           res = epsilon/2.d0
           goto 100
        else
           res = (F_g(0,Nrtotal)-F_g(0,Nrtotal-1))/dr(0) &
               - log(abs(F_g(0,Nrtotal)))/rr(0,Nrtotal)*F_g(0,Nrtotal)
        end if

!       Secant method:  Having found the difference for the two values of omega
!       we can use it to extrapolate linearly to the next best guess.

        if (.not.left) then

           left = .true.

           omega_old = dirac_omega
           dirac_omega = omega_right

        else

           if (.not.right) right = .true.

           aux = res - res_old

           if ((aux/=0.d0).and.(abs(aux)<1.d-1)) then
              omega_new = dirac_omega - res*(dirac_omega - omega_old)/aux
           else
              omega_new = 0.5d0*(omega_old + dirac_omega)
           end if

           omega_old = dirac_omega
           dirac_omega = omega_new

        end if

100     continue


!       **************************
!       ***   END ITERATIONS   ***
!       **************************

!       Output data to screen.

        write(*,"(A,I4,A,ES22.16,A,ES9.2)") ' Iteration: ',iter,'    Frequency: ',dirac_omega,'    Residual: ',res

     end do

!    Message in case we reached the maximum number of iterations.

     if ((iter>=maxiter).and.(maxiter/=1)) then
        print *
        print *, 'Maximum iteration number reached, initial data solver did not converge.'
        print *, 'Aborting! (subroutine idata_diracstar)'
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

           A_g(l-1,iaux)     = (9.d0*(A_g(l,i)+A_g(l,i+1))         - (A_g(l,i-1)+A_g(l,i+2)))/16.d0
           alpha_g(l-1,iaux) = (9.d0*(alpha_g(l,i)+alpha_g(l,i+1)) - (alpha_g(l,i-1)+alpha_g(l,i+2)))/16.d0
           F_g(l-1,iaux)     = (9.d0*(F_g(l,i)+F_g(l,i+1))         - (F_g(l,i-1)+F_g(l,i+2)))/16.d0
           G_g(l-1,iaux)     = (9.d0*(G_g(l,i)+G_g(l,i+1))         - (G_g(l,i-1)+G_g(l,i+2)))/16.d0

        end do

!       Fix ghost zones.

        do i=1,ghost
           A_g(l-1,1-i)     = + A_g(l-1,i)
           alpha_g(l-1,1-i) = + alpha_g(l-1,i)
           F_g(l-1,1-i)     = + F_g(l-1,i)
           G_g(l-1,1-i)     = - G_g(l-1,i)
        end do

     end do


!    *********************************
!    ***   RESCALE (ALPHA,OMEGA)   ***
!    *********************************

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

     omega_new = dirac_omega/alphafac

     if (rank==0) then
        write(*,'(A,E22.16)') ' Omega (not-rescaled) = ', dirac_omega
        write(*,'(A,E22.16)') ' Omega (rescaled)     = ', omega_new
        print *
     end if


!    ***********************************
!    ***   DIRAC STAR PERTURBATION   ***
!    ***********************************

!    NOT YET IMPLEMENTED.


! *************************************
! ***   FINISHED FINDING SOLUTION   ***
! *************************************

! When we get here the solution has been found,
! so we close the "if" statement for processor 0.

  end if


! ************************************************
! ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
! ************************************************

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! arrays with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  if (size==1) then

     A     = A_g
     alpha = alpha_g

     dirac_FR = F_g
     dirac_GI = G_g

  else

     call distribute(0,Nl-1,A,A_g)
     call distribute(0,Nl-1,alpha,alpha_g)

     call distribute(0,Nl-1,dirac_FR,F_g)
     call distribute(0,Nl-1,dirac_GI,F_g)

  end if

! Derivatives of A and alpha.

  diffvar => A

  do l=0,Nl-1
     D1_A(l,:) = diff1(l,+1)
  end do

  diffvar => alpha

  do l=0,Nl-1
     D1_alpha(l,:) = diff1(l,+1)
  end do


! **************************************************
! ***   IMAGINARY PART OF F AND REAL PART OF G   ***
! **************************************************

! Set the imaginary part of F and real partr of G to zero.

  dirac_FI = 0.d0
  dirac_GR = 0.d0


! ***************
! ***   END   ***
! ***************

  end subroutine idata_diracstar







! **********************************
! ***   RADIAL DERIVATIVE OF A   ***
! **********************************

! The radial derivative of A comes from the
! Hamiltonian constraint.

  function J1_DIR(A,alpha,F,G,rm)

  use param

  implicit none

  real(8) J1_DIR
  real(8) A,alpha,F,G,rm
  real(8) rho
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Energy density:
!
! rho  =  omega / (2 pi alpha) (F^2 + G^2)

  rho = dirac_omega/alpha*(F**2 + G**2)

! dA/dr = A [ (1-A)/r + 8 pi r A rho ]
!
! Notice that we already canceled a factor of 2*pi.

  J1_DIR = A*((1.d0 - A)/rm + 4.d0*rm*A*rho)

  end function J1_DIR







! **************************************
! ***   RADIAL DERIVATIVE OF ALPHA   ***
! **************************************

  function J2_DIR(A,alpha,F,G,rm)

  use param

  implicit none

  real(8) J2_DIR
  real(8) A,alpha,F,G,rm
  real(8) rho,SA
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Energy density:
!
! rho  =  omega / (2 pi alpha) (F^2 + G^2)

  rho = dirac_omega/(2.d0*smallpi*alpha)*(F**2 + G**2)

! SA = rho - 1/pi [ f g / r + m/2 (f^2 - g^2) ]

  SA = rho - (F*G/rm + 0.5d0*dirac_mass*(F**2 - G**2))/smallpi

! dalpha/dr = alpha [ (A-1)/2r + 4 pi r A SA ]

  J2_DIR = alpha*(0.5d0*(A-1.d0)/rm + 4.d0*smallpi*rm*A*SA)

  end function J2_DIR







! **********************************
! ***   RADIAL DERIVATIVE OF F   ***
! **********************************

  function J3_DIR(A,alpha,F,G,rm)

  use param

  implicit none

  real(8) J3_DIR
  real(8) A,alpha,F,G,rm,Dalpha
  real(8) rho,SA
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Energy density:
!
! rho  =  omega / (2 pi alpha) (F^2 + G^2)

  rho = dirac_omega/(2.d0*smallpi*alpha)*(F**2 + G**2)

! SA = rho - 1/pi [ f g / r + m/2 (f^2 - g^2) ]

  SA = rho - (F*G/rm + 0.5d0*dirac_mass*(F**2 - G**2))/smallpi

! dalpha/dr = alpha [ (A-1)/2r + 4 pi r A SA ]

  Dalpha = alpha*(0.5d0*(A-1.d0)/rm + 4.d0*smallpi*rm*A*SA)

!                                              1/2             1/2
! dF/dr  =  - F [ d alpha / 2 alpha  +  ( 1 - A   ) / r ]  -  A   G ( m + omega / alpha )
!                  r

  J3_DIR = - F*(Dalpha/(2.d0*alpha) + (1.d0 - sqrt(A))/rm) &
         - sqrt(A)*G*(dirac_mass + dirac_omega/alpha)

  end function J3_DIR







! **********************************
! ***   RADIAL DERIVATIVE OF G   ***
! **********************************

  function J4_DIR(A,alpha,F,G,rm)

  use param

  implicit none

  real(8) J4_DIR
  real(8) A,alpha,F,G,rm,Dalpha
  real(8) rho,SA
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Energy density:
!
! rho  =  omega / (2 pi alpha) (F^2 + G^2)

  rho = dirac_omega/(2.d0*smallpi*alpha)*(F**2 + G**2)

! SA = rho - 1/pi [ f g / r + m/2 (f^2 - g^2) ]

  SA = rho - (F*G/rm + 0.5d0*dirac_mass*(F**2 - G**2))/smallpi

! dalpha/dr = alpha [ (A-1)/2r + 4 pi r A SA ]

  Dalpha = alpha*(0.5d0*(A-1.d0)/rm + 4.d0*smallpi*rm*A*SA)

!                                              1/2             1/2
! dG/dr  =  - G [ d alpha / 2 alpha  +  ( 1 + A   ) / r ]  -  A   F ( m - omega / alpha )
!                  r

  J4_DIR = - G*(Dalpha/(2.d0*alpha) + (1.d0 + sqrt(A))/rm) &
         - sqrt(A)*F*(dirac_mass - dirac_omega/alpha)

  end function J4_DIR

