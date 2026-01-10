
  subroutine idata_ChargedDiracstar

! ****************************************************************
! ***   CHARGED DIRAC STAR INITIAL DATA IN POLAR AREAL GAUGE   ***
! ****************************************************************

! This subroutine calculates initial data for a charged
! Dirac star using a shooting method in the polar areal gauge.

! MAIN PROBLEM: The main problem with this routine is that
! it requires a very fine tuned initial guess since we are
! fighting a growing exponential.  This gets worse as
! we move the boundary further away. So the way to
! proceed is to put the boundary close in, find omega,
! and then move it slowly further out adjusting omega
! to high precision as we move out.

! Dirac stars are self-gravitating solutions of a Dirac
! spinor field such that the spacetime is static and the
! Dirac field has a harmonic time dependence.
!
! To obtain the initial data we assume that spacetime is
! static (K_ij=0), and also that the Dirac fields have
! the form:
!
! F(t,r) = f(r) exp(-i omega t)
! G(t,r) = i g(r) exp(-i omega t)
!
! In particular, at t=0 F is purely real while G is
! purely imaginary.
!
! Notice that with this ansatz the stress-energy tensor
! is time-independent so the metric can be static.
!
! The standard ansatz for the metric for Dirac stars is:
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
! with rho the total energy density given by the Dirac and
! Maxwell contributions:
!
!              2   2                               2
! rho  =  W [ f + g ] / (2 pi alpha)  +  A maxwellE / 8 pi
!
!
! where we have W := omega - q maxwellF, with maxwellF := alpha ePhi.
!
! For the lapse we use the polar slicing condition
! K_{theta,theta}=0, which implies:
!
!
! dalpha/dr  =  alpha [ (A - 1)/2r + 4 pi r A SA ]
!
!
! with SA now given by:
!
!                                      2   2
! SA  =  rho - 1/pi [ f g / r + m/2 ( f - g ) ]
!
! In this case there is no contribution from the charge.
! Notice that for a static solution this lapse condition should
! be equivalent to maximal slicing but is easier to solve.
!
! On the other hand, substituting our ansatz in the Dirac
! equations for (F,G) we find:
!
!                                              1/2             1/2
! dF/dr  =  - F [ d alpha / 2 alpha  +  ( 1 - A   ) / r ]  -  A   G ( m + W / alpha )
!                  r
!
!                                              1/2             1/2
! dG/dr  =  - G [ d alpha / 2 alpha  +  ( 1 + A   ) / r ]  -  A   F ( m - W / alpha )
!                  r
!
! with W the same as above.
!
! The equation for Maxwell electric field comes from the Gauss
! constraint, and takes the form:
!
!
! dmaxwellE/dr =  - E [ 2/r + (dA/dr) / 2A ] + 4 pi echarge
!
!
! with "echarge" the charge density of the Dirac field given by:
!
!                  2     2
! echarge  =  q ( F  +  G ) / 2 pi
!
!
! Finally, for the Maxwell scalar potential maxwellF := alpha ePhi
! we have:
!
!
! dmaxwellF/dr  =  - alpha A maxwellE
!
!
! For the boundary conditions at the origin we take:
!
! A(r=0)        = 1,           d A(r=0)     = 0
!                               r
!
! alpha(r=0)    = 1,           d alpha(r=0) = 0
!                               r
!
! diracF(r=0)   = dirac_f0,    d diracF = 0
!                               r
!
! diracG(r=0)   = 0,           d diracG(r=0) = f0 (W - m) / 3
!                               r
!
! maxwellF(r=0) = 0,           d maxwellF(r=0) = 0
!                               r
!
! maxwellE(r=0) = 0
!
! Notice that in the final solution we don't want alpha(r=0)=1,
! but rather alpha=1 at infinity.  But this is no problem as the
! slicing condition above is linear in alpha so we can always
! just rescale the lapse at the end, but in order not to affect
! the solution we must also rescale the final value of the frequency
! omega and maxwellF by the same factor.
!
! Also, we don't really want maxwellF(r=0)=0, but rather maxwellF=0
! at infinity. We can fix this by making a gauge transformation of
! the field once we found the solution that also changes the frequency
! (see comment in the code below for a detailed explanation).


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
  real(8) maxwellF0,maxwellE0           ! Initial values of Maxwell variables.

  real(8) A_rk,alpha_rk                 ! Runge-Kutta values of A and alpha.
  real(8) F_rk,G_rk                     ! Runge-Kutta values of Dirac variables.
  real(8) maxwellF_rk,maxwellE_rk       ! Runge-Kutta values of Maxwell variables.
  real(8) DF_rk,DG_rk                   ! Radial derivatives of F and G, for perturbations.
  real(8) rho_rk                        ! Energy density, for perturbations.

  real(8) k11,k12,k13,k14               ! Runge-Kutta sources for A.
  real(8) k21,k22,k23,k24               ! Runge-Kutta sources for alpha.
  real(8) k31,k32,k33,k34               ! Runge-Kutta sources for F.
  real(8) k41,k42,k43,k44               ! Runge-Kutta sources for G.
  real(8) k51,k52,k53,k54               ! Runge-Kutta sources for maxwellE.
  real(8) k61,k62,k63,k64               ! Runge-Kutta sources for maxwellF.

  real(8) J1_CDIR,J2_CDIR,J3_CDIR       ! Functions for sources of differential equations.
  real(8) J4_CDIR,J5_CDIR,J6_CDIR       ! Functions for sources of differential equations.
  real(8) J7_CDIR,J8_CDIR,J9_CDIR       ! Functions for sources of differential equations.

  real(8) res,res_old                   ! Residual.
  real(8) omega_new,omega_old,domega    ! Trial frequency and frequency interval.

  real(8) half,smallpi                  ! Numbers.
  real(8) rm,alphafac,Ffac,aux          ! Auxiliary quantities.
  real(8) :: epsilon = 1.d-8            ! Tolerance.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                        ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g,alpha_g               ! Radial metric and lapse global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: F_g,G_g                   ! Global arrays for Dirac functions (f,g).
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
     print *, 'Solving initial data for a charged Dirac star in polar-areal gauge using shooting method ...'
     print *
  end if

! Sanity check.

  if (spacetime=="minkowski") then
     print *, 'Charged Dirac star initial data is not compatible with a Minkowski background ...'
     print *, 'Aborting! (subroutine idata_chargeddiracstar)'
     print *
     call die
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

  maxwellF_g   = 0.d0
  maxwellE_g   = 0.d0
  maxwellPhi_g = 0.d0


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

              F_g(l,imin-1) = (9.d0*(F_g(l+1,Nrtotal-2)+F_g(l+1,Nrtotal-3)) &
                            - (F_g(l+1,Nrtotal-4)+F_g(l+1,Nrtotal-1)))/16.d0

              G_g(l,imin-1) = (9.d0*(G_g(l+1,Nrtotal-2)+G_g(l+1,Nrtotal-3)) &
                            - (G_g(l+1,Nrtotal-4)+G_g(l+1,Nrtotal-1)))/16.d0

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

!                Values of (F,G) at origin.

                 F0 = dirac_f0
                 G0 = 0.d0

!                Values of (maxwellF,maxwellE) at origin.

                 maxwellF0 = 0.d0
                 maxwellE0 = 0.d0

!             Grid spacing and values at previous grid point.

              else

                 delta = dr(l)
                 r0    = rr(l,i-1)

                 A0     = A_g(l,i-1)
                 alpha0 = alpha_g(l,i-1)

                 F0 = F_g(l,i-1)
                 G0 = G_g(l,i-1)

                 maxwellE0 = maxwellE_g(l,i-1)
                 maxwellF0 = maxwellF_g(l,i-1)

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

!                For maxwellE and maxwellF we have:  maxwellE' = maxwellF' = 0

                 k51 = 0.d0
                 k61 = 0.d0

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk     = A0
                 alpha_rk = alpha0

                 F_rk = F0
                 G_rk = G0

                 maxwellE_rk = maxwellE0
                 maxwellF_rk = maxwellF0

!                Sources.

                 k11 = delta*J1_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
                 k21 = delta*J2_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
                 k31 = delta*J3_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
                 k41 = delta*J4_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
                 k51 = delta*J5_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
                 k61 = delta*J6_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk     = A0     + half*k11
              alpha_rk = alpha0 + half*k21

              F_rk     = F0     + half*k31
              G_rk     = G0     + half*k41

              maxwellE_rk = maxwellE0 + half*k51
              maxwellF_rk = maxwellF0 + half*k61

!             Sources.

              k12 = delta*J1_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k22 = delta*J2_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k32 = delta*J3_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k42 = delta*J4_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k52 = delta*J5_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k62 = delta*J6_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk     = A0     + half*k12
              alpha_rk = alpha0 + half*k22

              F_rk     = F0     + half*k32
              G_rk     = G0     + half*k42

              maxwellE_rk = maxwellE0 + half*k52
              maxwellF_rk = maxwellF0 + half*k62

!             Sources.

              k13 = delta*J1_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k23 = delta*J2_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k33 = delta*J3_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k43 = delta*J4_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k53 = delta*J5_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k63 = delta*J6_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk     = A0     + k13
              alpha_rk = alpha0 + k23

              F_rk     = F0     + k33
              G_rk     = G0     + k43

              maxwellE_rk = maxwellE0 + k53
              maxwellF_rk = maxwellF0 + k63

!             Sources.

              k14 = delta*J1_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k24 = delta*J2_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k34 = delta*J3_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k44 = delta*J4_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k54 = delta*J5_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)
              k64 = delta*J6_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)

!             Advance variables to next grid point.

              A_g(l,i)     = A0     + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i) = alpha0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

              F_g(l,i)     = F0     + (k31 + 2.d0*(k32 + k33) + k34)/6.d0
              G_g(l,i)     = G0     + (k41 + 2.d0*(k42 + k43) + k44)/6.d0

              maxwellE_g(l,i) = maxwellE0 + (k51 + 2.d0*(k52 + k53) + k54)/6.d0
              maxwellF_g(l,i) = maxwellF0 + (k61 + 2.d0*(k62 + k63) + k64)/6.d0

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

              F_g(l,1-i) = + F_g(l,i)
              G_g(l,1-i) = - G_g(l,i)

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

        if (abs(F_g(0,Nrtotal))+abs(F_g(0,Nrtotal-1))<epsilon) then
           res = epsilon/2.d0
           goto 100
        else
           aux = dirac_omega - dirac_q*maxwellF_g(0,Nrtotal)
           aux = abs(dirac_mass**2 - (aux/alpha(0,Nrtotal))**2)
           res = (F_g(0,Nrtotal)-F_g(0,Nrtotal-1))/dr(0) + dsqrt(aux)*F_g(0,Nrtotal)
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

        write(*,"(A,I4,A,ES23.16,A,ES9.2)") ' Iteration: ',iter,'    Frequency: ',dirac_omega,'    Residual: ',res

     end do

!    Message in case we reached the maximum number of iterations.

     if ((iter>=maxiter).and.(maxiter/=1)) then
        print *
        print *, 'Maximum iteration number reached, initial data solver did not converge.'
        print *, 'Aborting! (subroutine idata_chargeddiracstar)'
        print *
        call die
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

           F_g(l-1,iaux) = (9.d0*(F_g(l,i)+F_g(l,i+1)) - (F_g(l,i-1)+F_g(l,i+2)))/16.d0
           G_g(l-1,iaux) = (9.d0*(G_g(l,i)+G_g(l,i+1)) - (G_g(l,i-1)+G_g(l,i+2)))/16.d0

           maxwellE_g(l-1,iaux) = (9.d0*(maxwellE_g(l,i)+maxwellE_g(l,i+1)) - (maxwellE_g(l,i-1)+maxwellE_g(l,i+2)))/16.d0
           maxwellF_g(l-1,iaux) = (9.d0*(maxwellF_g(l,i)+maxwellF_g(l,i+1)) - (maxwellF_g(l,i-1)+maxwellF_g(l,i+2)))/16.d0

        end do

!       Fix ghost zones.

        do i=1,ghost

           A_g(l-1,1-i)     = + A_g(l-1,i)
           alpha_g(l-1,1-i) = + alpha_g(l-1,i)

           F_g(l-1,1-i) = + F_g(l-1,i)
           G_g(l-1,1-i) = - G_g(l-1,i)

           maxwellF_g(l-1,1-i) = + maxwellF_g(l-1,i)
           maxwellE_g(l-1,1-i) = - maxwellE_g(l-1,i)

        end do

     end do


!    ******************************************
!    ***   RESCALE (ALPHA,OMEGA,maxwellF)   ***
!    ******************************************

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
        write(*,'(A,E23.16)') ' Omega (not-rescaled)      = ', dirac_omega
        write(*,'(A,E23.16)') ' Omega (rescaled)          = ', omega_new
     end if

!    Rescale maxwellF.

     maxwellF_g = maxwellF_g/alphafac


!    **********************************************
!    ***   OUTPUT GAUGE TRANSFORMED FREQUENCY   ***
!    **********************************************

!    Here we output the gauge transformed frequency.
!    But we don't really do the gauge transformation
!    since that should be done AFTER the perturbation
!    if there is one.

!    Figure out asymptotic value of F.

     if (order=="two") then
        Ffac = maxwellF_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
             *(3.d0*maxwellF_g(0,Nrtotal) - 4.d0*maxwellF_g(0,Nrtotal-1) &
             + maxwellF_g(0,Nrtotal-2))
     else
        Ffac = maxwellF_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
             *(25.d0*maxwellF_g(0,Nrtotal) - 48.d0*maxwellF_g(0,Nrtotal-1) &
             + 36.d0*maxwellF_g(0,Nrtotal-2) - 16.d0*maxwellF_g(0,Nrtotal-3) &
             + 3.d0*maxwellF_g(0,Nrtotal-4))/3.d0
     end if

!    Output transformed frequency.

     omega_new = omega_new + cproca_q*Ffac

     if (rank==0) then
        write(*,'(A,ES23.16)') ' Omega (gauge transformed) = ', omega_new
        print *
     end if


!    ***********************************
!    ***   DIRAC STAR PERTURBATION   ***
!    ***********************************

!    THIS SECTION IS NOT YET FINSIHED!!!  DO NOT USE IT!

!    Here we add a gaussian perturbation to the solution for the
!    real part of proca_F, or if we prefer to the real or imaginary
!    parts of proca_G.  Remember that in order to guarantee that
!    the momentum density is zero we must have proca_F purely
!    real and proca_G  purely imaginary (we could probably
!    consider more general perturbations, but at the moment we
!    only allow these ones).  The amplitudes of the perturbations
!    are rescaled with the maximum of dirac_FR and dirac_GI.
!
!    We use the same parameters to control the form of the gaussian
!    as in the routine idata_diracpulse.f90. We then solve again for
!    (A,alpha,maxwellF,maxwellE).
!
!    An important consideration is the fact that we need to calculate
!    the sources for (A,alpha) in a different way as above.  This is
!    because the energy density for the Dirac equation depends on the
!    time derivatives of (F,G).  For an unperturbed star we can calculate
!    those derivatives using the Dirac equation and the harmonic
!    ansatz.  But once we perturb the star this is not possible
!    since the harmonic ansatz is no longer valid.  What we do here
!    is take the functions (F,G) from the unperturbed solution,
!    calculate their spatial derivatives using finite differences and
!    use them calculate the energy density, which is then passed
!    directly to another version of the sources for (A,alpha).

     if ((diracgauss).and.(abs(diracFR_a0)+abs(diracGI_a0)>0.0d0)) then

!       Message to screen.

        print *, 'Charged Dirac star perturbation not yet implemented. Aborting!'
        !call die
        print *, 'Adding gaussian perturbation to charged Dirac star ...'
        print *

!       Sanity check. Remember diracF must be real and diracG must be imaginary.

        if (diracFI_a0/=0.d0) then
           print *, 'For a perturbation for a Dirac star we must have diracFI_a0=0 ...'
           print *, 'Aborting! (subroutine idata_diracpulse)'
           print *
           call die
        end if

        if (diracGR_a0/=0.d0) then
           print *, 'For a perturbation for a Dirac star we must have diracGR_a0=0 ...'
           print *, 'Aborting! (subroutine idata_diracpulse)'
           print *
           call die
        end if

!       Also, diracGI must remain zero at the origin.

        if (diracGI_r0==0.d0) then
           print *, 'For a perturbation for a Dirac star we must have diracGI_r0 non-zero ...'
           print *, 'Aborting! (subroutine idata_diracpulse)'
           print *
           call die
        end if

!       Rescale the amplitude of the perturbations.

        aux = maxval(abs(F_g))
        diracFR_a0 = aux*diracFR_a0

        aux = maxval(abs(G_g))
        diracGI_a0 = aux*diracGI_a0

!       Initialize again (A,alpha,maxwellF,maxwellE).

        A_g     = 1.d0
        alpha_g = 1.d0

        maxwellF_g = 0.d0
        maxwellE_g = 0.d0

!       Add perturbations to F_g (must be even) and G_g (must be odd).

        if (diracFR_r0==0.d0) then
           F_g = F_g + diracFR_a0*exp(-rr**2/diracFR_s0**2)
        else
           F_g = F_g + diracFR_a0 &
               *(exp(-(rr-diracFR_r0)**2/diracFR_s0**2) &
               + exp(-(rr+diracFR_r0)**2/diracFR_s0**2))
        end if

        G_g = G_g + diracGI_a0 &
            *(exp(-(rr-diracGI_r0)**2/diracGI_s0**2) &
            - exp(-(rr+diracGI_r0)**2/diracGI_s0**2))

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

              maxwellE_g(l-1,iaux) = (9.d0*(maxwellE_g(l,i)+maxwellE_g(l,i+1)) &
                                   - (maxwellE_g(l,i-1)+maxwellE_g(l,i+2)))/16.d0
              maxwellF_g(l-1,iaux) = (9.d0*(maxwellF_g(l,i)+maxwellF_g(l,i+1)) &
                                   - (maxwellF_g(l,i-1)+maxwellF_g(l,i+2)))/16.d0

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

!                Values of (F,G) at origin.

                 F0 = (9.d0*(F_g(l,0)+F_g(l,1)) - (F_g(l,-1)+F_g(l,2)))/16.d0
                 G0 = 0.d0

!                Values of (maxwellF,maxwellE) at origin.

                 maxwellF0 = 0.d0
                 maxwellE0 = 0.d0

              else

                 delta  = dr(l)
                 r0     = rr(l,i-1)

                 A0     = A_g(l,i-1)
                 alpha0 = alpha_g(l,i-1)

                 F0 = F_g(l,i-1)
                 G0 = G_g(l,i-1)

                 maxwellF0 = maxwellF_g(l,i-1)
                 maxwellE0 = maxwellE_g(l,i-1)

              end if

!             I) First Runge-Kutta step.

!             Sources at first grid point if we start
!             from the origin (for finer grid).

              if (i==1) then

!                At the origin we have:  A' = alpha' = 0.

                 k11 = 0.d0
                 k21 = 0.d0

!                For maxwellE and maxwellF we have:  maxwellE' = maxwellF' = 0

                 k51 = 0.d0
                 k61 = 0.d0

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk        = A0
                 alpha_rk    = alpha0

                 maxwellF_rk = maxwellF0
                 maxwellE_rk = maxwellE0

!                Fourth order derivatives of (F,G) at point i-1.

                 if (i==Nrtotal) then
                    DF_rk = (3.d0*F_g(l,i) + 10.d0*F_g(l,i-1) - 18.d0*F_g(l,i-2) &
                          + 6.d0*F_g(l,i-3) - F_g(l,i-4))/(12.d0*dr(l))
                    DG_rk = (3.d0*G_g(l,i) + 10.d0*G_g(l,i-1) - 18.d0*G_g(l,i-2) &
                          + 6.d0*G_g(l,i-3) - G_g(l,i-4))/(12.d0*dr(l))
                 else
                    DF_rk = (8.d0*(F_g(l,i) - F_g(l,i-2)) - (F_g(l,i+1) - F_g(l,i-3)))/(12.d0*dr(l))
                    DG_rk = (8.d0*(G_g(l,i) - G_g(l,i-2)) - (G_g(l,i+1) - G_g(l,i-3)))/(12.d0*dr(l))
                 end if

!                Sources.

                 rho_rk = 0.5d0/smallpi*((F_rk*DG_rk - G_rk*DF_rk)/sqrt(A_rk) &
                        + 2.d0*F_rk*G_rk/rm + dirac_mass*(F_rk**2 - G_rk**2)) &
                        + 0.125d0*A_rk*maxwellE_rk**2/smallpi

                 k11 = delta*J7_CDIR(A_rk,alpha_rk,F_rk,G_rk,rho_rk,rm)
                 k21 = delta*J8_CDIR(A_rk,alpha_rk,F_rk,G_rk,rho_rk,rm)

                 k51 = delta*J9_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rho_rk,rm)
                 k61 = delta*J6_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk     = A0     + half*k11
              alpha_rk = alpha0 + half*k21

              maxwellE_rk = maxwellE0 + half*k51
              maxwellF_rk = maxwellF0 + half*k61

              if (i==1) then  ! Linear interpolation for first and last points.
                 F_rk = 0.5d0*(F0 + F_g(l,i))
                 G_rk = 0.5d0*(G0 + G_g(l,i))
              else            ! Cubic interpolation for the rest.
                 F_rk = (9.d0*(F_g(l,i)+F_g(l,i-1)) - (F_g(l,i-2)+F_g(l,i+1)))/16.d0
                 G_rk = (9.d0*(G_g(l,i)+G_g(l,i-1)) - (G_g(l,i-2)+G_g(l,i+1)))/16.d0
              end if

!             Fourth order derivatives of (F,G) at point i-1/2.

              if (i==Nrtotal) then ! Second order at the boundary.
                 DF_rk = (F_g(l,i) - F_g(l,i-1))/dr(l)
                 DG_rk = (G_g(l,i) - G_g(l,i-1))/dr(l)
              else
                 DF_rk = (27.d0*(F_g(l,i)-F_g(l,i-1)) - (F_g(l,i+1) - F_g(l,i-2)))/(24.d0*dr(l))
                 DG_rk = (27.d0*(G_g(l,i)-G_g(l,i-1)) - (G_g(l,i+1) - G_g(l,i-2)))/(24.d0*dr(l))
              end if

!             Sources.

              rho_rk = 0.5d0/smallpi*((F_rk*DG_rk - G_rk*DF_rk)/sqrt(A_rk) &
                     + 2.d0*F_rk*G_rk/rm + dirac_mass*(F_rk**2 - G_rk**2)) &
                     + 0.125d0*A_rk*maxwellE_rk**2/smallpi

              k12 = delta*J7_CDIR(A_rk,alpha_rk,F_rk,G_rk,rho_rk,rm)
              k22 = delta*J8_CDIR(A_rk,alpha_rk,F_rk,G_rk,rho_rk,rm)

              k52 = delta*J9_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rho_rk,rm)
              k62 = delta*J6_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk     = A0     + half*k12
              alpha_rk = alpha0 + half*k22

              maxwellE_rk = maxwellE0 + half*k52
              maxwellF_rk = maxwellF0 + half*k62

!             Sources.

              rho_rk = 0.5d0/smallpi*((F_rk*DG_rk - G_rk*DF_rk)/sqrt(A_rk) &
                     + 2.d0*F_rk*G_rk/rm + dirac_mass*(F_rk**2 - G_rk**2)) &
                     + 0.125d0*A_rk*maxwellE_rk**2/smallpi

              k13 = delta*J7_CDIR(A_rk,alpha_rk,F_rk,G_rk,rho_rk,rm)
              k23 = delta*J8_CDIR(A_rk,alpha_rk,F_rk,G_rk,rho_rk,rm)

              k53 = delta*J9_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rho_rk,rm)
              k63 = delta*J6_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk     = A0     + k13
              alpha_rk = alpha0 + k23

              maxwellE_rk = maxwellE0 + k53
              maxwellF_rk = maxwellF0 + k63

              F_rk = F_g(l,i)
              G_rk = G_g(l,i)

!             Fourth order derivatives of (F,G) at point i.

              if (i==Nrtotal) then
                 DF_rk = (25.d0*F_g(l,i) - 48.d0*F_g(l,i-1) + 36.d0*F_g(l,i-2) &
                       - 16.d0*F_g(l,i-3) + 3.d0*F_g(l,i-4))/(12.d0*dr(l))
                 DG_rk = (25.d0*G_g(l,i) - 48.d0*G_g(l,i-1) + 36.d0*G_g(l,i-2) &
                       - 16.d0*G_g(l,i-3) + 3.d0*G_g(l,i-4))/(12.d0*dr(l))
              else if (i==Nrtotal-1) then
                 DF_rk = (3.d0*F_g(l,i+1) + 10.d0*F_g(l,i) - 18.d0*F_g(l,i-1) &
                        + 6.d0*F_g(l,i-2) - F_g(l,i-3))/(12.d0*dr(l))
                 DG_rk = (3.d0*G_g(l,i+1) + 10.d0*G_g(l,i) - 18.d0*G_g(l,i-1) &
                        + 6.d0*G_g(l,i-2) - G_g(l,i-3))/(12.d0*dr(l))
              else
                 DF_rk = (8.d0*(F_g(l,i+1) - F_g(l,i-1)) - (F_g(l,i+2) - F_g(l,i-2)))/(12.d0*dr(l))
                 DG_rk = (8.d0*(G_g(l,i+1) - G_g(l,i-1)) - (G_g(l,i+2) - G_g(l,i-2)))/(12.d0*dr(l))
              end if

!             Sources.

              rho_rk = 0.5d0/smallpi*((F_rk*DG_rk - G_rk*DF_rk)/sqrt(A_rk) &
                     + 2.d0*F_rk*G_rk/rm + dirac_mass*(F_rk**2 - G_rk**2)) &
                     + 0.125d0*A_rk*maxwellE_rk**2/smallpi

              k14 = delta*J7_CDIR(A_rk,alpha_rk,F_rk,G_rk,rho_rk,rm)
              k24 = delta*J8_CDIR(A_rk,alpha_rk,F_rk,G_rk,rho_rk,rm)

              k54 = delta*J9_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rho_rk,rm)
              k64 = delta*J6_CDIR(A_rk,alpha_rk,F_rk,G_rk,maxwellE_rk,maxwellF_rk,rm)

!             Advance variables to next grid point.

              A_g(l,i)     = A0     + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i) = alpha0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

              maxwellE_g(l,i) = maxwellE0 + (k51 + 2.d0*(k52 + k53) + k54)/6.d0
              maxwellF_g(l,i) = maxwellF0 + (k61 + 2.d0*(k62 + k63) + k64)/6.d0

           end do

!          Fix ghost zones.

           do i=1,ghost
              A_g(l,1-i)        = + A_g(l,i)
              alpha_g(l,1-i)    = + alpha_g(l,i)
              maxwellF_g(l,1-i) = + maxwellF_g(l,i)      
              maxwellE_g(l,1-i) = - maxwellE_g(l,i)
           end do

        end do

!       Restrict solution from fine to coarse grid.

        do l=Nl-1,1,-1

           do i=1,Nrtotal-ghost,2

              iaux = i/2 + 1
              rm = rr(l-1,iaux)

              A_g(l-1,iaux)     = (9.d0*(A_g(l,i)+A_g(l,i+1)) - (A_g(l,i-1)+A_g(l,i+2)))/16.d0
              alpha_g(l-1,iaux) = (9.d0*(alpha_g(l,i)+alpha_g(l,i+1)) - (alpha_g(l,i-1)+alpha_g(l,i+2)))/16.d0

              maxwellE_g(l-1,iaux) = (9.d0*(maxwellE_g(l,i)+maxwellE_g(l,i+1)) - (maxwellE_g(l,i-1)+maxwellE_g(l,i+2)))/16.d0
              maxwellF_g(l-1,iaux) = (9.d0*(maxwellF_g(l,i)+maxwellF_g(l,i+1)) - (maxwellF_g(l,i-1)+maxwellF_g(l,i+2)))/16.d0

           end do

!          Fix ghost zones.

           do i=1,ghost
              A_g(l-1,1-i)        = + A_g(l-1,i)
              alpha_g(l-1,1-i)    = + alpha_g(l-1,i)
              maxwellF_g(l-1,1-i) = + maxwellF_g(l-1,i)
              maxwellE_g(l-1,1-i) = - maxwellE_g(l-1,i)
           end do

        end do

!       Rescale lapse and maxwellF again.

        if (order=="two") then
           alphafac = alpha_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
               *(3.d0*alpha_g(0,Nrtotal) - 4.d0*alpha_g(0,Nrtotal-1) + alpha_g(0,Nrtotal-2))
        else
           alphafac = alpha_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
               *(25.d0*alpha_g(0,Nrtotal) - 48.d0*alpha_g(0,Nrtotal-1) &
               + 36.d0*alpha_g(0,Nrtotal-2) - 16.d0*alpha_g(0,Nrtotal-3) + 3.d0*alpha_g(0,Nrtotal-4))/3.d0
        end if

        alpha_g  = alpha_g/alphafac

        maxwellF_g = maxwellF_g/alphafac

!       Reconstruct maxwellPhi.

        maxwellPhi_g = maxwellF_g/alpha_g

     end if


!    ********************************
!    ***   GAUGE TRANSFORMATION   ***
!    ********************************

!    In a similar way to the lapse, we integrated taking maxwellF(r=0)=0,
!    but we really want maxwellF=0 at infinity.  We can fix this by making
!    a gauge transformation of the form:
!
!    maxwellF -> maxwellF - F_infty ,   omega ->  omega + q F_infty
!
!    with F_infty the asymptotic value of maxwellF.  This leaves the
!    combination (omega - q*maxwellF) unchanged. We find F_infty in
!    exactly the same way as we did it for the lapse above.
!
!    Notice that this in in fact a true gauge transformation.
!    The general form of the gauge transformation is:
!
!    a    ->  a   +  d  h
!     mu       mu     mu
!
!                  i h t
!    phi  ->  phi e
!
!              i h t
!    x  ->  X e
!
!    with "a_mu" the potential 1-form of the electromagnetic field
!    "phi" the complex scalar field, and "h" an arbitrary scalar function
!    of spacetime.  In this case we only have a_0 different from
!    zero, and we take h=kt, with k some constant. This reduces
!    to the transformation we wrote above for k=-F_infty.

     if (order=="two") then
        Ffac = maxwellF_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
             *(3.d0*maxwellF_g(0,Nrtotal) - 4.d0*maxwellF_g(0,Nrtotal-1) &
             + maxwellF_g(0,Nrtotal-2))
     else
        Ffac = maxwellF_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
             *(25.d0*maxwellF_g(0,Nrtotal) - 48.d0*maxwellF_g(0,Nrtotal-1) &
             + 36.d0*maxwellF_g(0,Nrtotal-2) - 16.d0*maxwellF_g(0,Nrtotal-3) &
             + 3.d0*maxwellF_g(0,Nrtotal-4))/3.d0
     end if

!    Gauge transformation.

     maxwellF_g = maxwellF_g - Ffac


!    ***************************************************
!    ***   RECONSTRUCT SCALAR POTENTIAL maxwellPhi   ***
!    ***************************************************

!    Reconstruct ePhi in terms of F:  ePhi = F/alpha.

     maxwellPhi_g = maxwellF_g/alpha_g


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

  if (size==1) then

     A     = A_g
     alpha = alpha_g

     dirac_FR = F_g
     dirac_GI = G_g

     ePhi     = maxwellPhi_g
     electric = maxwellE_g

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! arrays with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  else

     call distribute(0,Nl-1,A,A_g)
     call distribute(0,Nl-1,alpha,alpha_g)

     call distribute(0,Nl-1,dirac_FR,F_g)
     call distribute(0,Nl-1,dirac_GI,F_g)

     call distribute(0,Nl-1,ePhi,maxwellPhi_g)
     call distribute(0,Nl-1,electric,maxwellE_g)

  end if


! **************************************************
! ***   IMAGINARY PART OF F AND REAL PART OF G   ***
! **************************************************

! Set the imaginary part of F and real part of G to zero.

  dirac_FI = 0.d0
  dirac_GR = 0.d0

! Also set Maxwell vector potential to zero.

  eAr = 0.d0


! ***************
! ***   END   ***
! ***************

  end subroutine idata_ChargedDiracstar








! **********************************
! ***   RADIAL DERIVATIVE OF A   ***
! **********************************

! The radial derivative of A comes from the
! Hamiltonian constraint.

  function J1_CDIR(A,alpha,F,G,maxwellE,maxwellF,rm)

  use param

  implicit none

  real(8) J1_CDIR
  real(8) A,alpha,F,G,maxwellE,maxwellF,rm
  real(8) rho
  real(8) W,smallpi

! Numbers.

  smallpi = acos(-1.d0)

! W = omega - q maxwellF

  W = dirac_omega - dirac_q*maxwellF

! Energy density:
!
! rho  =  W / (2 pi alpha) (F^2 + G^2) + A maxwellE^2 / (8 pi)

  rho = W/(2.d0*smallpi*alpha)*(F**2 + G**2) + 0.125d0*A*maxwellE**2/smallpi

! dA/dr = A [ (1-A)/r + 8 pi r A rho ]

  J1_CDIR = A*((1.d0-A)/rm + 8.d0*smallpi*rm*A*rho)

  end function J1_CDIR







! **************************************
! ***   RADIAL DERIVATIVE OF alpha   ***
! **************************************

! The radial derivative of alpha comes from the
! polar-areal slicing condition.

  function J2_CDIR(A,alpha,F,G,maxwellE,maxwellF,rm)

  use param

  implicit none

  real(8) J2_CDIR
  real(8) A,alpha,F,G,maxwellE,maxwellF,rm
  real(8) rho,SA
  real(8) W,smallpi

! Numbers.

  smallpi = acos(-1.d0)

! W = omega - q maxwellF

  W = dirac_omega - dirac_q*maxwellF

! Energy density:
!
! rho  =  W / (2 pi alpha) (F^2 + G^2) + A maxwellE^2 / (8 pi)

  rho = W/(2.d0*smallpi*alpha)*(F**2 + G**2) + 0.125d0*A*maxwellE**2/smallpi

! SA = rho - 1/pi [ f g / r + m/2 (f^2 - g^2) ]

  SA = rho - (F*G/rm + 0.5d0*dirac_mass*(F**2 - G**2))/smallpi

! dalpha/dr = alpha [ (A-1)/2r + 4 pi r A SA ]

  J2_CDIR = alpha*(0.5d0*(A-1.d0)/rm + 4.d0*smallpi*rm*A*SA)

  end function J2_CDIR







! **********************************
! ***   RADIAL DERIVATIVE OF F   ***
! **********************************

! The radial derivative of F comes from the Dirac
! equation assuming a harmonic time dependence.

  function J3_CDIR(A,alpha,F,G,maxwellE,maxwellF,rm)

  use param

  implicit none

  real(8) J3_CDIR
  real(8) A,alpha,F,G,maxwellE,maxwellF,rm
  real(8) rho,SA,Dalpha
  real(8) W,smallpi

! Numbers.

  smallpi = acos(-1.d0)

! W = omega - q maxwellF

  W = dirac_omega - dirac_q*maxwellF

! Energy density:
!
! rho  =  W / (2 pi alpha) (F^2 + G^2) + A maxwellE^2 / (8 pi)

  rho = W/(2.d0*smallpi*alpha)*(F**2 + G**2) + 0.125d0*A*maxwellE**2/smallpi

! SA = rho - 1/pi [ f g / r + m/2 (f^2 - g^2) ]

  SA = rho - (F*G/rm + 0.5d0*dirac_mass*(F**2 - G**2))/smallpi

! dalpha/dr = alpha [ (A-1)/2r + 4 pi r A SA ]

  Dalpha = alpha*(0.5d0*(A-1.d0)/rm + 4.d0*smallpi*rm*A*SA)

!                                              1/2             1/2
! dF/dr  =  - F [ d alpha / 2 alpha  +  ( 1 - A   ) / r ]  -  A   G ( m + W / alpha )
!                  r

  J3_CDIR = - F*(Dalpha/(2.d0*alpha) + (1.d0 - sqrt(A))/rm) &
          - sqrt(A)*G*(dirac_mass + W/alpha)

  end function J3_CDIR







! **********************************
! ***   RADIAL DERIVATIVE OF G   ***
! **********************************

! The radial derivative of G comes from the Dirac
! equation assuming a harmonic time dependence.

  function J4_CDIR(A,alpha,F,G,maxwellE,maxwellF,rm)

  use param

  implicit none

  real(8) J4_CDIR
  real(8) A,alpha,F,G,maxwellE,maxwellF,rm
  real(8) rho,SA,Dalpha
  real(8) W,smallpi

! Numbers.

  smallpi = acos(-1.d0)

! W = omega - q maxwellF

  W = dirac_omega - dirac_q*maxwellF

! Energy density:
!
! rho  =  W / (2 pi alpha) (F^2 + G^2) + A maxwellE^2 / (8 pi)

  rho = W/(2.d0*smallpi*alpha)*(F**2 + G**2) + 0.125d0*A*maxwellE**2/smallpi

! SA = rho - 1/pi [ f g / r + m/2 (f^2 - g^2) ]

  SA = rho - (F*G/rm + 0.5d0*dirac_mass*(F**2 - G**2))/smallpi

! dalpha/dr = alpha [ (A-1)/2r + 4 pi r A SA ]

  Dalpha = alpha*(0.5d0*(A-1.d0)/rm + 4.d0*smallpi*rm*A*SA)

!                                              1/2             1/2
! dG/dr  =  - G [ d alpha / 2 alpha  +  ( 1 + A   ) / r ]  -  A   F ( m - W / alpha )
!                  r

  J4_CDIR = - G*(Dalpha/(2.d0*alpha) + (1.d0 + sqrt(A))/rm) &
          - sqrt(A)*F*(dirac_mass - W/alpha)

  end function J4_CDIR







! *****************************************
! ***   RADIAL DERIVATIVE OF maxwellE   ***
! *****************************************

! The radial derivative of maxwellE comes from the Gauss constraint
! and takes the form:
!
! dmaxwellE/dr =  - E [ 2/r + (dA/dr) / 2A ] + 4 pi echarge
!
! where in the second equality we used the Hamiltonian constraint
! to eliminate dA/dr, and with "echarge" the charge density of
! the Dirac field.

  function J5_CDIR(A,alpha,F,G,maxwellE,maxwellF,rm)

  use param

  implicit none

  real(8) J5_CDIR
  real(8) A,alpha,F,G,maxwellE,maxwellF,rm
  real(8) rho,DA,echarge
  real(8) W,smallpi

! Numbers.

  smallpi = acos(-1.d0)

! W = omega - q maxwellF

  W = dirac_omega - dirac_q*maxwellF

! Energy density:
!
! rho  =  W / (2 pi alpha) (F^2 + G^2) + A maxwellE^2 / (8 pi)

  rho = W/(2.d0*smallpi*alpha)*(F**2 + G**2) + 0.125d0*A*maxwellE**2/smallpi

! dA/dr = A [ (1-A)/r + 8 pi r A rho ]

  DA = A*((1.d0-A)/rm + 8.d0*smallpi*rm*A*rho)

! Charge density:
!                  2     2
! echarge  =  q ( F  +  G ) / 2 pi

  echarge = 0.5d0*dirac_q*(F**2 + G**2)/smallpi

! dmaxwellE/dr =  - E [ 2/r + (dA/dr) / 2A ] + 4 pi echarge

  J5_CDIR = - maxwellE*(0.5d0*DA/A + 2.d0/rm) + 4.d0*smallpi*echarge

  end function J5_CDIR







! *****************************************
! ***   RADIAL DERIVATIVE OF maxwellF   ***
! *****************************************

! The radial derivative of maxwellF comes from
! the Maxwell equations.

  function J6_CDIR(A,alpha,F,G,maxwellE,maxwellF,rm)

  use param

  implicit none

  real(8) J6_CDIR
  real(8) A,alpha,F,G,maxwellE,maxwellF,rm

! dmaxwellF/dr  =  - alpha A maxwellE

  J6_CDIR = - alpha*A*maxwellE

  end function J6_CDIR







! **********************************************
! ***   RADIAL DERIVATIVE OF A (VERSION 2)   ***
! **********************************************

! The radial derivative of A comes from the
! Hamiltonian constraint.
!
! This second version is for perturbed initial data,
! it requires previous knowledge of the energy density
! but does not assume the harmonic time dependence.

  function J7_CDIR(A,alpha,F,G,rho,rm)

  use param

  implicit none

  real(8) J7_CDIR
  real(8) A,alpha,F,G,rho,rm
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! dA/dr = A [ (1-A)/r + 8 pi r A rho ]

  J7_CDIR = A*((1.d0-A)/rm + 8.d0*smallpi*rm*A*rho)

  end function J7_CDIR








! **************************************************
! ***   RADIAL DERIVATIVE OF ALPHA (VERSION 2)   ***
! **************************************************

! The radial derivative of alpha comes from the
! polar-areal slicing condition.
!
! This second version is for perturbed initial data,
! it requires previous knowledge of the energy density
! but does not assume the harmonic time dependence.

  function J8_CDIR(A,alpha,F,G,rho,rm)

  use param

  implicit none

  real(8) J8_CDIR
  real(8) A,alpha,F,G,rho,rm
  real(8) SA
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! SA = rho - 1/pi [ f g / r + m/2 (f^2 - g^2) ]

  SA = rho - (F*G/rm + 0.5d0*dirac_mass*(F**2 - G**2))/smallpi

! dalpha/dr = alpha [ (A-1)/2r + 4 pi r A SA ]

  J8_CDIR = alpha*(0.5d0*(A-1.d0)/rm + 4.d0*smallpi*rm*A*SA)

  end function J8_CDIR







! *****************************************************
! ***   RADIAL DERIVATIVE OF maxwellE (VERSION 2)   ***
! *****************************************************

! The radial derivative of maxwellE comes from the Gauss constraint
! and takes the form:
!
! dmaxwellE/dr =  - E [ 2/r + (dA/dr) / 2A ] + 4 pi echarge
!
! where in the second equality we used the Hamiltonian constraint
! to eliminate dA/dr, and with "echarge" the charge density of
! the Dirac field.
!
! This second version is for perturbed initial data,
! it requires previous knowledge of the energy density
! but does not assume the harmonic time dependence.

  function J9_CDIR(A,alpha,F,G,maxwellE,maxwellF,rho,rm)

  use param

  implicit none

  real(8) J9_CDIR
  real(8) A,alpha,F,G,maxwellE,maxwellF,rm
  real(8) rho,DA,echarge
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! dA/dr = A [ (1-A)/r + 8 pi r A rho ]

  DA = A*((1.d0-A)/rm + 8.d0*smallpi*rm*A*rho)

! Charge density:
!                  2     2
! echarge  =  q ( F  +  G ) / 2 pi

  echarge = 0.5d0*dirac_q*(F**2 + G**2)/smallpi

! dmaxwellE/dr =  - E [ 2/r + (dA/dr) / 2A ] + 4 pi echarge

  J9_CDIR = - maxwellE*(0.5d0*DA/A + 2.d0/rm) + 4.d0*smallpi*echarge

  end function J9_CDIR

