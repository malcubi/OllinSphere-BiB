!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_ChargedBosonstar.f90,v 1.23 2023/03/01 19:20:52 malcubi Exp $

  subroutine idata_chargedboson

! *******************************************
! ***   CHARGED BOSON STAR INITIAL DATA   ***
! *******************************************

! Boson stars are solutions such that the spacetime is static
! and the complex scalar field has a harmonic dependence on time.
! This subroutine calculates initial data for a charged boson
! star using a shooting method in the polar-areal gauge (PA).

! MAIN PROBLEM: The main problem with this routine is that
! it requires a very fine tuned initial guess since we are
! fighting a growing exponential.  This gets worse as
! we move the boundary further away. So the way to
! proceed is to put the boundary close in, find omega,
! and then move it slowly further out adjusting omega
! to high precision as we move out.

! The initial data is almost the same as in the non-charged case,
! one only has to add the terms corresponding to the electric
! potential/field, and replace:
!
! omega/alpha  ->  omega/alpha - q ePhi
!
! where ePhi is the scalar (electric) potential and q is the
! electric charge parameter.  It is in fact easier to work
! with the rescaled potential F := alpha ePhi, so that the
! above substitution is:
!
! omega/alpha  ->  (omega - q F)/alpha
!
! Thus defined, F is just minus the time component of the
! potenrtial 1-form:  F = -a_0.
!
! To obtain the initial data we assume that spacetime is
! static (K_ij=0), and also that the complex scalar field
! has the form:
!
! Phi(t,r) = phi(r) exp(i omega t)
!
! Notice that with this ansatz the stress-energy tensor
! is time-independent so the metric can be static.
! The eletric field is also time-independent and real,
! so we are in electrostatics and we can take the
! vector potential equal to zero, eAr=0.
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
! Hamiltonian constraint and the Klein-Gordon equation to
! find (xi=dphi/dr):
!
!              /                      /            2              2       2      2    \          2 \
! d(A)/dr  = A | (1 - A)/r + 4 pi r A | 2 V  +  phi  (omega - q F) / alpha  +  xi / A | + r (A E)  |
!              \                      \                                               /            /
!
!                 /                                  \       /                       2       2 \
! d(xi)/dr = - xi | 2/r + d alpha / alpha - d A/(2A) |  +  A | V' - phi (omega - q F) / alpha  |
!                 \        r                 r       /       \                                 /
!
! In the above expressions V(phi) is the self-interaction
! potential of the scalar field and V'=dV/dphi. Also, E is
! the radial component of the electric field with index up.
!
! For the lapse we use the polar slicing condition
! K_{theta,theta}=0, which in this case implies:
!
!                   /                                             2 \
! dalpha/dr = alpha | (A - 1)/r + d A/(2A) - 8 pi r A V  - r (A E)  |
!                   \              r                                /
!
! Notice that for a static solution this should be equivalent
! to maximal slicing but is easier to solve.
!
! We also need to solve for the scalar potential F and the
! electric field E. The equation for F is obtained from the
! definition of E in terms of the potentials, essentially
! d(F)/dr = -E, but with metric corrections:
!
!
! d(F)/dr  = - alpha A E
!
!
! The equation for the electric field E comes from the Gauss
! constraint: D_i E^i = 4 pi rho_e, with rho_e the charge
! density, and takes the form:
!
!                 /                  \
! d(E)/dr  =  - E | 2/r  +  d A/(2A) |  +  4 pi rho_e
!                 \          r       /
!
! where the charge density is given by:
!
!              2
! rho_e = q phi  (omega -  q F) / alpha
!
!
! This routine takes as input parameters the value of the scalar
! field at the origin "boson_phi0", and two initial guesses to
! bracket the frequency: omega_left,omega_right. It then uses a
! shooting method with Runge-Kutta) to determine the value of the
! frequency "boson_omega" that is compatible with the asymptotically
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
! F(r=0)     = 0,     d F(r=0)  = 0
!                      r
!
! E(r=0)     = 0
!
! Notice that in the final solution we don't want alpha(r=0)=1,
! but rather alpha=1 at infinity.  But this is no problem as the
! slicing condition above is linear in alpha so we can always
! just rescale the lapse at the end, but in order not to affect
! the solution we must also rescale the final value of the frequency
! omega and the scalar electric potential F by the same factor.
!
! Also, we don't really want F(r=0)=0, but rather F=0 at infinity.
! We can fix this by making a gauge transformation of the field
! once we found the solution that also changes the frequency
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
  real(8) A0,alpha0                     ! Initial values of (A,alpha).
  real(8) phi0,xi0,piI0                 ! Initial values of (phi,xi).
  real(8) E0,F0                         ! Initial values of (E,F).
  real(8) DV0                           ! Initial value for derivative of potential.
  real(8) k11,k12,k13,k14               ! Runge-Kutta sources for A.
  real(8) k21,k22,k23,k24               ! Runge-Kutta sources for alpha.
  real(8) k31,k32,k33,k34               ! Runge-Kutta sources for phi.
  real(8) k41,k42,k43,k44               ! Runge-Kutta sources for xi.
  real(8) k51,k52,k53,k54               ! Runge-Kutta sources for F.
  real(8) k61,k62,k63,k64               ! Runge-Kutta sources for E.
  real(8) A_rk,alpha_rk                 ! Runge-Kutta values of (A,alpha).
  real(8) phi_rk,xi_rk,piI_rk           ! Runge-Kutta values (phi,xi).
  real(8) F_rk,E_rk                     ! Runge-Kutta values of (F,E).
  real(8) V_rk,K_rk                     ! Runge-Kutta values of (V,K) (potential and kinetic term).
  real(8) J1_CB,J2_CB,J3_CB             ! Functions for sources of differential equations.
  real(8) J4_CB,J5_CB,J6_CB             ! Functions for sources of differential equations.
  real(8) res,res_old                   ! Residual.
  real(8) omega_new,omega_old,domega    ! Trial frequency and frequency interval. 
  real(8) half,smallpi                  ! Numbers.
  real(8) rm,alphafac,Ffac,aux          ! Auxiliary quantities.
  real(8) :: epsilon = 1.d-8            ! Tolerance.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                 ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g,alpha_g        ! Radial metric and lapse global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: phi_g,xi_g,piI_g   ! Scalar field global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: E_g,F_g,ePhi_g     ! Electric field and scalar potential global arrays.


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  smallpi = acos(-1.d0)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for a charged boson star in polar-areal gauge using shooting method'
     print *
  end if

! Sanity check.

  if (spacetime=="minkowski") then
     print *, 'Charged boson star initial data is not compatible with a Minkowski background ...'
     print *, 'Aborting! (subroutine idata_ChargedBosonstar)'
     print *
     call die
  end if


! **************************************
! ***   SCALAR FIELD NORMALIZATION   ***
! **************************************

! Physical normalization.

  if (boson_factor=="physical") then

     if (rank==0) then
        write(*,'(A)') ' Using physical normalization at origin: phi(r=0) = phi0'
     end if

! Alternative normalization.

  else if (boson_factor=="harmonic") then

     boson_phi0 = boson_phi0/sqrt(4.d0*smallpi)

     if (rank==0) then
        write(*,'(A,E16.10)') ' Using harmonic normalization at origin: phi(r=0) = phi0/sqrt(4pi) = ', boson_phi0
     end if

  end if

! Rescale the tolerance with boson_phi0, so the final tolerance
! is measured with respect to this value.

  epsilon = epsilon*boson_phi0


! ********************************
! ***   CHARGE NORMALIZATION   ***
! ********************************

! Standard normalization.

  if (charge_factor=="standard") then

     if (rank==0) then
        write(*,'(A)') ' Using standard charge normalization'
        print *
     end if

! The normalization used by Jetzer has a factor of sqrt(2).
! (ver: Ph. Jetzer (Phys. Lett. B, vol. 227, p. 341 (1989).)

  else if (charge_factor=="jetzer") then

     complex_q = complex_q*sqrt(2.d0)

     if (rank==0) then
        write(*,'(A)') ' Using Jetzer charge normalization'
        print *
     end if

  end if


! ************************
! ***   SANITY CHECK   ***
! ************************

! Sanity check for potential.

  if ((complexpotential=="phi2").and.(complex_lambda/=0.d0)) then
     print *, 'You can not have a complexpotential=phi2 and complex_lambda different from zero.'
     print *, 'Aborting! (subroutine idata_ChargedBosonstar)'
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

  phi_g = 0.d0
  xi_g  = 0.d0

  F_g = 0.d0
  E_g = 0.d0


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
        print *, 'For charged boson star initial data we must have omega_right>omega_left.'
        print *, 'Aborting! (subroutine idata_ChargedBosonstar)'
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

     boson_omega = omega_left

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

!          Derivative of potential at the origin.

           aux = boson_phi0
           DV0 = complex_mass**2*aux + complex_lambda*aux**3

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

              phi_g(l,imin-1)   = (9.d0*(phi_g(l+1,Nrtotal-2)+phi_g(l+1,Nrtotal-3)) &
                                - (phi_g(l+1,Nrtotal-4)+phi_g(l+1,Nrtotal-1)))/16.d0
              xi_g(l,imin-1)    = (9.d0*(xi_g(l+1,Nrtotal-2)+xi_g(l+1,Nrtotal-3)) &
                                - (xi_g(l+1,Nrtotal-4)+xi_g(l+1,Nrtotal-1)))/16.d0

              F_g(l,imin-1)     = (9.d0*(F_g(l+1,Nrtotal-2)+F_g(l+1,Nrtotal-3)) &
                                - (F_g(l+1,Nrtotal-4)+F_g(l+1,Nrtotal-1)))/16.d0
              E_g(l,imin-1)     = (9.d0*(E_g(l+1,Nrtotal-2)+E_g(l+1,Nrtotal-3)) &
                                - (E_g(l+1,Nrtotal-4)+E_g(l+1,Nrtotal-1)))/16.d0

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

!                Values of (phi,xi) at origin.

                 phi0 = boson_phi0
                 xi0  = 0.d0

!                Values of (F,E) at origin.

                 F0 = 0.d0
                 E0 = 0.d0

!             Grid spacing and values at previous grid point.

              else

                 delta = dr(l)
                 r0    = rr(l,i-1)

                 A0     = A_g(l,i-1)
                 alpha0 = alpha_g(l,i-1)

                 phi0 = phi_g(l,i-1)
                 xi0  = xi_g(l,i-1)

                 F0 = F_g(l,i-1)
                 E0 = E_g(l,i-1)
 
              end if

!             I) First Runge-Kutta step.

!             Sources at first grid point if we start
!             from the origin (for finer grid).

              if (i==1) then

!                At the origin we have:  A' = alpha' = phi' = F' = 0.

                 k11 = 0.d0
                 k21 = 0.d0
                 k31 = 0.d0
                 k51 = 0.d0

!                Source for xi at origin. We assume that for small r
!                we have xi = k r and substitute in the equation for
!                dxi/dr at the origin. From that we solve for k to
!                find:
!
!                k  =  (V' - omega^2 phi0)/3

                 k41 = delta*(DV0 - boson_omega**2*boson_phi0)/3.d0

!                Source for E at origin.  We assume that for small r
!                we have E = k r and substitute in the equation for
!                dE/dr at the origin.  From that we solve for k to
!                find:
!
!                k = (4pi/3) q omega phi^2

                 k61 = delta*(4.d0*smallpi/3.d0*complex_q*boson_omega*boson_phi0**2)

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk     = A0
                 alpha_rk = alpha0

                 phi_rk = phi0
                 xi_rk  = xi0

                 F_rk   = F0
                 E_rk   = E0

!                Physical potential V and kinetic term K.

                 V_rk = 0.5d0*complex_mass**2*phi_rk**2 + 0.25d0*complex_lambda*phi_rk**4
                 K_rk = 0.5d0*(xi_rk**2/A_rk + phi_rk**2*((boson_omega - complex_q*F_rk)/alpha_rk)**2)

!                Sources.

                 k11 = delta*J1_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
                 k21 = delta*J2_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

                 k31 = delta*J3_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
                 k41 = delta*J4_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

                 k51 = delta*J5_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
                 k61 = delta*J6_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk     = A0     + half*k11
              alpha_rk = alpha0 + half*k21

              phi_rk   = phi0   + half*k31
              xi_rk    = xi0    + half*k41

              F_rk     = F0     + half*k51
              E_rk     = E0     + half*k61

!             Physical potential V and kinetic term K.

              V_rk = 0.5d0*complex_mass**2*phi_rk**2 + 0.25d0*complex_lambda*phi_rk**4
              K_rk = 0.5d0*(xi_rk**2/A_rk + phi_rk**2*((boson_omega - complex_q*F_rk)/alpha_rk)**2)

!             Sources.

              k12 = delta*J1_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k22 = delta*J2_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              k32 = delta*J3_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k42 = delta*J4_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              k52 = delta*J5_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k62 = delta*J6_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk     = A0     + half*k12
              alpha_rk = alpha0 + half*k22

              phi_rk   = phi0   + half*k32
              xi_rk    = xi0    + half*k42

              F_rk     = F0     + half*k52
              E_rk     = E0     + half*k62

!             Physical potential V and kinetic term K.

              V_rk = 0.5d0*complex_mass**2*phi_rk**2 + 0.25d0*complex_lambda*phi_rk**4
              K_rk = 0.5d0*(xi_rk**2/A_rk + phi_rk**2*((boson_omega - complex_q*F_rk)/alpha_rk)**2)

!             Sources.

              k13 = delta*J1_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k23 = delta*J2_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              k33 = delta*J3_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k43 = delta*J4_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              k53 = delta*J5_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k63 = delta*J6_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk     = A0     + k13
              alpha_rk = alpha0 + k23

              phi_rk   = phi0   + k33
              xi_rk    = xi0    + k43

              F_rk     = F0     + k53
              E_rk     = E0     + k63

!             Physical potential V and kinetic term K.

              V_rk = 0.5d0*complex_mass**2*phi_rk**2 + 0.25d0*complex_lambda*phi_rk**4
              K_rk = 0.5d0*(xi_rk**2/A_rk + phi_rk**2*((boson_omega - complex_q*F_rk)/alpha_rk)**2)

!             Sources.

              k14 = delta*J1_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k24 = delta*J2_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              k34 = delta*J3_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k44 = delta*J4_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              k54 = delta*J5_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k64 = delta*J6_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

!             Advance variables to next grid point.

              A_g(l,i)     = A0     + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i) = alpha0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

              phi_g(l,i) = phi0 + (k31 + 2.d0*(k32 + k33) + k34)/6.d0
              xi_g(l,i)  = xi0  + (k41 + 2.d0*(k42 + k43) + k44)/6.d0

              F_g(l,i) = F0 + (k51 + 2.d0*(k52 + k53) + k54)/6.d0
              E_g(l,i) = E0 + (k61 + 2.d0*(k62 + k63) + k64)/6.d0

!             Check if solution is blowing up. This helps to reduce
!             the need for a very fine tuned initial guess.

              if (abs(phi_g(l,i))>2.d0*boson_phi0) then
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

              if (abs(phi_g(l,i))+abs(xi_g(l,i-1))<epsilon) then
                 phi_g(l,i) = 0.d0
                 xi_g(l,i)  = 0.d0
              end if

           end do


!          ***********************
!          ***   GHOST ZONES   ***
!          ***********************

           do i=1,ghost

              A_g(l,1-i)     = + A_g(l,i)
              alpha_g(l,1-i) = + alpha_g(l,i)

              phi_g(l,1-i) = + phi_g(l,i)
              xi_g(l,1-i)  = - xi_g(l,i)

              F_g(l,1-i) = + F_g(l,i)
              E_g(l,1-i) = - E_g(l,i)

           end do

        end do


!       *******************************************
!       ***   FIND RESIDUAL AND CORRECT OMEGA   ***
!       *******************************************

!       Find the difference between the solution at the coarsest level
!       and the expected asymptotic behaviour. The solution should
!       decay exponetially far away as:
!
!       phi ~ exp(-k*r)
!
!       so that:
!
!       phi' + k phi = 0
!
!       Notice that far away we must have:  k = m**2 - (omega/alpha)**2

        res_old = res

        if (abs(phi_g(0,Nrtotal))+abs(phi_g(0,Nrtotal-1))<epsilon) then
           res = epsilon/2.d0
           goto 100
        else
           aux = complex_mass**2 - (boson_omega/alpha(0,Nrtotal))**2
           res = xi_g(0,Nrtotal) - log(abs(phi_g(0,Nrtotal)))/rr(0,Nrtotal)*phi_g(0,Nrtotal)
        end if

!       Secant method:  Having found the difference for the two values of omega
!       we can use it to extrapolate linearly to the next best guess.

        if (.not.left) then

           left = .true.

           omega_old = boson_omega
           boson_omega = omega_right

        else

           if (.not.right) right = .true.

           aux = res - res_old

           if ((aux/=0.d0).and.(abs(aux)<1.d-1)) then
              omega_new = boson_omega - res*(boson_omega - omega_old)/aux
           else
              omega_new = 0.5d0*(omega_old + boson_omega)
           end if

           omega_old = boson_omega
           boson_omega = omega_new

        end if

100     continue


!       **************************
!       ***   END ITERATIONS   ***
!       **************************

!       Output data to screen.

        write(*,"(A,I4,A,ES22.16,A,ES9.2)") ' Iteration: ',iter,'    Frequency: ',boson_omega,'    Residual: ',res

     end do

!    Message in case we reached the maximum number of iterations.

     if (iter>=maxiter) then
        print *
        print *, 'Maximum iteration number reached, initial data solver did not converge.'
        print *, 'Aborting! (subroutine idata_ChargedBosonstar)'
        print *
     else
        print *
        print *, 'Done!'
        print *
     end if

!    Output value of omega to screen.

     if (rank==0) then
        write(*,'(A,E22.16)') ' Omega (eigenvalue)        = ', boson_omega
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

           phi_g(l-1,iaux) = (9.d0*(phi_g(l,i)+phi_g(l,i+1)) - (phi_g(l,i-1)+phi_g(l,i+2)))/16.d0
           xi_g(l-1,iaux)  = (9.d0*(xi_g(l,i)+xi_g(l,i+1))   - (xi_g(l,i-1)+xi_g(l,i+2)))/16.d0

           F_g(l-1,iaux) = (9.d0*(F_g(l,i)+F_g(l,i+1)) - (F_g(l,i-1)+F_g(l,i+2)))/16.d0
           E_g(l-1,iaux) = (9.d0*(E_g(l,i)+E_g(l,i+1)) - (E_g(l,i-1)+E_g(l,i+2)))/16.d0

        end do

!       Fix ghost zones.

        do i=1,ghost

           A_g(l-1,1-i)     = + A_g(l-1,i)
           alpha_g(l-1,1-i) = + alpha_g(l-1,i)

           phi_g(l-1,1-i) = + phi_g(l-1,i)
           xi_g(l-1,1-i)  = - xi_g(l-1,i)

           F_g(l-1,1-i) = + F_g(l-1,i)
           E_g(l-1,1-i) = - E_g(l-1,i)

        end do

     end do


!    ******************************
!    ***   PUGLIESE FREQUENCY   ***
!    ******************************

!    This is here only in order to compare our frequency with the
!    one defined by Pugliese et al. (Phys. Rev. D 88, 024053 (2013)).
!    It must be found BEFORE we do any of our rescalings and gauge
!    transformations.
!
!    Their frequency is defined as:
!
!    w_p = 1 - q*F_infty/boson_omega
!
!    This value comes from first defining:
!
!    C(r) := (w - q*F(r))/w
!
!    with w the original frequency. Notice that we have C(r=0)=1.
!    The Pugliese frequency is then defined as the asymptotic
!    value of C(r).
!
!    However, in the paper by Pugliese et al. they seem to evaluate
!    this number at r~10 and not at infinity!!!

!    Figure out asymptotic value of F.

     if (.false.) then

        if (order=="two") then
           Ffac = F_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
                *(3.d0*F_g(0,Nrtotal) - 4.d0*F_g(0,Nrtotal-1) + F_g(0,Nrtotal-2))
        else
           Ffac = F_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
                *(25.d0*F_g(0,Nrtotal) - 48.d0*F_g(0,Nrtotal-1) &
                + 36.d0*F_g(0,Nrtotal-2) - 16.d0*F_g(0,Nrtotal-3) + 3.d0*F_g(0,Nrtotal-4))/3.d0
        end if

!       Output Pugliese frequency evaluated at the first grid point
!       on the base grid with r>10. This could be done better, but
!       it really is here only for comparisons.

        if (rank==0) then
           if (rr(0,Nrtotal)>10.d0) then
              i=1
              do while (rr(0,i)<10.d0)
                 i=i+1
              end do
              write(*,'(A,E22.16)') ' Omega (Pugliese r=10)     = ', 1.d0 - (complex_q/boson_omega)*F_g(0,i)
           end if
        end if

!       Output Pugliese frequency extrapolated to infinity.

        if (rank==0) then
           write(*,'(A,E22.16)') ' Omega (Pugliese r=infty)  = ', 1.d0 - (complex_q/boson_omega)*Ffac
           write(*,*)
        end if

     end if


!    *************************
!    ***   RESCALE LAPSE   ***
!    *************************

!    Since we integrated out assuming alpha(r=0)=1, we now
!    need to rescale the lapse to make sure that we have
!    alpha=1 at infinity. This can be done easily since
!    the equation for the lapse is linear. But notice that
!    we also need to rescale the final frequency omega,
!    and the electric potential F, by the same factor,
!    so that the terms of the form:
!
!    (omega - qF)/alpha
!
!    do not change.
!
!    In practice we just divide alpha, omega and F by the
!    asymptotic value of the lapse.  In order to obtain
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

!    Rescale omega and output its value to screen.

     omega_new = boson_omega/alphafac

     if (rank==0) then
        write(*,'(A,E22.16)') ' Omega (lapse rescaled)    = ', omega_new
     end if

!    Rescale electric potential F.

     F_g = F_g/alphafac


!    ********************************
!    ***   IMAGINARY PART OF pi   ***
!    ********************************

!    From the original ansatz, we take the time derivative of
!    the imaginary part equal to (omega/alpha)*phi.

     piI_g = (omega_new/alpha_g)*phi_g


!    *********************************************
!    ***  OUTPUT GAUGE TRANSFORMED FREQUENCY   ***
!    *********************************************

!    Here we output the gauge transformed frequency.
!    But we don't really do the gauge transformation
!    since that should be done AFTER the perturbation
!    if there is one.

!    Figure out asymtotic value of F.

     if (order=="two") then
        Ffac = F_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
             *(3.d0*F_g(0,Nrtotal) - 4.d0*F_g(0,Nrtotal-1) + F_g(0,Nrtotal-2))
     else
        Ffac = F_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
             *(25.d0*F_g(0,Nrtotal) - 48.d0*F_g(0,Nrtotal-1) &
             + 36.d0*F_g(0,Nrtotal-2) - 16.d0*F_g(0,Nrtotal-3) + 3.d0*F_g(0,Nrtotal-4))/3.d0
     end if

!    Output transformed frequency.

     omega_new = omega_new - complex_q*Ffac

     if (rank==0) then
        write(*,'(A,E22.16)') ' Omega (gauge transformed) = ', omega_new
        print *
     end if


!    ***********************************
!    ***   BOSON STAR PERTURBATION   ***
!    ***********************************

!    Here we add a gaussian perturbation to the solution
!    for the scalar field and solve again the hamiltonian
!    constraint for A, and the Gauss constraint for the
!    electric field E and scalar potential F. We also
!    solve again the polar slicing condition for the lapse.
!
!    The perturbation should be small, and since it must
!    be even we take the sum of two symmetric gaussians.
!
!    Notice that the perturbation must be done BEFORE
!    the final gauge transformation, otherwise it causes
!    a mess with the solution since F does not have the
!    correct behaviour.

     if ((bosongauss).and.((boson_phiR_a0/=0.d0).or.(boson_piI_a0/=0.d0))) then

        if (rank==0) then
           print *, 'Adding perturbation to initial data, recalculating (alpha,A,F,E) ...'
           print *
        end if

!       Rescale the amplitude of the gaussian with max(phi).

        aux = maxval(abs(phi_g))

        boson_phiR_a0 = aux*boson_phiR_a0
        boson_piI_a0  = aux*boson_piI_a0

!       Add perturbation to phi and xi.

        if (boson_phiR_r0==0.d0) then
           phi_g = phi_g + boson_phiR_a0*exp(-rr**2/boson_phiR_s0**2)
           xi_g  = xi_g - 2.d0*rr*boson_phiR_a0/boson_phiR_s0**2*exp(-rr**2/boson_phiR_s0**2)
        else
           phi_g = phi_g + boson_phiR_a0 &
                 *(exp(-(rr-boson_phiR_r0)**2/boson_phiR_s0**2) &
                 + exp(-(rr+boson_phiR_r0)**2/boson_phiR_s0**2))
           xi_g  = xi_g - 2.d0*boson_phiR_a0/boson_phiR_s0**2 &
                 *((rr-boson_phiR_r0)*exp(-(rr-boson_phiR_r0)**2/boson_phiR_s0**2) &
                 + (rr+boson_phiR_r0)*exp(-(rr+boson_phiR_r0)**2/boson_phiR_s0**2))
        end if

!       Add perturbation to piI. This perturbation is also
!       rescaled with (omega/alpha) to make it easier
!       to compare the perturbations in phiR and piI.

        if (boson_piI_r0==0.d0) then
           piI_g = piI_g + boson_piI_a0*(omega_new/alpha_g)*exp(-rr**2/boson_piI_s0**2)
        else
           piI_g = piI_g + boson_piI_a0*(omega_new/alpha_g) &
                 *(exp(-(rr-boson_piI_r0)**2/boson_piI_s0**2) &
                 + exp(-(rr+boson_piI_r0)**2/boson_piI_s0**2))
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

              F_g(l,imin-1) = (9.d0*(F_g(l+1,Nrtotal-2)+F_g(l+1,Nrtotal-3)) &
                            - (F_g(l+1,Nrtotal-4)+F_g(l+1,Nrtotal-1)))/16.d0
              E_g(l,imin-1) = (9.d0*(E_g(l+1,Nrtotal-2)+E_g(l+1,Nrtotal-3)) &
                            - (E_g(l+1,Nrtotal-4)+E_g(l+1,Nrtotal-1)))/16.d0

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

!                Values of (phi,xi,pi) at origin.

                 phi0 = (9.d0*(phi_g(l,0)+phi_g(l,1)) - (phi_g(l,-1)+phi_g(l,2)))/16.d0
                 xi0  = 0.d0
                 piI0 = (9.d0*(piI_g(l,0)+piI_g(l,1)) - (piI_g(l,-1)+piI_g(l,2)))/16.d0

!                Values of (F,E) at origin.

                 F0 = 0.d0
                 E0 = 0.d0

!             Grid spacing and values at previous grid point.

              else

                 delta = dr(l)
                 r0    = rr(l,i-1)

                 A0     = A_g(l,i-1)
                 alpha0 = alpha_g(l,i-1)

                 phi0 = phi_g(l,i-1)
                 xi0  = xi_g(l,i-1)
                 piI0 = piI_g(l,i-1)

                 F0 = F_g(l,i-1)
                 E0 = E_g(l,i-1)

              end if

!             I) First Runge-Kutta step.

              if (i==1) then

!                At the origin we have:  A' = alpha' = F' = 0.

                 k11 = 0.d0
                 k21 = 0.d0
                 k51 = 0.d0

!                Source for E at origin.  We assume that for small r
!                we have E = k r and substitute in the equation for
!                dE/dr at the origin.  From that we solve for k to
!                find:
!
!                k = (4pi/3) q phiR piI

                 k61 = delta*(4.d0*smallpi/3.d0*complex_q*phi0*piI0)

              else

                 rm = r0

                 A_rk     = A0
                 alpha_rk = alpha0

                 phi_rk = phi0
                 xi_rk  = xi0
                 piI_rk = piI0

                 F_rk = F0
                 E_rk = E0

!                Physical potential V and kinetic term K.

                 V_rk = 0.5d0*complex_mass**2*phi_rk**2 + 0.25d0*complex_lambda*phi_rk**4
                 K_rk = 0.5d0*(xi_rk**2/A_rk + (piI_rk - complex_q*F_rk*phi_rk/alpha_rk)**2)

!                Sources.

                 k11 = delta*J1_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
                 k21 = delta*J2_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

                 k51 = delta*J5_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
                 k61 = delta*J6_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              end if

!             II) Second Runge-Kutta step. We interpolate the values of (phi,xi).

              rm = r0 + half*delta

              A_rk     = A0     + half*k11
              alpha_rk = alpha0 + half*k21

              F_rk     = F0     + half*k51
              E_rk     = E0     + half*k61

              if (i==1) then  ! Linear interpolation for first point.
                 phi_rk = 0.5d0*(phi0 + phi_g(l,1))
                 xi_rk  = 0.5d0*(xi0 + xi_g(l,1))
                 piI_rk = 0.5d0*(piI0 + piI_g(l,1))
              else            ! Cubic interpolation for the rest.
                 phi_rk = (9.d0*(phi_g(l,i)+phi_g(l,i-1)) - (phi_g(l,i-2)+phi_g(l,i+1)))/16.d0
                 xi_rk  = (9.d0*(xi_g(l,i)+xi_g(l,i-1))   - (xi_g(l,i-2)+xi_g(l,i+1)))/16.d0
                 piI_rk = (9.d0*(piI_g(l,i)+piI_g(l,i-1)) - (piI_g(l,i-2)+piI_g(l,i+1)))/16.d0
              end if

!             Physical potential V and kinetic term K.

              V_rk  = 0.5d0*complex_mass**2*phi_rk**2 + 0.25d0*complex_lambda*phi_rk**4
              K_rk  = 0.5d0*(xi_rk**2/A_rk + (piI_rk - complex_q*F_rk*phi_rk/alpha_rk)**2)

!             Sources.

              k12 = delta*J1_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k22 = delta*J2_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              k52 = delta*J5_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k62 = delta*J6_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk     = A0     + half*k12
              alpha_rk = alpha0 + half*k22

              F_rk     = F0     + half*k52
              E_rk     = E0     + half*k62

!             Physical potential V and kinetic term K.

              V_rk = 0.5d0*complex_mass**2*phi_rk**2 + 0.25d0*complex_lambda*phi_rk**4
              K_rk = 0.5d0*(xi_rk**2/A_rk + (piI_rk - complex_q*F_rk*phi_rk/alpha_rk)**2)

!             Sources.

              k13 = delta*J1_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k23 = delta*J2_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              k53 = delta*J5_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k63 = delta*J6_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta
 
              A_rk     = A0     + k13
              alpha_rk = alpha0 + k23

              F_rk     = F0     + k53
              E_rk     = E0     + k63

              phi_rk = phi_g(l,i)
              xi_rk  = xi_g(l,i)
              piI_rk = piI_g(l,i)

!             Physical potential V and kinetic term K.

              V_rk = 0.5d0*complex_mass**2*phi_rk**2 + 0.25d0*complex_lambda*phi_rk**4
              K_rk = 0.5d0*(xi_rk**2/A_rk + (piI_rk - complex_q*F_rk*phi_rk/alpha_rk)**2)

!             Sources.

              k14 = delta*J1_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k24 = delta*J2_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

              k54 = delta*J5_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)
              k64 = delta*J6_CB(A_rk,alpha_rk,phi_rk,xi_rk,F_rk,E_rk,V_rk,K_rk,rm)

!             Advance variables to next grid point.

              A_g(l,i)     = A0     + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i) = alpha0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

              F_g(l,i) = F0 + (k51 + 2.d0*(k52 + k53) + k54)/6.d0
              E_g(l,i) = E0 + (k61 + 2.d0*(k62 + k63) + k64)/6.d0

           end do

!          Ghost zones.

           do i=1,ghost

              A_g(l,1-i)     = A_g(l,i)
              alpha_g(l,1-i) = alpha_g(l,i)

              F_g(l,1-i) = + F_g(l,i)
              E_g(l,1-i) = - E_g(l,i)

           end do

        end do

!       Restrict solution from fine to coarse grid.

        do l=Nl-1,1,-1

           do i=1,Nrtotal-ghost,2

              iaux = i/2 + 1
              rm = rr(l-1,iaux)

              A_g(l-1,iaux)     = (9.d0*(A_g(l,i)+A_g(l,i+1)) - (A_g(l,i-1)+A_g(l,i+2)))/16.d0
              alpha_g(l-1,iaux) = (9.d0*(alpha_g(l,i)+alpha_g(l,i+1)) - (alpha_g(l,i-1)+alpha_g(l,i+2)))/16.d0

              F_g(l-1,iaux) = (9.d0*(F_g(l,i)+F_g(l,i+1)) - (F_g(l,i-1)+F_g(l,i+2)))/16.d0
              E_g(l-1,iaux) = (9.d0*(E_g(l,i)+E_g(l,i+1)) - (E_g(l,i-1)+E_g(l,i+2)))/16.d0

           end do

!          Fix ghost zones.

           do i=1,ghost

              A_g(l-1,1-i)     = + A_g(l-1,i)
              alpha_g(l-1,1-i) = + alpha_g(l-1,i)

              F_g(l-1,1-i) = + F_g(l-1,i)
              E_g(l-1,1-i) = - E_g(l-1,i)

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

        alpha_g = alpha_g/alphafac
        F_g = F_g/alphafac

     end if


!    ********************************
!    ***   GAUGE TRANSFORMATION   ***
!    ********************************

!    In a similar way to the lapse, we integrated taking F(r=0)=0,
!    but we really want F=0 at infinity.  We can fix this by making
!    a gauge transformation of the form:
!
!    F -> F - F_infty ,    omega ->  omega - q F_infty
!
!    with F_infty the asymptotic value of F.  This leaves the
!    combination (omega - qF) unchanged. We find F_infty in
!    exactly the same way as we did it for the lapse above.
!
!    Notice that this in in fact a true gauge transformation.
!    The general form of the gauge transformation is:
!
!    a    ->  a   +  d  h
!     mu       mu     mu
!
!                 i h t
!    phi  -> phi e
!
!    with "a_mu" the potential 1-form of the electromagnetic field
!    "phi" the complex scalar field, and "h" an arbitrary scalar function
!    of spacetime.  In this case we only have a_0 different from
!    zero, and we take h=kt, with k some constant. This reduces
!    to the transformation we wrote above for k=-F_infty.

     if (order=="two") then
        Ffac = F_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
             *(3.d0*F_g(0,Nrtotal) - 4.d0*F_g(0,Nrtotal-1) + F_g(0,Nrtotal-2))
     else
        Ffac = F_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
             *(25.d0*F_g(0,Nrtotal) - 48.d0*F_g(0,Nrtotal-1) &
             + 36.d0*F_g(0,Nrtotal-2) - 16.d0*F_g(0,Nrtotal-3) + 3.d0*F_g(0,Nrtotal-4))/3.d0
     end if

!    Gauge transformation.

     F_g = F_g - Ffac

!    After the gauge transformation we need to correct pi.
!    Notice that since the gauge tranformation for piI
!    depends on phi, when we have a perturbation the
!    transformed piI will be modified even if we set
!    boson_piI_a0=0.

     piI_g = piI_g - complex_q*Ffac*phi_g/alpha_g


!    *********************************************
!    ***   RECONSTRUCT SCALAR POTENTIAL ePhi   ***
!    *********************************************

!    Reconstruct ePhi in terms of F:  ePhi = F/alpha.

     ePhi_g = F_g/alpha_g


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

     complex_phiR = phi_g
     complex_xiR  = xi_g
     complex_piI  = piI_g

     ePhi     = ePhi_g
     electric = E_g

  else

     call distribute(0,Nl-1,A,A_g)
     call distribute(0,Nl-1,alpha,alpha_g)

     call distribute(0,Nl-1,complex_phiR,phi_g)
     call distribute(0,Nl-1,complex_xiR,xi_g)
     call distribute(0,Nl-1,complex_piI,piI_g)

     call distribute(0,Nl-1,ePhi,ePhi_g)
     call distribute(0,Nl-1,electric,E_g)

  end if


! **********************************************************
! ***   IMAGINARY PART OF (phi,xi) AND REAL PART of pi   ***
! **********************************************************

! Set the imaginary part and its spatial derivative to zero.

  complex_phiI = 0.d0
  complex_xiI  = 0.d0

! Set time derivative of real part to zero.

  complex_piR  = 0.d0

! And make sure eAr is zero.

  eAr = 0.d0


! ***************
! ***   END   ***
! ***************

  end subroutine idata_chargedboson









! **********************************
! ***   RADIAL DERIVATIVE OF A   ***
! **********************************

! The radial derivative of A comes from the
! Hamiltonian constraint.

  function J1_CB(A,alpha,phi,xi,F,E,V,K,rm)

  use param

  implicit none

  real(8) J1_CB
  real(8) A,alpha,phi,xi,F,E
  real(8) V,K,rm,rho
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Energy density for scalar field.

  rho = K + V
  
!              /                      /            2              2       2      2    \          2 \
! d(A)/dr  = A | (1 - A)/r + 4 pi r A | 2 V  +  phi  (omega - q F) / alpha  +  xi / A | + r (A E)  |
!              \                      \                                               /            /

  J1_CB = A*((1.d0 - A)/rm + 8.d0*smallpi*rm*A*rho + rm*(A*E)**2)

  end function J1_CB







! **************************************
! ***   RADIAL DERIVATIVE OF ALPHA   ***
! **************************************

! The radial derivative of alpha comes from 
! the polar slicing condition:  dK_{theta,theta}/dt=0.

  function J2_CB(A,alpha,phi,xi,F,E,V,K,rm)

  use param

  implicit none

  real(8) J2_CB
  real(8) A,alpha,phi,xi,F,E
  real(8) V,K,rm,rho
  real(8) DA,smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Energy density for scalar field.

  rho = K + V

! Radial derivative of A (identical to function J1_CB).

  DA = A*((1.d0 - A)/rm + 8.d0*smallpi*rm*A*rho + rm*(A*E)**2)

!                   /                                             2 \
! dalpha/dr = alpha | (A - 1)/r + d A/(2A) - 8 pi r A V  - r (A E)  |
!                   \              r                                /

  J2_CB = alpha*((A - 1.d0)/rm + 0.5d0*DA/A - 8.d0*smallpi*rm*A*V - rm*(A*E)**2)

  end function J2_CB







! ************************************
! ***   RADIAL DERIVATIVE OF PHI   ***
! ************************************

! This is just the definition of xi:  xi=dphi/dr

  function J3_CB(A,alpha,phi,xi,F,E,V,K,rm)

  use param

  implicit none

  real(8) J3_CB
  real(8) A,alpha,phi,xi,F,E
  real(8) V,K,rm

! dphi/dr = xi.

  J3_CB = xi

  end function J3_CB







! ***********************************
! ***   RADIAL DERIVATIVE OF XI   ***
! ***********************************

! The radial derivative of xi comes from the
! Klein-Gordon equation.

  function J4_CB(A,alpha,phi,xi,F,E,V,K,rm)

  use param

  implicit none

  real(8) J4_CB
  real(8) A,alpha,phi,xi,F,E
  real(8) V,K,rm
  real(8) DV,rho
  real(8) DA,Dalpha,smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Energy density for scalar field.

  rho = K + V

! Potential derivative.

  DV = complex_mass**2*phi + complex_lambda*phi**3

! Radial derivative of A (identical to function J1_CB).

  DA = A*((1.d0 - A)/rm + 8.d0*smallpi*rm*A*rho + rm*(A*E)**2)

! Radial derivative of alpha (identical to function J2_CB).

  Dalpha = alpha*((A-1.d0)/rm + 0.5d0*DA/A - 8.d0*smallpi*rm*A*V - rm*(A*E)**2)

!                 /                                  \       /                       2       2 \
! d(xi)/dr = - xi | 2/r + d alpha / alpha - d A/(2A) |  +  A | V' - phi (omega - q F) / alpha  |
!                 \        r                 r       /       \                                 /

  J4_CB = - xi*(2.d0/rm + Dalpha/alpha - 0.5d0*DA/A) &
        + A*(DV - phi*(boson_omega - complex_q*F)**2/alpha**2)

  end function J4_CB







! **********************************
! ***   RADIAL DERIVATIVE OF F   ***
! **********************************

! The radial derivative of F comes from the
! definition of the electric field in terms
! of the scalar potential

  function J5_CB(A,alpha,phi,xi,F,E,V,K,rm)

  use param

  implicit none

  real(8) J5_CB
  real(8) A,alpha,phi,xi,F,E
  real(8) V,K,rm

! d(F)/dr  = - alpha A E

  J5_CB = - alpha*A*E

  end function J5_CB







! **********************************
! ***   RADIAL DERIVATIVE OF E   ***
! **********************************

! The radial derivative of E comes from
! the Gauss constraint.

  function J6_CB(A,alpha,phi,xi,F,E,V,K,rm)

  use param

  implicit none

  real(8) J6_CB
  real(8) A,alpha,phi,xi,F,E
  real(8) V,K,rm,rho,rhoe
  real(8) DA,smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Energy density for scalar field.

  rho = K + V

! Radial derivative of A (identical to function J1_CB).

  DA = A*((1.d0 - A)/rm + 8.d0*smallpi*rm*A*rho + rm*(A*E)**2)

!                 /                  \
! d(E)/dr  =  - E | 2/r  +  d A/(2A) |  +  4 pi rho_e
!                 \          r       /
!
! where the charge density is given by:
!
!              2
! rho_e = q phi  (omega -  q F) / alpha

  J6_CB = - E*(2.d0/rm + 0.5d0*DA/A) &
        + 4.d0*smallpi*complex_q*phi**2*(boson_omega - complex_q*F)/alpha

  end function J6_CB

