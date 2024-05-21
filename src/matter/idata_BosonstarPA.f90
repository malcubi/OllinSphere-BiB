!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_BosonstarPA.f90,v 1.119 2024/05/21 19:10:42 malcubi Exp $

  subroutine idata_BosonstarPA

! *******************************************************
! ***   BOSON STAR INITIAL DATA IN POLAR-AREAL GAUGE  ***
! *******************************************************

! Boson stars are solutions such that the spacetime is static
! and the complex scalar field has a harmonic dependence on time.
! This subroutine calculates initial data for a boson star
! in the polar-areal gauge (PA) using a shooting method.
! (There is in fact an option to use a relaxation method
! by setting boson_relax=.true. in the parameter file,
! but it is only second order accurate so the shooting
! method is better).

! MAIN PROBLEM: The main problem with this routine is that
! it requires a very fine tuned initial guess since we are
! fighting a growing exponential.  This gets worse as
! we move the boundary further away. So the way to
! proceed is to put the boundary close in, find omega,
! and then move it slowly further out adjusting omega
! to high precision as we move out.

! To obtain the initial data we assume that spacetime is
! static (K_ij=0), and also that the complex scalar field
! has the form:
!
! Phi(t,r) = phi(r) exp(i omega t)
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
! Hamiltonian constraint and the Klein-Gordon equation to
! find (xi=dphi/dr):
!
!             /                      /                      2    2      2    \ \
! dA/dr  =  A | (1 - A)/r + 4 pi r A | 2 V  +  (omega/alpha)  phi  +  xi / A | |
!             \                      \                                       / /
!
!               /                                \
! dxi/dr = - xi | 2/r + d alpha/alpha - d A/(2A) |
!               \        r               r       /
!
!             /                       2 \
!        +  A | V' - phi (omega/alpha)  |
!             \                         /
!
! In the above expressions xi=dphi/dr.  Also V(phi) is the
! self-interaction potential and V'=dV/dphi.
!
! For the lapse we use the polar slicing condition
! K_{theta,theta}=0, which implies:
!
!                   /                                   \
! dalpha/dr = alpha | (A - 1)/r + d A/(2A) - 8 pi r A V |
!                   \              r                    /
!
! Notice that for a static solution this should be equivalent
! to maximal slicing but is easier to solve.
!
! This routine takes as input parameter the value of the scalar
! field at the origin "boson_phi0", and uses a shooting method
! (with Runge-Kutta) to determine the value of the frequency
! "boson_omega" that is compatible with the asymptotically
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
! Notice that in the final solution we don't want alpha(r=0)=1,
! but rather alpha=1 at infinity.  But this is no problem as the
! slicing condition above is linear in alpha so we can always
! just rescale the lapse at the end, but in order not to affect
! the solution for A and phi we must also rescale the final value
! of the frequency omega by the same factor.
!
! ANGULAR MOMENTUM:  It is possible to have spherically symmetric
! configurations of a complex scalar field if we have a linear
! combination of modes (spherical harmonics) with a fixed azimuthal
! quantum number l and all possible values of the quantum number
! m=-l,...,+l, with the amplitude.
!
! In that case, the effect is simply to add a centrifugal
! potential term of the form:
!
! V_l  =  l(l+1) phi^2 / (2 r^2)
!
! But we also need to change the boundary conditions at the
! origin to:
!
! phi(r=0)  ~  phi0 r**l
!
! The sqrt(2*l+1) is a normalization factor, which can also
! be changed for convenience to:
!
! phi(r=0)  ~  phi0 r**l / sqrt(4 pi)
!
! For l>=2 it turns out that the fourth order Runge-Kutta
! is inconsistent at the origin, so in that case we solve
! instead for the rescaled scalar field phi_res = phi/r**l.

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
  real(8) phi0,xi0,piI0                 ! Initial values of variables.
  real(8) phires0,xires0,DV0            ! Initial values of variables.
  real(8) k11,k12,k13,k14               ! Runge-Kutta sources for A.
  real(8) k21,k22,k23,k24               ! Runge-Kutta sources for alpha.
  real(8) k31,k32,k33,k34               ! Runge-Kutta sources for phi.
  real(8) k41,k42,k43,k44               ! Runge-Kutta sources for xi.
  real(8) k51,k52,k53,k54               ! Runge-Kutta sources for phires.
  real(8) k61,k62,k63,k64               ! Runge-Kutta sources for xires.
  real(8) A_rk,alpha_rk                 ! Runge-Kutta values of variables.
  real(8) phi_rk,xi_rk,piI_rk           ! Runge-Kutta values of variables.
  real(8) phires_rk,xires_rk            ! Runge-Kutta values of variables.
  real(8) V_rk,VL_rk,K_rk               ! Runge-Kutta values of variables.
  real(8) J1_PA,J2_PA,J3_PA             ! Functions for sources of differential equations.
  real(8) J4_PA,J5_PA,J6_PA             ! Functions for sources of differential equations.
  real(8) V_PA,DV_PA                    ! Functions for scalar field potential and derivative.
  real(8) res,res_old                   ! Residual.
  real(8) omega_new,omega_old,domega    ! Trial frequency and frequency interval.
  real(8) half,smallpi                  ! Numbers.
  real(8) rm,alphafac,aux               ! Auxiliary quantities.
  real(8) :: epsilon = 1.d-8            ! Tolerance.
 
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                 ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g,alpha_g        ! Radial metric and lapse global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: phi_g,xi_g,piI_g   ! Scalar field global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: phires_g,xires_g   ! Rescaled scalar field arrays.


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  smallpi = acos(-1.d0)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for boson star in polar-areal gauge using shooting method ...'
     print *
  end if

! Sanity check.

  if (spacetime=="minkowski") then
     print *, 'Boson star initial data is not compatible with a Minkowski background ...'
     print *, 'Aborting! (subroutine idata_BosonstarPA)'
     print *
     call die
  end if


! *************************
! ***   NORMALIZATION   ***
! *************************

! Physical normalization.

  if (boson_factor=="physical") then

     if (rank==0) then
        if (complex_l==0) then
           write(*,'(A)') ' Using "physical" normalization at origin: phi(r=0) = phi0'
           print *
        else
           write(*,'(A)') ' Using "physical" normalization at origin: phi ~ phi0 r**l'
           print *
        end if
     end if

! Alternative normalization.

  else

     boson_phi0 = boson_phi0/sqrt(dble(2*complex_l+1))/sqrt(4.d0*smallpi)

     if (rank==0) then
        if (complex_l==0) then
           write(*,'(A,E16.10)') ' Using "harmonic" normalization at origin: phi(r=0) = phi0/sqrt(4pi) = ', boson_phi0
           print *
        else
           write(*,'(A,E16.10)') ' Using "harmonic" normalization at origin: phi/r**l ~ phi0/sqrt(4pi*(2l+1)) = ', boson_phi0
           print *
        end if
     end if

  end if

! Rescale the tolerance with boson_phi0, so the final tolerance
! is measured with respect to this value.

  epsilon = epsilon*boson_phi0


! ************************
! ***   SANITY CHECK   ***
! ************************

! Sanity check for potential.

  if ((complexpotential=="phi2").and.(complex_lambda/=0.d0)) then
     print *, 'You can not have a complexpotential=phi2 and complex_lambda different from zero.'
     print *, 'Aborting! (subroutine idata_BosonstarPA)'
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
  piI_g = 0.d0

  phires_g = 0.d0
  xires_g  = 0.d0


! **********************************
! ***   SOLVE WITH RELAXATION?   ***
! **********************************

! If we want to use the relaxation method from Numerical Recipes we call
! the routine "BosonstarPArelax" (last routine at the end of the file)
! and then return.
!
! Notice that at this point the relaxation routine does not allow for
! perturbations.  We could fix this later, but for the moment this is
! just for testing since the shooting algorithm works well and is fourth
! order, while the relaxation routine is only second order.

  if (boson_relax) then
     if (rank==0) then
        print *, 'Switching to relaxation method instead of shooting ...'
        print *
     end if
     call BosonstarPArelax(rr,A_g,alpha_g,phi_g,xi_g)
     return
  end if


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
        print *, 'For boson star initial data we must have omega_right>omega_left.'
        print *, 'Aborting! (subroutine idata_BosonstarPA)'
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
!          In fact we define:
!
!          DV := V'/r^l
!
!          and evaluate it at r=0, using the fact that:
!
!          phi ~ phi0*r^l

           aux = boson_phi0

           DV0 = DV_PA(aux)

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

              if (complex_l>2) then

                 phi_g(l,imin-1) = (9.d0*(phi_g(l+1,Nrtotal-2)+phi_g(l+1,Nrtotal-3)) &
                                 - (phi_g(l+1,Nrtotal-4)+phi_g(l+1,Nrtotal-1)))/16.d0
                 xi_g(l,imin-1)  = (9.d0*(xi_g(l+1,Nrtotal-2)+xi_g(l+1,Nrtotal-3)) &
                                 - (xi_g(l+1,Nrtotal-4)+xi_g(l+1,Nrtotal-1)))/16.d0

                 phires_g(l,imin-1) = phi_g(l,imin-1)/rr(l,imin-1)**complex_l
                 xires_g(l,imin-1)  = xi_g(l,imin-1)/rr(l,imin-1)**complex_l &
                                    - complex_l*phi_g(l,imin-1)/rr(l,imin-1)**(complex_l+1)

              else

                 phires_g(l,imin-1) = (9.d0*(phires_g(l+1,Nrtotal-2)+phires_g(l+1,Nrtotal-3)) &
                                    - (phires_g(l+1,Nrtotal-4)+phires_g(l+1,Nrtotal-1)))/16.d0
                 xires_g(l,imin-1)  = (9.d0*(xires_g(l+1,Nrtotal-2)+xires_g(l+1,Nrtotal-3)) &
                                    - (xires_g(l+1,Nrtotal-4)+xires_g(l+1,Nrtotal-1)))/16.d0

                 phi_g(l,imin-1) = phires_g(l,imin-1)*rr(l,imin-1)**complex_l
                 xi_g(l,imin-1)  = xires_g(l,imin-1)*rr(l,imin-1)**complex_l &
                                 + complex_l*phires_g(l,imin-1)*rr(l,imin-1)**(complex_l-1)

              end if

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

!                Values of (phires,xires) at origin.

                 phires0 = boson_phi0
                 xires0  = 0.d0

!                Values of (phi,xi) at origin.
 
                 if (complex_l==0) then
                    phi0 = phires0
                    xi0  = 0.d0
                 else if (complex_l==1) then
                    phi0 = 0.d0
                    xi0  = phires0
                 else
                    phi0 = 0.d0
                    xi0  = 0.d0
                 end if

!             Grid spacing and values at previous grid point.

              else

                 delta   = dr(l)
                 r0      = rr(l,i-1)

                 A0      = A_g(l,i-1)
                 alpha0  = alpha_g(l,i-1)

                 phi0    = phi_g(l,i-1)
                 xi0     = xi_g(l,i-1)

                 phires0 = phires_g(l,i-1)
                 xires0  = xires_g(l,i-1)

              end if

!             I) First Runge-Kutta step.

!             Sources at first grid point if we start
!             from the origin (for finer grid).

              if (i==1) then

!                At the origin we have:  A' = alpha' = 0.

                 k11 = 0.d0
                 k21 = 0.d0

!                For l=0 the first derivative of phi vanishes
!                at r=0, and its second derivative goes as:
!
!                phi'' = (V' - omega**2 phi0)/3

                 if (complex_l==0) then

                    k31 = 0.d0
                    K41 = delta*(DV0 - boson_omega**2*boson_phi0)/3.d0

!                For l=1 the first derivative of phi at r=0 is given by:
!
!                phi' = phi0
!
!                and the second derivative vanishes.

                 else if (complex_l==1) then

                    k31 = boson_phi0*delta
                    k41 = 0.d0

!                For l>1 we integrate phires instead of phi.
!                The first derivative of phires vanishes at r=0,
!                and its second derivative goes as:
!
!                phires'' = (V'/r**l - omega**2 phi0)/(2l+3)

                 else

                    k51 = 0.d0
                    k61 = delta*(DV0 - boson_omega**2*boson_phi0)/dble(2*complex_l+3)

                 end if

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk      = A0
                 alpha_rk  = alpha0

                 phi_rk    = phi0
                 xi_rk     = xi0
                 phires_rk = phires0
                 xires_rk  = xires0

!                Physical potential V, centrifugal potential VL, and kinetic term K.

                 V_rk = V_PA(phi_rk)

                 if (complex_l>0) then
                    VL_rk = 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk**2
                 endif

                 K_rk = 0.5d0*(xi_rk**2/A_rk + (phi_rk*boson_omega/alpha_rk)**2)

!                Sources.

                 k11 = delta*J1_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)
                 k21 = delta*J2_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)

                 k31 = delta*J3_PA(xi_rk)
                 k41 = delta*J4_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm)

                 k51 = delta*J5_PA(xires_rk)
                 k61 = delta*J6_PA(A_rk,alpha_rk,phires_rk,xires_rk,rm)

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk     = A0     + half*k11
              alpha_rk = alpha0 + half*k21

              if (complex_l<2) then
                 phi_rk    = phi0 + half*k31
                 xi_rk     = xi0  + half*k41
                 phires_rk = phi_rk/rm**complex_l
                 xires_rk  = xi_rk/rm**complex_l &
                           - complex_l*phi_rk/rm**(complex_l+1)
              else
                 phires_rk = phires0 + half*k51
                 xires_rk  = xires0  + half*k61
                 phi_rk    = phires_rk*rm**complex_l
                 xi_rk     = xires_rk*rm**complex_l &
                           + complex_l*phires_rk*rm**(complex_l-1)
              end if

!             Physical potential V, centrifugal potential VL, and kinetic term K.

              V_rk = V_PA(phi_rk)

              if (complex_l>0) then
                 VL_rk = 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk**2
              end if

              K_rk = 0.5d0*(xi_rk**2/A_rk + (phi_rk*boson_omega/alpha_rk)**2)

!             Sources.

              k12 = delta*J1_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)
              k22 = delta*J2_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)

              k32 = delta*J3_PA(xi_rk)
              k42 = delta*J4_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm)

              k52 = delta*J5_PA(xires_rk)
              k62 = delta*J6_PA(A_rk,alpha_rk,phires_rk,xires_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk     = A0     + half*k12
              alpha_rk = alpha0 + half*k22

              if (complex_l<2) then
                 phi_rk    = phi0 + half*k32
                 xi_rk     = xi0  + half*k42
                 phires_rk = phi_rk/rm**complex_l
                 xires_rk  = xi_rk/rm**complex_l &
                           - complex_l*phi_rk/rm**(complex_l+1)
              else
                 phires_rk = phires0 + half*k52
                 xires_rk  = xires0  + half*k62
                 phi_rk    = phires_rk*rm**complex_l
                 xi_rk     = xires_rk*rm**complex_l &
                           + complex_l*phires_rk*rm**(complex_l-1)
              end if

!             Physical potential V, centrifugal potential VL, and kinetic term K.

              V_rk = V_PA(phi_rk)

              if (complex_l>0) then
                 VL_rk = 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk**2
              end if

              K_rk = 0.5d0*(xi_rk**2/A_rk + (phi_rk*boson_omega/alpha_rk)**2)

!             Sources.

              k13 = delta*J1_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)
              k23 = delta*J2_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)

              k33 = delta*J3_PA(xi_rk)
              k43 = delta*J4_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm)

              k53 = delta*J5_PA(xires_rk)
              k63 = delta*J6_PA(A_rk,alpha_rk,phires_rk,xires_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk     = A0     + k13
              alpha_rk = alpha0 + k23

              if (complex_l<2) then
                 phi_rk    = phi0 + k33
                 xi_rk     = xi0  + k43
                 phires_rk = phi_rk/rm**complex_l
                 xires_rk  = xi_rk/rm**complex_l &
                           - complex_l*phi_rk/rm**(complex_l+1)
              else
                 phires_rk = phires0 + k53
                 xires_rk  = xires0  + k63
                 phi_rk    = phires_rk*rm**complex_l
                 xi_rk     = xires_rk*rm**complex_l &
                           + complex_l*phires_rk*rm**(complex_l-1)
              end if

!             Physical potential V, centrifugal potential VL, and kinetic term K.

              V_rk = V_PA(phi_rk)

              if (complex_l>0) then
                 VL_rk = 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk**2
              end if

              K_rk = 0.5d0*(xi_rk**2/A_rk + (phi_rk*boson_omega/alpha_rk)**2)

!             Sources.

              k14 = delta*J1_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)
              k24 = delta*J2_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)

              k34 = delta*J3_PA(xi_rk)
              k44 = delta*J4_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm)

              k54 = delta*J5_PA(xires_rk)
              k64 = delta*J6_PA(A_rk,alpha_rk,phires_rk,xires_rk,rm)

!             Advance variables to next grid point.

              A_g(l,i)     = A0     + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i) = alpha0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

              if (complex_l<2) then
                 phi_g(l,i)    = phi0   + (k31 + 2.d0*(k32 + k33) + k34)/6.d0
                 xi_g(l,i)     = xi0    + (k41 + 2.d0*(k42 + k43) + k44)/6.d0
                 phires_g(l,i) = phi_g(l,i)/rm**complex_l
                 xires_g(l,i)  = xi_g(l,i)/rm**complex_l &
                               - complex_l*phi_g(l,i)/rm**(complex_l+1)
              else
                 phires_g(l,i) = phires0 + (k51 + 2.d0*(k52 + k53) + k54)/6.d0
                 xires_g(l,i)  = xires0  + (k61 + 2.d0*(k62 + k63) + k64)/6.d0
                 phi_g(l,i)    = phires_g(l,i)*rm**complex_l
                 xi_g(l,i)     = xires_g(l,i)*rm**complex_l &
                               + complex_l*phires_g(l,i)*rm**(complex_l-1)
              end if

!             Check if solution is blowing up. This helps to reduce
!             the need for a very fine tuned initial guess.

              if (abs(phires_g(l,i))>2.d0*boson_phi0) then
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

              if (abs(phires_g(l,i))+abs(phires_g(l,i-1))<epsilon) then
                 phires_g(l,i) = 0.d0
                 xires_g(l,i) = 0.d0
                 phi_g(l,i) = 0.d0
                 xi_g(l,i) = 0.d0
              end if

           end do


!          ***********************
!          ***   GHOST ZONES   ***
!          ***********************

           do i=1,ghost

              A_g(l,1-i)      = + A_g(l,i)
              alpha_g(l,1-i)  = + alpha_g(l,i)

              phi_g(l,1-i)    = + (1.d0-2.d0*mod(complex_l,2))*phi_g(l,i)
              xi_g(l,1-i)     = - (1.d0-2.d0*mod(complex_l,2))*xi_g(l,i)

              phires_g(l,1-i) = + phires_g(l,i)
              xires_g(l,1-i)  = - xires_g(l,i)

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
!       Notice that far away we must have:  k = sqrt(m**2 - (omega/alpha)**2)

        res_old = res

        if (abs(phires_g(0,Nrtotal))+abs(phires_g(0,Nrtotal-1))<epsilon) then
           res = epsilon/2.d0
           goto 100
        else
           aux = abs(complex_mass**2 - (boson_omega/alpha_g(0,Nrtotal))**2)
           res = xi_g(0,Nrtotal) + dsqrt(aux)*phi_g(0,Nrtotal)
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

           if (complex_l<2) then

              phi_g(l-1,iaux) = (9.d0*(phi_g(l,i)+phi_g(l,i+1)) - (phi_g(l,i-1)+phi_g(l,i+2)))/16.d0
              xi_g(l-1,iaux)  = (9.d0*(xi_g(l,i)+xi_g(l,i+1))   - (xi_g(l,i-1)+xi_g(l,i+2)))/16.d0

              phires_g(l-1,iaux) = phi_g(l-1,iaux)/rm**complex_l
              xires_g(l-1,iaux)  = xi_g(l-1,iaux)/rm**complex_l &
                                 - complex_l*phi_g(l-1,iaux)/rm**(complex_l+1)

           else

              phires_g(l-1,iaux) = (9.d0*(phires_g(l,i)+phires_g(l,i+1)) - (phires_g(l,i-1)+phires_g(l,i+2)))/16.d0
              xires_g(l-1,iaux)  = (9.d0*(xires_g(l,i)+xires_g(l,i+1))   - (xires_g(l,i-1)+xires_g(l,i+2)))/16.d0

              phi_g(l-1,iaux) = phires_g(l-1,iaux)*rm**complex_l
              xi_g(l-1,iaux)  = xires_g(l-1,iaux)*rm**complex_l &
                              + complex_l*phires_g(l-1,iaux)*rm**(complex_l-1)

           end if

        end do

!       Fix ghost zones.

        do i=1,ghost

           A_g(l-1,1-i)      = + A_g(l-1,i)
           alpha_g(l-1,1-i)  = + alpha_g(l-1,i)

           phi_g(l-1,1-i)    = + (1.d0-2.d0*mod(complex_l,2))*phi_g(l-1,i)
           xi_g(l-1,1-i)     = - (1.d0-2.d0*mod(complex_l,2))*xi_g(l-1,i)

           phires_g(l-1,1-i) = + phires_g(l-1,i)
           xires_g(l-1,1-i)  = - xires_g(l-1,i)

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

     omega_new = boson_omega/alphafac

     write(*,'(A,E22.16)') ' Omega (not-rescaled) = ', boson_omega
     write(*,'(A,E22.16)') ' Omega (rescaled)     = ', omega_new
     print *


!    ********************************
!    ***   IMAGINARY PART OF pi   ***
!    ********************************

!    From the original ansatz, we take the time derivative of
!    the imaginary part equal to (omega/alpha)*phi.

     piI_g = (omega_new/alpha_g)*phi_g


!    ***********************************
!    ***   BOSON STAR PERTURBATION   ***
!    ***********************************

!    Here we add a gaussian perturbation to the solution
!    for the scalar field and solve again the hamiltonian
!    constraint for A and the polar slicing condition for
!    the lapse.
!
!    The perturbation should be small, and for complex_l
!    different from 0 we need to be sure that it respects
!    the correct behaviour at the origin, so we add it
!    to the rescaled field (the one divided with r**l).
!
!    We also rescale it with max(phi)/r0**l, to make sure that
!    at r=r0, the amplitude a0 really measures the size
!    of the perturbation in phi as a fraction of its maximum
!    value.
!
!    Also, since the perturbation must be even, we take
!    the sum of two symmetric gaussians.

     if ((bosongauss).and.((boson_phiR_a0/=0.d0).or.(boson_piI_a0/=0.d0))) then

!       Message to screen.

        print *, 'Adding gaussian perturbation to boson star ...'

!       Rescale the amplitude of the gaussian with max(phi)/r0**l.

        aux = maxval(abs(phi_g))/boson_phiR_r0**complex_l

        boson_phiR_a0 = aux*boson_phiR_a0
        boson_piI_a0  = aux*boson_piI_a0

!       Add perturbation to rescaled quantities phires and xires.

        if (boson_phiR_r0==0.d0) then
           phires_g = phires_g + boson_phiR_a0*exp(-rr**2/boson_phiR_s0**2)
           xires_g  = xires_g - 2.d0*rr*boson_phiR_a0/boson_phiR_s0**2*exp(-rr**2/boson_phiR_s0**2)
        else
           phires_g = phires_g + boson_phiR_a0 &
                    *(exp(-(rr-boson_phiR_r0)**2/boson_phiR_s0**2) &
                    + exp(-(rr+boson_phiR_r0)**2/boson_phiR_s0**2))
           xires_g  = xires_g - 2.d0*boson_phiR_a0/boson_phiR_s0**2 &
                   *((rr-boson_phiR_r0)*exp(-(rr-boson_phiR_r0)**2/boson_phiR_s0**2) &
                   + (rr+boson_phiR_r0)*exp(-(rr+boson_phiR_r0)**2/boson_phiR_s0**2))
        end if

!       Now find physical phi and xi.

        phi_g = phires_g*rr**complex_l
        xi_g  = xires_g*rr**complex_l + complex_l*phires_g*rr**(complex_l-1)

!       Add perturbation to piI. This perturbation is also rescaled
!       with (omega/alpha) to make it easier to compare the
!       perturbations in phiR and piI.

        if (boson_piI_r0==0.d0) then
           piI_g = piI_g + rr**complex_l*boson_piI_a0*(omega_new/alpha_g)*exp(-rr**2/boson_piI_s0**2)
        else
           piI_g = piI_g + rr**complex_l*boson_piI_a0*(omega_new/alpha_g) &
                 *(exp(-(rr-boson_piI_r0)**2/boson_piI_s0**2) &
                 + exp(-(rr+boson_piI_r0)**2/boson_piI_s0**2))
        end if

!       Solve again for (A,alpha).

        print *, 'Solving again for (alpha,A) ...'
        print *

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
           end if

!          Here we always use fourth order Runge-Kutta.

           do i=imin,Nrtotal

!             Grid spacing and values at first point
!             if we start from the origin (finer grid).

              if (i==1) then

                 delta = half*dr(l)
                 r0    = 0.d0

                 A0      = 1.d0
                 alpha0  = 1.d0

                 if (complex_l==0) then
                    phi0 = (9.d0*(phi_g(l,0)+phi_g(l,1)) - (phi_g(l,-1)+phi_g(l,2)))/16.d0
                    xi0  = 0.d0
                    piI0 = (9.d0*(piI_g(l,0)+piI_g(l,1)) - (piI_g(l,-1)+piI_g(l,2)))/16.d0
                 else if (complex_l==1) then
                    phi0 = 0.d0
                    xi0  = (9.d0*(xi_g(l,0)+xi_g(l,1)) - (xi_g(l,-1)+xi_g(l,2)))/16.d0
                    piI0 = 0.d0
                 else
                    phi0 = 0.d0
                    xi0  = 0.d0
                    piI0 = 0.d0
                 end if

!             Grid spacing and values at previous grid point.

              else

                 delta = dr(l)
                 r0    = rr(l,i-1)

                 A0     = A_g(l,i-1)
                 alpha0 = alpha_g(l,i-1)
                 phi0   = phi_g(l,i-1)
                 xi0    = xi_g(l,i-1)
                 piI0   = piI_g(l,i-1)

              end if

!             I) First Runge-Kutta step.

              if (i==1) then

                 k11 = 0.d0
                 k21 = 0.d0

              else

                 rm = r0

                 A_rk     = A0
                 alpha_rk = alpha0
                 phi_rk   = phi0
                 xi_rk    = xi0
                 piI_rk   = piI0

!                Physical potential V, centrifugal potential VL, and kinetic term K.

                 V_rk = V_PA(phi_rk)

                 if (complex_l>0) then
                    VL_rk = 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk**2
                 end if

                 K_rk  = 0.5d0*(xi_rk**2/A_rk + piI_rK**2)

!                Sources.

                 k11 = delta*J1_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)
                 k21 = delta*J2_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)

              end if

!             II) Second Runge-Kutta step. We interpolate the values of (phi,xi).

              rm = r0 + half*delta

              A_rk     = A0     + half*k11
              alpha_rk = alpha0 + half*k21

              if (i==1) then  ! Linear interpolation for first point.
                 phi_rk = 0.5d0*(phi0 + phi_g(l,1))
                 xi_rk  = 0.5d0*(xi0 + xi_g(l,1))
                 piI_rk = 0.5d0*(piI0 + piI_g(l,1))
              else            ! Cubic interpolation for the rest.
                 phi_rk = (9.d0*(phi_g(l,i)+phi_g(l,i-1)) - (phi_g(l,i-2)+phi_g(l,i+1)))/16.d0
                 xi_rk  = (9.d0*(xi_g(l,i)+xi_g(l,i-1))   - (xi_g(l,i-2)+xi_g(l,i+1)))/16.d0
                 piI_rk = (9.d0*(piI_g(l,i)+piI_g(l,i-1)) - (piI_g(l,i-2)+piI_g(l,i+1)))/16.d0
              end if

!             Physical potential V, centrifugal potential VL, and kinetic term K.

              V_rk = V_PA(phi_rk)

              if (complex_l>0) then
                 VL_rk = 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk**2
              end if

              K_rk  = 0.5d0*(xi_rk**2/A_rk + piI_rk**2)

!             Sources.

              k12 = delta*J1_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)
              k22 = delta*J2_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)

!             III) Third Runge-Kutta step.

              A_rk     = A0     + half*k12
              alpha_rk = alpha0 + half*k22

!             Physical potential V, centrifugal potential VL, and kinetic term K.

              V_rk = V_PA(phi_rk)

              if (complex_l>0) then
                 VL_rk = 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk**2
              end if

              K_rk  = 0.5d0*(xi_rk**2/A_rk + piI_rk**2)

!             Sources.

              k13 = delta*J1_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)
              k23 = delta*J2_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta
 
              A_rk     = A0     + k13
              alpha_rk = alpha0 + k23

              phi_rk = phi_g(l,i)
              xi_rk  = xi_g(l,i)
              piI_rk = piI_g(l,i)

!             Physical potential V, centrifugal potential VL, and kinetic term K.

              V_rk = V_PA(phi_rk)

              if (complex_l>0) then
                 VL_rk = 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk**2
              end if

              K_rk  = 0.5d0*(xi_rk**2/A_rk + piI_rk**2)

!             Sources.

              k14 = delta*J1_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)
              k24 = delta*J2_PA(A_rk,alpha_rk,phi_rk,xi_rk,rm,V_rk,VL_rk,K_rk)

!             Advance variables to next grid point.

              A_g(l,i)     = A0     + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i) = alpha0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

           end do

!          Ghost zones.

           do i=1,ghost
              A_g(l,1-i)     = A_g(l,i)
              alpha_g(l,1-i) = alpha_g(l,i)
           end do

        end do

!       Restrict solution from fine to coarse grid.

        do l=Nl-1,1,-1

           do i=1,Nrtotal-ghost,2

              iaux = i/2 + 1
              rm = rr(l-1,iaux)

              A_g(l-1,iaux)     = (9.d0*(A_g(l,i)+A_g(l,i+1)) - (A_g(l,i-1)+A_g(l,i+2)))/16.d0
              alpha_g(l-1,iaux) = (9.d0*(alpha_g(l,i)+alpha_g(l,i+1)) - (alpha_g(l,i-1)+alpha_g(l,i+2)))/16.d0

           end do

!          Fix ghost zones.

           do i=1,ghost
              A_g(l-1,1-i)     = A_g(l-1,i)
              alpha_g(l-1,1-i) = alpha_g(l-1,i)
           end do

        end do

!       Rescale lapse again.

        alphafac = alpha_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
            *(25.d0*alpha_g(0,Nrtotal) - 48.d0*alpha_g(0,Nrtotal-1) &
            + 36.d0*alpha_g(0,Nrtotal-2) - 16.d0*alpha_g(0,Nrtotal-3) + 3.d0*alpha_g(0,Nrtotal-4))/3.d0

        alpha_g = alpha_g/alphafac

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

  else

     call distribute(0,Nl-1,A,A_g)
     call distribute(0,Nl-1,alpha,alpha_g)

     call distribute(0,Nl-1,complex_phiR,phi_g)
     call distribute(0,Nl-1,complex_xiR,xi_g)
     call distribute(0,Nl-1,complex_piI,piI_g)

  end if


! **********************************************************
! ***   IMAGINARY PART OF (phi,xi) AND REAL PART OF pi   ***
! **********************************************************

! Set the imaginary part and its spatial derivative to zero.

  complex_phiI = 0.d0
  complex_xiI  = 0.d0

! Set time derivative of real part to zero.

  complex_piR  = 0.d0


! ***********************
! ***   DERIVATIVES   ***
! ***********************

! Derivatives of A and alpha.

  diffvar => A

  do l=0,Nl-1
     D1_A(l,:) = diff1(l,+1)
  end do

  diffvar => alpha

  do l=0,Nl-1
     D1_alpha(l,:) = diff1(l,+1)
  end do

! Auxiliary array to check normalization at origin, only for output.

  if (allocated(complex_phiaux)) then

     complex_phiaux = complex_phiR/r**complex_l

     if (boson_factor=="harmonic") then
        complex_phiaux = complex_phiaux*sqrt(4.d0*smallpi*dble(2*complex_l+1))
     end if

  end if

  200  continue


! ***************
! ***   END   ***
! ***************

  end subroutine idata_BosonstarPA







! **********************************
! ***   RADIAL DERIVATIVE OF A   ***
! **********************************

! The radial derivative of A comes from the
! Hamiltonian constraint.

  function J1_PA(A,alpha,phi,xi,rm,V,VL,K)

  use param

  implicit none

  real(8) J1_PA
  real(8) A,alpha,phi,xi,rm
  real(8) V,VL,K,rho
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Energy density.

  rho = K + V + VL

!           /                        \
! dA/dr = A | (1-A)/r + 8 pi r A rho |
!           \                        /
!
! Notice that rho=K+V.

  J1_PA = A*((1.d0 - A)/rm + 8.d0*smallpi*rm*A*rho)


! ***************
! ***   END   ***
! ***************

  end function J1_PA







! **************************************
! ***   RADIAL DERIVATIVE OF ALPHA   ***
! **************************************

! The radial derivative of alpha comes from
! the polar slicing condition:  dK_{theta,theta}/dt=0.

  function J2_PA(A,alpha,phi,xi,rm,V,VL,K)

  use param

  implicit none

  real(8) J2_PA
  real(8) A,alpha,phi,xi,rm
  real(8) V,VL,K
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

!                   /                                  \
! dalpha/dr = alpha | (A-1)/2r + 4 pi r A (K - V - VL) |
!                   \                                  /
!
! Notice that K-V is equal to SA (the radial stress), which is the source
! term in the polar slicing condition.

  J2_PA = alpha*(0.5d0*(A - 1.d0)/rm + 4.d0*smallpi*rm*A*(K-V-VL))


! ***************
! ***   END   ***
! ***************

  end function J2_PA







! ************************************
! ***   RADIAL DERIVATIVE OF PHI   ***
! ************************************

! This is just the definition of xi:  xi=dphi/dr

  function J3_PA(xi)

  use param

  implicit none

  real(8) J3_PA
  real(8) xi

! dphi/dr = xi.

  J3_PA = xi


! ***************
! ***   END   ***
! ***************

  end function J3_PA







! ***********************************
! ***   RADIAL DERIVATIVE OF XI   ***
! ***********************************

! The radial derivative of xi comes from the
! Klein-Gordon equation.

! This is specialized to the case of no perturbation,
! since for the perturbed solution we don't call this
! function again.

  function J4_PA(A,alpha,phi,xi,rm)

  use param

  implicit none

  real(8) J4_PA
  real(8) A,alpha,phi,xi,rm
  real(8) V,VL,DV,DVL
  real(8) V_PA,DV_PA
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Potential and derivative, including centrifugal terms.

  V  = V_PA(phi)
  DV = DV_PA(phi)

  VL  = 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi**2
  DVL = (dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi

!                 /                 2          \       /                         2     \
! dxi/dr = - xi/r | 1 + A  -  8 pi r  A (V+VL) |  +  A | (V'+VL') - (omega/alpha)  phi |
!                 \                            /       \                               /
!
! Notice that in the expression above, the term with V comes from the derivatives
! of (alpha,A), while the term with V' comes from the Klein-Gordon equation itself.
!
! For the case with non-zero angular momentum, both these terms should include
! contributions from the centrifugal effective potential.

  J4_PA = - xi/rm*(1.d0 + A - 8.d0*smallpi*rm**2*A*(V+VL)) &
        + A*((DV+DVL) - (boson_omega/alpha)**2*phi)


! ***************
! ***   END   ***
! ***************

  end function J4_PA







! ***************************************
! ***   RADIAL DERIVATIVE OF PHIRES   ***
! ***************************************

! This is just the definition of xi:  xi=dphi/dr

  function J5_PA(xires)

  use param

  implicit none

  real(8) J5_PA
  real(8) xires

! dphires/dr = xires.

  J5_PA = xires


! ***************
! ***   END   ***
! ***************

  end function J5_PA







! **************************************
! ***   RADIAL DERIVATIVE OF XIRES   ***
! **************************************

! The radial derivative of xi comes from the
! Klein-Gordon equation.

! This is specialized to the case of no perturbation,
! since for the perturbed solution we don't call this
! function again.

  function J6_PA(A,alpha,phires,xires,rm)

  use param

  implicit none

  real(8) J6_PA
  real(8) A,alpha,rm
  real(8) phires,xires
  real(8) V,VL,DV
  real(8) V_PA,DV_PA
  real(8) smallpi,aux

! Numbers.

  smallpi = acos(-1.d0)

! Potential and derivative, including centrifugal terms.

  aux = phires*rm**complex_l

  V  = V_PA(aux)
  DV = DV_PA(aux)

  VL = 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*aux**2

!                2              2
! dxires/dr  =  l (A-1) phires/r  -  (1+A+2l) xires/r
!
!
!            + 8 pi A (V+VL) (r xires + l*phires)
!
!                 /     l                 2        \
!            +  A | V'/r  -  (omega/alpha)  phires |
!                 \                                /
!
! Notice that in the expression above, the term with V comes from the derivatives
! of (alpha,A), while the term with V' comes from the Klein-Gordon equation itself.
!
! For the case with non-zero angular momentum the V term includes contributions
! from the centrifugal effective potential, while the V' terms does not since the
! centrifugal contribution has been used to make the term proportional to phires/r**2
! (the first one) manifestly regular.

  J6_PA = dble(complex_l)**2*(A - 1.d0)*phires/rm**2 &
        - (1.d0 + A + 2.d0*dble(complex_l))*xires/rm &
        + 8.d0*smallpi*A*(V+VL)*(rm*xires + dble(complex_l)*phires) &
        + A*(DV/rm**complex_l - (boson_omega/alpha)**2*phires)


! ***************
! ***   END   ***
! ***************

  end function J6_PA







! **********************************
! ***   SCALAR FIELD POTENTIAL   ***
! **********************************

! Scalar field potential.

  function V_PA(phi)

  use param

  implicit none

  real(8) V_PA
  real(8) phi

! By default set potential to zero.

  V_PA = 0.d0

! Quadratic potential.

  if (complexpotential=="phi2") then

     V_PA = 0.5d0*complex_mass**2*phi**2

! Phi^4 potential.

  else if (complexpotential=="phi4") then

     V_PA = 0.5d0*complex_mass**2*phi**2 + 0.25d0*complex_lambda*phi**4

! Unknown potential.

  else

     print *, 'Unknown potential type ...'
     print *, 'Aborting! (subroutine idata_BosonstarPA)'
     print *
     call die

  end if


! ***************
! ***   END   ***
! ***************

  end function V_PA







! *********************************************
! ***   SCALAR FIELD POTENTIAL DERIVATIVE   ***
! *********************************************

! Scalar field potential.

  function DV_PA(phi)

  use param

  implicit none

  real(8) DV_PA
  real(8) phi

! By default set potential to zero.

  DV_PA = 0.d0

! Quadratic potential.

  if (complexpotential=="phi2") then

     DV_PA = complex_mass**2*phi

! Phi^4 potential.

  else if (complexpotential=="phi4") then

     DV_PA = complex_mass**2*phi + complex_lambda*phi**3

! Unknown potential.

  else

     print *, 'Unknown potential type ...'
     print *, 'Aborting! (subroutine idata_BosonstarPA)'
     print *
     call die

  end if


! ***************
! ***   END   ***
! ***************

  end function DV_PA







  subroutine BosonstarPArelax(rr,A_g,alpha_g,phi_g,xi_g)

! *****************************************************************
! ***   SOLVE THE EIGENVALUE PROBLEM WITH A RELAXATION METHOD   ***
! *****************************************************************

! Here we solve the eigenvalue problem using the routine "solvde"
! from Numerical Recipes, which solves the equations with a relaxation
! that uses Newton's method.
!
! The routine "solvde" calls the subroutine "difeq", where one must define
! the equations to solve and the derivatives for Newton's method.  Have a look
! at it, it is in the directory "base".
!
! The routine "solvde" is adapated to first order systems, so here I take:
!
! y1  =  A
! y2  =  alpha
! y3  =  phi
! y4  =  dphi/dr  =  xi
! y5  =  omega (eigenvalue)
!
! This section is just for testing, as it is only second order.
! If you want to use it you must set "boson_relax = .true." in the
! paremeter file.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Include Numerical Recipes modules.

  use nrtype
  use nr, only : solvde

! Extra variables.

  implicit none

  integer i,i0

  real(8) r0,delta,ddr
  real(8) alphafac,omega_new
  real(8) aux

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                 ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g,alpha_g        ! Radial metric and lapse global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: phi_g,xi_g,piI_g   ! Scalar field global arrays.

! Variables passed to solvde (I4B and DP are Numerical Recipes types).

  integer(I4B) itmax,nb,l
  integer(I4B), dimension(5) :: indexv

  real(DP) conv,slowc
  real(DP), dimension(5) :: scalv
  real(DP), dimension(1:Nrtotal) :: x    ! Radial position for solvde.
  real(DP), dimension(5,1:Nrtotal) :: y  ! Solution vector for solvde.


! **************************************************
! ***   CALL RELAXATION ROUTINE ON COARSE GRID   ***
! **************************************************

! Only processor 0 solves the system of equation.

  if (rank==0) then

!    Message to screen.

     print *, 'Solving the system using the "solvde" routine from Numerical Recipes ...'
     print *

!    We have four boundary conditions at left hand side:
!
!    A(r=0)     =  1
!    alpha(r=0) =  1
!    phi(r=0)   =  phi0
!    xi(r=0)    =  0

     nb = 4

!    Set "solvde" parameters (slowc,scalv) to 1 (the default).

     slowc = 1.d0
     scalv = 1.d0

!    Maximum number of iterations and tolerance.

     itmax = 100
     conv  = 1.d-10

!    Trivial ordering.

     indexv(1) = 1
     indexv(2) = 2
     indexv(3) = 3
     indexv(4) = 4
     indexv(5) = 5

!    Initialize (alpha,A) to 1.

     y(1,:) = 1.d0
     y(2,:) = 1.d0

!    Initial guess for (phi,xi).  For the moment a simple
!    gaussian. The width of the gaussian must be inverse to the
!    amplitude phi0. The expresion below works reasonably well.
!
!    I will improve this later and allow one to read a previous
!    solution as initial guess.

     aux = complexR_s0/boson_phi0

     y(3,:) = boson_phi0*exp(-rr(0,:)**2/aux**2)
     y(4,:) = - 2.d0*rr(0,:)*y(3,:)/aux**2

!    Initialize eigenvalue to average between omega_left and omega_right.

     y(5,:) = 0.5d0*(omega_left+omega_right)

!    Copy radial position from rr to x.

     do i=1,Nrtotal
        x(i) = rr(0,i)
     end do

!    Call "solvde" routine.

     call solvde(itmax,conv,slowc,scalv,indexv,nb,x,y,0)
     print *

!    Copy solution to arrays (A,alpha,phi,xi).

     do i=1,Nrtotal
        A_g(0,i)     = y(1,i)
        alpha_g(0,i) = y(2,i)
        phi_g(0,i)   = y(3,i)
        xi_g(0,i)    = y(4,i)
     end do

!    The new frequency must be the value of y(5,:),
!    which should be constant.

     boson_omega = y(5,1)

!    Ghost zones.

     do i=1,ghost
        A_g(0,1-i)     = + A_g(0,i)
        alpha_g(0,1-i) = + alpha_g(0,i)
        phi_g(0,1-i)   = + phi_g(0,i)
        xi_g(0,1-i)    = - xi_g(0,i)
     end do


!    *********************************
!    ***   RESCALE (ALPHA,OMEGA)   ***
!    *********************************

!    Normalize lapse and omega.  Since we integrated initially
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

     omega_new = boson_omega/alphafac

     write(*,'(A,E22.16)') ' Omega (not-rescaled) = ', boson_omega
     write(*,'(A,E22.16)') ' Omega (rescaled)     = ', omega_new
     print *


!    ************************
!    ***   FINNER GRIDS   ***
!    ************************

!    At this point we only have the solution on the coarse
!    grid.  If we have finner grids we need to find the
!    solution on those grids.
!
!    Here we just interpolate.  This is not particularly
!    accurate, but it is easy and works. It would be much
!    better to solve on the finner grids using Dirichlet
!    boundary conditions from the coarse grid.  I might
!    try this later.

     if (Nl>1) then

        do l=1,Nl-1

!          Grid spacing.

           ddr = dr(l-1)

!          Interpolation.

           do i=1,Nrtotal

!             Position to interpolate.

              r0 = rr(l,i)

!             Find grid point to the left of r0.

              i0 = int((r0-rr(l-1,0))/ddr)
              delta = r0 - rr(l-1,i0)
              aux = delta/ddr

!             Cubic interpolation.

              A_g(l,i) = A_g(l-1,i0) - aux*(A_g(l-1,i0+2) - 6.d0*A_g(l-1,i0+1) &
                 + 3.d0*A_g(l-1,i0) + 2.d0*A_g(l-1,i0-1))/6.d0 &
                 + aux**2*(3.d0*A_g(l-1,i0+1) - 6.d0*A_g(l-1,i0) + 3.d0*A_g(l-1,i0-1))/6.d0 &
                 + aux**3*(A_g(l-1,i0+2) - 3.d0*A_g(l-1,i0+1) + 3.d0*A_g(l-1,i0) - A_g(l-1,i0-1))/6.d0

              alpha_g(l,i) = alpha_g(l-1,i0) + aux*(alpha_g(l-1,i0+1) - alpha_g(l-1,i0))

              alpha_g(l,i) = alpha_g(l-1,i0) - aux*(alpha_g(l-1,i0+2) - 6.d0*alpha_g(l-1,i0+1) &
                 + 3.d0*alpha_g(l-1,i0) + 2.d0*alpha_g(l-1,i0-1))/6.d0 &
                 + aux**2*(3.d0*alpha_g(l-1,i0+1) - 6.d0*alpha_g(l-1,i0) + 3.d0*alpha_g(l-1,i0-1))/6.d0 &
                 + aux**3*(alpha_g(l-1,i0+2) - 3.d0*alpha_g(l-1,i0+1) + 3.d0*alpha_g(l-1,i0) - alpha_g(l-1,i0-1))/6.d0

              phi_g(l,i) = phi_g(l-1,i0) - aux*(phi_g(l-1,i0+2) - 6.d0*phi_g(l-1,i0+1) &
                 + 3.d0*phi_g(l-1,i0) + 2.d0*phi_g(l-1,i0-1))/6.d0 &
                 + aux**2*(3.d0*phi_g(l-1,i0+1) - 6.d0*phi_g(l-1,i0) + 3.d0*phi_g(l-1,i0-1))/6.d0 &
                 + aux**3*(phi_g(l-1,i0+2) - 3.d0*phi_g(l-1,i0+1) + 3.d0*phi_g(l-1,i0) - phi_g(l-1,i0-1))/6.d0

              xi_g(l,i) = xi_g(l-1,i0) - aux*(xi_g(l-1,i0+2) - 6.d0*xi_g(l-1,i0+1) &
                 + 3.d0*xi_g(l-1,i0) + 2.d0*xi_g(l-1,i0-1))/6.d0 &
                 + aux**2*(3.d0*xi_g(l-1,i0+1) - 6.d0*xi_g(l-1,i0) + 3.d0*xi_g(l-1,i0-1))/6.d0 &
                 + aux**3*(xi_g(l-1,i0+2) - 3.d0*xi_g(l-1,i0+1) + 3.d0*xi_g(l-1,i0) - xi_g(l-1,i0-1))/6.d0

           end do

!          Ghost zones.

           do i=1,ghost
              A_g(l,1-i)     = + A_g(l,i)
              alpha_g(l,1-i) = + alpha_g(l,i)
              phi_g(l,1-i)   = + phi_g(l,i)
              xi_g(l,1-i)    = - xi_g(l,i)
           end do

        end do

     end if


!    ********************************
!    ***   IMAGINARY PART OF pi   ***
!    ********************************

!    From the original ansatz, we take the time derivative of
!    the imaginary part equal to (omega/alpha)*phi.

     piI_g = (omega_new/alpha_g)*phi_g


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

  else

     call distribute(0,Nl-1,A,A_g)
     call distribute(0,Nl-1,alpha,alpha_g)

     call distribute(0,Nl-1,complex_phiR,phi_g)
     call distribute(0,Nl-1,complex_xiR,xi_g)
     call distribute(0,Nl-1,complex_piI,piI_g)

  end if


! **********************************************************
! ***   IMAGINARY PART OF (phi,xi) AND REAL PART OF pi   ***
! **********************************************************

! Set the imaginary part of phi and its spatial derivative to zero.

  complex_phiI = 0.d0
  complex_xiI  = 0.d0

! Set time derivative of the real part of phi to zero.

  complex_piR  = 0.d0


! ***************
! ***   END   ***
! ***************

  end subroutine BosonstarPArelax


