
  subroutine idata_TOVstar

! *********************************
! ***   TOV STAR INITIAL DATA   ***
! *********************************

! This subroutine calculates initial data for a TOV star.
! In order to find the initial data we need to solve for
! the mass function m(r) (which gives us the radial metric),
! and the pressure of the star.
!
! To obtain the initial data we assume that spacetime is
! static (K_ij=0).  For the spatial metric we take:
!
!   2          2   2        2     2      2
! ds  = - alpha  dt  +  A dr  +  r dOmega
!
! In other words, we are assuming psi=B=1, beta=0, and
! the radial coordinate is the area radius.
!
! We substitute the above metric and this ansatz into the
! Hamiltonian constraint to find an equation for the radial
! mertric A(r) of the form:
!
!                      
! dA/dr  =  A [ (1 - A)/r + 8 pi r A rho ]
!
!
! with "rho" the total energy density, which for a fluid is
! given by:
!
! rho  =  rho0 ( 1 + e )
!
! where now "rho0" is the rest mass energy density (fluid_rho)
! and "e" is the specific internal energy (fluid_e).
!
! On the other hand, from the TOV equations we find the
! following equation for the pressure "p":
!
!                                 2
! dp/dr  =  - (rho + p) ( m(r) / r  +  4 pi r p ) / ( 1 - 2m(r) / r )
!
!
! with "m(r)" the mass function defined as:
!
! m(r)  =  (r/2) ( 1 - 1/A )
!
! In order to solve the system we still need an equation of state.
! Here we assume we have a polytropic,equation of state:
!
! p  =  Kappa rho0^gamma
!
! which in turn implies that the specific internal energy is:
!
! e  =  [ Kappa/(gamma-1) ] rho0^(gamma-1)
!
! Using this equation of state, we can substitute in the equation
! for the pressure "p" above to find:
!
!
! d rho0  =  - rho0 [ rho0^(1-gamma) / Kappa  +  gamma / (1 - gamma) ]
!  r
!                   2
!         * ( m / r  +  4 pi r Kappa rho0^gamma ) / [ gamma ( 1 - 2 m / r) ]
!
!
! In principle we should also solve for the lapse, but we really
! don't need to since it does not couple to the two equations above,
! so that we can just solve the maximal slicing equation at the
! end.  This is much better since this way the lapse remains
! smooth at the surface of the star.
!
! SOME SPECIAL TEST CASES:
!
! 1) gamma = 2, kappa = 100, rho0 = 0.00128
!
!    This is a TOV star with a stiff equation of state (gamma=2),
!    and should result in a star with mass M=1.4 and radius R=9.58
!    in code units.  If we assume that the units are such that
!    c=C=Msol=1, then the radius is R=14.14 km.
!
! 2) gamma=4/3, kappa=0.4654, rho0=1.d-10 (very small)
!
!    This is a TOV star with an ultra-relativistic fluid (gamma=4/3)
!    but newtonian gravity (very small rho0, so that alpha~A~1).
!    It is the Chandrasekhar limit, so it corresponds to a mass
!    of M=1.44 (in units such that c=G=Msol=1).  The mass should
!    not depend of rho0 as long as it is very small, but the
!    radius does depend on rho0 (for rho0=1.d-10) it should be
!    R~5718 or R~8446 km (it can vary a little depending resolution).
!    The value of kappa=0.4654 is obtained from the corresponding
!    Lane-Emden solution taking cG=Msol=1).

! NOTE FOR PARALLEL RUNS:  The initial data is not really
! solved in parallel.  It is in fact solved only on processor
! zero on a full size array, and then it is distributed
! among the other processors.  This is slow, but works.

! Include modules.

  use procinfo
  use param
  use arrays
  use derivatives

! Extra variables.

  implicit none

  logical foundsurf                     ! Did we find the star's surface?

  integer i,l,iter                      ! Counters.
  integer imin                          ! Leftmost grid point.
  integer iaux                          ! Auxiliary quantity.
  integer irad                          ! Position of star's surface.

  real(8) r0,delta                      ! Local radius and grid spacing.
  real(8) A0,rho0                       ! Initial values of variables.
  real(8) k11,k12,k13,k14               ! Runge-Kutta sources for A.
  real(8) k21,k22,k23,k24               ! Runge-Kutta sources for rho0.
  real(8) A_rk,rho0_rk                  ! Runge-Kutta values of variables.
  real(8) J1_TOV,J2_TOV                 ! Functions for sources of differential equations.
  real(8) rm,aux                        ! Auxiliary quantities.
  real(8) half,smallpi                  ! Numbers.
  real(8) rhoatmos,Eatmos,patmos        ! Atmosphere values.
  real(8) c,G,Msol,Rsol

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr            ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g           ! Radial metric global array.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rho0_g        ! Rest mass density global array.


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  smallpi = acos(-1.d0)

! Speed of light in SI.

  c = 2.99792458d8

! Newton's constant in SI.

  G = 6.6743d-11

! Solar mass in SI.

  Msol = 1.988920d30

! Solar radius in meters.

  Rsol = G*Msol/c**2
  !print *, Rsol


! ******************************
! ***   SURFACE PARAMETERS   ***
! ******************************

! Initialize flag and surface position.

  foundsurf = .false.
  irad = 0


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, "Solving initial data for a TOV star ..."
  end if


! ************************
! ***   SANITY CHECK   ***
! ************************

! Sanity check: gamma must be larger than 1,
! and kappa must be positive.

  if (fluid_gamma<=1.0) then
     print *
     print *, "The adiabatic index gamma must be larger than 1."
     print *, "Aborting in 'idata_TOV.f90' ..."
     print *
     call die
  end if

  if (fluid_kappa<0.0) then
     print *
     print *, "The polytropic constant kappa must be non-negative."
     print *, "Aborting in 'idata_TOV.f90' ..."
     print *
     call die
  end if

! TOV_rad should be zero at this point.

  if (TOV_rad/=0.d0) then
     print *
     print *, "The value of TOV_rad should not be set in the parameter file."
     print *, "Aborting in 'idata_TOV.f90' ..."
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
  
! We solve the system of equations using fourth
! order Runge--Kutta.

! For the moment only processor 0 solves for the
! initial data.  The solution is later distributed
! among processors.

! Initialize arrays.

  A_g = 1.d0
  rho0_g = 0.d0


! *********************************************
! ***   ONLY PROCESSOR 0 SOLVES THE ODE's   ***
! *********************************************

  if (rank==0) then


!    *********************************
!    ***   LOOP OVER GRID LEVELS   ***
!    *********************************

!    We solve from fine to coarse grid.

     do l=Nl-1,0,-1

!       Find initial point. Only the finest grid
!       integrates from the origin.

        if (l==Nl-1) then
           imin = 1
        else
           imin = Nrtotal/2
        end if

!       For coarse grids we interpolate the initial point
!       (which is half way through the grid) from the last
!       points of the higher resolution grid.  Here we use
!       cubic interpolation.

        if (l<Nl-1) then
           A_g(l,imin-1)= (9.d0*(A_g(l+1,Nrtotal-2)+A_g(l+1,Nrtotal-3)) &
                        - (A_g(l+1,Nrtotal-4)+A_g(l+1,Nrtotal-1)))/16.d0
           rho0_g(l,imin-1) = (9.d0*(rho0_g(l+1,Nrtotal-2)+rho0_g(l+1,Nrtotal-3)) &
                            - (rho0_g(l+1,Nrtotal-4)+rho0_g(l+1,Nrtotal-1)))/16.d0
        end if


!       ************************************
!       ***   FOURTH ORDER RUNGE-KUTTA   ***
!       ************************************

        do i=imin,Nrtotal

!          Grid spacing and values at first point
!          if we start from the origin (finer grid).

           if (i==1) then

!             For the first point we use dr/2.

              delta = half*dr(l)
              r0    = 0.d0

!             Values of (A,rho) at origin.

              A0   = 1.d0
              rho0 = TOV_rho0

!          Grid spacing and values at previous grid point.

           else

              delta = dr(l)
              r0    = rr(l,i-1)

              A0   = A_g(l,i-1)
              rho0 = rho0_g(l,i-1)

           end if

!          I) First Runge-Kutta step.

!          Sources at first grid point if we start
!          from the origin (for finer grid).

           if (i==1) then

!             At the origin we have A'=rho'=0.

              k11 = 0.d0
              k21 = 0.d0

!          Sources at previous grid point.

           else

              rm = r0

              A_rk = A0
              rho0_rk = rho0

!             Sources.

              k11 = delta*J1_TOV(A_rk,rho0_rk,rm)

              if (.not.foundsurf) then
                 k21 = delta*J2_TOV(A_rk,rho0_rk,rm)
              else
                 k21 = 0.d0
              end if

           end if

!          II) Second Runge-Kutta step.

           rm = r0 + half*delta

           A_rk    = A0   + half*k11
           rho0_rk = rho0 + half*k21

!          Sources.

           k12 = delta*J1_TOV(A_rk,rho0_rk,rm)

           if (.not.foundsurf) then
              k22 = delta*J2_TOV(A_rk,rho0_rk,rm)
           else
              k22 = 0.d0
           end if

!          III) Third Runge-Kutta step.

           A_rk    = A0   + half*k12
           rho0_rk = rho0 + half*k22

!          Sources.

           k13 = delta*J1_TOV(A_rk,rho0_rk,rm)

           if (.not.foundsurf) then
              k23 = delta*J2_TOV(A_rk,rho0_rk,rm)
           else
              k23 = 0.d0
           end if

!          IV) Fourth Runge-Kutta step.

           rm = r0 + delta

           A_rk    = A0   + k13
           rho0_rk = rho0 + k23

!          Sources.

           k14 = delta*J1_TOV(A_rk,rho0_rk,rm)

           if (.not.foundsurf) then
              k24 = delta*J2_TOV(A_rk,rho0_rk,rm)
           else
              k24 = 0.d0
           end if

!          Advance variables to next grid point.

           A_g(l,i) = A0 + (k11 + 2.d0*(k12 + k13) + k14)/6.d0

           if (.not.foundsurf) then
              rho0_g(l,i) = rho0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0
           else
              rho0_g(l,i) = 0.d0
           end if

!          If the density becomes negative we
!          have reached the surface of the star.

           if (rho0_g(l,i)<=0.d0) then

              if (.not.foundsurf) then

!                Set flag to true and save position
!                of last grid point inside star.

                 foundsurf = .true.

                 irad = i-1                

!                Find star radius using linear extrapolation from
!                previous two points.

                 TOV_rad = rr(l,irad-1) - dr(l)*rho0_g(l,irad-1)/(rho0_g(l,irad) - rho0_g(l,irad-1))

                 write(*,'(A,E19.12)') ' Radius of TOV star = ',TOV_rad

!                We also output the radius in km assuming our units are
!                such that c=G=Msol=1.  For this we need to multiply
!                the result with G*Msol/c^2 ~ 1.477.

                 write(*,'(A,E19.12)') ' Radius of TOV star in km (taking c=G=Msol=1) = ',TOV_rad*Rsol/1.d3
                 print *

!                Set fluid density equal to zero.

                 rho0_g(l,i) = 0.d0

              end if

           end if

        end do


!       ***********************
!       ***   GHOST ZONES   ***
!       ***********************

        do i=1,ghost
           A_g(l,1-i) = A_g(l,i)
           rho0_g(l,1-i) = rho0_g(l,i)
        end do

     end do

!    ****************************************
!    ***   WARNING IF SURFACE NOT FOUND   ***
!    ****************************************

!    Check if we didn't find the surface.

     if (.not.foundsurf) then
        print *
        print *, 'WARNING: Did not reach surface of star, consider moving the outer boundary further away.'
        print *, '         Setting radius to edge of grid ...'
        print *
        TOV_rad = rr(0,Nrtotal)
     end if


!    ************************
!    ***   PERTURBATION   ***
!    ************************

!    When desired, here we add a small gaussian in the
!    density as a perturbation and solve the Hamiltonian
!    constraint again.

     if (TOVgauss) then

!       Message to screen.

        print *, 'Adding gaussian perturbation to TOV star ...'

!       Rescale the amplitude of the gaussian with max(rho).

        aux = maxval(rho0_g)
        TOV_a0 = aux*TOV_a0

!       Add gaussian perturbation to density.

        if (TOV_r0==0.d0) then
           do i=1-ghost,Nrtotal
              rho0_g(:,i) = rho0_g(:,i) + TOV_a0*exp(-rr(:,i)**2/TOV_s0**2)
           end do
        else
           do i=1-ghost,Nrtotal
              rho0_g(:,i) = rho0_g(:,i) + TOV_a0 &
                          *(exp(-(rr(:,i)-TOV_r0)**2/TOV_s0**2) &
                          + exp(-(rr(:,i)+TOV_r0)**2/TOV_s0**2))
           end do
        end if

!       Solve again Hamiltonian constraint.

        print *, 'Solving hamiltonian constraint again'
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
           end if

!          Fourth order Runge-Kutta.

           do i=imin,Nrtotal

!             Grid spacing and values at first point
!             if we start from the origin (finer grid).

              if (i==1) then

                 delta = half*dr(l)
                 r0 = 0.d0

                 A0 = 1.d0
                 rho0 = (9.d0*(rho0_g(l,0)+rho0_g(l,1)) - (rho0_g(l,-1)+rho0_g(l,2)))/16.d0

!             Grid spacing and values at previous grid point.

              else

                 delta = dr(l)
                 r0 = rr(l,i-1)

                 A0 = A_g(l,i-1)
                 rho0 = rho0_g(l,i-1)

              end if

!             I) First Runge-Kutta step.

!             Sources at first grid point if we start
!             from the origin (for finer grid).

              if (i==1) then

!                At the origin we have A'=0

                 k11 = 0.d0

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk = A0
                 rho0_rk = rho0

!                Sources.

                 k11 = delta*J1_TOV(A_rk,rho0_rk,rm)

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk = A0 + half*k11

              if (i==1) then ! Linear interpolation for first point.
                 rho0_rk = 0.5d0*(rho0 + rho0_g(l,1))
              else           ! Cubic interpolation for the rest.
                 rho0_rk = (9.d0*(rho0_g(l,i)+rho0_g(l,i-1)) - (rho0_g(l,i-2)+rho0_g(l,i+1)))/16.d0
              end if

!             Sources.

              k12 = delta*J1_TOV(A_rk,rho0_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk = A0 + half*k12

!             Sources.

              k13 = delta*J1_TOV(A_rk,rho0_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk = A0 + k13
              rho0_rk = rho0_g(l,i)

!             Sources.

              k14 = delta*J1_TOV(A_rk,rho0_rk,rm)

!             Advance to next grid point.

              A_g(l,i) = A0 + (k11 + 2.d0*(k12 + k13) + k14)/6.d0

           end do

!          Ghost zones.

           do i=1,ghost
              A_g(l,1-i) = A_g(l,i)
           end do

        end do

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

           A_g(l-1,iaux) = (9.d0*(A_g(l,i)+A_g(l,i+1)) - (A_g(l,i-1)+A_g(l,i+2)))/16.d0
           rho0_g(l-1,iaux) = (9.d0*(rho0_g(l,i)+rho0_g(l,i+1)) - (rho0_g(l,i-1)+rho0_g(l,i+2)))/16.d0

        end do

!       Fix ghost zones.

        do i=1,ghost
           A_g(l-1,1-i) = A_g(l-1,i)
           rho0_g(l-1,1-i) = rho0_g(l-1,i)
        end do

     end do


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
     A = A_g
     fluid_rho = rho0_g
  else
     call distribute(0,Nl-1,A,A_g)
     call distribute(0,Nl-1,fluid_rho,rho0_g)
  end if


! ******************************************
! ***   FIND ALL OTHER FLUID VARIABLES   ***
! ******************************************

! Set atmosphere.

  rhoatmos = TOV_rho0*fluid_atmos/10.d0

  do l=0,Nl-1
     do i=1-ghost,Nr
        if (fluid_rho(l,i)<=TOV_rho0*fluid_atmos) then
           fluid_rho(l,i) = rhoatmos
        end if
     end do
  end do

! Find pressure and internal energy.

  if (fluid_EOS=="ideal") then
     fluid_p = fluid_kappa*fluid_rho**fluid_gamma
     fluid_e = fluid_kappa*fluid_rho**(fluid_gamma-1.d0)/(fluid_gamma-1.d0)
  end if

! Find enthalpy: h = 1 + e + p/rho.

  fluid_h = 1.d0 + fluid_e + fluid_p/fluid_rho

! Set fluid velocity to zero.

  fluid_v = 0.d0

! Set Lorentz factor to one.

  fluid_W = 1.d0

! Conserved quantities.

  fluid_cD = fluid_rho
  fluid_cE = fluid_rho*fluid_e
  fluid_cS = 0.d0

! Initialize artificial viscosity to zero.

  fluid_q = 0.d0

! Initialize speed of sound for a polytropic equation of state:
!
!   2
! vs  =  gamma (gamma - 1) e / (1 + gamma e)

  fluid_vs = sqrt(abs(fluid_gamma*(fluid_gamma-1.d0)*fluid_e/(1.d0+fluid_gamma*fluid_e)))



! **********************
! ***   FIND LAPSE   ***
! **********************

! Now solve for the lapse using maximal slicing.
! But we first need to calculate all auxiliary
! variables on all grid levels, which include
! the stress-energy tensor.

  do l=0,Nl-1
     call auxiliary(l)
  end do

  call alphamaximal(0,Nl-1,"robin",1.d0)


! ***************
! ***   END   ***
! ***************

  end subroutine idata_TOVstar







! **********************************
! ***   RADIAL DERIVATIVE OF A   ***
! **********************************

! The radial derivative of A comes from the
! Hamiltonian constraint.

  function J1_TOV(A,rho0,rm)

  use param

  implicit none

  real(8) J1_TOV
  real(8) A,rho0,rm
  real(8) e,rhotot
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Specific internal energy.

  e = fluid_kappa/(fluid_gamma - 1.d0)*abs(rho0)**(fluid_gamma-1.d0)

! Total energy density.

  rhotot = rho0*(1.d0 + e)

!           /                        \
! dA/dr = A | (1-A)/r + 8 pi r A rho |
!           \                        /
!
! Notice that rho=K+V.

  J1_TOV = A*((1.d0-A)/rm + 8.d0*smallpi*rm*A*rhotot)

  end function J1_TOV







! *************************************
! ***   RADIAL DERIVATIVE OF rho0   ***
! *************************************

! The radial derivative of rho0 comes from
! the TOV equations.

  function J2_TOV(A,rho0,rm)

  use param

  implicit none

  real(8) J2_TOV
  real(8) A,rho0,rm
  real(8) mass
  real(8) smallpi

! Numbers.

  smallpi = acos(-1.d0)

! Mass function:  m(r) = (r/2) ( 1 - 1/A )

  mass = 0.5d0*rm*(1.d0 - 1.d0/A)

!
! d rho0  =  - rho0 [ rho0^(1-gamma) / Kappa  +  gamma / (1 - gamma) ]
!  r
!                   2
!         * ( m / r  +  4 pi r Kappa rho0^gamma ) / [ gamma ( 1 - 2 m / r) ]

  J2_TOV = - abs(rho0)*(abs(rho0)**(1.d0-fluid_gamma)/fluid_kappa + fluid_gamma/(fluid_gamma - 1.d0)) &
         *(mass/rm**2 + 4.d0*smallpi*rm*fluid_kappa*abs(rho0)**fluid_gamma) &
         /(fluid_gamma*(1.d0 - 2.d0*mass/rm))

  end function J2_TOV
