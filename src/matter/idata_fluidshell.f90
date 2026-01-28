
  subroutine idata_fluidshell

! ************************************
! ***   FLUID SHELL INITIAL DATA   ***
! ************************************

! This initial data corresponds to a gaussian shell of fluid
! initially at rest.  The lapse, shift and conformal metric
! and extrinsic curvature are left as in Minkowski, and
! the Hamiltonian constraint is solved for the conformal
! factor.
!
! Notice that the initial shell is assumed to be at rest,
! as otherwise we would have a non-zero momentum density
! and we would need to solve the coupled momentum and
! hamiltonian constraints.
!
! Notice also that here we assume that we have a polytropic
! equation of state:
!
! p = K*rho0^gamma
!
! which in turn implies that the specific internal energy is:
!
! e = K/(gamma-1)*rho0^(gamma-1)
!
! We will assume a conformally flat metric and
! solve the Hamiltonian constraint for the conformal
! factor, which in this case takes the form:
!
!  2                                5
! d psi  +  (2/r) d psi  +  2 pi psi rho  =  0
!  r               r
!
! with rho the energy density of the fluid at rest, which
! is given by:
!
! rho  =  rho0*h - p  =  rho0*(1 + e)
!
! IMPORTANT: Since the equation is non-linear, we solve it
! by integrating from the origin using a fourth order
! Runge-Kutta method. We later need to rescale the solution
! to have the correct asymptotic value of psi. This means
! in particular that the final value of rho(r=0) will no
! longer correspond to the value of the parameter fluid_a0.
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
  use radialfunctions

! Extra variables.

  implicit none

  integer i,l,imin
  integer j,maxiter

  real(8) r0,delta                      ! Local radius and grid spacing.
  real(8) psi0,Dpsi0                    ! Initial values of variables.
  real(8) k11,k12,k13,k14               ! Runge-Kutta sources for psi.
  real(8) k21,k22,k23,k24               ! Runge-Kutta sources for Dpsi.
  real(8) psi_rk,Dpsi_rk,rho_rk         ! Runge-Kutta values of variables.
  real(8) rm,aux,small,res
  real(8) one,half,smallpi

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr            ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: psi_g,Dpsi_g  ! Global solution and derivative.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rho_g         ! Global energy density.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: mass_g        ! Global rest-mass density.


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0
  one  = 1.d0

  smallpi = acos(-1.d0)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for fluid shell ...'
  end if

! Sanity check: gamma must be larger than 1,
! and kappa must be positive.

  if (fluid_gamma<=1.d0) then
     print *
     print *, "The adiabatic index gamma must be larger than 1."
     print *, 'Aborting! (subroutine idata_fluidshell)'
     print *
     call die
  end if

  if (fluid_kappa<0.d0) then
     print *
     print *, "The polytropic constant kappa must be non-negative."
     print *, 'Aborting! (subroutine idata_fluidshell)'
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


! *****************************************
! ***   SET UP SHELL IN FLUID DENSITY   ***
! *****************************************

! Set up shell in fluid density. Remember that the
! rest-mass density must be even.

! Gaussian profile.

  if (fluidprofile=="gaussian") then

     if (fluid_r0==0.d0) then
        fluid_rho = gaussian(fluid_a0,0.d0,fluid_s0)
     else
        fluid_rho = gaussian(fluid_a0,+fluid_r0,fluid_s0) &
                  + gaussian(fluid_a0,-fluid_r0,fluid_s0)
     end if

! Smooth top-hat profile.

  else if (fluidprofile=="tophat") then

     if (fluid_r0==0.d0) then
        fluid_rho = tophat(fluid_a0,0.d0,fluid_s0,fluid_t0)
     else
        fluid_rho = tophat(fluid_a0,+fluid_r0,fluid_s0,fluid_t0) &
                  + tophat(fluid_a0,-fluid_r0,fluid_s0,fluid_t0)
     end if

  end if

! Add contribution from internal energy.

  if (fluid_EOS=="ideal") then
     rho_g = fluid_rho + fluid_kappa/(fluid_gamma-1.d0)*fluid_rho**fluid_gamma
  else
     print *
     print *, 'At the moment fluidshell initial data is only compatible'
     print *, 'with an ideal gas (polytropic) equation of state'
     print *, 'Aborting! (subroutine idata_fluidshell)'
     print *
     call die
  end if


! ***********************
! ***   INTEGRATION   ***
! ***********************

! Initialize arrays.

  psi_g  = 1.d0
  Dpsi_g = 0.d0


! *********************************************
! ***   ONLY PROCESSOR 0 SOLVES THE ODE's   ***
! *********************************************

  if (rank==0) then


!    *********************************
!    ***   LOOP OVER GRID LEVELS   ***
!    *********************************

!    We solve from fine to coarse grid.

     do l=Nl-1,0,-1


!       *******************************
!       ***   FIND INITIAL VALUES   ***
!       *******************************

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
           psi_g(l,imin-1)  = (9.d0*(psi_g(l+1,Nrtotal-2)+psi_g(l+1,Nrtotal-3)) &
                            - (psi_g(l+1,Nrtotal-4)+psi_g(l+1,Nrtotal-1)))/16.d0
           Dpsi_g(l,imin-1) = (9.d0*(Dpsi_g(l+1,Nrtotal-2)+Dpsi_g(l+1,Nrtotal-3)) &
                            - (Dpsi_g(l+1,Nrtotal-4)+Dpsi_g(l+1,Nrtotal-1)))/16.d0
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

!             Values of (psi,Dpsi) at origin.

              psi0  = 1.d0
              Dpsi0 = 0.d0

!             Grid spacing and values at previous grid point.

           else

              delta = dr(l)
              r0    = rr(l,i-1)

              psi0  = psi_g(l,i-1)
              Dpsi0 = Dpsi_g(l,i-1)

           end if

!          I) First Runge-Kutta step.

!          Sources at first grid point if we start
!          from the origin (for finer grid).

           if (i==1) then

!             At the origin we have:  psi' = 0.

              k11 = 0.d0

!             For the derivative of Dpsi at the first grid point we
!             use the following aproximation:  we assume that close
!             to the origin  Dpsi = a*r, and substitute this in the
!             Hamiltonian constraint.  Doing this we find that:
!
!             a  = - (2 pi / 3) rho(r=0)

              if (order=="two") then
                 rho_rk = 0.5d0*(rho_g(l,0) + rho_g(l,1))
              else
                 rho_rk = (9.d0*(rho_g(l,0) + rho_g(l,1)) - (rho_g(l,-1) + rho_g(l,2)))/16.d0
              end if

              k21 = - delta*(2.d0*smallpi/3.d0)*rho_rk

!          Sources at previous grid point.

           else

              rm = r0

              rho_rk  = rho_g(l,i-1)

              psi_rk  = psi0
              Dpsi_rk = Dpsi0

              k11 = delta*Dpsi_rk
              k21 = - 2.d0*delta*(Dpsi_rk/rm + smallpi*psi_rk**5*rho_rk)

           end if

!          II) Second Runge-Kutta step.

           rm = r0 + half*delta

           rho_rk = (9.d0*(rho_g(l,i-1) + rho_g(l,i)) - (rho_g(l,i-2) + rho_g(l,i+1)))/16.d0

           psi_rk  = psi0  + half*k11
           Dpsi_rk = Dpsi0 + half*k21

           k12 = delta*Dpsi_rk
           k22 = - 2.d0*delta*(Dpsi_rk/rm + smallpi*psi_rk**5*rho_rk)

!          III) Third Runge-Kutta step.

           psi_rk  = psi0  + half*k12
           Dpsi_rk = Dpsi0 + half*k22

           k13 = delta*Dpsi_rk
           k23 = - 2.d0*delta*(Dpsi_rk/rm + smallpi*psi_rk**5*rho_rk)

!          IV) Fourth Runge-Kutta step.

           rm = r0 + delta

           rho_rk  = rho_g(l,i)

           psi_rk  = psi0  + k13
           Dpsi_rk = Dpsi0 + k23

           k14 = delta*Dpsi_rk
           k24 = - 2.d0*delta*(Dpsi_rk/rm + smallpi*psi_rk**5*rho_rk)

!          Advance variables to next grid point.

           psi_g(l,i)  = psi0   + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
           Dpsi_g(l,i) = Dpsi0  + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

        end do


!       ***********************
!       ***   GHOST ZONES   ***
!       ***********************

        do i=1,ghost
           psi_g(l,1-i)  = + psi_g(l,i)
           Dpsi_g(l,1-i) = - Dpsi_g(l,i)
        end do

     end do


!    ************************************
!    ***   RESCALE CONFORMAL FACTOR   ***
!    ************************************

!    Figure out the asymptotic value of the conformal factor and use it
!    to rescale the central value so that the new asymptotic value is 1.
!    To find the asymptotic value I assume we have a fall-off of the form:
!
!    u = u0 + k/r   =>   du/dr = (u0-u)/r
!
!    Solving for u0 we find:
!
!    u0 = u + r (du/dr)
!
!    I then approximate the derivative at the boundary using one-sided
!    differences.

     aux = psi_g(0,Nrtotal) + rr(0,Nrtotal)*Dpsi_g(0,Nrtotal)
     psi_g = psi_g/aux

!    Notice that since the Hamiltonian constraint is not linear
!    in psi we can't just rescale psi without also rescaling
!    the energy density rho.  In fact, if we change psi -> psi/aux
!    must take rho -> rho*aux**4 to keep the Hamiltonian constraint
!    satisfied.

     rho_g = rho_g*aux**4


!    ***************************************
!    ***   RECALCULATE FLUID VARIABLES   ***
!    ***************************************

!    Since we have rescaled rho above, we must now rescale also
!    the fluid density fluid_rho0.  However, in order to do
!    this for an ideal equation of state we must invert a
!    trascendental equation:
!
!    rho = rho0 + kappa/(gamma-1)*rho0**gamma
!
!    with rho the energy density and rho0 the rest-mass density.
!    Notice that here we use rho=rho_g and rho0=mass_g.
!
!    To invert the equation we use Newton-Rapson's method. 
!    We define a function:
!
!    F(rho0) := rho0 + kappa/(gamma-1)*rho0**gamma - rho
!
!    whose derivative is:
!
!    F'(rho0) = 1 + gamma*kappa/(gamma-1)*rho0**(gamma-1)
!
!    The method then takes an initial guess rho0=rho, and
!    iterates the following equation to some tolerance:
!
!    rho0_new = rho0_old - F(rho0_old)/F'(rho0_old)

     if (fluid_EOS=="ideal") then

!       Initialize.

        small = 1.d-10
        maxiter = 50

        aux = fluid_kappa/(fluid_gamma-1.d0)

!       Loop over grid levels.

        do l=Nl-1,0,-1

!          Loop over grid points.

           do i=1,Nrtotal

!             Initialize counter and residue.

              j = 0
              res = 1.d0

!             Trial value for fluid_rho0.

              mass_g(l,i) = rho_g(l,i)

!             Start iterations for Newton-Rapson's method.

              do while ((abs(res)>small).and.(j<maxiter))

!                Increment counter.

                 j = j+1

!                Calculate residual.

                 res = (mass_g(l,i) + aux*mass_g(l,i)**fluid_gamma - rho_g(l,i)) &
                     /(1.d0 + fluid_gamma*aux*mass_g(l,i)**(fluid_gamma-1.d0))

!                Update value.

                 mass_g(l,i) = mass_g(l,i) - res

              end do

           end do

!          If we exceeded the maximum number of iterations send message to screen.

           if (j==maxiter) then
              print *
              print *, 'Maximum iteration number reached in idata_fluidshell.f90 at grid point: ',i
              print *
           end if

!          Ghost points.

           do i=1,ghost
              mass_g(l,1-i) = mass_g(l,i)
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
           psi_g(l-1,i/2+1)  = (9.d0*(psi_g(l,i)+psi_g(l,i+1))   - (psi_g(l,i-1)+psi_g(l,i+2)))/16.d0
           Dpsi_g(l-1,i/2+1) = (9.d0*(Dpsi_g(l,i)+Dpsi_g(l,i+1)) - (Dpsi_g(l,i-1)+Dpsi_g(l,i+2)))/16.d0
           mass_g(l-1,i/2+1) = (9.d0*(mass_g(l,i)+mass_g(l,i+1)) - (mass_g(l,i-1)+mass_g(l,i+2)))/16.d0
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
! array with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  if (size==1) then
     psi = psi_g
     fluid_rho = mass_g
  else
     call distribute(0,Nl-1,psi,psi_g)
     call distribute(0,Nl-1,fluid_rho,mass_g)
  end if


! ******************************************
! ***   FIND ALL OTHER FLUID VARIABLES   ***
! ******************************************

! Find pressure and internal energy.

  if (fluid_EOS=="ideal") then
     fluid_p = fluid_kappa*fluid_rho**fluid_gamma
     fluid_e = fluid_kappa/(fluid_gamma-1.d0)*fluid_rho**(fluid_gamma-1.d0)
  end if

! Find enthalpy: h = 1 + e + p/rho0

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


! ***************************************
! ***   FIND CONFORMAL FUNCTION phi   ***
! ***************************************

! Find phi and chi.

  phi = dlog(psi)

  if (chimethod) then
     chi = one/psi**chipower
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine idata_fluidshell

