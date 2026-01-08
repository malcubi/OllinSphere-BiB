
  subroutine idata_BosonstarCF

! *************************************************************
! ***   BOSON STAR INITIAL DATA IN CONFORMALLY FLAT GAUGE   ***
! *************************************************************

! Boson stars are solutions such that the spacetime is static
! and the complex scalar field has a harmonic dependence on time.
! This subroutine calculates initial data for a boson star
! in the conformally flat gauge using a shooting method.

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
! The standard ansatz for the metric in this case is:
!
!   2          2   2       4    2     2      2
! ds  = - alpha  dt  +  psi ( dr  +  r dOmega )
!
! that is, we are asuming a conformally flat spatial metric.
!
! We substitute the above metric and this ansatz into the
! Hamiltonian constraint and the Klein-Gordon equation to
! find (xi=dphi/dr):
!
!  2                               5
! d psi  = - (2/r) d psi - 2 pi psi rho
!  r                r
!
!  2               /                                   \
! d phi  = - d phi | 2/r + d alpha/alpha + 2 d psi/psi |
!  r          r    \        r                 r        /
!
!              4 /                       2 \
!        +  psi  | V' - phi (omega/alpha)  |
!                \                         /
!
! Where V(phi) is the self-interaction potential, V'=dV/dphi,
! and where the energy density rho is given by:
!
!               /        2     4                         2       \
! rho  =  (1/2) | (d phi) / psi  +  ( omega phi / alpha )  + 2 V |
!               \   r                                            /
!
! For the lapse we use maximal slicing, which implies:
!
!  2                     /                 \                   4
! d alpha = - 2 d alpha  | 1/r + d psi/psi |  +  4 pi alpha psi (rho + S)
!  r             r       \        r        /
!
! where now:
!                                      2
! rho + S  =  2 [ ( omega phi / alpha )  -  V ]
!
! This routine takes as input parameter the value of the scalar
! field at the origin "boson_phi0", and uses a shooting method
! (with Runge-Kutta) to determine the value of the frequency
! "boson_omega" that is compatible with the asymptotically
! decaying solution.
!
! For the boundary conditions at the origin we take:
!
! psi(r=0)   = 1,     d psi(r=0)   = 0
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
! just rescale the lapse at the end. But in order not to affect
! the solution for A and phi we must also rescale the final value
! of the frequency omega by the same factor.
!
! We also don't want psi(r=0)=1, but rather psi=1 at infinity.
! Unfortunately, the equation for psi is non-linear in the
! source term, so here we can't just rescale the solution.
! However, once we have a solution, we can rescale it as:
!
! r -> k r ,       psi ->  psi/sqrt(k)
!
! and we still have a solution (we can check that this is
! true for all equations above, since the density "rho"
! for a scalar field is invariant with this change).
!
! In principle we could just do this change of radial
! coordinate and interpolate, but we would loose precision.
! What I do instead is figure out the asymptotic value of psi
! in the solution with psi(r=0)=1, and use it to find the factor
! needed to rescale psi so that it goes to 1 at infinity.  I then
! rescale the central value of psi with this same factor and solve
! again the whole system. Doing this we should recover the same
! frequency (up to truncation error).
!
! PERTURBATIONS:  For the moment this rountine does not solve
! for perturbed initial data.
!
! ANGULAR MOMENTUM:  This routine does not yet allow for the
! case of non-zero angular momentum (complex_l must be zero).


! NOTE FOR PARALLEL RUNS:  The initial data is not really
! solved in parallel.  It is in fact solved only on processor
! zero on a full size array, and then it is distributed
! among the other processors.  This is slow, but works.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  logical :: left  = .false.            ! This flag is true if the left value of omega is well behaved.
  logical :: right = .false.            ! This flag is true if the right value of omega is well behaved.

  integer i,l,iter                      ! Counters.
  integer imin                          ! Leftmost grid point.
  integer :: maxiter = 1000             ! Maximum number of iterations.

  real(8) r0,delta                      ! Local radius and grid spacing.
  real(8) psi0,Dpsi0                    ! Initial values of variables.
  real(8) alpha0,Dalpha0                ! Initial values of variables.
  real(8) phi0,xi0                      ! Initial values of variables.
  real(8) V0,DV0                        ! Potential and derivative at origin.
  real(8) k11,k12,k13,k14               ! Runge-Kutta sources for psi.
  real(8) k21,k22,k23,k24               ! Runge-Kutta sources for Dpsi.
  real(8) k31,k32,k33,k34               ! Runge-Kutta sources for alpha.
  real(8) k41,k42,k43,k44               ! Runge-Kutta sources for Dalpha.
  real(8) k51,k52,k53,k54               ! Runge-Kutta sources for phi.
  real(8) k61,k62,k63,k64               ! Runge-Kutta sources for xi.
  real(8) psi_rk,Dpsi_rk                ! Runge-Kutta values of variables.
  real(8) alpha_rk,Dalpha_rk            ! Runge-Kutta values of variables.
  real(8) phi_rk,xi_rk                  ! Runge-Kutta values of variables.
  real(8) J1_CF,J2_CF,J3_CF             ! Functions for sources of differential equations.
  real(8) J4_CF,J5_CF,J6_CF             ! Functions for sources of differential equations.
  real(8) res,res_old                   ! Residual.
  real(8) omega_new,omega_old,domega    ! Trial frequency and frequency interval.
  real(8) one,half,smallpi              ! Numbers.
  real(8) rm,aux,psiorigin              ! Auxiliary quantities.
  real(8) :: epsilon = 1.d-8            ! Tolerance.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                 ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: psi_g,Dpsi_g       ! Conformal factor and derivative global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: alpha_g,Dalpha_g   ! Lapse and derivative global arrays.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: phi_g,xi_g         ! Scalar field global arrays.


! *******************
! ***   NUMBERS   ***
! *******************

  one  = 1.d0
  half = 0.5d0

  smallpi = acos(-1.d0)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for a boson star in conformally flat gauge using shooting method ...'
     print *
  end if


! *************************
! ***   NORMALIZATION   ***
! *************************

! Normalization.  For the rotating case the phi0 is rescaled by
! sqrt(2*l+1).  An additional factor of 1/sqrt(4*pi) can also be
! included (this is Olvier Sarbach's convention).

  if (boson_factor=="physical") then

     boson_phi0 = boson_phi0*sqrt(dble(2*complex_l+1))

     if (rank==0) then
        if (complex_l==0) then
           write(*,'(A,E20.14)') ' Using physical normalization for phi0 ...'
           print *
        else
           write(*,'(A,E20.14)') ' Rescaled phi0 (phi0*sqrt(2l+1)): ', boson_phi0
           print *
        end if
     end if

  else

     boson_phi0 = boson_phi0*sqrt(dble(2*complex_l+1)/(4.d0*smallpi))

     if (rank==0) then
        if (complex_l==0) then
           write(*,'(A,E20.14)') ' Physical phi0 (phi0/sqrt(4pi)): ', boson_phi0
           print *
        else
           write(*,'(A,E20.14)') ' Rescaled phi0 (phi0*sqrt((2l+1)/4pi)): ', boson_phi0
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
     print *, 'Aborting! (subroutine idata_BosonstarCF)'
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
! integration is done with either second or
! fourth order Runge-Kutta.

! For the moment only processor 0 solves for the
! initial data.  The solution is later distributed
! among processors.

! Initialize value of psi at origin.

  psiorigin = 1.d0

! Initialize arrays.

  psi_g  = 1.d0
  Dpsi_g = 0.d0

  alpha_g  = 1.d0
  Dalpha_g = 0.d0

  phi_g = 0.d0
  xi_g  = 0.d0


! *********************************************
! ***   ONLY PROCESSOR 0 SOLVES THE ODE's   ***
! *********************************************

  if (rank==0) then


!    **********************
!    ***   INITIALIZE   ***
!    **********************

!    Find initial domega.

     domega = (omega_right - omega_left)/20.d0

!    Sanity check.

     if (domega<0.d0) then
        print *
        print *, 'For boson star initial data we must have omega_right>omega_left.'
        print *, 'Aborting! (subroutine idata_BosonstarCF)'
        print *
        call die
     end if

!    If omega_right=omega_left it means we want a fixed frequency with NO iterations.

     if (domega==0.d0) then
        maxiter = 1
     end if


!    ****************************
!    ***   BEGIN ITERATIONS   ***
!    ****************************

!    This continue statement is here for when we solve again with
!    a corrected value of the conformal factor at the origin.

     100  continue

!    Initialize residual.

     res = 1.d0

!    Start iterations.

     boson_omega = omega_left

     iter = 0

     do while ((abs(res).gt.epsilon).and.(iter.lt.maxiter))

        iter = iter + 1


!       *********************************
!       ***   LOOP OVER GRID LEVELS   ***
!       *********************************

!       We solve from fine to coarse grid.

        do l=Nl-1,0,-1


!          ****************************************
!          ***   FIND INITIAL VALUES AT ORIGIN  ***
!          ****************************************

!          Potential and its derivative at the origin.

           V0  = 0.5d0*complex_mass**2*boson_phi0**2 + 0.25d0*complex_lambda*boson_phi0**4
           DV0 = complex_mass**2*boson_phi0 + complex_lambda*boson_phi0**3

!          Find initial point.  Only the finest grid
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

              psi_g(l,imin-1)  = (9.d0*(psi_g(l+1,Nrtotal-2)+psi_g(l+1,Nrtotal-3)) &
                               - (psi_g(l+1,Nrtotal-4)+psi_g(l+1,Nrtotal-1)))/16.d0
              Dpsi_g(l,imin-1) = (9.d0*(Dpsi_g(l+1,Nrtotal-2)+Dpsi_g(l+1,Nrtotal-3)) &
                               - (Dpsi_g(l+1,Nrtotal-4)+Dpsi_g(l+1,Nrtotal-1)))/16.d0

              alpha_g(l,imin-1)  = (9.d0*(alpha_g(l+1,Nrtotal-2)+alpha_g(l+1,Nrtotal-3)) &
                                 - (alpha_g(l+1,Nrtotal-4)+alpha_g(l+1,Nrtotal-1)))/16.d0
              Dalpha_g(l,imin-1) = (9.d0*(Dalpha_g(l+1,Nrtotal-2)+Dalpha_g(l+1,Nrtotal-3)) &
                                 - (Dalpha_g(l+1,Nrtotal-4)+Dalpha_g(l+1,Nrtotal-1)))/16.d0

              phi_g(l,imin-1) = (9.d0*(phi_g(l+1,Nrtotal-2)+phi_g(l+1,Nrtotal-3)) &
                              - (phi_g(l+1,Nrtotal-4)+phi_g(l+1,Nrtotal-1)))/16.d0
              xi_g(l,imin-1)  = (9.d0*(xi_g(l+1,Nrtotal-2)+xi_g(l+1,Nrtotal-3)) &
                              - (xi_g(l+1,Nrtotal-4)+xi_g(l+1,Nrtotal-1)))/16.d0

           end if


!          ************************************
!          ***   SECOND ORDER RUNGE-KUTTA   ***
!          ************************************

           if (order=='two') then

              do i=imin,Nrtotal

!                Grid spacing and values at first point
!                if we start from the origin (finer grid).

                 if (i==1) then

!                   For the first point we use dr/2.

                    delta = half*dr(l)
                    r0    = 0.d0

!                   Values of (psi,alpha) and derivatives at r=0.

                    psi0  = psiorigin
                    Dpsi0 = 0.d0

                    alpha0  = 1.d0
                    Dalpha0 = 0.d0

!                   For zero angular momentum (complex_l=0) the
!                   scalar field goes as phi ~ phi0 + O(r^2)
!                   at the origin.

                    if (complex_l==0) then

                       phi0 = boson_phi0
                       xi0  = 0.d0

!                   For complex_l=1, the scalar field goes as
!                   phi ~ phi0*r + 0(r^3) at the origin.

                    else if (complex_l==1) then

                       phi0 = 0.d0
                       xi0  = boson_phi0

!                   For non-zero angular momentum the scalar field
!                   goes as phi ~ phi0*r^l + O(r^(l+2)) at the origin.

                    else

                       phi0 = 0.d0
                       xi0  = 0.d0
 
                    end if

!                Grid spacing and current values if we are not at the origin.

                 else

                    delta  = dr(l)
                    r0     = rr(l,i-1)

                    psi0  = psi_g(l,i-1)
                    Dpsi0 = Dpsi_g(l,i-1)

                    alpha0  = alpha_g(l,i-1)
                    Dalpha0 = Dalpha_g(l,i-1)

                    phi0 = phi_g(l,i-1)
                    xi0  = xi_g(l,i-1)

                 end if

!                First Runge-Kutta step.

!                Sources at first grid point if we start
!                from the origin (for finer grid).

                 if (i==1) then

!                   At the origin we have: Dpsi=Dalpha=xi=0.
          
                    k11 = 0.d0
                    k31 = 0.d0
                    k51 = 0.d0

!                   For zero angular momentum we also have:
!
!                    2                      5                       2
!                   d psi   ~ - (pi/3)  psi0 ( (omega phi0 / alpha0) + 2 V0 )
!                    r
!
!                    2                      4                       2
!                   d alpha ~ + (8pi/3) psi0 ( (omega phi0 / alpha0)  - 2 V0)
!                    r
!
!                                       4                            2
!                   d xi    ~ (1/3) psi0 ( V' - phi0 (omega / alpha0) )
!                    r

                    if (complex_l==0) then

                       k21 = - smallpi/3.d0*psi0**5*((boson_omega*boson_phi0/alpha0)**2 + 2*V0)*delta
                       k41 = 8.d0*smallpi/3.d0*psi0**4*((boson_omega*boson_phi0/alpha0)**2 - 2*V0)*delta
                       k61 = psi0**4/3.d0*(DV0 - boson_phi0*(boson_omega/alpha0)**2)*delta

                    else

                       print *
                       print *, 'Non-zero values of complex_l not yet implemented.'
                       print *, 'Aborting! (subroutine idata_BosonstarCF)'
                       print *
                       call die

                    end if

!                Sources if we are not at the origin.

                 else

                    rm = r0

                    psi_rk  = psi0
                    Dpsi_rk = Dpsi0

                    alpha_rk = alpha0
                    Dalpha_rk = Dalpha0

                    phi_rk = phi0
                    xi_rk  = xi0

                    k11 = delta*J1_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                    k21 = delta*J2_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                    k31 = delta*J3_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                    k41 = delta*J4_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                    k51 = delta*J5_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                    k61 = delta*J6_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

                 end if

!                Second Runge-Kutta step.

                 rm = r0 + half*delta

                 psi_rk  = psi0  + half*k11
                 Dpsi_rk = Dpsi0 + half*k21

                 alpha_rk = alpha0   + half*k31
                 Dalpha_rk = Dalpha0 + half*k41

                 phi_rk = phi0 + half*k51
                 xi_rk  = xi0  + half*k61

                 k12 = delta*J1_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k22 = delta*J2_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k32 = delta*J3_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k42 = delta*J4_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k52 = delta*J5_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k62 = delta*J6_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

!                Advance variables.

                 psi_g(l,i)  = psi0  + k12
                 Dpsi_g(l,i) = Dpsi0 + k22

                 alpha_g(l,i)  = alpha0  + k32
                 Dalpha_g(l,i) = Dalpha0 + k42

                 phi_g(l,i) = phi0 + k52
                 xi_g(l,i)  = xi0  + k62

!                Check if solution is blowing up. This helps to reduce
!                the need for a very fine tuned initial guess.

                 if (abs(phi_g(l,i))>boson_phi0) then
                    if (.not.left) then
                       omega_left = omega_left + domega
                       boson_omega = omega_left
                       goto 200
                    else if (.not.right) then
                       omega_right = omega_right - domega
                       boson_omega = omega_right
                       goto 200
                    end if
                 end if

              end do


!          ************************************
!          ***   FOURTH ORDER RUNGE-KUTTA   ***
!          ************************************

           else

              do i=imin,Nrtotal

!                Grid spacing and values at first point
!                if we start from the origin (finer grid).

                 if (i==1) then

!                   For the first point we use dr/2.

                    delta = half*dr(l)
                    r0    = 0.d0

!                   Values of psi alpha and derivatives at r=0.

                    psi0  = psiorigin
                    Dpsi0 = 0.d0

                    alpha0  = 1.d0
                    Dalpha0 = 0.d0

!                   For zero angular momentum (complex_l=0) the
!                   scalar field goes as phi ~ phi0 + O(r^2)
!                   at the origin.

                    if (complex_l==0) then

                       phi0 = boson_phi0
                       xi0  = 0.d0

!                   For complex_l=1, the scalar field goes as
!                   phi ~ phi0*r + 0(r^3) at the origin.

                    else if (complex_l==1) then

                       phi0 = 0.d0
                       xi0  = boson_phi0

!                   For non-zero angular momentum the scalar field
!                   goes as phi ~ phi0*r^l + O(r^(l+2)) at the origin.

                    else

                       phi0 = 0.d0
                       xi0  = 0.d0
 
                    end if

!                Grid spacing and current values if we are not at the origin.

                 else

                    delta  = dr(l)
                    r0     = rr(l,i-1)

                    psi0  = psi_g(l,i-1)
                    Dpsi0 = Dpsi_g(l,i-1)

                    alpha0  = alpha_g(l,i-1)
                    Dalpha0 = Dalpha_g(l,i-1)

                    phi0 = phi_g(l,i-1)
                    xi0  = xi_g(l,i-1)

                 end if

!                First Runge-Kutta step.

!                Sources at first grid point if we start
!                from the origin (for finer grid).

                 if (i==1) then

!                   At the origin we have: Dpsi=Dalpha=xi=0.
          
                    k11 = 0.d0
                    k31 = 0.d0
                    k51 = 0.d0

!                   We also have:
!
!                    2                      5                       2
!                   d psi   ~ - (pi/3)  psi0 ( (omega phi0 / alpha0) + 2 V0 )
!                    r
!
!                    2                      4                       2
!                   d alpha ~ + (8pi/3) psi0 ( (omega phi0 / alpha0)  - 2 V0)
!                    r
!
!                                       4                            2
!                   d xi    ~ (1/3) psi0 ( V' - phi0 (omega / alpha0) )
!                    r

                    if (complex_l==0) then

                       k21 = - smallpi/3.d0*psi0**5*((boson_omega*boson_phi0/alpha0)**2 + 2*V0)*delta
                       k41 = 8.d0*smallpi/3.d0*psi0**4*((boson_omega*boson_phi0/alpha0)**2 - 2*V0)*delta
                       k61 = psi0**4/3.d0*(DV0 - boson_phi0*(boson_omega/alpha0)**2)*delta

                    else

                       print *
                       print *, 'Non-zero values of complex_l not yet implemented.'
                       print *, 'Aborting! (subroutine idata_BosonstarCF)'
                       print *
                       call die

                    end if

!                Sources if we are not at the origin.

                 else

                    rm = r0

                    psi_rk  = psi0
                    Dpsi_rk = Dpsi0

                    alpha_rk = alpha0
                    Dalpha_rk = Dalpha0

                    phi_rk = phi0
                    xi_rk  = xi0

                    k11 = delta*J1_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                    k21 = delta*J2_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                    k31 = delta*J3_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                    k41 = delta*J4_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                    k51 = delta*J5_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                    k61 = delta*J6_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

                 end if

!                Second Runge-Kutta step.

                 rm = r0 + half*delta

                 psi_rk  = psi0  + half*k11
                 Dpsi_rk = Dpsi0 + half*k21

                 alpha_rk = alpha0   + half*k31
                 Dalpha_rk = Dalpha0 + half*k41

                 phi_rk = phi0 + half*k51
                 xi_rk  = xi0  + half*k61

                 k12 = delta*J1_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k22 = delta*J2_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k32 = delta*J3_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k42 = delta*J4_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k52 = delta*J5_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k62 = delta*J6_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

!                Third Runge-Kutta step.

                 psi_rk  = psi0  + half*k12
                 Dpsi_rk = Dpsi0 + half*k22

                 alpha_rk = alpha0   + half*k32
                 Dalpha_rk = Dalpha0 + half*k42

                 phi_rk = phi0 + half*k52
                 xi_rk  = xi0  + half*k62

                 k13 = delta*J1_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k23 = delta*J2_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k33 = delta*J3_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k43 = delta*J4_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k53 = delta*J5_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k63 = delta*J6_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

!                Fourth Runge-Kutta step.

                 rm = r0 + delta

                 psi_rk  = psi0  + k13
                 Dpsi_rk = Dpsi0 + k23

                 alpha_rk = alpha0   + k33
                 Dalpha_rk = Dalpha0 + k43

                 phi_rk = phi0 + k53
                 xi_rk  = xi0  + k63

                 k14 = delta*J1_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k24 = delta*J2_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k34 = delta*J3_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k44 = delta*J4_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k54 = delta*J5_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)
                 k64 = delta*J6_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

!                Advance variables.

                 psi_g(l,i)    = psi0    + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
                 Dpsi_g(l,i)   = Dpsi0   + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

                 alpha_g(l,i)  = alpha0  + (k31 + 2.d0*(k32 + k33) + k34)/6.d0
                 Dalpha_g(l,i) = Dalpha0 + (k41 + 2.d0*(k42 + k43) + k44)/6.d0

                 phi_g(l,i)    = phi0    + (k51 + 2.d0*(k52 + k53) + k54)/6.d0
                 xi_g(l,i)     = xi0     + (k61 + 2.d0*(k62 + k63) + k64)/6.d0

!                Check if solution is blowing up. This helps to reduce
!                the need for a very fine tuned initial guess.

                 if (abs(phi_g(l,i))>boson_phi0) then
                    if (.not.left) then
                       omega_left = omega_left + domega
                       boson_omega = omega_left
                       goto 200
                    else if (.not.right) then
                       omega_right = omega_right - domega
                       boson_omega = omega_right
                       goto 200
                    end if
                 end if

              end do

           end if

        end do


!       ***********************
!       ***   GHOST ZONES   ***
!       ***********************

!       Ghost zones only for fine grid.

        do i=1,ghost

           psi_g(Nl-1,1-i)  = + psi_g(Nl-1,i)
           Dpsi_g(Nl-1,1-i) = - Dpsi_g(Nl-1,i)

           alpha_g(Nl-1,1-i)  = + alpha_g(Nl-1,i)
           Dalpha_g(Nl-1,1-i) = - Dalpha_g(Nl-1,i)

           phi_g(Nl-1,1-i) = + phi_g(Nl-1,i)
           xi_g(Nl-1,1-i)  = - xi_g(Nl-1,i)

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

        if (abs(phi_g(0,Nrtotal))+abs(phi_g(0,Nrtotal-1))<epsilon) then
           res = epsilon/2.d0
           goto 200
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
           omega_new = boson_omega - res*(boson_omega - omega_old)/aux

           omega_old = boson_omega
           boson_omega = omega_new

        end if

200     continue


!       **************************
!       ***   END ITERATIONS   ***
!       **************************

!       Output data to screen.

        write(*,"(A,I4,A,ES20.14,A,ES9.2)") ' Iteration: ',iter,'    Frequency: ',boson_omega,'    Residual: ',res

     end do

!    Message in case we reached the maximum iteration number.

     if (iter>=maxiter) then
        print *
        print *, 'Maximum iteration number reached in idata_bosonstar.f90.'
        print *, 'Initial data solver did not converge, aborting!'
        print *
        call die
     end if


!    ************************************
!    ***   RESCALE CONFORMAL FACTOR   ***
!    ************************************

!    Figure out the asymptotic value of the conformal factor and use it
!    to rescale the central value so that the new asymptotic value is 1.
!    We then solve the whole system again.  This is not terribly
!    efficient, but seems to work fine.
!
!    To find the asymptotic value I assume we have a fall off of the form:
!
!    u = u0 + k/r   =>   du/dr = (u0-u)/r
!
!    Solving for u0 we find:
!
!    u0 = u + r (du/dr)
!
!    I then approximate the derivative at the boundary using one-sided
!    differences.

     if (order=="two") then
        aux = psi_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
            *(3.d0*psi_g(0,Nrtotal) - 4.d0*psi_g(0,Nrtotal-1) + psi_g(0,Nrtotal-2))
     else
        aux = psi_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
            *(25.d0*psi_g(0,Nrtotal) - 48.d0*psi_g(0,Nrtotal-1) &
            + 36.d0*psi_g(0,Nrtotal-2) - 16.d0*psi_g(0,Nrtotal-3) + 3.d0*psi_g(0,Nrtotal-4))/3.d0
     end if

     print *
     print *, 'Asymptotic value of conformal factor: ',aux
     print *, 'Omega: ',boson_omega

     if (psiorigin==1.d0) then
        print *
        print *, 'Restarting iterations with rescaled conformal factor at origin ...'
        psiorigin = 1.d0/aux
        goto 100
     end if

!    Print success message.

     print *
     print *, 'Done!'
     print *


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

           alpha_g(l-1,i/2+1)  = (9.d0*(alpha_g(l,i)+alpha_g(l,i+1))   - (alpha_g(l,i-1)+alpha_g(l,i+2)))/16.d0
           Dalpha_g(l-1,i/2+1) = (9.d0*(Dalpha_g(l,i)+Dalpha_g(l,i+1)) - (Dalpha_g(l,i-1)+Dalpha_g(l,i+2)))/16.d0

           phi_g(l-1,i/2+1) = (9.d0*(phi_g(l,i)+phi_g(l,i+1)) - (phi_g(l,i-1)+phi_g(l,i+2)))/16.d0
           xi_g(l-1,i/2+1)  = (9.d0*(xi_g(l,i)+xi_g(l,i+1))   - (xi_g(l,i-1)+xi_g(l,i+2)))/16.d0

        end do

!       Fix ghost zones.

        do i=1,ghost

           psi_g(l-1,1-i)  = + psi_g(l-1,i)
           Dpsi_g(l-1,1-i) = - Dpsi_g(l-1,i)

           alpha_g(l-1,1-i)  = + alpha_g(l-1,i)
           Dalpha_g(l-1,1-i) = - Dalpha_g(l-1,i)

           phi_g(l-1,1-i) = + phi_g(l-1,i)
           xi_g(l-1,1-i)  = - xi_g(l-1,i)

        end do

     end do


!    *************************
!    ***   RESCALE LAPSE   ***
!    *************************

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
        aux = alpha_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
            *(3.d0*alpha_g(0,Nrtotal) - 4.d0*alpha_g(0,Nrtotal-1) + alpha_g(0,Nrtotal-2))
     else
        aux = alpha_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
            *(25.d0*alpha_g(0,Nrtotal) - 48.d0*alpha_g(0,Nrtotal-1) &
            + 36.d0*alpha_g(0,Nrtotal-2) - 16.d0*alpha_g(0,Nrtotal-3) + 3.d0*alpha_g(0,Nrtotal-4))/3.d0
     end if

     alpha_g = alpha_g/aux
     omega_new = boson_omega/aux

!    Write value of omega to screen.

     if (rank==0) then
        write(*,'(A,E22.16)') ' Omega (not-rescaled) = ', boson_omega
        write(*,'(A,E22.16)') ' Omega (rescaled)     = ', omega_new
        print *
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

! Broadcast the new value of the frequency to all processors.

  call MPI_BCAST(boson_omega,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! arrays with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  if (size==1) then
     psi    = psi_g
     D1_psi = Dpsi_g
     alpha    = alpha_g
     D1_alpha = Dalpha_g
     complex_phiR = phi_g
     complex_xiR  = xi_g
  else
     call distribute(0,Nl-1,psi,psi_g)
     call distribute(0,Nl-1,D1_psi,Dpsi_g)
     call distribute(0,Nl-1,alpha,alpha_g)
     call distribute(0,Nl-1,D1_alpha,Dalpha_g)
     call distribute(0,Nl-1,complex_phiR,phi_g)
     call distribute(0,Nl-1,complex_xiR,xi_g)
  end if


! **********************************************************************
! ***   IMAGINARY PART AND TIME DERIVATIVE OF COMPLEX SCALAR FIELD   ***
! **********************************************************************

! At t=0 we set the imaginary part of the complex scalar field
! and its spatial derivatve to 0.

  complex_phiI = 0.d0
  complex_xiI  = 0.d0
 
! Set time derivative of complex field.  From the ansatz we know
! that the time derivative is purely imaginary.

  complex_piR = 0.d0
  complex_piI = omega_new*complex_phiR/alpha


! ***************************************
! ***   FIND CONFORMAL FUNCTION phi   ***
! ***************************************

! Find phi and derivative.

  phi  = dlog(psi)
  D1_phi = D1_psi/psi

! Find chi and derviative.

  if (chimethod) then
     chi  = one/psi**chipower
     D1_chi = - dble(chipower)*D1_psi/psi**3
  end if

! Find (psi2,psi4).

  psi2 = psi**2
  psi4 = psi**4


! ***************
! ***   END   ***
! ***************

  end subroutine idata_BosonstarCF







! ************************************
! ***   RADIAL DERIVATIVE OF PSI   ***
! ************************************

! This is just the definition:  Dpsi=dpsi/dr

  function J1_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

  implicit none

  real(8) J1_CF,rm
  real(8) psi_rk,Dpsi_rk
  real(8) alpha_rk,Dalpha_rk
  real(8) xi_rk,phi_rk

! dpsi/dr = Dpsi.

  J1_CF = Dpsi_rk

  end function J1_CF







! *************************************
! ***   RADIAL DERIVATIVE OF DPSI   ***
! *************************************

  function J2_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

! This comes from the Hamiltonian constraint.

  use param

  implicit none

  real(8) J2_CF,rm
  real(8) psi_rk,Dpsi_rk
  real(8) alpha_rk,Dalpha_rk
  real(8) xi_rk,phi_rk
  real(8) smallpi,V

! Numbers.

  smallpi = acos(-1.d0)

! Potential.

  V = 0.5d0*complex_mass**2*phi_rk**2 + 0.25d0*complex_lambda*phi_rk**4

! Angular momentum.

  if (complex_l>0) then
     V = V + 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk**2
  end if

!                            /    5                            2           2 \
! d Dpsi = - (2/r) Dpsi - pi | psi ( 2 V + (omega phi / alpha )  ) + psi xi  |
!  r                         \                                               /

  J2_CF = - 2.d0/rm*Dpsi_rk - smallpi*(psi_rk**5*(2.d0*V &
        + (boson_omega*phi_rk/alpha_rk)**2) + psi_rk*xi_rk**2)

  end function J2_CF








! **************************************
! ***   RADIAL DERIVATIVE OF ALPHA   ***
! **************************************

! This is just the definition:  Dalpha=dalpha/dr

  function J3_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

  implicit none

  real(8) J3_CF,rm
  real(8) psi_rk,Dpsi_rk
  real(8) alpha_rk,Dalpha_rk
  real(8) xi_rk,phi_rk

! dalpha/dr = Dalpha.

  J3_CF = Dalpha_rk

  end function J3_CF







! ***************************************
! ***   RADIAL DERIVATIVE OF DALPHA   ***
! ***************************************

  function J4_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

! This comes from maximal slicing.

  use param

  implicit none

  real(8) J4_CF,rm
  real(8) psi_rk,Dpsi_rk
  real(8) alpha_rk,Dalpha_rk
  real(8) xi_rk,phi_rk
  real(8) smallpi,V

! Numbers.

  smallpi = acos(-1.d0)

! Potential.

  V = 0.5d0*complex_mass**2*phi_rk**2 + 0.25d0*complex_lambda*phi_rk**4

! Angular momentum.

  if (complex_l>0) then
     V = V + 0.5d0*(dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk**2
  end if

!                       /                 \                 4 /                    2      \
! d Dalpha = - 2 Dalpha | 1/r + d psi/psi | + 8 pi alpha psi  | (omega phi / alpha)  -  V |
!  r                    \        r        /                   \                           /

  J4_CF = - 2.d0*Dalpha_rk*(1.d0/rm + Dpsi_rk/psi_rk) &
        + 8.d0*smallpi*alpha_rk*psi_rk**4*((boson_omega*phi_rk/alpha_rk)**2 - V)

  end function J4_CF








! ************************************
! ***   RADIAL DERIVATIVE OF PHI   ***
! ************************************

! This is just the definition:  xi=dphi/dr

  function J5_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

  implicit none

  real(8) J5_CF,rm
  real(8) psi_rk,Dpsi_rk
  real(8) alpha_rk,Dalpha_rk
  real(8) xi_rk,phi_rk

! dphi/dr = xi.

  J5_CF = xi_rk

  end function J5_CF







! ***********************************
! ***   RADIAL DERIVATIVE OF XI   ***
! ***********************************

  function J6_CF(psi_rk,Dpsi_rk,alpha_rk,Dalpha_rk,phi_rk,xi_rk,rm)

! This comes from the Klein-Gordon equation.

  use param

  implicit none

  real(8) J6_CF,rm
  real(8) psi_rk,Dpsi_rk
  real(8) alpha_rk,Dalpha_rk
  real(8) xi_rk,phi_rk
  real(8) DV

! Derivative of potential.

  DV = complex_mass**2*phi_rk + complex_lambda*phi_rk**3

! Angular momentum.

  if (complex_l>0) then
     DV = DV + (dble(complex_l)*(dble(complex_l)+1.d0)/rm**2)*phi_rk
  end if

!             /                                   \        4 /                       2 \
! d xi = - xi | 2/r + d alpha/alpha + 2 d psi/psi |  +  psi  | V' - phi (omega/alpha)  |
!  r          \        r                 r        /          \                         /

  J6_CF = - xi_rk*(2.d0/rm + Dalpha_rk/alpha_rk + 2.d0*Dpsi_rk/psi_rk) &
        + psi_rk**4*(DV - phi_rk*(boson_omega/alpha_rk)**2)

  end function J6_CF







