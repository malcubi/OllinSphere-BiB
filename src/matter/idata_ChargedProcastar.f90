!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_ChargedProcastar.f90,v 1.6 2024/08/28 17:27:18 malcubi Exp $

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
!
! To obtain the initial data we assume that spacetime is
! static (K_ij=0), and also that the complex scalar and
! vector potentials (Phi,A) have the form:
!
! procaPhi(t,r) = phi(r) exp(-i omega t)
! procaA(t,r)   = i a(r) exp(-i omega t)
!
! In particular, at t=0 procaPhi is purely real while procaA is
! purely imaginary.
!
! Notice that with this ansatz the stress-energy tensor
! is time-independent so the metric can be static.
! The eletric field is also time-independent and real,
! so we are in electrostatics and we can take the Maxwell
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
! We substitute the above ansatz into the Hamiltonian
! constraint to find:
!
! dA/dr  =  A [ (1 - A)/r + 8 pi r A rho ]
!
! where the energy density is now given by:
!
!                              2           2        2                   2          2  
! rho = + 1/(8 pi) { A ( procaE  + maxwellE  )  +  m  [ (procaF / alpha)  +  procaA / A ] }
!
!
! with procaE the Proca electric field, maxwellE the Maxwell
! electric field, and procaF = alpha*procaPhi.
! Notice that this is very similar to the uncharged case, except
! for the contribution from the Maxwell electric field.
!
! For the lapse we use the polar slicing condition
! K_{theta,theta}=0, which when combined with the
! Hamiltonian constraint above takes the form:
!
! dalpha/dr  =  alpha [ (A - 1)/(2r) + 4 pi r A SA ]
!
! where now:
!
!                              2           2        2                   2          2
! SA  = - 1/(8 pi) { A ( procaE  + maxwellE  )  -  m  [ (procaF / alpha)  +  procaA / A ] }
!
!
! Again, this is the same as the uncharged case except for
! the contribution from maxwellE. In the above equations
! procaE is found explicitely to be:
!
!             2
! procaE = - m  alpha procaA/ W A
!
!
! where:  W := omega + q maxwellF
!
! Notice that procaE is defined in such a way that it already
! includes the charge contributions.  This means in particular
! that the Gauss constraint does not have explicit terms
! proportional to the charge (they have been absorved in the
! definition of procaE).
!
! On the other hand, from the Proca field equations we find 
! the following equations for procaF and procaA:
!
!                                     2
! dprocaF/dr = procaA W [ (m alpha/ W) - 1 ]
!
!                                2
! dprocaA/dr =  W procaF A /alpha  -  procaA [ (A+1)/r + 4 pi r A (SA - rho) ]
!
!                 2                   2
!            + q A maxwellE procaE / m
!
! In particular notice that we now have:
!
!                        2          2
! SA - rho  = - A (procaE + maxwellE ) / 4 pi
!
!
! The equation for Maxwell electric field comes from the Gauss
! constraint, and takes the form:
!
!                            2      2
! dmaxwellE/dr =  - q alpha m procaA / (A W)  -  maxwellE [ (5-A) / 2r + 4 pi r A rho ]
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
! A(r=0)        = 1,     d A(r=0)     = 0
!                         r
!
! alpha(r=0)    = 1,     d alpha(r=0) = 0
!                         r
!
! procaF(r=0)   = proca_phi0
!     
!
! procaA(r=0)   = 0
!
!
! maxwellF(r=0) = 0
!
!
! maxwellE(r=0) = 0

! Notice that in the final solution we don't want alpha(r=0)=1,
! but rather alpha=1 at infinity.  But this is no problem as the
! slicing condition above is linear in alpha so we can always
! just rescale the lapse at the end, but in order not to affect
! the solution we must also rescale the final value of the frequency
! omega and both procaF and maxwellF by the same factor.
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

  logical :: left  = .false.             ! This flag is true if the left  value of omega is well behaved.
  logical :: right = .false.             ! This flag is true if the right value of omega is well behaved.

  integer i,l,iter                       ! Counters.
  integer imin                           ! Leftmost grid point.
  integer iaux                           ! Auxiliary quantity.
  integer :: maxiter = 200               ! Maximum number of iterations.

  real(8) r0,delta                       ! Local radius and grid spacing.
  real(8) A0,alpha0                      ! Initial values of (A,alpha).
  real(8) procaF0,procaA0,procaE0        ! Initial values of Proca variables.
  real(8) maxwellF0,maxwellE0            ! Initial values of Maxwell variables. 

  real(8) k11,k12,k13,k14                ! Runge-Kutta sources for A.
  real(8) k21,k22,k23,k24                ! Runge-Kutta sources for alpha.
  real(8) k31,k32,k33,k34                ! Runge-Kutta sources for procaF.
  real(8) k41,k42,k43,k44                ! Runge-Kutta sources for procaA.
  real(8) k51,k52,k53,k54                ! Runge-Kutta sources for maxwellE.
  real(8) k61,k62,k63,k64                ! Runge-Kutta sources for maxwellF.
  real(8) k71,k72,k73,k74                ! Runge-Kutta sources for procaE.

  real(8) A_rk,alpha_rk                  ! Runge-Kutta values of (A,alpha).
  real(8) procaF_rk,procaA_rk,procaE_rk  ! Runge-Kutta values of Proca variables.
  real(8) maxwellF_rk,maxwellE_rk        ! Runge-Kutta values of Maxwell variables.
  real(8) cps_omega_rk                   ! Modified frequency value 
  real(8) J1_CPS,J2_CPS,J3_CPS           ! Functions for sources of differential equations.
  real(8) J4_CPS,J5_CPS,J6_CPS,J7_CPS    ! Functions for sources of differential equations.

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
        print *
     end if

! Alternative normalization.

  else

     proca_phi0  = proca_phi0 /sqrt(4.d0*smallpi)

     if (rank==0) then
        write(*,'(A,ES23.16)') ' Using harmonic normalization at origin: procaPhi(r=0) = proca_phi0 /sqrt(4pi) = ', proca_phi0
        print *
     end if

  end if

! Rescale the tolerance with proca_phi0, so the final tolerance
! is measured with respect to this value.

  epsilon = epsilon*proca_phi0


! ********************************
! ***   CHARGE NORMALIZATION   ***
! ********************************

! Standard normalization.

  if (charge_factor=="standard") then

     if (rank==0) then
        write(*,'(A)') ' Using standard charge normalization'
        print *
     end if

  else

     if (rank==0) then
        write(*,'(A)') ' Using DIFFERENT charge normalization'
        print *
     end if

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
                 
!                For maxwellE and maxwellF we have:  maxwellE' = maxwellF' = 0

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

        write(*,"(A,I4,A,ES23.16,A,ES9.2)") ' Iteration: ',iter,'    Frequency: ',proca_omega,'    Residual: ',res

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

           maxwellE_g(l-1,iaux) = (9.d0*(maxwellE_g(l,i)+maxwellE_g(l,i+1)) - (maxwellE_g(l,i-1)+maxwellE_g(l,i+2)))/16.d0
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
        write(*,'(A,ES23.16)') ' Omega (not-rescaled)      = ', proca_omega
        write(*,'(A,ES23.16)') ' Omega (rescaled)          = ', omega_new
     end if

!    Rescale procaF and maxwellF.

     procaF_g   = procaF_g/alphafac
     maxwellF_g = maxwellF_g/alphafac


!    *******************************************
!    ***   RECONSTRUCT procaPhi AND procaE   ***
!    *******************************************

!    Having found the solution, the scalar potential
!    is recovered as: phi = F/alpha.

     maxwellPhi_g = maxwellF_g/alpha_g
     procaPhi_g   = procaF_g/alpha_g

!    The Proca electric field is given by (for an unperturbed star):
! 
!    E = - m**2 alpha procaA / [A (omega + q*maxwellF)]

     procaE_g = - cproca_mass**2*alpha_g*procaA_g/(A_g*(omega_new + cproca_q*maxwellF_g))


!    *********************************************
!    ***  OUTPUT GAUGE TRANSFORMED FREQUENCY   ***
!    *********************************************

!    Here we output the gauge transformed frequency.
!    But we don't really do the gauge transformation
!    since that should be done AFTER the perturbation
!    if there is one.

!    Figure out asymtotic value of F.

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
!    ***   PROCA STAR PERTURBATION   ***
!    ***********************************

!    Here we add a gaussian perturbation to the solution for the
!    Proca scalar potential procaF, leaving procaA unchanged.
!    We then solve again for (A,alpha,maxwellF,maxwellE,procaE).
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

!       Initialize again (A,alpha,maxwellF,maxwellE,procaE).

        A_g      = 1.d0
        alpha_g  = 1.d0

        maxwellF_g = 0.d0
        maxwellE_g = 0.d0

        procaE_g = 0.d0

!       Rescale back procaF (otherwise it won't be consistent any more).

        procaF_g = procaF_g*alphafac

!       Add perturbation to proca_Phi.

        if (proca_phiR_r0==0.d0) then
           procaF_g = procaF_g + proca_PhiR_a0*exp(-rr**2/proca_phiR_s0**2)
        else
           procaF_g = procaF_g + proca_PhiR_a0 &
                    *(exp(-(rr-proca_PhiR_r0)**2/proca_PhiR_s0**2) &
                    + exp(-(rr+proca_PhiR_r0)**2/proca_PhiR_s0**2))
        end if

!       Loop over grid levels. We solve from fine to coarse grid.

        do l=Nl-1,0,-1

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

!                Values of (maxwellF,maxwellE) at origin.

                 maxwellF0 = 0.d0
                 maxwellE0 = 0.d0

!                Values of (procaF,procaA,procaE) at origin.

                 procaF0 = (9.d0*(procaF_g(l,0)+procaF_g(l,1)) - (procaF_g(l,-1)+procaF_g(l,2)))/16.d0
                 procaA0 = 0.d0
                 procaE0 = 0.d0

!             Grid spacing and values at previous grid point.

              else

                 delta = dr(l)
                 r0    = rr(l,i-1)

                 A0        = A_g(l,i-1)
                 alpha0    = alpha_g(l,i-1)

                 maxwellF0 = maxwellF_g(l,i-1)
                 maxwellE0 = maxwellE_g(l,i-1)

                 procaF0   = procaF_g(l,i-1)
                 procaA0   = procaA_g(l,i-1)
                 procaE0   = procaE_g(l,i-1)

              end if

!             I) First Runge-Kutta step.

              if (i==1) then

!                At the origin we have:  A' = alpha' = 0.

                 k11 = 0.d0
                 k21 = 0.d0

!                For maxwellE and maxwellF we have:  maxwellE' = maxwellF' = 0

                 k51 = 0.d0
                 k61 = 0.d0

!                If we take E ~ k r with k constant close to
!                the origin we find from the Gauss constraint:
!
!                k = - m^2 Phi(r=0) / 3.

                 k71 = - delta*cproca_mass**2*procaF0/3.d0

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk        = A0
                 alpha_rk    = alpha0

                 maxwellF_rk = maxwellF0
                 maxwellE_rk = maxwellE0

                 procaF_rk   = procaF0
                 procaA_rk   = procaA0
                 procaE_rk   = procaE0

!                Sources.

                 k11 = delta*J1_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
                 k21 = delta*J2_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
                 k51 = delta*J5_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
                 k61 = delta*J6_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
                 k71 = delta*J7_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk        = A0        + half*k11
              alpha_rk    = alpha0    + half*k21
              maxwellE_rk = maxwellE0 + half*k51
              maxwellF_rk = maxwellF0 + half*k61
              procaE_rk   = procaE0   + half*k71

              if (i==1) then  ! Linear interpolation for first point.
                 procaF_rk = 0.5d0*(procaF0 + procaF_g(l,1))
                 procaA_rk = 0.5d0*(procaA0 + procaA_g(l,1))
              else            ! Cubic interpolation for the rest.
                 procaF_rk = (9.d0*(procaF_g(l,i)+procaF_g(l,i-1)) - (procaF_g(l,i-2)+procaF_g(l,i+1)))/16.d0
                 procaA_rk = (9.d0*(procaA_g(l,i)+procaA_g(l,i-1)) - (procaA_g(l,i-2)+procaA_g(l,i+1)))/16.d0
              end if

!             Sources.

              k12 = delta*J1_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k22 = delta*J2_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k52 = delta*J5_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k62 = delta*J6_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k72 = delta*J7_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)

!             III) Third Runge-Kutta step.

              A_rk        = A0        + half*k12
              alpha_rk    = alpha0    + half*k22
              maxwellE_rk = maxwellE0 + half*k52
              maxwellF_rk = maxwellF0 + half*k62
              procaE_rk   = procaE0   + half*k72

!             Sources.

              k13 = delta*J1_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k23 = delta*J2_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k53 = delta*J5_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k63 = delta*J6_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k73 = delta*J7_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk        = A0        + k13
              alpha_rk    = alpha0    + k23
              maxwellE_rk = maxwellE0 + k53
              maxwellF_rk = maxwellF0 + k63
              procaE_rk   = procaE0   + k73

              procaF_rk = procaF_g(l,i)
              procaA_rk = procaA_g(l,i)

!             Sources.

              k14 = delta*J1_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k24 = delta*J2_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k54 = delta*J5_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k64 = delta*J6_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)
              k74 = delta*J7_CPS(A_rk,alpha_rk,procaF_rk,procaA_rk,maxwellE_rk,maxwellF_rk,procaE_rk,rm)

!             Advance variables to next grid point.

              A_g(l,i)        = A0        + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i)    = alpha0    + (k21 + 2.d0*(k22 + k23) + k24)/6.d0
              maxwellE_g(l,i) = maxwellE0 + (k51 + 2.d0*(k52 + k53) + k54)/6.d0
              maxwellF_g(l,i) = maxwellF0 + (k61 + 2.d0*(k62 + k63) + k64)/6.d0
              procaE_g(l,i)   = procaE0   + (k71 + 2.d0*(k72 + k73) + k74)/6.d0

           end do

!          Fix ghost zones.

           do i=1,ghost
              A_g(l,1-i)        = + A_g(l,i)
              alpha_g(l,1-i)    = + alpha_g(l,i)
              maxwellF_g(l,1-i) = + maxwellF_g(l,i)      
              maxwellE_g(l,1-i) = - maxwellE_g(l,i)
              procaE_g(l,1-i)   = - procaE_g(l,i)
           end do

        end do

!       Restrict solution from fine to coarse grid.

        do l=Nl-1,1,-1

           do i=1,Nrtotal-ghost,2

              iaux = i/2 + 1
              rm = rr(l-1,iaux)

              A_g(l-1,iaux)        = (9.d0*(A_g(l,i)+A_g(l,i+1)) - (A_g(l,i-1)+A_g(l,i+2)))/16.d0
              alpha_g(l-1,iaux)    = (9.d0*(alpha_g(l,i)+alpha_g(l,i+1)) - (alpha_g(l,i-1)+alpha_g(l,i+2)))/16.d0

              maxwellE_g(l-1,iaux) = (9.d0*(maxwellE_g(l,i)+maxwellE_g(l,i+1)) - (maxwellE_g(l,i-1)+maxwellE_g(l,i+2)))/16.d0
              maxwellF_g(l-1,iaux) = (9.d0*(maxwellF_g(l,i)+maxwellF_g(l,i+1)) - (maxwellF_g(l,i-1)+maxwellF_g(l,i+2)))/16.d0

              procaE_g(l-1,iaux)   = (9.d0*(procaE_g(l,i)+procaE_g(l,i+1)) - (procaE_g(l,i-1)+procaE_g(l,i+2)))/16.d0

           end do

!          Fix ghost zones.

           do i=1,ghost
              A_g(l-1,1-i)        = + A_g(l-1,i)
              alpha_g(l-1,1-i)    = + alpha_g(l-1,i)
              maxwellF_g(l-1,1-i) = + maxwellF_g(l-1,i)
              maxwellE_g(l-1,1-i) = - maxwellE_g(l-1,i)
              procaE_g(l-1,1-i)   = - procaE_g(l-1,i)
           end do

        end do

!       Rescale lapse, maxwellF and procaF again.

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
        procaF_g = procaF_g/alphafac

!       Reconstruct maxwellPhi and procaPhi.

        maxwellPhi_g = maxwellF_g/alpha_g
        procaPhi_g   = procaF_g/alpha_g

     end if


!    ********************************
!    ***   GAUGE TRANSFORMATION   ***
!    ********************************

!    In a similar way to the lapse, we integrated taking maxwellF(r=0)=0,
!    but we really want maxwellF=0 at infinity.  We can fix this by making
!    a gauge transformation of the form:
!
!    Fm -> Fm - F_infty ,   omega ->  omega + q F_infty
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

! For the energy density we have:
!
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

! For SA we have:
!
!                              2           2        2                   2          2
! SA  = - 1/(8 pi) { A ( procaE  + maxwellE  )  -  m  [ (procaF / alpha)  +  procaA / A ] }
!
! Notice that we don't divide by 8*pi since it cancels. 

  SA = - A*(procaE**2 + maxwellE**2) + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A)

! dalpha/dr  =  alpha [ (A - 1)/(2r) + 4 pi r A SA ]
! 
! Notice that we don't multiply the last term with 4*pi*A since
! it cancels (we are left only with 1/2 which we factor out).

  J2_CPS = 0.5d0*alpha*((A - 1.d0)/rm + rm*A*SA)

  end function J2_CPS







! ***************************************
! ***   RADIAL DERIVATIVE OF procaF   ***
! ***************************************

! The radial derivative of procaF comes from the
! Proca equations.

  function J3_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J3_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm
  real(8) W

! dprocaF/dr = procaA W [ (m alpha/ W)**2 - 1 ]
!
! where:  W  :=  omega + q maxwellF

  W = proca_omega + cproca_q*maxwellF

  J3_CPS = procaA*W*((cproca_mass*alpha/W)**2 - 1.d0)

  end function J3_CPS







! ***************************************
! ***   RADIAL DERIVATIVE OF procaA   ***
! ***************************************

! The radial derivative of procaA from the
! Proca equations.

  function J4_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J4_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm
  real(8) aux,W
  
!                        2          2
! SA - rho  = - A (procaE + maxwellE ) / 4 pi
! 
! Notice that we don't divide by 4*pi since it cancels.
 
  aux = - A*(procaE**2 + maxwellE**2)
              
! dprocaA/dr =  W procaF A /alpha**2  -  procaA [ (A+1)/r + 4 pi r A (SA - rho) ]
!
!                 2                   2
!            + q A maxwellE procaE / m
!
!
! where:  W  :=  omega + q maxwellF

  W = proca_omega + cproca_q*maxwellF

  J4_CPS = W*procaF*A/alpha**2 - procaA*((A + 1.d0)/rm + rm*A*aux) &
         + cproca_q*A**2*maxwellE*procaE/cproca_mass**2

  end function J4_CPS










! *****************************************
! ***   RADIAL DERIVATIVE OF maxwellE   ***
! *****************************************

! The radial derivative of maxwellE comes from
! the Maxwell equations (the Gauss constraint).

  function J5_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J5_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm
  real(8) rho,ecurrent

! The radial derivative of maxwellE comes from the Gauss constraint
! and takes the form:
!
! dmaxwellE/dr =  - E [ 2/r + (dA/dr) / 2A ] + 4 pi ecurrent
!
!              =  - maxwellE [ (5-A) / 2r + 4 pi r A rho ] + 4 pi ecurrent
!
! where in the second equality we used the Hamiltonian constraint
! to eliminate dA/dr, and with "ecurrent" the electric current of
! the Proca field.
!
! For the energy density we have:
!
!                              2           2        2                   2          2  
! rho = + 1/(8 pi) { A ( procaE  + maxwellE  )  +  m  [ (procaF / alpha)  +  procaA / A ] }
!
!
! And for the electric current we have:
!
! ecurrent = (q / 4 pi) procaA procaE

  rho = A*(procaE**2 + maxwellE**2) + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A)

  ecurrent = cproca_q*procaA*procaE

! Notice that we don't divide by the factors of 8*pi and 4*pi
! in the last expressions since they cancel.

  J5_CPS = - maxwellE*((5.d0-A)/(2.d0*rm) + 0.5d0*rm*A*rho) + ecurrent

  end function J5_CPS







! *****************************************
! ***   RADIAL DERIVATIVE OF maxwellF   ***
! *****************************************

! The radial derivative of maxwellF comes from
! the Maxwell equations.

  function J6_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J6_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm

! dmaxwellF/dr  =  - alpha A maxwellE

  J6_CPS = - alpha*A*maxwellE

  end function J6_CPS







! ***************************************
! ***   RADIAL DERIVATIVE OF procaE   ***
! ***************************************

! The radial derivative of procaE comes from the
! Proca equations.

  function J7_CPS(A,alpha,procaF,procaA,maxwellE,maxwellF,procaE,rm)

  use param
  
  implicit none
  
  real(8) J7_CPS
  real(8) A,alpha,procaF,procaA,procaE,maxwellE,maxwellF,rm
  real(8) rho,aux

! We use the fact that:
!
! dA/dr  =  A [ (1 - A)/r + 8 pi r A rho ]
!
! with rho given by:
!                              2           2        2                   2          2  
! rho = + 1/(8 pi) { A ( procaE  + maxwellE  )  +  m  [ (procaF / alpha)  +  procaA / A ] }
!
! Notice that we don't divide by 8*pi since it cancels.

  rho = A*(procaE**2 + maxwellE**2) + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A)
  aux = (1.d0-A)/rm + rm*A*rho

! Proca Gauss constraint:
!
!                                                  2
! dprocE/dr  =  - procaE [ 2/r + (1/2A) dA/dr ] - m  procaF/alpha

  J7_CPS = - procaE*(2.d0/rm + 0.5d0*aux) - cproca_mass**2*procaF/alpha

  end function J7_CPS



