!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_LProcastar.f90,v 1.4 2024/02/02 17:32:56 malcubi Exp $

  subroutine idata_LProcastar

! *************************************
! ***   L-PROCA STAR INITIAL DATA   ***
! *************************************

! This subroutine calculates initial data for a l-Proca star
! using a shooting method in the polar/areal gauge.


! MAIN PROBLEM: THIS ROUTINE DOES NOT WORK WELL IN PRACTICE,
! SINCE WE HAVE FOUND THAT L-PROCA STARS REQUIRE US TO FIX
! TWO PARAMETERS, SO A SIMPLE SHOOTING FAILS.
!
! CLAUDIO LAZARTE HAS WRITTEN A SPECTRAL CODE THAT SOLVES THIS
! PROBLEM, BUT IT IS STILL NOT ADDED TO THE MAIN OLLINSPHERE CODE.


! L-Proca stars are self-gravitating solutions of a collection
! of odd number 2l+1 of complex massive vector fields (Proca fields) 
! such that the spacetime is static and spherically symmetric.  
! Those Proca fields have the same harmonic time dependence and the 
! same radial profile, and their angular dependence given
! by spherical harmonics Y^{lm} with a fixed azimuthal quantum
! number l and all possible values of the quantum number m=-l,...,+l.
! (for details see the Master thesis of Claudio Lazarte, or the
! paper in preparation).
!
! To obtain the initial data we assume that spacetime is
! static (K_ij=0), and also that the complex scalar and
! vector potentials (phi,A) for each (lm) have the form:
!
! Phi  =  phi(r) Y^{lm} exp(-i omega t) 
! A    =  i a(r) Y^{lm} exp(-i omega t) 
!  
! Where a(r) := procaA is the radial component of the vector potential.
! These two functions are enough to describe l-Proca stars with l=0, 
! i.e. a "single field" Proca star.  However, when l is not 0 the 
! vector potential has non-zero angular components, and both the polar 
! and the azimutal component are in terms of a function:
!
! B  =  i b(r) dY^{lm} exp(-i omega t) 
!
! where dY are the angular derivatives of the spherical harmonics.
! In those cases, the effect in the equations is adding centrifugal
! terms proportional to (Phi,A,B,G) with the factor l(l+1)/r^2. 
!    
! Notice that with this ansatz the total stress-energy tensor is both
! time-independent and spherically symmetric.  Also, at t=0 the
! scalar potential phi is purely real, while the components of
! the vector potential A and B are purely imaginary.
!
! The standard ansatz for the metric for l-Proca stars is:
!
!   2          2   2        2     2      2
! ds  = - alpha  dt  +  A dr  +  r dOmega
!
! In other words, we are assuming psi=B=1, beta=0, and the radial
! coordinate is the area radius.
!
! It is in fact much better to work with the quantity F := alpha*phi
! instead of phi itself, and then reconstruct phi at the end.
! With the above ansatz, and using F, the l-Proca evolution equations
! reduce to:
!
!                                     2                      2                  2 
! dF/dr  =  a omega[ ( m alpha/omega )  -  1 ] + l(l+1) (alpha/omega) (a - G) /r 
!
!                            2                 2  
! da/dr  =  omega F A / alpha +  l(l+1) A b / r  -  a [ (A+1)/r + 4 pi r A (SA - rho) ]
!
!     
! db/dr  =  G 
!
!              2                2           2
! dG/dr  =  [ m  - (omega/alpha)  + l(l+1)/r  ] A b  - 2a/r  - G [ (A-1)/r + 4 pi r A (SA - rho) ]   
!              
!
! The last two equations are valid only for l different from 0.
! Furthermore, the energy density rho is given by:
!
!                           2     2           2    2                  2    2 2                         2       2
! rho  =    1 / (8 pi) { A E  +  m [ (F/alpha)  + a / A ]  +  l(l+1)/r  [ m b + ( (omega b + F)/alpha ) + (a-G) /A ] }
!
!
! and the radial stress SA is given by:
!
!                             2     2           2    2                 2     2 2                         2       2 
! SA   =    1 / (8 pi) { - A E  +  m [ (F/alpha)  + a / A ] +  l(l+1)/r  [ -m b + ( (omega b + F)/alpha ) + (a-G) /A ] }
!
!
! with radial contravariant component of electric field :
!
!                                2                     2
! E  =  -  alpha /(omega A) [  m  a  + l(l+1) (a-G) / r  ]
!
!
! In particular:
!
!                               2           2 2   2
! SA - rho  = - 1 / (4 pi) { A E  + l(l+1) m b / r  }
!
!
! On the other hand, the Hamiltonian constraint takes the form:
!
!                      
! dA/dr  =  A [ (1 - A)/r + 8 pi r A rho ]
!
!
! For the lapse we use the polar slicing condition K_{theta,theta}=0,
! which when combined with the hamiltonian constraint above takes the form:
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
! The boundary conditions at the origin for (alpha,A) are:
!
! A(r=0)     = 1,     d A(r=0)     = 0
!                      r
!
! alpha(r=0) = 1,     d alpha(r=0) = 0
!                      r
!
! And for Proca fields we have:
!
! F(r=0) ~ phi0 r**l
!
! a(r=0) ~ a0 r**(l+1)
!
! b(r=0) ~ b0 r**(l+2)
!      
! it is noted that for l>=1 it turns out that the fourth order 
! Runge-Kutta is inconsistent at the origin due to the behavior of
! the fields, so in that case we solve instead for the rescaled fields  
!
! F_res = F/r**l,
! a_res = a/r**(l+1),  
! b_res = b/r**(l+2), and     
! G_res = G/r**(l+1). 
!    
! Therefore, the boundary conditions for these reescaled fields are
!      
! F_res(r=0) = phi0,                               d F_res(r=0) = 0
!                                                   r     
!
! a_res(r=0) = omega phi0/(2l+3),                  d a_res(r=0) = 0
!                                                   r     
!
! b_res(r=0) = - omega phi0/[(2l+3)(l+1)],         d b_res(r=0) = 0
!                                                   r 
!                                                                           
! G_res(r=0) = - (l+2) omega phi0/[(2l+3)(l+1)],   d G_res(r=0) = 0     
!                                                   r   
!         
!
! Notice that in the final solution we don't want alpha(r=0)=1,
! but rather alpha=1 at infinity.  But this is no problem as the
! slicing condition above is linear in alpha so we can always
! just rescale the lapse at the end, but in order not to affect
! the solution for (A,phi,a,b) we must also rescale the final value
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
  real(8) procaFres0,procaAres0         ! Initial values of reescaled variables.
  real(8) procaBres0,procaGres0         ! Initial values of reescaled variables.  
  real(8) k11,k12,k13,k14               ! Runge-Kutta sources for A.
  real(8) k21,k22,k23,k24               ! Runge-Kutta sources for alpha.
  real(8) k31,k32,k33,k34               ! Runge-Kutta sources for procaFres.
  real(8) k41,k42,k43,k44               ! Runge-Kutta sources for procaAres.
  real(8) k51,k52,k53,k54               ! Runge-Kutta sources for procaBres.
  real(8) k61,k62,k63,k64               ! Runge-Kutta sources for procaGres.
  real(8) A_rk,alpha_rk                 ! Runge-Kutta values of variables.
  real(8) procaFres_rk,procaAres_rk     ! Runge-Kutta values of reescaled variables.
  real(8) procaBres_rk,procaGres_rk     ! Runge-Kutta values of reescales variables.
  real(8) J1_lproca,J2_lproca           ! Functions for sources of differential equations.
  real(8) J3_lproca,J4_lproca           ! Functions for sources of differential equations.
  real(8) J5_lproca,J6_lproca           ! Functions for sources of differential equations.
  real(8) res_old,res                   ! Residual.
  real(8) omega_new,omega_old,domega    ! Trial frequency and frequency interval. 
  real(8) half,smallpi                  ! Numbers.
  real(8) rm,alphafac,aux               ! Auxiliary quantities.
  real(8) :: epsilon = 1.d-8            ! Tolerance.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                      ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g,alpha_g             ! Global arrays for radial metric and lapse.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaPhi_g,procaA_g     ! Global arrays for potentials.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaB_g,procaG_g       ! Global arrays for potentials.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaF_g                ! Global array for F = alpha*phi.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaE_g                ! Global array for radial electric field.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaXi_g               ! Global array for angular electric field.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaFres_g,procaAres_g ! Global arrays for reescaled potentials.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: procaBres_g,procaGres_g ! Global arrays for reescaled potentials.


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  smallpi = acos(-1.d0)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for a l-Proca star in the polar/areal gauge using a shooting method ...'
     print *
  end if


! *************************
! ***   NORMALIZATION   ***
! *************************

! Physical normalization.

  if (proca_factor=="physical") then

     if (rank==0) then
        if (cproca_l==0) then
           write(*,'(A)') ' Using "physical" normalization at origin: phi(r=0) = phi0'
           print *
        else
           write(*,'(A)') ' Using "physical" normalization at origin: phi ~ phi0 r**l'
           print *
     end if
  end if

! Harmonic normalization.

  else

     proca_phi0 = proca_phi0/sqrt(dble(2*cproca_l+1))/sqrt(4.d0*smallpi)

     if (rank==0) then
        if (cproca_l==0) then
           write(*,'(A,E16.10)') ' Using "harmonic" normalization at origin: phi(r=0) = phi0/sqrt(4pi) = ', proca_phi0
           print *
        else
           write(*,'(A,E16.10)') ' Using "harmonic" normalization at origin: phi/r**l ~ phi0/sqrt(4pi*(2l+1)) = ', proca_phi0
           print *
        end if
     end if

  end if

! Rescale the tolerance with proca_phi0, so the final
! tolerance is measured with respect to this value.

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

  procaPhi_g = 0.d0
  procaF_g   = 0.d0
  procaA_g   = 0.d0
  procaB_g   = 0.d0
  procaG_g   = 0.d0
  procaE_g   = 0.d0
  procaXi_g  = 0.d0

  procaFres_g = 0.d0
  procaAres_g = 0.d0
  procaBres_g = 0.d0
  procaGres_g = 0.d0


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
        print *, 'For l-Proca star initial data we must have omega_right>omega_left.'
        print *, 'Aborting! (subroutine idata_Lprocastar)'
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

              procaFres_g(l,imin-1) = (9.d0*(procaFres_g(l+1,Nrtotal-2)+procaFres_g(l+1,Nrtotal-3)) &
                                    - (procaFres_g(l+1,Nrtotal-4)+procaFres_g(l+1,Nrtotal-1)))/16.d0
              procaAres_g(l,imin-1) = (9.d0*(procaAres_g(l+1,Nrtotal-2)+procaAres_g(l+1,Nrtotal-3)) &
                                    - (procaAres_g(l+1,Nrtotal-4)+procaAres_g(l+1,Nrtotal-1)))/16.d0

              procaBres_g(l,imin-1) = 0              
              procaGres_g(l,imin-1) = 0

              if (cproca_l/=0) then 
                 procaBres_g(l,imin-1) = (9.d0*(procaBres_g(l+1,Nrtotal-2)+procaBres_g(l+1,Nrtotal-3)) &
                                       - (procaBres_g(l+1,Nrtotal-4)+procaBres_g(l+1,Nrtotal-1)))/16.d0
                 procaGres_g(l,imin-1) = (9.d0*(procaGres_g(l+1,Nrtotal-2)+procaGres_g(l+1,Nrtotal-3)) &
                                       - (procaGres_g(l+1,Nrtotal-4)+procaGres_g(l+1,Nrtotal-1)))/16.d0
              end if               

              procaF_g(l,imin-1) = procaFres_g(l,imin-1)*rr(l,imin-1)**cproca_l
              procaA_g(l,imin-1) = procaAres_g(l,imin-1)*rr(l,imin-1)**(cproca_l+1)
              procaB_g(l,imin-1) = procaBres_g(l,imin-1)*rr(l,imin-1)**(cproca_l+2)
              procaG_g(l,imin-1) = procaGres_g(l,imin-1)*rr(l,imin-1)**(cproca_l+1) 

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

!                Values of (procaFres,procaAres,procaBres,procaGres) at origin.

                 procaFres0 = proca_phi0
                 procaAres0 = proca_omega*proca_phi0/dble(2*cproca_l+3)
                 procaBres0 = 0
                 procaGres0 = 0

                 if (cproca_l/=0) then   
                    procaBres0 = - proca_omega*proca_phi0/dble(2*cproca_l+3)/dble(cproca_l+1)
                    procaGres0 = - dble(cproca_l+2)*proca_omega*proca_phi0/dble(2*cproca_l+3)/dble(cproca_l+1)
                 end if 

!             Grid spacing and values at previous grid point.

              else

                 delta = dr(l)
                 r0    = rr(l,i-1)

                 A0     = A_g(l,i-1)
                 alpha0 = alpha_g(l,i-1)

                 procaFres0 = procaFres_g(l,i-1)
                 procaAres0 = procaAres_g(l,i-1)
                 procaBres0 = 0
                 procaGres0 = 0

                 if (cproca_l/=0) then
                    procaBres0 = procaBres_g(l,i-1)
                    procaGres0 = procaGres_g(l,i-1)
                 end if

              end if

!             I) First Runge-Kutta step.

!             Sources at first grid point if we start
!             from the origin (for finer grid).

              if (i==1) then

!                At the origin we have:  A' = alpha' = procaFres' = procaAres' = 0.

                 k11 = 0.d0
                 k21 = 0.d0
                 k31 = 0.d0
                 k41 = 0.d0

                 if (cproca_l/=0) then

!                   At the origin we have:  procaBres' = procaGres'= 0

                    k51 = 0.d0
                    k61 = 0.d0

                 end if

!             Sources at previous grid point.

              else

                 rm = r0

                 A_rk     = A0
                 alpha_rk = alpha0

                 procaFres_rk = procaFres0
                 procaAres_rk = procaAres0

                 procaBres_rk = 0
                 procaGres_rk = 0

                 if (cproca_l/=0) then
                    procaBres_rk = procaBres0
                    procaGres_rk = procaGres0
                 end if

!                Sources.

                 k11 = delta*J1_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
                 k21 = delta*J2_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
                 k31 = delta*J3_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
                 k41 = delta*J4_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)

                 if(cproca_l/=0) then 
                    k51 = delta*J5_lproca(procaBres_rk,procaGres_rk,rm)
                    k61 = delta*J6_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
                 end if 

              end if

!             II) Second Runge-Kutta step.

              rm = r0 + half*delta

              A_rk     = A0     + half*k11
              alpha_rk = alpha0 + half*k21

              procaFres_rk = procaFres0 + half*k31
              procaAres_rk = procaAres0 + half*k41

              procaBres_rk = 0
              procaGres_rk = 0

              if (cproca_l/=0) then
                  procaBres_rk = procaBres0 + half*k51
                  procaGres_rk = procaGres0 + half*k61
              end if

!             Sources.

              k12 = delta*J1_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              k22 = delta*J2_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              k32 = delta*J3_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              k42 = delta*J4_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)

              if (cproca_l/=0) then 
                  k52 = delta*J5_lproca(procaBres_rk,procaGres_rk,rm)
                  k62 = delta*J6_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              end if 

!             III) Third Runge-Kutta step.

              A_rk     = A0     + half*k12
              alpha_rk = alpha0 + half*k22

              procaFres_rk = procaFres0 + half*k32
              procaAres_rk = procaAres0 + half*k42

              procaBres_rk = 0
              procaGres_rk = 0 

              if (cproca_l/=0) then 
                  procaBres_rk = procaBres0 + half*k52
                  procaGres_rk = procaGres0 + half*k62
              end if  

!             Sources.

              k13 = delta*J1_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              k23 = delta*J2_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              k33 = delta*J3_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              k43 = delta*J4_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)

              if (cproca_l/=0) then 
                  k53 = delta*J5_lproca(procaBres_rk,procaGres_rk,rm)
                  k63 = delta*J6_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              end if 

!             IV) Fourth Runge-Kutta step.

              rm = r0 + delta

              A_rk     = A0     + k13
              alpha_rk = alpha0 + k23

              procaFres_rk = procaFres0 + k33
              procaAres_rk = procaAres0 + k43

              procaBres_rk = 0
              procaGres_rk = 0 

              if (cproca_l/=0) then 
                  procaBres_rk = procaBres0 + k53
                  procaGres_rk = procaGres0 + k63
              end if  

!             Sources.

              k14 = delta*J1_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              k24 = delta*J2_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              k34 = delta*J3_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              k44 = delta*J4_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)

              if (cproca_l/=0) then 
                  k54 = delta*J5_lproca(procaBres_rk,procaGres_rk,rm)
                  k64 = delta*J6_lproca(A_rk,alpha_rk,procaFres_rk,procaAres_rk,procaBres_rk,procaGres_rk,rm)
              end if 

!             Advance variables to next grid point.

              A_g(l,i)     = A0     + (k11 + 2.d0*(k12 + k13) + k14)/6.d0
              alpha_g(l,i) = alpha0 + (k21 + 2.d0*(k22 + k23) + k24)/6.d0

              procaFres_g(l,i) = procaFres0 + (k31 + 2.d0*(k32 + k33) + k34)/6.d0
              procaAres_g(l,i) = procaAres0 + (k41 + 2.d0*(k42 + k43) + k44)/6.d0

              procaBres_g(l,i) = 0
              procaGres_g(l,i) = 0

              if (cproca_l/=0) then 
                  procaBres_g(l,i) = procaBres0 + (k51 + 2.d0*(k52 + k53) + k54)/6.d0
                  procaGres_g(l,i) = procaGres0 + (k61 + 2.d0*(k62 + k63) + k64)/6.d0
              end if 

!             Reconstruct physical potentials.

              procaF_g(l,i) = procaFres_g(l,i)*rr(l,i)**cproca_l
              procaA_g(l,i) = procaAres_g(l,i)*rr(l,i)**(cproca_l+1)
              procaB_g(l,i) = procaBres_g(l,i)*rr(l,i)**(cproca_l+2)
              procaG_g(l,i) = procaGres_g(l,i)*rr(l,i)**(cproca_l+1)

!             Check if solution is blowing up. This helps to reduce
!             the need for a very fine tuned initial guess.

              if (abs(procaFres_g(l,i))>2.d0*proca_phi0) then
                 if (.not.left) then
                    omega_left = omega_left + domega
                    proca_omega = omega_left
                    goto 100
                 else if (.not.right) then
                    omega_right = omega_right - domega
                    proca_omega = omega_right
                    goto 100
                 end if
              end if

!             Check if solution is already very small.

              if (abs(procaFres_g(l,i))+abs(procaFres_g(l,i-1))<epsilon) then
                 procaFres_g(l,i) = 0.d0
                 procaF_g(l,i) = 0.d0
              end if

           end do


!          ***********************
!          ***   GHOST ZONES   ***
!          ***********************

           do i=1,ghost

              A_g(l,1-i)     = + A_g(l,i)
              alpha_g(l,1-i) = + alpha_g(l,i)

              procaF_g(l,1-i) = +(1.d0-2.d0*mod(cproca_l,2))*procaF_g(l,i)
              procaA_g(l,1-i) = +(1.d0-2.d0*mod(cproca_l+1,2))*procaA_g(l,i)
              procaB_g(l,1-i) = +(1.d0-2.d0*mod(cproca_l,2))*procaB_g(l,i)
              procaG_g(l,1-i) = +(1.d0-2.d0*mod(cproca_l+1,2))*procaG_g(l,i)

              procaFres_g(l,1-i) = + procaFres_g(l,i)
              procaAres_g(l,1-i) = + procaAres_g(l,i)
              procaBres_g(l,i-1) = + procaBres_g(l,i)
              procaGres_g(l,i-1) = + procaGres_g(l,i) 

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

        if (abs(procaFres_g(0,Nrtotal)) + abs(procaFres_g(0,Nrtotal-1))<epsilon) then
           res = epsilon/2.d0
           goto 100
        else
           aux = abs(cproca_mass**2 - (boson_omega/alpha_g(0,Nrtotal))**2)
           res = (procaF_g(0,Nrtotal) - procaF_g(0,Nrtotal-1))/dr(0) + dsqrt(aux)*procaF_g(0,Nrtotal)
           !res = (procaB_g(0,Nrtotal) - procaB_g(0,Nrtotal-1))/dr(0) + dsqrt(aux)*procaB_g(0,Nrtotal)
        end if

!       Secant method:  Having found the difference for the two values of omega
!       we can use it to extrapolate linearly to the next best guess.

        if (.not.left) then

           left = .true.

           omega_old   = proca_omega
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
        print *, 'Aborting! (subroutine idata_LProcastar)'
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

           procaFres_g(l-1,iaux) = (9.d0*(procaFres_g(l,i)+procaFres_g(l,i+1)) - (procaFres_g(l,i-1)+procaFres_g(l,i+2)))/16.d0
           procaAres_g(l-1,iaux) = (9.d0*(procaAres_g(l,i)+procaAres_g(l,i+1)) - (procaAres_g(l,i-1)+procaAres_g(l,i+2)))/16.d0

           procaBres_g(l-1,iaux) = 0
           procaGres_g(l-1,iaux) = 0

           if (cproca_l/=0) then
              procaBres_g(l-1,iaux) = (9.d0*(procaBres_g(l,i)+procaBres_g(l,i+1)) - (procaBres_g(l,i-1)+procaBres_g(l,i+2)))/16.d0
              procaGres_g(l-1,iaux) = (9.d0*(procaGres_g(l,i)+procaGres_g(l,i+1)) - (procaGres_g(l,i-1)+procaGres_g(l,i+2)))/16.d0
           end if 

           procaF_g(l-1,iaux) = procaFres_g(l-1,iaux)*rr(l-1,iaux)**cproca_l
           procaA_g(l-1,iaux) = procaAres_g(l-1,iaux)*rr(l-1,iaux)**(cproca_l+1)
           procaB_g(l-1,iaux) = procaBres_g(l-1,iaux)*rr(l-1,iaux)**(cproca_l+2)
           procaG_g(l-1,iaux) = procaGres_g(l-1,iaux)*rr(l-1,iaux)**(cproca_l+1) 

        end do

!       Fix ghost zones.

        do i=1,ghost

           A_g(l-1,1-i)      = + A_g(l-1,i)
           alpha_g(l-1,1-i)  = + alpha_g(l-1,i)

           procaF_g(l-1,1-i) = +(1.d0-2.d0*mod(cproca_l,2))*procaF_g(l-1,i)
           procaA_g(l-1,1-i) = +(1.d0-2.d0*mod(cproca_l+1,2))*procaA_g(l-1,i)
           procaB_g(l-1,1-i) = +(1.d0-2.d0*mod(cproca_l,2))*procaB_g(l-1,i)
           procaG_g(l-1,1-i) = +(1.d0-2.d0*mod(cproca_l+1,2))*procaG_g(l-1,i)

           procaFres_g(l-1,1-i) = + procaFres_g(l-1,i)
           procaAres_g(l-1,1-i) = + procaAres_g(l-1,i)
           procaBres_g(l-1,i-1) = + procaBres_g(l-1,i)
           procaGres_g(l-1,i-1) = + procaGres_g(l-1,i)

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
        write(*,'(A,I1)' )    ' l-Procastar parameter  = ', cproca_l
        write(*,'(A,E22.16)') ' Omega (not-rescaled)   = ', proca_omega
        write(*,'(A,E22.16)') ' Omega (rescaled)       = ', omega_new
        print *
     end if

!    Rescale procaF.

     procaFres_g = procaFres_g/alphafac

     procaF_g = procaF_g/alphafac


!    ****************************************************
!    ***   RECONSTRUCT procaPhi, procaE and procaXi   ***
!    ****************************************************

!    Having found the solution, the scalar potential
!    is recovered as: phi = F/alpha.

     procaPhi_g = procaF_g/alpha_g

!    Radial component (index up) of electric field is given by:
!
!    E = - (omega*procaA + dF/dr)/(alpha*A)
!
!      = - alpha [ m**2 procaA + l(l+1)(procaA -procaG)/r**2] / (A omega)

     procaE_g = - alpha_g*(cproca_mass**2*procaA_g &
              + dble(cproca_l)*(dble(cproca_l)+1.d0)*(procaA_g-procaG_g)/rr**2)/(A_g*omega_new)

!    Angular component (index down) of electric field is given by:
!
!    Xi = - (omega*procaB + alpha*procaPhi)/alpha 
    
     procaXi_g = - (omega_new*procaB_g + alpha_g*procaPhi_g)/alpha_g


!    *************************************
!    ***   L-PROCA STAR PERTURBATION   ***
!    *************************************
!     
!    Comming soon! ...


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
     cprocaB_I   = procaB_g
     cprocaG_I   = procaG_g

     cprocaE_R   = procaE_g
     cprocaXi_R  = procaXi_g

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! arrays with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  else

     call distribute(0,Nl-1,A,A_g)
     call distribute(0,Nl-1,alpha,alpha_g)

     call distribute(0,Nl-1,cprocaPhi_R,procaPhi_g)
     call distribute(0,Nl-1,cprocaA_I,procaA_g)
     call distribute(0,Nl-1,cprocaB_I,procaB_g)
     call distribute(0,Nl-1,cprocaG_I,procaG_g)

     call distribute(0,Nl-1,cprocaE_R,procaE_g)
     call distribute(0,Nl-1,cprocaXi_R,procaXi_g)

  end if


! ********************************************************
! ***   IMAGINARY PARTS OF (procaPhi,procaE,procaXi)   ***
! ***   AND REAL PARTS  OF (procaA,procaB,procaG)      ***
! ********************************************************

! Set imaginary parts of (procaPhi,procaE,procaXi) to zero.

  cprocaPhi_I = 0.d0
  cprocaE_I   = 0.d0
  cprocaXi_I  = 0.d0

! Set real part of (procaA,procaB,procaG) to zero.

  cprocaA_R = 0.d0
  cprocaB_R = 0.d0
  cprocaG_R = 0.d0

! Find procaL.

  cprocaL_R = 0.d0
  cprocaL_I= (cprocaA_I - cprocaG_I)/r**2


! ***************
! ***   END   ***
! ***************

  end subroutine idata_LProcastar








! **********************************
! ***   RADIAL DERIVATIVE OF A   ***
! **********************************

! The radial derivative of A comes from the
! Hamiltonian constraint.

  function J1_lproca(A,alpha,procaFres,procaAres,procaBres,procaGres,rm)

  use param

  implicit none

  real(8) J1_lproca
  real(8) A,alpha,procaFres,procaAres,procaBres,procaGres,rm
  real(8) procaF,procaA,procaB,procaG,procaE
  real(8) rho

! Non-rescaled variables.

  procaF = procaFres*rm**cproca_l
  procaA = procaAres*rm**(cproca_l+1)
  procaB = procaBres*rm**(cproca_l+2)
  procaG = procaGres*rm**(cproca_l+1)

! Radial electric field:
!
! E  = - alpha [ m**2 a + l(l+1)(a-G)/r**2] / (A omega)

  procaE = - alpha/(A*proca_omega)*(cproca_mass**2*procaA &
         + dble(cproca_l)*(dble(cproca_l)+1.d0)*(procaA-procaG)/rm**2)

! Energy density:
!                        2     2           2    2                 2 2                        2        2
! rho  = + 1/(8 pi) { A E  +  m [ (F/alpha)  + a / A ] + l(l+1)[ m b  + ([omega b + F]/alpha)  + (a-G) / A ]/ r**2 }
!
! Notice that we don't divide by 8*pi since it cancels.

  rho = A*procaE**2 + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A) &
      + dble(cproca_l)*(dble(cproca_l)+1.d0)*((cproca_mass*procaB)**2 &
      + ((proca_omega*procaB + procaF)/alpha)**2 + (procaA-procaG)**2/A)/rm**2

! Radial derivative of A:
!
! dA/dr  =  A [ (1 - A)/r + 8 pi r A rho ]
!
! Notice that we don't multiply the last term with 8*pi since it cancels.

  J1_lproca = A*((1.d0-A)/rm + rm*A*rho)

  end function J1_lproca








! **************************************
! ***   RADIAL DERIVATIVE OF ALPHA   ***
! **************************************

! The radial derivative of alpha comes from
! the polar slicing condition:  dK_{theta,theta}/dt=0.

  function J2_lproca(A,alpha,procaFres,procaAres,procaBres,procaGres,rm)

  use param

  implicit none

  real(8) J2_lproca
  real(8) A,alpha,procaFres,procaAres,procaBres,procaGres,rm
  real(8) procaF,procaA,procaB,procaG,procaE
  real(8) SA

! Non-rescaled variables.

  procaF = procaFres*rm**cproca_l
  procaA = procaAres*rm**(cproca_l+1)
  procaB = procaBres*rm**(cproca_l+2)
  procaG = procaGres*rm**(cproca_l+1)

! Radial electric field:
!
! E  = - alpha [ m**2 a + l(l+1)(a-G)/r**2] / (A omega)

  procaE = - alpha/(A*proca_omega)*(cproca_mass**2*procaA &
         + dble(cproca_l)*(dble(cproca_l)+1.d0)*(procaA-procaG)/rm**2)

! Radial stress tensor component:
!
!                         2     2           2    2                   2 2                        2        2
! SA  = + 1/(8 pi) { - A E  +  m [ (F/alpha)  + a / A ] + l(l+1)[ - m b  + ([omega b + F]/alpha)  + (a-G) / A ]/ r**2 }
!
! Notice that we don't divide by 8*pi since it cancels.

  SA = - A*procaE**2 + cproca_mass**2*((procaF/alpha)**2 + procaA**2/A) &
     + dble(cproca_l)*(dble(cproca_l)+1.d0)*(((proca_omega*procaB + procaF)/alpha)**2 &
     - cproca_mass**2*procaB**2 + (procaA-procaG)**2/A)/rm**2

! Radial derivative of alpha:
!
! dalpha/dr  =  alpha [ (A - 1)/(2r) + 4 pi r A SA ]
!
! Notice that we don't multiply the last term with 4*pi*A since
! it cancels (we are left only with 1/2).

  J2_lproca = alpha*(0.5d0*(A-1.d0)/rm + 0.5d0*rm*A*SA)

  end function J2_lproca








! ******************************************
! ***   RADIAL DERIVATIVE OF procaFres   ***
! ******************************************

! The radial derivative of procaFres from the
! Proca evolution equations.

  function J3_lproca(A,alpha,procaFres,procaAres,procaBres,procaGres,rm)

  use param

  implicit none

  real(8) J3_lproca
  real(8) A,alpha,procaFres,procaAres,procaBres,procaGres,rm

! Radial derivative of procaFres:
!
! dF_res/dr  =  ( -l F_res + l(l+1)(alpha**2 /omega)[ a_res - G_res ] ) / r  +  r omega a_res [ (m alpha/omega)**2 - 1 ]
!

  J3_lproca = (-dble(cproca_l)*procaFres & 
                + dble(cproca_l)*(dble(cproca_l)+1.d0)*(alpha**2/proca_omega)*(procaAres-procaGres))/rm &  
                + rm*proca_omega*procaAres*((cproca_mass*alpha/proca_omega)**2 - 1.d0)

  end function J3_lproca








! ******************************************
! ***   RADIAL DERIVATIVE OF procaAres   ***
! ******************************************

! The radial derivative of procaAres from the Proca
! evolution equations.

  function J4_lproca(A,alpha,procaFres,procaAres,procaBres,procaGres,rm)

  use param

  implicit none

  real(8) J4_lproca
  real(8) A,alpha,procaFres,procaAres,procaBres,procaGres,rm
  real(8) procaF,procaA,procaB,procaG,procaE
  real(8) aux

! Non-rescaled variables.

  procaF = procaFres*rm**cproca_l
  procaA = procaAres*rm**(cproca_l+1)
  procaB = procaBres*rm**(cproca_l+2)
  procaG = procaGres*rm**(cproca_l+1) 

! Radial electric field:
!
! E  = - alpha [ m**2 a + l(l+1)(a-G)/r**2] / (A omega)

  procaE = - alpha*(cproca_mass**2*procaA &
         + dble(cproca_l)*(dble(cproca_l)+1.d0)*(procaA-procaG)/rm**2)/(A*proca_omega)

!                            2           2 2   2
! SA - rho  = - (1/4pi) { A E  + l(l+1) m b / r }
!
! Notice that we don't divide by 4*pi since it cancels.

  aux = - A*procaE**2 - dble(cproca_l)*(dble(cproca_l)+1.d0)*cproca_mass**2*procaB**2/rm**2

!
! da_res/dr  = [ (omega A/alpha**2) F_res - (l+2+A) a_res  + l(l+1) A b_res ]/r - 4pi r a_res  A (SA-rho)
!
! Notice that we don't multiply the last term with 4*pi since
! it cancels.

  J4_lproca = (proca_omega*procaFres*A/alpha**2 - (dble(cproca_l)+2.d0+A)*procaAres &
            + dble(cproca_l)*(dble(cproca_l)+1.d0)*A*procaBres)/rm - rm*procaAres*A*aux

  end function J4_lproca







! ******************************************
! ***   RADIAL DERIVATIVE OF procaBres   ***
! ******************************************

! The radial derivative of procaBres. This is just the definition of G_res.

  function J5_lproca(procaBres,procaGres,rm)

  use param

  implicit none

  real(8) J5_lproca
  real(8) procaBres,procaGres,rm

! d b_res/dr = (G_res - (l+2) b_res)/r

  J5_lproca = (procaGres - (dble(cproca_l)+2.d0)*procaBres)/rm

  end function J5_lproca






! ******************************************
! ***   RADIAL DERIVATIVE OF procaGres   ***
! ******************************************

! The radial derivative of procaGres from the Proca
! evolution equations.

  function J6_lproca(A,alpha,procaFres,procaAres,procaBres,procaGres,rm)

  use param

  implicit none

  real(8) J6_lproca
  real(8) A,alpha,procaFres,procaAres,procaBres,procaGres,rm
  real(8) procaF,procaA,procaB,procaG,procaE
  real(8) aux

! Non-rescaled variables.

  procaF = procaFres*rm**cproca_l
  procaA = procaAres*rm**(cproca_l+1)
  procaB = procaBres*rm**(cproca_l+2)
  procaG = procaGres*rm**(cproca_l+1)

! Radial electric field:
!
! E  = - alpha [ m**2 a + l(l+1)(a-G)/r**2] / (A omega)

  procaE = - alpha*(cproca_mass**2*procaA &
         + dble(cproca_l)*(dble(cproca_l)+1.d0)*(procaA-procaG)/rm**2)/(A*proca_omega)
 
!                            2           2 2   2
! SA - rho  = - (1/4pi) { A E  + l(l+1) m b / r }
!
! Notice that we don't divide by 4*pi since it cancels.

  aux = - A*procaE**2 - dble(cproca_l)*(dble(cproca_l)+1.d0)*cproca_mass**2*procaB**2/rm**2

!
! dG_res/dr  =  [ -2 a_res - (l+A) G_res + l(l+1) A b_res ] / r
!
!                       2               2
!               + r[ ( m - (omega/alpha) ) A b_res - 4pi A G_res (SA-rho) ]
!
! Notice that we don't multiply the last term with 4*pi since
! it cancels.

  J6_lproca = (- 2*procaAres - (dble(cproca_l)+A)*procaGres + dble(cproca_l)*(dble(cproca_l)+1.d0)*A*procaBres)/rm  &
            + rm*A*((cproca_mass**2-(proca_omega/alpha)**2)*procaBres - procaGres*aux)

  end function J6_lproca

