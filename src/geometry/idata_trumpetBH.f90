!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/idata_trumpetBH.f90,v 1.23 2022/07/12 23:58:58 malcubi Exp $

  subroutine idata_trumpetBH

! *******************************************
! ***   TRUMPET BLACK HOLE INITIAL DATA   ***
! *******************************************

! This program solves for "trumpet" black hole
! initial data.
!
! Then original equation for the conformal factor
! has the form:
!
!  2                          8          7  6
! d psi  +  2 d psi / 2  +  k1  / ( 4 psi  r  )  =  0
!  r           r
!
! with k1 some constant. Notice that this equation has
! an exact solution of the form:
!
! psi = k1/sqrt(r)
!
! Unfortunately, this solution does not satisfy the
! asymptotic boundary condition, which requires that
! psi goes to 1 far away.
!
! In order to solve the equation with the correct asymptotic
! condition, we first define:
!
! u  :=  sqrt(r) psi
!
! The equation to solve now becomes:
!
!  2                        2      8   7
! d u  +  d u / r  +  ( 1/4r ) [ k1 / u  -  u ]  =  0
!  r       r
!
! Notice that u=k1 solves the equation.
!
! One could now try to solve this equation integrating
! outwards from r=0, and using the boundary conditions
! u(0) = k1, u'(0) = 0. Unfortunately, one then only recovers
! the exact solution above, which does not satisfy the asymptotic
! condition.
!
! One can instead ask that, for small r, u behaves as:
!
! u ~ k1 + k2 r^n
!
! When we substitute this in the equation we find that the
! constant k2 is arbitrary, and the power "n" must be equal
! to sqrt(2).
!
! This has the disadvantage that for small r the second derivative
! of u diverges (u is only C1 at r=0). But we can use the expansion
! above to start the integration at some small value of r.  Here
! we set initial data at r=dr/2.
!
! Notice also that the differential equation for psi has the
! property that if psi is a solution for k1=1, then psi/c
! will be a solution for k1=c.  So in fact we only need
! to solve the equation with k1=1, and then rescale the
! solution in order to ensure that psi goes to 1 at infinity.
!
! Tke constant k2, on the other hand, remains free, and it
! can be adjusted in a shooting method to obtain a desired
! value for the black hole mass using the fact that far away:
!
! psi  ~  1 + M/2r
!
! Finally, for the full trumpet solution we must fix the
! lapse, shift and extrinsic curvature to:
!
! alpha  =  1 + 2 r d psi / psi  =  2 r d u / u
!                    r                   r
!
! beta   =  Cb / (r**2  psi**6)
!
!
! KTA    =  Ck / (r**3  psi**6)
!
! with Cb and Ck constants given by:
!
! Cb  =  + 1/sqrt(3) k1**4
!
! Ck  =  - 2 Cb
!
! Notice also that this initial data is conformally flat (A=B=1),
! and maximal (trK=0).

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays
  use derivatives

! Declare variables.

  implicit none

  integer i,l,iter               ! Counters.
  integer imin                   ! Leftmost grid point.
  integer iaux                   ! Auxiliary quantity.
  integer :: maxiter = 100       ! Maximum number of iterations.

  real(8) hdr                    ! dr/2.
  real(8) kk1,kk2,kk2_old        ! Coefficients of u expansion for small r:  u~kk1+kk2*sqrt(r)
  real(8) u_rk,su1,su2,su3,su4   ! Runge-Kutta coefficients for u.
  real(8) v_rk,sv1,sv2,sv3,sv4   ! Runge-Kutta coefficients for v.
  real(8) psifar                 ! Asymptotic value of psi.
  real(8) mass                   ! Black hole mass.
  real(8) res,res_old            ! Residual.
  real(8) Cb,Ck                  ! Constants for expressions for beta and KA.
  real(8) rm,delta,half,aux      ! Auxiliary quantities.
  real(8) :: epsilon = 1.d-10    ! Tolerance.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr      ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: u       ! Solution of equation.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: v       ! du/dr.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: psi_g   ! Conformal factor.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: alpha_g ! Lapse.


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for a trumpet black hole ...'
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


! *****************************
! ***   START INTEGRATION   ***
! *****************************

! Give initial values of (k1,k2).

  kk1 = 1.0d0
  kk2 = 0.1d0

! Initialize u and v.

  u = kk1
  v = 0.d0


! *********************************************
! ***   ONLY PROCESSOR 0 SOLVES THE ODE's   ***
! *********************************************

  if (rank==0) then


!    ****************************
!    ***   BEGIN ITERATIONS   ***
!    ****************************

!    Initialize residual.

     res = 1.d0

!    Iterations for shooting method.

     iter = 0

     do while ((abs(res).gt.epsilon).and.(iter.lt.maxiter))

        iter = iter + 1


!       *********************************
!       ***   LOOP OVER GRID LEVELS   ***
!       *********************************

!       We solve from fine to coarse grid.

        do l=Nl-1,0,-1

!          Find initial point. Only the finest grid
!          integrates from the origin.

           if (l==Nl-1) then

              imin = 2

!             First grid point.  For this we use the fact that
!             for small r the solution goes as:
!
!             u  ~  k1 + k2 r^sqrt(2)
!
!             so that:
!
!             v  ~ k2 sqrt(2) r^(sqrt(2)-1)
!
!             with "a" a constant that will be used to find the
!             correct asymptotic behaviour.
!
!             Also, remember that the first grid point is at r=dr/2.

              aux = dsqrt(2.d0)
              hdr = half*dr(l)

              u(l,1) = kk1 + kk2*hdr**aux
              v(l,1) = kk2*aux*hdr**(aux-1.d0)

           else

              imin = Nrtotal/2

!             For coarse grids we interpolate the initial point
!             (which is half way through the grid) from the last
!             points of the higher resolution grid.  Here we use
!             cubic interpolation.

              u(l,imin-1) = (9.d0*(u(l+1,Nrtotal-2)+u(l+1,Nrtotal-3)) &
                                - (u(l+1,Nrtotal-4)+u(l+1,Nrtotal-1)))/16.d0

              v(l,imin-1) = (9.d0*(v(l+1,Nrtotal-2)+v(l+1,Nrtotal-3)) &
                                - (v(l+1,Nrtotal-4)+v(l+1,Nrtotal-1)))/16.d0

           end if


!          **********************************
!          ***   RUNGE-KUTTA INTEGRATION  ***
!          **********************************

           do i=imin,Nrtotal

!             First Runge-Kutta step.

              su1 = v(l,i-1)
              sv1 = - v(l,i-1)/rr(l,i-1) + 0.25d0*(u(l,i-1) - kk1**8/u(l,i-1)**7)/rr(l,i-1)**2

              u_rk = u(l,i-1) + half*dr(l)*su1
              v_rk = v(l,i-1) + half*dr(l)*sv1

!             Second Runge-Kutta step.

              rm = rr(l,i-1) + half*dr(l)

              su2 = v_rk
              sv2 = - v_rk/rm + 0.25d0*(u_rk - kk1**8/u_rk**7)/rm**2

              u_rk = u(l,i-1) + half*dr(l)*su2
              v_rk = v(l,i-1) + half*dr(l)*sv2

!             Third Runge-Kutta step.

              su3 = v_rk
              sv3 = - v_rk/rm + 0.25d0*(u_rk - kk1**8/u_rk**7)/rm**2

              u_rk = u(l,i-1) + dr(l)*su3
              v_rk = v(l,i-1) + dr(l)*sv3

!             Fourth Runge-Kutta step.

              rm = rr(l,i)

              su4 = v_rk
              sv4 = - v_rk/rm + 0.25d0*(u_rk - kk1**8/u_rk**7)/rm**2

!             Advance variables.

              !u(l,i) = u(l,i-1) + dr(l)*su2     ! Second order Runge-Kutta
              !v(l,i) = v(l,i-1) + dr(l)*sv2     ! Second order Runge-Kutta

              u(l,i) = u(l,i-1) + dr(l)*(su1 + 2.d0*(su2 + su3) + su4)/6.d0
              v(l,i) = v(l,i-1) + dr(l)*(sv1 + 2.d0*(sv2 + sv3) + sv4)/6.d0

           end do

!          Ghost zones. Notice that u inherits the parity from psi.
!          On the other hand, v does not have well-defined parity
!          since it goes as r**(sqrt(2)-1).  This does not matter
!          too much as it is never used, but giving it positive
!          parity makes for nicer plots!

           do i=1,ghost
              u(l,1-i) = +u(l,i)
              v(l,1-i) = +v(l,i)
           end do

        end do


!       ************************************************
!       ***   FIND CONFORMAL FACTOR AND RESCALE IT   ***
!       ************************************************

!       Conformal factor.

        psi_g = u/sqrt(abs(rr))

!       Find asymptotic value of psi. We obtain this by assuming
!       that far away we have:  psi ~ psifar + const/r, which
!       implies that:
!
!       psifar  =  psi + r dpsi/dr
!
!       We calculate dalpha/dr at the boundary with one sided
!       differences and solve for alpha0.

        psifar = 1.d0

        if (order=="two") then
           psifar = psi_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
               *(3.d0*psi_g(0,Nrtotal) - 4.d0*psi_g(0,Nrtotal-1) + psi_g(0,Nrtotal-2))
        else
           psifar = psi_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
               *(25.d0*psi_g(0,Nrtotal) - 48.d0*psi_g(0,Nrtotal-1) &
               + 36.d0*psi_g(0,Nrtotal-2) - 16.d0*psi_g(0,Nrtotal-3) + 3.d0*psi_g(0,Nrtotal-4))/3.d0
        end if

!       Rescale phi.

        psi_g = psi_g/psifar


!       **********************
!       ***   FIND LAPSE   ***
!       **********************

!       The lapse is given by:  alpha = 2 r v / u.

        alpha_g = 2.d0*rr*v/u


!       *************************************************************
!       ***   FIND ADM MASS, RESIDUAL, AND CORRECT k2 IF NEEDED   ***
!       *************************************************************

!       Remember that for Schwarzschild  asymptotically we
!       must have psi ~ 1 + M/2r, so that the mass is:
!
!       M ~ 2 r (psi-1)
!
!       In order to find the correct mass, we use Newton's method.

        mass = 2.d0*rr(0,Nrtotal)*(psi_g(0,Nrtotal) - 1.d0)

        res_old = res
        res = BHmass - mass

        if (iter==1) then
           delta = 0.1d0
        else
           delta = - res*(kk2-kk2_old)/(res-res_old)
        end if

        kk2_old = kk2
        kk2 = kk2 + delta


!       **************************
!       ***   END ITERATIONS   ***
!       **************************

!       Output data to screen.

        write(*,"(A,I4,A,ES10.4,A,ES9.2)") ' Iteration: ',iter,'    Mass: ',mass,'    Residual: ',res

     end do

!    If we succeded rescale k1 with psifar.

     if (iter<maxiter) then

!       Rescale k1.

        kk1 = kk1/psifar

!       Message to screen.

        print *
        print *, 'Found solution!'
        print *


!    If we failed abort.

     else

        print *
        print *, 'Maximum iteration number reached in idata_trumpetBH.f90.'
        print *, 'Initial data solver did not converge.'
        print *, 'Aborting! (subroutine idata_trumpetBH.f90)'
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

!       Loop over grid.

        do i=1,Nrtotal-ghost,2

           iaux = i/2 + 1

           rm = rr(l-1,iaux)

           u(l-1,iaux) = (9.d0*(u(l,i)+u(l,i+1)) - (u(l,i-1)+u(l,i+2)))/16.d0
           v(l-1,iaux) = (9.d0*(v(l,i)+v(l,i+1)) - (v(l,i-1)+v(l,i+2)))/16.d0

           psi_g(l-1,iaux) = (9.d0*(psi_g(l,i)+psi_g(l,i+1)) - (psi_g(l,i-1)+psi_g(l,i+2)))/16.d0
           alpha_g(l-1,iaux) = (9.d0*(alpha_g(l,i)+alpha_g(l,i+1)) - (alpha_g(l,i-1)+alpha_g(l,i+2)))/16.d0

        end do

!       Fix ghost zones.

        do i=1,ghost
           u(l-1,1-i) = +u(l-1,i)
           v(l-1,1-i) = +v(l-1,i)
           psi_g(l-1,1-i)   = + psi_g(l-1,i)
           alpha_g(l-1,1-i) = + alpha_g(l-1,i)
        end do

     end do


! *************************************
! ***   FINISHED FINDING SOLUTION   ***
! *************************************

! When we get here either the solution has been found,
! so we close the "if" statement for processor 0.

  end if


! ************************************************
! ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
! ************************************************

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! arrays with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

! Single processor run.

  if (size==1) then

!    Conformal factor and lapse.

     psi = psi_g
     alpha = alpha_g

!    Variables for output.

     if (allocated(trumpet_u)) then
        trumpet_u = u
     end if

     if (allocated(trumpet_v)) then
        trumpet_v = v
     end if

! Parallel run.

  else

!    Conformal factor and lapse.

     call distribute(0,Nl-1,psi,psi_g)
     call distribute(0,Nl-1,alpha,alpha_g)

!    Variables for output.

     if (allocated(trumpet_u)) then
        call distribute(0,Nl-1,trumpet_u,u)
     end if

     if (allocated(trumpet_v)) then
        call distribute(0,Nl-1,trumpet_v,v)
     end if

  end if


! *****************************************************
! ***   RECONSTRUCT SHIFT AND EXTRINSIC CURVATURE   ***
! *****************************************************

! Shift.

  Cb = 1.d0/dsqrt(3.d0)*kk1**4

  beta = Cb/r**2/psi**6

! Extrinsic curvature.

  Ck = - 2.d0*Cb

  KTA = Ck/abs(r)**3/psi**6
  KTB = - 0.5d0*KTA

  Klambda = (KTA - KTB)/r**2

  if (regular2) then
     Klambda2 = Klambda/dexp(4.0d0*phi)
  end if

! Quantities related to the conformal factor.

  psi2 = psi**2
  psi4 = psi**4

  phi = dlog(psi)
  chi = 1.d0/psi**chipower

! Solve maximal slicing for lapse. This is not really
! necessary so it is commented out.  If you use it beware:
! It must be called with the Dirichlet boundary condition
! since the Robin condition turns out to be incompatible
! with the static initial data.

! call alphamaximal("conformal",1.d0)


! ***************
! ***   END   ***
! ***************

  end subroutine idata_trumpetBH


