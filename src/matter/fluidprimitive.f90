
  subroutine fluidprimitive(l)

! *********************************************
! ***   RECOVER FLUID PRIMITIVE VARIABLES   ***
! *********************************************

! This routine recovers the fluid primitive variables (rho0,p,v)
! starting from the conserved quantities (D,E,S).
!
! Since the system of equations we need to invert is coupled and
! non-linear, we do the inversion using Newton-Raphson's method.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  logical :: firstcall

  integer :: i,l

  real(8) :: aux

  data firstcall / .true. /


! ***************************************
! ***   RECOVER PRIMITIVE VARIABLES   ***
! ***************************************

! On firstcall rescale fluid_atmos with maximum of fluid_rho.

  if (firstcall) then
     fluid_atmos = maxval(fluid_rho)*fluid_atmos
     call MPI_Allreduce(fluid_atmos,aux,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
     fluid_atmos = aux
     firstcall = .false.
  end if

! For initial data we don't need to recover the primitive
! variables, they should be fixed by the initial data routines. 
! We also don't add the atmosphere at t=0.

  if (t(l)>0.d0) then
     call recover(l)
  end if


! ************************
! ***   GHOST POINTS   ***
! ************************

  if (rank==0) then
     do i=1,ghost

        fluid_rho(l,1-i) = + fluid_rho(l,i)
        fluid_p(l,1-i)   = + fluid_p(l,i)
        fluid_e(l,1-i)   = + fluid_e(l,i)
        fluid_csi(l,1-i) = + fluid_csi(l,i)

        fluid_u(l,1-i)   = - fluid_u(l,i)
        fluid_v(l,1-i)   = - fluid_v(l,i)
        fluid_W(l,1-i)   = + fluid_W(l,i)

     end do
  end if


! ********************
! ***   ENTHALPY   ***
! ********************

! h  =  1 + e + (p+q)/rho
!
! Avoid divisions by zero!

  do i=1-ghost,Nr
     if (fluid_rho(l,i)==0.d0) then
        fluid_h(l,i) = 0.d0
     else
        fluid_h(l,i) = 1.d0 + fluid_e(l,i) + (fluid_p(l,i) + fluid_q(l,i))/fluid_rho(l,i)
     end if
  end do


! **************************
! ***   SPEED OF SOUND   ***
! **************************

  if (fluid_EOS/="none") then

!    Ideal gas equation of state.

     if (fluid_EOS=="ideal") then

!       Speed of sound for an ideal gas equation of state:
!
!       vs^2  = gamma (gamma - 1) e / ( 1 + gamma e)
!
!             =  p gamma (gamma - 1) / (p gamma + rho0 (gamma - 1))
!
!       Avoid divisions by zero!

        do i=1-ghost,Nr
           if (fluid_rho(l,i)==0.d0) then
              fluid_vs(l,i) = 0.d0
           else
              fluid_vs(l,i) = sqrt(abs(fluid_p(l,i)*fluid_gamma*(fluid_gamma-1.d0) &
                            /(fluid_p(l,i)*fluid_gamma + fluid_rho(l,i)*(fluid_gamma-1.d0))))
          end if
        end do

!    ADD IF STATEMENTS FOR NEW EQUATIONS OF STATE HERE.

!    Unknown equation of state.

     else

        print *
        print *, 'Unknown equation of state'
        print *, 'Aborting! (subroutine fluidprimitive, speed of sound)'
        print *
        call die

     end if

! If we don't have an equation of state we just use the polytropic relation
! in order to avoid problems in other routines if vs is not defined.  But
! really this is not a good idea.

  else

     do i=1-ghost,Nr
        if (fluid_rho(l,i)==0.d0) then
           fluid_vs(l,i) = 0.d0
        else
           fluid_vs(l,i) = sqrt(abs(fluid_gamma*fluid_p(l,i)/(fluid_rho(l,i)*fluid_h(l,i))))
        end if
     end do

  end if

! Calculate Mach number.

  fluid_Mach(l,:) = fluid_v(l,:)/fluid_vs(l,:)


! *********************************
! ***   CHARACTERISTIC SPEEDS   ***
! *********************************

! Characteristic speeds (in spherical symmetry, see page 262 of my book):
!
! v_+  =  [ alpha/(1 - A psi**4 fluid_v**2 fluid_vs**2) ] [ fluid_v ( 1 - fluid_vs**2) + aux ]
!
! v_-  =  [ alpha/(1 - A psi**4 fluid_v**2 fluid_vs**2) ] [ fluid_v ( 1 - fluid_vs**2) - aux ]
!
! with:
!
! aux  =  fluid_vs (1 - A psi**4 fluid_v**2) / ( A^(1/2) psi**2 )

  do i=1-ghost,Nr

     aux = fluid_vs(l,i)*abs(1.d0 - A(l,i)*exp(4.d0*phi(l,i))*fluid_v(l,i)**2) &
         /sqrt(A(l,i))/exp(2.d0*phi(l,i))

     fluid_vcp(l,i) = alpha(l,i)/(1.d0 - A(l,i)*exp(4.d0*phi(l,i))*fluid_v(l,i)**2*fluid_vs(l,i)**2) &
                    *(fluid_v(l,i)*(1.d0 - fluid_vs(l,i)**2) + aux)
     fluid_vcm(l,i) = alpha(l,i)/(1.d0 - A(l,i)*exp(4.d0*phi(l,i))*fluid_v(l,i)**2*fluid_vs(l,i)**2) &
                    *(fluid_v(l,i)*(1.d0 - fluid_vs(l,i)**2) - aux)

     if (shift/="none") then
        fluid_vcp(l,i) = - beta(l,i) + fluid_vcp(l,i)
        fluid_vcm(l,i) = - beta(l,i) + fluid_vcm(l,i)
     end if

  end do


! ********************************
! ***   ARTIFICIAL VISCOSITY   ***
! ********************************

! Update artificial viscosity:
!
! 1) If the divergence of S is positive then set
! fluid_q to zero.
!
! 2) If the divergence is negative compute the artificial
! viscosity contribution. The artificial viscosity has
! two contributions and has the following general form:
!
! q  =  dr abs(div.S) ( q1 vs + q2 dr abs(div.S) / rho )
!
! with q1 and q2 positive constants, and vs the speed
! of sound. Notice that the term with q1 introduces a first
! order error, so one must keep q1 small.  The term with q2
! introduces a second order error. For rho we use the ADM
! energy density: rho = E + D.

  if ((fluid_q1/=0.d0).or.(fluid_q2/=0.d0)) then

     do i=1-ghost,Nr

        aux = dr(l)*D1_fluid_cS(l,i)/(A(l,i)*exp(4.d0*phi(l,i)))

        if (aux>=0.d0) then
           fluid_q(l,i) = 0.d0
        else
           fluid_q(l,i) = abs(aux)*(fluid_q1*fluid_vs(l,i) &
                        + fluid_q2*abs(aux)/(fluid_cE(l,i) + fluid_cD(l,i)))
        end if

     end do

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine fluidprimitive










  subroutine recover(l)

! Here we recover primitive variables using different methods.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  integer :: i,j,l
  integer :: maxiter = 1000

  real(8) :: rhoatmos,Eatmos,patmos
  real(8) :: W,W2,res,aux
  real(8) :: p1,p2,f1,f2
  real(8) :: csi1,csi2,S2

  real(8) :: epsilon = 1.d-10 ! Tolerance
  real(8) :: Wmax = 1.d5      ! Maximum allowed Lorentz factor


! ********************************
! ***   FIX ATMOSPHERE SCALE   ***
! ********************************

! Find atmosphere values.  We set both "rho0" and "e" to a small value,
! and obtain the pressure from the polytropic equation of state.
!
! Notice that we can't just set e=p=0, since in that case the speed of
! sound vanishes and the routine that calculates the fluid sources fails.

  rhoatmos = fluid_atmos/2.d0

  patmos = fluid_kappa*rhoatmos**fluid_gamma
  Eatmos = fluid_kappa*rhoatmos**(fluid_gamma-1.d0)/(fluid_gamma-1.d0)


! *********************************
! ***   LOOP OVER GRID POINTS   ***
! *********************************

  do i=1,Nr


!    **********************
!    ***   ATMOSPHERE   ***
!    **********************

!    Atmosphere: If the density is too low we run the risk of dividing
!    by a very small quantity in order to recover the fluid speed,
!    which will introduce very large round-off errors and can cause
!    the code to crash.
!
!    In order to avoid this, we set the conserved mass density to a small
!    value and the momentum density to zero. For the conserved energy
!    density we take E = rho0*e.

     if (fluid_cD(l,i)<=fluid_atmos) then

!       Conserved quantities (D,E,S).

        fluid_cD(l,i) = rhoatmos
        fluid_cE(l,i) = rhoatmos*Eatmos
        fluid_cS(l,i) = 0.d0

!       Fluid speed and Lorentz factor.

        fluid_v(l,i) = 0.d0
        fluid_u(l,i) = 0.d0
        fluid_W(l,i) = 1.d0

!       Primitive variables.

        fluid_rho(l,i) = rhoatmos
        fluid_p(l,i) = patmos
        fluid_e(l,i) = Eatmos

!       Enthalpy and csi.

        fluid_h(l,i) = 1.d0 + fluid_e(l,i) + fluid_p(l,i)/fluid_rho(l,i)
        fluid_csi(l,i) = fluid_rho(l,i)*fluid_h(l,i)*fluid_W(l,i)

        goto 10

     end if


!    *************************
!    ***   DIRECT METHOD   ***
!    *************************

!    The idea of the "direct" method is to choose a trial value of
!    the pressure (the value at the old time level).  From this we
!    calculate first the fluid velocity from:
!
!     r                             4 phi
!    v   =  S  / (E + D + p) / ( A e     )
!            r
!
!    where the geometric factor is there in order to raise the index
!    (v has index up, while S has index down), and where "q" is the
!    contribution to the pressure coming from the artificial viscosity.
!
!    We then calculate the Lorentz factor W from:
!
!                             4 phi  2
!    W  =  1 / sqrt( 1  -  A e      v  )
!
!
!    Notice that we must always have:  A exp(4 phi) v^2 < 1.
!
!    The rest-mass energy density is then obtained as:
!
!
!    rho  =  D / W
!
!
!    and the specific internal energy is:
!
!                                          2
!    e  =  ( E + D (1 - W) + (p + q) (1 - W ) ) / D W
!
!
!    The above expression can be found from the definitions of
!    E and the enthalpy h.
!
!    Having found these variables, we find the difference between
!    the trial value of the pressure and the one obtained from
!    the equation of state, and iterate until this difference
!    is smaller than a tolerance set by "epsilon".

     if (fluid_primitive=="direct") then

!       Initialize counter and residue.

        j = 0
        res = 1.d0

!       Start iterations for Newton-Raphson's method.

        do while ((abs(res)>epsilon).and.(j<maxiter))

!          Increment counter.

           j = j + 1

!          Set initial guess for pressure to old value.

           if (j==1) p1 = fluid_p(l,i)

!          Find value of v using old value of pressure.

           if (abs(fluid_cS(l,i))==0.d0) then
              fluid_v(l,i) = 0.d0
           else
              fluid_v(l,i) = fluid_cS(l,i)/(fluid_cE(l,i) + fluid_cD(l,i) + p1 + fluid_q(l,i)) &
                           /(A(l,i)*exp(4.d0*phi(l,i)))
           end if

!          Find Lorentz factor W, and check that v^2 is less than 1.

           aux = 1.d0 - A(l,i)*exp(4.d0*phi(l,i))*fluid_v(l,i)**2

           if (aux<=0.d0) then
              W = Wmax
              if (fluid_cS(l,i)>0.d0) then
                 fluid_v(l,i) = + sqrt((1.d0 - 1.d0/W**2))/sqrt(abs(A(l,i)))/exp(2.d0*phi(l,i))
              else
                 fluid_v(l,i) = - sqrt((1.d0 - 1.d0/W**2))/sqrt(abs(A(l,i)))/exp(2.d0*phi(l,i))
              end if
           else
              W = 1.d0/sqrt(abs(aux))
           end if

!          Find rest-mass density: rho = D/W.

           fluid_rho(l,i) = fluid_cD(l,i)/W

!          Find specific internal energy.

           aux = fluid_cD(l,i)*(W - 1.d0) + (p1 + fluid_q(l,i))*(W**2 - 1.d0)
           fluid_e(l,i) = (fluid_cE(l,i) - aux)/(fluid_cD(l,i)*W)

!          Calculate difference between trial value of the pressure
!          and its predicted value from the equation of state, and
!          update the value of the pressure for the next iteration.
!
!          If there is NO equation of state (for example if you forced
!          some specific form of rho for the initial data and solved
!          the TOV equations), then we just set the difference to 0.

!          Ideal gas equation of state: p = (1-gamma) rho e

           if (fluid_EOS=="ideal") then

              f1 = (fluid_gamma-1.d0)*fluid_rho(l,i)*fluid_e(l,i) - p1

              if (j==1) then
                 aux = (fluid_gamma-1.d0)*fluid_rho(l,i)*fluid_e(l,i)
              else
                 aux = p1 - f1*(p2-p1)/(f2-f1)
              end if

              p2 = p1
              f2 = f1

              p1 = aux
              if (p1<=patmos) p1=patmos

!             Calculate residual.

              res = abs(p2 - p1)

!          ADD IF STATEMENTS FOR NEW EQUATIONS OF STATE HERE.

           else

              print *
              print *, 'Unknown equation of state'
              print *, 'Aborting! (subroutine fluidprimitive)'
              print *
              call die

           end if

        end do

!       If we exceeded the maximum number of iterations
!       send message to screen.

        if (j==maxiter) then
           print *
           print *, 'Maximum iteration number reached in fluidprimitive.f90 at point: i = ',i,', r = ',r(l,i)
           print *, 'Residue = ',real(res)
           print *
           call die
        end if

!       Save last value of p.

        fluid_p(l,i) = p1

!       Recalculate fluid speed and Lorentz number using new value of pressure.

        fluid_v(l,i) = fluid_cS(l,i)/(fluid_cE(l,i) + fluid_cD(l,i) + fluid_p(l,i) + fluid_q(l,i)) &
                     /(A(l,i)*exp(4.d0*phi(l,i)))

        aux = 1.d0 - A(l,i)*exp(4.d0*phi(l,i))*fluid_v(l,i)**2

        if (aux<=0.d0) then
           W = Wmax
           print *, "Warning: Fluid speed larger than light speed at point ",i," radius",r(l,i)," time",t(l)
           print *, aux,fluid_v(l,i),fluid_cS(l,i),fluid_cD(l,i),fluid_cE(l,i)
           print *
           if (fluid_cS(l,i)>0.d0) then
              fluid_v(l,i) = + sqrt((1.d0 - 1.d0/W**2))/sqrt(abs(A(l,i)))/exp(2.d0*phi(l,i))
           else
              fluid_v(l,i) = - sqrt((1.d0 - 1.d0/W**2))/sqrt(abs(A(l,i)))/exp(2.d0*phi(l,i))
           end if
        else
           W = 1.d0/sqrt(abs(aux))
        end if

        fluid_W(l,i) = W

        fluid_u(l,i) = fluid_W(l,i)*fluid_v(l,i)

!       Recalculate rho and e.

        fluid_rho(l,i) = fluid_cD(l,i)/fluid_W(l,i)

        aux = fluid_cD(l,i)*(fluid_W(l,i) - 1.d0) + (p1 + fluid_q(l,i))*(fluid_W(l,i)**2 - 1.d0)

        if (fluid_cE(l,i)>aux) then
           fluid_e(l,i) = (fluid_cE(l,i) - aux)/(fluid_cD(l,i)*W)
        else
           fluid_cE(l,i) = aux
           fluid_e(l,i)  = 0.d0
        end if

        fluid_e(l,i) = (fluid_cE(l,i) - aux)/(fluid_cD(l,i)*fluid_W(l,i))

!       Recalculate enthalpy and csi.

        fluid_h(l,i) = 1.d0 + fluid_e(l,i) + (fluid_p(l,i) + fluid_q(l,i))/fluid_rho(l,i)
        fluid_csi(l,i) = fluid_rho(l,i)*fluid_h(l,i)*fluid_W(l,i)


!    **********************
!    ***   CSI METHOD   ***
!    **********************

!    With this method we solve for csi=rho*h*W instead
!    of the pressure.  We choose a trial value of csi
!    (we use the value at the last time level), and
!    calculate the 4-velocity from:
!
!                                  4
!    u^r  =  ( S_r / csi ) / (A psi )
!
!
!    Next we can calculate the Lorentz factor from:
!
!                        1/2                 4      2  1/2
!    W = ( 1 + u_r u^r  )    =    ( 1 + A psi  (u^r)  )
!               
!
!    Notice that this guarantees that W is real.
!    We can now calculate the fluid speed as:
!
!    v^r  =  u^r / W
!
!    with this procedure we will always have v^2 < 1.
!
!    To proceed we caculate the energy density rho,
!    the enthalpy h and the pressure p as:
!
!    rho = D / W
!
!    h = csi / D
!
!    p = csi W - E - D
!
!    We can now use the fact that the enthalpy for
!    a perfect fluid can also be calculated as:
!
!    h' = ( 1 + gamma/(gamma-1) (rho/p) )
!
!    Wew then calculate the residue h-h', and iterate
!    until this difference reaches the desired tolerance.

     else if (fluid_primitive=="csi") then

!       Initialize counter and residue.

        j = 0
        res = 1.d0

!       Start iterations for Newton-Raphson's method.

        do while ((abs(res)>epsilon).and.(j<maxiter))

!          Increment counter.

           j = j + 1

!          Set initial guess for csi to old value.

           if (j==1) csi1 = fluid_csi(l,i)

!          Find 4-velocity from u^r = (S_r/csi) / (A psi^4)

           fluid_u(l,i) = fluid_cS(l,i)/csi1/(A(l,i)*exp(4.d0*phi(l,i)))

!          Find Lorentz factor from W = sqrt(1 + u^2).

           W2 = 1.d0 + A(l,i)*exp(4.d0*phi(l,i))*fluid_u(l,i)**2

           fluid_W(l,i) = sqrt(W2)

!          Find fluid speed.

           fluid_v(l,i) = fluid_u(l,i)/fluid_W(l,i)

!          Find rest-mass density: rho = D/W.

           fluid_rho(l,i) = fluid_cD(l,i)/fluid_W(l,i)

!          Find enthalpy.

           fluid_h(l,i) = csi1/fluid_cD(l,i)

!          Find pressure from p = csi*W - E - D.

           fluid_p(l,i) = csi1*fluid_W(l,i) - fluid_cE(l,i) - fluid_cD(l,i)

!          Calculate difference between trial value of the enthalpy and
!          its predicted value from equation of state:
!
!          h = ( 1 + gamma/(gamma-1) (rho/p) )

           f1 = (1.d0 + fluid_gamma/(fluid_gamma-1.d0)*fluid_p(l,i)/fluid_rho(l,i)) - fluid_h(l,i)

           if (j==1) then
              aux = (1.d0 + fluid_gamma/(fluid_gamma-1.d0)*fluid_p(l,i)/fluid_rho(l,i))
           else
              aux = csi1 - f1*(csi2-csi1)/(f2-f1)
           end if

           csi2 = csi1
           f2 = f1

           csi1 = aux

!          Calculate residual.

           res = abs(csi2 - csi1)

        end do

!       If we exceeded the maximum number of iterations
!       send message to screen.

        if (j==maxiter) then
           print *
           print *, 'Maximum iteration number reached in fluidprimitive.f90 at point: i = ',i,', r = ',r(l,i)
           print *, 'Residue = ',real(res)
           print *
           call die
        end if

!       Save last value of csi.

        fluid_csi(l,i) = csi1

!       Recalculate fluid speed and Lorentz number.

        fluid_u(l,i) = fluid_cS(l,i)/fluid_csi(l,i)/(A(l,i)*exp(4.d0*phi(l,i)))
        fluid_W(l,i) = sqrt(1.d0 + A(l,i)*exp(4.d0*phi(l,i))*fluid_u(l,i)**2)
        fluid_v(l,i) = fluid_u(l,i)/fluid_W(l,i)

!       Recalculate (rho,p,e,h).

        fluid_rho(l,i) = fluid_cD(l,i)/fluid_W(l,i)
        fluid_p(l,i)   = fluid_csi(l,i)*fluid_W(l,i) - fluid_cE(l,i) - fluid_cD(l,i)
        fluid_e(l,i)   = fluid_p(l,i)/(fluid_gamma - 1.d0)/fluid_rho(l,i)
        fluid_h(l,i)   = fluid_csi(l,i)/fluid_cD(l,i)

     end if


!    ****************************************
!    ***   END OF LOOP OVER GRID POINTS   ***
!    ****************************************

     10 continue

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine recover





