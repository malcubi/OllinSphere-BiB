
  subroutine fluidprimitive(l)

! *********************************************
! ***   RECOVER FLUID PRIMITIVE VARIABLES   ***
! *********************************************

! This routine recovers the fluid primitive variables (rho0,p,v)
! starting from the conserved quantities (D,E,S). The different
! quantities are related through:
!
! D   = rho0 / W
!
! S_i = rho0 h W^2 v_i
!
! E   = rho0 h W^2 - p - rho0 W
!
! with W the Lorentz factor:
!
! W = 1 / ( 1 - v_i v^i)
!
! and h the enthalpy:
!
! h = 1 + e + p /rho0
!
! where e is the specific internal energy which can be recovered
! from the definition of h as:
!
! e  =  h - 1 - p / rho0 
!
!    = [ E + rho0 W (1-W) + p (1-W^2) ] / rho W^2
!
!    = [ E + D (1-W) + p (1-W^2) ] / D W
!
! Since the system of equations above is coupled and non-linear,
! we do the inversion using Newton-Raphson's method.

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
! Avoid possible divisions by zero!

  do i=1-ghost,Nr
     if (fluid_rho(l,i)==0.d0) then
        fluid_h(l,i) = 0.d0
     else
        fluid_h(l,i) = 1.d0 + fluid_e(l,i) + fluid_p(l,i)/fluid_rho(l,i)
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

  logical :: atmos

  integer :: i,j,l
  integer :: maxiter = 100

  real(8) :: rhoatmos,Eatmos,patmos
  real(8) :: u1,u2,v1,v2,W1,W2
  real(8) :: rho1,rho2,e1,e2,h1,h2
  real(8) :: p1,p2,dp
  real(8) :: f1,f2,df
  real(8) :: csi1,csi2,dcsi
  real(8) :: aux

  real(8) :: epsilon = 1.d-12  ! Tolerance
  real(8) :: Wmax = 1.d5       ! Maximum allowed Lorentz factor


! ********************************
! ***   FIX ATMOSPHERE SCALE   ***
! ********************************

! Find atmosphere values.  We set both "rho0" and "e" to a small value,
! and obtain the pressure from the polytropic equation of state.
!
! Notice that we can't just set e=p=0, since in that case the speed of
! sound vanishes and the routine that calculates the fluid sources fails.

  rhoatmos = fluid_atmos/10.d0

  patmos = fluid_kappa*rhoatmos**fluid_gamma
  Eatmos = fluid_kappa*rhoatmos**(fluid_gamma-1.d0)/(fluid_gamma-1.d0)


! *********************************
! ***   LOOP OVER GRID POINTS   ***
! *********************************

  do i=1,Nr

     atmos = .false.


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

     if ((fluid_cD(l,i)<=fluid_atmos).or.(fluid_cE(l,i)<=rhoatmos*Eatmos)) then

        atmos = .true.

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

!    The idea of the "direct" method is to solve for the pressure p.
!
!    The way we proceed is as follows:
!
!    1) We take an initial value of p=p1 (the
!       value at the previous time level).
!
!    2) Using this value of p, and the new values of
!       the conserved quantities (D,E,S), we calculate
!       the velocity v, Lorentz factor W, density rho,
!       and specific internal energy e from:
!
!       v1 = S  / (E + D + p1) / ( A psi**4 )
!
!       W1 = 1 / sqrt( 1  -  A psi**4 v1**2 )
!
!       rho1 = D / W
!
!       e1  =  [ E + D (1 - W1) + p1 (1 - W1**2 ) ] / D W1

!       The factors of A psi**4 are there since the code
!       uses S with lower index, but u and v with upper index.
!
!    3) We calculate the difference between the pressure
!       p1 and the one we obtain from the equation of state:
!
!       f1 = p(EOS) - p1

!    4) We consider a increment in p: p2 = p1 + delta.
!       We then do the same calculation all over again in order
!       to obtain a new value of f=f2. This is because for
!       Newton-Raphson we need to calculate the derivative
!       of f, and we do it numerically:
!
!       f' = (f2-f1)/(p2-p1)
!
!    5) We update p1:
!
!       p1 = p1 - f1/h'  = p1 - f1 (p2-p1)/(f2-f1)
!
!    6) If abs(f1) is below the tolerance we end the iterations.

     if (fluid_primitive=="direct") then

!       Initial guess for pressure.

        p1 = fluid_p(l,i)

!       Start iterations for Newton-Raphson's method.

        do j=1,maxiter

!          Calculate f1 for p1.

           v1 = fluid_cS(l,i)/(fluid_cE(l,i) + fluid_cD(l,i) + p1)/(A(l,i)*exp(4.d0*phi(l,i)))

           aux = 1.d0 - A(l,i)*exp(4.d0*phi(l,i))*v1**2

           if (aux<=0.d0) then
              W1 = Wmax
           else
              W1 = 1.d0/sqrt(aux)
           end if

           rho1 = fluid_cD(l,i)/W1

           aux = fluid_cD(l,i)*(W1 - 1.d0) + p1*(W1**2 - 1.d0)
           e1 = (fluid_cE(l,i) - aux)/(fluid_cD(l,i)*W1)

           if (fluid_EOS=="ideal") then
              f1 = (fluid_gamma-1.d0)*rho1*e1 - p1
           else
              print *
              print *, 'Unknown equation of state'
              print *, 'Aborting! (subroutine fluidprimitive)'
              print *
              call die
           end if

!          Calculate f1 for p2.

           p2 = p1 + 1.d-8

           v2 = fluid_cS(l,i)/(fluid_cE(l,i) + fluid_cD(l,i) + p2)/(A(l,i)*exp(4.d0*phi(l,i)))

           aux = 1.d0 - A(l,i)*exp(4.d0*phi(l,i))*v2**2

           if (aux<=0.d0) then
              W2 = Wmax
           else
              W2 = 1.d0/sqrt(aux)
           end if

           rho2 = fluid_cD(l,i)/W2

           aux = fluid_cD(l,i)*(W2 - 1.d0) + p2*(W2**2 - 1.d0)
           e2 = (fluid_cE(l,i) - aux)/(fluid_cD(l,i)*W2)

           if (fluid_EOS=="ideal") then
              f2 = (fluid_gamma-1.d0)*rho2*e2 - p2
           else
              print *
              print *, 'Unknown equation of state'
              print *, 'Aborting! (subroutine fluidprimitive)'
              print *
              call die
           end if

!          Find dcsi and df.

           dp = p2 - p1
           df = f2 - f1

!          Update cs1.

           p1 = p1 - f1*dp/df

!          Did we reach the desired tolerance?

           if (abs(f1)<epsilon) exit

        end do

!       Save last value of p.

        fluid_p(l,i) = p1

!       Recalculate fluid speed and Lorentz number using new value of pressure.

        fluid_v(l,i) = fluid_cS(l,i)/(fluid_cE(l,i) + fluid_cD(l,i) + fluid_p(l,i)) &
                     /(A(l,i)*exp(4.d0*phi(l,i)))

        aux = 1.d0 - A(l,i)*exp(4.d0*phi(l,i))*fluid_v(l,i)**2

        if (aux<=0.d0) then
           atmos = .true.
        else
           fluid_W(l,i) = 1.d0/sqrt(aux)
        end if

        fluid_u(l,i) = fluid_W(l,i)*fluid_v(l,i)

!       Recalculate rho and e.

        fluid_rho(l,i) = fluid_cD(l,i)/fluid_W(l,i)

        aux = fluid_cD(l,i)*(fluid_W(l,i) - 1.d0) + p1*(fluid_W(l,i)**2 - 1.d0)
        fluid_e(l,i) = (fluid_cE(l,i) - aux)/(fluid_cD(l,i)*fluid_W(l,i))

!       Recalculate enthalpy and csi.

        fluid_h(l,i) = 1.d0 + fluid_e(l,i) + fluid_p(l,i)/fluid_rho(l,i)
        fluid_csi(l,i) = fluid_rho(l,i)*fluid_h(l,i)*fluid_W(l,i)


!    **********************
!    ***   CSI METHOD   ***
!    **********************

!    With this method we solve for:
!
!    csi = rho h W
!
!    The way we proceed is as follows:
!
!    1) We take an initial value of csi=csi1 (the
!       value at the previous time level).
!
!    2) Using this value of csi, and the new values of
!       the conserved quantities (D,E,S), we calculate
!       the 4-velocity u, Lorentz factor W, density rho,
!       pressure p, and enthalpy h from:
!
!       u1 = S / csi1 / A psi**4
!
!       W1 = sqrt(1 + A psi**4 u1**2)
!
!       rho1 = D/W1
!
!       p1 = csi1 W1 - E - D
!
!       h1 = csi1 / D
!
!       The factors of A psi**4 are there since the code
!       uses S with lower index, but u and v with upper index.
!
!    3) We calculate the difference between the enthalpy
!       h we just calculated and the one we obtain from
!       the equation of state:
!
!       f1 = h1(EOS) - h1(csi)
!
!    4) We consider a increment in csi:  csi2 = csi1 + delta.
!       We then do the same calculation all over again in order
!       to obtain a new value of f=f2. This is because for
!       Newton-Raphson we need to calculate the derivative
!       of f, and we do it numerically:
!
!       f' = (f2-f1)/(csi2-csi1)
!
!    5) We update csi1:
!
!       csi1 = csi1 - f1/h'  = csi1 - f1 (csi2-csi1)/(f2-f1)
!
!    6) If abs(f1) is below the tolerance we end the iterations.

     else if (fluid_primitive=="csi") then

!       Initial guess.

        csi1 = fluid_csi(l,i)

!       Start iterations for Newton-Raphson's method.

        do j=1,maxiter

!          Calculate f1 for csi1.

           u1 = fluid_cS(l,i)/csi1/(A(l,i)*exp(4.d0*phi(l,i)))

           W1 = sqrt(1.d0 + A(l,i)*exp(4.d0*phi(l,i))*u1**2)

           rho1 = fluid_cD(l,i)/W1
           p1   = csi1*W1 - fluid_cE(l,i) - fluid_cD(l,i)
           h1   = csi1/fluid_cD(l,i)

           f1 = (1.d0 + fluid_gamma/(fluid_gamma-1.d0)*p1/rho1) - h1

!          Calculate f2 for csi2.

           csi2 = csi1 + 1.d-8

           u2 = fluid_cS(l,i)/csi1/(A(l,i)*exp(4.d0*phi(l,i)))

           W2 = sqrt(1.d0 + A(l,i)*exp(4.d0*phi(l,i))*u2**2)

           rho2 = fluid_cD(l,i)/W2
           h2   = csi2/fluid_cD(l,i)
           p2   = csi2*W2 - fluid_cE(l,i) - fluid_cD(l,i)

           f2 = (1.d0 + fluid_gamma/(fluid_gamma-1.d0)*p2/rho2) - h2

!          Find dcsi and df.

           dcsi = csi2 - csi1
           df = f2 - f1

!          Update cs1.

           csi1 = csi1 - f1*dcsi/df

!          Did we reach the desired tolerance?

           if (abs(f1)<epsilon) exit

        end do

!       Save last value of csi.

        fluid_csi(l,i) = csi1

!       Recalculate fluid 4-velocity, Lorentz factor and fluid speed.

        fluid_u(l,i) = fluid_cS(l,i)/fluid_csi(l,i)/(A(l,i)*exp(4.d0*phi(l,i)))
        fluid_W(l,i) = sqrt(1.d0 + A(l,i)*exp(4.d0*phi(l,i))*fluid_u(l,i)**2)
        fluid_v(l,i) = fluid_u(l,i)/fluid_W(l,i)

!       Recalculate (rho,p,h).

        fluid_rho(l,i) = fluid_cD(l,i)/fluid_W(l,i)
        fluid_p(l,i)   = fluid_csi(l,i)*fluid_W(l,i) - fluid_cE(l,i) - fluid_cD(l,i)
        fluid_h(l,i)   = fluid_csi(l,i)/fluid_cD(l,i)

!       Calculate the specific internal energy e.

        aux = fluid_cD(l,i)*(fluid_W(l,i) - 1.d0) + p1*(fluid_W(l,i)**2 - 1.d0)
        fluid_e(l,i) = (fluid_cE(l,i) - aux)/(fluid_cD(l,i)*fluid_W(l,i))

     end if


!    ************************
!    ***   DID WE FAIL?   ***
!    ************************

!    Check if we are below atmosphere levels.

     if ((fluid_rho(l,i)<rhoatmos).or.(fluid_p(l,i)<patmos).or.(fluid_e(l,i)<eatmos)) atmos=.true.

!    If we exceeded the maximum number of iterations or failed
!    in other ways we set everything to the atmosphere.

     if ((j==maxiter).or.atmos) then

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





