!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/fluidprimitive.f90,v 1.11 2023/03/24 18:00:08 malcubi Exp $

  subroutine fluidprimitive(l)

! *********************************************
! ***   RECOVER FLUID PRIMITIVE VARIABLES   ***
! *********************************************

! This routine recovers the fluid primitive variables (rho0,p,v)
! starting from the conserved quantities (D,E,S).
!
! Since the system of equations we need to invert is coupled and
! non-linear, we do the inversion using Newton-Raphson's method.
!
! The idea is to choose a trial value of the pressure (the value
! at the old time level).  From this we calculate first the fluid
! velocity from:
!
!  r                                4 phi
! v   =  S / (E + D + p + q) / ( A e     )
!         r
!
! where the geometric factor is there in order to raise the index
! (v has index up, while S has index down), and where "q" is the
! contribution to the pressure coming from the artificial viscosity.
!
! We then calculate the Lorentz factor W from:
!
!                          4 phi  2
! W  =  1 / sqrt( 1  -  A e      v  )
!
!
! Notice that we must always have:  A exp(4 phi) v^2 < 1.
!
! The rest-mass energy density is then obtained as:
!
! rho  =  D / W
!
! and the specific internal energy is:
!
!                                       2
! e  =  ( E + D (1 - W) + (p + q) (1 - W ) )
!
!
! The above expression can be found from the definitions of
! E and the enthalpy h.
!
! Finally, we also calculate the enthalpy as:
!
! h = 1 + e + (p + q)/rho
!
! and the speed of sound from the equation of state.
!
! Having found these variables, we find the difference between
! the trial value of the pressure and the one obtained from
! the equation of state, and iterate until this difference
! becomes zero.

! Include modules.

  use param
  use arrays
  use derivatives

! Extra variables.

  implicit none

  integer i,j,l,maxiter

  real(8) p1,p2,f1,f2
  real(8) rhoatmos,Eatmos,patmos
  real(8) W,res,aux
  real(8) :: epsilon = 1.d-10


! ***************************************
! ***   RECOVER PRIMITIVE VARIABLES   ***
! ***************************************

  maxiter = 500

! Find atmosphere values (from polytropic equation of state).

  rhoatmos = fluid_atmos

  patmos = fluid_kappa*rhoatmos**fluid_gamma
  Eatmos = fluid_kappa/(fluid_gamma-1.0)*rhoatmos**(fluid_gamma-1.0)

! Loop over grid points.

  do i=1,Nr

!    Atmosphere: If the density is too low we run the risk of dividing
!    by a very small quantity to recover the fluid speed, which will
!    introduce very large round-off errors.  In order to avoid this
!    we set the conserved mass density D to a small value and the
!    momentum density S to 0. For the conserved energy density we
!    use the polytropic equation of state, which implies:
!
!    E  =  rho0 e  =  kappa rho0**gamma / (gamma - 1)
!
!    We then set all other variables to their values consistent with
!    this and we jump out of here.

     if ((fluid_cD(l,i)<=1.5d0*rhoatmos).or.(fluid_p(l,i)<=1.5d0*patmos)) then

!       Primitive variables.

        fluid_rho(l,i) = rhoatmos

        fluid_p(l,i) = patmos
        fluid_e(l,i) = Eatmos

!       Fluid speed and Lorentz factor.

        fluid_v(l,i)  = 0.d0
        fluid_W(l,i)  = 1.d0

!       Enthalpy.

        fluid_h(l,i) = 1.d0 + fluid_e(l,i) + fluid_p(l,i)/fluid_rho(l,i)

!       Speed of sound (assuming a polytropic relation):  vs^2 = gamma p / rho0 h

        fluid_vs(l,i) = sqrt(abs(fluid_gamma*fluid_p(l,i)/(fluid_rho(l,i)*fluid_h(l,i))))

!       Conserved quantities (D,E,S).

        fluid_cD(l,i) = fluid_rho(l,i)
        fluid_cE(l,i) = fluid_rho(l,i)*fluid_e(l,i)
        fluid_cS(l,i) = 0.d0

        goto 10

     end if

!    Initialize counter and residue.

     j = 0
     res = 1.d0

!    Set initial guess for pressure to old value.

     p1 = fluid_p(l,i)

!    Start iterations for Newton-Raphson's method.

     do while ((abs(res)>epsilon).and.(j<maxiter))

!       Increment counter.

        j = j + 1

!       Find value of v using old value of pressure.

        if (fluid_cS(l,i)==0.d0) then
           fluid_v(l,i) = 0.d0
        else
           fluid_v(l,i) = fluid_cS(l,i)/(fluid_cE(l,i) + fluid_cD(l,i) + p1 + fluid_q(l,i)) &
                        /(A(l,i)*exp(4.d0*phi(l,i)))
        end if

!       Find Lorentz factor W, and check that v^2 is less than 1.

        aux = 1.d0 - A(l,i)*exp(4.d0*phi(l,i))*fluid_v(l,i)**2

        if (aux<=0.d0) then
           W = 1.d10
        !   print *, "Fluid speed is larger than 1"
        else
           W = 1.d0/sqrt(abs(aux))
        end if

        W = 1.d0/sqrt(abs(aux))
        fluid_W(l,i) = W

!       Find specific internal energy (it should be positive).

        aux = fluid_cD(l,i)*(W - 1.d0) + (p1 + fluid_q(l,i))*(W**2 - 1.d0)

        if (fluid_cE(l,i)>aux) then
           fluid_e(l,i) = (fluid_cE(l,i) - aux)/(fluid_cD(l,i)*W)
        else
           fluid_cE(l,i) = aux
           fluid_e(l,i)  = 0.d0
        end if

!       Find rest-mass density.

        fluid_rho(l,i) = fluid_cD(l,i)/W

!       Find enthalpy.

        fluid_h(l,i) = 1.d0 + fluid_e(l,i) + (p1 + fluid_q(l,i))/fluid_rho(l,i)

!       Find speed of sound (at the moment only for ideal gas equation of state):
!
!       vs^2  = gamma (gamma-1) e / ( 1 + gamma e)  =  p gamma (gamma - 1) / (p gamma + rho0 (gamma - 1))

        if (fluid_EOS=="ideal") then

           fluid_vs(l,i) = sqrt(abs(p1*fluid_gamma*(fluid_gamma-1.d0) &
                         /(p1*fluid_gamma + fluid_rho(l,i)*(fluid_gamma - 1.d0))))

!       If we don't have an equation of state we just use the polytropic relation
!       in order to avoid problems in other routines if vs is not defined.  But
!       really this is not a good idea.

        else

           fluid_vs(l,i) = sqrt(abs(fluid_gamma*fluid_p(l,i)/(fluid_rho(l,i)*fluid_h(l,i))))

        end if

!       Update artificial viscosity:
!
!       1) If the divergence of S is positive then set
!       fluid_q to zero.
!
!       2) If the divergence is negative compute the artificial
!       viscosity contribution. The artificial viscosity has
!       two contributions and has the following general form:
!
!       q  =  dr abs(div.S) ( q1 vs + q2 dr abs(div.S) / rho )
!
!       with q1 and q2 positive constants, and vs the speed
!       of sound. Notice that the term with q1 introduces a first
!       order error, so one must keep q1 small.  The term with q2
!       introduces a second order error. For rho we use the ADM
!       energy density: rho = E + D.

        if ((fluid_q1/=0.d0).and.(fluid_q2/=0.d0)) then
           aux = dr(l)*D1_fluid_cS(l,i)/(A(l,i)*exp(4.d0*phi(l,i)))
           if (aux>=0.d0) then
              fluid_q(l,i) = 0.d0
           else
              fluid_q(l,i) = abs(aux)*(fluid_q1*fluid_vs(l,i) &
                           + fluid_q2*abs(aux)/(fluid_cE(l,i) + fluid_cD(l,i)))
           end if
        else
           fluid_q(l,i) = 0.d0
        end if

!       Calculate difference between trial value of pressure
!       and predicted value from the equation of state.
!
!       If there is NO equation of state (for example if you forced
!       some specific form of rho for the initial data and solved
!       the TOV equations), then we just set the difference to 0.

        if (fluid_EOS/="none") then

           if (fluid_EOS=="ideal") then
              f1 = (fluid_gamma - 1.d0)*fluid_rho(l,i)*fluid_e(l,i) - p1
           else
              print *
              print *, 'Unknown equation of state'
              print *, 'Aborting! (subroutine fluidprimitive)'
              print *
              call die
           end if

        else

           f1 = 0.d0

        end if
         
!       Update the value of the pressure for the next iteration.

        if (fluid_EOS/="none") then

           if (j==1) then
              if (fluid_EOS=="ideal") then
                 aux = (fluid_gamma - 1.d0)*fluid_rho(l,i)*fluid_e(l,i)
              else
                 aux = fluid_p(l,i)
              end if
           else
              aux = p1 - f1*(p2-p1)/(f2-f1)
           end if

           p2 = p1
           f2 = f1

           p1 = aux
           if (p1<=patmos) p1=patmos

!          Calculate residual.

           res = abs(p2 - p1)

!       No equation of state.

        else

           res = 0.d0

        end if

     end do

!    Save last value of p1 as the new pressure.

     fluid_p(l,i) = p1

!    If we exceeded the maximum number of iterations send message to screen.

     if (j==maxiter) then
        print *
        print *, 'Maximum iteration number reached in fluidprimitive.f90 at point: i = ',i,', r = ',r(l,i)
        print *, 'Residue = ',real(res)
        print *
        call die
     end if

!    Simple continue for the goto above (yes, a goto, so what?).

10   continue

!    Calculate Mach number.

     fluid_Mach(l,i) = fluid_v(l,i)/fluid_vs(l,i)

!    Characteristic speeds.

     aux = fluid_vs(l,i)/sqrt(A(l,i))*abs(1.d0 - A(l,i)*fluid_v(l,i)**2)

     fluid_vcp(l,i) = alpha(l,i)/(1.d0 - A(l,i)*fluid_v(l,i)**2*fluid_vs(l,i)**2) &
                    *(fluid_v(l,i)*(1.d0 - fluid_vs(l,i)**2) + aux)
     fluid_vcm(l,i) = alpha(l,i)/(1.d0 - A(l,i)*fluid_v(l,i)**2*fluid_vs(l,i)**2) &
                    *(fluid_v(l,i)*(1.d0 - fluid_vs(l,i)**2) - aux)

  end do

! Ghost points.

  do i=1,ghost
     fluid_rho(l,1-i) = + fluid_rho(l,i)
     fluid_p(l,1-i)   = + fluid_p(l,i)
     fluid_e(l,1-i)   = + fluid_e(l,i)
     fluid_h(l,1-i)   = + fluid_h(l,i)
     fluid_v(l,1-i)   = - fluid_v(l,i)
     fluid_W(l,1-i)   = + fluid_W(l,i)
     fluid_q(l,1-i)   = + fluid_q(l,i)
     fluid_vs(l,1-i)  = + fluid_vs(l,i)
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine fluidprimitive

