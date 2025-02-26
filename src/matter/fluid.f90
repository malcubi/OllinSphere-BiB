!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/fluid.f90,v 1.16 2025/02/26 20:11:26 malcubi Exp $

  subroutine sources_fluid(l)

! *****************************
! ***   SOURCES FOR FLUID   ***
! *****************************

! This routine calculates the sources for the relativistic Euler
! equations for a perfect fluid, written in conservative form.
!
! These equations take the general form:
!
!                                                                         1/2
! d D  =  beta d D  -  d [ alpha v D ] +  alpha trK D  -  alpha v D d ln g
!  t            r       r                                            r
!
!
! d E  =  beta d E  -  d [ alpha v (E + p) ]  +  alpha trK (E + p)
!  t            r       r  
!                              2    4 phi
!      +  (E + p + D) [ alpha v  A e     (KTA + trK/3) - v d alpha ]
!                                                           r
!                               1/2
!      -  alpha v (E + p) d ln g
!                          r
!
! d S  =  beta d S  +  S d beta  -  d [ alpha ( v S + p ) ]  +  alpha trK S
!  t            r         r          r
!
!      -  alpha v S [ d B / B  +  4 d phi  +  2/r ]  -  (E + D) d alpha
!                      r             r                           r
!
!
! where g is the determinant of the spatial metric and:
!
!       1/2
! d ln g   =  d g / 2g  =  d A / 2A + d B / B + 6 d phi + 2/r
!  r           r            r          r           r
!
! Notice that when we use artificial viscosity we need
! to add an extra contribution to the pressure:
!
! p  ->  p + q
!
! We use a flux conservative formulation, so that the
! sources are calculated as:
!
! source = ( F      -  F      ) / dr
!             i+1/2     i-1/2
!
! The routine can use several different differencing
! methods for the flux terms.
!
! Notice that for notation we in fact use:
!
! F      =  F(i-1)
!  i-1/2
!
! F      =  F(i)
!  i+1/2
!
! In the routine we first calculate the flux terms, and
! at the end we add all other source tems (including
! terms with Christoffels and shift terms).

! Include modules.

  use param
  use arrays
  use derivatives
  use derivadvect

! Extra variables.

  implicit none

  integer i,l

  real(8) one,half,third
  real(8) slope1,slope2,slopelim
  real(8) vpp,vmm
  real(8) idr,aux

  real(8) flux_D(1-ghost:Nr),flux_DL(1-ghost:Nr),flux_DR(1-ghost:Nr)
  real(8) flux_E(1-ghost:Nr),flux_EL(1-ghost:Nr),flux_ER(1-ghost:Nr)
  real(8) flux_S(1-ghost:Nr),flux_SL(1-ghost:Nr),flux_SR(1-ghost:Nr)

  real(8) D_L(1-ghost:Nr),D_R(1-ghost:Nr)
  real(8) E_L(1-ghost:Nr),E_R(1-ghost:Nr)
  real(8) S_L(1-ghost:Nr),S_R(1-ghost:Nr)

  real(8) vp_L(1-ghost:Nr),vp_R(1-ghost:Nr)
  real(8) vm_L(1-ghost:Nr),vm_R(1-ghost:Nr)

  character(len=20) limiter


! *******************
! ***   NUMBERS   ***
! *******************

  one   = 1.d0
  half  = 0.5d0
  third = 1.d0/3.d0

  idr = 1.d0/dr(l)


! ****************************
! ***   CALCULATE FLUXES   ***
! ****************************

! The flux terms at the cell interfaces (between grid points)
! can be calculated using several different methods.  At the
! moment all these methods are either first or second order.
!
! 1) Center:   Standard centered average (second order).
!
! 2) Upwind:   Copy from the adjacent grid point from the upwind
!              side (first order).
!
! 3) Limiter:  The advection term in the flux is either interpolated
!              (averaged) or extrapolated from the upwind side using
!              a minmod limiter (second order).
!
! 4) LLF:      Local Lax-Friedrichs.  First order scheme that should
!              work reasonably well.
!
! 5) HLLE:     Approximate Riemann solver of Harten, Lax, van Leer
!              and Einfeldt.

! Initialize fluxes to zero.

  flux_D = 0.d0
  flux_E = 0.d0
  flux_S = 0.d0


! **********************
! ***   LLF METHOD   ***
! **********************

! Local Lax-Friedrichs method.  This is only first order.

  if (fluid_method=="llf") then

!    NOT YET IMPLEMENTED


! **************************************
! ***   HLLE METHOD (SECOND ORDER)   ***
! **************************************

! Here I use an HLLE flux calculated from the reconstructed
! primitive variables at the grid interfaces. For this we
! first need to do a linear reconstruction of the fluid
! variables at cell interfaces using both left and right
! sided extrapolations.
!
! This method is more complicated, and last I checked the
! method "limiter" above seems to be more accurate. I need
! to check if I don't have anything wrong here.

  else if (fluid_method=="hlle") then

     limiter = fluid_limiter

!    Reconstruct conserved quantities at cell interfaces.

     call reconstruct(fluid_cD(l,:),D_L(:),D_R(:),limiter,+1)
     call reconstruct(fluid_cE(l,:),E_L(:),E_R(:),limiter,+1)
     call reconstruct(fluid_cS(l,:),S_L(:),S_R(:),limiter,-1)

!    Reconstruct fluxes at cell interfaces.

     flux_D = alpha(l,:)*fluid_v(l,:)*fluid_cD(l,:)
     call reconstruct(flux_D,flux_DL,flux_DR,limiter,-1)

     flux_E = alpha(l,:)*fluid_v(l,:)*(fluid_cE(l,:) + fluid_p(l,:) + fluid_q(l,:))
     call reconstruct(flux_E,flux_EL,flux_ER,limiter,-1)

     flux_S = alpha(l,:)*(fluid_v(l,:)*fluid_cS(l,:) + fluid_p(l,:))
     call reconstruct(flux_S,flux_SL,flux_SR,limiter,+1)

!    Reconstruct characteristic speeds at cell boundaries.

     if (fluid_usesoundspeed) then
        call reconstruct(fluid_vcp(l,:),vp_L(:),vp_R(:),limiter,+1)
        call reconstruct(fluid_vcm(l,:),vm_L(:),vm_R(:),limiter,+1)
     end if

!    Find fluxes.

     do i=0,Nr-1

!       Find maximum and minimum speeds using the
!       characteristic speeds coming from the speed
!       of sound.

        if (fluid_usesoundspeed) then

           vpp = max(0.d0,vp_L(i),vp_R(i))
           vmm = min(0.d0,vm_L(i),vm_R(i))

           if (vpp-vmm==0.d0) then
              print *
              print *, 'Problem in HLLE: Maximum and minimum speeds (vmm,vpp) are equal: ',vmm,vpp
              print *, 'at point: ', i
              print *, 'Aborting! (subroutine fluid.f90)'
              print *
              call die
           end if

!       If we don't know the speed of sound, we can use
!       the local speed of light instead.  This works just
!       fine, though it is less accurate (more dissipative).

        else

           aux = 0.5d0*(alpha(l,i  )/sqrt(abs(A(l,i  )))/psi(l,i  )**2 &
                      + alpha(l,i+1)/sqrt(abs(A(l,i+1)))/psi(l,i+1)**2)

           if (shift=="none") then
              vpp = max(0.d0,+aux)
              vmm = min(0.d0,-aux)
           else
              vpp = max(0.d0,-0.5d0*(beta(l,i)+beta(l,i+1))+aux)
              vmm = min(0.d0,-0.5d0*(beta(l,i)+beta(l,i+1))-aux)
           end if

        end if

!       Calculate HLLE fluxes.

        flux_D(i) = (vpp*flux_DL(i) - vmm*flux_DR(i) + vpp*vmm*(D_R(i) - D_L(i)))/(vpp - vmm)
        flux_E(i) = (vpp*flux_EL(i) - vmm*flux_ER(i) + vpp*vmm*(E_R(i) - E_L(i)))/(vpp - vmm)
        flux_S(i) = (vpp*flux_SL(i) - vmm*flux_SR(i) + vpp*vmm*(S_R(i) - S_L(i)))/(vpp - vmm)

     end do


! **************************
! ***   UNKNOWN METHOD   ***
! **************************

! This should never happen, but just in case.

  else

     print *, 'Unknown fluid_method ...'
     print *, 'Aborting! (subroutine fluid.f90)'
     print *
     call die

  end if


! *******************************************
! ***   SOURCE CONTRIBUTION FROM FLUXES   ***
! *******************************************

  do i=1,Nr-1
     sfluid_cD(l,i) = - idr*(flux_D(i) - flux_D(i-1))
     sfluid_cE(l,i) = - idr*(flux_E(i) - flux_E(i-1))
     sfluid_cS(l,i) = - idr*(flux_S(i) - flux_S(i-1))
  end do


! **************************************
! ***   ADD ALL OTHER SOURCE TERMS   ***
! **************************************

! Here we add all other source terms, including terms related to
! the extrinsic curvature, the Christoffel symbols and the shift.

! Sources for D.  The terms we are still missing are:
!
!                                     1/2
!  +  alpha trK D  -  alpha v D d ln g
!                                r
!
! with:
!
!       1/2
! d ln g   =  d g / 2g  =  d A / 2A + d B / B + 6 d phi + 2/r
!  r           r            r          r           r

  sfluid_cD(l,:) = sfluid_cD(l,:) + alpha(l,:)*trK(l,:)*fluid_cD(l,:) &
                 - alpha(l,:)*fluid_v(l,:)*fluid_cD(l,:) &
                 *(half*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 6.d0*D1_phi(l,:) + 2.d0/r(l,:))

! Shift terms for D (scalar):
!
!   +  beta d D
!            r

  if (shift/="none") then
     sfluid_cD(l,:) = sfluid_cD(l,:) + beta(l,:)*DA_fluid_cD(l,:)
  end if

! Sources for E.  The terms we are still missing are:
!
!   +  alpha trK (E + p)
!
!                                4 phi
!   +  (E + p + D) [ alpha v  A e      (KTA + trK/3) - v d alpha ]
!                                                         r
!                            1/2
!   -  alpha v (E + p) d ln g
!                       r
!
! with:
!
!       1/2
! d ln g   =  d g / 2g  =  d A / 2A + d B / B + 6 d phi + 2/r
!  r           r            r          r           r

  sfluid_cE(l,:) = sfluid_cE(l,:) + alpha(l,:)*trK(l,:)*(fluid_cE(l,:) + fluid_p(l,:) + fluid_q(l,:)) &
                 + (fluid_cE(l,:) + fluid_cD(l,:) + fluid_p(l,:) + fluid_q(l,:)) &
                 *(alpha(l,:)*fluid_v(l,:)**2*A(l,:)*exp(4.d0*phi(l,:)) &
                 *(KTA(l,:) + third*trK(l,:)) - fluid_v(l,:)*D1_alpha(l,:)) &
                 - alpha(l,:)*fluid_v(l,:)*(fluid_cE(l,:) + fluid_p(l,:) + fluid_q(l,:)) &
                 *(half*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 6.d0*D1_phi(l,:) + 2.d0/r(l,:))

! Shift terms for E (scalar):
!
!   +  beta d E
!            r

  if (shift/="none") then
     sfluid_cE(l,:) = sfluid_cE(l,:) + beta(l,:)*DA_fluid_cE(l,:)
  end if

! Sources for S.  The terms we are still missing are:
!
!   +  alpha trK S
!
!
!   -  alpha v S [ d B / B  +  4 d phi  +  2/r ]
!                   r             r
!
!   -  (E + D + p) d alpha  -  alpha d p
!                   r                 r

  sfluid_cS(l,:) = sfluid_cS(l,:) + alpha(l,:)*trK(l,:)*fluid_cS(l,:) &
                 - alpha(l,:)*fluid_v(l,:)*fluid_cS(l,:) &
                 *(D1_B(l,:)/B(l,:) + 4.d0*D1_phi(l,:) + 2.d0/r(l,:)) &
                 - (fluid_cE(l,:) + fluid_cD(l,:))*D1_alpha(l,:)

! Shift terms for S (1-form):
!
!   +  beta d S  +  S d beta
!            r         r

  if (shift/="none") then
     sfluid_cS(l,:) = sfluid_cS(l,:) + beta(l,:)*DA_fluid_cS(l,:) + fluid_cS(l,:)*D1_beta(l,:)
  end if


! *************************************
! ***   SOURCES AT OUTER BOUNDARY   ***
! *************************************
 
! For the moment I just copy the sources from one point in.

  sfluid_cD(l,Nr) = sfluid_cD(l,Nr-1)
  sfluid_cE(l,Nr) = sfluid_cE(l,Nr-1)
  sfluid_cS(l,Nr) = sfluid_cS(l,Nr-1)


! ************************
! ***   GHOST POINTS   ***
! ************************

! Ghost points using symmetries.

  do i=1,ghost
     sfluid_cD(l,1-i) = + sfluid_cD(l,i)
     sfluid_cE(l,1-i) = + sfluid_cE(l,i)
     sfluid_cS(l,1-i) = - sfluid_cS(l,i)
  end do


! ***********************
! ***   DISSIPATION   ***
! ***********************

! This is here just for testing simple methods. 
! In practice it seems that adding extra
! dissipation to the fluid methods makes things
! go unstable.

  if (fluiddiss/=0.d0) then

     dissipvar => fluid_cD
     sourcevar => sfluid_cD
     call dissipation(l,+1,fluiddiss)

     dissipvar => fluid_cE
     sourcevar => sfluid_cE
     call dissipation(l,+1,fluiddiss)

     dissipvar => fluid_cS
     sourcevar => sfluid_cS
     call dissipation(l,-1,fluiddiss)

  end if


! ***********************************
! ***   COSMOLOGICAL BACKGROUND   ***
! ***********************************

! Sources for cosmological background when needed.


! ***************
! ***   END   ***
! ***************

  end subroutine sources_fluid








  subroutine reconstruct(var,varl,varr,limiter,sym)

! ****************************************************
! ***   LINEAR RECONSTRUCTION AT CELL INTERFACES   ***
! ****************************************************

! This routine finds a linear reconstruction of the array
! "var" at cell interfaces. The reconstruction is done
! with either a centered interpolation or a one sided
! extrapolation using a slope limiter.
!
! We do both left and right sided reconstructions, and store
! them in the arrays "varl" and "varr".
!
! The slow limiters work by reconstructing a function at the
! half interval between grid points using the slopes to the
! left, center and right:
!
!     slope1 = var(i  ) - var(i-1)
!     slope2 = var(i+1) - var(i  )
!     slope3 = var(i+2) - var(i+1)
!
! For each case we need to do a reconstruction from the left
! and from the right, as this is needed for the Riemann solver.
! For the left-sided reconstruction we take:
!
!      varl(i)  = var(i) + 0.5d0*slopelim
!
! and for the right-sided reconstruction:
!
!     varr(i)  = var(i+1) - 0.5d0*slopelim
!
! where "slopelim" is an approximation to the slope obtained
! be using the ratio of the center and left slopes (for varl),
! and center and right slopes (for varr).
!
! For the left-sided reconstruction we take:
!
!     slopelim = phi*slope2
!
! and for the right-sided reconstruction:
!
!     slopelim = phi*slope2
!
! where "phi" is a function of the ratios (r=slope1/slope2)
! and (r=slope2/slope3) respectively.
!
! The function "phi" is called the "limiter", and can be chosen
! in a variety of different ways, but in order to give us a TVD
! method it must satisfy the following conditions:
!
!     r <=  phi(r) <= 2r      for:   0 <= r <= 1/2
!     r >=  phi(r) <= 1       for:   1/2 <= r <= 1
!     phi(1) = 1
!     1 <=  phi(r) <= r       for:   1 <= r <= 2
!     1 <=  phi(r) <= 2       for    r > 2
 
! Include modules.

  use param

! Extra variables.

  implicit none

  integer i,sym

  real(8) slope1,slope2,slope3,slopelim
  real(8) ratio,phi,philimiter
  real(8) var(1-ghost:Nr),varl(1-ghost:Nr),varr(1-ghost:Nr)

  character(len=*) limiter


! ***********************************************
! ***   FIND LEFT AND RIGHT RECONSTRUCTIONS   ***
! ***********************************************

! Sanity check.

  if (abs(sym)/=1) then
     print *
     print *, 'The symmetry of a variable must be +-1, aborting ...'
     print *
     call die
  end if


! ***************************
! ***   INTERIOR POINTS   ***
! ***************************

  do i=1,Nr-2

!    Find slopes.

     slope1 = var(i  ) - var(i-1)
     slope2 = var(i+1) - var(i  )
     slope3 = var(i+2) - var(i+1)

!    Left side reconstruction.

     if (slope1*slope2>0.d0) then
        ratio = slope1/slope2
        phi = philimiter(limiter,ratio)
     else
        phi = 0.d0
     end if

     slopelim = phi*slope2
     varl(i)  = var(i) + 0.5d0*slopelim

!    Right side reconstruction.

     if (slope2*slope3>0.d0) then
        ratio = slope2/slope3
        phi = philimiter(limiter,ratio)
     else
        phi = 0.d0
     end if

     slopelim = phi*slope3
     varr(i)  = var(i+1) - 0.5d0*slopelim

  end do


! ******************
! ***   ORIGIN   ***
! ******************

! Find slopes.  Notice that for "slope1", which
! corresponds to var(0) - var(-1), we use the symmetry
! of the corresponding variable, so that:
!
! var(0) - var(-1) = +-(var(1) - var(2))
!
! with the + sign for even variables and the
! - for odd variables.

  slope1 = sym*(var(1) - var(2))
  slope2 = var(1) - var(0)
  slope3 = var(2) - var(1)

! Left side reconstruction.

  if (slope1*slope2>0.d0) then
     ratio = slope1/slope2
     phi = philimiter(limiter,ratio)
  else
     phi = 0.d0
  end if

  slopelim = phi*slope2
  varl(0) = var(0) + 0.5d0*slopelim

! Right side reconstruction.

  if (slope2*slope3>0.d0) then
     ratio = slope2/slope3
     phi = philimiter(limiter,ratio)
  else
     phi = 0.d0
  end if

  slopelim = phi*slope3
  varr(0) = var(1) - 0.5d0*slopelim


! **************************
! ***   OUTER BOUNDARY   ***
! **************************

! Find slopes.  In this case there is no slope3.

  slope1 = var(Nr-1) - var(Nr-2)
  slope2 = var(Nr  ) - var(Nr-1)

! Left side reconstruction.

  if (slope1*slope2>0.d0) then
     ratio = slope1/slope2
     phi = philimiter(limiter,ratio)
  else
     phi = 0.d0
  end if

  slopelim = phi*slope2
  varl(Nr-1) = var(Nr-1) + 0.5d0*slopelim

! Right side reconstruction.

  varr(Nr-1) = var(Nr  ) - 0.5d0*slope2


! ***************
! ***   END   ***
! ***************

  end subroutine reconstruct









! ******************************
! ***   DIFFERENT LIMITERS   ***
! ******************************

 function philimiter(limiter,ratio)

 real(8) philimiter
 real(8) ratio
 real(8) aux1,aux2,beta

 character(len=20) limiter

! 1) minmod limiter.

  if (limiter=="minmod") then

     philimiter = min(1.d0,ratio)

! 2) vanleer limiter.

  else if (limiter=="vanleer") then

     philimiter = 2.d0*ratio/(1.d0+ratio)

! 3) superbee limiter.

  else if (limiter=="superbee") then

     aux1 = min(1.d0,2.d0*ratio)
     aux2 = min(2.d0,ratio)
     philimiter = max(aux1,aux2)

! 4) monotonized central (MC) limiter.

  else if (limiter=="mc") then

     philimiter = min(2.d0*ratio,0.5d0*(1.d0+ratio),2.d0)

! 5) koren limiter.

  else if (limiter=="koren") then

     philimiter = min(2.d0*ratio,(1.d0+2.d0*ratio)/3.d0,2.d0)

! 6) ospre limiter.

  else if (limiter=="ospre") then

     philimiter = 1.5d0*(ratio + ratio**2)/(1.d0 + ratio + ratio**2)

! 7) sweby limiter (beta=1.5).  Notice that superbee is the same
!    as sweby with beta=2.

  else if (limiter=="sweby") then

     beta = 1.5d0

     aux1 = min(1.d0,beta*ratio)
     aux2 = min (beta,ratio)
     philimiter = max(aux1,aux2)
 
  end if


! ***************
! ***   END   ***
! ***************

 end function philimiter
