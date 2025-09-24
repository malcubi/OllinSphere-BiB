!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/fluid.f90,v 1.20 2025/09/24 23:57:40 malcubi Exp $

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

  use procinfo
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

  real(8) flux_D(1-ghost:Nrmax),flux_DL(1-ghost:Nrmax),flux_DR(1-ghost:Nrmax)
  real(8) flux_E(1-ghost:Nrmax),flux_EL(1-ghost:Nrmax),flux_ER(1-ghost:Nrmax)
  real(8) flux_S(1-ghost:Nrmax),flux_SL(1-ghost:Nrmax),flux_SR(1-ghost:Nrmax)

  real(8) D_L(1-ghost:Nrmax),D_R(1-ghost:Nrmax)
  real(8) E_L(1-ghost:Nrmax),E_R(1-ghost:Nrmax)
  real(8) S_L(1-ghost:Nrmax),S_R(1-ghost:Nrmax)

  real(8) vp_L(1-ghost:Nrmax),vp_R(1-ghost:Nrmax)
  real(8) vm_L(1-ghost:Nrmax),vm_R(1-ghost:Nrmax)

  character(len=20) limiter


! *******************
! ***   NUMBERS   ***
! *******************

  one   = 1.d0
  half  = 0.5d0
  third = 1.d0/3.d0

  idr = 1.d0/dr(l)


! ********************************
! ***   NO EQUATION OF STATE   ***
! ********************************

! If we don't have an equation of state we force
! the method to use the speed of light instead
! of the speed of sound.

  if (fluid_EOS=="none") then
     fluid_usesoundspeed = .false.
  end if


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


! ****************************
! ***   CALCULATE FLUXES   ***
! ****************************

! Calculate fluxes at grid points.

  flux_D(:) = alpha(l,:)*fluid_v(l,:)*fluid_cD(l,:) &
            *(sqrt(A(l,:))*B(l,:)*exp(6.d0*phi(l,:)))

  flux_E(:) = alpha(l,:)*fluid_v(l,:)*(fluid_cE(l,:) + fluid_p(l,:) + fluid_q(l,:)) &
            *(sqrt(A(l,:))*B(l,:)*exp(6.d0*phi(l,:)))

  flux_S(:) = alpha(l,:)*(fluid_v(l,:)*fluid_cS(l,:) + fluid_p(l,:))


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

  else if (fluid_method=="hlle") then

     limiter = fluid_limiter

!    Reconstruct conserved quantities at cell interfaces.

     call reconstruct(fluid_cD(l,:),D_L(:),D_R(:),limiter,+1)
     call reconstruct(fluid_cE(l,:),E_L(:),E_R(:),limiter,+1)
     call reconstruct(fluid_cS(l,:),S_L(:),S_R(:),limiter,-1)

!    Reconstruct fluxes at cell interfaces.

     call reconstruct(flux_D,flux_DL,flux_DR,limiter,-1)
     call reconstruct(flux_E,flux_EL,flux_ER,limiter,-1)
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
!       fine, though it is more dissipative.

        else

           aux = 0.5d0*(alpha(l,i  )/sqrt(abs(A(l,i  )))/psi(l,i  )**2 &
                      + alpha(l,i+1)/sqrt(abs(A(l,i+1)))/psi(l,i+1)**2)

           if (shift=="none") then
              vpp = max(0.d0,+aux)
              vmm = min(0.d0,-aux)
           else
              vpp = max(0.d0,-0.5d0*(beta(l,i)+beta(l,i+1)) + aux)
              vmm = min(0.d0,-0.5d0*(beta(l,i)+beta(l,i+1)) - aux)
           end if

        end if

!       Calculate HLLE fluxes.  Notice that, even though we use the same
!       names for flux arrays, they will now represent the values
!       at cell interfaces "i+1" instrad of at the grid points "i".

        flux_D(i) = (vpp*flux_DL(i) - vmm*flux_DR(i) + vpp*vmm*(D_R(i) - D_L(i)))/(vpp - vmm)
        flux_E(i) = (vpp*flux_EL(i) - vmm*flux_ER(i) + vpp*vmm*(E_R(i) - E_L(i)))/(vpp - vmm)
        flux_S(i) = (vpp*flux_SL(i) - vmm*flux_SR(i) + vpp*vmm*(S_R(i) - S_L(i)))/(vpp - vmm)

     end do


! ************************
! ***   METHOD = MP5   ***
! ************************

! This is the MP5 method coded by Matthew Smith.

  else if (fluid_method=="mp5") then

     call mp5fluxes(flux_D,flux_E,flux_S,l)


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
     !sfluid_cD(l,i) = - idr*(flux_D(i) - flux_D(i-1))
     !sfluid_cE(l,i) = - idr*(flux_E(i) - flux_E(i-1))
     sfluid_cD(l,i) = - idr*(flux_D(i) - flux_D(i-1))/(sqrt(A(l,i))*B(l,i)*exp(6.d0*phi(l,i)))
     sfluid_cE(l,i) = - idr*(flux_E(i) - flux_E(i-1))/(sqrt(A(l,i))*B(l,i)*exp(6.d0*phi(l,i)))
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

!  sfluid_cD(l,:) = sfluid_cD(l,:) + alpha(l,:)*trK(l,:)*fluid_cD(l,:) &
!                 - alpha(l,:)*fluid_v(l,:)*fluid_cD(l,:) &
!                 *(half*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 6.d0*D1_phi(l,:) + 2.d0/r(l,:))

  sfluid_cD(l,:) = sfluid_cD(l,:) + alpha(l,:)*trK(l,:)*fluid_cD(l,:) &
                 - alpha(l,:)*fluid_v(l,:)*fluid_cD(l,:)*(2.d0/r(l,:))

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
                 !*(half*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 6.d0*D1_phi(l,:) + 2.d0/r(l,:))
                 *(2.d0/r(l,:))

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

  if (rank==0) then
     do i=1,ghost
        sfluid_cD(l,1-i) = + sfluid_cD(l,i)
        sfluid_cE(l,1-i) = + sfluid_cE(l,i)
        sfluid_cS(l,1-i) = - sfluid_cS(l,i)
     end do
  end if


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
! Not yet implemented.

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
! them in the arrays "varl" and "varr".  Notice that, while
! var(i) donotes the value of the original variable at the
! grid point "i", varl(i) and varr(i) denote the left and
! right sided reconstructions of the function evaluated at
! the cell interface correponding to the point at "i+1/2".
!
! The slope limiters work by reconstructing a function at the
! cell interface between grid points using the slopes to the
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
!     varl(i)  = var(i)   + 0.5d0*slopelim
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

  use procinfo
  use param

! Extra variables.

  implicit none

  integer i,sym

  real(8) slope1,slope2,slope3,slopelim
  real(8) ratio,phi,philimiter
  real(8) var(1-ghost:Nrmax),varl(1-ghost:Nrmax),varr(1-ghost:Nrmax)

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

  do i=0,Nr-2

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








! ***********************************
! ***   MP5 ROUTINES START HERE   ***
! ***********************************

! These routines were coded by Matthew Smith.
!
! I am not sure how all this works, I need to check it (miguel).


! **********************
! ***   MP5 FLUXES   ***
! **********************

  subroutine mp5fluxes(FD,FE,FS,l)

  use procinfo
  use arrays
  use param

  implicit none

  integer i,l

  real(8) vllf,aux

! Numerical fluxes.

  real(8) FD(1-ghost:Nrmax)
  real(8) FE(1-ghost:Nrmax)
  real(8) FS(1-ghost:Nrmax)

! Left and right conserved variables.

  real(8) D_L(1-ghost:Nrmax),D_R(1-ghost:Nrmax)
  real(8) E_L(1-ghost:Nrmax),E_R(1-ghost:Nrmax)
  real(8) S_L(1-ghost:Nrmax),S_R(1-ghost:Nrmax)
   
! Left and right fluxes.

  real(8) FD_L(1-ghost:Nrmax),FD_R(1-ghost:Nrmax)
  real(8) FE_L(1-ghost:Nrmax),FE_R(1-ghost:Nrmax)
  real(8) FS_L(1-ghost:Nrmax),FS_R(1-ghost:Nrmax)

! Left and right speeds.

  real(8) vp_L(1-ghost:Nrmax),vp_R(1-ghost:Nrmax)
  real(8) vm_L(1-ghost:Nrmax),vm_R(1-ghost:Nrmax)


! ****************************
! ***   CALCULATE FLUXES   ***
! ****************************

! Calculate fluxes at grid points.

  FD = alpha(l,:)*fluid_v(l,:)*fluid_cD(l,:) &
     *(sqrt(A(l,:))*B(l,:)*exp(6.d0*phi(l,:)))

  FE = alpha(l,:)*fluid_v(l,:)*(fluid_cE(l,:) + fluid_p(l,:) + fluid_q(l,:)) &
     *(sqrt(A(l,:))*B(l,:)*exp(6.d0*phi(l,:)))

  FS = alpha(l,:)*(fluid_v(l,:)*fluid_cS(l,:) + fluid_p(l,:))

! Reconstruct conserved quantities at cell interfaces.

  call mp5reconstruct(fluid_cD(l,:),D_L,D_R,+1)
  call mp5reconstruct(fluid_cE(l,:),E_L,E_R,+1)
  call mp5reconstruct(fluid_cS(l,:),S_L,S_R,-1)

! Reconstruct fluxes at cell interfaces.

  call mp5reconstruct(FD,FD_L,FD_R,-1)
  call mp5reconstruct(FE,FE_L,FE_R,-1)
  call mp5reconstruct(FS,FS_L,FS_R,+1)


! ***************************************
! ***   LOCAL LAX-FRIEDRICHS SPEEDS   ***
! ***************************************

! Reconstruct characteristic speeds at cell boundaries.

  if (fluid_usesoundspeed) then
     call mp5reconstruct(fluid_vcp(l,:),vp_L(:),vp_R(:),+1)
     call mp5reconstruct(fluid_vcm(l,:),vm_L(:),vm_R(:),+1)
  end if

! Find local Lax-Fredrichs speeds.

  do i=0,Nr-1

!    Use speed of sound.

     if (fluid_usesoundspeed) then

        vllf = max(abs(vp_L(i)),abs(vp_R(i)),abs(vm_L(i)),abs(vm_R(i)))

!    Use speed of light.

     else

        aux = 0.5d0*(alpha(l,i  )/sqrt(abs(A(l,i  )))/psi(l,i  )**2 &
                   + alpha(l,i+1)/sqrt(abs(A(l,i+1)))/psi(l,i+1)**2)
 
        if (shift=="none") then
           vllf = abs(aux)
        else
           vllf = max(abs(-0.5d0*(beta(l,i)+beta(l,i+1)) + aux), &
                      abs(-0.5d0*(beta(l,i)+beta(l,i+1)) - aux))
        end if

     end if

!    Calculate fluxes according to LLF interface values.

     FD(i) = 0.5d0*(FD_L(i) + FD_R(i)) - 0.5d0*vllf*(D_R(i) - D_L(i))
     FE(i) = 0.5d0*(FE_L(i) + FE_R(i)) - 0.5d0*vllf*(E_R(i) - E_L(i))
     FS(i) = 0.5d0*(FS_L(i) + FS_R(i)) - 0.5d0*vllf*(S_R(i) - S_L(i))

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine mp5fluxes







! ***************************
! ***   MP5 RECONSTRUCT   ***
! ***************************

! Performs a left and right reconstruction of a variable
! according to the mp5interface function.

  subroutine mp5reconstruct(var,varl,varr,sym)

  use procinfo
  use param
  
  implicit none

  integer i,sym

  real(8) mp5interface
  real(8) aux

  real(8) var(1-ghost:Nrmax),varl(1-ghost:Nrmax),varr(1-ghost:Nrmax)

  !! For the origin, we should be able to apply the same reconstruction
  !! as for the rest of the interior, provided there are sufficient ghost
  !! cells (test runs are using 3 ghost cells).
  !!
  !! It is questionable to reconstruct up to to index Nr-1,
  !! since our array indices only go up to Nr. 
  !! It appears however to not lead directly to crashes.
  !!
  !! For the sake of correctness, it is best to run the left- and right-
  !! reconstructions separately owing to their differing stencils.

  !! left reconstruction

  do i=0,Nr-2 
     varl(i) = mp5interface(var(i-2),var(i-1),var(i),var(i+1),var(i+2))
  end do

  ! varl(0) = mp5interface(sym*var(3), sym*var(2), sym*var(1), var(1), var(2))

  !! right reconstruction

  do i=0,Nr-3
     varr(i) = mp5interface(var(i+3),var(i+2),var(i+1),var(i),var(i-1))
  end do

  ! varr(0) = mp5interface(var(3), var(2), var(1), sym*var(1), sym*var(2))

  ! print *, var(0), "    ", var(1), "    ", sym
  ! print *, var(2), "    ", var(-1)
  !! uncomment these and everything is fine
  ! varl(0) = var(0) + sym*0.5d0*(var(2) - var(1))
  ! varr(0) = var(1) - 0.5d0*(var(2) - var(1))
  ! varl(0) = var(0) - 0.5d0*(var(0) - var(-1))
  ! varr(0) = var(1) - 0.5d0*(var(1) - var(0))
  ! varl(0) = sym*varr(0)
  ! if (sym == 1) then
     ! aux = var(2) - var(1)
     ! varr(0) = 0.0d0 !var(1)
     ! varl(0) = 0.0d0 !var(0)
  ! else
     ! varr(0) = 0.0d0
     ! varl(0) = 0.0d0
  ! endif

  !! We need a way to reconstruct at index Nr-1 and
  !! we also need the right-reconstruction at Nr-2.
  !! Choose to do non-limited reconstruction here,
  !! and replace with something more sophisticated later.
  !!
  !! This could be a source of error, but the first NaN
  !! appear far away from the outer boundary in tests.

  varl(Nr-1) = var(Nr-1) ! + 0.5d0*(var(Nr) - var(Nr-1))

  varr(Nr-2) = var(Nr-1) !var(Nr-1) - 0.5d0*(var(Nr-1) - var(Nr-2))
  varr(Nr-1) = var(Nr)   !var(Nr) - 0.5d0*(var(Nr) - var(Nr-1))
    
  end subroutine mp5reconstruct







! *************************
! ***   MP5 INTERFACE   ***
! *************************

  function mp5interface(um2, um1, u, up1, up2) result(mp5val)

  implicit none

  real(8) :: um2,um1,u,up1,up2
  real(8) :: minmod2,minmod4
  real(8) :: mp5val
  real(8) :: vnorm,vl,vmp
  real(8) :: djm1,dj,djp1,dm4jph,dm4jmh
  real(8) :: vul,vav,vmd,vlc,vmin,vmax
  real(8) :: mp5alpha,eps

  mp5alpha = 4.0d0
  eps = 1.0d-10

  vnorm = sqrt(um2*um2 + um1*um1 + u*u + up1*up1 + up2*up2)

! Not sure about this formula ...

  vl = (2.0d0*um2 - 13.0d0*um1 + 47.0d0*u + 27.0d0*up1 - 3.0d0*up2)/60.0d0
  ! vl = (3.0d0*um2 - 20.0d0*um1 + 90.0d0*u + 60.0d0*up1 - 5.0d0*up2)/128.d0

  vmp = u + minmod2(up1-u,mp5alpha*(u-um1))

  if ((vl-u)*(vl-vmp)<eps) then

     mp5val = vl

  else

     djm1 = um2 - 2.0d0*um1 + u
     dj   = um1 - 2.0d0*u + up1
     djp1 = u - 2.0d0*up1 + up2

     dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
     dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)

     vul = u + mp5alpha*(u - um1)
     vav = 0.5d0*(u + up1)
     vmd = vav - 0.5d0*dm4jph
     vlc = u + 0.5d0*(u - um1) + 4.0d0*dm4jmh/3.0d0

     vmin = max(min(u, up1, vmd), min(u, vul, vlc))
     vmax = min(max(u, up1, vmd), max(u, vul, vlc))

     mp5val = vl + minmod2(vmin - vl,vmax - vl)

  end if


! ***************
! ***   END   ***
! ***************

  end function mp5interface








! *******************
! ***   MINMOD2   ***
! *******************

  function minmod2(a,b)

  implicit none

  real(8) :: a,b
  real(8) :: minmod2
  real(8) :: signof

  minmod2 = 0.5d0*(signof(a) + signof(b))*min(abs(a),abs(b))

  end function minmod2







! *******************
! ***   MINMOD4   ***
! *******************

  function minmod4(a,b,c,d)

  implicit none

  real(8) :: a,b,c,d
  real(8) :: minmod4
  real(8) :: signof

  minmod4 = 0.125d0*(signof(a) + signof(b))*abs(signof(a) + signof(c)) &
          *abs(signof(a) + signof(d))*min(abs(a),abs(b),abs(c),abs(d))

  end function minmod4







! *******************
! ***   SIGN OF   ***
! *******************

  function signof(a)

  implicit none

  real(8) :: a
  real(8) :: signof

  if (a>0.0d0) then
     signof = +1.0d0
  else
     signof = -1.0d0
  end if

  end function signof

