!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/fluid.f90,v 1.14 2024/12/04 16:54:44 malcubi Exp $

  subroutine sources_fluid(l)

! *****************************
! ***   SOURCES FOR FLUID   ***
! *****************************

! This routine calculates the sources for the relativistic Euler
! equations for a perfect fluid, written in conservative form.
!
! These equations take the general form:
!
!                                                          1/2
! d D  =  beta d D  -  d [ alpha v D ]  -  alpha v D d ln g
!  t            r       r                             r
!
!      +  alpha trK D
!
!                                                                      1/2
! d E  =  beta d E  -  d [ alpha v (E + p) ]  -  alpha v (E + p) d ln g
!  t            r       r                                         r
!                              2    4 phi
!      +  (E + p + D) [ alpha v  A e     (KTA + trK/3) - v d alpha ]
!                                                           r
!
!      +  alpha trK (E + p)
!
!
! d S  =  beta d S  +  S d beta  -  d [ alpha ( v S + p ) ]
!  t            r         r          r
!
!      -  alpha v S [ d B / B  +  4 d phi  +  2/r ]
!                      r             r
!
!      -  (E + D) d alpha  +  alpha trK S
!                  r
!
! where g is the determinant of the spatial metric and:
!
!       1/2
! d ln g   = d g / 2g = d A / 2A + d B / B + 6 d phi + 2/r
!  r          r          r          r           r
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


! ******************************************
! ***   METHOD = CENTER (SECOND ORDER)   ***
! ******************************************

! The fluxes at cell interfaces are just averaged over
! adjacent grid points.  This is not very stable!

  if (fluid_method=="center") then

!    NOT YET IMPLEMENTED.

     print *, 'fluid_method=center not yet implemented.'
     print *, 'Aborting! (subroutine fluid.f90)'
     print *

     call die


! *****************************************
! ***   METHOD = UPWIND (FIRST ORDER)   ***
! *****************************************

! The fluxes at cell interfaces are copied from
! the upwind side.
!
! Notice that this is not a "true" upwind since
! I don't decompose into eigenfunctions, and I
! just base the direction in which the derivatives
! are calculated on the fluid speed and not the
! characteristic speeds.
!
! This method is not very good, it is only first
! order accurate, and also it tends to go unstable
! at the origin and develop high frequency noise
! when the fluid speed approaches zero.  It is only
! here as a quick and dirty way to test things.

  else if (fluid_method=="upwind") then

     do i=0,Nr-1

        aux = half*(fluid_v(l,i) + fluid_v(l,i+1))

!       Fluxes for D.

        flux_D(i) = half*(alpha(l,i  )*fluid_cD(l,i  )*(one + sign(one,aux)) &
                        + alpha(l,i+1)*fluid_cD(l,i+1)*(one - sign(one,aux)))*aux

!       Fluxes for E.

        flux_E(i) = half*(alpha(l,i  )*fluid_cE(l,i  )*(one + sign(one,aux)) &
                        + alpha(l,i+1)*fluid_cE(l,i+1)*(one - sign(one,aux)))*aux &
                  + half*(alpha(l,i  )*fluid_v(l,i  )*(fluid_p(l,i  ) + fluid_q(l,i  )) &
                        + alpha(l,i+1)*fluid_v(l,i+1)*(fluid_p(l,i+1) + fluid_q(l,i+1))) 

!       Fluxes for S.

        flux_S(i) = half*(alpha(l,i  )*fluid_cS(l,i  )*(one + sign(one,aux)) &
                        + alpha(l,i+1)*fluid_cS(l,i+1)*(one - sign(one,aux)))*aux &
                  + half*(alpha(l,i  )*(fluid_p(l,i  ) + fluid_q(l,i  )) &
                        + alpha(l,i+1)*(fluid_p(l,i+1) + fluid_q(l,i+1)))

     end do


! *******************************************
! ***   METHOD = LIMITER (SECOND ORDER)   ***
! *******************************************

! The fluxes at cell interfaces are either interpolated (averaged)
! or extrapolated from the upwind side using a slope limiter.
!
! Just as the "uwpind" method above, we only take into account the
! fluid speed and not the characteristic speeds.  This method is
! second order, but it is not very accurate, nor very stable.
! It is only here as a quick and dirty way to test things.

  else if (fluid_method=="limiter") then

     do i=0,Nr-1

        aux = half*(fluid_v(l,i) + fluid_v(l,i+1))

!       Fluxes for D.

        if (aux>=0.d0) then

           slope1 = alpha(l,i+1)*fluid_cD(l,i+1) - alpha(l,i  )*fluid_cD(l,i  )
           slope2 = alpha(l,i  )*fluid_cD(l,i  ) - alpha(l,i-1)*fluid_cD(l,i-1)

           if (slope1*slope2>0.d0) then
              if (fluid_limiter=="minmod") then
                 if (abs(slope1)<abs(slope2)) then
                    slopelim = slope1
                 else
                    slopelim = slope2
                 end if
              else if (fluid_limiter=="vanleer") then
                 slopelim = 2.d0*slope1*slope2/(slope1 + slope2)
              end if
           else
              slopelim = 0.d0
           end if

           flux_D(i) = aux*(alpha(l,i  )*fluid_cD(l,i  ) + half*slopelim)

        else

           slope1 = alpha(l,i+1)*fluid_cD(l,i+1) - alpha(l,i  )*fluid_cD(l,i  )
           slope2 = alpha(l,i+2)*fluid_cD(l,i+2) - alpha(l,i+1)*fluid_cD(l,i+1)

           if (slope1*slope2>0.d0) then
              if (fluid_limiter=="minmod") then
                 if (abs(slope1)<abs(slope2)) then
                    slopelim = slope1
                 else
                    slopelim = slope2
                 end if
              else if (fluid_limiter=="vanleer") then
                 slopelim = 2.d0*slope1*slope2/(slope1 + slope2)
              end if
           else
              slopelim = 0.d0
           end if

           flux_D(i) = aux*(alpha(l,i+1)*fluid_cD(l,i+1) - half*slopelim)

        end if

!       Fluxes for E.

        if (aux>=0.d0) then

           slope1 = alpha(l,i+1)*(fluid_cE(l,i+1) + fluid_p(l,i+1) + fluid_q(l,i+1)) &
                  - alpha(l,i  )*(fluid_cE(l,i  ) + fluid_p(l,i  ) + fluid_q(l,i  ))
           slope2 = alpha(l,i  )*(fluid_cE(l,i  ) + fluid_p(l,i  ) + fluid_q(l,i  )) &
                  - alpha(l,i-1)*(fluid_cE(l,i-1) + fluid_p(l,i-1) + fluid_q(l,i-1))

           if (slope1*slope2>0.d0) then
              if (fluid_limiter=="minmod") then
                 if (abs(slope1)<abs(slope2)) then
                    slopelim = slope1
                 else
                    slopelim = slope2
                 end if
              else if (fluid_limiter=="vanleer") then
                 slopelim = 2.d0*slope1*slope2/(slope1 + slope2)
              end if
           else
              slopelim = 0.d0
           end if

           flux_E(i) = aux*(alpha(l,i  )*(fluid_cE(l,i  ) + fluid_p(l,i  ) + fluid_q(l,i  )) &
                     + half*slopelim)

        else

           slope1 = alpha(l,i+1)*(fluid_cE(l,i+1) + fluid_p(l,i+1) + fluid_q(l,i+1)) &
                  - alpha(l,i  )*(fluid_cE(l,i  ) + fluid_p(l,i  ) + fluid_q(l,i  ))
           slope2 = alpha(l,i+2)*(fluid_cE(l,i+2) + fluid_p(l,i+2) + fluid_q(l,i+2)) &
                  - alpha(l,i+1)*(fluid_cE(l,i+1) + fluid_p(l,i+1) + fluid_q(l,i+1))

           if (slope1*slope2>0.d0) then
              if (fluid_limiter=="minmod") then
                 if (abs(slope1)<abs(slope2)) then
                    slopelim = slope1
                 else
                    slopelim = slope2
                 end if
              else if (fluid_limiter=="vanleer") then
                 slopelim = 2.d0*slope1*slope2/(slope1 + slope2)
              end if
           else
              slopelim = 0.d0
           end if

           flux_E(i) = aux*(alpha(l,i+1)*(fluid_cE(l,i+1) + fluid_p(l,i+1) + fluid_q(l,i+1)) &
                     - half*slopelim)

        end if

!       Fluxes for S.

        if (aux>=0.d0) then

           slope1 = alpha(l,i+1)*fluid_cS(l,i+1) - alpha(l,i  )*fluid_cS(l,i  )
           slope2 = alpha(l,i  )*fluid_cS(l,i  ) - alpha(l,i-1)*fluid_cS(l,i-1)

           if (slope1*slope2>0.d0) then
              if (fluid_limiter=="minmod") then
                 if (abs(slope1)<abs(slope2)) then
                    slopelim = slope1
                 else
                    slopelim = slope2
                 end if
              else if (fluid_limiter=="vanleer") then
                 slopelim = 2.d0*slope1*slope2/(slope1 + slope2)
              end if
           else
              slopelim = 0.d0
           end if

           flux_S(i) = aux*(alpha(l,i)*fluid_cS(l,i) + half*slopelim)

        else

           slope1 = alpha(l,i+1)*fluid_cS(l,i+1) - alpha(l,i  )*fluid_cS(l,i  )
           slope2 = alpha(l,i+2)*fluid_cS(l,i+2) - alpha(l,i+1)*fluid_cS(l,i+1)

           if (slope1*slope2>0.d0) then
              if (fluid_limiter=="minmod") then
                 if (abs(slope1)<abs(slope2)) then
                    slopelim = slope1
                 else
                    slopelim = slope2
                 end if
              else if (fluid_limiter=="vanleer") then
                 slopelim = 2.d0*slope1*slope2/(slope1 + slope2)
              end if
           else
              slopelim = 0.d0
           end if

           flux_S(i) = aux*(alpha(l,i+1)*fluid_cS(l,i+1) - half*slopelim)

        end if

!       Add contribution from pressure (it is independent of the speed v).

        flux_S(i) = flux_S(i) &
                  + half*(alpha(l,i  )*(fluid_p(l,i  ) + fluid_q(l,i  )) &
                        + alpha(l,i+1)*(fluid_p(l,i+1) + fluid_q(l,i+1)))

     end do


! **********************
! ***   LLF METHOD   ***
! **********************

! Local Lax-Friedrichs method.  This is only first order.

  else if (fluid_method=="llf") then

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

     flux_S = alpha(l,:)*(fluid_v(l,:)*fluid_cS(l,:) + fluid_p(l,:) + fluid_q(l,:))
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
     sfluid_cD(l,i) = - (flux_D(i) - flux_D(i-1))*idr
     sfluid_cE(l,i) = - (flux_E(i) - flux_E(i-1))*idr
     sfluid_cS(l,i) = - (flux_S(i) - flux_S(i-1))*idr
  end do


! **************************************
! ***   ADD ALL OTHER SOURCE TERMS   ***
! **************************************

! Here we add all other source terms, including terms related to
! the extrinsic curvature, the Christoffel symbols and the shift.

! Sources for D.

  sfluid_cD(l,:) = sfluid_cD(l,:) + alpha(l,:)*trK(l,:)*fluid_cD(l,:) &
                 - alpha(l,:)*fluid_v(l,:)*fluid_cD(l,:) &
                 *(half*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 6.d0*D1_phi(l,:) + 2.d0/r(l,:))

  if (shift/="none") then
     sfluid_cD(l,:) = sfluid_cD(l,:) + beta(l,:)*DA_fluid_cD(l,:)
  end if

! Sources for E.

  sfluid_cE(l,:) = sfluid_cE(l,:) + (fluid_cE(l,:) + fluid_cD(l,:) + fluid_p(l,:) + fluid_q(l,:)) &
                 *(alpha(l,:)*fluid_v(l,:)**2*A(l,:)*exp(4.d0*phi(l,:)) &
                 *(KTA(l,:) + third*trK(l,:)) - fluid_v(l,:)*D1_alpha(l,:)) &
                 + alpha(l,:)*trK(l,:)*(fluid_cE(l,:) + fluid_p(l,:) + fluid_q(l,:)) &
                 - alpha(l,:)*fluid_v(l,:)*(fluid_cE(l,:) + fluid_p(l,:) + fluid_q(l,:)) &
                 *(half*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 6.d0*D1_phi(l,:) + 2.d0/r(l,:))

  if (shift/="none") then
     sfluid_cE(l,:) = sfluid_cE(l,:) + beta(l,:)*DA_fluid_cE(l,:)
  end if

! Sources for S.

  sfluid_cS(l,:) = sfluid_cS(l,:) + alpha(l,:)*trK(l,:)*fluid_cS(l,:) &
                 - (fluid_cE(l,:) + fluid_cD(l,:))*D1_alpha(l,:) &
                 - alpha(l,:)*fluid_v(l,:)*fluid_cS(l,:) &
                 *(D1_B(l,:)/B(l,:) + 4.d0*D1_phi(l,:) + 2.d0/r(l,:))

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
! extrapolation using a slope limiter (minmod or vanleer).  
! We do both left and right sided reconstructions, and store
! them in the arrays "varl" and "varr".

! Include modules.

  use param

! Extra variables.

  implicit none

  integer i,sym

  real(8) slope1,slope2,slope3,slopelim
  real(8) var(1-ghost:Nr),varl(1-ghost:Nr),varr(1-ghost:Nr)

  character(len=*) limiter


! ***********************************************
! ***   FIND LEFT AND RIGHT RECONSTRUCTIONS   ***
! ***********************************************

! Sanity check.

  if ((limiter/="minmod").and.(limiter/="vanleer")) then
     print *
     print *, 'Unknown limiter type in fluid.f90, aborting ...'
     print *
     call die
  end if

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
        if (limiter=="minmod") then
           if (abs(slope1)<abs(slope2)) then
              slopelim = slope1
           else
              slopelim = slope2
           end if
        else if (limiter=="vanleer") then
           slopelim = 2.d0*slope1*slope2/(slope1 + slope2)
        end if
     else
        slopelim = 0.d0
     end if

     varl(i) = var(i) + 0.5d0*slopelim

!    Right side reconstruction.

     if (slope2*slope3>0.d0) then
        if (limiter=="minmod") then
           if (abs(slope2)<abs(slope3)) then
              slopelim = slope2
           else
              slopelim = slope3
           end if
        else if (limiter=="vanleer") then
           slopelim = 2.d0*slope2*slope3/(slope2 + slope3)
        end if
     else
        slopelim = 0.d0
     end if

     varr(i) = var(i+1) - 0.5d0*slopelim

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
     if (limiter=="minmod") then
        if (abs(slope1)<abs(slope2)) then
           slopelim = slope1
        else
           slopelim = slope2
        end if
     else if (limiter=="vanleer") then
        slopelim = 2.d0*slope1*slope2/(slope1 + slope2)
     end if
  else
     slopelim = 0.d0
  end if

  varl(0) = var(0) + 0.5d0*slopelim

! Right side reconstruction.

  if (slope2*slope3>0.d0) then
     if (limiter=="minmod") then
        if (abs(slope2)<abs(slope3)) then
           slopelim = slope2
        else
           slopelim = slope3
        end if
     else if (limiter=="vanleer") then
        slopelim = 2.d0*slope2*slope3/(slope2 + slope3)
     end if
  else
     slopelim = 0.d0
  end if

  varr(0) = var(1) - 0.5d0*slopelim


! **************************
! ***   OUTER BOUNDARY   ***
! **************************

! Find slopes.  In this case there is no slope3.

  slope1 = var(Nr-1) - var(Nr-2)
  slope2 = var(Nr  ) - var(Nr-1)

! Left side reconstruction.

  if (slope1*slope2>0.d0) then
     if (limiter=="minmod") then
        if (abs(slope1)<abs(slope2)) then
           slopelim = slope1
        else
           slopelim = slope2
        end if
     else if (limiter=="vanleer") then
        slopelim = 2.d0*slope1*slope2/(slope1 + slope2)
     end if
  else
     slopelim = 0.d0
  end if

  varl(Nr-1) = var(Nr-1) + 0.5d0*slopelim

! Right side reconstruction.

  varr(Nr-1) = var(Nr  ) - 0.5d0*slope2


! ***************
! ***   END   ***
! ***************

  end subroutine reconstruct


