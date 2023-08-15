!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/dust.f90,v 1.4 2022/07/05 16:58:16 malcubi Exp $

  subroutine sources_dust(l)

! ****************************
! ***   SOURCES FOR DUST   ***
! ****************************

! This routine calculates the sources for the relativistic
! Euler equations for dust (p=0), written in conservative form.
!
! These equations take the general form:
!
!                                                          1/2
! d D  =  beta d D  -  d [ alpha v D ]  -  alpha v D d ln g
!  t            r       r                             r
!
!      +  alpha trK D
!
!                                                          1/2
! d E  =  beta d E  -  d [ alpha v E ]  -  alpha v E d ln g
!  t            r       r                             r
!                          2    4 phi
!      +  (E + D) [ alpha v  A e     (KTA + trK/3) - v d alpha ]
!                                                       r
!
!      +  alpha trK E
!
!
! d S  =  beta d S  +  S d beta  -  d [ alpha v S ]
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
! and with W the Lorentz factor.  Notice that for dust
! the conserved energy E is essentially only the kinetic
! energy of the fluid elements.
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
  real(8) idr,aux

  real(8) flux_D(0:Nr),flux_E(0:Nr),flux_S(0:Nr)


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

! Initialize fluxes to zero.

  flux_D = 0.d0
  flux_E = 0.d0
  flux_S = 0.d0


! ******************************************
! ***   METHOD = CENTER (SECOND ORDER)   ***
! ******************************************

! The fluxes at cell interfaces are just averaged over
! adjacent grid points.  This is not very stable!

  if (dust_method=="center") then

     do i=0,Nr-1

!       Fluxes for D.

        flux_D(i) = half*(alpha(l,i  )*dust_v(l,i  )*dust_cD(l,i  ) &
                        + alpha(l,i+1)*dust_v(l,i+1)*dust_cD(l,i+1))

!       Fluxes for E.

        flux_E(i) = half*(alpha(l,i  )*dust_v(l,i  )*dust_cE(l,i  ) &
                        + alpha(l,i+1)*dust_v(l,i+1)*dust_cE(l,i+1))

!       Fluxes for S.

        flux_S(i) = half*(alpha(l,i  )*dust_v(l,i  )*dust_cS(l,i  ) &
                        + alpha(l,i+1)*dust_v(l,i+1)*dust_cS(l,i+1))

     end do


! *****************************************
! ***   METHOD = UPWIND (FIRST ORDER)   ***
! *****************************************

! The fluxes at cell interfaces are copied from
! the upwind side.

  else if (dust_method=="upwind") then

     do i=0,Nr-1

        aux = half*(dust_v(l,i) + dust_v(l,i+1))

!       Fluxes for D.

        flux_D(i) = half*(alpha(l,i  )*dust_cD(l,i  )*(one + sign(one,aux)) &
                        + alpha(l,i+1)*dust_cD(l,i+1)*(one - sign(one,aux)))*aux

!       Fluxes for E.

        flux_E(i) = half*(alpha(l,i  )*dust_cE(l,i  )*(one + sign(one,aux)) &
                        + alpha(l,i+1)*dust_cE(l,i+1)*(one - sign(one,aux)))*aux
!       Fluxes for S.

        flux_S(i) = half*(alpha(l,i  )*dust_cS(l,i  )*(one + sign(one,aux)) &
                        + alpha(l,i+1)*dust_cS(l,i+1)*(one - sign(one,aux)))*aux

     end do


! *******************************************
! ***   METHOD = LIMITER (SECOND ORDER)   ***
! *******************************************

! The fluxes at cell interfaces are either interpolated (averaged)
! or extrapolated from the upwind side using a minmod limiter.

  else if (dust_method=="limiter") then

     do i=1,Nr-2

        aux = half*(dust_v(l,i) + dust_v(l,i+1))

        if (aux>=0.d0) then

!          Fluxes for D.

           slope1 = alpha(l,i+1)*dust_cD(l,i+1) - alpha(l,i  )*dust_cD(l,i  )
           slope2 = alpha(l,i  )*dust_cD(l,i  ) - alpha(l,i-1)*dust_cD(l,i-1)

           if (slope1*slope2>0.d0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.d0
           end if

           flux_D(i) = aux*(alpha(l,i)*dust_cD(l,i) + half*slopelim)

!          Fluxes for E.

           slope1 = alpha(l,i+1)*dust_cE(l,i+1) - alpha(l,i  )*dust_cE(l,i  )
           slope2 = alpha(l,i  )*dust_cE(l,i  ) - alpha(l,i-1)*dust_cE(l,i-1)

           if (slope1*slope2>0.d0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.d0
           end if

           flux_E(i) = aux*(alpha(l,i)*dust_cE(l,i) + half*slopelim)

!          Fluxes for S.

           slope1 = alpha(l,i+1)*dust_cS(l,i+1) - alpha(l,i  )*dust_cS(l,i  )
           slope2 = alpha(l,i  )*dust_cS(l,i  ) - alpha(l,i-1)*dust_cS(l,i-1)

           if (slope1*slope2>0.d0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.d0
           end if

           flux_S(i) = aux*(alpha(l,i)*dust_cS(l,i) + half*slopelim)

        else

!          Fluxes for D.

           slope1 = alpha(l,i+1)*dust_cD(l,i+1) - alpha(l,i  )*dust_cD(l,i  )
           slope2 = alpha(l,i+2)*dust_cD(l,i+2) - alpha(l,i+1)*dust_cD(l,i+1)

           if (slope1*slope2>0.d0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.d0
           end if

           flux_D(i) = aux*(alpha(l,i+1)*dust_cD(l,i+1) - half*slopelim)

!          Fluxes for E.

           slope1 = alpha(l,i+1)*dust_cE(l,i+1) - alpha(l,i  )*dust_cE(l,i  )
           slope2 = alpha(l,i+2)*dust_cE(l,i+2) - alpha(l,i+1)*dust_cE(l,i+1)

           if (slope1*slope2>0.d0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.d0
           end if

           flux_E(i) = aux*(alpha(l,i+1)*dust_cE(l,i+1) - half*slopelim)

!          Fluxes for S.

           slope1 = alpha(l,i+1)*dust_cS(l,i+1) - alpha(l,i  )*dust_cS(l,i  )
           slope2 = alpha(l,i+2)*dust_cS(l,i+2) - alpha(l,i+1)*dust_cS(l,i+1)

           if (slope1*slope2>0.d0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.d0
           end if

           flux_S(i) = aux*(alpha(l,i+1)*dust_cS(l,i+1) - half*slopelim)

        end if

     end do


! **************************
! ***   UNKNOWN METHOD   ***
! **************************

! This should never happen, but just in case.

  else

     print *, 'Unknown dust_method ...'
     print *, 'Aborting! (subroutine dust.f90)'
     print *
     call die

  end if


! *******************************************
! ***   SOURCE CONTRIBUTION FROM FLUXES   ***
! *******************************************

  do i=1,Nr-1
     sdust_cD(l,i) = - (flux_D(i) - flux_D(i-1))*idr
     sdust_cE(l,i) = - (flux_E(i) - flux_E(i-1))*idr
     sdust_cS(l,i) = - (flux_S(i) - flux_S(i-1))*idr
  end do


! **************************************
! ***   ADD ALL OTHER SOURCE TERMS   ***
! **************************************

! Here we add all other source terms, including terms related to
! the extrinsic curvature, the Christoffel symbols and the shift.

! Sources for D.

  sdust_cD(l,:) = sdust_cD(l,:) + alpha(l,:)*trK(l,:)*dust_cD(l,:) &
                - alpha(l,:)*dust_v(l,:)*dust_cD(l,:) &
                *(half*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 6.d0*D1_phi(l,:) + 2.d0/r(l,:))

  if (shift/="none") then
     sdust_cD(l,:) = sdust_cD(l,:) + beta(l,:)*DA_dust_cD(l,:)
  end if

! Sources for E.

  sdust_cE(l,:) = sdust_cE(l,:) + alpha(l,:)*trK(l,:)*dust_cE(l,:) &
                + (dust_cE(l,:) + dust_cD(l,:)) &
                *(alpha(l,:)*dust_v(l,:)**2*A(l,:)*psi4(l,:)*(KTA(l,:) + third*trK(l,:)) &
                - dust_v(l,:)*D1_alpha(l,:)) - alpha(l,:)*dust_v(l,:)*dust_cE(l,:) &
                *(half*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) + 6.d0*D1_phi(l,:) + 2.d0/r(l,:))

  if (shift/="none") then
     sdust_cE(l,:) = sdust_cE(l,:) + beta(l,:)*DA_dust_cE(l,:)
  end if

! Sources for S.

  sdust_cS(l,:) = sdust_cS(l,:) + alpha(l,:)*trK(l,:)*dust_cS(l,:) &
                - (dust_cE(l,:) + dust_cD(l,:))*D1_alpha(l,:) &
                - alpha(l,:)*dust_v(l,:)*dust_cS(l,:) &
                *(D1_B(l,:)/B(l,:) + 4.d0*D1_phi(l,:) + 2.d0/r(l,:))

  if (shift/="none") then
     sdust_cS(l,:) = sdust_cS(l,:) + beta(l,:)*DA_dust_cS(l,:) &
                   + dust_cS(l,:)*D1_beta(l,:)
  end if


! *************************************
! ***   SOURCES AT OUTER BOUNDARY   ***
! *************************************

! For the moment I just set the sources at the
! boundaries to zero.  For fine grid this should
! be corrected when we restrict from the coarse
! grids.

  sdust_cD(l,Nr) = 0.d0
  sdust_cE(l,Nr) = 0.d0
  sdust_cS(l,Nr) = 0.d0


! ************************
! ***   GHOST POINTS   ***
! ************************

! Ghost points using symmetries.

  do i=1,ghost
     sdust_cD(l,1-i) = + sdust_cD(l,i)
     sdust_cE(l,1-i) = + sdust_cE(l,i)
     sdust_cS(l,1-i) = - sdust_cS(l,i)
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine sources_dust


