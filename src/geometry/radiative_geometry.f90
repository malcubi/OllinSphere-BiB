!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/radiative_geometry.f90,v 1.49 2023/10/04 19:42:33 malcubi Exp $

  subroutine radiative_geometry(l)

! ****************************************
! ***   RADIATIVE BOUNDARY CONDITION   ***
! ****************************************

! Radiative boundaries for geometry.
!
! The routine radiative_boundary applies outgoing wave boundary
! conditions to the geometric variables (trK,KTA,Klambda,dtbeta).
! These boundary conditions take into account the characteristic
! structure, but not the constraints.
!
! Notice that the boundary conditions here are applied at the level
! of the source terms, and not directly to the dynamical quantities
! themselves.
!
! For the metric variables (alpha,beta,phi,A,B,lambda), which are
! declared with the attribute NOBOUND, we calculate the sources
! all the way up to the boundary using one-sided differences.
!
! We treat "Deltar" just as the metric variables, calculating the
! sources all the way to the boundaries, which seems to work fine
! when the shift is not dynamical.  For the Gammadriver2 shift
! we apply an outgoing wave boundary condition to "dbeta".
!
! NOTE: Even though the grid level l is passed as an argument,
! the routine is really only ever called for l=0.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer l

  real(8) zero,half,third
  real(8) vl,va,vs
  real(8) aux


! **********************************
! ***   DO WE WANT TO BE HERE?   ***
! **********************************

  if (spacetime/="dynamic") return


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  half  = 0.5d0
  third = 1.d0/3.d0


! **********************************
! ***   ONLY PHYSICAL BOUNDARY   ***
! **********************************

! Only the processor at the physical boundary applies boundary conditions.

  if ((rank==size-1).and.(l==0)) then


!    ***********************
!    ***   EIGENSPEEDS   ***
!    ***********************

!    Find asymptotic eigenspeeds:
!
!    vl = speed of light
!    va = slicing gauge speed
!    vs = shift gauge speed.

!    Notice that here we assume that we are using the
!    Lagrangian formulation.  For the Eulerian formulation
!    the characteristic analysis becomes much more difficult.
!
!    Also, to calculate these speeds we approximate alpha=A=psi=1
!    far away, so the boundary conditions should really be applied
!    quite far away.

     vl = 1.d0
     va = dsqrt(gauge_f)
     vs = dsqrt(4.d0*third*drivercsi)


!    ****************************
!    ***   BOUNDARY FOR TRK   ***
!    ****************************

!    Radiative boundary for trK (this is only needed for
!    Bona-Masso type slicing conditions).
!
!    In order to impose the boundary condition we assume that
!    far away trK behaves as an outgoing spherical wave:
!
!    trK  ~  g(r-va*t)/r
!
!    with "va" the slicing gauge speed.
!
!    This can be shown to imply:
!
!    strK  ~  - va (D1_trK + trK/r)

     if (slicing=="maximal") then

        strK(l,Nr) = zero

     else

        if (.not.cosmic_run) then

!          For asymptotically flat spacetimes I take the gauge
!          speed at the boundary equal to sqrt(gauge_f), that is,
!          I assume the metric is essentially Minkowski far away.

           strK(l,Nr) = - va*(D1_trK(l,Nr) + trK(l,Nr)/r(l,Nr))

        else

!          For cosmological spacetimes, the gauge speed at the
!          boundary is given in terms of the background metric.
!          We also assume that trK evolves as the background
!          plus an outgoing wave.

           if (index(adjustl(slicing),"cosmo")==1) then
              aux = sqrt(gauge_f)*cosmobg_alpha(l)/cosmobg_a(l)
           else
              aux = sqrt(cosmobg_falpha(l))/cosmobg_a(l)
           end if

           strK(l,Nr) = scosmobg_trK(l) - aux*(D1_trK(l,Nr) + (trK(l,Nr) - cosmobg_trK(l))/r(l,Nr))

!          For cosmological spacetimes the above condition causes
!          causes a drift at late times that can become large for
!          very long evolutions.  Here I add a term to correct this,
!          but it spoils the outgoing wave condition, so I only
!          turn it on after the main pulse has gone through the
!          boundary.

           if (t(0)>sqrt(gauge_f)*r(0,Nr)) then
              strK(l,Nr) = strK(l,Nr) - 0.01d0*(trK(l,Nr) - cosmobg_trK(l))
           end if

        end if

     end if


!    ****************************
!    ***   BOUNDARY FOR KTA   ***
!    ****************************

!    Radiative boundary for KTA.  Notice that this condition works
!    reasonably well, but it will violate the constraints at the
!    boundary.
!
!    Now, for the eigenfields that move at the speed of
!    light we have:
!
!    wlp + wlm  ~  KTA  -  2/3  trK
!
!    If we now assume that far away these eigenfields 
!    also behave as an outgoing spherical wave, but
!    with vl, we will have:
!
!    wlp + wlm  ~  h(r-vl*t)/r
!
!    Both these assumptions together can be shown to imply that:
!
!    sKTA  ~  - vl (D1_KTA + KTA/r) + 2/3 (1-vl/va) strK
!
!    Do notice that if va=vl then the last term vanishes.

     if (boundtype/="constraint") then

        if (.not.cosmic_run) then

!          For asymptotically flat spacetimes I take the speed
!          of light at the boundary equal to 1, that is, I assume
!          the metric is essentially Minkowski far away.

           sKTA(l,Nr) = - vl*(D1_KTA(l,Nr) + KTA(l,Nr)/r(l,Nr)) &
                      + 2.d0*third*(1.d0 - 1.d0/va)*strK(l,Nr)

        else

!          For cosmological spacetimes, the speed of light at the
!          boundary is given in terms of the background metric.

           aux = cosmobg_alpha(l)/cosmobg_a(l)

           if (index(adjustl(slicing),"cosmo")==1) then
              sKTA(l,Nr) = - aux*(D1_KTA(l,Nr) + KTA(l,Nr)/r(l,Nr)) &
                         + 2.d0*third*(1.d0 - 1.d0/va)*(strK(l,Nr) - scosmobg_trK(l))
           else
              sKTA(l,Nr) = - aux*(D1_KTA(l,Nr) + KTA(l,Nr)/r(l,Nr)) &
                         + 2.d0*third*(1.d0 - cosmobg_alpha(l)/sqrt(cosmobg_falpha(l)))*(strK(l,Nr)-scosmobg_trK(l))
           end if

        end if

!       Having found the source for KTA at the boundary,
!       we now recover the source for Klambda.

        if (.not.nolambda) then
           sKlambda(l,Nr) = 1.5d0*sKTA(l,Nr)/r(l,Nr)**2
           if (regular2) then
              sKlambda2(l,Nr) = (sKlambda(l,Nr) - 4.d0*Klambda(l,Nr)*sphi(l,Nr))/psi4(l,Nr)
           end if
        end if

     end if


!    ************************************
!    ***   BOUNDARY FOR DELTA/DTBETA  ***
!    ************************************

!    Radiative boundary for dtbeta.
!
!    At the moment, we only apply a boundary condition
!    for the Gammadriver shifts.

     if ((shift(1:11)=="Gammadriver").and.(drivercsi/=0.d0)) then

!       Gammadriver1.

        if (shift=="Gammadriver1") then

!          This is the same as for Gammadriver2 below, but we apply the condition
!          to Deltar instead of dtbeta (remember that dtbeta ~ csi Deltar).

           if (.not.cosmic_run) then
              sDeltar(l,Nr) = - vs*(D1_Deltar(l,Nr) + Deltar(l,Nr)/r(l,Nr)) &
                 - vs**2*va/(vs + va)*(D1_trK(l,Nr) + trK(l,Nr)/r(l,Nr))/drivercsi
           else
              sDeltar(l,Nr) = - vs*(D1_Deltar(l,Nr) + Deltar(l,Nr)/r(l,Nr)) &
                 - vs**2*va/(vs + va)*(D1_trK(l,Nr) + (trK(l,Nr)-cosmobg_trK(l))/r(l,Nr))/drivercsi
           end if

!       Gammadriver2.

        else if (shift=="Gammadriver2") then

!          Assuming that the spacetime is asymptotically flat at the
!          boundary (and the shift is small), the sum of ingoing
!          and outgoing shift eigenfields at the boundary is:
!
!          wbp  +  wbm  ~  dtbeta  +  Q  D1_alpha
!
!          with:
!
!          Q  =  vs**2 / ( vs**2 - va**2 )
!
!          with "vs" the shift eigenspeed and "va" the slicing eigenspeed.
!          Assuming now that far away we have:
!
!          wbp  +  wbm  ~  h(r-vs*t)/r
!
!          we find that:
!
!          sdtbeta  ~  - vs (D1_dtbeta + dtbeta/r)
!
!                   +  Q (vs - va) (D2_alpha + D1_alpha/r)
!
!                   =  - vs (D1_dtbeta + dtbeta/r)
!
!                   -  vs**2 / (vs + va) (D2_alpha + D1_alpha/r)
!
!          This expression works fine in many cases, but causes a drift
!          when the final lapse is non-trivial. We can find a better
!          expression  if we assume that the incoming slicing eigenfield
!          is very small, which implies:  D1_alpha ~ sqrt(gauge_f) trK.
!
!          Notice that for travelling waves both expressions should give
!          similar results, but for stationary non-trivial solutions the
!          spatial derivatives of the lapse will not vanish, while trK will,
!          so the second expression should produce a much smaller late-time
!          drift (that in fact should vanish for dtbeta~1/r).

!          sdtbeta(l,Nr) = - vs*(D1_dtbeta(l,Nr) + dtbeta(l,Nr)/r(l,Nr)) &
!             - v2**2/(vs + va)*(D2_alpha(l,Nr) + D1_alpha(l,Nr)/r(l,Nr))

           if (.not.cosmic_run) then
              sdtbeta(l,Nr) = - vs*(D1_dtbeta(l,Nr) + dtbeta(l,Nr)/r(l,Nr)) &
                 - vs**2*va/(vs + va)*(D1_trK(l,Nr) + trK(l,Nr)/r(l,Nr))
           else
              sdtbeta(l,Nr) = - vs*(D1_dtbeta(l,Nr) + dtbeta(l,Nr)/r(l,Nr)) &
                 - vs**2*va/(vs + va)*(D1_trK(l,Nr) + (trK(l,Nr)-cosmobg_trK(l))/r(l,Nr))
           end if

!       Gammadriver3.

        else if (shift=="Gammadriver3") then

!          I still need to explain this ...

!          sdtbeta(l,Nr) = - vs*(D1_dtbeta(l,Nr) + dtbeta(l,Nr)/r(l,Nr)) &
!             - (vs-va)*(D2_alpha(l,Nr) + D1_alpha(l,Nr)/r(l,Nr))

           sdtbeta(l,Nr) = - vs*(D1_dtbeta(l,Nr) + dtbeta(l,Nr)/r(l,Nr)) &
              - va*(vs - va)*(D1_trK(l,Nr) + trK(l,Nr)/r(l,Nr))

        end if

     end if


!    **********************************
!    ***   BOUNDARY FOR Z4C THETA   ***
!    **********************************

!    For the Z4c formulation we need to apply a
!    boundary condition to z4theta.
!
!    In order to impose the boundary condition we assume that
!    far away trK behaves as an outgoing spherical wave:
!
!    Theta  ~  g(r-t)/r
!
!    This can be shown to imply:
!
!    sTheta  ~  - (D1_Theta + Theta/r)

     if (formulation=="z4c") then
        sz4theta(l,Nr) = - (D1_z4theta(l,Nr) + z4theta(l,Nr)/r(l,Nr))
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine radiative_geometry

