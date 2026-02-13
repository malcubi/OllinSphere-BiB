!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/sources_shift.f90,v 1.38 2026/02/10 20:05:44 malcubi Exp $

  subroutine sources_shift(l)

! *****************************
! ***   SOURCES FOR SHIFT   ***
! *****************************

! This routine calculates the sources for the shift.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer l

  real(8) sigma
  real(8) zero,two,half,third


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  two   = 2.d0
  half  = 0.5d0
  third = 1.d0/3.d0

  if (bssnflavor=='eulerian') then
     sigma = 0.d0
  else if (bssnflavor=='lagrangian') then
     sigma = 1.d0
  end if


! ************************
! ***   STATIC SHIFT   ***
! ************************

  if ((shift=="zero").or.(shift=="static")) then

     sbeta(l,:)   = zero
     sdtbeta(l,:) = zero


! *************************************************
! ***   GAMMADRIVER0:  PARABOLIC GAMMA DRIVER   ***
! *************************************************

! This is the parabolic Gamma driver:
!
! d beta  =  csi d Deltar
!  t              t

  else if (shift=="Gammadriver0") then

!    Source for shift.

     sbeta(l,:) = drivercsi*sDeltar(l,:)

!    We are not evolving dtbeta.

     sdtbeta(l,:) = 0.d0


! **************************************************************
! ***   GAMMADRIVER1:  FIRST ORDER HYPERBOLIC GAMMA DRIVER   ***
! **************************************************************

! This is the first order hyperbolic Gamma driver:
!
! d beta  =  beta d beta  +  csi Deltar
!  t               r
!
! Where "csi" is the parameter that fixes the wave speed
! (csi=0.75 corresponds to a speed of 1 in the asymptotic region).
!
! I have introduced a modifications that helps avoid
! coordinate shocks:
!
!
! d beta  =  beta d beta  +  csi fdriver Deltar
!  t               r
!
! In this case "fdriver" is a dynamical function that evolves through:
!
! d fdriver =  beta d fdriver  -  (2/3) fdriver DIV.beta
!  t                 r

  else if ((shift=="Gammadriver1").or.(shift=="Gammadrivershock1")) then

     if (shift=="Gammadriver1") then

!       Source for shift.

        if (.not.cosmic_run) then

           sbeta(l,:) = drivercsi*Deltar(l,:)

           if (driveradv) then
              sbeta(l,:) = sbeta(l,:) + beta(l,:)*DA_beta(l,:)
           end if

        else

           sbeta(l,:) = drivercsi*Deltar(l,:)/cosmobg_a(l)**2

        end if

        if (driverD0) then
           if (.not.cosmic_run) then
              sbeta(l,:) = sbeta(l,:) - drivercsi*Deltar0(l,:)
           else
              sbeta(l,:) = sbeta(l,:) - drivercsi*Deltar0(l,:)/cosmobg_a(l)**2
           end if
        end if

     else if (shift=="Gammadrivershock1") then

!       Source for fdriver.

        sfdriver(l,:) = - two*third*fdriver(l,:)*DIV_beta(l,:)

        if (driveradv) then
           sfdriver(l,:) = sfdriver(l,:) + beta(l,:)*DA_fdriver(l,:)
        end if

!       Source for shift.

        sbeta(l,:) = drivercsi*fdriver(l,:)*Deltar(l,:)

        if (driveradv) then
           sbeta(l,:) = sbeta(l,:) + beta(l,:)*DA_beta(l,:)
        end if

        if (driverD0) then
           sbeta(l,:) = sbeta(l,:) - drivercsi*fdriver(l,:)*Deltar0(l,:)
        end if

     end if

!    Damping term.

     sbeta(l,:) = sbeta(l,:) - drivereta*beta(l,:)

!    We are not evolving dtbeta.

     sdtbeta(l,:) = zero


! ***************************************************************
! ***   GAMMADRIVER2:  SECOND ORDER HYPERBOLIC GAMMA DRIVER   ***
! ***************************************************************

! This is the standard second order hyperbolic Gamma driver used by
! everybody written in first order form by defining an extra variable
! "dtbeta" as the time derivative of beta:
!
! d beta   =  beta d beta  +  dtbeta
!  t                r
!
! d dtbeta =  csi d Deltar  -  eta d beta
!  t               t                t
!
! As before, "csi" is the parameter that fixes the wave speed
! (csi=0.75 corresponds to a speed of 1 in the asymptotic region),
! and "eta" is a damping parameter.  Notice that eta has dimensions
! of 1/distance, so that its value has to be changed depending
! on the scale of the problem.
!
! I have introduced a modification that helps avoid
! coordinate shocks:
!
!
! d dtbeta  =  csi d (fdriver Deltar)  -  eta d beta
!  t                t                          t
!
! For this case "fdriver" is a dynamical function that evolves through:
!
! d fdriver =  beta d fdriver  -  (2/3) fdriver DIV.beta
!  t                 r

  else if ((shift=="Gammadriver2").or.(shift=="Gammadrivershock2")) then

!    Source for shift.

     sbeta(l,:) = dtbeta(l,:) 

     if (driveradv) then
        sbeta(l,:) = sbeta(l,:) + beta(l,:)*DA_beta(l,:)
     end if

!    Source for dtbeta.

     if (shift=="Gammadriver2") then

!       Standard run.

        if (.not.cosmic_run) then

           sdtbeta(l,:) = drivercsi*sDeltar(l,:)

           if (driveradv) then
              sdtbeta(l,:) = sdtbeta(l,:) + beta(l,:)*(DA_dtbeta(l,:)-DA_Deltar(l,:))
           end if

!       Cosmic run.

        else

           sdtbeta(l,:) = drivercsi*sDeltar(l,:)/cosmobg_a(l)**2

        end if

     else if (shift=="Gammadrivershock2") then

!       Source for fdriver.

        sfdriver(l,:) = - two*third*fdriver(l,:)*DIV_beta(l,:)

        if (driveradv) then
           sfdriver(l,:) = sfdriver(l,:) + beta(l,:)*DA_fdriver(l,:)
        end if

!       Source for dtbeta.

        sdtbeta(l,:) = drivercsi*(fdriver(l,:)*sDeltar(l,:) + sfdriver(l,:)*Deltar(l,:))

        if (driveradv) then
           sdtbeta(l,:) = sdtbeta(l,:) + beta(l,:)*(DA_dtbeta(l,:)-DA_Deltar(l,:))
        end if

        if (driverD0) then
           sdtbeta(l,:) = sdtbeta(l,:) - drivercsi*sfdriver(l,:)*Deltar0(l,:)
        end if

     end if

!    Damping term.

     sdtbeta(l,:) = sdtbeta(l,:) - drivereta*dtbeta(l,:)


! ****************************************************
! ***   NEW SECOND ORDER HYPERBOLIC GAMMA DRIVER   ***
! ****************************************************

! EXPERIMENTAL GAMMADRIVER. NOT IMPLEMENTED.

  else if (shift=="Gammadriver3") then


! ************************
! ***   SANITY CHECK   ***
! ************************

  else

     if (rank==0) then
        print *
        print *, "sources_shift.f90: Unknown shift condition."
        print *
     end if

     call die

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_shift

