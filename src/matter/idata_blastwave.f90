!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_blastwave.f90,v 1.1 2022/08/08 17:59:46 malcubi Exp $

  subroutine idata_blastwave

! ***********************************
! ***   BLAST WAVE INITIAL DATA   ***
! ***********************************

! Blast wave initial data.

! This initial data corresponds to a simple blast wave
! where we have initially a in two regions with higher
! constant density and pressure inside a certain radius,
! and lower constant density and pressure on the outside.
!
! This ninitial data requires a static Minkowski spacetime
! and as such it violates the constrait equations.

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer i,l


! ************************
! ***   SANITY CHECK   ***
! ************************

  if (spacetime=="dynamic") then
     print *, 'Blastwave initial data needs spacetime=background ...'
     print *, 'Aborting! (subroutine idata_blastwave)'
     print *
     call die
  end if


! *******************************
! ***   SET UP INITIAL DATA   ***
! *******************************

! Set up initial density and pressure.

  do l=0,Nl-1
     do i=1-ghost,Nrtotal

!       Inner regions.

        if (r(l,i)<=blast_R) then

           fluid_rho(l,i) = blast_rhol
           fluid_p(l,i)   = blast_pl

!       Outer region.

        else

           fluid_rho(l,i) = blast_rhor
           fluid_p(l,i)   = blast_pr

        end if

     end do
  end do

! Find specific internal energy.

  if (fluid_EOS=="ideal") then
     fluid_e = fluid_p/fluid_rho/(fluid_gamma - 1.d0)
  end if

! Find enthalpy: h = 1 + e + p/rho.

  fluid_h = 1.d0 + fluid_e + fluid_p/fluid_rho

! Set fluid velocity to zero.

  fluid_v = 0.d0

! Set Lorentz factor to one.

  fluid_W = 1.d0

! Conserved quantities.

  fluid_cD = fluid_rho
  fluid_cE = fluid_rho*fluid_e
  fluid_cS = 0.d0

! Initialize artificial viscosity to zero.

  fluid_q = 0.d0


! ***************
! ***   END   ***
! ***************

  end subroutine idata_blastwave


