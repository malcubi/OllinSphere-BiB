!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/Kodama_mass.f90,v 1.4 2025/09/24 17:25:32 malcubi Exp $

  subroutine Kodamamass

! *************************************************
! ***   CALCULATION OF INTEGRATED KODAMA MASS   ***
! *************************************************
!
! Note: For unstable matter configurations is not clear
! how to calculate the mass of the distribution.
!
! The Kodama mass is a quasi-local conserved energy in a 
! spherically symmetric space time. It is defined with a
! vector K^A:
!
! K^A  =  eps^{AB} d_B R
!
! where R is the areal radius of the 2-sphere at a given (t,r),
! and eps^AB is the Levi-Civita tensor on the 2-sphere.
! Using K^A, we can define the vector S^mu
!
! S^mu  =  T^{mu nu} K_nu
!
! with T^{mu nu} the stress-energy tensor. S^mu is a conserved
! 4-flux and satisfies the conservation law:
!
! d_mu (sqrt(-g) S^mu)  =  0
!
! We can define the conserved mass (Kodama mass) in a sphere 
! of constant radius r at a given time t as follows:
!
!                  r
!                  /
! M(t,r)  =  4 pi  | r^2 alpha exp(6 phi) A^{1/2} B S^t dr
!                  /
!                  0
!
! where S^t is defined as:
!
! S^t  =  rho/alpha (B/A)^1/2 r (2 Dr_phi + Dr_B/2B + 1/r)
!      -  JA/alpha  (B/A)^1/2 r (K/3 + KTB)

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use integrals

! Extra variables.

  implicit none

  logical firstcall  ! Is this the first call?

  integer i,l,i0

  real(8) r0,interp
  real(8) delta,m0,m1
  real(8) R95,R99,MTOT,CMAX
  real(8) aux
  real(8) half,one,two,third,smallpi

  real(8), dimension (0:Nl-1,1-ghost:Nrmax) :: St,Sr ! Auxiliary vectors


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0

  one = 1.d0
  two = 2.d0
  third = 1.d0/3.d0

  smallpi = acos(-one)


! **********************
! ***   INITIALIZE   ***
! **********************

  Kodama_mass = 0.d0


! **************************
! ***   INTEGRATE MASS   ***
! **************************

! Find integrand.

  St = sqrt(B/A)*r/alpha*(rho*(two*D1_phi + half*D1_B/B + one/r) - JA*(third*trK + KTB))

  auxarray = 4.d0*smallpi*r**2*alpha*exp(6.d0*phi)*sqrt(A)*B*St

! Integrate.

  intvar => auxarray

  do l=Nl-1,0,-1
     Kodama_mass(l,:) = integral(l)
  end do

! Restrict integral.

  if (Nl>1) then
     intvar => Kodama_mass
     call restrictintegral
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine Kodamamass








  subroutine KEnergy_flux(l)

! *********************************************
! ***   CALCULATION OF KODAMA ENERGY FLUX   ***
! *********************************************

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer i,l

  real(8) zero,sixth,third
  real(8) half,one,two,smallpi
  real(8) aux,r0,interp


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  half  = 0.5d0
  third = 1.d0/3.d0

  one   = 1.d0
  two   = 2.d0

  smallpi = acos(-one)


! ***********************
! ***   ENERGY FLUX   ***
! ***********************

  auxarray(l,:) = sqrt(B(l,:)/A(l,:))*r(l,:)*(JA(l,:)/(psi4(l,:)*A(l,:))*(two*D1_phi(l,:) &
                + half*D1_B(l,:)/B(l,:)+one/r(l,:)) - SAA(l,:)*(third*trK(l,:) + KTB(l,:)))

  if (shift/="none") then
     auxarray(l,:) = auxarray(l,:) - rho(l,:)*beta(l,:)/alpha(l,:)*sqrt(B(l,:)/A(l,:))*r(l,:) &
                   *(two*D1_phi(l,:) + half*D1_B(l,:)/B(l,:) + one/r(l,:)) &
                   + rho(l,:)*beta(l,:)**2/alpha(l,:)**2*psi4(l,:)*sqrt(A(l,:)*B(l,:))*r(l,:)*(third*trK(l,:) + KTB(l,:)) &
                   - JA(l,:)*beta(l,:)/alpha(l,:)*sqrt(B(l,:)/A(l,:))*r(l,:)*(third*trK(l,:) + KTB(l,:))
  end if

  auxarray(l,:) = 4.d0*smallpi*alpha(l,:)*psi(l,:)**6*sqrt(A(l,:))*B(l,:)*r(l,:)**2*auxarray(l,:)

  sP_Kodama(l,:) = auxarray(l,:)


! ***************
! ***   END   ***
! ***************

  end subroutine KEnergy_flux
