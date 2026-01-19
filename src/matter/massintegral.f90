
  subroutine massintegral

! ******************************************
! ***   CALCULATION OF INTEGRATED MASS   ***
! ******************************************

! This is the integrated mass.  In order to find it we start
! from the fact that in areal coordinates the spatial metric
! can be written as:
!
!   2        2                         2       2
! dl  =  dr_a / (1 - 2m(r_a)/r_a) + r_a  dOmega
!
!
! with r_a the "areal" (or Schwarzschild) radius, dOmega^2 the
! standard solid angle element, and m(r_a) the so-called "mass function".
!
! In those coordinates the spatial Ricci scalar reduces to the simple
! form:
!
!             2
! R  =  (4/r_a ) dm/dr_a
!
!
! The Hamiltonian constraint then implies:
!
!                2 /                        2        2           2  \
! dm/dr_a  =  r_a  | 4 pi rho  + (1/4) ( KTA  + 2 KTB   - 2/3 trK ) |
!                  \                                                /
!
! The mass function can then be integrated to find:
!
!       /                          2        2           2       2
! m  =  | [ 4 pi rho  + (1/4) ( KTA  + 2 KTB   - 2/3 trK ) ] r_a  dr_a
!       /
!
! Notice that if the spacetime reduces to Schwarzschild far away, then
! the asymptotic value of "m" will in fact correspond to the ADM mass.
!
! Also, for a static spacetime the extrinsic curvature vanishes, and
! the mass function "m" reduces to the Newtonian expression:
!
!       /             2
! m  =  | 4 pi rho r_a  dr_a
!       /
!
! On the other hand, since the code does not work in areal coordinates,
! we have in general the following transformation:
!
! r_a  =  r psi**2 sqrt(B)
!
! with "r" the standard radial coordinate for the code. This implies:
!
!    2          2    6  3/2
! r_a dr_a  =  r  psi  B   ( 1 + r d B / (2 B)  +  2 r d psi / psi ) dr
!                                   r                   r
!
! Notice that this assumes that the spacetime is REGULAR everywhere,
! particularly at the origin, so it will not work for eternal black holes
! such as Schwarzschild or Reissner-Nordstrom.  It will also fail if the
! spacetime develops a singularity during evolution.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use integrals

! Extra variables.

  implicit none

  logical contains

  integer i,l,i0,p
  integer status(MPI_STATUS_SIZE)

  real(8) r0,interp
  real(8) delta,m0,m1
  real(8) R95,R99,MTOT,CMAX
  real(8) aux
  real(8) half,one,two,third,smallpi


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

  mass_int = 0.d0


! **************************
! ***   INTEGRATE MASS   ***
! **************************

! Find integrand.

  auxarray = r**2*psi**6*B**1.5d0*(one + r*(half*D1_B/B + two*D1_phi)) &
           *(4.d0*smallpi*rho + 0.25d0*(KTA**2 + two*KTB**2 - two*third*trK**2))

! Integrate.

  intvar => auxarray

  do l=0,Nl-1
     mass_int(l,:) = integral(l)
  end do

! Restrict integral.

  if (Nl>1) then
     intvar => mass_int
     call restrictintegral
  end if


! *****************************
! ***   OUTPUT TOTAL MASS   ***
! *****************************

! At t=0 output total mass.  Notice that for cases
! when the density decays slowly, such as when we
! have an electric fiels, this number does not mean
! much as it converges very slowly, so we don't output it.

  if (t(0)==0.d0) then

     MTOT = mass_int(0,Nr)

     if (size==1) then
        write(*,'(A,ES23.16)') ' Total integrated mass M_int = ',MTOT
        print *
     else
        if (rank==0) then
           p = size-1
           call MPI_RECV(MTOT,1,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           write(*,'(A,ES23.16)') ' Total integrated mass M_int = ',MTOT
           print *
        else if (rank==size-1) then
           call MPI_SEND(MTOT,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end if
     end if

  end if


! ***********************
! ***   COMPACTNESS   ***
! ***********************

! The local compactness of an object is defined as:
!
! C(r)  :=  M(r)/r_area
!
! Notice that I don't include a factor of 2, so for
! a black hole we expect C ~ 0.5.

  compactness = abs(mass_int)/r_area

! At t=0 output maximum of compactness function.

  CMAX = maxval(compactness)

  if ((t(0)==0.d0).and.(.not.contains(mattertype,"electric"))) then
     if (rank==0) then
        write(*,'(A,ES23.16)') ' Maximum of compactness function (M(r)/r_area) = ',CMAX
        print *
     end if
  end if


! *****************************
! ***   M/R FOR TOV STARS   ***
! *****************************

! For TOV stars we output M/R at t=0.

  if ((t(0)==0.d0).and.(idata=="TOVstar")) then
     if (rank==0) then
        write(*,'(A,ES23.16)') ' Total compactness M/R = ',MTOT/TOV_rad
        print *
     end if
  end if


! *************************
! ***   R99 AND M/R99   ***
! *************************

! Calculate the radius for which we have 99% of the
! total integrated mass. We do this only on the coarse
! grid and use linear interpolation.
!
! Notice that R99 really only makes sense for matter
! distributions that have no definite radius, and
! for which the energy density decays very rapidly.
!
! We only do this calculation at t=0, and for certain types
! of initial data such as boson stars, Proca stars, etc.
!
! In particular, when we have a electric field, like in the
! case of a charged boson star, it really makes no sense
! to find R99 associated to the integtated mass since the
! electric field decays only as 1/r, so that the energy
! density decays very slowly.

  if ((t(0)==0.d0).and. &
     ((idata=="bosonstar").or.(idata=="procastar").or.(idata=="l-procastar").or.(idata=='diracstar') &
        .or.(idata=="chargedboson").or.(idata=="chargedproca").or.(idata=="chargeddirac"))) then

     R99 = 0.d0

!    Single processor run.

     if (size==1) then

        do i=1,Nr
           if ((abs(mass_int(0,i-1))<0.99d0*abs(MTOT)).and.(abs(mass_int(0,i))>=0.99d0*abs(MTOT))) then
              R99 = r(0,i-1) + dr(0)*(0.99d0*MTOT-mass_int(0,i-1))/(mass_int(0,i)-mass_int(0,i-1))
           end if
        end do

!    Multiple processor run. Broadcast the total mass from
!    processor size-1 (the boundary) to all other processors.

     else

        call MPI_BCAST(MTOT,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)

!       Now each processor tries to find R99.

        if (rank==0) then
           i0 = 1
        else
           i0 = 1+ghost
        end if

        do i=i0,Nr
           if ((mass_int(0,i-1)<0.99d0*MTOT).and.(mass_int(0,i)>=0.99d0*MTOT)) then
              R99 = r(0,i-1) + dr(0)*(0.99d0*MTOT-mass_int(0,i-1))/(mass_int(0,i)-mass_int(0,i-1))
           end if
        end do

!       We assume that only one processor found R99 (not necessarily the same one).
!       We then find the maximum values of R99 across all processors.
!
!       Notice that if the integrated mass grows monotonically, then only one processor must
!       have found the a value for R99, for the others that value should remain zero.
!       But, if for some reason the integrated mass is not monotone (numerical error, there
!       is a black hole and the integrated mass makes no sense, we have a non-minimally coupled
!       scalar field, etc.), then we keep the largest value.

        call MPI_Allreduce(R99,aux,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
        R99 = aux

     end if

!    Processor 0 writes result to screen.

     if (rank==0) then
        write(*,'(A,ES23.16)') ' Effective radius R99 from M_int(r) = ',R99
        write(*,'(A,ES23.16)') ' Effective compactness M/R99        = ',MTOT/R99
        print *
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine massintegral
