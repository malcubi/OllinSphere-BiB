!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_TOVcomplex.f90,v 1.3 2025/09/24 17:31:02 malcubi Exp $

  subroutine idata_TOVcomplex

! *******************************************
! ***   TOV STAR + COMPLEX SCALAR FIELD   ***
! ***       (FERMION-BOSON SYSTEM)        ***
! *******************************************

! This initial data is for a TOV star plus a complex
! scalar field pulse.  It is a fermion-boson system,
! but NOT a fermion-boson star since the complex
! scalar field is not assumed to be stationary.
!
! The idea is to first solve for the TOV star and later
! add a scalar field pulse and solve the Hamiltonian
! constraint again.
!
! For the initial data I use the areal coordinate
! gauge in order to be consistent with the TOV
! star data.  That is, I set B=psi=1 and solve
! for the radial metric coefficient A.
!
! The Hamiltonian constraint then has the form:
!
! dA/dr = A [ (1 - A) / r + 8 pi r A rho ]
!
! where now the density "rho" must include the
! contributions from the fluid and scalar field.
!
! rho  =  rhoF  +  rhoS
!
! with:
!
! rhoF  =  fluid_cE  +  fluid_cD
!
! We separate the scalar field contribution in two:
! one term that depends on the metric function A,
! and one that does not:
!
! rhoS  =  rhoX/A  +  rhoP
!
! with:
!
! rhoX  =  (xiR**2 + xiI**2) / 2
!
! rhoP  =  (piR**2 + piI**2) / 2  +  V
!
! Finally, we solve the maximal slicing condition again.
!
!
! NOTE FOR PARALLEL RUNS:  The initial data is not really
! solved in parallel.  It is in fact solved only on processor
! zero on a full size array, and then it is distributed
! among the other processors.  This is slow, but works.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays
  use derivatives
  use radialfunctions

! Extra variables.

  implicit none

  integer :: i,l,p
  integer :: i0,imin,imax,iaux,Naux
  integer :: status(MPI_STATUS_SIZE)

  real(8) :: r0,delta,rm
  real(8) :: A0,A_rk
  real(8) :: rho_rk,rhoF_rk,rhoX_rk,rhoP_rk
  real(8) :: s1,s2,s3,s4
  real(8) :: half,one,smallpi

! Global arrays.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr,A_g
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rhoF_g,rhoX_g,rhoP_g

! Local arrays.

  real(8), dimension (0:Nl-1,1-ghost:Nrmax) :: rhoF,rhoX,rhoP

! rr       Radial coordinate global array.
! A_g      Radial metric global array.
! rho*     Local and global contributions to energy density.


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0
  one  = 1.d0

  smallpi = acos(-1.d0)


! **************************************
! ***   FIRST SOLVE FOR A TOV STAR   ***
! **************************************

  call idata_TOVstar


! ************************************************
! ***   NOW ADD A COMPLEX SCALAR FIELD PULSE   ***
! ************************************************

! Message to screen.

  if (rank==0) then
     print *, 'Adding a complex scalar field pulse ...'
     print *
  end if

! Remember that the complex scalar field must be even.

! Gaussian profle.

  if (complexprofile=="gaussian") then

     if (complexR_r0==0.d0) then
        complex_phiR = gaussian(complexR_a0,0.d0,complexR_s0)
     else
        complex_phiR = gaussian(complexR_a0,+complexR_r0,complexR_s0) &
                     + gaussian(complexR_a0,-complexR_r0,complexR_s0)
     end if

     if (complexI_r0==0.d0) then
        complex_phiI = gaussian(complexI_a0,0.d0,complexI_s0)
     else
        complex_phiI = gaussian(complexI_a0,+complexI_r0,complexI_s0) &
                     + gaussian(complexI_a0,-complexI_r0,complexI_s0)
     end if

! Smooth top-hat profile.

  else if (complexprofile=="tophat") then

     if (complexR_r0==0.d0) then
        complex_phiR = tophat(complexR_a0,0.d0,complexR_s0,complexR_t0)
     else
        complex_phiR = tophat(complexR_a0,+complexR_r0,complexR_s0,complexR_t0) &
                     + tophat(complexR_a0,-complexR_r0,complexR_s0,complexR_t0)
     end if

     if (complexI_r0==0.d0) then
        complex_phiI = tophat(complexI_a0,0.d0,complexI_s0,complexI_t0)
     else
        complex_phiI = tophat(complexI_a0,+complexI_r0,complexI_s0,complexI_t0) &
                     + tophat(complexI_a0,-complexI_r0,complexI_s0,complexI_t0)
     end if

  end if

! Time derivatives.

  if (k_parameter==0.d0) then

!    The initial time derivatives are set to zero.

     complex_piR = 0.d0
     complex_piI = 0.d0

  else if (k_parameter/=0.d0) then

!    When k_parameter is not zero we set up an initially
!    purely real field.

     complex_phiI = 0.d0

!    The initial time derivative is harmonic with frequency k_parameter.
!    This guarantees that the initial data still satisfies the momentum
!    constraint trivially.

     complex_piR = 0.d0
     complex_piI = k_parameter*complex_phiR

  end if

! The spatial derivatives are calculated with finite differences.

  diffvar => complex_phiR

  do l=0,Nl-1
     complex_xiR(l,:) = diff1(l,+1)
  end do

  diffvar => complex_phiI

  do l=0,Nl-1
     complex_xiI(l,:) = diff1(l,+1)
  end do


! ******************************************
! ***   COMPLEX SCALAR FIELD POTENTIAL   ***
! ******************************************

! Find complex field potential.

  do l=0,Nl-1
     call potential(l)
  end do


! ****************************************
! ***   ENERGY DENSITY CONTRIBUTIONS   ***
! ****************************************

! Fluid contribution:
!
! rho_fluid  =  fluid_cE + fluid_cD
!
! Take care not to include the artificial atmosphere.

  do l=0,Nl-1
     do i=1-ghost,Nr
        if ((fluid_cD(l,i)>fluid_atmos).and.(fluid_rho(l,i)>fluid_atmos)) then
           rhoF(l,i) = fluid_cE(l,i) + fluid_cD(l,i)
        else
           rhoF(l,i) = 0.d0
        end if
     end do
  end do

! Complex scalar field contribution. We separate
! this in two contributions depending on the power
! of A, so that we have the full scalar field
! contribution given by:
!
! rho_scalar = rhoX/A + rhoP

  rhoX = 0.5d0*(complex_xiR**2 + complex_xiI**2)
  rhoP = 0.5d0*(complex_piR**2 + complex_piI**2) + complex_V


! *****************************************
! ***   FILL IN GLOBAL DENSITY ARRAYS   ***
! *****************************************

! Array size.

  Naux = (Nrmax + ghost)

! Procesor 0.

  if (rank==0) then

     if (size==1) then
        imax = Nrl(0)
     else
        imax = Nrl(0) - ghost
     end if

     do i=1-ghost,imax
        rhoF_g(:,i) = rhoF(:,i)
        rhoX_g(:,i) = rhoX(:,i)
        rhoP_g(:,i) = rhoP(:,i)
     end do

!    Iterate over other processors to receive data.

     i0 = imax

     do p=1,size-1

!       Receive local density arrays (rhoF,rhoX,rhoP) from other
!       processors and copy them into global arrays.

           do l=0,Nl-1
              call MPI_RECV(rhoF(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(rhoX(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(rhoP(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           end do

           if (p==size-1) then
              imax = Nrl(p)
           else
              imax = Nrl(p) - ghost
           end if

           do i=1,imax
              rhoF_g(:,i+i0) = rhoF(:,i)
              rhoX_g(:,i+i0) = rhoX(:,i)
              rhoP_g(:,i+i0) = rhoP(:,i)
           end do

        i0 = i0 + imax

     end do

! Other processors send local density arrays to processor 0.

  else

     do l=0,Nl-1
        call MPI_SEND(rhoF(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(rhoX(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(rhoP(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
     end do

  end if


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

  do l=0,Nl-1
     do i=1-ghost,Nrtotal
        rr(l,i) = (dble(Nmin(rank) + i) - half)*dr(l)
     end do
  end do


! ***********************
! ***   INTEGRATION   ***
! ***********************
  
! We solve the system of equations using fourth
! order Runge--Kutta.

! For the moment only processor 0 solves for the
! initial data.  The solution is later distributed
! among processors.

! Initialize array for metric coefficient A.

  A_g = 1.d0

  
! *******************************************
! ***   ONLY PROCESSOR 0 SOLVES THE ODE   ***
! *******************************************

  if (rank==0) then


!    ****************************************
!    ***   SOLVE HAMILTONIAN CONSTRAINT   ***
!    ****************************************

!    The Hamiltonian constraint is solved for A
!    using a fourth order Runge-Kutta method.

!    Loop over grid levels. We solve from fine to coarse grid.

     do l=Nl-1,0,-1

!       Find initial point. Only the finest grid
!       integrates from the origin.

        if (l==Nl-1) then
           imin = 1
        else
           imin = Nrtotal/2
        end if

!       For coarse grids we interpolate the initial point.

        if (l<Nl-1) then
           A_g(l,imin-1) = (9.d0*(A_g(l+1,Nrtotal-2)+A_g(l+1,Nrtotal-3)) &
                         - (A_g(l+1,Nrtotal-4)+A_g(l+1,Nrtotal-1)))/16.d0
        end if

!       Fourth order Runge-Kutta.

        do i=imin,Nrtotal

!          Grid spacing and values at first point
!          if we start from the origin (finer grid).

           if (i==1) then

              delta = half*dr(l)
              r0 = 0.d0

              A0 = 1.d0

              rhoF_rk = (9.d0*(rhoF_g(l,0)+rhoF_g(l,1)) - (rhoF_g(l,-1)+rhoF_g(l,2)))/16.d0
              rhoX_rk = (9.d0*(rhoX_g(l,0)+rhoX_g(l,1)) - (rhoX_g(l,-1)+rhoX_g(l,2)))/16.d0
              rhoP_rk = (9.d0*(rhoP_g(l,0)+rhoP_g(l,1)) - (rhoP_g(l,-1)+rhoP_g(l,2)))/16.d0

!          Grid spacing and values at previous grid point.

           else

              delta = dr(l)
              r0 = rr(l,i-1)

              A0 = A_g(l,i-1)

              rhoF_rk = rhoF_g(l,i-1)
              rhoX_rk = rhoX_g(l,i-1)
              rhoP_rk = rhoP_g(l,i-1)

           end if

!          I) First Runge-Kutta step.

!          Sources at first grid point if we start
!          from the origin (for finer grid).

           if (i==1) then

!             At the origin we have A'=0.

              s1 = 0.d0

!          Sources at previous grid point.

           else

              rm = r0
              A_rk = A0

              rho_rk = rhoF_rk + (rhoX_rk/A_rk + rhoP_rk)
              s1 = delta*A_rk*((1.d0-A_rk)/rm + 8.d0*smallpi*rm*A_rk*rho_rk)

           end if

!          II) Second Runge-Kutta step.

           rm = r0 + half*delta
           A_rk = A0 + half*s1

           if (i==1) then
              rhoF_rk = 0.5d0*(rhoF_rk + rhoF_g(l,1))
              rhoX_rk = 0.5d0*(rhoX_rk + rhoX_g(l,1))
              rhoP_rk = 0.5d0*(rhoP_rk + rhoP_g(l,1))
           else
              rhoF_rk = (9.d0*(rhoF_g(l,i)+rhoF_g(l,i-1)) - (rhoF_g(l,i-2)+rhoF_g(l,i+1)))/16.d0
              rhoX_rk = (9.d0*(rhoX_g(l,i)+rhoX_g(l,i-1)) - (rhoX_g(l,i-2)+rhoX_g(l,i+1)))/16.d0
              rhoP_rk = (9.d0*(rhoP_g(l,i)+rhoP_g(l,i-1)) - (rhoP_g(l,i-2)+rhoP_g(l,i+1)))/16.d0
           end if

           rho_rk = rhoF_rk + (rhoX_rk/A_rk + rhoP_rk)
           s2 = delta*A_rk*((1.d0-A_rk)/rm + 8.d0*smallpi*rm*A_rk*rho_rk)

!          III) Third Runge-Kutta step.

           A_rk = A0 + half*s2

           s3 = delta*A_rk*((1.d0-A_rk)/rm + 8.d0*smallpi*rm*A_rk*rho_rk)

!          IV) Fourth Runge-Kutta step.

           rm = r0 + delta
           A_rk = A0 + s3

           rhoF_rk = rhoF_g(l,i)
           rhoX_rk = rhoX_g(l,i)
           rhoP_rk = rhoP_g(l,i)

           rho_rk = rhoF_rk + (rhoX_rk/A_rk + rhoP_rk)
           s4 = delta*A_rk*((1.d0-A_rk)/rm + 8.d0*smallpi*rm*A_rk*rho_rk)

!          Advance to next grid point.

           A_g(l,i) = A0 + (s1 + 2.d0*(s2 + s3) + s4)/6.d0

        end do

!       Ghost zones.

        do i=1,ghost
           A_g(l,1-i) = A_g(l,i)
        end do

     end do


!    ************************************
!    ***   RESTRICT TO COARSE GRIDS   ***
!    ************************************

!    Restrict solution from fine to coarse grid.
!    We don't call the subroutine "restrict"
!    since here we are running only on processor 0.
!    We use cubic interpolation.

     do l=Nl-1,1,-1

        do i=1,Nrtotal-ghost,2
           iaux = i/2 + 1
           rm = rr(l-1,iaux)
           A_g(l-1,iaux) = (9.d0*(A_g(l,i)+A_g(l,i+1)) - (A_g(l,i-1)+A_g(l,i+2)))/16.d0
        end do

!       Fix ghost zones.

        do i=1,ghost
           A_g(l-1,1-i) = A_g(l-1,i)
        end do

     end do


! *************************************
! ***   FINISHED FINDING SOLUTION   ***
! *************************************

! When we get here the solution has been found,
! so we close the "if" statement for processor 0.

  end if


! ************************************************
! ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
! ************************************************

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! array with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  if (size==1) then
     A = A_g
  else
     call distribute(0,Nl-1,A,A_g)
  end if


! ***************************
! ***   FIND LAPSE AGAIN  ***
! ***************************

! Now solve for the lapse using maximal slicing.
! But we first need to calculate all auxiliary
! variables on all grid levels, which include
! the stress-energy tensor.

  do l=0,Nl-1
     call auxiliary(l)
  end do

  call alphamaximal(0,Nl-1,"robin",1.d0)


! ***************
! ***   END   ***
! ***************

  end subroutine idata_TOVcomplex

