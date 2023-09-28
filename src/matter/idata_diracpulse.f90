!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_diracpulse.f90,v 1.2 2023/08/22 17:59:32 malcubi Exp $

  subroutine idata_diracpulse

! ************************************
! ***   DIRAC PULSE INITIAL DATA   ***
! ************************************

! This initial data leaves the lapse, shift, conformal
! spatial metric and extrinsic curvature as in Minkowski
! (which have already been set up in "initial").
!
! We then give a simple initial profile in the Dirac field
! and solve the Hamiltonian constraint for the conformal
! factor.
!
! Notice that the initial pulse is assumed to be time-symmetric,
! as otherwise we would have a non-zero energy flux and we would
! need to solve the coupled momentum and Hamiltonian constraints.
!
! The Hamiltonian constraint that we need to solve has the form:
!
!  __2                5
!  \/ psi  +  2 pi psi rho  =  0
!
!
! with rho the energy density of the Dirac field:
!
!
! rho  =  1/(2 pi) [ m ( FR**2 + FI**2 - GR**2 - GI**2 )
!
!      +  2/(r psi**2 sqrt(B)) ( FR*GI - FI*GR )
!
!      +  1/(psi**2 sqrt(A)) ( FR*dGI/dr - FI*dGR/dr + GR*dFI/dr - DI*dFR/dr ) ]
!
!
! (assuming that the initial pulse is time-symmetric), and
! where the Laplacian of psi is given in general by:
!
! __2          2
! \/ psi  = [ d psi  +  (2/r - d A/2A + d B/B) d psi ] / A
!              r                r        r      r
!
! In most cases we will have A=B=1, so the Laplacian simplifies.
! But this won't be the case for a transformed radial coordinate.
! Notice also that due to the form of the energy density, the
! Hamiltonian constraint is always non-linear in psi.
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

  integer i,l,p,iter
  integer lmin,lmax
  integer i0,imax,Naux
  integer status(MPI_STATUS_SIZE)
  integer :: maxiter = 100

  real(8) :: one,two
  real(8) :: res,aux
  real(8) :: epsilon = 1.d-10

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: u,u_old
  real(8), dimension (0:Nl-1,1-ghost:Nrmax)   :: CL0,CL1,CL2,CLA,fL1,fL2
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: C0,C1,C2,CA,f1,f2

! u             Global solution of ODE.
! CL0,C0        Local and global coefficient of first derivative in ODE.
! CL1,C1        Local and global coefficient of linear term in ODE.
! CL2,C2        Local and global source term in ODE.
! f12,f2,CA     Auxiliary global arrays.
! fL1,fL2,CLA   Auxiliary local arrays.


! *******************
! ***   NUMBERS   ***
! *******************

  one = 1.d0
  two = 2.d0


! ***********************************************
! ***   SET UP INITIAL PULSE IN DIRAC FIELD   ***
! ***********************************************

! Initial pulse.  Remember that F must be even and G must be odd.

  if (diracGR_r0*diracGI_r0==0.d0) then
     print *, 'For "diracpulse" initial data diracGR_r0 and diracGI_r0 must both be non-zero ...'
     print *, 'Aborting! (subroutine idata_diracpulse)'
     print *
     call die
  end if

! Gaussian profile.

  if (diracprofile=="gaussian") then

!    FR and FI.

     if (diracFR_r0==0.d0) then
        dirac_FR = gaussian(diracFR_a0,0.d0,diracFR_s0)
     else
        dirac_FR = gaussian(diracFR_a0,+diracFR_r0,diracFR_s0) &
                 + gaussian(diracFR_a0,-diracFR_r0,diracFR_s0)
     end if

     if (diracFI_r0==0.d0) then
        dirac_FI = gaussian(diracFI_a0,0.d0,diracFI_s0)
     else
        dirac_FI = gaussian(diracFI_a0,+diracFI_r0,diracFI_s0) &
                 + gaussian(diracFI_a0,-diracFI_r0,diracFI_s0)
     end if

!    GR and GI.

     dirac_GR = gaussian(diracGR_a0,+diracGR_r0,diracGR_s0) &
              - gaussian(diracGR_a0,-diracGR_r0,diracGR_s0)

     dirac_GI = gaussian(diracGI_a0,+diracGI_r0,diracGI_s0) &
              - gaussian(diracGI_a0,-diracGI_r0,diracGI_s0)

! Smooth top-hat profile.

  else if (diracprofile=="tophat") then

!   FR and FI.

     if (diracFR_r0==0.d0) then
        dirac_FR = tophat(diracFR_a0,0.d0,diracFR_s0,diracFR_t0)
     else
        dirac_FR = tophat(diracFR_a0,+diracFR_r0,diracFR_s0,diracFR_t0) &
                 + tophat(diracFR_a0,-diracFR_r0,diracFR_s0,diracFR_t0)
     end if

     if (diracFI_r0==0.d0) then
        dirac_FI = tophat(diracFI_a0,0.d0,diracFI_s0,diracFI_t0)
     else
        dirac_FI = tophat(diracFI_a0,+diracFI_r0,diracFI_s0,diracFI_t0) &
                 + tophat(diracFI_a0,-diracFI_r0,diracFI_s0,diracFI_t0)
     end if

!   GR and GI.

    dirac_GR = tophat(diracGR_a0,+diracGR_r0,diracGR_s0,diracGR_t0) &
             - tophat(diracGR_a0,-diracGR_r0,diracGR_s0,diracGR_t0)

    dirac_GI = tophat(diracGI_a0,+diracGI_r0,diracGI_s0,diracGI_t0) &
             - tophat(diracGI_a0,-diracGI_r0,diracGI_s0,diracGI_t0)

  end if

! Spatial derivatives with finite differences.

  diffvar => dirac_FR
  do l=0,Nl-1
     D1_dirac_FR(l,:) = diff1(l,+1)
  end do

  diffvar => dirac_FI
  do l=0,Nl-1
     D1_dirac_FI(l,:) = diff1(l,+1)
  end do

  diffvar => dirac_GR
  do l=0,Nl-1
     D1_dirac_GR(l,:) = diff1(l,-1)
  end do

  diffvar => dirac_GI
  do l=0,Nl-1
     D1_dirac_GI(l,:) = diff1(l,-1)
  end do

! If we are running on Minkowski spacetime just return.

  if (spacetime=="minkowski") return


! ************************
! ***   SANITY CHECK   ***
! ************************

! Here we check if the initial data does guarantee
! that the momentum density JA is zero.

  if (diractype==1) then

!    For type 1 initial data both F and G must be either
!    purely real of purely imaginary.

     if (diracFR_a0*diracFI_a0/=0.d0) then
        print *, 'For diractype=1 the function F must be either purely real of purely imaginary ...'
        print *, 'Aborting! (subroutine idata_diracpulse)'
        print *
        call die
     end if

     if (diracGR_a0*diracGI_a0/=0.d0) then
        print *, 'For diractype=1 the function G must be either purely real of purely imaginary ...'
        print *, 'Aborting! (subroutine idata_diracpulse)'
        print *
        call die
     end if

  else if (diractype==2) then

     print *, 'diractype=2 initial data not yet working ...'
     print *, 'Aborting! (subroutine idata_diracpulse)'
     print *
     call die

  else if (diractype==3) then

     print *, 'diractype=3 initial data not yet working ...'
     print *, 'Aborting! (subroutine idata_diracpulse)'
     print *
     call die

  else

     print *, 'Unkown type of diracpulse initial data ...'
     print *, 'Aborting! (subroutine idata_diracpulse)'
     print *
     call die

  end if


! ****************************************
! ***   SOLVE HAMILTONIAN CONSTRAINT   ***
! ****************************************

! Here we solve the Hamiltonian constraint for
! the conformal factor psi.

! Message to screen.

  if (rank==0) then
     print *, 'Solving initial data for dirac pulse ...'
     print *
  end if

! Array size.

  Naux = (Nrmax + ghost)

! Initialize array u.

  u = 1.d0

! Calculate energy density. Notice that we don't 
! use the factor 1/2pi since it cancels:
!
! rho  =  m ( FR**2 + FI**2 - GR**2 - GI**2 )
!
!      +  2/(r sqrt(B) psi**2) [ FR GI - FI GR ]
!
!      +  1/(sqrt(A) psi**2) [ FR dGI/dr - FI dGR/dr + GR dFI/dr - GI dFR/dr ]
!
! Remember also that for the initial guess we have psi=1.
!
! We define f1 and f2 since they have different powers of psi,
! and we multiply them with a the A that comes from the laplacian.

  fL1 = dirac_mass*A*(dirac_FR**2 + dirac_FI**2 - dirac_GR**2 - dirac_GI**2)
  fL2 = 2.d0*A/(r*sqrt(B))*(dirac_FR*dirac_GI - dirac_FI*dirac_GR) &
      +(dirac_FR*D1_dirac_GI - dirac_FI*D1_dirac_GR &
      + dirac_GR*D1_dirac_FI - dirac_GI*D1_dirac_FR)/sqrt(A)

  rho = f1 + f2

! Fill in local coefficients of linearized equation.
! Remember that C0 is the coefficient of dpsi/dr,
! C1 is the coefficient of psi, and C2 is the
! right hand side.

  if (newr) then
     CL0 = two/r - 0.5d0*D1_A/A + D1_B/B
     CLA = 0.25d0/r**2*(one - A/B) - 0.25d0*D1_A/(r*A) + 0.75d0*D1_B/(r*B) &  ! This is the Ricci scalar, it should
         - 0.125d0*D1_A*D1_B/(A*B) - 0.0625d0*D1_B**2/B**2 + 0.25d0*D2_B/B    ! be zero for conformally flat data.
  else
     CL0 = two/r
     CLA = 0.d0
  end if

  CL1 = rho + CLA
  CL2 = 0.d0

! Processor zero receives data and solves ODE.

  if (rank==0) then

     if (size==1) then
        imax = Nrl(0)
     else
        imax = Nrl(0) - ghost
     end if

     do i=1-ghost,imax
        C0(:,i) = CL0(:,i)
        C1(:,i) = CL1(:,i)
        CA(:,i) = CLA(:,i)
        f1(:,l) = fL1(:,i)
        f2(:,l) = fL2(:,i)
     end do

!    Iterate over other processors to receive data.

     i0 = imax

     do p=1,size-1

!       Receive local coefficients (CL0,CL1) from other processors
!       and copy them into large arrays (C0,C1).

        do l=0,Nl-1
           call MPI_RECV(CL0(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(CL1(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(CLA(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(fL1(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(fL2(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        end do

        if (p==size-1) then
           imax = Nrl(p)
        else
           imax = Nrl(p) - ghost
        end if

        do i=1,imax
           C0(:,i+i0) = CL0(:,i)
           C1(:,i+i0) = CL1(:,i)
        end do

        i0 = i0 + imax

     end do

!    Solve for initial guess.

     lmin = 0; lmax = Nl-1
     call invertmatrix(lmin,lmax,1.d0,u,C0,C1,C2,+1,"robin")

!    Initialize residual.

     res = 1.d0

!    Begin iterations.

     iter = 0

     do while ((res.gt.epsilon).and.(iter.lt.maxiter))

        iter = iter + 1

!       Save old solution.

        u_old = u

!       Correct linear term and source term.   For this I am using
!       a variation of Newton's method.  We first notice that the
!       equation has the form of a linear term A*u plus a non-linear
!       purely diagonal term F(u):
!
!       A*u  +  F(u)  =  0
!
!       We can then expand in Taylor the non-linear term to first
!       order around a previous guess u_old to find:
!
!       A*u  +  F(u_old)  +  F'(u_old)*(u - u_old)  =  0
!
!       and express the iterations as:
!
!       [A + F'(u_old)]*u_new = F'(u_old)*u_old - F(u_old)
!
!       The effect is then just to add to the diagonal of the
!       original matrix A the term F'(u_old), and introduce a
!       diagonal source term on the right hand side. This has
!       the advantage that I don't need to calculate A*u_old.
!
!       In our case we have:
!
!                    5         3
!       F(u)  =  f1 u   +  f2 u
!
!
!       with:
!
!       f1  =  m ( |F|^2 - |G|^2 )
!
!       f2  =  2/(r*sqrt(B)) [ FR*GI - FI*GR ]
!           +  1/sqrt(A) [ FR*dGI/dr - FI*dGR/dr + GR*dFI/dr - GI*dFR/dr ]
!
!       From this we find:
!
!                       4          2
!       F'(u)  =  5 f1 u  +  3 f2 u
!
!
!       Notice that this then imples:
!
!                                                5              3
!       F'(u_old)*u_old - F(u_old)  =  4 f1 u_old  +  2 f2 u_old

        if (newr) then
           C1 = 5.d0*fL1*u_old**4 + 3.d0*fL2*u_old**2 + CA
           C2 = 4.d0*fL1*u_old**5 + 2.d0*fL2*u_old**3
        else
           C1 = 5.d0*fL1*u_old**4 + 3.d0*fL2*u_old**2
           C2 = 4.d0*fL1*u_old**5 + 2.d0*fL2*u_old**3
        end if

!       Call matrix inversion.

        lmin = 0; lmax = Nl-1
        call invertmatrix(lmin,lmax,1.d0,u,C0,C1,C2,+1,"robin")

!       Find residual (across all grid levels).

        res = 0.d0

        do l=0,Nl-1
           do i=1,Nrtotal
              aux = abs(u(l,i)-u_old(l,i))
              if (aux>res) res = aux
           end do
        end do

        write(*,"(A,I0,A,ES8.2)") ' Iteration: ',iter,'    Residual: ',res

     end do

!    Message in case we reached the maximum iteration number.

     if (iter>=maxiter) then
        print *
        print *, 'Maximum iteration number reached in idata_scalarpulse.f90.'
        print *, 'Initial data solver did not converge, aborting!'
        print *
     else
        print *
        print *, 'Done!'
        print *
     end if

! Other processors send local coefficients (CL0,CL1) to processor 0.

  else

     do l=0,Nl-1
        call MPI_SEND(CL0(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(CL1(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(CLA(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(fL1(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(fL2(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
     end do

  end if

! ************************************************
! ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
! ************************************************

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! array with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  if (size==1) then
     psi = u
  else
     call distribute(0,Nl-1,psi,u)
  end if


! ***************************************
! ***   FIND CONFORMAL FUNCTION phi   ***
! ***************************************

! Find phi.

  phi  = dlog(psi)
  D1_phi = D1_psi/psi

  if (chimethod) then
     chi  = one/psi**chipower
     D1_chi = - dble(chipower)*D1_psi/psi**3
  end if

! Find psi2 and psi4.

  psi2 = psi**2
  psi4 = psi**4


! ***************
! ***   END   ***
! ***************

  end subroutine idata_diracpulse
