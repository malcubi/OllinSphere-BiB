!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_scalarpulse.f90,v 1.70 2024/03/11 16:52:16 malcubi Exp $

  subroutine idata_scalarpulse

! *************************************
! ***   SCALAR PULSE INITIAL DATA   ***
! *************************************

! This initial data leaves the lapse, shift, conformal
! spatial metric and extrinsic curvature as in Minkowski
! (which have already been set up in "initial").
!
! We then give a simple initial profile in the scalar field
! and solve the Hamiltonian constraint for the conformal
! factor.
!
! Notice that the initial pulse is assumed to be time-symmetric,
! (i.e. pi=0), as otherwise we would have a non-zero energy
! flux and we would need to solve the coupled momentum and
! Hamiltonian constraints.
!
! The Hamiltonian constraint that we need to solve has the form:
!
!  __2                5
!  \/ psi  +  2 pi psi rho  =  0
!
!
! with rho the energy density of the scalar field:
!
!           2          4
! rho  =  xi / (2 A psi )  +  V
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
!
! Notice also that when the potential vanishes, the Hamiltonian
! constraint becomes linear in psi.
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

! Include Numerical Recipes modules.

  use nrtype
  use nr, only : solvde

! Extra variables.

  implicit none

  integer i,l,p,iter
  integer lmin,lmax
  integer i0,imax,Naux
  integer status(MPI_STATUS_SIZE)
  integer :: maxiter = 1000

  real(8) half,one,two,smallpi
  real(8) aux,res
  real(8) :: epsilon = 1.d-10

  real(8), dimension(0:Nl-1,1-ghost:Nrmax)   :: CL0,CL1,CL2,VL
  real(8), dimension(0:Nl-1,1-ghost:Nrtotal) :: C0,C1,C2,CA,VG
  real(8), dimension(0:Nl-1,1-ghost:Nrtotal) :: u,u_old

! u        Global solution of ODE.
! CL0,C0   Local and global coefficient of first derivative in ODE.
! CL1,C1   Local and global coefficient of linear term in ODE.
! CL2,C2   Local and global source term in ODE.
! CA       Auxiliary global array.
! VL,VG    Local and global potential.


! Now declare arrays for using the Numerical Recipes routine "solvde".
! The declarations I4B and DP are nrtypes.

  integer(I4B) nb
  integer(I4B), dimension(2) :: indexv

  real(DP) slowc
  real(DP), dimension(2) :: scalv
  real(DP), dimension(1:Nrtotal) :: x    ! Radial position for solvde.
  real(DP), dimension(5,1:Nrtotal) :: y  ! Solution vector for solvde.


! *******************
! ***   NUMBERS   ***
! *******************

  half = 0.5d0
  one  = 1.d0
  two  = 2.d0

  smallpi = acos(-1.d0)


! ************************************************
! ***   SET UP INITIAL PULSE IN SCALAR FIELD   ***
! ************************************************

! Initial pulse.  Remember that the scalar field must be even.

! Gaussian profile.

  if (scalarprofile=="gaussian") then

     if (scalar_r0==0.d0) then
        scalar_phi = gaussian(scalar_a0,0.d0,scalar_s0)
     else
        scalar_phi = gaussian(scalar_a0,+scalar_r0,scalar_s0) &
                   + gaussian(scalar_a0,-scalar_r0,scalar_s0)
     end if

! r2 Gaussian profile.

  else if (scalarprofile=="r2gaussian") then

     if (scalar_r0==0.d0) then
        scalar_phi = r2gaussian(scalar_a0,0.d0,scalar_s0)
     else
        scalar_phi = r2gaussian(scalar_a0,+scalar_r0,scalar_s0) &
                   + r2gaussian(scalar_a0,-scalar_r0,scalar_s0)
     end if

! Smooth top-hat profile.

  else if (scalarprofile=="tophat") then

     if (scalar_r0==0.d0) then
        scalar_phi = tophat(scalar_a0,0.d0,scalar_s0,scalar_t0)
     else
        scalar_phi = tophat(scalar_a0,+scalar_r0,scalar_s0,scalar_t0) &
                   + tophat(scalar_a0,-scalar_r0,scalar_s0,scalar_t0)
     end if

  end if

! The spatial derivative is calculated with finite differences.

  diffvar => scalar_phi

  do l=0,Nl-1
     scalar_xi(l,:) = diff1(l,+1)
  end do

! The initial time derivative is set to 0.

  scalar_pi = 0.d0

! If we are running on Minkowski spacetime just return.

  if (spacetime=="minkowski") return


! **************************
! ***   ZERO POTENTIAL   ***
! **************************

! In the case of zero potential the Hamiltonian
! constraint is linear in psi, so we can just solve
! it directly by matrix inversion.

  if ((scalarpotential=="none").and.(.not.scalar_relax)) then

!    Message to screen.

     if (rank==0) then
        print *, 'Solving initial data for a scalar pulse with zero potential ...'
        print *
     end if

!    Array size.

     Naux = (Nrmax + ghost)

!    Fill in (local) coefficients of linear equation.
!    Remember that C0 is the coefficient of dpsi/dr,
!    C1 is the coefficient of psi, and C2 is the
!    right hand side.

     if (newr) then
        CL0 = two/r - 0.5d0*D1_A/A + D1_B/B
        CL1 = smallpi*scalar_xi**2 &
            + 0.25d0/r**2*(one - A/B) - 0.25d0*D1_A/(r*A) + 0.75d0*D1_B/(r*B) &  ! This is the Ricci scalar, it should
            - 0.125d0*D1_A*D1_B/(A*B) - 0.0625d0*D1_B**2/B**2 + 0.25d0*D2_B/B    ! be zero for conformally flat data.
     else
        CL0 = two/r
        CL1 = smallpi*scalar_xi**2
     end if

     CL2 = 0.d0

!    Processor zero receives data and solves ODE.

     if (rank==0) then

        if (size==1) then
           imax = Nrl(0)
        else
           imax = Nrl(0) - ghost
        end if

        do i=1-ghost,imax
           C0(:,i) = CL0(:,i)
           C1(:,i) = CL1(:,i)
        end do

!       Iterate over other processors to receive data.

        i0 = imax

        do p=1,size-1

!          Receive local coefficients (CL0,CL1) from other processors
!          and copy them into large arrays (C0,C1).

           do l=0,Nl-1
              call MPI_RECV(CL0(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(CL1(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
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

!       Call ODE solver. Notice that in this case C2 is always 0.

        C2 = 0.d0

        lmin = 0; lmax = Nl-1

        call invertmatrix(lmin,lmax,1.d0,u,C0,C1,C2,+1,"robin")

!    Other processors send local coefficients (CL0,CL1) to processor 0.

     else

        do l=0,Nl-1
           call MPI_SEND(CL0(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(CL1(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end do

     end if

!    Message to screen.

     if (rank==0) then
        print *, 'Done!'
        print *
     end if


! ******************************
! ***   NON-ZERO POTENTIAL   ***
! ******************************

! When the potential is not zero the equation to solve
! is no longer linear. We solve this equation using
! Newton's iterative method.

  else if (.not.scalar_relax) then

!    Message to screen.

     if (rank==0) then
        print *, 'Solving initial data for scalar pulse with non-zero potential ...'
        print *
     end if

!    Find scalar field potential.

     do l=0,Nl-1
        call potential(l)
     end do

!    Array size.

     Naux = (Nrmax + ghost)

!    Fill in local coefficients of linearized equation.

     if (newr) then
        CL0 = two/r - 0.5d0*D1_A/A + D1_B/B
        CL1 = smallpi*scalar_xi**2 &
            + 0.25d0/r**2*(one - A/B) - 0.25d0*D1_A/(r*A) + 0.75d0*D1_B/(r*B) &  ! This is the Ricci scalar, it should
            - 0.125d0*D1_A*D1_B/(A*B) - 0.0625d0*D1_B**2/B**2 + 0.25d0*D2_B/B    ! be zero for conformally flat data.
     else
        CL0 = two/r
        CL1 = smallpi*scalar_xi**2
     end if

!    Fill in local potential.

     VL = scalar_V

!    Processor zero receives data and solves ODE.

     if (rank==0) then

        if (size==1) then
           imax = Nrl(0)
        else
           imax = Nrl(0) - ghost
        end if

        do i=1-ghost,imax
           C0(:,i) = CL0(:,i)
           C1(:,i) = CL1(:,i)
           VG(:,i) = VL (:,i)
        end do

!       Iterate over other processors to receive data.

        i0 = imax

        do p=1,size-1

!          Receive local coefficients (CL0,CL1) from other processors
!          and copy them into large arrays (C0,C1).

           do l=0,Nl-1
              call MPI_RECV(CL0(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(CL1(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(VL(l,:) ,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           end do

           if (p==size-1) then
              imax = Nrl(p)
           else
              imax = Nrl(p) - ghost
           end if

           do i=1,imax
              C0(:,i+i0) = CL0(:,i)
              C1(:,i+i0) = CL1(:,i)
              VG(:,i+i0) = VL (:,i)
           end do

           i0 = i0 + imax

        end do

!       For initial guess solve with zero potential.

        C2 = 0.d0

        lmin = 0; lmax = Nl-1
        call invertmatrix(lmin,lmax,1.d0,u,C0,C1,C2,+1,"robin")

!       Initialize residual.

        res = 1.d0

!       Begin iterations.

        iter = 0

        do while ((res>epsilon).and.(iter<maxiter))

           iter = iter + 1

!          Save old solution.

           u_old = u

!          Correct linear term and source term.  For this I am using
!          a variation of Newton's method.  We first notice that the
!          equation has the form of a linear term A*u plus a non-linear
!          purely diagonal term F(u):
!
!          A*u  +  F(u)  =  0
!
!          We can then expand in Taylor the non-linear term to first
!          order around a previous guess u_old to find:
!
!          A*u  +  F(u_old)  +  F'(u_old)*(u - u_old)  =  0
!
!          and express the iterations as:
!
!          [A + F'(u_old)]*u_new = F'(u_old)*u_old - F(u_old)
!
!          The effect is then just to add to the diagonal of the
!          original matrix A the term F'(u_old), and introduce a
!          diagonal source term on the right hand side. This has
!          the advantage that I don't need to calculate A*u_old.
!
!          In our case F(u) = 2*pi*V*u**5, so that:  F'(u) = 10*pi*V*u**4.
!          Notice also that:  F'(u_old)*u_old - F(u_old)  =  8*pi*V*u_old**5

           if (newr) then
              CA = C1 + 10.d0*smallpi*A*VG*u_old**4
              C2 = 8.d0*smallpi*A*VG*u_old**5
           else
              CA = C1 + 10.d0*smallpi*VG*u_old**4
              C2 = 8.d0*smallpi*VG*u_old**5
           end if

!          Call matrix inversion.

           lmin = 0; lmax = Nl-1
           call invertmatrix(lmin,lmax,1.d0,u,C0,CA,C2,+1,"robin")

!          Find residual (across all grid levels).

           res = 0.d0

           do l=0,Nl-1
              do i=1,Nrtotal
                 aux = abs(u(l,i)-u_old(l,i))
                 if (aux>res) res = aux
              end do
           end do

           write(*,"(A,I0,A,ES8.2)") ' Iteration: ',iter,'    Residual: ',res

        end do

!       Message in case we reached the maximum iteration number.

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

!    Other processors send local coefficients (CL0,CL1) to processor 0.

     else

        do l=0,Nl-1
           call MPI_SEND(CL0(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(CL1(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(VL(l,:) ,Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end do

     end if


! **********************
! ***   RELAXATION   ***
! **********************

! Here I use the routine "solvde" from Numerical Recipes,
! which solves the equation with a relaxation that uses
! Newton's method.
!
! The routine "solvde" calls the subroutine "difeq",
! where one must define the equations to solve and
! the derivatives for Newton's method.  Have a look
! at it, it is in the directory "base".
!
! This section is just for testing, as it is only second order.
! If you want to use it you must set "scalar_relax = .true."
! in the parameter file.

  else

!    Message to screen.

     if (rank==0) then
        print *, 'Solving initial data for a scalar pulse using the "solvde" routine from Numerical Recipes ...'
        print *
     end if

!    Find scalar field potential.

     do l=0,Nl-1
        call potential(l)
     end do

!    One boundary condition at left hand side (dpsi/dr=0).

     nb = 1

!    Set "solvde" parameters (slowc,scalv) to 1 (the default).

     slowc = 1.d0
     scalv = 1.d0

!    Trivial ordering.

     indexv(1) = 1
     indexv(2) = 2

!    The routine "solvde" is adapated to first order systems,
!    so here I take:
!
!    y1 = dpsi/dr    (initialized to 0)
!    y2 = psi        (initialized to 1)
!
!    Notice that the first index of y(j,i) refers to the unkown
!    variables, and the second one to the grid point.
!
!    The reason for taking y1 as the derivative is the fact
!    that the "solvde" routine expects the boundary conditions
!    at the left hand side to correspond to the first variables
!    in the array y. Here we only have 1 boundary condition at
!    the origin, dpsi/dr=0, so the first array must correspond
!    to dpsi/dr.
!
!    We solve from coarse to fine grid.  On the coarse grid we
!    use the correct boundary condition at the last grid point,
!    but on the fine grids we use Dirichlet boundary conditions
!    interpolating results from the coarser grids. Since we need
!    to know the current grid level I pass it to the Numerical Recipes
!    routines as the last parameter.

     do l=0,Nl-1

!       Initialize y1 and y2.

        y(1,:) = 0.d0
        y(2,:) = 1.d0

!       Call "solvde" routine.

        write(*,"(A,I0)") ' Level ',l
        print *

        call solvde(maxiter,epsilon,slowc,scalv,indexv,nb,x,y,l)

        print *

!       Copy solution for conformal factor onto array u.

        do i=1,Nrtotal
           u(l,i) = y(2,i)
           psi(l,i) = u(l,i)
        end do

!       Ghost zones.

        do i=1,ghost
           u(l,1-i) = u(l,i)
        end do

     end do

!    Restrict solution from fine to coarse grid.
!    We don't call the subroutine "restrict"
!    since here we are running only on processor 0.
!    We use cubic interpolation.

     do l=Nl-1,1,-1

!       Substract constant difference at edge of fine grid.

        i0 = Nrtotal/2 - 1

        aux = u(l-1,i0) - (9.d0*(u(l,2*i0)+u(l,2*i0-1)) - (u(l,2*i0-2)+u(l,2*i0+1)))/16.d0

        u(l-1,:) = u(l-1,:) - aux

!       Interpolate.

        do i=1,Nrtotal-ghost,2
           u(l-1,i/2+1) = (9.d0*(u(l,i)+u(l,i+1)) - (u(l,i-1)+u(l,i+2)))/16.d0
        end do

!       Fix symmetries.

        do i=1,ghost
           u(l-1,1-i) = u(l-1,i)
        end do

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

! Find phi and chi.

  phi = dlog(psi)

  if (chimethod) then
     chi = one/psi**chipower
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine idata_scalarpulse
