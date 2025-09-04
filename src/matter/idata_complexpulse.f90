!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_complexpulse.f90,v 1.40 2025/09/04 16:04:58 malcubi Exp $

  subroutine idata_complexpulse

! **************************************
! ***   COMPLEX PULSE INITIAL DATA   ***
! **************************************

! Adapted from scalarpulse initial data.
!
! This initial data leaves the lapse, shift, conformal
! spatial metric and extrinsic curvature as in Minkowski
! (which have already been set up in "initial".)
!
! We then give a simple initial profile in the scalar field
! and solve the Hamiltonian constraint for the conformal
! factor.
!
! The Hamiltonian constraint that we need to solve has the form:
!
!  __2                5
!  \/ psi  +  2 pi psi rho  =  0
!
!
! with rho the energy density of the complex scalar field:
!
!            2       2            4
! rho  =  (xi_r  + xi_i) / 2 A psi  +  V
!
!
! (assuming that the initial pulse is time-symmetric).
! Notice that when the time derivatives are not zero
! the density becomes instead:
!
!              2        2        2       2          4
! rho  =  [ (pi_r  +  pi_i) + (xi_r  + xi_i) / A psi  ] / 2  +  V
!
! But in that case we need to take care in order to satisfy
! the momentum constraints.
!
! The Laplacian of psi is given in general by:
!
! __2           2
! \/ psi  =  [ d psi  +  (2/r - d A/2A + d B/B) d psi ] / A
!               r                r        r      r
!
! In most cases we will have A=B=1, so the Laplacian simplifies.
! But this won't be the case for a transformed radial coordinate.
!
! Notice that the initial pulse is assumed to be time-symmetric,
! (i.e. pi=0), or else purely real with purely imaginary time
! derivative, as otherwise we would have a non-zero energy
! flux and we would need to solve the coupled momentum and
! hamiltonian constraints.
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

! Extra variables.

  implicit none

  logical :: contains

  integer :: i,l,p,iter
  integer :: lmin,lmax
  integer :: i0,imax,Naux
  integer :: status(MPI_STATUS_SIZE)
  integer :: maxiter = 100

  real(8) :: one,two,smallpi
  real(8) :: aux,res
  real(8) :: epsilon = 1.d-10

  real(8), dimension (0:Nl-1,1-ghost:Nrmax)   :: CL0,CL1,CL2,VL
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: C0,C1,C2,CA,VG
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: u,u_old

! u        Global solution of ODE.
! CL0,C0   Local and global coefficient of first derivative in ODE.
! CL1,C1   Local and global coefficient of linear term in ODE.
! CL2,C2   Local and global source term in ODE.


! *******************
! ***   NUMBERS   ***
! *******************

  one = 1.d0
  two = 2.d0

  smallpi = acos(-1.d0)


! *********************************************************
! ***   SET UP GAUSSIAN PULSE IN COMPLEX SCALAR FIELD   ***
! *********************************************************

! Initial pulse.  Remember that the complex scalar field must be even.
!
! Also, if complex_l is non-zero, then the field must behave as r^l
! close to the origin. This means that the gaussian profile must be
! far from the origin, and the top-hat profile can't be used.

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

  else if ((k_parameter/=0.d0).and.(complex_mass>0.d0)) then

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

! If we are running on Minkowski spacetime just return.

  if (spacetime=="minkowski") return


! *************************************************
! ***   INITIALIZE ELECTRIC FIELD WHEN NEEDED   ***
! *************************************************

! This initial data can also be used for the case of a charged
! complex field with zero initial charge density.  As long as
! the real and imaginary parts are different and overlap in some
! region we will have a non-zero initial current density.
! This initial data is described in Section 3.1:  J.M. Torres
! and M. Alcubierre, Gen.Rel.Grav. 46:1773 (2014).

  if (contains(mattertype,"electric")) then
     ePhi = 0.d0
     eAr  = 0.d0
     electric = 0.d0
  end if


! ***************************************************
! ***   ZERO POTENTIAL AND ZERO TIME DERIVATIVE   ***
! ***************************************************

! In the case of zero potential and zero time derivative
! the Hamiltonian  constraint is linear in psi, so we can
! just solve it directly by matrix inversion.
!
! Notice that here we are assuming not only that there
! is no potential, but also that the time derivatives
! of the scalar field are zero, as otherwise the equation
! would no longer be linear.
!
! For l-complex fields we need to add an extra term
! to the left hand side of the above equation:
!
! l(l+1) |phi|^2 / r^2

  if ((complexpotential=="none").and.(k_parameter==0.d0)) then

!    Message to screen.

     if (rank==0) then
        if (complex_l==0) then
           print *, 'Solving initial data for a complex pulse with zero potential'
        else
           write(*,'(A,I2)') ' Solving initial data for a complex pulse with zero potential and complex_l =',complex_l
        end if
        print *
     end if

!    Array size.

     Naux = (Nrmax + ghost)

!    Fill in local coefficients of linearized equation.

     if (newr) then
        CL0 = two/r - 0.5d0*D1_A/A + D1_B/B
        CL1 = smallpi*(complex_xiR**2 + complex_xiI**2) &
            + 0.25d0/r**2*(one - A/B) - 0.25d0*D1_A/(r*A) + 0.75d0*D1_B/(r*B) &  ! This is the Ricci scalar, it should
            - 0.125d0*D1_A*D1_B/(A*B) - 0.0625d0*D1_B**2/B**2 + 0.25d0*D2_B/B    ! be zero for conformally flat data.
     else
        CL0 = two/r
        CL1 = smallpi*(complex_xiR**2 + complex_xiI**2)
     end if

     CL2 = 0.d0

     if (complex_l/=0) then
        CL1 = CL1 + smallpi*complex_l*(complex_l+1)*(complex_phiR**2 + complex_phiI**2)/r**2
     end if

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
!
! For l-complex fields we need to add an extra term
! to the left hand side of the above equation:
!
! l(l+1) |phi|^2 / r^2
!
! Notice also that when k_parameter is not zero we need
! to take into account the contribution of the time
! dervatives. The easiest way to do this is simply
! add the corresponding term always.  In most cases
! it will just be zero.

  else

!    Message to screen.

     if (complex_l==0) then
        if (rank==0) then
           print *, 'Solving initial data for a complex pulse with non-zero potential'
           print *
        end if
     else if (complex_lambda/=0.d0) then
        if (rank==0) then
           print *, 'A complex scalar field with non-zero complex_l and non-zero complex_lambda is inconsistent'
           print *, 'Aborting (subroutine idata_complexpulse.f90)'
        end if
        call die
     else
        if (rank==0) then
           write(*,'(A,I2)') ' Solving initial data for a massive complex pulse with complex_l =',complex_l
        end if
     end if

!    Find complex field potential.

     do l=0,Nl-1
        call potential(l)
     end do

!    Array size.

     Naux = (Nrmax + ghost)

!    Fill in local coefficients of linearized equation.

     if (newr) then
        CL0 = two/r - 0.5d0*D1_A/A + D1_B/B
        CL1 = smallpi*(complex_xiR**2 + complex_xiI**2) &
            + 0.25d0/r**2*(one - A/B) - 0.25d0*D1_A/(r*A) + 0.75d0*D1_B/(r*B) &  ! This is the Ricci scalar, it should
            - 0.125d0*D1_A*D1_B/(A*B) - 0.0625d0*D1_B**2/B**2 + 0.25d0*D2_B/B    ! be zero for conformally flat data.
     else
        CL0 = two/r
        CL1 = smallpi*(complex_xiR**2 + complex_xiI**2)
     end if

     CL2 = 0.d0

     if (complex_l/=0) then
        CL1 = CL1 + smallpi*complex_l*(complex_l+1)*(complex_phiR**2 + complex_phiI**2)/r**2
     end if

!    Fill in local potential (plus time derivative terms).

     VL = complex_V + 0.5d0*(complex_piR**2 + complex_piI**2)

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

!       For initial guess solve taking psi=1 on the right-hand side.

        C2 = - 2.d0*smallpi*VG

        lmin = 0; lmax = Nl-1
        call invertmatrix(lmin,lmax,1.d0,u,C0,C1,C2,+1,"robin")

!       Initialize residual.

        res = 1.d0

!       Begin iterations.

        iter = 0

        do while ((res.gt.epsilon).and.(iter.lt.maxiter))

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
           print *, 'Maximum iteration number reached in idata_complexpulse.f90.'
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

! Derivative of psi.

  diffvar => psi

  do l=0,Nl-1
     D1_psi(l,:) = diff1(l,+1)
  end do


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

  end subroutine idata_complexpulse

