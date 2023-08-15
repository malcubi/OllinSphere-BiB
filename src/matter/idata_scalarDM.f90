!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_scalarDM.f90,v 1.23 2021/11/29 19:14:02 malcubi Exp $

  subroutine idata_scalarDM

! ****************************************
! ***   SCALAR DARK MATTER COSMOLOGY   ***
! ****************************************

! This routine first sets up a dark matter scalar
! field cosmological background. It can also add a perturbation.
!
! For the perturbation we assume that the momentum constraint
! is trivially satisfied, and solve the Hamiltonian
! constraint for the conformal factor.
!
! NOTE FOR PARALLEL RUNS:  The initial data is not really
! solved in parallel.  It is in fact solved only on processor
! zero on a full size array, and then it is distributed
! among the other processors.  This is slow, but works.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use derivatives
  use radialfunctions

! Extra variables.

  implicit none

  integer i,l,p,iter
  integer lmin,lmax
  integer i0,imax,Naux
  integer status(MPI_STATUS_SIZE)
  integer :: maxiter = 100

  real(8) zero,one,two,smallpi
  real(8) aux,res
  real(8) :: epsilon = 1.d-10

  real(8), dimension (0:Nl-1,1-ghost:Nrmax)   :: CL0,CL1,CL2,VL
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: C0,C1,C2,CA,VG
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: u,u_old

! u        Global solution of ODE.
! CL0,C0   Local and global coefficient of first derivative in ODE.
! CL1,C1   Local and global coefficient of linear term in ODE.
! CL2,C2   Local and global source term in ODE.
! CA       Auxiliary global array.
! VL,VG    Local and global potential.


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  one  = 1.d0
  two  = 2.d0

  smallpi = acos(-one)


! **********************
! ***   BACKGROUND   ***
! **********************

! For the background we set up a uniform scalar field.
! We then set up a uniform trK in different ways
! depending on the form of the potential.
!
! Notice that if we want an expanding universe
! we must use a negative value for trK.

  scalar_phi = scalar_bg_phi0
  scalar_xi  = zero

! I) Zero potential.

  if (scalarpotential=="none") then

!    For the case of zero potential we must take the time
!    derivative of the scalar field different from zero,
!    as otherwise we end up with a trivial solution (Minkowski).
!    The Hamiltonian constraint in this case reduces to:
!
!    trK  =  - sqrt(24 pi rho)  =  - sqrt(12 pi scalar_pi**2)

     scalar_pi = scalar_bg_pi0
     trK = - sqrt(12.d0*smallpi)*abs(scalar_bg_pi0)

! II) Non-zero potential.

  else

!    For the case of non-zero potential we can take the
!    time derivative of the scalar field equal to zero.
!    The Hamiltonian constraint then reduces to:
!
!    trK  =  - sqrt(24 pi rho)  =  - sqrt(24 pi V)
!
!    with V the scalar field potential.

     scalar_pi = zero

     do l=0,Nl-1
        call potential(l)
     end do

     trK = - dsqrt(24.d0*smallpi*scalar_V)

  end if

! Cosmological background quantities.

  if (cosmic_run) then

!    Geometry.

     cosmobg_trK = trK(0,Nr)
     cosmobg_H = - cosmobg_trK/3.d0

!    Scalar field.

     cosmobg_scalar_phi = scalar_bg_phi0

     if (scalarpotential=="none") then
        cosmobg_scalar_pi = scalar_bg_pi0
     else
        cosmobg_scalar_pi = zero
     end if

  end if


! ************************
! ***   PERTURBATION   ***
! ************************

! Add a perturbation to the background.

  if (scalar_bg_pert) then

!    Gaussian profile.

     if (scalarprofile=="gaussian") then

        if (scalar_r0==0.d0) then
           scalar_phi = scalar_phi + gaussian(scalar_a0,0.d0,scalar_s0)
        else
           scalar_phi = scalar_phi + gaussian(scalar_a0,+scalar_r0,scalar_s0) &
                                   + gaussian(scalar_a0,-scalar_r0,scalar_s0)
        end if

!    Smooth top-hat profile.

     else if (scalarprofile=="tophat") then

        if (scalar_r0==0.d0) then
           scalar_phi = scalar_phi + tophat(scalar_a0,0.d0,scalar_s0,scalar_t0)
        else
           scalar_phi = scalar_phi + tophat(scalar_a0,+scalar_r0,scalar_s0,scalar_t0) &
                                   + tophat(scalar_a0,-scalar_r0,scalar_s0,scalar_t0)
        end if

     end if

!    The spatial derivative is calculated with finite differences.

     diffvar => scalar_phi

     do l=0,Nl-1
        scalar_xi(l,:) = diff1(l,+1)
     end do


!    *****************************
!    ***   SOLVE HAMILTONIAN   ***
!    *****************************

!    Now we need to solve the Hamiltonian constraint for the
!    conformal factor, which in this case takes the form:
!
!    __2                          2          5
!    \/ psi  +  ( 2 pi rho  -  trK / 12 ) psi  =  0
!
!
!    with rho the energy density of the scalar field given by:
!
!                2      2       4
!    rho  =  ( Pi  +  Xi / A psi ) / 2  +  V
!
!
!    Notice that if we assume that trK is uniform and solves the
!    background, then the Hamiltonian constraint reduces to:
!
!    __2              2                    5
!    \/ psi  +  pi  Xi  psi  +  2 pi dV psi  =  0
!
!
!    with dV the difference between the full potential and the
!    background potential V0:
!
!    dV  =  V  -  V0
!
!    In the expressions above, the Laplacian of psi is given by:
!
!    __2           2
!    \/ psi  =  [ d psi  +  (2/r - d A/2A + d B/B) d psi ] / A
!                  r                r        r      r
!
!    In most cases we will have A=B=1, so the Laplacian simplifies.
!    But this won't be the case for a transformed radial coordinate.
!
!    Notice that here we are assuming that there is no perturbation
!    to the time derivative.
!
!    Notice also that V0 = trK**2/(24 pi), so that:
!
!                                 2
!    dV  =  V  -  V0  =  V  -  trK / (24 pi)


!    **************************
!    ***   ZERO POTENTIAL   ***
!    **************************

!    This case is actually more complicated since we have
!    a non-trivial background time derivative, and with the
!    perturbation we also have a non-trivial spatial derivative.
!    This implies that the momentum constraint is non-trivial
!    and also needs to be solved. I will then leave this for later.

     if (scalarpotential=="none") then

        if (rank==0) then
           print *, 'Cosmological perturbation for a scalar field with zero potential not yet implemented.'
           print *, 'Aborting!  (subroutine idata_scalarDM)'
        end if

        call die


!    ******************************
!    ***   NON-ZERO POTENTIAL   ***
!    ******************************

!    In this case the time derivative is zero so the momentum
!    constraint is trivial, so we just solve the Hamiltonian.

     else

!       Message to screen.

        if (rank==0) then
           print *, 'Solving initial data for perturbed scalar dark matter pulse with non-zero potential ...'
           print *
         end if

!       Find scalar field potential.

        do l=0,Nl-1
           call potential(l)
        end do

!       Array size.

        Naux = (Nrmax + ghost)

!       Fill in local coefficients of linearized equation.

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

!       Fill in local potential (remember to subtract background value).

        VL = scalar_V - trK**2/(24.d0*smallpi)

!       Processor zero receives data and solves ODE.

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

!          Iterate over other processors to receive data.

           i0 = imax

           do p=1,size-1

!             Receive local coefficients (CL0,CL1) from other processors
!             and copy them into large arrays (C0,C1).

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

!          For initial guess solve with zero potential.

           C2 = 0.d0

           lmin = 0; lmax = Nl-1
           call invertmatrix(lmin,lmax,1.d0,u,C0,C1,C2,+1,"robin")

!          Initialize residual.

           res = 1.d0

!          Begin iterations.

           iter = 0

           do while ((res.gt.epsilon).and.(iter.lt.maxiter))

              iter = iter + 1

!             Save old solution.

              u_old = u

!             Correct linear term and source term.  For this I am using
!             a variation of Newton's method.  We first notice that the
!             equation has the form of a linear term A*u plus a non-linear
!             purely diagonal term F(u):
!
!             A*u  +  F(u)  =  0
!
!             We can then expand in Taylor the non-linear term to first
!             order around a previous guess u_old to find:
!
!             A*u  +  F(u_old)  +  F'(u_old)*(u - u_old)  =  0
!
!             and express the iterations as:
!
!             [A + F'(u_old)]*u_new = F'(u_old)*u_old - F(u_old)
!
!             The effect is then just to add to the diagonal of the
!             original matrix A the term F'(u_old), and introduce a
!             diagonal source term on the right hand side. This has
!             the advantage that I don't need to calculate A*u_old.
!
!             In our case F(u) = 2*pi*V*u**5, so that:  F'(u) = 10*pi*V*u**4.
!             Notice also that:  F'(u_old)*u_old - F(u_old)  =  8*pi*V*u_old**5

              if (newr) then
                 CA = C1 + 10.d0*smallpi*A*VG*u_old**4
                 C2 = 8.d0*smallpi*A*VG*u_old**5
              else
                 CA = C1 + 10.d0*smallpi*VG*u_old**4
                 C2 = 8.d0*smallpi*VG*u_old**5
              end if

!             Call matrix inversion.

              lmin = 0; lmax = Nl-1
              call invertmatrix(lmin,lmax,1.d0,u,C0,CA,C2,+1,"robin")

!             Find residual (across all grid levels).

              res = 0.d0

              do l=0,Nl-1
                 do i=1,Nrtotal
                    aux = abs(u(l,i)-u_old(l,i))
                    if (aux>res) res = aux
                 end do
              end do

              write(*,"(A,I0,A,ES8.2)") ' Iteration: ',iter,'    Residual: ',res

           end do

!          Message in case we reached the maximum iteration number.

           if (iter>=maxiter) then
              print *
              print *, 'Maximum iteration number reached in idata_scalarDM.f90.'
              print *, 'Initial data solver did not converge, aborting!'
              print *
           else
              print *
              print *, 'Done!'
              print *
           end if

!       Other processors send local coefficients (CL0,CL1) to processor 0.

        else

           do l=0,Nl-1
              call MPI_SEND(CL0(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
              call MPI_SEND(CL1(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
              call MPI_SEND(VL(l,:) ,Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           end do

        end if

     end if


!    ************************************************
!    ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
!    ************************************************

!    For parallel runs, when we get here the solution
!    is known only on processor zero for the full size
!    array with dimensions Nrtotal.  We must now distribute
!    the solution among all other processors.

     if (size==1) then
        psi = u
     else
        call distribute(0,Nl-1,psi,u)
     end if

!    Derivative of psi.

     diffvar => psi

     do l=0,Nl-1
        D1_psi(l,:) = diff1(l,+1)
     end do


!    ***************************************
!    ***   FIND CONFORMAL FUNCTION phi   ***
!    ***************************************

!    Find phi.

     phi  = dlog(psi)
     D1_phi = D1_psi/psi

     if (chimethod) then
        chi  = one/psi**chipower
        D1_chi = - dble(chipower)*D1_psi/psi**chipower
     end if

!    Find psi2 and psi4.

     psi2 = psi**2
     psi4 = psi**4

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine idata_scalarDM
