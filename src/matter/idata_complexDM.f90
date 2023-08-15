!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_complexDM.f90,v 1.27 2021/11/29 19:14:02 malcubi Exp $

  subroutine idata_complexDM

! ************************************************
! ***   COMPLEX SCALAR DARK MATTER COSMOLOGY   ***
! ************************************************

! This routine first sets up a dark matter complex scalar
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

  real(8) zero,half,one,two,smallpi
  real(8) aux,res
  real(8) :: epsilon = 1.d-10

  real(8), dimension (0:Nl-1,1-ghost:Nrmax)   :: CL0,CL1,CL2,VL,V0
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
  half = 0.5d0
  one  = 1.d0
  two  = 2.d0

  smallpi = acos(-one)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for complex scalar field dark matter cosmology ...'
     print *
  end if


! **********************
! ***   BACKGROUND   ***
! **********************

! For the background we set up a uniform field.
! We then set up a uniform trK in different ways
! depending on the form of the potential.
!
! Notice that if we want an expanding universe
! we must use a negative value for trK.

  complex_phiR = complex_bg_phiR0
  complex_phiI = complex_bg_phiI0

  complex_xiR = zero
  complex_xiI = zero

! For the complex scalar field we assume a harmonic
! time dependence, with a frequency equal to the mass.
!
! Since we are assuming a harmonic time dependence, the
! mass can't be zero.  If the potential is zero, then
! the mass works as an ad-hoc frequency.

  if ((complexpotential=="none").and.(complex_mass==0.d0)) then
     if (rank==0) then
        print *, 'For a cosmological perturbation of a complex scalar field with zero potential'
        print *, 'the parameter complex_mass is used as an ad-hoc frequency and must be non-zero,'
        print *, 'as otherwise the background has zero energy density.'
        print *, 'Aborting!  (subroutine idata_complexDM)'
     end if
  end if

  complex_piR = - complex_mass*complex_bg_phiI0
  complex_piI = + complex_mass*complex_bg_phiR0

! Find background potential and save it in V0.

  do l=0,Nl-1
     call potential(l)
  end do

  V0 = complex_V

! The Hamiltonian constraint has the form:
!
!                                                  2      2
! trK  =  - sqrt(24 pi rho)  =  - sqrt[ 24 pi ( (Pi_r + Pi_i)/2 + V ) ]
!
! with V the complex field potential.

  trK = - sqrt(24.d0*smallpi*(half*(complex_piR**2 + complex_piI**2) + V0))

! Cosmological background quantities.

  if (cosmic_run) then

!    Geometry.

     cosmobg_trK = trK(0,Nr)
     cosmobg_H = - cosmobg_trK/3.d0

!    Complex field.

     cosmobg_complex_phiR = complex_bg_phiR0
     cosmobg_complex_phiI = complex_bg_phiI0

!    Time derivative.

     cosmobg_complex_piR = - complex_mass*complex_bg_phiI0
     cosmobg_complex_piI = + complex_mass*complex_bg_phiR0

  end if


! ************************
! ***   SANITY CHECK   ***
! ************************

! Check that the momentum constraint is trivially satisfied.

  if (complex_bg_phiR0*complex_bg_phiI0/=zero) then

     if (rank==0) then
        print *, 'For a cosmological perturbation of a complex scalar field'
        print *, 'either the real or the imaginary background values must be zero.'
        print *, 'Aborting!  (subroutine idata_complexDM)'
     end if

     call die

  else if ((complex_bg_phiR0/=zero).and.(complexI_a0/=zero)) then

     if (rank==0) then
        print *, 'For a cosmological perturbation of a complex scalar field with'
        print *, 'real background, and trivially satisfied momentum constraint,'
        print *, 'the imaginary part of the perturbation must be zero.'
        print *, 'Aborting!  (subroutine idata_complexDM)'
     end if

     call die

  else if ((complex_bg_phiI0/=zero).and.(complexR_a0/=zero)) then

     if (rank==0) then
        print *, 'For a cosmological perturbation of a complex scalar field with'
        print *, 'imaginary background, and trivially satisfied momentum constraint,'
        print *, 'the real part of the perturbation must be zero.'
        print *, 'Aborting!  (subroutine idata_complexDM)'
     end if

     call die

  end if


! ************************
! ***   PERTURBATION   ***
! ************************

! Add a perturbation to the background.

  if (complex_bg_pert) then

!    Gaussian profile.

     if (complexprofile=="gaussian") then

        if (complexR_r0==0.d0) then
           complex_phiR = complex_phiR + gaussian(complexR_a0,0.d0,complexR_s0)
        else
           complex_phiR = complex_phiR + gaussian(complexR_a0,+complexR_r0,complexR_s0) &
                                       + gaussian(complexR_a0,-complexR_r0,complexR_s0)
        end if

        if (complexI_r0==0.d0) then
           complex_phiI = complex_phiI + gaussian(complexI_a0,0.d0,complexI_s0)
        else
           complex_phiI = complex_phiI + gaussian(complexI_a0,+complexI_r0,complexI_s0) &
                                       + gaussian(complexI_a0,-complexI_r0,complexI_s0)
        end if

!    Smooth top-hat profile.

     else if (complexprofile=="tophat") then

        if (complexR_r0==0.d0) then
           complex_phiR = complex_phiR + tophat(complexR_a0,0.d0,complexR_s0,complexR_t0)
        else
           complex_phiR = complex_phiR + tophat(complexR_a0,+complexR_r0,complexR_s0,complexR_t0) &
                                       + tophat(complexR_a0,-complexR_r0,complexR_s0,complexR_t0)
        end if

        if (complexI_r0==0.d0) then
           complex_phiI = complex_phiI + tophat(complexI_a0,0.d0,complexI_s0,complexI_t0)
        else
           complex_phiI = complex_phiI + tophat(complexI_a0,+complexI_r0,complexI_s0,complexI_t0) &
                                       + tophat(complexI_a0,-complexI_r0,complexI_s0,complexI_t0)
        end if

     end if

!    The spatial derivatives are calculated with finite differences.

     diffvar => complex_phiR

     do l=0,Nl-1
        complex_xiR(l,:) = diff1(l,+1)
     end do

     diffvar => complex_phiI

     do l=0,Nl-1
        complex_xiI(l,:) = diff1(l,+1)
     end do

!    Harmonic time dependence.

     complex_piR = - complex_mass*complex_phiI
     complex_piI = + complex_mass*complex_phiR


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
!    with "rho" the energy density of the scalar field given by:
!
!                 2        2         2        2          4
!    rho  =  [  Pi_r  +  Pi_i  +  (Xi_r  +  Xi_i) / A psi  ] / 2  +  V 
!
!
!    Notice that if we assume that trK is uniform and solves the
!    background, then the Hamiltonian constraint simplifies to:
!
!    __ 2              2        2                                       5
!    \/ psi  +  pi ( Xi_r  +  Xi_i ) psi / A  +  2 pi ( dV  +  dPi ) psi  =  0
! 
!
!    with dV the difference between the full potential and the
!    background potential V0:
!
!    dV  :=  V  -  V0
!
!    And where we have defined:
!
!                      2      2       2      2
!    dPi :=  (1/2) [ Pi_r + Pi_i - (Pi_r + Pi_i)          ]
!                                               background
!
!    which using the harmonic time dependence becomes:
!
!              2         2       2        2       2
!    dPi :=  (m /2) [ Phi_r + Phi_i - (Phi_r + Phi_i)          ]
!                                                   background
!
!    In the expressions above, the Laplacian of psi is given by:
!
!    __2           2
!    \/ psi  =  [ d psi  +  (2/r - d A/2A + d B/B) d psi ] / A
!                  r                r        r      r
!
!    In most cases we will have A=B=1, so the Laplacian simplifies.
!    But this won't be the case for a transformed radial coordinate.

!    Message to screen.

     if (rank==0) then
        print *, 'Solving initial data for a perturbed complex dark matter pulse ...'
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
        CL1 = smallpi*(complex_xiR**2 + complex_xiI**2) &
            + 0.25d0/r**2*(one - A/B) - 0.25d0*D1_A/(r*A) + 0.75d0*D1_B/(r*B) &  ! This is the Ricci scalar, it should
            - 0.125d0*D1_A*D1_B/(A*B) - 0.0625d0*D1_B**2/B**2 + 0.25d0*D2_B/B    ! be zero for conformally flat data.
     else
        CL0 = two/r
        CL1 = smallpi*(complex_xiR**2 + complex_xiI**2)
     end if

     CL2 = 0.d0

!    Fill in local potential. Remember to subtract the background
!    and add contributions from time derivatives of the field.

     VL = (complex_V - V0) &
        + half*(complex_piR**2 + complex_piI**2 - complex_mass**2*(complex_bg_phiR0**2 + complex_bg_phiI0**2))

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

!       For initial guess solve taking phi=1 on the right-hand side.

        C2 = - 2.d0*smallpi*VG

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
!          Notice also that:  F'(u_old)*u_old - F(u_old) = 8*pi*V*u_old**5

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
           print *, 'Maximum iteration number reached in idata_complexDM.f90.'
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


!    ***************************************
!    ***   FIND CONFORMAL FUNCTION phi   ***
!    ***************************************

!    Derivative of psi.

     diffvar => psi

     do l=0,Nl-1
        D1_psi(l,:) = diff1(l,+1)
     end do

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

  end subroutine idata_complexDM
