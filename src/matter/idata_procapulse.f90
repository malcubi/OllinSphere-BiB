
  subroutine idata_procapulse

! ************************************
! ***   PROCA PULSE INITIAL DATA   ***
! ************************************

! This initial data leaves the lapse, shift, conformal
! spatial metric and extrinsic curvature as in Minkowski
! (which have already been set up in "initial").
!
! For this initial data we need to solve first the Gauss
! constraint, which for a Proca field has the form:
!
!    i     2
! D E  +  m  Phi  =  0
!  i
!
! with E^i the electric Proca field and Phi the scalar
! potential.
!
! We do this by assuming we have a known conformal
! metric, and rescale the electric field and scalar
! potential as:
!
!  i    ^i      6                  ^       6
! E  =  E  / psi           Phi  = Phi / psi
!
! We then assume that the conformal electric field
! is minus the gradient of an auxiliary potential V,
! so the Gauss constraint becomes:
!
! __2        2  ^
! \/  V  =  m  Phi
!
! where we take the conformal space Laplacian.
! 
! We can then choose the conformal scalar potential Phi
! freely and solve for V.  Once we have V we find its
! gradient to obtain the conformal electric field:
!
! ^r       r
! E  =  - d V  =  - d V / A
!                    r
!
! We then assume that we have time-symmetric initial
! data so that the extrinsic curvature and momentum flux
! are initially zero. This implies that the vector
! potential also vanishes.
!
! We subsitute the solution of the Gauss constraint
! into the Hamiltonian constraint and solve for the
! conformal factor:
!
! __2                5
! \/ psi  +  2 pi psi rho  =  0
!
!
! with the energy density of the Proca field given by:
!
!                       r     2    2
! rho  =  1/(8 pi) [ E E  +  m  Phi ]
!                     r
!
!                      ^ 2     8     2  ^ 2      12
!      =  1/(8 pi) [ A E  / psi  +  m  Phi  / psi  ]
!
! where now E := E^r, so that:
!
!         5                   ^ 2     3     2  ^ 2      7
! 2 pi psi  rho  =  (1/4) [ A E  / psi  +  m  Phi  / psi  ]
!
!
! Both in the Gauss and Hamiltonian constraints the Laplacian
! operator has the general form:
!
! __2           2
! \/ psi  =  [ d psi  +  (2/r - d A/2A + d B/B) d psi ] / A
!               r                r        r      r
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

  real(8) idr
  real(8) aux,res
  real(8) zero,one,two
  real(8) :: epsilon = 1.d-10

  real(8), dimension (0:Nl-1,1-ghost:Nrmax)   :: CL0,CL1,CL2,CLA,AL
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: C0,C1,C2
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: AG,VG,EG,PhiG
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: PsiG,Psi_old

! VG       Global auxiliary potential.
! EG       Global electric field.
! PhiG     Global scalar potential.

! PsiG     Global conformal factor.
! AL,AG    Local and global radial metric coefficient.

! CL0,C0   Local and global coefficient of first derivative in ODE.
! CL1,C1   Local and global coefficient of linear term in ODE.
! CL2,C2   Local and global source term in ODE.
! CLA      Auxiliary local array.


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0 
  one  = 1.d0
  two  = 2.d0


! *****************************
! ***   INITIALIZE ARRAYS   ***
! *****************************

  VG   = zero
  EG   = zero
  PhiG = zero

  AG   = one
  PsiG = one


! **************************************************************
! ***   SET UP INITIAL PULSE IN CONFORMAL SCALAR POTENTIAL   ***
! **************************************************************

! Message to screen.

  if (rank==0) then
     print *, 'Solving initial data for a Proca pulse ...'
  end if

! Initial pulse.  The scalar potential must be even.

! Gaussian profile.

  if (procaprofile=="gaussian") then

     if (proca_r0==0.d0) then
        procaPhi = gaussian(proca_a0,0.d0,proca_s0)
     else
        procaPhi = gaussian(proca_a0,+proca_r0,proca_s0) &
                 + gaussian(proca_a0,-proca_r0,proca_s0)
     end if

! Smooth top-hat profile.

  else if (procaprofile=="tophat") then

     if (proca_r0==0.d0) then
        procaPhi = tophat(proca_a0,0.d0,proca_s0,proca_t0)
     else
        procaPhi = tophat(proca_a0,+proca_r0,proca_s0,proca_t0) &
                 + tophat(proca_a0,-proca_r0,proca_s0,proca_t0)
     end if

  end if

! Remember that the vector potential is zero.

  procaA = 0.d0


! *************************************
! ***   FILL IN GLOBAL ARRAY PhiG   ***
! *************************************

! Single processor run.

  if (size==1) then

     phiG = procaPhi

! Parallel run.

  else

     Naux = (Nrmax + ghost)

     if (rank==0) then

        if (size==1) then
           imax = Nrl(0)
        else
           imax = Nrl(0) - ghost
        end if

        do i=1-ghost,imax
           PhiG(:,i) = procaPhi(:,i)
        end do

!       Iterate over other processors to receive data.

        i0 = imax

        do p=1,size-1

           do l=0,Nl-1
              call MPI_RECV(CLA(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           end do

           if (p==size-1) then
              imax = Nrl(p)
           else
              imax = Nrl(p) - ghost
           end if

           do i=1,imax
              PhiG(:,i+i0) = CLA(:,i)
           end do

           i0 = i0 + imax

        end do

!    Other processors send data.

     else

        do l=0,Nl-1
           call MPI_SEND(procaPhi(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end do

     end if

  end if


! ********************************************
! ***   GLOBAL RADIAL METRIC COEFFICIENT   ***
! ********************************************

! In case the initial conformal metric is not flat,
! we need to copy the global metric array A.

! Copy A into AL.

  AL = A

! Processor 0.

  if (rank==0) then

     if (size==1) then
        imax = Nrl(0)
     else
        imax = Nrl(0) - ghost
     end if

     do i=1-ghost,imax
        AG(:,i) = AL(:,i)
     end do

!    Iterate over other processors to receive data.

     i0 = imax

     do p=1,size-1

!       Receive local coefficients from other processors
!       and copy them into large arrays.

        do l=0,Nl-1
           call MPI_RECV(AL(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        end do

        if (p==size-1) then
           imax = Nrl(p)
        else
           imax = Nrl(p) - ghost
        end if

        do i=1,imax
           AG(:,i+i0) = AL(:,i)
        end do

        i0 = i0 + imax

     end do

! Other processors send local coefficients to processor 0.

  else

     do l=0,Nl-1
        call MPI_SEND(AL(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
     end do

  end if


! **********************************
! ***   SOLVE GAUSS CONSTRAINT   ***
! **********************************

! The Gauss constraint takes the form:
!
!  2                     2  ^
! d  V  + (2/r) d V  =  m  Phi
!  r             r
!
! We solve this by direct matrix inversion.

  if (rank==0) then
     print *, 'Solving Gauss constraint ...'
  end if

! Array size.

  Naux = (Nrmax + ghost)

! Fill in (local) coefficients of linear equation.  We take
! into account the possibility that A and B are non-trivial.

  CL0 = two/r - 0.5d0*D1_A/A + D1_B/B
  CL1 = 0.d0
  CL2 = A*proca_mass**2*procaPhi

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
        C2(:,i) = CL2(:,i)
     end do

!    Iterate over other processors to receive data.

     i0 = imax

     do p=1,size-1

!       Receive local coefficients from other processors
!       and copy them into large arrays.

        do l=0,Nl-1
           call MPI_RECV(CL0(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(CL1(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(CL2(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        end do

        if (p==size-1) then
           imax = Nrl(p)
        else
           imax = Nrl(p) - ghost
        end if

        do i=1,imax
           C0(:,i+i0) = CL0(:,i)
           C1(:,i+i0) = CL1(:,i)
           C2(:,i+i0) = CL2(:,i)
        end do

        i0 = i0 + imax

     end do

!    Call ODE solver. Notice that in this case C1 is always 0.

     lmin = 0; lmax = Nl-1

     call invertmatrix(lmin,lmax,0.d0,VG,C0,C1,C2,+1,"robin")

! Other processors send local coefficientes (CL0,CL1,CL2) to processor 0.

  else

     do l=0,Nl-1
        call MPI_SEND(CL0(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(CL1(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(CL2(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
     end do

  end if


! *****************************************
! ***   FIND CONFORMAL ELECTRIC FIELD   ***
! *****************************************

! The conformal electric field is minus the gradient
! of the auxiliary potential V.   We take into account
! the possibility that A is non-trivial when we raise
! the index to find E^r.

  do l=0,Nl-1

     idr = 1.d0/dr(l)

     if (order=="two") then

!       Interior points.

        do i=1,Nrtotal-1
           EG(l,i) = - 0.5d0*idr*(VG(l,i+1) - VG(l,i-1))/AG(l,i)
        end do

!       Boundary with one-sided differences.

        EG(l,Nrtotal) = - 0.5d0*idr*(3.d0*VG(l,Nrtotal) - 4.d0*VG(l,Nrtotal-1) + VG(l,Nrtotal-2))/AG(l,Nrtotal)

     else

!       Interior points.

        do i=1,Nrtotal-2
           EG(l,i) = - idr*(8.d0*(VG(l,i+1) - VG(l,i-1)) - (VG(l,i+2) - VG(l,i-2)))/12.d0/AG(l,i)
        end do

!       Boundary. We use sixth-order left-sided differences
!       to improve the boundary behavior.

        EG(l,Nrtotal-1) = - idr*(1.d0/6.d0*VG(l,Nrtotal) + 77.d0/60.d0*VG(l,Nrtotal-1) &
           - 5.d0/2.d0*VG(l,Nrtotal-2) + 5.d0/3.d0*VG(l,Nrtotal-3) - 5.d0/6.d0*VG(l,Nrtotal-4) &
           + 0.25d0*VG(l,Nrtotal-5) - 1.d0/30.d0*VG(l,Nrtotal-6))/AG(l,Nrtotal-1)

        EG(l,Nrtotal  ) = - idr*(49.d0/20.d0*VG(l,Nrtotal) - 6.d0*VG(l,Nrtotal-1) &
           + 7.5d0*VG(l,Nrtotal-2) - 20.d0/3.d0*VG(l,Nrtotal-3) + 15.d0/4.d0*VG(l,Nrtotal-4) &
           - 6.d0/5.d0*VG(l,Nrtotal-5) + 1.d0/6.d0*VG(l,Nrtotal-6))/AG(l,Nrtotal)

     end if

!    Ghost points at origin.

     do i=1,ghost
        EG(l,1-i) = - EG(l,i)
     end do

  end do


! ****************************************
! ***   SOLVE HAMILTONIAN CONSTRAINT   ***
! ****************************************

! The Hamiltonian constraint has the form:
!
! __2                  5
! \/ psi  =  - 2 pi psi rho
!
! where the Laplacian of psi is:
!
! __2           2
! \/ psi  =  [ d psi  +  (2/r - d A/2A + d B/B) d psi ] / A
!               r                r        r      r
!
! and the energy denisty term is given by:
!
!         5                     ^ 2     3     2  ^ 2      7
! 2 pi psi  rho  =  (1/4) [ A  E  / psi  +  m  Phi  / psi  ]
!
! This equation is clearly not linear in psi, we solve it
! using Newton's iterative method.

  if (rank==0) then
     print *, 'Solving Hamiltonian constraint ...'
     print *
  end if

! Array size.

  Naux = (Nrmax + ghost)

! Fill in local coefficient CL0 of linearized equation.
! The other coefficients are calculated globally later.
! We take into account the possibility that A and B
! are non-trivial.

  CL0 = two/r - 0.5d0*D1_A/A + D1_B/B

! Processor zero receives data and solves ODE.

  if (rank==0) then

     if (size==1) then
        imax = Nrl(0)
     else
        imax = Nrl(0) - ghost
     end if

     do i=1-ghost,imax
        C0(:,i) = CL0(:,i)
     end do

!    Iterate over other processors to receive data.

     i0 = imax

     do p=1,size-1

!       Receive local coefficient CL0 from other processors
!       and copy it into large array.

        do l=0,Nl-1
           call MPI_RECV(CL0(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        end do

        if (p==size-1) then
           imax = Nrl(p)
        else
           imax = Nrl(p) - ghost
        end if

        do i=1,imax
           C0(:,i+i0) = CL0(:,i)
        end do

        i0 = i0 + imax

     end do

!    Solve for initial guess taking psi=1.

     C1 = 0.d0
     C2 = - 0.25d0*AG*(A*EG**2 + proca_mass**2*PhiG**2)

     lmin = 0; lmax = Nl-1
     call invertmatrix(lmin,lmax,1.d0,PsiG,C0,C1,C2,+1,"robin")

!    Initialize residual.

     res = 1.d0

!    Begin iterations.

     iter = 0

     do while ((res.gt.epsilon).and.(iter.lt.maxiter))

        iter = iter + 1

!       Save old solution.

        Psi_old = PsiG

!       Correct linear source term.  For this I am using a
!       variation of Newton's method.  We first notice that
!       the equation has the form of a linear term A*u plus
!       a non-linear purely diagonal term F(u):
!
!       A*u  =  F(u)
!
!       We can then expand in Taylor the non-linear term to
!       first order around a previous guess u_old to find:
!
!       A*u  =  F(u_old)  +  F'(u_old)*(u - u_old)
!
!       and express the iterations as:
!
!       [A - F'(u_old)]*u_new = F(u_old) - F'(u_old)*u_old
!
!       The effect is then just to add to the diagonal of the
!       original matrix A the term -F'(u_old), and introduce a
!       diagonal source term on the right hand side. This has
!       the advantage that I don't need to calculate A*u_old.
!
!       In our case we have:
!
!       F(u)  =  - 1/4 ( E^2 / u^3  +  m^2 Phi^2 / u^7 )
!
!             =  - 1/4 ( E^2  +  m^2 Phi^2 / u^4 ) / u^3
!
!       F'(u) =  + 1/4 ( 3 E^2 / u^4  +  7 m^2 Phi^2 / u^8 )
!
!             =  1/(4 u^4) ( 3 E^2  +  7 m^2 Phi^2 / u^4 )
!
!       Notice also that:
!
!       F'(u_old)*u_old  =  1/(4 u_old^3) ( 3 E^2 + 7 m^2 Phi^2 / u_old^4 )
!
!       so that:
!
!       F(u_old) - F'(u:old)*u_old  =  - 1/u_old^3 ( E^2 + 2 m^2 Phi^2 / u_old^4 )

        C1 = - 0.25d0*AG*(3.d0*AG*EG**2 + 7.d0*proca_mass**2*PhiG**2/Psi_old**4)/Psi_old**4
        C2 = - AG*(AG*EG**2 + 2.d0*proca_mass**2*PhiG**2/Psi_old**4)/Psi_old**3

!       Call matrix inversion.

        lmin = 0; lmax = Nl-1
        call invertmatrix(lmin,lmax,1.d0,PsiG,C0,C1,C2,+1,"robin")

!       Find residual (across all grid levels).

        res = 0.d0

        do l=0,Nl-1
           do i=1,Nrtotal
              aux = abs(PsiG(l,i)-Psi_old(l,i))
              if (aux>res) res = aux
           end do
        end do

        write(*,"(A,I0,A,ES8.2)") ' Iteration: ',iter,'    Residual: ',res

     end do

!    Message in case we reached the maximum iteration number.

     if (iter>=maxiter) then
        print *
        print *, 'Maximum iteration number reached in idata_procapulse.f90.'
        print *, 'Initial data solver did not converge, aborting!'
        print *
     else
        print *
        print *, 'Done!'
        print *
     end if

! Other processors send local coefficients to processor 0.

  else

     do l=0,Nl-1
        call MPI_SEND(CL0(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
     end do

  end if


! ************************************************
! ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
! ************************************************

  if (size==1) then
     procaV = VG
     procaE = EG
     psi    = PsiG
  else
     call distribute(0,Nl-1,procaV,VG)
     call distribute(0,Nl-1,procaE,EG)
     call distribute(0,Nl-1,psi,PsiG)
  end if

! Now rescale electric field and scalar potential.

  procaE   = procaE/psi**6
  procaPhi = procaPhi/psi**6


! ***************************************
! ***   FIND CONFORMAL FUNCTION phi   ***
! ***************************************

! Find derivative of psi.

  diffvar => psi

  do l=0,Nl-1
     D1_psi(l,:) = diff1(l,+1)
  end do

! Find phi.

  phi = dlog(psi)
  D1_phi = D1_psi/psi

  if (chimethod) then
     chi = one/psi**chipower
     D1_chi = - dble(chipower)*D1_psi/psi**chipower
  end if

! Find psi2 and psi4.

  psi2 = psi**2
  psi4 = psi**4


! ***************
! ***   END   ***
! ***************

  end subroutine idata_procapulse
