
  subroutine alphamaximal(lmin,lmax,maxbound,alpha0)

! ***************************
! ***   MAXIMAL SLICING   ***
! ***************************

! This routine solves the elliptic equation:
!
! __2                  ij
! \/ alpha  =  alpha  K  K   +  4 pi alpha ( rho + S )
!                         ij
!
! where the Laplacian of alpha is given in general by:
!
! __2            2                                                            4
! \/ alpha  = [ d alpha  -  d alpha (d A/2A - d B/B - 2 d phi - 2/r) ] / A psi
!                r           r        r        r         r
!
! It is important to notice that, since this routine is
! called BEFORE we calculate all auxiliary geometric
! quantities, here one should not assume that any of
! those quantities has been correctly defined yet.
!
! The boundary conditions at the moment can be:
!       robin      (alpha0 is the value at infinity)
!       dirichlet  (alpha0 is the boundary value)
!
! Notice that we might be solving at an intermediate
! time level where not all grids coincide. Or we can be
! solving at a coarse level when the fine levels have not
! yet reached that point. That is the purpose of the
! parameters "lmin" and "lmax".
!
! WARNING: The current maximal slicing routine uses
! a direct matrix inversion algorithm for band-diagonal
! matrices.  The problem with this is that it requires
! an outer boundary condition.  At the moment I have
! included Dirichlet and Robin boundary conditions.
! This, though in principle correct, might not be ideal
! in some cases. For example, for trumpet black hole initial
! data the Robin boundary condition turns out to be
! incompatible with the rest of the solution (it gave me
! a big headache and took me several days to discover this).
!
! It might be a good idea to code a Runge-Kutta maximal
! solver that starts from the origin.  I will do that
! later.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays
  use derivatives

! Extra variables.

  implicit none

  logical contains

  integer lmin,lmax,l,p,i
  integer i0,imax,Naux
  integer status(MPI_STATUS_SIZE)

  real(8) smallpi,alpha0
  real(8) aux

  real(8), dimension (0:Nl-1,1-ghost:Nrmax)   :: CL0,CL1
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: u,C0,C1,C2

  character maxbound*(*),bound*10


! *******************
! ***   NUMBERS   ***
! *******************

  smallpi = acos(-1.d0)


! **********************************************************
! ***   MAKE SURE WE HAVE CORRECT GEOMETRIC QUANTITIES   ***
! **********************************************************

! Conformal factor.

  if (chimethod) then
     do l=lmin,lmax
        phi(l,:) = - log(abs(chi(l,:)))/dble(chipower)
     end do
  end if

! Derivatives of conformal factor.

  diffvar => phi

  if (dorigin=="onesided") then
     do l=lmin,lmax
        D1_phi(l,:) = diff1(l,+1,.true.)
     end do
  else
     do l=lmin,lmax
        D1_phi(l,:) = diff1(l,+1)
     end do
  end if

! Derivatives of spatial metric.

  diffvar => A

  do l=lmin,lmax
     D1_A(l,:) = diff1(l,+1)
  end do

  diffvar => B

  do l=lmin,lmax
     D1_B(l,:) = diff1(l,+1)
  end do

! Square of extrinsic curvature. Remember that KTB = - KTA/2,
! and also that for maximal slicing trK=0.

  if (noKTA) then
     do l=lmin,lmax
        KTA(l,:) = 2.d0*r(l,:)**2*Klambda(l,:)/3.d0
     end do
  end if

  do l=lmin,lmax
     K2(l,:) = 1.5d0*KTA(l,:)**2
  end do


! **************************
! ***   SOLVE EQUATION   ***
! **************************

! Array size.

  Naux = (Nrmax + ghost)

! Fill in (local) coefficients of linear equation.

  do l=lmin,lmax
     CL0(l,:) = - 0.5d0*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) &
              + 2.d0*D1_phi(l,:) + 2.d0/r(l,:)
     CL1(l,:) = - A(l,:)*psi4(l,:)*K2(l,:)
  end do

  if (mattertype /= "vacuum") then
     if (contains(mattertype,"nonmin")) then
        do l=lmin,lmax
           CL1(l,:) = CL1(l,:) - 0.5d0/nonmin_f(l,:)*A(l,:)*psi4(l,:)*(rho(l,:) + trS(l,:))
        end do
     else
        do l=lmin,lmax
           CL1(l,:) = CL1(l,:) - 4.d0*smallpi*A(l,:)*psi4(l,:)*(rho(l,:) + trS(l,:))
        end do
     end if
  end if

! Processor zero.

  if (rank==0) then

     if (size==1) then
        imax = Nrl(0)
     else
        imax = Nrl(0) - ghost
     end if

     do i=1-ghost,imax
        do l=lmin,lmax
           C0(l,i) = CL0(l,i)
           C1(l,i) = CL1(l,i)
        end do
     end do

!    Iterate over other processors.

     i0 = imax

     do p=1,size-1

!       Receive data from other processors.

        do l=lmin,lmax
           call MPI_RECV(CL0(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(CL1(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        end do

        if (p==size-1) then
           imax = Nrl(p)
        else
           imax = Nrl(p) - ghost
        end if

        do l=lmin,lmax
           do i=1,imax
              C0(l,i+i0) = CL0(l,i)
              C1(l,i+i0) = CL1(l,i)
           end do
        end do

        i0 = i0 + imax

     end do

! Other processors send data to processor 0.

  else

     do l=lmin,lmax
        call MPI_SEND(CL0(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(CL1(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
     end do

  end if

! Call ODE solver (only processor 0).
! Notice that C2 is 0 for maximal slicing.

  if (rank==0) then

!    Set C2=0.

     C2 = 0.d0

!    Boundary value. At the moment only allow Robin
!    and Dirichlet.

     if ((maxbound=="robin").or.(maxbound=="dirichlet")) then

        bound = maxbound
        aux = alpha0

     else if (maxbound=="conformal") then

!       In this case the outer boundary value is fixed
!       to:  alpha = alpha0 + 2r dphi/dr.
!       This is the correct value for conformally flat data.

        bound = 'dirichlet'
        aux = alpha0 + 2.d0*r(0,Nr)*D1_phi(0,Nr)

     end if

!    Call matrix inversion. The solution is saved in "u".

     call invertmatrix(lmin,lmax,aux,u,C0,C1,C2,+1,bound)

  end if


! ************************************************
! ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
! ************************************************

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! array with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  if (size==1) then
     do l=lmin,lmax
        alpha(l,:) = u(l,:)
     end do
  else
     call distribute(lmin,lmax,alpha,u)
  end if


! *****************************
! ***   LAPSE DERIVATIVES   ***
! *****************************

! Now we need to calculate the derivatives
! of the lapse.

  do l=lmin,lmax

!    Derivatives of lapse.

     diffvar => alpha
     D1_alpha(l,:) = diff1(l,+1)
     D2_alpha(l,:) = diff2(l,+1)

!                                   4
!    DD_alphar = d ( d alpha / r psi )
!                 r   r

     auxarray(l,:) = D1_alpha(l,:)/r(l,:)/psi4(l,:)

     diffvar => auxarray
     DD_alphar(l,:) = diff1(l,+1)

!    Second covariant derivative of lapse.

     Dcov2_alpha(l,:) = 1.d0/(A(l,:)*psi4(l,:))*(D2_alpha(l,:) &
                      - D1_alpha(l,:)*(0.5d0*D1_A(l,:)/A(l,:) + 2.d0*D1_phi(l,:)))

!    Laplacian of lapse.

     Lapla_alpha(l,:) = 1.d0/(A(l,:)*psi4(l,:))*(D2_alpha(l,:) &
                      - D1_alpha(l,:)*(0.5d0*D1_A(l,:)/A(l,:) - D1_B(l,:)/B(l,:) &
                      - 2.d0*D1_phi(l,:) - 2.d0/r(l,:)))

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine alphamaximal

