!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/sources_lapse.f90,v 1.35 2023/03/13 17:47:51 malcubi Exp $

  subroutine sources_lapse(l,dtw)

! *****************************
! ***   SOURCES FOR LAPSE   ***
! *****************************

! This routine calculates the sources for the lapse.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use derivatives

! Extra variables.

  implicit none

  logical contains
  logical evolvelapse

  integer i,l,p
  integer i0,imax,Naux
  integer status(MPI_STATUS_SIZE)

  real(8) r0,dtw
  real(8) interp
  real(8) aux0,aux1,aux2

  real(8), dimension (0:Nl-1,1-ghost:Nrmax)   :: CL0,CL1,CL2
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: u,C0,C1,C2

  character bound*10


! *****************************
! ***   STATIC OR MAXIMAL   ***
! *****************************

! Sources is set to zero for static or maximal slicing.

  if ((slicing=="static").or.(slicing=="maximal")) then

!    Set flag for lapse evolution to false.

     evolvelapse = .false.

!    Set lapse source to 0.

     salpha(l,:) = 0.d0

!    Cosmological background lapse.

     if (cosmic_run) then
        scosmobg_alpha(l) = 0.d0
     end if


! **********************************
! ***   FIRST ORDER BONA-MASSO   ***
! **********************************

! First order Bona-Masso type slicings:
!
! d alpha  = - alpha^2 f(alpha) trK  +  beta d alpha
!  t                                          r

  else if ((slicing=="harmonic").or.(slicing=="1+log").or. &
           (slicing=="shockavoid").or.(slicing=="alphaminus2")) then

!    Set flag for lapse evolution to true.

     evolvelapse = .true.

!    Source for cosmological background lapse.

     if (cosmic_run) then
       scosmobg_alpha(l) = - cosmobg_falpha(l)*cosmobg_trK(l)
     end if

!    First the term with no shift.

     salpha(l,:) = - falpha(l,:)*trK(l,:)

!    Shift term.

     if (shift/="none") then
        salpha(l,:) = salpha(l,:) + beta(l,:)*DA_alpha(l,:)
     end if

!    Add extra term for the case of a non-minimally coupled
!    scalar field.

     if (contains(mattertype,"nonmin")) then
        salpha(l,:) = salpha(l,:) + alpha(l,:)**2*nonmin_theta*(nonmin_fp(l,:)/nonmin_f(l,:))*nonmin_pi(l,:)
     end if


! ***********************************
! ***   COSMOLOGICAL BONA-MASSO   ***
! ***********************************

! For cosmological runs we can use a synchronous background lapse,
! or a conformal background lapse, and evolve the "dynamical" part
! using some other Bona-Masso type condition.

  else if (index(adjustl(slicing),"cosmo")/=0) then

     if (cosmic_run) then

!       Source for cosmological background lapse.

        scosmobg_alpha(l) = - cosmobg_falpha(l)*cosmobg_trK(l)

!       Source for dynamical part.

        salpha(l,:) = scosmobg_alpha(l) - falpha(l,:)*(trK(l,:) - cosmobg_trK(l))

!       Shift term.

        if (shift/="none") then
           salpha(l,:) = salpha(l,:) + beta(l,:)*DA_alpha(l,:)
        end if

     else

        if (rank==0) then
           print *
           print *, 'Slicings of type "cosmo" must have cosmic_run=.true.'
           print *, 'Aborting! (subroutine sources_lapse.f90)'
           print *
        end if

        call die

     end if


! ************************
! ***   SANITY CHECK   ***
! ************************

  else

     if (rank==0) then
        print *
        print *, 'Unknown slicing condition.'
        print *, 'Aborting! (subroutine sources_lapse.f90)'
        print *
     end if

     call die

  end if


! ****************************
! ***   REAL DISSIPATION   ***
! ****************************

! Here we add a real (i.e. not numerical) dissipation term
! to the Bona-Masso slicings, which means that we now have
! a parabolic slicing condition instead of hyperbolic.
! Because of this we need to use an implicit numerical
! method for updating the lapse.
!
! The evolution equation now takes the form:
!
!                                     __2
!    d alpha  =  salpha  +  lapsediss \/ alpha
!     t
!
! and the equation we are in fact inverting becomes:
!
!   __2     n+1         n+1
!   \/ alpha    -  alpha   / (lapsediss dt)  =
!
!               n             n
!       - (alpha  +  dt salpha ) / (lapsediss dt)
!
! where "lapsediss" is the dissipation coefficient, "salpha"
! is the Bona-Masso source calculated above, and where the
! super-indices n and n+1 indicate current and next time step.
!
! Remember that the Laplacian of alpha is given in general by:
!
! __2            2                                                            4
! \/ alpha  = [ d alpha  -  d alpha (d A/2A - d B/B - 2 d phi - 2/r) ] / A psi
!                r           r        r        r         r
!
! Notice that we have calculated the dissipation term on the
! top time level for stability, i.e. we use an implicit scheme,
! which requires a matrix inversion.  The whole thing is only
! first order accurate, but this should be OK since the higher
! accuracy is later obtained in the Runge-Kutta iterations.

  if ((lapsediss/=0.d0).and.evolvelapse) then

!    Array size.

     Naux = (Nrmax + ghost)

!    Auxiliary variable, essentially the inverse of
!    the dissipation parameter.  The A*psi4 comes from
!    the coefficient of D2_alpha in the laplacian.

     auxarray(l,:) = A(l,:)*psi4(l,:)/(lapsediss*dtw)

!    Here we make sure that we kill the dissipation close 
!    to the physical outer boundary (at 80% of the distance).

     auxarray(l,:) = 2.d0*auxarray(l,:)/(1.d0 - tanh(r(l,:)-0.8d0*rbound))

!    Fill in (local) coefficients of linear equation.  I have tested
!    slightly different ways to solve the equation.  The difference
!    seems to be quite small.

     CL0(l,:) = - 0.5d0*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) &
              + 2.d0*D1_phi(l,:) + 2.d0/r(l,:)
     CL1(l,:) = - auxarray(l,:)
     CL2(l,:) = - auxarray(l,:)*(alpha_p(l,:) + dtw*salpha(l,:))

     !CL0(l,:) = 0.d0
     !CL1(l,:) = - auxarray(l,:)
     !CL2(l,:) = - auxarray(l,:)*(alpha_p(l,:) + dtw*salpha(l,:)) &
     !         + (0.5d0*D1_A(l,:)/A(l,:) - D1_B(l,:)/B(l,:) &
     !         - 2.d0*D1_phi(l,:) - 2.d0/r(l,:))*D1_alpha(l,:)

     !CL0(l,:) = - 0.5d0*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:) &
     !         + 2.d0*D1_phi(l,:) + 2.d0/r(l,:)
     !CL1(l,:) = - 2.d0*auxarray(l,:)
     !CL2(l,:) = - 2.d0*auxarray(l,:)*(alpha_p(l,:) + dtw*salpha(l,:)) &
     !         + (0.5d0*D1_A(l,:)/A(l,:) - D1_B(l,:)/B(l,:) &
     !         - 2.d0*D1_phi(l,:) - 2.d0/r(l,:))*D1_alpha(l,:) &
     !         - D2_alpha(l,:)

!    Processor zero receives data and solves ODE.

     if (rank==0) then

        if (size==1) then
           imax = Nrl(0)
        else
           imax = Nrl(0) - ghost
        end if

        do i=1-ghost,imax
           C0(l,i) = CL0(l,i)
           C1(l,i) = CL1(l,i)
           C2(l,i) = CL2(l,i)
        end do

!       Iterate over other processors.

        i0 = imax

        do p=1,size-1

!          Receive data from other processors.

           call MPI_RECV(CL0(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(CL1(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(CL2(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

           if (p==size-1) then
              imax = Nrl(p)
           else
              imax = Nrl(p) - ghost
           end if

           do i=1,imax
              C0(l,i+i0) = CL0(l,i)
              C1(l,i+i0) = CL1(l,i)
              C2(l,i+i0) = CL2(l,i)
           end do

           i0 = i0 + imax

        end do

!    Other processors send data to processor 0.

     else

        call MPI_SEND(CL0(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(CL1(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(CL2(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

     end if

!    Boundary. We use a Dirichlet boundary condition with
!    the boundary value given by what it would have been without
!    dissipation (notice that the dissipation is switched
!    off close to the boundary).  But we need to broadcast
!    the value at the physical boundary to other processors
!    (in particular processor 0).

!    Coarse grid: physical boundary.

     if (l==0) then

        aux0 = alpha_p(0,Nr) + dtw*salpha(0,Nr)

!    Fine grids: interpolate from coarser grid. This code
!    is bassically copied from the one in "onestep.f90".

     else

        r0 = (dble(Nrtotal)-0.5d0)*dr(l)

        interpvar => alpha
        aux1 = interp(l-1,r0,.false.)

        call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

        aux0 = dtw/(t(l-1)-t(l))*aux2 + (t(l)+dtw-t(l-1))/(t(l)-t(l-1))*alpha_p(l,Nr)

     end if

     call MPI_BCAST(aux0,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)

!    Call ODE solver (only processor 0).

     if (rank==0) then

!       Boundary condition.

        bound = 'dirichlet'

!       Call matrix inversion. The solution is saved in "u".

        call invertmatrix(l,l,aux0,u,C0,C1,C2,+1,bound)

     end if

!    Distribute solution among processors. The distributed
!    solution is saved in auxarray.

     if (size==1) then
        auxarray(l,:) = u(l,:)
     else
        call distribute(l,l,auxarray,u)
     end if

!    Now set the source term to the difference between the solution
!    above and the old value of the lapse.  The idea is not to
!    change the value of "alpha" until the update step in order
!    to avoid contaminating other sources that still need to be
!    calculated.

     salpha(l,:) = (auxarray(l,:) - alpha_p(l,:))/dtw

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine sources_lapse
