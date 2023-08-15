!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_nonminpulse.f90,v 1.22 2021/11/29 19:14:02 malcubi Exp $

  subroutine idata_nonminpulse

! *************************************
! ***   NONMIN PULSE INITIAL DATA   ***
! *************************************

! This initial data leaves the lapse, shift, conformal
! spatial metric and extrinsic curvature as in Minkowski
! (which have already been set up in "initial").
!
! We then give a simple initial profile in the nonminimally coupled
! scalar field and solve the Hamiltonian constraint for the conformal
! factor.
!
! Notice that the initial pulse is assumed to be time-symmetric,
! (i.e. pi=0), as otherwise we would have a non-zero energy
! flux and we would need to solve the coupled momentum and
! hamiltonian constraints.
!
! The Hamiltonian constraint that we need to solve has the form:
!
!  __2               5
!  \/ psi  +  psi rho / 4f =  0
!
!
! with rho the energy density of the non-minimally coupled scalar field:
!
!           2                          4
! rho  =  xi  ( 1 + 2 fpp ) / ( 2 A psi )  +  V
!
!                                                            4
!      +  fp ( d xi  +  2 xi d psi / psi  +  2 xi / r ) / psi
!               r             r
!
! (assuming that the initial pulse is time-symmetric), and
! where the Laplacian of psi given in general by:
!
! __2           2
! \/ psi  =  [ d psi  +  (2/r - d A/2A + d B/B) d psi ] / A
!               r                r        r      r
!
! In most cases we will have A=B=1, so the Laplacian simplifies,
! but this won't be the case for a transformed radial coordinate.
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

  integer i,l,p
  integer lmin,lmax
  integer i0,imax,Naux
  integer status(MPI_STATUS_SIZE)

  real(8) one,two,smallpi

  real(8), dimension (0:Nl-1,1-ghost:Nrmax)   :: CL0,CL1,CL2
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: C0,C1,C2
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: u

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


! ************************************************
! ***   SET UP INITIAL PULSE IN NONMIN FIELD   ***
! ************************************************

! Initial pulse.  Remember that the nonmin field must be even.

! Gaussian profile.

  if (nonminprofile=="gaussian") then

     if (nonmin_r0==0.d0) then
        nonmin_phi = gaussian(nonmin_a0,0.d0,nonmin_s0)
     else
        nonmin_phi = gaussian(nonmin_a0,+nonmin_r0,nonmin_s0) &
                   + gaussian(nonmin_a0,-nonmin_r0,nonmin_s0)
     end if

! r2 Gaussian  profile.

  else if (nonminprofile=="r2gaussian") then

     if (nonmin_r0==0.d0) then
        nonmin_phi = r2gaussian(nonmin_a0,0.d0,nonmin_s0)
     else
        nonmin_phi = r2gaussian(nonmin_a0,+nonmin_r0,nonmin_s0) &
                   + r2gaussian(nonmin_a0,-nonmin_r0,nonmin_s0)
     end if

! Smooth top-hat profile.

  else if (nonminprofile=="tophat") then

     if (nonmin_r0==0.d0) then
        nonmin_phi = tophat(nonmin_a0,0.d0,nonmin_s0,nonmin_t0)
     else
        nonmin_phi = tophat(nonmin_a0,+nonmin_r0,nonmin_s0,nonmin_t0) &
                   + tophat(nonmin_a0,-nonmin_r0,nonmin_s0,nonmin_t0)
     end if

  end if

! The spatial derivative is calculated with finite differences.

  diffvar => nonmin_phi
 
  do l=0,Nl-1
     nonmin_xi(l,:) = diff1(l,+1)
  end do

! We also need the derivative of xi.

  diffvar => nonmin_xi

  do l=0,Nl-1
     D1_nonmin_xi(l,:) = diff1(l,-1)
  end do

! The initial time derivative is set to 0.

  nonmin_pi = 0.d0

! Non-minimal coupling.

  call nonmincoupling


! **************************
! ***   ZERO POTENTIAL   ***
! **************************

! In the case of zero potential the Hamiltonian
! constraint is in fact linear in psi, so we can
! just solve it directly by matrix inversion.
!
! The equation takes the final form:
!
!  2
! d psi  +  C0 d psi  +  C1 psi  =  0
!  r            r
!
! with:
!
! C0  =  2/r  +  (1/2f) fp xi
!
!                   2
! C1  =  (1/8f) [ xi ( 1 + 2 fpp)  +  2 fp ( d xi + 2 xi / r) ]
!                                             r

  if (nonminpotential=="none") then

!    Message to screen.

     if (rank==0) then
        print *, 'Solving initial data for a nonmin pulse with zero potential ...'
        print *
     end if

!    Array size.

     Naux = (Nrmax + ghost)

!    Fill in (local) coefficients of linear equation.
!    Remember that C0 is the coefficient of dpsi/dr,
!    C1 is the coefficient of psi, and C2 is the
!    right hand side.

     if (newr) then
        CL0 = two/r + 0.5d0/nonmin_f*nonmin_fp*nonmin_xi - 0.5d0*D1_A/A + D1_B/B
        CL1 = 0.25d0/r**2*(one - A/B) + 0.25/r*(3.d0*D1_B/B - D1_A/A + 2.d0*nonmin_fp/nonmin_f*nonmin_xi)  &
            - 0.125d0*D1_A/A*(D1_B/B + nonmin_fp*nonmin_xi/nonmin_f) - 0.0625*(D1_B**2 - 4.d0*B*D2_B)/B**2 &
            + one/nonmin_f*(0.125*(one + 2.d0*nonmin_fpp)*nonmin_xi**2 + 0.25*nonmin_fp*(D1_B*nonmin_xi/B + D1_nonmin_xi))
     else
        CL0 = two/r + 0.5d0/nonmin_f*nonmin_fp*nonmin_xi
        CL1 = (0.125d0/nonmin_f)*(nonmin_xi**2*(1.d0 + 2.d0*nonmin_fpp) &
            + 2.d0*nonmin_fp*(D1_nonmin_xi + 2.d0*nonmin_xi/r))
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

        if (newr) then
           call invertmatrix(lmin,lmax,1.d0,u,C0,C1,C2,+1,"robin")
        else
           call invertmatrix(lmin,lmax,1.d0,u,C0,C1,C2,+1,"robin")
        end if

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
! is no longer linear.

  else

!    At the moment this is not yet implemented.

     if (rank==0) then
        print *
        print *, 'Initial data with a non-zero potential for the nonmin field'
        print *, 'is not yet implemented.  Choose nonminpotential="none".'
        print *, '(subroutine idata_nonminpulse.f90)'
        print *
     end if

     call die

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
     D1_chi = - dble(chipower)*D1_psi/psi**chipower
  end if

! Find psi2 and psi4.

  psi2 = psi**2
  psi4 = psi**4


! ***************
! ***   END   ***
! ***************

end subroutine idata_nonminpulse

