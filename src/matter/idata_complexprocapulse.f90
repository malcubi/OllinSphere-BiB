
  subroutine idata_complexprocapulse

! ********************************************
! ***   COMPLEX PROCA PULSE INITIAL DATA   ***
! ********************************************

! This initial data leaves the lapse, shift, conformal
! spatial metric and extrinsic curvature as in Minkowski
! (which have already been set up in "initial").
!
! For this initial data we need to solve first the Gauss
! constraint (real and imaginary parts), which for a Proca
! field has the form:
!
!    i     2
! D E  +  m  Phi  =  0
!  i
!
! with E^i the electric Proca field and Phi the scalar
! potential.

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
  real(8) zero,one,two

  real(8), dimension (0:Nl-1,1-ghost:Nrmax)   :: CL0,CL1,CL2,CLA,AL
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: C0,C1,C2,CA
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: AG,VG_R,VG_I,EG_R,EG_I
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: PhiG_R,PhiG_I
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: PsiG

! VG       Global auxiliary potentials (real and imaginary).
! EG       Global electric fields (real and imaginary).
! PhiG     Global scalar potential (real and imaginary).

! PsiG     Global conformal factor.
! AL,AG    Local and global radial metric coefficient.

! CL0,C0   Local and global coefficient of first derivative in ODE.
! CL1,C1   Local and global coefficient of linear term in ODE.
! CL2,C2   Local and global source term in ODE.
! CLA,CA   Auxiliary local and global array.


  print *, 'Complex Proca pulse initial data not yet implemented!'
  call die


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0 
  one  = 1.d0
  two  = 2.d0


! **************************************************************
! ***   SET UP INITIAL PULSE IN CONFORMAL SCALAR POTENTIAL   ***
! **************************************************************

! Message to screen.

  if (rank==0) then
     print *, 'Solving initial data for a Proca pulse ...'
  end if

! Initial pulse.  The scalar potential must be even.

  if (procaprofile=="gaussian") then

     if (proca_PhiR_r0==0.d0) then
        cprocaPhi_R = gaussian(proca_PhiR_a0,0.d0,proca_PhiR_s0)
     else
        cprocaPhi_R = gaussian(proca_PhiR_a0,+proca_PhiR_r0,proca_PhiR_s0) &
                    + gaussian(proca_PhiR_a0,-proca_PhiR_r0,proca_PhiR_s0)
     end if

     if (proca_PhiI_r0==0.d0) then
        cprocaPhi_I = gaussian(proca_PhiI_a0,0.d0,proca_PhiI_s0)
     else
        cprocaPhi_I = gaussian(proca_PhiI_a0,+proca_PhiI_r0,proca_PhiI_s0) &
                    + gaussian(proca_PhiI_a0,-proca_PhiI_r0,proca_PhiI_s0)
     end if

! Smooth top-hat profile.

  else if (procaprofile=="tophat") then

     if (proca_PhiR_r0==0.d0) then
        cprocaPhi_R = tophat(proca_PhiR_a0,0.d0,proca_PhiR_s0,proca_PhiR_t0)
     else
        cprocaPhi_R = tophat(proca_PhiR_a0,+proca_PhiR_r0,proca_PhiR_s0,proca_PhiR_t0) &
                    + tophat(proca_PhiR_a0,-proca_PhiR_r0,proca_PhiR_s0,proca_PhiR_t0)
     end if

     if (proca_PhiI_r0==0.d0) then
        cprocaPhi_I = tophat(proca_PhiI_a0,0.d0,proca_PhiI_s0,proca_PhiI_t0)
     else
        cprocaPhi_I = tophat(proca_PhiI_a0,+proca_PhiI_r0,proca_PhiI_s0,proca_PhiI_t0) &
                    + tophat(proca_PhiI_a0,-proca_PhiI_r0,proca_PhiI_s0,proca_PhiI_t0)
     end if

  end if

! Remember that the vector potential is zero.

  cprocaA_R = 0.d0
  cprocaA_I = 0.d0


! **************************************
! ***   FILL IN GLOBAL ARRAYS PhiG   ***
! **************************************

! Single processor run.

  if (size==1) then

     phiG_R = cprocaPhi_R
     phiG_I = cprocaPhi_I

! Parallel run.

  else

     print *, 'Complex Proca pulse initial data for pareller runs not yet implemented!'
     call die

     Naux = (Nrmax + ghost)

     if (rank==0) then

     else

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

     end do

! Other processors send local coefficients to processor 0.

  else

     print *, 'Complex Proca pulse initial data for pareller runs not yet implemented!'
     call die

  end if


! **********************************
! ***   SOLVE GAUSS CONSTRAINT   ***
! **********************************

! The Gauss constraints take the form:
!
!  2                     2  ^
! d  V  + (2/r) d V  =  m  Phi
!  r             r
!
! for both the real and imaginary parts.
! We solve this by direct matrix inversion.

  if (rank==0) then
     print *, 'Solving Gauss constraints ...'
  end if

! Array size.

  Naux = (Nrmax + ghost)

! Fill in (local) coefficients of linear equation.  We take
! into account the possibility that A and B are non-trivial.

  CL0 = two/r - 0.5d0*D1_A/A + D1_B/B
  CL1 = 0.d0

! Real part.

  CL2 = A*cproca_mass**2*cprocaPhi_R

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

     end do

!    Call ODE solver. Notice that in this case C1 is always 0.

     lmin = 0; lmax = Nl-1

     call invertmatrix(lmin,lmax,0.d0,VG_R,C0,C1,C2,+1,"robin")

  else

     print *, 'Complex Proca pulse initial data for pareller runs not yet implemented!'
     call die

  end if

! Imaginary part.

  CL2 = A*cproca_mass**2*cprocaPhi_I

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

     end do

!    Call ODE solver. Notice that in this case C1 is always 0.

     lmin = 0; lmax = Nl-1

     call invertmatrix(lmin,lmax,0.d0,VG_I,C0,C1,C2,+1,"robin")

  else

     print *, 'Complex Proca pulse initial data for pareller runs not yet implemented!'
     call die

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
           EG_R(l,i) = - 0.5d0*idr*(VG_R(l,i+1) - VG_R(l,i-1))/AG(l,i)
           EG_I(l,i) = - 0.5d0*idr*(VG_I(l,i+1) - VG_I(l,i-1))/AG(l,i)
        end do

!       Boundary with one-sided differences.

        EG_R(l,Nrtotal) = - 0.5d0*idr*(3.d0*VG_R(l,Nrtotal) &
           - 4.d0*VG_R(l,Nrtotal-1) + VG_R(l,Nrtotal-2))/AG(l,Nrtotal)

        EG_I(l,Nrtotal) = - 0.5d0*idr*(3.d0*VG_I(l,Nrtotal) &
           - 4.d0*VG_I(l,Nrtotal-1) + VG_I(l,Nrtotal-2))/AG(l,Nrtotal)

     else

!       Interior points.

        do i=1,Nrtotal-2
           EG_R(l,i) = - idr*(8.d0*(VG_R(l,i+1) - VG_R(l,i-1)) &
                     - (VG_R(l,i+2) - VG_R(l,i-2)))/12.d0/AG(l,i)
           EG_I(l,i) = - idr*(8.d0*(VG_I(l,i+1) - VG_I(l,i-1)) &
                     - (VG_I(l,i+2) - VG_I(l,i-2)))/12.d0/AG(l,i)
        end do

!       Boundary. We use sixth-order left-sided differences
!       to improve the boundary behavior.

        EG_R(l,Nrtotal-1) = - idr*(1.d0/6.d0*VG_R(l,Nrtotal) + 77.d0/60.d0*VG_R(l,Nrtotal-1) &
           - 5.d0/2.d0*VG_R(l,Nrtotal-2) + 5.d0/3.d0*VG_R(l,Nrtotal-3) - 5.d0/6.d0*VG_R(l,Nrtotal-4) &
           + 0.25d0*VG_R(l,Nrtotal-5) - 1.d0/30.d0*VG_R(l,Nrtotal-6))/AG(l,Nrtotal-1)

        EG_I(l,Nrtotal-1) = - idr*(1.d0/6.d0*VG_I(l,Nrtotal) + 77.d0/60.d0*VG_I(l,Nrtotal-1) &
           - 5.d0/2.d0*VG_I(l,Nrtotal-2) + 5.d0/3.d0*VG_I(l,Nrtotal-3) - 5.d0/6.d0*VG_I(l,Nrtotal-4) &
           + 0.25d0*VG_I(l,Nrtotal-5) - 1.d0/30.d0*VG_I(l,Nrtotal-6))/AG(l,Nrtotal-1)

        EG_R(l,Nrtotal  ) = - idr*(49.d0/20.d0*VG_R(l,Nrtotal) - 6.d0*VG_R(l,Nrtotal-1) &
           + 7.5d0*VG_R(l,Nrtotal-2) - 20.d0/3.d0*VG_R(l,Nrtotal-3) + 15.d0/4.d0*VG_R(l,Nrtotal-4) &
           - 6.d0/5.d0*VG_R(l,Nrtotal-5) + 1.d0/6.d0*VG_R(l,Nrtotal-6))/AG(l,Nrtotal)

        EG_I(l,Nrtotal  ) = - idr*(49.d0/20.d0*VG_I(l,Nrtotal) - 6.d0*VG_I(l,Nrtotal-1) &
           + 7.5d0*VG_I(l,Nrtotal-2) - 20.d0/3.d0*VG_I(l,Nrtotal-3) + 15.d0/4.d0*VG_I(l,Nrtotal-4) &
           - 6.d0/5.d0*VG_I(l,Nrtotal-5) + 1.d0/6.d0*VG_I(l,Nrtotal-6))/AG(l,Nrtotal)

     end if

!    Ghost points at origin.

     do i=1,ghost
        EG_R(l,1-i) = - EG_R(l,i)
        EG_I(l,1-i) = - EG_I(l,i)
     end do

  end do


! ****************************************
! ***   SOLVE HAMILTONIAN CONSTRAINT   ***
! ****************************************


! ************************************************
! ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
! ************************************************

  if (size==1) then

     cprocaV_R = VG_R
     cprocaV_I = VG_I

     cprocaE_R = EG_R
     cprocaE_I = EG_I

     psi = PsiG

  else

     print *, 'Complex Proca pulse initial data for pareller runs not yet implemented!'
     call die

  end if

! Now rescale electric field and scalar potential.

  cprocaE_R   = cprocaE_R/psi**6
  cprocaE_I   = cprocaE_I/psi**6

  cprocaPhi_R = cprocaPhi_R/psi**6
  cprocaPhi_I = cprocaPhi_I/psi**6


! ***************************************
! ***   FIND CONFORMAL FUNCTION phi   ***
! ***************************************

! Find derivative of psi.

  diffvar => psi

  do l=0,Nl-1
     D1_psi(l,:) = diff1(l,+1)
  end do

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

  end subroutine idata_complexprocapulse
