!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/idata_complexDMmom.f90,v 1.27 2022/07/12 23:59:18 malcubi Exp $

  subroutine idata_complexDMmom

! ************************************************
! ***   COMPLEX SCALAR DARK MATTER COSMOLOGY   ***
! ************************************************

! This routine first sets up a dark matter complex scalar field
! cosmological background.  It can also add a perturbation.
!
! The difference between this routine and the one in found in
! idata_complexDM.f90 is that in that case we assumed initial
! data such that the momentum constraint is trivially satisfied.
! Here, however, we want to set up a perturbation such that
! the initial momentum flux is non-zero, so we need to solve
! the coupled hamiltonian and momentum constraints.
!
! We will assume conformally flat initial data.  We will also assume
! that trK only has a background constant value. In that case the
! Hamiltonian and momentum constraints take the form:
!
! __ 2            5           2          2
! \/ psi  =  - psi  [ (3/2 KTA  - 2/3 trK ) / 8 + 2 pi rho ]
!
!
!
! d KTA  =  - 3 KTA [ 2 (d psi / psi)  +  d B / 2B  +  1/r ]  +  8 pi JA
!  r                      r                r
!
! with rho and JA the energy and momentum density of the complex
! scalar field:
!
!             2      2        2     2         4
! rho =  [ PiR  + PiI  +  (XiR  + XiI )/ A psi  ] / 2  +  V(phi)
!
!
! JA  =  - PiR XiR - PiI XiI
!
! In the Hamiltonian constraint the Laplacian of psi is given by:
!
! __2           2
! \/ psi  =  [ d psi  +  (2/r - d A/2A + d B/B) d psi ] / A
!               r                r        r      r
!
! In most cases we will have A=B=1, so the Laplacian simplifies.
! But this won't be the case for a transformed radial coordinate
! (notice that we still assume that the conformal Ricci vanishes,
! otherwise we will need to add the corresponding terms).
!
! For the boundary conditions at the origin I take:
!
! psi(r=0)  =  psi0   d psi(r=0)  =  0      (psi ~ psi0 + a r^2)
!                      r
!
! KTA(r=0)  =  0      d KTA(r=0)  =  0      (KTA ~ b r^2)
!                      r
!
! Two important comments should be made about these conditions:
!
! I)  At the origin KTA goes as r^2, which means that a second
!     order integration won't be able to capture this correctly
!     (since close to the origin the error is of the same order
!     as KTA). In fact, I have found that even with fourth order
!     Runge-Kutta, Klambda (essentially KTA/r^2) does not converge
!     at the origin.  This is probably because asking for KTA~r^2
!     requires to have both KTA and its first derivative equal
!     to 0 at r=0, and this is overdetermined since the momentum
!     constraint is a first order differential equation.
!
!     But fortunately, if we rewrite the momentum constraint in
!     terms of Khat := KTA/r, then this problem seems to go away:
!
!     d Khat  =  - 2 Khat [ 2/r  +  3 d psi / psi  + (3/4) d B/B ]
!      r                               r                    r
!
!             +  8 pi JA / r
!
!
!     But we only get second order convergence at the origin even
!     with the fourth order scheme (I am not sure why).
!
! II) I first solve the coupled constraints taking psi0=1. But
!     we don't really want psi0=1, but rather psi=1 at infinity.
!     Unfortunately, the equation for psi is non-linear in the source
!     term, so we can't just multiply psi with a constant to fix this.
!     However, once we have found a solution, we can rescale it as:
!
!     r -> k r ,       psi -> psi/sqrt(k)
!
!     and we will still have a solution (we can check that this
!     is true for all equations above, since the density "rho"
!     for a scalar field is invariant with this change).
!
!     But in order for this rescaling to work we must rescale the
!     initial profile width multipling all length scales with k.
!
!     This means that the initial profle will no longer correspond
!     to what was set up in the parameter file (in some cases
!     this change can be quite large).
!
!
! NOTE FOR PARALLEL RUNS:  The initial data is not really
! solved in parallel.  It is in fact solved only on processor
! zero on a full size array, and then it is distributed
! among the other processors.  This is slow, but works.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical minimize                       ! Find solution with original parameters.

  integer l                              ! Counter

  real(8) H0                             ! Hubble parameter.
  real(8) R_r0,R_s0,R_t0,R_r1,R_s1,R_t1  ! Initial real profile parameters.
  real(8) I_r0,I_s0,I_t0,I_r1,I_s1,I_t1  ! Initial imaginary profile parameters.
  real(8) ax,bx,cx,fa,fb,fc,xmin,fmin    ! Parameters for minimization.
  real(8) zero,half,one,two,third,smallpi
  real(8) aux

  common / pert / H0,minimize
  common / initialprofile / R_r0,R_s0,R_t0,R_r1,R_s1,R_t1,I_r0,I_s0,I_t0,I_r1,I_s1,I_t1


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  half  = 0.5d0
  one   = 1.d0
  two   = 2.d0
  third = 1.d0/3.d0

  smallpi = acos(-one)


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  if (rank==0) then
     print *, 'Solving initial data for complex scalar field dark matter cosmology'
     print *, 'with non-trivial momentum constraint ...'
     print *
  end if


! ************************
! ***   SANITY CHECK   ***
! ************************

! Since we are assuming a harmonic time dependence for the
! background, the mass can't be zero.

  if (complex_mass==0.d0) then
     if (rank==0) then
        print *, 'For a cosmological perturbation of a complex scalar field'
        print *, 'the parameter complex_mass must be non-zero, as otherwise'
        print *, 'the background has zero energy density.'
        print *, 'Aborting!  (subroutine idata_complexDMmom)'
        print *
     end if
     call die
  end if

! For initial data corresponding to a growing mode the imaginary
! part of the background scalar field must be zero. Also, the
! real part of the perturbation must be zero (this is explained
! below).

  if (complexDM_type=='growing') then

     if (complex_bg_phiI0/=0.d0) then
        if (rank==0) then
           print *, 'For a cosmological perturbation of a complex with complexDM_type=growing'
           print *, 'the imaginary part of the background scalar field must be zero.'
           print *, 'Aborting!  (subroutine idata_complexDMmom)'
           print *
        end if
        call die
     end if

     if (complexR_a0/=0.d0) then
        if (rank==0) then
           print *, 'For a cosmological perturbation of a complex with complexDM_type=growing'
           print *, 'the real part of the pertrubartion must be zero.'
           print *, 'Aborting!  (subroutine idata_complexDMmom)'
           print *
        end if
        call die
     end if

  end if


! ***********************************
! ***   COSMOLOGICAL BACKGROUND   ***
! ***********************************

! For the background we set up a uniform field.

  complex_phiR = complex_bg_phiR0
  complex_phiI = complex_bg_phiI0

  complex_xiR = zero
  complex_xiI = zero

! For the background we also assume a harmonic time
! dependence, with a frequency equal to the mass.

  complex_piR = - complex_mass*complex_bg_phiI0
  complex_piI = + complex_mass*complex_bg_phiR0

! Find background potential.

  do l=0,Nl-1
     call potential(l)
  end do

! The Hamiltonian constraint has the form:
!                                                  2      2
! trK  =  - sqrt(24 pi rho)  =  - sqrt[ 24 pi ( (Pi_r + Pi_i)/2 + V ) ]
!
! with V the complex field potential. The minus sign
! corresponds to an expanding Universe.

  trK = - sqrt(24.d0*smallpi*(half*(complex_piR**2 + complex_piI**2) + complex_V))

! Cosmological background quantities.

  if (cosmic_run) then

!    Geometry.

     cosmobg_trK = trK(0,Nr)

     H0 = - trK(0,Nr)/3.d0
     cosmobg_H = H0

!    Complex field.

     cosmobg_complex_phiR = complex_bg_phiR0
     cosmobg_complex_phiI = complex_bg_phiI0

!    Time derivative.

     cosmobg_complex_piR = - complex_mass*complex_bg_phiI0
     cosmobg_complex_piI = + complex_mass*complex_bg_phiR0

  end if


! ************************
! ***   PERTURBATION   ***
! ************************

! Save original profile parameters.

  R_r0 = complexR_r0
  R_s0 = complexR_s0
  R_t0 = complexR_t0

  I_s0 = complexI_s0
  I_s0 = complexI_s0
  I_t0 = complexI_t0

  if (complexprofile(1:4)=="comp") then

     R_r1 = complexR_r1
     R_s1 = complexR_s1
     R_t1 = complexR_t1

     I_s1 = complexI_s1
     I_s1 = complexI_s1
     I_t1 = complexI_t1

  end if

! Check if we want to rescale the initial profile and
! redefine initial profile parameters.

  if (rescaledata/=1.d0) then

     call rescaleprofile(rescaledata,.true.)

     R_r0 = complexR_r0
     R_s0 = complexR_s0
     R_t0 = complexR_t0

     I_s0 = complexI_s0
     I_s0 = complexI_s0
     I_t0 = complexI_t0

     if (complexprofile(1:4)=="comp") then

        R_r1 = complexR_r1
        R_s1 = complexR_s1
        R_t1 = complexR_t1

        I_s1 = complexI_s1
        I_s1 = complexI_s1
        I_t1 = complexI_t1

     end if

  end if

! Add a perturbation to the background scalar field.
! The logical flag "minimize" calls a minimization
! algorithm to try to find a solution with the
! original profile parameters.  Is is hard wired
! to ".true.", but can easily be changed.

  minimize = .true.
  !minimize = .false.

  if (complex_bg_pert) then

!    Solve using the given profile parameters, and accept
!    that the final solution will rescale the initial
!    profile with some factor.

     if (.not.minimize) then

        call complexDMmom_pert

!    Use a minimization algorithm to try to find a
!    solution compatible with the original parameters.

     else

!       Bracket solutiion that minimizes the distance
!       to the original parameters.  For this we call
!       the Numrrical Recipes routine "mnbrak".

        ax = 1.0d0
        bx = 1.1d0
        cx = 0.d0

        call mnbrak(ax,bx,cx,fa,fb,fc)
        !print *,ax,bx,cx,fa,fb,fc

!       Now call the minimization routine.  For this we
!       use the Numerical Recipes routine "golden".

        call golden (ax,bx,cx,xmin,fmin)
        !print *,xmin,fmin

!       Check if the mimimun is small enough.

        if (fmin>1.d-5) then
           if (rank==0) then
              print *, 'WARNING:  a solution with the original profile parameters'
              print *, 'could not be found (it probably does not exist).'
              print *, 'Aborting!  (subroutine idata_complexDMmom)'
              print *
           end if
           call die
        end if

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine idata_complexDMmom












  subroutine complexDMmom_pert

! ************************
! ***   PERTURBATION   ***
! ************************

! This routine solves the coupled constraints.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use derivatives
  use radialfunctions

! Extra variables.

  implicit none

  logical rescale                        ! Flag for rescaling conformal factor.
  logical minimize                       ! Find solution with original parameters.

  integer l,i,p                          ! Counters.
  integer imin,imax,i0,Naux              ! Local grid limits.
  integer status(MPI_STATUS_SIZE)

  real(8) H0                             ! Hubble parameter.
  real(8) psiorigin,psiasymp             ! Values of psi at origin and infinity.
  real(8) r0,rm,delta                    ! Local position and step size.
  real(8) psi0,Dpsi0,KTA0,Khat0          ! Local values of variables.
  real(8) k11,k12,k13,k14                ! Runge-Kutta data.
  real(8) k21,k22,k23,k24                ! Runge-Kutta data.
  real(8) k31,k32,k33,k34                ! Runge-Kutta data.
  real(8) psi_rk,Dpsi_rk                 ! Local values of psi and Dpsi for Runge-Kutta.
  real(8) KTA_rk,Khat_rk,trK_rk          ! Local values of KTA, Khat and trK for Runge-Kutta.
  real(8) A_rk,B_rk,DA_rk,DB_rk          ! Local values of (A,B) and derivatives for Runge-Kutta.
  real(8) pi2_rk,xi2_rk                  ! Local values of Pi^2 and Xi^2 for Runge-Kutta.
  real(8) V_rk,JA_rk,rho_rk              ! Local potential, momentum and energy density for Runge-Kutta.
  real(8) zero,half,one,two,third,smallpi
  real(8) aux

  real(8), dimension (0:Nl-1,1-ghost:Nrmax) :: phiR_l,phiI_l       ! Local scalar field.
  real(8), dimension (0:Nl-1,1-ghost:Nrmax) :: xiR_l,xiI_l         ! Local derivative of scalar field.
  real(8), dimension (0:Nl-1,1-ghost:Nrmax) :: piR_l,piI_l         ! Local time derivative of scalar field.
  real(8), dimension (0:Nl-1,1-ghost:Nrmax) :: V_l                 ! Local potential.
  real(8), dimension (0:Nl-1,1-ghost:Nrmax) :: trK_l               ! Local trK.
  real(8), dimension (0:Nl-1,1-ghost:Nrmax) :: A_l,B_l,DA_l,DB_l   ! Local (A,B) and derivatives.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr                ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: KTA_g,Khat_g      ! Global KTA and Khat.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: trK_g             ! Global trK.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g,B_g,DA_g,DB_g ! Global (A,B)
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: psi_g,Dpsi_g      ! Global conformal factor.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: phiR_g,phiI_g     ! Global scalar field.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: xiR_g,xiI_g       ! Global derivative of scalar field.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: piR_g,piI_g       ! Global time derivative of scalar field.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: V_g               ! Global potential.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: Xi2_g,Pi2_g       ! Global squares of Xi and Pi.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: JA_g              ! Global momentum density.

  common / pert / H0,minimize


! *******************
! ***   NUMBERS   ***
! *******************

  zero  = 0.d0
  half  = 0.5d0
  one   = 1.d0
  two   = 2.d0
  third = 1.d0/3.d0

  smallpi = acos(-one)


! **************************************
! ***   SOLVE FOR THE PERTURBATION   ***
! **************************************

! Initialize value of psi at origin.  We need to do this here,
! because once we rescale the value of "psiorigin" at the end
! of the routine we must recalculate the perturbation.

  rescale = .false.
  psiorigin = 1.d0

! Continue statement for when we recalculate the perturbation
! with a rescaled conformal factor and gaussian width.

  100 continue

! Gaussian profile.

  if (complexprofile=="gaussian") then

     if (complexR_r0==0.d0) then
        complex_phiR = complex_bg_phiR0 + gaussian(complexR_a0,0.d0,complexR_s0)
     else
        complex_phiR = complex_bg_phiR0 + gaussian(complexR_a0,+complexR_r0,complexR_s0) &
                                        + gaussian(complexR_a0,-complexR_r0,complexR_s0)
     end if

     if (complexI_r0==0.d0) then
        complex_phiI = complex_bg_phiI0 + gaussian(complexI_a0,0.d0,complexI_s0)
     else
        complex_phiI = complex_bg_phiI0 + gaussian(complexI_a0,+complexI_r0,complexI_s0) &
                                        + gaussian(complexI_a0,-complexI_r0,complexI_s0)
     end if

! Smooth top-hat profile.

  else if (complexprofile=="tophat") then

     if (complexR_r0==0.d0) then
        complex_phiR = complex_bg_phiR0 + tophat(complexR_a0,0.d0,complexR_s0,complexR_t0)
     else
        complex_phiR = complex_bg_phiR0 + tophat(complexR_a0,+complexR_r0,complexR_s0,complexR_t0) &
                                        + tophat(complexR_a0,-complexR_r0,complexR_s0,complexR_t0)
     end if

     if (complexI_r0==0.d0) then
        complex_phiI = complex_bg_phiI0 + tophat(complexI_a0,0.d0,complexI_s0,complexI_t0)
     else
        complex_phiI = complex_bg_phiI0 + tophat(complexI_a0,+complexI_r0,complexI_s0,complexI_t0) &
                                        + tophat(complexI_a0,-complexI_r0,complexI_s0,complexI_t0)
     end if

! Compensated gaussian profile. Here we add two gaussian functions centered at r=0.

  else if (complexprofile=="comp-gaussian") then

     complex_phiR = complex_bg_phiR0 + gaussian(complexR_a0,0.d0,complexR_s0) &
                                     + gaussian(complexR_a1,0.d0,complexR_s1)

     complex_phiI = complex_bg_phiI0 + gaussian(complexI_a0,0.d0,complexI_s0) &
                                     + gaussian(complexI_a1,0.d0,complexI_s1)

! Compensated top-hat profile. Here we add two top-hat functions centered at r=0.

  else if (complexprofile=="comp-tophat") then

     complex_phiR = complex_bg_phiR0 + tophat(complexR_a0,0.d0,complexR_s0,complexR_t0) &
                                     + tophat(complexR_a1,0.d0,complexR_s1,complexR_t1)

     complex_phiI = complex_bg_phiI0 + tophat(complexI_a0,0.d0,complexI_s0,complexI_t0) &
                                     + tophat(complexI_a1,0.d0,complexI_s1,complexI_t1)

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

! Time derivatives.

  if (complexDM_type=='harmonicMOM') then

!    In this case we assume a purely harmonic time
!    dependence for the initial time derivative.

     complex_piR = - complex_mass*complex_phiI
     complex_piI = + complex_mass*complex_phiR

  else if (complexDM_type=='growing') then

!    In this case the initial perturbation corresponds to a
!    purely growing mode (see notes by J.C Hidalgo and D. Núñez).
!    For the perturbed scalar field we have:
!
!                            1/2
!    dphi_R  =  13/15 (3/8pi)   (H0/m) F(r)  =  13/15 phi0 F(r)
!
!                           1/2
!    dphi_I  = - 4/15 (3/8pi)   F(r)  = - 4/13 (m/H0) dphi_R
!
!    where phi0 is the initial background (purely real) scalar field,
!    and F(r) is the perturbation profile (called R^(2)(r) in the notes).
!    Since we assume that H0/m is small, it turns out that the imaginary
!    part of the perturbation dominates, so it is more natural to write:
!
!    dphi_I  =  G(r)
!
!    dphi_R  =  - 13/4 (H0/m) G(r)
!
!    with G(r) the perturbation profile (gaussian). We then assume that
!    the imaginary part of perturbation was already set up above and
!    rewrite the real part. Notice also that complex_phiI only contains
!    the perturbed part (the background part is zero).

     complexR_s0  = complexI_s0
     complex_phiR = complex_bg_phiR0 - 13.d0/4.d0*(H0/complex_mass)*complex_phiI

     diffvar => complex_phiR

     do l=0,Nl-1
        complex_xiR(l,:) = diff1(l,+1)
     end do

!    For the time derivatives we use the expression from the notes:
!
!                                     1/2
!    Pi  =  i m phi  -  (m/15) (3/8pi)   (H0/m) [ (13/2) (H0/m)  + 4 i ] F(r)
!
!    Here we are ignoring a term in the notes that is not really
!    consistent.  With the above definitions this becomes:
!
!    Pi  =  i m phi  +  (H0/4) G(r) [ (13/2) (H0/m)  +  4 i ]
!
!    Notice that this includes both background and perturbation.
!    From this we find:
!
!    Pi_R  =  m G(r)  [ (13/8) (H0/m)^2 - 1 ]  =  m Phi_I [ (13/8) (H0/m)^2 - 1 ]
!
!    Pi_I  =  m [ Phi0 - 9/4 (H0/m) G(r) ]

     complex_piR = complex_mass*((13.d0/8.d0)*(H0/complex_mass)**2 - one)*complex_phiI
     complex_piI = complex_mass*(complex_bg_phiR0 - (9.d0/4.d0)*(H0/complex_mass)*complex_phiI)

  end if

! Calculate potential again after setting up perturbation.

  do l=0,Nl-1
     call potential(l)
  end do


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

  do l=0,Nl-1
     do i=1-ghost,Nrtotal
        rr(l,i) = (dble(Nmin(rank) + i) - half)*dr(l)
     end do
  end do


! ********************************************
! ***   COPY LOCAL DATA TO GLOBAL ARRAYS   ***
! ********************************************

! Array size.

  Naux = (Nrmax + ghost)

! For single processor runs just copy the data.

  if (size==1) then

!    Scalar field and derivatives.

     phiR_g = complex_phiR
     phiI_g = complex_phiI

     xiR_g = complex_xiR
     xiI_g = complex_xiI

     piR_g = complex_piR
     piI_g = complex_piI

!    Potential.

     V_g = complex_V

!    Trace of extrinsic curvature.

     trK_g = trK

!    Metric coefficients and derivatives.

     A_g = A
     B_g = B

     DA_g = D1_A
     DB_g = D1_B

! For parallel runs processor 0 must receive
! data from other processors.

  else

     if (rank==0) then

        if (size==1) then
           imax = Nrl(0)
        else
           imax = Nrl(0) - ghost
        end if

!       Copy local data from proc 0 to global arrays.

        do i=1-ghost,imax

           phiR_g(:,i) = complex_phiR(:,i)
           phiI_g(:,i) = complex_phiI(:,i)

           xiR_g(:,i) = complex_xiR(:,i)
           xiI_g(:,i) = complex_xiI(:,i)

           piR_g(:,i) = complex_piR(:,i)
           piI_g(:,i) = complex_piI(:,i)

           V_g(:,i)   = complex_V(:,i)
           trK_g(:,i) = trK(:,i)

           A_g(:,i) = A(:,i)
           B_g(:,i) = B(:,i)

           DA_g(:,i) = D1_A(:,i)
           DB_g(:,i) = D1_B(:,i)

        end do

!       Iterate over other processors to receive data.

        i0 = imax

        do p=1,size-1

!          Receive local data from other processors
!          and copy them into global arrays.

           do l=0,Nl-1

              call MPI_RECV(phiR_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(phiI_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

              call MPI_RECV(xiR_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(xiI_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

              call MPI_RECV(piR_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(piI_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

              call MPI_RECV(V_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(trK_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

              call MPI_RECV(A_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(B_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

              call MPI_RECV(DA_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(DB_l(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

           end do

           if (p==size-1) then
              imax = Nrl(p)
           else
              imax = Nrl(p) - ghost
           end if

           do i=1,imax

              phiR_g(:,i+i0) = phiR_l(:,i)
              phiI_g(:,i+i0) = phiI_l(:,i)

              xiR_g(:,i+i0) = xiR_l(:,i)
              xiI_g(:,i+i0) = xiI_l(:,i)

              piR_g(:,i+i0) = piR_l(:,i)
              piI_g(:,i+i0) = piI_l(:,i)

              V_g(:,i+i0)   = V_l(:,i)
              trK_g(:,i+i0) = trK_l(:,i)

              A_g(:,i+i0) = A_l(:,i)
              B_g(:,i+i0) = B_l(:,i)

              DA_g(:,i+i0) = DA_l(:,i)
              DB_g(:,i+i0) = DB_l(:,i)

           end do

           i0 = i0 + imax

        end do

!    Other processors send local data to processor 0.

     else

        do l=0,Nl-1

           call MPI_SEND(complex_phiR(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(complex_phiI(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

           call MPI_SEND(complex_xiR(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(complex_xiI(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

           call MPI_SEND(complex_piR(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(complex_piI(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

           call MPI_SEND(complex_V(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(trK(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

           call MPI_SEND(A(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(B(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

           call MPI_SEND(D1_A(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(D1_B(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

        end do

     end if

  end if

! Squares of Xi and Pi.

  Xi2_g = XiR_g**2 + XiI_g**2
  Pi2_g = PiR_g**2 + PiI_g**2

! Momentum density.

  JA_g = - piR_g*xiR_g - piI_g*xiI_g


! *************************************
! ***   SOLVE COUPLED CONSTRAINTS   ***
! *************************************

! We solve the coupled constraints using fourth order Runge-Kutta.
! Only processor 0 solves the initial data.

! Initialize arrays.

  psi_g  = psiorigin
  Dpsi_g = 0.d0

  KTA_g  = 0.d0
  Khat_g = 0.d0


! *********************************************
! ***   ONLY PROCESSOR 0 SOLVES THE ODE's   ***
! *********************************************

  if (rank==0) then

!   Loop over grid levels, we solve from fine to coarse grid.

     do l=Nl-1,0,-1


!       *******************************
!       ***   FIND INITIAL VALUES   ***
!       *******************************

!       Background value of trK. We are assuming that trK
!       is uniform (constant throughout the grid).

        trK_rk = trK_g(l,Nr)

!       Find initial point. Only the finest grid
!       integrates from the origin.

        if (l==Nl-1) then
           imin = 2
        else
           imin = Nrtotal/2
        end if

!       For coarse grids we interpolate the initial point
!       (which is half way through the grid) from the last
!       points of the higher resolution grid.  Here we use
!       cubic interpolation.

        if (l<Nl-1) then
           psi_g(l,imin-1)  = (9.d0*(psi_g(l+1,Nrtotal-2)+psi_g(l+1,Nrtotal-3)) &
                            - (psi_g(l+1,Nrtotal-4)+psi_g(l+1,Nrtotal-1)))/16.d0
           Dpsi_g(l,imin-1) = (9.d0*(Dpsi_g(l+1,Nrtotal-2)+Dpsi_g(l+1,Nrtotal-3)) &
                            - (Dpsi_g(l+1,Nrtotal-4)+Dpsi_g(l+1,Nrtotal-1)))/16.d0
           Khat_g(l,imin-1) = (9.d0*(Khat_g(l+1,Nrtotal-2)+Khat_g(l+1,Nrtotal-3)) &
                            - (Khat_g(l+1,Nrtotal-4)+Khat_g(l+1,Nrtotal-1)))/16.d0
        end if


!       ************************************
!       ***   SECOND ORDER RUNGE-KUTTA   ***
!       ************************************

!       Remember that we integrate the equation for Khat:=KTA/r
!       instead of KTA, and recover KTA at the end.

        if (order=='two') then

!          First grid point.  We assume that close to the origin
!          we have:
!                               2
!          psi   ~  psi0  +  a r    =>   Dpsi  ~ 2 a r
!
!          Khat  ~  b r
!
!          Substituting in the constraints we find that:
!
!                     2                           5
!          a  =  [ trK / 72  -  (pi/3) rho0 ] psi0
!
!          b  =  (8pi/5) (JA/r)
!
!          We then use these expansions at the very first grid point.

           if (l==Nl-1) then

!             Energy density at r=0. We use linear interpolation.

              pi2_rk = half*(pi2_g(l,0) + pi2_g(l,1))
              V_rk   = half*(V_g(l,0) + V_g(l,1))

              A_rk   = half*(A_g(l,0) + A_g(l,1))

              rho_rk = half*pi2_rk + V_rk

!             Variables at i=1.

              aux = A_rk*(trK_rk**2/72.d0 - third*smallpi*rho_rk)*psiorigin**5

              psi_g(l,1)  = psiorigin + aux*rr(l,1)**2
              Dpsi_g(l,1) = two*aux*rr(l,1)

              Khat_g(l,1) = 8.d0*smallpi/5.d0*JA_g(l,1)

           end if

!          Points away from the origin.

           do i=imin,Nrtotal

              r0    = rr(l,i-1)
              delta = dr(l)

              psi0  = psi_g(l,i-1)
              Dpsi0 = Dpsi_g(l,i-1)
              Khat0 = Khat_g(l,i-1)

!             I) First Runge-Kutta step.

!             Variables at point i-1.

              rm = r0

              psi_rk  = psi0
              Dpsi_rk = Dpsi0
              Khat_rk = Khat0

!             Value of (xi,pi,V,JA,A,B) at point i-1.

              xi2_rk = xi2_g(l,i-1)
              pi2_rk = pi2_g(l,i-1)

              V_rk   = V_g(l,i-1)
              JA_rk  = JA_g(l,i-1)

              A_rk = A_g(l,i-1)
              B_rk = B_g(l,i-1)

              DA_rk = DA_g(l,i-1)
              DB_rk = DB_g(l,i-1)

!             Energy density.

              rho_rk = half*(pi2_rk + xi2_rk/A_rk/psi_rk**4) + V_rk

!             Sources for (psi,Dpsi,Khat).

              k11 = Dpsi_rk
              k21 = - (two/rm - 0.5d0*DA_rk/A_rk + DB_rk/B_rk)*Dpsi_rk &
                  - A_rk*psi_rk**5*((1.5d0*(Khat_rk*rm)**2 - two*third*trK_rk**2)/8.d0 + two*smallpi*rho_rk)
              k31 = - 2.d0*(Khat_rk*(2.d0/rm + 0.75d0*DB_rk/B_rk + 3.d0*Dpsi_rk/psi_rk) - 4.d0*smallpi*JA_rk/rm)

!             II) Second Runge-Kutta step.

!             Variables at mid point.

              rm = r0 + half*delta

              psi_rk  = psi0  + half*delta*k11
              Dpsi_rk = Dpsi0 + half*delta*k21
              Khat_rk = Khat0 + half*delta*k31

!             Interpolated value of (xi,pi,V,JA,A,B) at mid point.

              xi2_rk = half*(xi2_g(l,i-1) + xi2_g(l,i))
              pi2_rk = half*(pi2_g(l,i-1) + pi2_g(l,i))

              V_rk   = half*(V_g(l,i-1) + V_g(l,i))
              JA_rk  = half*(JA_g(l,i-1) + JA_g(l,i))

              A_rk  = half*(A_g(l,i-1) + A_g(l,i))
              B_rk  = half*(B_g(l,i-1) + B_g(l,i))

              DA_rk  = half*(DA_g(l,i-1) + DA_g(l,i))
              DB_rk  = half*(DB_g(l,i-1) + DB_g(l,i))

!             Energy density at mid point.

              rho_rk = half*(pi2_rk + xi2_rk/A_rk/psi_rk**4) + V_rk

!             Sources for (psi,Dpsi,Khat).

              k12 = Dpsi_rk
              k22 = - (two/rm - 0.5d0*DA_rk/A_rk + DB_rk/B_rk)*Dpsi_rk &
                  - A_rk*psi_rk**5*((1.5d0*(Khat_rk*rm)**2 - two*third*trK_rk**2)/8.d0 + two*smallpi*rho_rk)
              k32 = - 2.d0*(Khat_rk*(2.d0/rm + 0.75d0*DB_rk/B_rk + 3.d0*Dpsi_rk/psi_rk) - 4.d0*smallpi*JA_rk/rm)

!             III) Advance variables.

              psi_g(l,i)  = psi0  + delta*k12
              Dpsi_g(l,i) = Dpsi0 + delta*k22
              Khat_g(l,i) = Khat0 + delta*k32

           end do

!          Recover KTA from Khat.

           KTA_g = rr*Khat_g


!       ************************************
!       ***   FOURTH ORDER RUNGE-KUTTA   ***
!       ************************************

!       Remember that we integrate the equation for Khat:=KTA/r
!       instead of KTA, and recover KTA at the end.

        else

!          First grid point.  We assume that close to the origin
!          we have:
!                               2
!          psi   ~  psi0  +  a r    =>   Dpsi  ~ 2 a r
!
!          Khat  ~  b r
!
!          Substituting in the constraints we find that:
!
!                     2                           5
!          a  =  [ trK / 72  -  (pi/3) rho0 ] psi0
!
!          b  =  (8pi/5) (JA/r)
!
!          We then use these expansions at the very first grid point.

           if (l==Nl-1) then

!             Energy density at r=0.  We use cubic interpolation.

              pi2_rk = (9.d0*(pi2_g(l,0) + pi2_g(l,1)) &
                     - (pi2_g(l,-1) + pi2_g(l,2)))/16.d0
              V_rk   = (9.d0*(V_g(l,0) + V_g(l,1)) &
                     - (V_g(l,-1) + V_g(l,2)))/16.d0

              A_rk   = (9.d0*(A_g(l,0) + A_g(l,1)) &
                     - (A_g(l,-1) + A_g(l,2)))/16.d0

              rho_rk = half*pi2_rk + V_rk

!             Variables at i=1.

              aux = A_rk*(trK_rk**2/72.d0 - third*smallpi*rho_rk)*psiorigin**5

              psi_g(l,1)  = psiorigin + aux*rr(l,1)**2
              Dpsi_g(l,1) = two*aux*rr(l,1)

              Khat_g(l,1) = (8.d0*smallpi/5.d0)*JA_g(l,1)

           end if

!          Points away from the origin.

           do i=imin,Nrtotal

              r0    = rr(l,i-1)
              delta = dr(l)

              psi0  = psi_g(l,i-1)
              Dpsi0 = Dpsi_g(l,i-1)
              Khat0 = Khat_g(l,i-1)

!             I) First Runge-Kutta step.

!             Variables at point i-1.

              rm = r0

              psi_rk  = psi0
              Dpsi_rk = Dpsi0
              Khat_rk = Khat0

!             Value of (xi,pi,V,JA) at point i-1.

              xi2_rk = xi2_g(l,i-1)
              pi2_rk = pi2_g(l,i-1)

              V_rk  = V_g(l,i-1)
              JA_rk = JA_g(l,i-1)

              A_rk = A_g(l,i-1)
              B_rk = B_g(l,i-1)

              DA_rk = DA_g(l,i-1)
              DB_rk = DB_g(l,i-1)

!             Energy density.

              rho_rk = half*(pi2_rk + xi2_rk/A_rk/psi_rk**4) + V_rk

!             Sources for (psi,Dpsi,Khat).

              k11 = Dpsi_rk
              k21 = - (two/rm - 0.5d0*DA_rk/A_rk + DB_rk/B_rk)*Dpsi_rk &
                  - A_rk*psi_rk**5*((1.5d0*(Khat_rk*rm)**2 - two*third*trK_rk**2)/8.d0 + two*smallpi*rho_rk)
              k31 = - 2.d0*(Khat_rk*(2.d0/rm + 0.75d0*DB_rk/B_rk + 3.d0*Dpsi_rk/psi_rk) - 4.d0*smallpi*JA_rk/rm)

!             II) Second Runge-Kutta step.

!             Variables at mid point.

              rm = r0 + half*delta

              psi_rk  = psi0  + half*delta*k11
              Dpsi_rk = Dpsi0 + half*delta*k21
              Khat_rk = Khat0 + half*delta*k31

!             Interpolated value of (xi,pi,V,JA) at mid point
!             (careful with the last point).

              if (i==Nrtotal) then  ! Second order at the moment, can fix it later.

                 xi2_rk = half*(xi2_g(l,i-1) + xi2_g(l,i))
                 pi2_rk = half*(pi2_g(l,i-1) + pi2_g(l,i))

                 V_rk  = half*(V_g(l,i-1) + V_g(l,i))
                 JA_rk = half*(JA_g(l,i-1) + JA_g(l,i))

                 A_rk = half*(A_g(l,i-1) + A_g(l,i))
                 B_rk = half*(B_g(l,i-1) + B_g(l,i))

                 DA_rk = half*(DA_g(l,i-1) + DA_g(l,i))
                 DB_rk = half*(DB_g(l,i-1) + DB_g(l,i))

              else

                 xi2_rk = (9.d0*(xi2_g(l,i-1) + xi2_g(l,i)) &
                        - (xi2_g(l,i-2) + xi2_g(l,i+1)))/16.d0
                 pi2_rk = (9.d0*(pi2_g(l,i-1) + pi2_g(l,i)) &
                        - (pi2_g(l,i-2) + pi2_g(l,i+1)))/16.d0

                 V_rk  = (9.d0*(V_g(l,i-1) + V_g(l,i)) &
                       - (V_g(l,i-2) + V_g(l,i+1)))/16.d0
                 JA_rk = (9.d0*(JA_g(l,i-1) + JA_g(l,i)) &
                       - (JA_g(l,i-2) + JA_g(l,i+1)))/16.d0

                 A_rk = (9.d0*(A_g(l,i-1) + A_g(l,i)) &
                      - (A_g(l,i-2) + A_g(l,i+1)))/16.d0
                 B_rk = (9.d0*(B_g(l,i-1) + B_g(l,i)) &
                      - (B_g(l,i-2) + B_g(l,i+1)))/16.d0

                 DA_rk = (9.d0*(DA_g(l,i-1) + DA_g(l,i)) &
                       - (DA_g(l,i-2) + DA_g(l,i+1)))/16.d0
                 DB_rk = (9.d0*(DB_g(l,i-1) + DB_g(l,i)) &
                       - (DB_g(l,i-2) + DB_g(l,i+1)))/16.d0

              end if

!             Energy density at mid point.

              rho_rk = half*(pi2_rk + xi2_rk/A_rk/psi_rk**4) + V_rk

!             Sources for (psi,Dpsi,Khat).

              k12 = Dpsi_rk
              k22 = - (two/rm - 0.5d0*DA_rk/A_rk + DB_rk/B_rk)*Dpsi_rk & 
                  - A_rk*psi_rk**5*((1.5d0*(Khat_rk*rm)**2 - two*third*trK_rk**2)/8.d0 + two*smallpi*rho_rk)
              k32 = - 2.d0*(Khat_rk*(2.d0/rm + 0.75d0*DB_rk/B_rk + 3.d0*Dpsi_rk/psi_rk) - 4.d0*smallpi*JA_rk/rm)

!             III) Third Runge-Kutta step.

!             New variables at mid point.

              psi_rk  = psi0  + half*delta*k12
              Dpsi_rk = Dpsi0 + half*delta*k22
              Khat_rk = Khat0 + half*delta*k32

!             New energy density at mid point.

              rho_rk = half*(pi2_rk + xi2_rk/A_rk/psi_rk**4) + V_rk

!             Sources for (psi,Dpsi,Khat).

              k13 = Dpsi_rk
              k23 = - (two/rm - 0.5d0*DA_rk/A_rk + DB_rk/B_rk)*Dpsi_rk &
                  - A_rk*psi_rk**5*((1.5d0*(Khat_rk*rm)**2 - two*third*trK_rk**2)/8.d0 + two*smallpi*rho_rk)
              k33 = - 2.d0*(Khat_rk*(2.d0/rm + 0.75d0*DB_rk/B_rk + 3.d0*Dpsi_rk/psi_rk) - 4.d0*smallpi*JA_rk/rm)

!             IV) Fourth Runge-Kutta step.

!             Variables at point i.

              rm = r0 + delta

              psi_rk  = psi0  + delta*k13
              Dpsi_rk = Dpsi0 + delta*k23
              Khat_rk = Khat0 + delta*k33

!             Value of (xi,pi,V,JA) at point i.

              xi2_rk = xi2_g(l,i)
              pi2_rk = pi2_g(l,i)

              V_rk  = V_g(l,i)
              JA_rk = JA_g(l,i)

              A_rk = A_g(l,i)
              B_rk = B_g(l,i)

              DA_rk = DA_g(l,i)
              DB_rk = DB_g(l,i)

!             Energy density at point i.

              rho_rk = half*(pi2_rk + xi2_rk/A_rk/psi_rk**4) + V_rk

!             Sources for (psi,Dpsi,Khat).

              k14 = Dpsi_rk
              k24 = - (two/rm - 0.5d0*DA_rk/A_rk + DB_rk/B_rk)*Dpsi_rk &
                  - A_rk*psi_rk**5*((1.5d0*(Khat_rk*rm)**2 - two*third*trK_rk**2)/8.d0 + two*smallpi*rho_rk)
              k34 = - 2.d0*(Khat_rk*(2.d0/rm + 0.75d0*DB_rk/B_rk + 3.d0*Dpsi_rk/psi_rk) - 4.d0*smallpi*JA_rk/rm)

!             V) Advance variables.

              psi_g(l,i)  = psi0  + delta*(k11 + 2.d0*(k12 + k13) + k14)/6.d0
              Dpsi_g(l,i) = Dpsi0 + delta*(k21 + 2.d0*(k22 + k23) + k24)/6.d0
              Khat_g(l,i) = Khat0 + delta*(k31 + 2.d0*(k32 + k33) + k34)/6.d0

           end do

!          Recover KTA from Khat.

           KTA_g = rr*Khat_g

        end if

     end do


!    ***********************
!    ***   GHOST ZONES   ***
!    ***********************

!    Ghost zones only for fine grid.

     do i=1,ghost

        psi_g(Nl-1,1-i)  = + psi_g(Nl-1,i)
        Dpsi_g(Nl-1,1-i) = - Dpsi_g(Nl-1,i)

        Khat_g(Nl-1,1-i) = + Khat_g(Nl-1,i)
        KTA_g(Nl-1,1-i)  = + KTA_g(Nl-1,i)

     end do


!    ****************************************
!    ***   FIND ASYMPTOTIC VALUE OF PSI   ***
!    ****************************************

!    We need to find the asymtptotic value of psi in order to rescale it.
!    To find this value we assume we have a fall-off of the form:
!
!    u = u0 + k/r   =>   du/dr = (u0-u)/r
!
!    Solving for u0 we find:
!
!    u0 = u + r (du/dr)
!
!    We then approximate the derivative at the boundary using one-sided
!    differences.

     if (order=="two") then
        psiasymp = psi_g(0,Nrtotal) + rr(0,Nrtotal)*0.5d0/dr(0) &
            *(3.d0*psi_g(0,Nrtotal) - 4.d0*psi_g(0,Nrtotal-1) + psi_g(0,Nrtotal-2))
     else
        psiasymp = psi_g(0,Nrtotal) + rr(0,Nrtotal)*0.25d0/dr(0) &
            *(25.d0*psi_g(0,Nrtotal) - 48.d0*psi_g(0,Nrtotal-1) &
            + 36.d0*psi_g(0,Nrtotal-2) - 16.d0*psi_g(0,Nrtotal-3) + 3.d0*psi_g(0,Nrtotal-4))/3.d0
     end if

     if (.not.minimize) then
        print *, 'Asymptotic value of conformal factor: ',psiasymp
        if (rescale) print *
     end if


!    ************************************
!    ***   RESTRICT TO COARSE GRIDS   ***
!    ************************************

!    Restrict solution from fine to coarse grid.
!    We don't call the subroutine "restrict"
!    since here we are running only on processor 0.
!    We use cubic interpolation.

     do l=Nl-1,1,-1

        do i=1,Nrtotal-ghost,2
           psi_g(l-1,i/2+1)  = (9.d0*(psi_g(l,i)+psi_g(l,i+1))   - (psi_g(l,i-1)+psi_g(l,i+2)))/16.d0
           Dpsi_g(l-1,i/2+1) = (9.d0*(Dpsi_g(l,i)+Dpsi_g(l,i+1)) - (Dpsi_g(l,i-1)+Dpsi_g(l,i+2)))/16.d0
           KTA_g(l-1,i/2+1)  = (9.d0*(KTA_g(l,i)+KTA_g(l,i+1))   - (KTA_g(l,i-1)+KTA_g(l,i+2)))/16.d0
           Khat_g(l-1,i/2+1) = (9.d0*(Khat_g(l,i)+Khat_g(l,i+1)) - (Khat_g(l,i-1)+Khat_g(l,i+2)))/16.d0
        end do

!       Fix ghost zones.

        do i=1,ghost
           psi_g(l-1,1-i)  = + psi_g(l-1,i)
           Dpsi_g(l-1,1-i) = - Dpsi_g(l-1,i)
           KTA_g(l-1,1-i)  = + KTA_g(l-1,i)
           Khat_g(l-1,1-i) = + Khat_g(l-1,i)
        end do

     end do


! *************************************
! ***   FINISHED FINDING SOLUTION   ***
! *************************************

! When we get here the solution has been found,
! so we close the "if" statement for processor 0.

  end if


! ************************************
! ***   RESCALE CONFORMAL FACTOR   ***
! ************************************

! Broadcast asymptotic value of psi to all processors.

  call MPI_BCAST(psiasymp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! Rescale psi and initial profile width as explained in the comment
! at the top of the routine.

  if (.not.rescale) then

!    Rescale solution.

     rescale = .true.
     psiorigin = 1.d0/psiasymp

     aux = psiasymp**2
     call rescaleprofile(aux,.false.)

!    Messages to screen.

     if (.not.minimize) then

        print *, 'Restarting integration with rescaled conformal factor at origin'

        if (complexprofile=="gaussian") then

           print *, 'Gaussian parameters changed to: '

           if (complexDM_type/='growing') then
              print *, 'complexR_r0 = ',complexR_r0
              print *, 'complexR_s0 = ',complexR_s0
           end if

           print *, 'complexI_r0 = ',complexI_r0
           print *, 'complexI_s0 = ',complexI_s0

        else if (complexprofile=="tophat") then

           print *, 'Top-hat parameters changed to: '

           if (complexDM_type/='growing') then
              print *, 'complexR_r0 = ',complexR_r0
              print *, 'complexR_s0 = ',complexR_s0
              print *, 'complexR_t0 = ',complexR_t0
           end if

           print *, 'complexI_r0 = ',complexI_r0
           print *, 'complexI_s0 = ',complexI_s0
           print *, 'complexI_t0 = ',complexI_t0

        else if (complexprofile=="comp-gaussian") then

           if (complexDM_type/='growing') then
              print *, 'complexR_r0 = ',complexR_r0
              print *, 'complexR_s0 = ',complexR_s0
              print *, 'complexR_r1 = ',complexR_r1
              print *, 'complexR_s1 = ',complexR_s1
           end if

           print *, 'complexI_r0 = ',complexI_r0
           print *, 'complexI_s0 = ',complexI_s0
           print *, 'complexI_r1 = ',complexI_r1
           print *, 'complexI_s1 = ',complexI_s1

        else if (complexprofile=="comp-tophat") then

           print *, 'Top-hat parameters changed to: '

           if (complexDM_type/='growing') then
              print *, 'complexR_r0 = ',complexR_r0
              print *, 'complexR_s0 = ',complexR_s0
              print *, 'complexR_t0 = ',complexR_t0
              print *, 'complexR_r1 = ',complexR_r1
              print *, 'complexR_s1 = ',complexR_s1
              print *, 'complexR_t1 = ',complexR_t1
           end if

           print *, 'complexI_r0 = ',complexI_r0
           print *, 'complexI_s0 = ',complexI_s0
           print *, 'complexI_t0 = ',complexI_t0
           print *, 'complexI_r1 = ',complexI_r1
           print *, 'complexI_s1 = ',complexI_s1
           print *, 'complexI_t1 = ',complexI_t1

        end if

     end if

     goto 100

  end if


! ************************************************
! ***   DISTRIBUTE SOLUTION AMONG PROCESSORS   ***
! ************************************************

! For parallel runs, when we get here the solution
! is known only on processor zero for the full size
! array with dimensions Nrtotal.  We must now distribute
! the solution among all other processors.

  if (size==1) then
     psi = psi_g
     KTA = KTA_g
  else
     call distribute(0,Nl-1,psi,psi_g)
     call distribute(0,Nl-1,KTA,KTA_g)
  end if

! Remember that KTB = - KTA/2.

  KTB = - half*KTA


! ***************************************
! ***   FIND CONFORMAL FUNCTION phi   ***
! ***************************************

! Derivative of psi.

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

  end subroutine complexDMmom_pert









  subroutine rescaleprofile(factor,initial)

! ***************************
! ***   RESCALE PROFILE   ***
! ***************************

! Rescale initial profile parameters.

  use param

  implicit none

  logical initial

  real(8) factor
  real(8) R_r0,R_s0,R_t0,R_r1,R_s1,R_t1
  real(8) I_r0,I_s0,I_t0,I_r1,I_s1,I_t1

  common / initialprofile / R_r0,R_s0,R_t0,R_r1,R_s1,R_t1,I_r0,I_s0,I_t0,I_r1,I_s1,I_t1

! Rescale from original profile.

  if (initial) then

     complexR_r0 = R_r0*factor
     complexR_s0 = R_s0*factor
     complexR_t0 = R_t0*factor

     complexI_r0 = I_r0*factor
     complexI_s0 = I_s0*factor
     complexI_t0 = I_t0*factor

     if (complexprofile(1:4)=="comp") then

        complexR_r1 = R_r1*factor
        complexR_s1 = R_s1*factor
        complexR_t1 = R_t1*factor

        complexI_r1 = I_r1*factor
        complexI_s1 = I_s1*factor
        complexI_t1 = I_t1*factor

     end if

! Rescale from current profile.

  else

     complexR_r0 = complexR_r0*factor
     complexR_s0 = complexR_s0*factor
     complexR_t0 = complexR_t0*factor

     complexI_r0 = complexI_r0*factor
     complexI_s0 = complexI_s0*factor
     complexI_t0 = complexI_t0*factor

     if (complexprofile(1:4)=="comp") then

        complexR_r1 = complexR_r1*factor
        complexR_s1 = complexR_s1*factor
        complexR_t1 = complexR_t1*factor

        complexI_r1 = complexI_r1*factor
        complexI_s1 = complexI_s1*factor
        complexI_t1 = complexI_t1*factor

     end if

  end if

  end subroutine rescaleprofile









  real(8) function func(x)

! ********************************
! ***   FUNCTION TO MINIMIZE   ***
! ********************************

! Here we define a function that measures the
! distance between the latest profile and the
! original profile parameters.

  use param

  implicit none

  real(8) x
  real(8) R_r0,R_s0,R_t0,R_r1,R_s1,R_t1
  real(8) I_r0,I_s0,I_t0,I_r1,I_s1,I_t1

  common / initialprofile / R_r0,R_s0,R_t0,R_r1,R_s1,R_t1,I_r0,I_s0,I_t0,I_r1,I_s1,I_t1

! Rescale the initial profile by a factor x from the
! original parameters.

  call rescaleprofile(x,.true.)
  call complexDMmom_pert

! Define the distance, for this we simply compare the
! parameters s0. But beware, for the "growing" type
! the real part is ignored.

  if (complexDM_type/='growing') then
     func = (R_s0 - complexR_s0)**2 + (I_s0 - complexI_s0)**2
  else
     func = (I_s0 - complexI_s0)**2
  end if

! Simple quadratic function, for testing.

! func = 1.0 + (x-2.d0)**2


! ***************
! ***   END   ***
! ***************

  end function func










  subroutine mnbrak(ax,bx,cx,fa,fb,fc)

! *****************************************
! ***   MNBRAK FROM NUMERICAL RECIPES   ***
! *****************************************

! Given a function "func", and distincrt initial points ax and bx,
! this routine searches in the downhill direction and returns new
! points (ax,bx,cx) that bracket the minimun of the function.
! Also returned are the function values at the three points (fa,fb,fc).

  implicit none

  real(8) ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
  real(8) dum,fu,q,r,u,ulim,cu

  external func
  parameter (GOLD=1.618034,GLIMIT=100.0,TINY=1.d-20)

! Evaluate function at ax and bx.

  fa = func(ax)
  fb = func(bx)

! Swithc roles of ax and bx so that we go in the downhill direction.
 
  if (fb.gt.fa) then
    dum = ax
    ax = bx
    bx = dum
    dum = fb
    fb = fa
    fa = dum
  end if

! First guess for cx.

  cx = bx + GOLD*(bx-ax)
  fc = func(cx)

! Keep returning here until we bracket.

  1 continue

  if (fb.gt.fc) then

!   Compute u by parabolic extrapolation from a,b,c.
!   TINY is used to prevent possible divisions by zero.

    r = (bx-ax)*(fb-fc)
    q = (bx-cx)*(fb-fa)
    u = bx - ((bx-cx)*q - (bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))

!   We won't go further than this.

    ulim = bx + GLIMIT*(cx-bx)

!   Check if parabolic u is between b and c.
 
    if ((bx-u)*(u-cx).gt.0.d0) then

       fu=func(u)

!      Got a minimum between b and c.

       if (fu.lt.fc) then

          ax = bx
          fa = fb
          bx = u
          fb = fu
          return

!      Got a minimum between a and u.

       else if (fu.gt.fb) then

          cx = u
          fc = fu
          return

       end if

!      Parabolic fit was no use.  Use dafaut magnification.

       u = cx + GOLD*(cx-bx)
       fu = func(u)

!   Parabolic fit is between c and its allowed limit.

    else if ((cx-cu)*(u-ulim).gt.0.d0) then

       fu = func(u)

       if (fu.lt.fc) then
          bx = cx
          cx = u
          u = cx + GOLD*(cx-bx)
          fb = fc
          fc = fu
          fu = func(u)
       end if

!   Limit parabolic u to maximum allowed value.

    else if ((u-ulim)*(ulim-cx).ge.0.d0) then

       u=ulim
       fu = func(u)

!   Reject parabolic u, use default magnification.

    else

       u = cx + GOLD*(cx-bx)
       fu = func(u)

    end if

!   Eliminate oldest point and continue
    ax = bx
    bx = cx
    cx = u
    fa = fb
    fb = fc
    fc = fu
    goto 1

  end if


! ***************
! ***   END   ***
! ***************

  return
  end









  subroutine golden(ax,bx,cx,xmin,fmin)

! *****************************************
! ***   GOLDEN FROM NUMERICAL RECIPES   ***
! *****************************************

! Given a function "func", and a bracketing triplet of abcissas (ax,bx,cx),
! this routine performs a golen section search for the minimum, isolating
! it to a fractional precision "tol". The location of the minimum is returned
! as "xmin", and the minimum function value as "fmin".

  implicit none

  real(8) ax,bx,cx,xmin,fmin
  real(8) tol,func,R,C
  real(8) f1,f2,x0,x1,x2,x3

  external func
  parameter(tol=1.d-10,R=0.61803399d0,C=1.d0-R)

! At any given point we keep rack of fou points: x0,x1,x2,x3.

  x0 = ax
  x3 = cx

! Make x0 to x1 the smallest segment and fill in
! the new point to be tried.

  if (abs(cx-bx).gt.abs(bx-ax)) then
     x1 = bx
     x2 = bx + C*(cx-bx)
  else
     x2 = bx
     x1 = bx - C*(bx-ax)
  end if

! Initial function evaluations.

  f1 = func(x1)
  f2 = func(x2)

! We keep returning here.

  1 continue

  if (abs(x3-x0).gt.tol*(abs(x1)+abs(x2))) then

     if (f2.lt.f1) then
        x0 = x1
        x1 = x2
        x2 = R*x1 + C*x3
        f1 = f2
        f2 = func(x2)
     else
        x3 = x2
        x2 = x1
        x1 = R*x2 + C*x0
        f2 = f1
        f1 = func(x1)
     end if

     goto 1

  end if

! We are done. Output the best of the two current value.

  if (f1.lt.f2) then
     fmin = f1
     xmin = x1
  else
     fmin = f2
     xmin = x2
  end if


! ***************
! ***   END   ***
! ***************

  return
  end




