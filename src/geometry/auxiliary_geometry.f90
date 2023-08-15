!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/auxiliary_geometry.f90,v 1.54 2023/03/09 23:29:22 malcubi Exp $

  subroutine auxiliary_geometry(l)

! *******************************************
! ***   AUXLIARY VARIABLES FOR GEOMETRY   ***
! *******************************************

! Include modules.

  use param
  use arrays
  use derivatives
  use derivadvect
  use procinfo

! Extra variables.

  implicit none

  integer i,l

  real(8) zero,half,third,one,two


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  one  = 1.d0
  two  = 2.d0

  half  = 0.5d0
  third = 1.d0/3.d0


! ****************************
! ***   CONFORMAL FACTOR   ***
! ****************************

! Calculate (phi,chi) and their derivatives. Notice that at t=0
! chi might not be defined, so we always need to define it.

  if (.not.chimethod) then
     chi(l,:) = exp(-dble(chipower)*phi(l,:))
  else
     phi(l,:) = - log(abs(chi(l,:)))/dble(chipower)
  end if

! Derivatives.

  diffvar => phi

  if (dorigin=="onesided") then
     D1_phi(l,:) = diff1(l,+1,.true.)
     D2_phi(l,:) = diff2(l,+1,.true.)
  else
     D1_phi(l,:) = diff1(l,+1)
     D2_phi(l,:) = diff2(l,+1)
  end if

  diffvar => chi
  D1_chi(l,:) = diff1(l,+1)
  D2_chi(l,:) = diff2(l,+1)

! Find psi and derivatives.

  psi(l,:)  = dexp(phi(l,:))
  psi2(l,:) = psi(l,:)**2
  psi4(l,:) = psi(l,:)**4

  diffvar => psi

  if (dorigin=="onesided") then
     D1_psi(l,:) = diff1(l,+1,.true.)
     D2_psi(l,:) = diff2(l,+1,.true.)
  else
     D1_psi(l,:) = diff1(l,+1)
     D2_psi(l,:) = diff2(l,+1)
  end if

  !D1_psi(l,:) = psi(l,:)*D1_phi(l,:)
  !D2_psi(l,:) = psi(l,:)*(D1_phi(l,:)**2 + D2_phi(l,:))

!                            2
! DD_phir = d ( d phi / r psi  )
!            r   r

  auxarray(l,:) = D1_phi(l,:)/r(l,:)/psi2(l,:)

  diffvar => auxarray
  DD_phir(l,:) = diff1(l,+1)


! *****************
! ***   LAPSE   ***
! *****************

! Derivatives of lapse.

  diffvar => alpha
  D1_alpha(l,:) = diff1(l,+1)
  D2_alpha(l,:) = diff2(l,+1)

!                                4
! DD_alphar = d ( d alpha / r psi )
!              r   r

  auxarray(l,:) = D1_alpha(l,:)/r(l,:)/psi4(l,:)

  diffvar => auxarray
  DD_alphar(l,:) = diff1(l,+1)


! **************************************
! ***   BONA-MASSO GAUGE FUNCTIONS   ***
! **************************************

! Here we calculate the Bona-Masso gauge function
! which is used for the slicing condition, and
! also to calculate gauge speeds. Notice that I
! actually define falpha as:
!
! falpha = alpha**2 f(alpha)
!
! This is in order to avoid divisions with
! alpha for small alpha.

! GEODESIC OR MAXIMAL.

  if ((slicing=="static").or.(slicing=="maximal")) then

!    In this case we just set falpha it to zero.

     falpha(l,:) = 0.d0

     if (cosmic_run) then
        cosmobg_falpha(l) = 0.d0
     end if

! HARMONIC FAMILY.

  else if (index(adjustl(slicing),"harmonic")==1) then

!    In this case the gauge function is:
!
!    f  =  gauge_f   =>  falpha = gauge_f*alpha**2
!
!    True harmonic slicing corresponds to the case when
!    the constant gauge_f is equal to 1.
!
!    The conformal foliation corresponds to alpha=psi**2,
!    which implies f=1/3.

     falpha(l,:) = gauge_f*alpha(l,:)**2

     if (cosmic_run) then
        cosmobg_falpha(l) = gauge_f*cosmobg_alpha(l)**2
     end if

! 1+LOG FAMILY.

  else if (index(adjustl(slicing),"1+log")==1) then

!    In this case the gauge function is:
!
!    f  =  gauge_f/alpha  =>  falpha = gauge_f*alpha
!
!    Standard 1+log slicing corresponds to the case when
!    the constant gauge_f is equal to 2.

     falpha(l,:) = gauge_f*alpha(l,:)

     if (cosmic_run) then
        cosmobg_falpha(l) = gauge_f*cosmobg_alpha(l)
     end if

! SHOCK AVOIDING FAMILY.

  else if (index(adjustl(slicing),"shockavoid")==1) then

!    In this case the gauge function is:
!
!    f = 1 + k/alpha^2  =>  falpha = alpha**2 + k
!
!    with k=constant. The whole family avoids shocks.
!    (See: Class.Quant.Grav. 20 (2003) 607-624; gr-qc/0210050)
!
!    Below we take k=gauge_f-1.  This is just to guarantee
!    that for alpha=1 we always have falpha=gauge_f (as in
!    the harmonic and 1+log cases above).

     falpha(l,:) = alpha(l,:)**2 + (gauge_f-1.d0)

     if (cosmic_run) then
        cosmobg_falpha(l) = cosmobg_alpha(l)**2 + (gauge_f-1.d0)
     end if

! TESTING NEW CONDITIONS.

  else if (index(adjustl(slicing),"alphaminus2")==1) then

     falpha(l,:) = gauge_f

! COSMOLOGICAL.

  else if (index(adjustl(slicing),"cosmo")==1) then

!    In this case the background lapse evolves in a
!    different way to the perturbed lapse.
!
!    For "comosync" we use a background synchronous lapse (alpha=1):
!
!    f_bg = 0
!
!    while for "cosmocf" we use a background conformal lapse:
!
!    f_bg = 1/3  =>  falpha_bg = alpha_bg/3
!
!    The perturbed part of the lapse can evolve with any of the
!    Bona-Masso standard types.

     if (cosmic_run) then

!       Background gauge function.

        if (index(adjustl(slicing),"cosmosync")==1) then
           cosmobg_falpha(l) = zero
        else if (index(adjustl(slicing),"cosmocf")==1) then
           cosmobg_falpha(l) = third*cosmobg_alpha(l)**2
        end if

!       Perturbed gauge function:

        if (index(adjustl(slicing),"harmonic")>0) then
           falpha(l,:) = gauge_f*alpha(l,:)**2
        else if (index(adjustl(slicing),"1+log")>0) then
           falpha(l,:) = gauge_f*alpha(l,:)
        else if (index(adjustl(slicing),"shockavoid")>0) then
           falpha(l,:) = alpha(l,:)**2 + (gauge_f-1.d0)
        end if

     else

        if (rank==0) then
           print *
           print *, 'Slicings of type "cosmo" must have cosmic_run=.true.'
           print *, 'Aborting! (subroutine auxiliary_geometry.f90)'
           print *
        end if

        call die

     end if

  end if


! **************************
! ***   SPATIAL METRIC   ***
! **************************

! For Lagrangian evolutions, find B from A.

  if ((t(l)>0.d0).and.(bssnflavor=="lagrangian")) then
     B(l,:) = sqrt(abs(AB2(l,:)/A(l,:)))
  end if

! Derivatives of conformal spatial metric A.

  diffvar => A
  D1_A(l,:) = diff1(l,+1)
  D2_A(l,:) = diff2(l,+1)

! Derivatives of conformal spatial metric B.

  diffvar => B
  D1_B(l,:) = diff1(l,+1)
  D2_B(l,:) = diff2(l,+1)

! DD_Ar = d ( d A / r )
!          r   r

  auxarray(l,:) = D1_A(l,:)/r(l,:)
  diffvar => auxarray
  DD_Ar(l,:) = diff1(l,+1)

! DD_Br = d ( d B / r )
!          r   r

  auxarray(l,:) = D1_B(l,:)/r(l,:)
  diffvar => auxarray
  DD_Br(l,:) = diff1(l,+1)

! Physical spatial metric.

  if (allocated(APHYS)) then
     APHYS(l,:) = psi4(l,:)*A(l,:)
  end if

  if (allocated(BPHYS)) then
     BPHYS(l,:) = psi4(l,:)*B(l,:)
  end if

! Determinant of conformal metric (without the r**2): A*B**2.
! For lagrangian evolutions it should be time independent,
! while for eulerian evolutions it should change in time.
!
! Notice that for lagrangian evolutions the code only calculates
! this quantity for the first time step, since it should remain
! static. If you want to test if it really does you can just
! comment the if statement.

  if ((t(l)==0.d0).or.(bssnflavor=="eulerian")) then

     AB2(l,:) = A(l,:)*B(l,:)**2

     diffvar => AB2
     D1_AB2(l,:) = diff1(l,+1)
     D2_AB2(l,:) = diff2(l,+1)

  end if


! *****************
! ***   SHIFT   ***
! *****************

  if (shift/="none") then

!    Derivatives of shift.

     diffvar => beta

     D1_beta(l,:) = diff1(l,-1)
     D2_beta(l,:) = diff2(l,-1)
     DA_beta(l,:) = diffadv(l,-1)

!    Derivatives of dtbeta.

     diffvar => dtbeta

     D1_dtbeta(l,:) = diff1(l,-1)
     DA_dtbeta(l,:) = diffadv(l,-1)

!    Derivatives of fdriver.

     diffvar => fdriver

     D1_fdriver(l,:) = diff1(l,+1)
     DA_fdriver(l,:) = diffadv(l,+1)

!    Conformal divergence of shift. Notice that since this
!    requires derivatives of the metric it must be calculated
!    after those derivatives have been computed.

!    DIV_beta(l,:) = D1_beta(l,:) + beta(l,:)*(half*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:)) &
!                  + two*beta(l,:)/r(l,:)

     DIV_beta(l,:) = D1_beta(l,:) + half*beta(l,:)*D1_AB2(l,:)/AB2(l,:) &
                   + two*beta(l,:)/r(l,:)

!    DD_beta = d ( beta / r )
!               r

     auxarray(l,:) = beta(l,:)/r(l,:)
     diffvar => auxarray
     DD_beta(l,:) = diff1(l,+1)

!    Derivative of conformal divergence of shift.

!    D1_DIV_beta(l,:) = D2_beta(l,:) + D1_beta(l,:)*(half*D1_A(l,:)/A(l,:) + D1_B(l,:)/B(l,:)) &
!                     + beta(l,:)*(half/A(l,:)*(D2_A(l,:) - D1_A(l,:)**2/A(l,:)) &
!                     + one/B(l,:)*(D2_B(l,:) - D1_B(l,:)**2/B(l,:))) + two*DD_beta(l,:)

     D1_DIV_beta(l,:) = D2_beta(l,:) + half*D1_beta(l,:)*D1_AB2(l,:)/AB2(l,:) &
                      + half*beta(l,:)*(D2_AB2(l,:)/AB2(l,:) - (D1_AB2(l,:)/AB2(l,:))**2) &
                      + two*DD_beta(l,:)

  end if


! *******************************
! ***   EXTRINSIC CURVATURE   ***
! *******************************

! Notice that in fact KTA and Klambda are proportional
! to each other through the relation below, which one
! might choose to use or ignore.  By default the code
! does not evolve KTA independently which seems to
! improve accuracy.

  if (noKTA) then
     if (.not.regular2) then
        KTA(l,:) = 2.d0*third*r(l,:)**2*Klambda(l,:)
     else
        KTA(l,:) = 2.d0*third*r(l,:)**2*Klambda2(l,:)*psi(l,:)**lambdapower
     end if
  end if

! Angular traceless-conformal extrinsic curvature.

  KTB(l,:) = - half*KTA(l,:)

! For the components of the physical extrinsic curvature
! we must undo the conformal rescaling and add the trace.

  KAPHYS(l,:) = KTA(l,:) + third*trK(l,:)
  KBPHYS(l,:) = KTB(l,:) + third*trK(l,:)

! Square of extrinsic curvature.

  K2(l,:) = KTA(l,:)**2 + two*KTB(l,:)**2 + trK(l,:)**2/3.d0

! Derivatives of trK.

  diffvar => trK
  D1_trK(l,:) = diff1(l,+1)

! Derivatives of KTA.

  diffvar => KTA
  D1_KTA(l,:) = diff1(l,+1)

! Derivatives of KTB.

  diffvar => KTB
  D1_KTB(l,:) = diff1(l,+1)


! *******************
! ***   LAMBDAS   ***
! *******************

! If we don't want to introduce the regularization variables
! lambda and Klambda, then here we calculate them in terms
! of metric and extrinsic curvature.

  if (nolambda) then

     if (noKTA) then
        print *
        print *, 'Having "nolambda" and "noKTA" both true is incompatible.'
        print *
        call die
     end if

     lambda(l,:)  = (one - A(l,:)/B(l,:))/r(l,:)**2
     Klambda(l,:) = 1.5d0*KTA(l,:)/r(l,:)**2

     lambda2(l,:)  =  lambda(l,:)/psi(l,:)**lambdapower
     Klambda2(l,:) = Klambda(l,:)/psi(l,:)**lambdapower

  end if

! Puncture regularization.

  if (regular2) then
     lambda(l,:)  =  lambda2(l,:)*psi(l,:)**lambdapower
     Klambda(l,:) = Klambda2(l,:)*psi(l,:)**lambdapower
  end if

! Derivatives of lambda.

  diffvar => lambda
  D1_lambda(l,:) = diff1(l,+1)
  D2_lambda(l,:) = diff2(l,+1)

! Derivatives of Klambda.

  diffvar => Klambda
  D1_Klambda(l,:) = diff1(l,+1)

! Derivatives of lambda2 and Klambda2.

  if (regular2) then

     diffvar => lambda2
     D1_lambda2(l,:) = diff1(l,+1)
     D2_lambda2(l,:) = diff2(l,+1)

     diffvar => Klambda2
     D1_Klambda2(l,:) = diff1(l,+1)

  end if


! **********************
! ***   BSSN Delta   ***
! **********************

! Delta in terms of its definition.

!$OMP PARALLEL DO SCHEDULE(GUIDED)
  do i=1-ghost,Nrmax
     DeltaAB(l,i) = (half*D1_A(l,i)/A(l,i) - D1_B(l,i)/B(l,i) &
                  - two*r(l,i)*lambda(l,i))/A(l,i)
  end do
!$OMP END PARALLEL DO

! If we don't want to evolve the BSSN variable Delta, then
! here we calculate it in terms of derivatives of the metric.

  if (noDeltar) then
     Deltar(l,:) = DeltaAB(l,:)
  end if

! Derivative of Deltar.

  diffvar => Deltar
  D1_Deltar(l,:) = diff1(l,-1)

!$OMP PARALLEL DO SCHEDULE(GUIDED)
  do i=1-ghost,Nrmax
     D1_DeltaAB(l,i) = - D1_A(l,i)/A(l,i)*DeltaAB(l,i) + one/A(l,i) &
                     *(half*(D2_A(l,i)/A(l,i) - (D1_A(l,i)/A(l,i))**2) &
                     - (D2_B(l,i)/B(l,i) - (D1_B(l,i)/B(l,i))**2) &
                     - two*(lambda(l,i) + r(l,i)*D1_lambda(l,i)))
  end do
!$OMP END PARALLEL DO

! DD_Deltar = d ( Deltar / r )
!              r

  auxarray(l,:) = Deltar(l,:)/r(l,:)
  diffvar => auxarray
  DD_Deltar(l,:) = diff1(l,+1)


! ********************
! ***   Z4 THETA   ***
! ********************

  if (formulation=="z4c") then

!    Derivatives of z4theta.

     diffvar => z4theta
     D1_z4theta(l,:) = diff1(l,+1)

  end if


! *********************************
! ***   ADVECTIVE DERIVATIVES   ***
! *********************************

  if (shift/="none") then

!    Conformal factor.

     diffvar => phi

     if (dorigin=="onesided") then
        DA_phi(l,:) = diffadv(l,+1,.true.)
     else
        DA_phi(l,:) = diffadv(l,+1)
     end if

     diffvar => chi
     DA_chi(l,:) = diffadv(l,+1)

!    A.

     diffvar => A
     DA_A(l,:) = diffadv(l,+1)

!    B.

     diffvar => B
     DA_B(l,:) = diffadv(l,+1)

!    alpha.

     diffvar => alpha
     DA_alpha(l,:) = diffadv(l,+1)

!    trK.

     diffvar => trK
     DA_trK(l,:) = diffadv(l,+1)

!    KTA.

     diffvar => KTA
     DA_KTA(l,:) = diffadv(l,+1)

!    KTB.

     diffvar => KTB
     DA_KTB(l,:) = diffadv(l,+1)

!    lambda.

     diffvar => lambda
     DA_lambda(l,:) = diffadv(l,+1)

!    Klambda.

     diffvar => Klambda
     DA_Klambda(l,:) = diffadv(l,+1)

!    Deltar.

     diffvar => Deltar
     DA_Deltar(l,:) = diffadv(l,-1)

!    z4theta.

     if (formulation=="z4c") then
        diffvar => z4theta
        DA_z4theta(l,:) = diffadv(l,+1)
     end if

!    Puncture regularization.

     if (regular2) then

        diffvar => lambda2
        DA_lambda2(l,:) = diffadv(l,+1)

        diffvar => Klambda2
        DA_Klambda2(l,:) = diffadv(l,+1)

     end if

  end if


! *****************************************
! ***   COVARIANT DERIVATIVES OF LAPSE  ***
! *****************************************

! Second covariant derivative of lapse with
! one index up and one down: D^r D_r alpha.

!$OMP PARALLEL DO SCHEDULE(GUIDED)

  do i=1-ghost,Nrmax

     Dcov2_alpha(l,i) = one/(A(l,i)*psi4(l,i))*(D2_alpha(l,i) &
                      - D1_alpha(l,i)*(half*D1_A(l,i)/A(l,i) + two*D1_phi(l,i)))

!    Laplacian of lapse.

     Lapla_alpha(l,i) = one/(A(l,i)*psi4(l,i))*(D2_alpha(l,i) &
                      - D1_alpha(l,i)*(half*D1_A(l,i)/A(l,i) - D1_B(l,i)/B(l,i) &
                      - two*D1_phi(l,i) - two/r(l,i)))

  end do

!$OMP END PARALLEL DO


! ************************
! ***   RICCI TENSOR   ***
! ************************

!$OMP PARALLEL DO SCHEDULE(GUIDED)

  do i=1-ghost,Nrmax

!    Mixed radial component of physical Ricci tensor R^r_r, found with MAPLE.

     RICA(l,i) = - one/(A(l,i)*psi4(l,i))*(half*D2_A(l,i)/A(l,i) - A(l,i)*D1_Deltar(l,i) &
               - 0.75d0*(D1_A(l,i)/A(l,i))**2 + half*(D1_B(l,i)/B(l,i))**2 - half*Deltar(l,i)*D1_A(l,i) &
               + D1_A(l,i)/r(l,i)/B(l,i) + two*lambda(l,i)*(one + r(l,i)*D1_B(l,i)/B(l,i)) &
               + two*(two*D2_phi(l,i) - D1_phi(l,i)*(D1_A(l,i)/A(l,i) - D1_B(l,i)/B(l,i) - two/r(l,i))))

!    Scalar curvature (trace of Ricci), found with MAPLE.

     RSCAL(l,i) = - one/(A(l,i)*psi4(l,i))*(half*D2_A(l,i)/A(l,i) + D2_B(l,i)/B(l,i) &
                - A(l,i)*D1_Deltar(l,i) - (D1_A(l,i)/A(l,i))**2 + half*(D1_B(l,i)/B(l,i))**2 &
                + two/r(l,i)/B(l,i)*(3.d0 - A(l,i)/B(l,i))*D1_B(l,i) + 4.d0*lambda(l,i))

     RSCAL(l,i) = RSCAL(l,i) - 8.d0/(A(l,i)*psi4(l,i))*(D2_phi(l,i) + D1_phi(l,i)**2 &
                - D1_phi(l,i)*(half*D1_A(l,i)/A(l,i) - D1_B(l,i)/B(l,i) - two/r(l,i)))

!    The term below is equivalent to the second term above, I keep it here
!    for testing.

!    RSCAL(l,i) = RSCAL(l,i) - 8.d0/(A(l,i)*psi(l,i)**5)*(D2_psi(l,i) &
!               - D1_psi(l,i)*(half*D1_A(l,i)/A(l,i) - D1_B(l,i)/B(l,i) - two/r(l,i)))

  end do

!$OMP END PARALLEL DO


! ***************
! ***   END   ***
! ***************

  end subroutine auxiliary_geometry

