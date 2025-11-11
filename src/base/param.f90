
! ****************************************
! ***   PARAMETERS ARE DECLARED HERE   ***
! ****************************************

  module param

!               VERY IMPORTANT  (PLEASE READ)
!
! Declare only one parameter per line.  Even if FORTRAN allows
! one to declare several parameters separated by commas, this
! file will be processed by a PERL script that does NOT understand
! that.
!
! Use always the format:
!
!    type :: name = value
!
! The following variations are permitted:  Character-type parameters
! that are allowed to receive multiple values at the same time
! (separated by commas) should be declared as:
!
!    character :: name = value     ! multiple
!
! Also, a range can be defined for character-type parameters as:
!
!    character :: name = value     ! range = (value1,value2,...,valuen)
!
! The range is not compulsory.  If it is not there, any value
! is allowed (for example, directory names).
!
! REMEMBER:  All parameters must be initialized.  The initial
! value should basically correspond to the code NOT doing
! anything weird.  That is, initialize to Minkowski, static
! slicing, no shift, vacuum, and all special features turned off.


! **************************
! ***   PARAMETER FILE   ***
! **************************

! The parameter file is piped to the code at run time,
! so it is not really itself a "parameter", but it is
! convenient to declare it here since its name will
! be known to all other routines (particularly) the
! checkpoint routine.

  character(1000) :: parfile = ""


! ****************
! ***   GRID   ***
! ****************

! dr0:           Spatial interval for base (coarsest) grid.
! rbound:        Position of outer boundary.
!
! Nrtotal:       Total number of grid points.
! Nr:            Local number of grid points.
! Nl:            Number of refinement levels.
! ghost:         Number of ghost zones.
! intorder:      Order of interpolation (1,3).
!
! nproctot:      Total number of processors used in the run.  It is not really a parameter
!                and I use it only for checkpointing (should not be set in parameter file).
!
! Notice that the number of ghost zones should
! not be given in the parameter file since it
! is controlled by the parameter "order".
!
! Also, the local number of grid points Nr should
! not be given either in the parameter file.

  real(8) :: dr0  = 1.d0
  real(8) :: rbound = 0.0    ! Should not be set in parameter file!

  integer :: Nrtotal = 10
  integer :: Nr = 1
  integer :: Nl = 1

  integer :: ghost = 0
  integer :: intorder = 3

  integer :: nproctot = 1    ! Should not be set in parameter file!


! *************************
! ***   TIME STEPPING   ***
! *************************

! dt0:           Time step for base (coarsest) grid.
!
! adjuststep:    Do we ajust the time step using the CFL condition?
! adjuststepmax  Do we allow the time step to increase from its original value?
!
! dtfac:         Courant parameter (dtfac=dt/dr).
! Nt:            Total number of time steps.

  logical :: adjuststep = .false.
  logical :: adjuststepmax = .true.

  real(8) :: dt0 = 0.d0
  real(8) :: dtfac = 0.5d0

  integer :: Nt = 10


! *********************************
! ***   OUTPUT AND CHECKPOINT   ***
! *********************************

! directory:       Name of directory for output.
! checkpointfile:  Name of checkpoint directory.
!
! An important point to remember is that, because of the
! way Fortran works, this string has to be defined with a
! fixed length.  This means that before using it to identify
! a directory, it must first be trimmed to its true size.

  character(100) :: directory = "output"
  character(100) :: checkpointfile = "checkpoint"

! Output parameters.
!
! checkpoint:    Do we output checkpoint files?
! checkpointinitial:  Do we checkpoint initial data?
!
! closefiles:    Do we close output files each time we write to them?
! Ninfo:         How often do we output information to screen?
! Noutput0D:     How often do we do 0D output?
! Noutput1D:     How often do we do 1D output?
!
! outvars0D:     Variables that need 0D output (a list separated by commas).
! outvars1D:     Variables that need 1D output (a list separated by commas).
!
! checkvars:     Variables that need to be output for a checkpoint (should not be set in parameter file).
!                It is initialized with the coordinate array (r).  The evolving arrays are later
!                added to the list at run time when they are allocated.
!
! nvars0D:       Number of variables with 0D output (should not be set in parameter file).
! nvars1D:       Number of variables with 1D output (should not be set in parameter file).
!
! commenttype:   Type of comment lines on files (xgraph,gnuplot).
! norm_rmin:     Minimum value of r for calculating norms.
! output_r1:     Point for 0D output.

  logical :: checkpointinitial = .true.
  logical :: checkpoint = .false.

  logical :: closefiles = .true.

  integer :: Ninfo = 1
  integer :: Noutput0D = 1
  integer :: Noutput1D = 1
  integer :: Ncheckpoint = 1000
  integer :: nvars0D = 0     ! Not a real parameter, DO NOT set in parameter file.
  integer :: nvars1D = 0     ! Not a real parameter, DO NOT set in parameter file.

  real(8) :: norm_rmin = 0.d0
  real(8) :: output_r1 = 0.d0

  character(1000) :: outvars0D = "alpha"      ! multiple
  character(1000) :: outvars1D = "alpha"      ! multiple
  character(1000) :: checkvars = "r"          ! multiple

  character(1000) :: commenttype = "gnuplot"  ! range = (xgraph,gnuplot)


! **********************
! ***   BOUNDARIES   ***
! **********************

! boundary:      Type of boundary condition.
! dorigin:       Derivatives at origin: centered or one-sided

  character(30) :: boundtype = "radiative"    ! range = (none,static,flat,radiative,constraint)
  character(30) :: dorigin = "centered"       ! range = (centered,onesided)


! *******************
! ***   SLICING   ***
! *******************

! slicing:       Type of slicing condition.
!
! ilapse:        Type of initial lapse:
!                   none (default):  This parameter is ignored and we can set up the initial lapse elsewhere.
!                   one:             Initial lapse equal to 1.
!                   isotropic:       Isotropic lapse for Schwarzschild and Reissner-Norstrom data.
!                   psiminus2:       alpha = 1/psi**2
!                   psiminus4:       alpha = 1/psi**4
!                   psi2:            alpha = psi**2
!                   psi4:            alpha = psi**4
!                   maximal:         Initialize lapse to maximal slicing (assuming trK=0).
!
! maximalbound:  Boundary condition for maximal slicing.
! gauge_f:       Coefficient for Bona-Masso type slicing condition (positive).
!
! lapsepert:     Perturbation to initial lapse.
! lapse_a0:      Amplitude of initial profile.
! lapse_r0:      Center of initial profile.
! lapse_s0:      Width of initial profile.
! lapse_t0:      Transition with of initial profile (for top-hat profile).
!
! lapseeta:      Coefficient of damping term for second order Bona-Masso slicing conditions.
!
! lapsediss:     Coefficient of real dissipation for evolution-type slicing
!                conditions.  This is NOT numerical dissipation, but rather
!                corresponds to adding a finite parabolic term to the Bona-Masso
!                condition in order to avoid gauge shocks.

  character(30) :: slicing = "harmonic"    ! range = (static,maximal,harmonic,1+log,shockavoid,alphaminus2,cosmocf-harmonic,cosmocf-1+log,cosmocf-shockavoid,cosmosync-harmonic,cosmosync-1+log,cosmosync-shockavoid)

  character(30) :: ilapse  = "none"        ! range = (none,one,isotropic,psiminus2,psiminus4,psi2,psi4,maximal)
  character(30) :: maximalbound = "robin"  ! range = (robin,dirichlet,conformal)

  real(8) :: gauge_f = 1.d0

  character(30) :: lapsepert = "none"      ! range = (none,gaussian,tophat)
  real(8) :: lapse_a0 = 0.d0
  real(8) :: lapse_r0 = 0.d0
  real(8) :: lapse_s0 = 1.d0
  real(8) :: lapse_t0 = 1.d0

  real(8) :: lapseeta  = 0.d0
  real(8) :: lapsediss = 0.d0


! *****************
! ***   SHIFT   ***
! *****************

! shift:         Type of shift condition.
!
! driverD0:      Do we add the term Delta0 to driver equations?
! drivercsi:     Coefficient of Deltar for Gammadriver.
! drivereta:     Coefficient of damping term for Gammadriver.
!
! shiftpert:     Perturbation to initial shift.
! shift_a0:      Amplitude of initial profile.
! shift_r0:      Center of initial profile.
! shift_s0:      Width of initial profile.
! shift_t0:      Transition width of initial profile (for top-hat profile).

  character(30) :: shift = "none"          ! range = (none,zero,static,Gammadriver0,Gammadriver1,Gammadriver2,Gammadriver3,Gammadrivershock1,Gammadrivershock2)

  logical :: driverD0 = .false.

  real(8) :: drivercsi = 0.75d0
  real(8) :: drivereta = 0.d0

  character(30) :: shiftpert = "none"      ! range = (none,gaussian,tophat)
  real(8) :: shift_a0 = 0.d0
  real(8) :: shift_r0 = 0.d0
  real(8) :: shift_s0 = 1.d0
  real(8) :: shift_t0 = 1.d0


! ************************
! ***   INITIAL DATA   ***
! ************************

! idata:         Type of initial data.
!
! rescaledata:   Rescale initial profile by the factor. Useful if we want to
!                rescale the initial profiles by this factor.

  real(8) :: rescaledata = 1.d0

  character(30) :: idata = "minkowski"  ! range = (checkpoint,minkowski,schwarzschild,schwarzschildKS,trumpetBH,reissnernordstrom,desitter,scalarpulse,ghostpulse,nonminpulse,complexpulse,complexghostpulse,bosonstar,chargedboson,procapulse,procastar,l-procastar,chargedproca,diracpulse,diracstar,dustshell,fluidshell,TOVstar,TOVcomplex,blastwave,scalarDM,complexDM,ghostwormhole,duststep)


! *********************
! ***   EVOLUTION   ***
! *********************

! spacetime:     Dynamic spacetime, background (static) spacetime, or just Minkowski
!
! integrator:    Time integration method.
! order:         Order of spatial differencing.
! icniter:       Number of iterations for icn.
!
! formulation:   Formulation of evolution equations (bssn,z4c).
! bssnflavor:    Evolution of volume element in BSSN (Eulerian/Lagrangian).
! eta:           Multiple of momentum constraint in source of Deltar.
!
! chimethod:     Use chi=1/psi**n=exp(-n*phi) in evolutions instead of phi (for BH evolutions).
! chipower:      Power of 1/psi in the definition of chi.
!
! nolambda:      Do not use regularization variables.
! noDeltar:      Do not use BSSN Delta variable and calculate it from derivatives of the metric.
! noKTA:         Do not use KTA as an independent variable and calculate it from Klambda.
! regular2:      Puncture regularization based on lambda2.
! lambdapower:   Power of psi in the definition of lambda2.
!
! kappa1:        Constraint dissipation coefficient for Z4c.
! kappa2:        Constraint dissipation coefficient for Z4c.

  character(30) :: spacetime   = "dynamic"    ! range = (dynamic,background,minkowski)
  character(30) :: formulation = "bssn"       ! range = (bssn,z4c)
  character(30) :: bssnflavor  = "lagrangian" ! range = (eulerian,lagrangian)
  character(30) :: integrator  = "rk4"        ! range = (icn,rk4)
  character(30) :: order       = "four"       ! range = (two,four,six,eight)

  logical :: nolambda  = .false.              ! By default we evolve the regularization variables.
  logical :: noDeltar  = .false.              ! By default we evolve the variable Deltar.
  logical :: noKTA     = .true.               ! By default we obtain KTA from Klambda.
  logical :: chimethod = .true.               ! By default we evolve chi instead of phi.
  logical :: regular2  = .true.               ! By default we evolve Klambda2 instead of Klamba.

  integer :: icniter = 3

  integer :: chipower = 4
  integer :: lambdapower = 4

  real(8) :: eta = 2.0d0

  real(8) :: kappa1 = 1.d0
  real(8) :: kappa2 = 1.d0

! Kreiss-Oliger dissipation coefficients.  It is better if they are independent.

  real(8) :: geodiss    = 0.01d0      ! Dissipation for geometric variables.
  real(8) :: scalardiss = 0.01d0      ! Dissipation for scalar field.
  real(8) :: nonmindiss = 0.01d0      ! Dissipation for nonmin field.
  real(8) :: elecdiss   = 0.01d0      ! Dissipation for electric field.
  real(8) :: procadiss  = 0.01d0      ! Dissipation for Proca field.
  real(8) :: diracdiss  = 0.01d0      ! Dissipation for Dirac field.
  real(8) :: fluiddiss  = 0.0d0       ! Dissipation for fluids (better zero).


! ********************
! ***   HORIZONS   ***
! ********************

! ahfind:        Do we want to find an apparent horizon? 
! ahfind_every:  How often do we look for horizons?
! ahafter:       Start looking for horizons after a given time.

  logical :: ahfound = .false.
  logical :: ahfind = .false.

  integer :: ahfind_every = 1

  real(8) :: ahafter = 0.d0
  real(8) :: ahmass = 0.d0


! ******************
! ***   MATTER   ***
! ******************

! IMPORTANT:  The different matter types must have very different names,
!             and one name must NOT be contained in another.  That is,
!             don't use names such as "scalar" and "scalar-ghost" for
!             different matter types.  The reason for this is that the
!             code just checks for a given string to be contained in
!             "mattertype", and in such a case it will be contained
!             in both cases.
!
! mattertype (type of matter):
!
!    vacuum       =  No matter.
!    cosmo        =  Cosmological constant.
!    scalar       =  Real scalar field.
!    ghost        =  Real ghost scalar field.
!    complex      =  Complex scalar field.
!    complexghost =  Complex ghost scalar field.
!    nonmin       =  Non-minimally coupled scalar field.
!    electric     =  Electric field (there is no magnetic field in spherical symmetry).
!    proca        =  Proca field (complex massive EM field).
!    complexproca =  Complex Proca field.
!    dirac        =  Dirac field.
!    fluid        =  Perfect fluid type matter.
!    dust         =  Dust (fluid with zero pressure).

  character(1000) :: mattertype = "vacuum"  ! multiple, range=(vacuum,cosmo,scalar,ghost,complex,complexghost,nonmin,electric,proca,complexproca,dirac,fluid,dust)


! *********************
! ***   COSMOLOGY   ***
! *********************

! cosmic_run:    Flag for assigning storage for cosmological variables.

  logical :: cosmic_run = .false.

! lambda_cosmo:  Cosmological constant.

  real(8) :: lambda_cosmo = 0.d0


! *********************************
! ***   BLACK HOLE SPACETIMES   ***
! *********************************

! BHmass:        Black hole mass.
! BHcharge:      Black hole charge.

  real(8) :: BHmass = 1.d0
  real(8) :: BHcharge = 0.d0


! *****************************
! ***   REAL SCALAR FIELD   ***
! *****************************

! scalarprofile:     Form of initial profile.
! scalar_a0:         Amplitude of initial profile.
! scalar_r0:         Center of initial profile.
! scalar_s0:         Width of initial profile.
! scalar_t0:         Transition width of initial profile (for top-hat profile).
!
! scalarpotential:   Type of scalar field potential.
! scalar_mass:       Scalar field mass parameter.
! scalar_lambda:     Coefficient of phi**4 term in potential.
!
! scalarmethod:      Finite difference method.
! scalar_relax:      Do we use relaxation method for initial data?

  real(8) :: scalar_a0 = 0.d0
  real(8) :: scalar_r0 = 0.d0
  real(8) :: scalar_s0 = 1.d0
  real(8) :: scalar_t0 = 1.d0

  real(8) :: scalar_mass = 0.d0
  real(8) :: scalar_lambda = 0.d0

  character(1000) :: scalarprofile   = "gaussian"   ! range=(gaussian,r2gaussian,tophat)
  character(1000) :: scalarpotential = "none"       ! range=(none,phi2,phi4)
  character(1000) :: scalarmethod    = "second"     ! range=(first,second)

  logical :: scalar_relax = .false.

! Cosmological runs:
!
! scalar_bg_phi0:    Amplitude of background scalar field.
! scalar_bg_pi0:     Amplitude of background time derivative of scalar field.
!
! scalar_bg_pert:    Do we add a perturbation to the background field?

  logical :: scalar_bg_pert = .false.

  real(8) :: scalar_bg_phi0 = 0.d0
  real(8) :: scalar_bg_pi0  = 0.d0


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

! complexprofile:    Form of initial profile.
!
! complexR_a0:       Amplitude of initial profile.
! complexR_r0:       Center of initial profile.
! complexR_s0:       Width of initial profile.
! complexR_t0:       Transition width of initial profile (for top-hat profile).
! complexR_a1:       Amplitude of second profile for compensated cases.
! complexR_r1:       Center of second profile for compensated cases.
! complexR_s1:       Width of second profile for compensated cases.
! complexR_t1        Transition width of second profile for compensated cases.
!
! complexI_a0:       Amplitude of initial profile.
! complexI_r0:       Center of initial profile.
! complexI_s0:       Width of initial profile.
! complexI_t0:       Transition width of initial profile (for top-hat profile).
! complexI_a1:       Amplitude of second profile for compensated cases.
! complexI_r1:       Center of second profile for compensated cases.
! complexI_s1:       Width of second profile for compensated cases.
! complexI_t1        Transition width of second profile for compensated cases.
!
! complexpotential:  Type of complex field potential.
! complex_mass:      Complex field mass parameter.
! complex_lambda:    Coefficient of phi**4 term in potential.
!
! complex_l:         Total angular momentum for "rotating" case (azimuthal quantum number).
! complex_lmax:      Maximum value of l for the case of multiple l's.
!
! complexmethod:     Finite difference method.
!
! complexDM_type:    Type of initial data for complex dark matter. Allowed types are:
!                       harmonicTS:    Time symmetric initial data with harmonic time dependence.
!                       harmonicMOM:   Harmonic time dependence but with with non-zero momentum density.
!                       growing:       Purely growing mode.

  real(8) :: complexR_a0 = 0.d0
  real(8) :: complexR_r0 = 0.d0
  real(8) :: complexR_s0 = 1.d0
  real(8) :: complexR_t0 = 1.d0

  real(8) :: complexR_a1 = 0.d0
  real(8) :: complexR_r1 = 0.d0
  real(8) :: complexR_s1 = 1.d0
  real(8) :: complexR_t1 = 1.d0

  real(8) :: complexI_a0 = 0.d0
  real(8) :: complexI_r0 = 0.d0
  real(8) :: complexI_s0 = 1.d0
  real(8) :: complexI_t0 = 1.d0

  real(8) :: complexI_a1 = 0.d0
  real(8) :: complexI_r1 = 0.d0
  real(8) :: complexI_s1 = 1.d0
  real(8) :: complexI_t1 = 1.d0

  real(8) :: complex_mass = 0.d0
  real(8) :: complex_lambda = 0.d0

  integer :: complex_l    = 0
  integer :: complex_lmax = 0

  real(8) :: k_parameter = 0.d0

  character(1000) :: complexprofile   = "gaussian"    ! range=(gaussian,tophat,comp-gaussian,comp-tophat)
  character(1000) :: complexpotential = "none"        ! range=(none,phi2,phi4)
  character(1000) :: complexmethod    = "second"      ! range=(first,second)
  character(1000) :: complexDM_type   = "harmonicTS"  ! range=(harmonicTS,harmonicMOM,growing)

! Boson star initial data:
!
! boson_phi0:        Boson star parameter for value at origin (see boson_factor).
! boson_factor:      Normalization factor:
!                         * physical:  phi(r<<1) ~ phi0 r**l
!                         * harmonic:  phi(r<<1) ~ phi0/sqrt(4pi*(2l+1)) r**l
! boson_omega:       Frequency of the boson star.
! omega_left:        Left  frequency guess.
! omega_right:       Right frequency guess.
!
! boson_gauge:       Gauge for boson star initial data (PA=polar-areal,CF=conformally flat)   
! boson_relax:       Do we use relaxation method for boson star initial data?

  real(8) :: boson_phi0  = 0.d0

  real(8) :: boson_omega = 0.d0
  real(8) :: omega_left    = 0.d0
  real(8) :: omega_right   = 0.d0

  character(1000) :: boson_factor = "harmonic" ! range=(physical,harmonic)
  character(1000) :: boson_gauge  = "PA"       ! range=(PA,CF)

  logical :: boson_relax = .false.

! Boson star perturbation.
!
! bosongauss:        Do we add an initial profile to the boson star solution?
! boson_a0:          Amplitude of initial profile.
! boson_r0:          Center of initial profile.
! boson_s0:          Width of initial profile.

  logical :: bosongauss = .false.

  real(8) :: boson_phiR_a0 = 0.d0
  real(8) :: boson_phiR_r0 = 0.d0
  real(8) :: boson_phiR_s0 = 1.d0
  real(8) :: boson_piI_a0  = 0.d0
  real(8) :: boson_piI_r0  = 0.d0
  real(8) :: boson_piI_s0  = 1.d0

! Cosmological runs:
!
! complex_phiR0:    Amplitude of the real part of background.
! complex_phiI0:    Amplitude of the imaginary part of background.
!
! complex_bg_pert:  Add a perturbation to the background complex scalar field?

  logical :: complex_bg_pert = .false.

  real(8) :: complex_bg_phiR0 = 0.d0
  real(8) :: complex_bg_phiI0 = 0.d0

! Charged complex scalar field:
!
! complex_initialcharged:   Charged initial data?
! complex_q:                Charge of complex scalar field.
! charge_factor:            Charge normalization for boson stars.

  logical :: complex_initialcharged = .false.

  real(8) :: complex_q = 0.d0
  character(1000) :: charge_factor = "standard" ! range=(standard,jetzer)


! ******************************
! ***   GHOST SCALAR FIELD   ***
! ******************************

! ghostprofile:      Form of initial profile.
! ghost_a0:          Amplitude of initial profile.
! ghost_r0:          Center of initial profile.
! ghost_s0:          Width of initial profile.
! ghost_t0:          Transition width of initial profile (for top-hat profile).
!
! ghostpotential:    Type of ghost field potential.
! ghost_mass:        Complex field mass parameter.
! ghost_lambda:      Coefficient of phi**4 term in potential.
!
! ghostmethod:       Finite difference method.

  real(8) :: ghost_a0 = 0.d0
  real(8) :: ghost_r0 = 0.d0
  real(8) :: ghost_s0 = 1.d0
  real(8) :: ghost_t0 = 1.d0

  real(8) :: ghost_mass = 0.d0
  real(8) :: ghost_lambda = 0.d0

  character(1000) :: ghostprofile   = "gaussian"   ! range=(gaussian,tophat)
  character(1000) :: ghostpotential = "none"       ! range=(none,phi2,phi4)
  character(1000) :: ghostmethod    = "second"     ! range=(first,second)


! ***************************************
! ***   COMPLEX GHOST SCALAR FIELD    ***
! ***************************************

! complexghostprofile:    Form of initial profile.
!
! complexghostR_a0:       Amplitude of initial profile.
! complexghostR_r0:       Center of initial profile.
! complexghostR_s0:       Width of initial profile.
! complexghostR_t0:       Transition width of initial profile (for top-hat profile).
!
! complexghostI_a0:       Amplitude of initial profile.
! complexghostI_r0:       Center of initial profile.
! complexghostI_s0:       Width of initial profile.
! complexghostI_t0:       Transition width of initial profile (for top-hat profile).
!
! complexghostpotential:  Type of complex field potential.
! complexghost_mass:      Complex field mass parameter.
! complexghost_lambda:    Coefficient of phi**4 term in potential.
! complexghost_l:         Total angular momentum for "rotating" case (azimuthal quantum number).
! complexghost_lmax:      Maximum valur of l for the case of multiple l's.
!
! complexghostmethod:     Finite difference method.

  real(8) :: complexghostR_a0 = 0.d0
  real(8) :: complexghostR_r0 = 0.d0
  real(8) :: complexghostR_s0 = 1.d0
  real(8) :: complexghostR_t0 = 1.d0

  real(8) :: complexghostI_a0 = 0.d0
  real(8) :: complexghostI_r0 = 0.d0
  real(8) :: complexghostI_s0 = 1.d0
  real(8) :: complexghostI_t0 = 1.d0

  real(8) :: complexghost_mass = 0.d0
  real(8) :: complexghost_lambda = 0.d0

  integer :: complexghost_l    = 0
  integer :: complexghost_lmax = 0

  character(1000) :: complexghostprofile   = "gaussian" ! range=(gaussian,tophat)
  character(1000) :: complexghostpotential = "none"     ! range=(none,phi2,phi4)
  character(1000) :: complexghostmethod    = "second"   ! range=(first,second)


! **********************************************
! ***   NON-MINIMALLY COUPLED SCALAR FIELD   ***
! **********************************************

! nonminprofile:     Form of initial profile.
! nonmin_a0:         Amplitude of initial profile.
! nonmin_r0:         Center of initial profile.
! nonmin_s0:         Width of initial profile.
! nonmin_t0:         Transition width of initial profile (for top-hat profile).
!
! nonmin_fxi:        Coefficient of quadratic term in nonminimal coupling.
! nonmin_theta:      Parameter for slicing condition.
! nonminf:           Form of coupling function f(phi).
! nonmin_phi0:       For quadratic coupling we have:  f = 1/(8 pi)  +  xi*(phi-phi0)^2
! nonminpotential:   Type of potential for nonmin field.
!
! nonminmethod:      Finite difference method.
! nonmin_pulse:      Pulse shape.

  real(8) :: nonmin_a0 = 0.d0
  real(8) :: nonmin_r0 = 0.d0
  real(8) :: nonmin_s0 = 1.d0
  real(8) :: nonmin_t0 = 1.d0

  real(8) :: nonmin_fxi   = 1.d0
  real(8) :: nonmin_theta = 1.d0
  real(8) :: nonmin_phi0  = 0.d0

  character(1000) :: nonminf         = "trivial"  ! range=(trivial,quadratic,BransDicke)
  character(1000) :: nonminprofile   = "gaussian" ! range=(gaussian,r2gaussian,tophat)
  character(1000) :: nonminpotential = "none"     ! range=(none)
  character(1000) :: nonminmethod    = "second"   ! range=(first,second)
  character(1000) :: nonmin_pulse    = "gaussian" ! range=(gaussian,tanh)


! ***********************
! ***   PROCA FIELD   ***
! ***********************

! proca_mass:        Mass parameter for Proca field.
!
! procaprofile:      Form of initial profile.
! proca_a0:          Amplitude of initial profile.
! proca_r0:          Center of initial profile.
! proca_s0:          Width of initial profile.
! proca_t0:          Transition width of initial profile (for top-hat profile).

  real(8) :: proca_mass = 1.d0

  real(8) :: proca_a0 = 0.d0
  real(8) :: proca_r0 = 0.d0
  real(8) :: proca_s0 = 1.d0
  real(8) :: proca_t0 = 1.d0

  character(1000) :: procaprofile = "gaussian"    ! range=(gaussian,tophat)


! *******************************
! ***   COMPLEX PROCA FIELD   ***
! *******************************

! cproca_mass:       Mass parameter for complex Proca field.
! complex_q:         Charge of complex Proca field.
! cproca_l:          Total angular momentum for "rotating" case (azimuthal quantum number).

  real(8) :: cproca_mass = 1.d0
  real(8) :: cproca_q = 0.d0

  integer :: cproca_l    = 0

! Proca star initial data (we use the same parameters omega_right
! and omega_left declared above for boson stars):
!
! proca_phi0:        Proca star parameter for value at origin (see proca_factor).
! proca_factor:      Normalization factor:
!                         * physical:  phi(r<<1) ~ phi0 r**l
!                         * harmonic:  phi(r<<1) ~ phi0/sqrt(4pi*(2l+1)) r**l
! proca_omega:       Frequency of the boson star.

  real(8) :: proca_phi0  = 0.d0
  real(8) :: proca_omega = 0.d0

  character(1000) :: proca_factor = "harmonic"    ! range=(physical,harmonic)

! Proca star perturbation.
!
! procagauss:        Do we add a initial profile to the proca star solution?
! proca_a0:          Amplitude of initial profile.
! proca_r0:          Center of initial profile.
! proca_s0:          Width of initial profile.

  logical :: procagauss = .false.

  real(8) :: proca_PhiR_a0 = 0.d0
  real(8) :: proca_PhiR_r0 = 0.d0
  real(8) :: proca_PhiR_s0 = 1.d0

  real(8) :: proca_AI_a0 = 0.d0
  real(8) :: proca_AI_r0 = 1.d0
  real(8) :: proca_AI_s0 = 1.d0


! ***********************
! ***   DIRAC FIELD   ***
! ***********************

! dirac_mass:        Mass parameter for Dirac field.

  real(8) :: dirac_mass = 1.d0

! diracprofile:      Form of initial profile.
!
! dirac*_a0:         Amplitude of initial profile.
! dirac*_r0:         Center of initial profile.
! dirac*_s0:         Width of initial profile.
! dirac*_t0:         Transition width of initial profile (for top-hat profile).
!
! diractype:         Type of initial data (to guarantee time symmetry).
!                        1  Both F and G must be purely real or purely imaginary.
!                        2  FI = k1 FR and GI = k2 GR, with (k1,k2) constants.
!                        3  GR = +- FR/k,  GI = -+ k FI, or GR = +- FI/k,  GI = -+ k FR, with k constant.
! dirac_k:           Constant por type 3 initial data.

  character(1000) :: diracprofile   = "gaussian"   ! range=(gaussian,tophat)

  integer :: diractype = 1

  real(8) :: dirac_k = 1.d0

! Profile parameters for (GR,GI).

  real(8) :: diracFR_a0 = 0.d0
  real(8) :: diracFR_r0 = 0.d0
  real(8) :: diracFR_s0 = 1.d0
  real(8) :: diracFR_t0 = 1.d0

  real(8) :: diracFI_a0 = 0.d0
  real(8) :: diracFI_r0 = 0.d0
  real(8) :: diracFI_s0 = 1.d0
  real(8) :: diracFI_t0 = 1.d0

! Profile parameters for (GR,GI). Notice that diracGR_r0 and diracGI_r0
! must both be non-zero since G must be odd.

  real(8) :: diracGR_a0 = 0.d0
  real(8) :: diracGR_r0 = 1.d0
  real(8) :: diracGR_s0 = 1.d0
  real(8) :: diracGR_t0 = 1.d0

  real(8) :: diracGI_a0 = 0.d0
  real(8) :: diracGI_r0 = 1.d0
  real(8) :: diracGI_s0 = 1.d0
  real(8) :: diracGI_t0 = 1.d0

! Dirac star initial data (we use the same parameters omega_right
! and omega_left declared above for boson stars):
!
! dirac_f0:          Dirac star parameter for value of function f at origin.
! dirac_omega:       Frequency of the Dirac star.

  real(8) :: dirac_f0    = 0.d0
  real(8) :: dirac_omega = 0.d0

! Dirac star perturbation.
!
! diracgauss:        Do we add a initial profile to the boson star solution?

  logical :: diracgauss = .false.


! ****************
! ***   DUST   ***
! ****************

! dustprofile:       Form of initial profile.
! dust_a0:           Amplitude of initial profile.
! dust_r0:           Center of initial profile.
! dust_s0:           Width of initial profile.
! dust_t0:           Transition width of initial profile (for top-hat profile).
!
! dust_atmos:        Size of artificial atmosphere for dust.
!
! dust_method:       Method for dust integration.

  real(8) :: dust_a0 = 0.d0
  real(8) :: dust_r0 = 0.d0
  real(8) :: dust_s0 = 1.d0
  real(8) :: dust_t0 = 1.d0

  real(8) :: dust_atmos = 1.d-6 ! This is rescaled with the maximum density at t=0

  character(1000) :: dustprofile = "gaussian"  ! range=(gaussian,tophat)
  character(1000) :: dust_method = "limiter"   ! range=(center,upwind,limiter,mp5)


! *****************
! ***   FLUID   ***
! *****************

! fluidprofile:      Form of initial profile.
! fluid_a0:          Amplitude of initial profile.
! fluid_r0:          Center of initial profile.
! fluid_s0:          Width of initial profile.
! fluid_t0:          Transition width of initial profile (for top-hat profile).
!
! fluid_EOS:         Fluid's equation of state.
!
! fluid_N:           Polytropic index.
! fluid_gamma:       Adiabatic index.
! fluid_kappa:       Polytropic constant.
!
! fluid_atmos:       Size of artificial atmosphere.
!
! fluid_q1:          Coefficient of linear artificial viscosity.
! fluid_q2:          Coefficient of quadratic artificial viscosity.
!
! fluid_method:      Method for fluid integration.
! fluid_limiter:     Type of limiter for reconstruction.
!
! Notice that one should have gamma = 1 + 1/N.  The default values
! use here are compatible with that:  N=3/2, gamma=5/3.
!
! The parameter N is NEVER used in the evolutions, just in initial
! data.  However, one should be very careful: I am not quite sure if
! it makes sense to use values of N in the initial data that are not
! consistent with the value of gamma used in the evolution.

  real(8) :: fluid_a0 = 0.d0
  real(8) :: fluid_r0 = 0.d0
  real(8) :: fluid_s0 = 1.d0
  real(8) :: fluid_t0 = 1.d0

  real(8) :: fluid_N = 1.5d0                    ! 3/2, consistent with gamma=5/3.
  real(8) :: fluid_gamma = 1.6666666666666666d0 ! 5/3 for a monoatomic non-relativistic gas.
  real(8) :: fluid_kappa = 1.d0


! Don't set the artificial atmosphere to less than 1.d-8 since
! this is the level of roundoff error in fluidprimitive.f90
! (it involves square roots).

  real(8) :: fluid_atmos = 1.d-8

! Set artificial viscosity to zero by default.

  real(8) :: fluid_q1 = 0.d0
  real(8) :: fluid_q2 = 0.d0

  character(1000) :: fluidprofile  = "gaussian" ! range=(gaussian,tophat)
  character(1000) :: fluid_EOS     = "ideal"    ! range=(none,ideal)
  character(1000) :: fluid_method  = "hlle"     ! range=(llf,hlle,mp5)
  character(1000) :: fluid_limiter = "mc"       ! range=(minmod,vanleer,superbee,mc,koren,ospre,sweby)

! Do we use the speed of sound for calculating 
! fluxed in the different fluid methods?
! Using the sound speed is more accurate, but for
! cases when the density changes by many orders
! of magnitude (like stars) it is less stable.
! If the flag is set to false we use the speed
! of light instead to calculate fluxes

  logical :: fluid_usesoundspeed = .false.

! TOV initial data:
!
! TOV_rho0:          Central density of the star.
! TOV_rad:           Radius of TOV star. This is not really a free parameter
!                    and should not be fixed in the parameter file.

  real(8) :: TOV_rho0 = 0.001d0
  real(8) :: TOV_rad  = 0.d0

! TOV star perturbation.
!
! TOV_gauss:         Do we add an initial profile to the TOV star solution?
! TOV_a0:            Amplitude of initial profile.
! TOV_r0:            Center of initial profile.
! TOV_s0:            Width of initial profile.

  logical :: TOVgauss = .false.

  real(8) :: TOV_a0 = 0.d0
  real(8) :: TOV_r0 = 0.d0
  real(8) :: TOV_s0 = 1.d0

! Blast wave initial data:
!
! blast_R            Initial radius of inner region.
!
! blast_rhol         Density to the left of R.
! blast_rhor         Density to the right of R.
!
! blast_pl           Density to the left of R.
! blast_pr           Density to the right of R.

  real(8) :: blast_R = 1.d0

  real(8) :: blast_rhol = 0.d0
  real(8) :: blast_rhor = 0.d0

  real(8) :: blast_pl = 0.d0
  real(8) :: blast_pr = 0.d0


! ******************************************
! ***   COMPACTIFIED RADIAL COORDINATE   ***
! ******************************************

! newr:              Do we use a transformed radial coordinate?
!
! rtype:             Type of radial coordinate transformation.
!
! drinf:             Derivative of new radius at infinity.
! beta_transf:       Beta parameter for radial transformation.
! gamma_transf:      Gamma parameter for radial transformation.
! delta_transf:      Delta parameter for radial transformation.
!
! rtmin:             r_trans position of first grid point (on finer grid).
! rtmax:             r_trans position of last grid point (on coarsest grid).
!
! integralmethod     Method used to integrate the transformation equation.
! Ncheb              Number of Chebyshev polynomials in expansion.

  logical :: newr = .false.

  character(1000) :: rtype = "sigmoid" ! range=(sigmoid,choptuik,smoothstep,polynomal,sinh)

  real(8) :: drinf        = 1.d0
  real(8) :: beta_transf  = -1.d0
  real(8) :: gamma_transf = 1.d0
  real(8) :: delta_transf = 3.d0

! Integration method.

  character(1000) :: integralmethod = "quadrature" ! range=(quadrature,fourth)

  integer :: NCheb = 10

  real(8) :: rtmin = 0.d0    ! Not a real parameter, DO NOT set in parameter file.
  real(8) :: rtmax = 0.d0    ! Not a real parameter, DO NOT set in parameter file.


! ****************************************
! ***   TRACKING ANALYTIC SPACETIMES   ***
! ****************************************

! Minkowski.

  logical :: TrackMinkowski = .false.

! Schwarzschild.

  logical :: TrackSchwarzschild = .false.


! ***************
! ***   END   ***
! ***************

  end module param

