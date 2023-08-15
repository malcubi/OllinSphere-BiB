!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/potential.f90,v 1.20 2023/03/02 18:43:10 malcubi Exp $

  subroutine potential(l)

! *****************************************************
! ***   POTENTIAL FOR THE DIFFERENT SCALAR FIELDS   ***
! *****************************************************

! This routine calculates the self interaction potential and its derivative
! for the different types of scalar field: scalar,ghost,complex.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  logical contains

  integer l


! ************************
! ***   SCALAR FIELD   ***
! ************************

  if (contains(mattertype,"scalar")) then

!    **************************
!    ***   ZERO POTENTIAL   ***
!    **************************

!    By default, set to zero.

     if (scalarpotential=="none") then

        scalar_V(l,:)  = 0.d0
        scalar_VP(l,:) = 0.d0

        if (cosmic_run) then
           cosmobg_scalar_V(l)  = 0.d0
           cosmobg_scalar_VP(l) = 0.d0
        end if


!    ***************************
!    ***   PHI^2 POTENTIAL   ***
!    ***************************

!    A quadratic potential indicates a mass term
!    for the scalar field and has the form:
!
!    V  =  m^2 phi^2 / 2

     else if (scalarpotential=="phi2") then

        scalar_V(l,:)  = 0.5d0*scalar_mass**2*scalar_phi(l,:)**2
        scalar_VP(l,:) = scalar_mass**2*scalar_phi(l,:)

        if (cosmic_run) then
           cosmobg_scalar_V(l)  = 0.5d0*scalar_mass**2*cosmobg_scalar_phi(l)**2
           cosmobg_scalar_VP(l) = scalar_mass**2*cosmobg_scalar_phi(l)
        end if


!    ***************************
!    ***   PHI^4 POTENTIAL   ***
!    ***************************

!    A phi^4 potential indicates a self interaction of the scalar
!    field, its strenght is characterized by a dimensionless coupling
!    constant lambda and has the form:
!             
!    V  =  lambda phi^4 / 4

     else if (scalarpotential=="phi4") then

        scalar_V(l,:)  = 0.5d0*scalar_mass**2*scalar_phi(l,:)**2 &
                       + 0.25d0*scalar_lambda*scalar_phi(l,:)**4
        scalar_VP(l,:) = scalar_mass**2*scalar_phi(l,:) &
                       + scalar_lambda*scalar_phi(l,:)**3

        if (cosmic_run) then
           cosmobg_scalar_V(l)  = 0.5d0*scalar_mass**2*cosmobg_scalar_phi(l)**2 &
                                + 0.25d0*scalar_lambda*cosmobg_scalar_phi(l)**4
           cosmobg_scalar_VP(l) = scalar_mass**2*cosmobg_scalar_phi(l) &
                                + scalar_lambda*cosmobg_scalar_phi(l)**3
        end if


!    *******************
!    ***   UNKNOWN   ***
!    *******************

     else

        print *, 'Unknown potential type ...'
        print *, 'Aborting! (subroutine potential)'
        print *
        call die

     end if

  end if


! ********************************
! ***   COMPLEX SCALAR FIELD   ***
! ********************************

  if (contains(mattertype,"complex")) then

!    **************************
!    ***   ZERO POTENTIAL   ***
!    **************************

!    By default, set to zero.

     if (complexpotential=="none") then

        complex_V(l,:)   = 0.d0
        complex_VPR(l,:) = 0.d0
        complex_VPI(l,:) = 0.d0

        if (cosmic_run) then
           cosmobg_complex_V(l)   = 0.d0
           cosmobg_complex_VPR(l) = 0.d0
           cosmobg_complex_VPI(l) = 0.d0
        end if


!    ***************************
!    ***   PHI^2 POTENTIAL   ***
!    ***************************

!    A quadratic potential indicates a mass term
!    for the scalar field and has the form:
!             
!    V  =  m^2 |phi|^2 / 2

     else if (complexpotential=="phi2") then

        complex_V(l,:) = 0.5d0*complex_mass**2*(complex_phiR(l,:)**2 + complex_phiI(l,:)**2)

        complex_VPR(l,:) = complex_mass**2*complex_phiR(l,:)
        complex_VPI(l,:) = complex_mass**2*complex_phiI(l,:)

        if (cosmic_run) then
           cosmobg_complex_V(l) = 0.5d0*complex_mass**2*(cosmobg_complex_phiR(l)**2 + cosmobg_complex_phiI(l)**2)
           cosmobg_complex_VPR(l) = complex_mass**2*cosmobg_complex_phiR(l)
           cosmobg_complex_VPI(l) = complex_mass**2*cosmobg_complex_phiI(l)
        end if


!    ***************************
!    ***   PHI^4 POTENTIAL   ***
!    ***************************

!    A phi^4 potential indicates a self interaction of the scalar
!    field, its strenght is characterized by a dimensionless coupling
!    constant lambda and has the form:
!             
!    V  =  lambda |phi|^4 / 4

     else if (complexpotential=="phi4") then

        complex_V(l,:) = 0.5d0*complex_mass**2*(complex_phiR(l,:)**2 + complex_phiI(l,:)**2) &
                       + 0.25d0*complex_lambda*(complex_phiR(l,:)**2 + complex_phiI(l,:)**2)**2

        complex_VPR(l,:) = complex_mass**2*complex_phiR(l,:) &
                         + complex_lambda*complex_phiR(l,:)*(complex_phiR(l,:)**2 + complex_phiI(l,:)**2)
        complex_VPI(l,:) = complex_mass**2*complex_phiI(l,:) &
                         + complex_lambda*complex_phiI(l,:)*(complex_phiR(l,:)**2 + complex_phiI(l,:)**2)

        if (cosmic_run) then
           cosmobg_complex_V(l) = 0.5d0*complex_mass**2*(cosmobg_complex_phiR(l)**2 + cosmobg_complex_phiI(l)**2) &
              + 0.25d0*complex_lambda*(cosmobg_complex_phiR(l)**2 + cosmobg_complex_phiI(l)**2)**2
           cosmobg_complex_VPR(l) = complex_mass**2*cosmobg_complex_phiR(l) &
              + complex_lambda*cosmobg_complex_phiR(l)*(cosmobg_complex_phiR(l)**2 + cosmobg_complex_phiI(l)**2)
           cosmobg_complex_VPI(l) = complex_mass**2*cosmobg_complex_phiI(l) &
              + complex_lambda*cosmobg_complex_phiI(l)*(cosmobg_complex_phiR(l)**2 + cosmobg_complex_phiI(l)**2)
        end if


!    *******************
!    ***   UNKNOWN   ***
!    *******************

     else

        print *, 'Unknown potential type ...'
        print *, 'Aborting! (subroutine potential)'
        print *
        call die

     end if

  end if


! ************************
! ***   NONMIN FIELD   ***
! ************************

  if (contains(mattertype,"nonmin")) then

!    **************************
!    ***   ZERO POTENTIAL   ***
!    **************************

!    For the moment the nomin field is always
!    assumed to have zero potential.

     nonmin_V(l,:)  = 0.d0
     nonmin_VP(l,:) = 0.d0

  end if


! ******************************
! ***   GHOST SCALAR FIELD   ***
! ******************************

! The convention of the cpde is such that for the ghost
! scalar field we take the potential terms with OPPOSITE
! sign to the standard ones.
!
! The potential then enters with the usual sign in the
! stress-energy tensor, but since here it is defined as
! negative it has the opposite effect.
!
! The potential enters with the opposite sign to the usual
! one in the Klein-Gordon equation, but again since it is
! negative it results in the standard Klein-Gordon equation.

  if (contains(mattertype,"ghost")) then

!    **************************
!    ***   ZERO POTENTIAL   ***
!    **************************

!    By default, set to zero.

     if (ghostpotential=="none") then

        ghost_V(l,:)  = 0.d0
        ghost_VP(l,:) = 0.d0


!    ***************************
!    ***   PHI^2 POTENTIAL   ***
!    ***************************

!    A quadratic potential indicates a mass term
!    for the scalar field and has the form:
!             
!    V  = - m^2 |phi|^2 / 2

     else if (ghostpotential=="phi2") then

        ghost_V(l,:)  = - 0.5d0*ghost_mass**2*ghost_phi(l,:)**2 
        ghost_VP(l,:) = - ghost_mass**2*ghost_phi(l,:)

        if (cosmic_run) then

        end if


!    ***************************
!    ***   PHI^4 POTENTIAL   ***
!    ***************************

!    A phi^4 potential indicates a self interaction of the scalar
!    field, its strenght is characterized by a dimensionless coupling
!    constant lambda and has the form:
!             
!    V  = - lambda |phi|^4 / 4

     else if (ghostpotential=="phi4") then

        ghost_V(l,:) = - 0.5d0*ghost_mass**2*ghost_phi(l,:)**2 &
                     - 0.25d0*ghost_lambda*ghost_phi(l,:)**4

        ghost_VP(l,:) = - ghost_mass**2*ghost_phi(l,:) - ghost_lambda*ghost_phi(l,:)**3

        if (cosmic_run) then

        end if

!    *******************
!    ***   UNKNOWN   ***
!    *******************

     else

        print *, 'Unknown potential type ...'
        print *, 'Aborting! (subroutine potential)'
        print *
        call die

     end if

  end if


! **************************************
! ***   COMPLEX GHOST SCALAR FIELD   ***
! **************************************

! The convention of the cpde is such that for the ghost
! scalar field we take the potential terms with OPPOSITE
! sign to the standard ones.
!
! The potential then enters with the usual sign in the
! stress-energy tensor, but since here it is defined as
! negative it has the opposite effect.
!
! The potential enters with the opposite sign to the usual
! one in the Klein-Gordon equation, but again since it is
! negative it results in the standard Klein-Gordon equation.

  if (contains(mattertype,"complexghost")) then

!    **************************
!    ***   ZERO POTENTIAL   ***
!    **************************

!    By default, set to zero.

     if (complexghostpotential=="none") then

        complexghost_V(l,:)   = 0.d0
        complexghost_VPR(l,:) = 0.d0
        complexghost_VPI(l,:) = 0.d0

        if (cosmic_run) then

        end if


!    ***************************
!    ***   PHI^2 POTENTIAL   ***
!    ***************************

!    A quadratic potential indicates a mass term
!    for the scalar field and has the form:
!             
!    V  = - m^2 |phi|^2 / 2

     else if (complexghostpotential=="phi2") then

        complexghost_V(l,:) = - 0.5d0*complexghost_mass**2*(complexghost_phiR(l,:)**2 + complexghost_phiI(l,:)**2)

        complexghost_VPR(l,:) = - complexghost_mass**2*complexghost_phiR(l,:)
        complexghost_VPI(l,:) = - complexghost_mass**2*complexghost_phiI(l,:)

        if (cosmic_run) then

        end if


!    ***************************
!    ***   PHI^4 POTENTIAL   ***
!    ***************************

!    A phi^4 potential indicates a self interaction of the scalar
!    field, its strenght is characterized by a dimensionless coupling
!    constant lambda and has the form:
!             
!    V  = - lambda |phi|^4 / 4

     else if (complexghostpotential=="phi4") then

        complexghost_V(l,:) = - 0.5d0*complexghost_mass**2*(complexghost_phiR(l,:)**2 + complexghost_phiI(l,:)**2) &
           - 0.25d0*complexghost_lambda*(complexghost_phiR(l,:)**2 + complexghost_phiI(l,:)**2)**2

        complexghost_VPR(l,:) = - complexghost_mass**2*complexghost_phiR(l,:) &
           - complexghost_lambda*complexghost_phiR(l,:)*(complexghost_phiR(l,:)**2 + complexghost_phiI(l,:)**2)
        complexghost_VPI(l,:) = - complexghost_mass**2*complexghost_phiI(l,:) &
           - complexghost_lambda*complexghost_phiI(l,:)*(complexghost_phiR(l,:)**2 + complexghost_phiI(l,:)**2)

        if (cosmic_run) then

        end if


!    *******************
!    ***   UNKNOWN   ***
!    *******************

     else

        print *, 'Unknown potential type ...'
        print *, 'Aborting! (subroutine potential)'
        print *
        call die

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine potential


