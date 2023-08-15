!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/main.f90,v 1.21 2022/07/12 23:58:01 malcubi Exp $

! ********************************
! ***   PROGRAM  OLLINSPHERE   ***
! ********************************

  program main

! Load MPI module.  If we compile without MPI this will
! just load a dummy module.

  use mpi

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  integer p,Naux

  character(100) message,string1,string2

! The variables introduced above are:
!
! parfile       Name of parameter file.
!
! message       Function to format messages to screen (see functions.f90)
! string*       Auxiliary strings.


! **************************
! ***   INITIALIZE MPI   ***
! **************************

! When compiling without MPI this call does nothing at all.

  call MPI_INIT(ierr)


! ********************************************************
! ***   FIND OUT NUMBER OF PROCESSORS AND LOCAL RANK   ***
! ********************************************************

! When compiling without MPI these calls just return size=1 and rank=0.

  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  nproctot = size


! ***************************
! ***   INITIAL MESSAGE   ***
! ***************************

  if (rank==0) then

     print *
     print *, '******************************************'
     print *, '******************************************'
     print *, '***                                    ***'
     print *, '***            OLLINSPHERE             ***'
     print *, '***                                    ***'
     print *, '***   Evolving Einstein''s equations    ***'
     print *, '***       in spherical symmetry        ***'
     print *, '***     using the BSSN formulation     ***'
     print *, '***                                    ***'    
     print *, '******************************************'
     print *, '******************************************'
     print *

  end if


! *********************************
! ***   CALL PARAMETER PARSER   ***
! *********************************

! Get name of parameter file.

  call getarg(1,parfile)

  if (parfile==" ") then
     if (rank==0) then
        print *, 'Missing parfile name.'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if

! Call parser.

  call parse(rank,parfile)


! *************************
! ***   SANITY CHECKS   ***
! *************************

! Various sanity checks.  Might add more later.

! Ninfo must be positive.

  if (Ninfo<=0) then
     if (rank==0) then
        print *, 'Ninfo must be positive'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if

! Noutput must be non-negative.  Notice that 0 is allowed,
! but this is a special case for testing purposes only.
! When the corresponding parameter Noutput is equal to 0,
! we do output al ALL time steps for each box and level,
! and not only at the coarse time steps.

  if (Noutput0D<0) then
     if (rank==0) then
        print *, 'Noutput0D must be positive'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if

  if (Noutput1D<0) then
     if (rank==0) then
        print *, 'Noutput1D must be positive'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if

! Ncheckpoint must be positive.

  if (Ncheckpoint<=0) then
     if (rank==0) then
        print *, 'Ncheckpoint must be positive'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  end if


! **************************
! ***   OTHER MESSAGES   ***
! **************************

  if (rank==0) then

     print *, '******************************************'
     print *, '******************************************'
     print *, '***                                    ***'

!    Slicing message.

     string1 = "slicing:"
     string2 = trim(adjustl(slicing))

     print *, message(string1,string2)

!    Shift message.

     string1 = "shift:"
     string2 = trim(adjustl(shift))

     print *, message(string1,string2)

!    Type of matter.

     print *, '***                                    ***'

     string1 = "matter type:"
     string2 = trim(adjustl(mattertype))

     print *, message(string1,string2)

!    Initial data message.

     print *, '***                                    ***'

     string1 = "initial data:"
     string2 = trim(adjustl(idata))

     print *, message(string1,string2)

!    Evolution method message.

     print *, '***                                    ***'

     string1 = "evolution method:"
     string2 = trim(adjustl(integrator))

     print *, message(string1,string2)

     string1 = "order of integration:"
     string2 = trim(adjustl(order))

     print *, message(string1,string2)

!    Boundary type message.

     string1 = "boundary condition:"
     string2 = trim(adjustl(boundtype))

     print *, message(string1,string2)

!    Regularization message.

     print *, '***                                    ***'

     string1 = "regularization of origin:"

     if (nolambda) then
        string2 = "off"
     else
        string2 = "on"
     end if

     print *, message(string1,string2)

!    Output directory message.

     print *, '***                                    ***'

     string1 = "output  directory:"
     string2 = " "
     print *, message(string1,string2)

     string1 = trim(adjustl(directory))
     string2 = " "

     print *, message(string1,string2)

!    End message.

     print *, '***                                    ***'   
     print *, '******************************************'
     print *, '******************************************'
     print *

  end if


! ***********************************
! ***   CREATE OUTPUT DIRECTORY   ***
! ***********************************

! Create output directory and copy the parameter
! file inside it.  But in case it already existed,
! remove it first.

  if (rank==0) then
     call system('rm -rf '//trim(directory))
     call system('mkdir -p '//trim(directory))
     call system('cp '//trim(parfile)//' '//trim(directory))
  end if


! ***********************
! ***   GHOST ZONES   ***
! ***********************

! Sanity check.

  if (ghost/=0) then

     if (rank==0) then
        print *
        print *, 'The value of "ghost" should not be fixed in the parameter file.'
        print *, 'Aborting! (subroutine main)'
        print *
     end if

     call die

  else

!    For second order I use 2 ghost zones and for
!    fourth order 3 ghost zones.  This might seem
!    one more ghost zone than needed, but one must
!    remember that for dissipation and one-sided
!    derivatives we need an extra ghost zone.

     if (order=="two") then
        ghost = 2
     else if (order=="four") then
        ghost = 3
     else if (order=="six") then
        ghost = 4
     else if (order=="eight") then
        ghost = 5
     end if

  end if

! Output number of ghost zones.

  if (rank==0) then
     write (*,'(A,I0)') ' Ghost zones = ',ghost
     print *
  end if


! *************************************
! ***   MESSAGE ABOUT GRID LEVELS   ***
! *************************************

! Output number of grid levels.

  if (Nl<1) then
     if (rank==0) then
        print *, 'The number of grid levels l must be at least 1.'
        print *, 'Aborting! (subroutine main)'
        print *
     end if
     call die
  else
     if (rank==0) then
        if (Nl==1) then
           write (*,'(A)')    ' No refinement levels, only base grid'
        else
           write (*,'(A,I0)') ' Total grid levels = ',Nl
        end if
        print *
     end if
  end if


! ***************************************
! ***   OUTPUT NUMBER OF PROCESSORS   ***
! ***************************************

! Output total number of processors.

  if (rank==0) then
     if (size==1) then
        print *, 'Running on 1 processor'
        print *
     else
        write(*,'(A,I0,A)') ' Running on ',size,' processors'
        print *
     end if
  end if


! ********************************
! ***   DOMAIN DECOMPOSITION   ***
! ********************************

! Allocate arrays (Nmim,Nmax,Nrl,rleft,rright).

  allocate(Nmin(0:size-1),Nmax(0:size-1),Nrl(0:size-1))
  allocate(rleft(0:size-1,0:Nl-1),rright(0:size-1,0:Nl-1))

! The idea here is that each processor owns grid points
! with:
!
! i = 1,...,Nr-ghost      (interior points)
!
! Boundary points, or points belonging to a different
! proccesor that need to be synchronized, are points
! such that:
!
! i = 1-ghost,...,0    (lower boundary)
! i = Nr-ghost,...,Nr  (upper boundary)
!
! For example, for ghost=1, the boundary points are i=0
! and i=Nr.  For ghost=2, the boundary points are i=-1,0
! and i=Nr-1,Nr.

! Sanity check.

  if (Nrtotal<2+ghost) then

     if (rank==0) then
        print *
        print *, 'The value of Nrtotal should be greater than 2 plus the number of ghost zones.'
        print *, 'Aborting! (subroutine main)'
        print *
     end if

     call die

  else if (Nr/=1) then

     if (rank==0) then
        print *, 'WARNING: The value of Nr should not be fixed'
        print *, 'in the parameter file, use Nrtotal instead.'
        print *
        print *, 'Setting  Nrtotal = Nr'
        print *
     end if

     Nrtotal = Nr

  end if

! Integer division to figure out local number of grid points.

  Naux = Nrtotal/size

  do p=0,size-1

     Nmax(p) = (p+1)*Naux
     Nmin(p) = p*Naux

!    Take ghost zones into account, and make sure
!    that the last processor always has Nmax=Ntotal.

     if (p==size-1) then
        Nmax(p) = Nrtotal
     else
        Nmax(p) = Nmax(p) + ghost
     end if

!    Now find local number of grid points in each direction.

     Nrl(p) = Nmax(p) - Nmin(p)

  end do

! print *, rank, Nmin(rank),Nmax(rank),Nrl(rank)

! Copy local number of grid points to Nr.

  Nr = Nrl(rank)


! ***************************
! ***   ALLOCATE ARRAYS   ***
! ***************************

! Find out largest value of Nr.  This is important
! since all arrays will have to be of the same size
! in order to make communications easier.

  Nrmax = 0

  do p=0,size-1
     if (Nrl(p)>Nrmax) Nrmax=Nrl(p)
  end do

! Allocate arrays and initialize to zero.

  call allocatearrays('on')


! *********************
! ***   EVOLUTION   ***
! *********************

  call evolve


! *****************************
! ***   DEALLOCATE ARRAYS   ***
! *****************************

! Deallocate arrays.

  call allocatearrays('off')

! Deallocate arrays with processor information.

  deallocate(Nmin,Nmax,Nrl)
  deallocate(rleft,rright)


! ************************
! ***   FINALIZE MPI   ***
! ************************

! When compiling without MPI this call does nothing at all.

  call MPI_FINALIZE(ierr)


! ***************
! ***   END   ***
! ***************

  if (rank==0) then
     print *,'------------------------------'
     print *
     print *, 'PROGRAM OLLINSPHERE HAS FINISHED'
     print *
     print *, 'Have a nice day!'
     print *
     print *
     print *
  end if

  end program main

