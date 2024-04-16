!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/checkpoint.f90,v 1.21 2024/04/16 20:33:32 malcubi Exp $

  subroutine checkpointsave

! **********************
! ***   CHECKPOINT   ***
! **********************

! This routine saves checkpoint files for a restart.

! Include modules.

  use param
  use arrays
  use procinfo
  use mpi

! Extra variables.

  implicit none

  logical firstcall

  integer i,j,l
  integer nvars
  integer, allocatable, dimension (:) :: commas

  character(1)   aa
  character(2)   filen
  character(20)  filestatus
  character(100) outdir,outtime

  character(50), allocatable, dimension (:) :: outvars

  data firstcall / .true. /

  save firstcall,nvars,outvars


! ***************************************
! ***   CREATE CHECKPOINT DIRECTORY   ***
! ***************************************

! Create a checkpoint directory and copy the parameter
! file inside it.

  write(outtime,"(F12.5)") t(0)
  outdir = trim(directory)//'/checkpoint'//'_t='//adjustl(trim(outtime))

  if (rank==0) then
     call system('mkdir -p '//outdir)
     call system('cp '//trim(directory)//'/'//trim(parfile)//' '//outdir)
  end if

! wait here until all processors reach this point.  We don't want
! processors saving anything until the directory has been created.

  call MPI_Barrier(MPI_COMM_WORLD,ierr)


! ***********************************
! ***   IS THIS THE FIRST CALL?   ***
! ***********************************

! On first call figure out what needs output.

  if (firstcall) then

     firstcall = .false.

!    Find out length of string "checkvars".

     l = len_trim(checkvars)

!    Find out how many variables need output.

     nvars = 1

     do i=1,l
        aa = checkvars(i:i)
        if (aa==",") then
           nvars = nvars + 1
        end if
     end do

!    Identify comma positions.

     allocate(commas(0:nvars))

     j = 0
     commas(0) = 0
     commas(nvars) = l+1

     do i=1,l
        aa = checkvars(i:i)
        if (aa==",") then
           j = j+1
           commas(j)=i
        end if
     end do

!    Now read variable names, and eliminate spaces.

     allocate(outvars(1:nvars))

     do i=1,nvars
        outvars(i) = checkvars(commas(i-1)+1:commas(i)-1)
        outvars(i) = trim(adjustl(outvars(i)))
        !print *, 'Array: ',outvars(i)
     end do

  end if


! ***********************************
! ***   SAVE CRUCIAL PARAMETERS   ***
! ***********************************

! Some parameters are crucial for a restart:
! grid structure, matter content,
! etcetera, and they need to be the same.
!
! It is not enough to re-read the old parameter
! file since many may have had their default
! values and that won't show up in the parfile.
! So I save them here to the file "checkparam".

  if (rank==0) then

!    Open file.

     open(1,file=trim(outdir)//'/checkparam',form='formatted',status='replace')

!    Save number of processors used in run.

     write(1,*) 'nproctot = ',nproctot

!    Save grid spacing.

     write(1,*) 'dr0 = ',dr0

!    Save grid size.

     write(1,*) 'Nrtotal = ',Nrtotal

!    Save number of levels.

     write(1,*) 'Nl = ',Nl

!    Save matter content.

     write(1,*) 'mattertype = ',trim(mattertype)

!    Is this a cosmological run?

     write(1,*) 'cosmic_run = ',cosmic_run

!    Close file.

     close(1)

  end if


! *********************
! ***   SAVE DATA   ***
! *********************

! For a checkpoint we always replace the old files.

  filestatus = 'replace'

! Save name of parfile and list of variables
! that have output.

  if (rank==0) then

     open(1,file=trim(outdir)//'/checklist',form='formatted',status='replace')

     write(1,*) trim(parfile)

     do i=1,nvars
        write(1,*) trim(outvars(i))
     end do

     close(1)

  end if

! Loop over grid levels.

  do l=0,Nl-1

!    Append grid level to file name.

     if (l<10) then
        write(filen,'(i1)') l
     else
        write(filen,'(i2)') l
     end if

!    Save all variables that require output. Notice that
!    we call the subroutine "save1Dvariable" that is
!    found in the file "save1Ddata.f90".

     do i=1,nvars

!       Grab array.

        call grabarray(trim(outvars(i)))

!       Save 0D arrays.

        if (associated(savevar0D).and.(.not.associated(savevar))) then

           call save0Dvariable(outdir,trim(outvars(i)),i,s(0),t(0),filestatus,-1)

!       Save standard 1D arrays.

        else

           call save1Dvariable(outdir,trim(outvars(i)),i,s(0),filestatus,-1)

!          Save internal boundaries, but not for coordinates (only processor at boundary).

           if (outvars(i)=='r') cycle

           if (rank==size-1) then

              open(1,file=trim(outdir)//'/'//trim(outvars(i))//'_bound'//trim(filen),form='formatted',status='replace')

              write(1,'(3ES24.16)') t(l),t1(l),t2(l)
              write(1,*)

              do j=0,ghost-1
                 write(1,'(3ES24.16)') grabvar_bound(l,j,1),grabvar_bound(l,j,2),grabvar_bound(l,j,3)
              end do

              close(1)

           end if

        end if

     end do

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine checkpointsave









  subroutine checkpointrestart

! *****************************
! ***   CHECKPOINT RESTART  ***
! *****************************

! This routine restarts from checkpoint files.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical flag
  logical cosmic_run_new

  integer i,j,l,nvars
  integer nproctot_new,Nr_new,Nl_new

  real(8) dr_new,tdata

  character(2)    filen
  character(1000) matter_new,parfile_old,line
  character(1000) checkparam

  character(50) outvars(1:1000)


! *****************
! ***   START   ***
! *****************

! Message to screen.

  if (rank==0) then
     print *, 'Starting from checkpoint!'
     print *
  end if

! Check if the checkpoint directory exists.
! It turns out that although GFORTRAN can
! use INQUIRE to check if a directory exists,
! the INTEL compiler can't.  So now I check
! to see if the file "." exists inside the
! given directory (it should always be there).

  inquire(FILE='./'//trim(checkpointfile)//'/.',EXIST=flag)

  if (.not.flag) then
     if (rank==0) then
        print *, 'Directory "',trim(checkpointfile),'" does not exist.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if


! ************************
! ***   SANITY CHECK   ***
! ************************

! We need to make sure that the values of certain
! crucial parameters are the same as in the old run,
! or the new run will not work.
!
! Notice that many parameters can in fact change
! and it won't matter, but parameters related to
! grid structure, matter content, etcetera, need
! to be the same.
!
! The first step it to save the values that those parameters
! have after reading the new parameter file.

  nproctot_new = nproctot

  dr_new = dr0

  Nr_new = Nrtotal

  Nl_new = Nl

  matter_new = mattertype

  cosmic_run_new = cosmic_run

! Now parse the file "checkparam".

  checkparam = trim(checkpointfile)//'/checkparam'

  call parse(rank,trim(checkparam))

! Check number of processors.

  if (nproctot_new/=nproctot) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint we must use the same number of processors.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check grid spacing.

  if (dr0/=dr_new) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the grid spacing must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check grid size.

  if (Nrtotal/=Nr_new) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the grid size must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check number of levels.

  if (Nl/=Nl_new) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the number of levels must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check matter content.

  if (matter_new/=mattertype) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the matter type must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if

! Check if this is a cosmological run.

  if (cosmic_run_new.neqv.cosmic_run) then
     if (rank==0) then
        print *, 'When restarting from a checkpoint the value of cosmic_run must be the same.'
        print *, 'Aborting (subroutine checkpoint.f90)'
        print *
     end if
     call die
  end if


! *******************************
! ***   READ VARIABLE NAMES   ***
! *******************************

  open(1,file=trim(checkpointfile)//'/checklist',form='formatted',status='old')

  read(1,*) parfile_old

  i = 1

  do while(.true.)
     read(1,*,end=100) outvars(i)
     i=i+1
  end do

  100 nvars = i-1

  close(1)


! *********************
! ***   READ DATA   ***
! *********************

! Loop over grid variables.

  do i=1,nvars

!    Grab array.

     call grabarray(trim(outvars(i)))

!    Read standard 1D arrays.

     if (associated(savevar)) then

        call read1Ddata(outvars(i),checkpointfile)

!    Read 0D arrays.

     else

!       Loop over grid levels.

        do l=0,Nl-1

!          Append grid level to file name.

           if (l<10) then
             write(filen,'(i1)') l
           else
              write(filen,'(i2)') l
           end if

!          Open file.

           open(1,file=trim(checkpointfile)//'/'//trim(outvars(i))//trim(filen)//'.tl',form='formatted',status='old')

!          Read data.  Remember that the first line is a comment.

           read(1,'(A)') line
           read(1,*) tdata,savevar0D(l)

!          Close file.

           close(1)

        end do

     end if

  end do


! ******************************
! ***   READ BOUNDARY DATA   ***
! ******************************

! Loop over grid levels.

  do l=0,Nl-1

!    Append grid level to file name.

     if (l<10) then
        write(filen,'(i1)') l
     else
        write(filen,'(i2)') l
     end if

!    Loop over grid variables.

     do i=1,nvars

        call grabarray(trim(outvars(i)))

        continue

!       Don't read boundary data for coordinates (it doesn't exist).

        if (outvars(i)=='r') cycle

!       Read boundary data.

        if (associated(savevar)) then

!          Open file.

           open(1,file=trim(checkpointfile)//'/'//trim(outvars(i))//'_bound'//trim(filen),form='formatted',status='old')

!          Read time data (remember there is a blank line after it).

           read(1,*) t(l),t1(l),t2(l)
           read(1,*)

!          Read boundary data.

           do j=0,ghost-1
              read(1,*) grabvar_bound(l,j,1),grabvar_bound(l,j,2),grabvar_bound(l,j,3)
           end do

!          Close file.

           close(1)

        end if

     end do

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine checkpointrestart

