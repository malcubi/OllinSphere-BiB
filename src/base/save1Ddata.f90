!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/save1Ddata.f90,v 1.30 2021/06/21 16:37:06 malcubi Exp $

  subroutine save1Ddata

! ********************************
! ***   SAVE 1D DATA TO FILE   ***
! ********************************

! This routine saves 1D data, that is the
! whole spatial arrays at a given time.

! Include modules.

  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical firstcall

  integer i,j,l
  integer maxfiles
  integer, allocatable, dimension (:) :: commas

  character(1)  aa
  character(20) filestatus

  data firstcall / .true. /

  save firstcall


! ***********************************
! ***   IS THIS THE FIRST CALL?   ***
! ***********************************

! On first call, replace file and figure
! out what needs output.

  if (firstcall) then

     firstcall = .false.

!    File status.

     filestatus = 'replace'

!    Find out length of string "outvars1D".

     l = len_trim(outvars1D)

!    Find out how many variables need output.

     nvars1D = 1

     do i=1,l
        aa = outvars1D(i:i)
        if (aa==",") then
           nvars1D = nvars1D + 1
        end if
     end do

!    Identify comma positions.

     allocate(commas(0:nvars1D))

     j = 0
     commas(0) = 0
     commas(nvars1D) = l+1

     do i=1,l
        aa = outvars1D(i:i)
        if (aa==",") then
           j = j+1
           commas(j)=i
        end if
     end do

!    Now read variable names, and eliminate spaces.

     allocate(outvars1Darray(1:nvars1D))

     do i=1,nvars1D
        outvars1Darray(i) = outvars1D(commas(i-1)+1:commas(i)-1)
        outvars1Darray(i) = trim(adjustl(outvars1Darray(i)))
     end do

!    Check if any name is repeated.

     if (rank==0) then
        do i=1,nvars1D
           do j=1,nvars1D
              if ((i.ne.j).and.(trim(outvars1Darray(i))==trim(outvars1Darray(j)))) then
                 print *
                 print *, 'Error in parfile, array name repeated in outvars1D: ',outvars1Darray(i)
                 print *, 'Aborting! (subroutine save1Ddata)'
                 print *
                 call die
              end if
           end do
        end do
     end if

!    If we are not closing the files, check that
!    we do not have too many open at once.

     if ((.not.closefiles).and.(rank==0)) then

        call system('ulimit -n > .maxfiles')
        open(1,file='.maxfiles')
        read(1,*) maxfiles
        close(1)
        call system('/bin/rm .maxfiles')

        if ((5*Nl*nvars0D+(Nl+1)*nvars1D)>maxfiles) then
           write(*,'(A,I0,A)') ' You can not have more than ',maxfiles,' files open at once.'
           print *, 'Consider increasing this by adding e.g. "ulimit -n 1000" to you .bash_profile'
           print *, 'Aborting! (subroutine save1Ddata)'
           print *
           call die
        end if

     end if

! Not first call.

  else

     filestatus = 'old'

  end if


! *********************
! ***   SAVE DATA   ***
! *********************

! Save all variables that require output.

  do i=1,nvars1D
     call grabarray(trim(outvars1Darray(i)))
     call save1Dvariable(directory,trim(outvars1Darray(i)),i,s(0),filestatus,-1)
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine save1Ddata









  subroutine save1Dvariable(outdir,varname,nvar,step,filestatus,ll)

! This subroutine outputs 1D files.
!
! The parameters passed to this routine are:
!
! outdir      =  Name of output directory.
! varname     =  A string with the name of the variable to be saved.
! nvar        =  Number of the variable to be saved in teh list outvars1Darray.
!                This is needed because we keep a file open for each variable.
! step        =  Time step. This is needed so that we known when to close the
!                open files (when step=Nt).
! filestatus  =  "replace" or "old".
! ll          =  If ll is less than 0 we output all time levels, otherwise we
!                only output the time lebel ll.

! Load modules.

  use mpi
  use param
  use arrays
  use procinfo

! Declare variables.

  implicit none

  logical savemerged

  integer i,l,p
  integer ll,lmin,lmax
  integer imin,imax,Naux
  integer nvar,step
  integer unit1,unit2
  integer status(MPI_STATUS_SIZE)

  real(8) tiny,aux
  real(8), dimension (1-ghost:Nrmax) :: var
 
  character(len=5)  filen
  character(len=10) form
  character(len=*)  varname,outdir,filestatus


! ************************
! ***   SANITY CHECK   ***
! ************************

  if (associated(savevar0D).and.(.not.associated(savevar))) then
     print *, 'ZEROD declared variables can not have 1D output: ',varname
     print *, 'Aborting! (subroutine save1Ddata)'
     print *
     call die
  end if


! *************************
! ***   OUTPUT FORMAT   ***
! *************************

! For checkpoint we need output with many significant figures,
! while for plotting we don't.  The easiest way to ensure this
! is saving 15 significant figures whenever the file is replaced,
! and saving fewer when it is only appended to.

  if (filestatus=='replace') then
     form = "(2ES24.16)"
  else
     form = "(2ES16.8)"
  end if


! *************************************
! ***   LOOP OVER ALL GRID LEVELS   ***
! *************************************

  tiny = 1.d-40

! Do we save merged data?

  savemerged = .false.

! Find total number of points.

  Naux = Nrmax + ghost

! Open merged file.

  if (savemerged.and.(rank==0).and.(ll<0)) then

!    Define output unit 2.

     unit2 = 10000 + (nvar-1)

!    Open file for merged data.

     if (filestatus=='replace') then
        open(unit2,file=trim(outdir)//'/'//varname//'-m.rl',form='formatted', &
        status='replace')
     else if (closefiles) then
        open(unit2,file=trim(outdir)//'/'//varname//'-m.rl',form='formatted', &
        status='old',position='append')
     end if

     if (commenttype=='xgraph') then
        write(unit2,"(A8,ES18.10)") '"Time = ',t(0)
     else
        write(unit2,"(A8,ES18.10)") '#Time = ',t(0)
     end if

  end if

! Loop over all grid levels.

  if (ll<0) then
     lmin = 0
     lmax = Nl-1
  else
     lmin = ll
     lmax = ll
  end if

  do l=lmin,lmax

!    Processor 0 does output.

     if (rank==0) then

        if (l<10) then
           write(filen,'(i1)') l
        else
           write(filen,'(i2)') l
        end if

!       Define output unit 1 (for different grid levels).

        unit1 = 11000 + l + Nl*(nvar-1)

!       Open file for level l.

        if (filestatus=='replace') then
           open(unit1,file=trim(outdir)//'/'//varname//trim(filen)//'.rl',form='formatted', &
           status='replace')
        else if (closefiles) then
           open(unit1,file=trim(outdir)//'/'//varname//trim(filen)//'.rl',form='formatted', &
           status='old',position='append')
        end if

!       Write current time.

        if (commenttype=='xgraph') then
           write(unit1,"(A8,ES18.10)") '"Time = ',t(l)
        else
           write(unit1,"(A8,ES18.10)") '#Time = ',t(l)
        end if

!       Save data from processor 0.

        if (size==1) then
           imax = Nrl(0)
        else
           imax = Nrl(0) - ghost
        end if

        do i=1-ghost,imax
           write(unit1,form) r(l,i),savevar(l,i)+tiny
        end do

!       Save data from processor 0 to merged file.

        if ((savemerged).and.(ll<0)) then

           if (l==Nl-1) then
              imin = 1-ghost
           else
              imin = int((rright(size-1,l+1)-rleft(0,l))/dr(l))
           end if

           do i=imin,imax
              write(unit2,"(2ES18.10)") r(l,i),savevar(l,i)+tiny
           end do

        end if

!       Iterate over other processors.

        do p=1,size-1

!          Receive and save data from other processors.

           call MPI_RECV(var,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

           if (p==size-1) then
              imax = Nrl(p)
           else
              imax = Nrl(p) - ghost
           end if

           aux = dr(l)*(Nmin(p)-Nmin(0))

           do i=1,imax
              write(unit1,form) r(l,i)+aux,var(i)+tiny ! + 0.01*p
           end do

!          Save data from processor p to merged file.

           if ((savemerged).and.(ll<0)) then

              if ((l==Nl-1).or.(rleft(p,l)>rright(size-1,l+1))) then
                 imin = 1
              else
                 imin = int((rright(size-1,l+1)-rleft(p,l))/dr(l))
              end if

              do i=imin,imax
                 write(unit2,form) r(l,i)+aux,var(i)+tiny ! + 0.01*p
              end do

           end if

        end do

!       Leave two blank spaces before next time.
!       The reason to leave two spaces is that 'gnuplot' asks
!       for two spaces to distinguish different records.

        write(unit1,*)
        write(unit1,*)

!       Close files.

        if ((closefiles).or.(step==Nt)) then
           close(unit1)
        end if

!    Other processors send data to processor 0.

     else

        call MPI_SEND(savevar(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

     end if

  end do

! Close file for merged data.

  if (savemerged.and.(rank==0).and.(ll<0)) then

     write(unit2,*)
     write(unit2,*)

     if ((closefiles).or.(step==Nt)) then
        close(unit2)
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine save1Dvariable
