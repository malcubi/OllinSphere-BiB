!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/save0Ddata.f90,v 1.26 2023/03/09 01:58:02 malcubi Exp $

  subroutine save0Ddata

! ********************************
! ***   SAVE 0D DATA TO FILE   ***
! ********************************

! This routine saves "0D" data to files. By 0D data I mean reduced
! quantities as functions of time (minimum, maximum, norms, etc).

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

!    Find out length of string "outvars0D".

     l = len_trim(outvars0D)

!    Find out how many variables need output.

     nvars0D = 1

     do i=1,l
        aa = outvars0D(i:i)
        if (aa==",") then
           nvars0D = nvars0D + 1
        end if
     end do

!    Identify comma positions.

     allocate(commas(0:nvars0D))

     j = 0
     commas(0) = 0
     commas(nvars0D) = l+1

     do i=1,l
        aa = outvars0D(i:i)
        if (aa==",") then
           j = j+1
           commas(j)=i
        end if
     end do

!    Now read variable names, and eliminate spaces.

     allocate(outvars0Darray(1:nvars0D))

     do i=1,nvars0D
        outvars0Darray(i) = outvars0D(commas(i-1)+1:commas(i)-1)
        outvars0Darray(i) = trim(adjustl(outvars0Darray(i)))
     end do

!    Check if any name is repeated.

     if (rank==0) then
        do i=1,nvars0D
           do j=1,nvars0D
              if ((i.ne.j).and.(trim(outvars0Darray(i))==trim(outvars0Darray(j)))) then
                 print *, 'Error in parfile, array name repeated in outvars0D: ',outvars0Darray(i)
                 print *, 'Aborting! (subroutine save0Ddata)'
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

        if (5*Nl*Nvars0D>maxfiles) then
           write(*,'(A,I0,A)') ' You can not have more than ',maxfiles,' files open at once.'
           print *, 'Consider increasing this by adding e.g. "ulimit -n 1000" to you .bash_profile'
           print *, 'Aborting! (subroutine save0Ddata)'
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

  do i=1,nvars0D
     call grabarray(trim(outvars0Darray(i)))
     call save0Dvariable(directory,trim(outvars0Darray(i)),i,s(0),t(0),filestatus,-1)
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine save0Ddata








  subroutine save0Dvariable(outdir,varname,nvar,step,time,filestatus,ll)

! This subroutine calculates minimum, maximum, norm 1, norm 2
! and total variation of a variable and outputs them.
!
!
! The parameters passed to this routine are:
!
! varname     =  A string with the name of the variable to be saved.
! outdir      =  Name of output directory.
! nvar        =  Number of the variable to be saved in the list outvars1Darray.
!                This is needed because we keep a file open for each variable.
! step        =  Time step. This is needed so that we known when to close the
!                open files (when step=Nt).
! filestatus  =  "replace" or "old".
! ll          =  If ll is negative we output all time levels, otherwise we
!                only output the time level ll.
!
! The quantities saved are:
!
! vmax:    Maximum value of the variable over the grid.
! vmin:    Minimum value of the variable over the grid.
!
! nm1:     Maximum absolute value.
! nm2:     Root mean square (rms).
!
! tvar:    Total variation.
!
! The routine also outputs the value of the function at the origin,
! and its value at the outer boundary.

! Load modules.

  use mpi
  use param
  use arrays
  use procinfo

! Declare variables.

  implicit none

  integer i,l,p
  integer ll,lmin,lmax
  integer imin,imax
  integer nvar,step
  integer unit
  integer bestl

  integer status(MPI_STATUS_SIZE)

  real(8) foundnan
  real(8) time,interp
  real(8) tiny,aux

  real(8), dimension(0:Nl-1) :: vmax,gmax,vmin,gmin
  real(8), dimension(0:Nl-1) :: nm1,gnm1,nm2,gnm2
  real(8), dimension(0:Nl-1) :: tvar,gtvar
  real(8), dimension(0:Nl-1) :: NN,NNtot

  character(len=1) comment
  character(len=5) filen
  character(len=*) varname,outdir,filestatus


! ************************
! ***   COMMENT TYPE   ***
! ************************

! Decide on comment type.

  if (commenttype=='xgraph') then
     comment = '"'
  else
     comment = '#'
  end if


! ******************************
! ***   FIND LMIN AND LMAX   ***
! ******************************

! Output for all time levels.

  if (ll<0) then

     lmin = 0
     lmax = Nl-1

! Output for a single time level.

  else

     lmin = ll
     lmax = ll

  end if


! *********************
! ***   0D ARRAYS   ***
! *********************

! For arrays declared as ZEROD we just output the value.

  if (associated(savevar0D).and.(.not.associated(savevar))) then

     if (rank==0) then
        do l=lmax,lmin,-1

           if (l<10) then
              write(filen,'(i1)') l
           else
              write(filen,'(i2)') l
           end if

           if (filestatus == 'replace') then
              open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.tl',form='formatted', &
              status='replace')
              write(1,*) comment//varname//'_max.tl'
           else
              open(1,file=trim(outdir)//'/'//varname//trim(filen)//'.tl',form='formatted', &
              status='old',position='append')
           end if

           write(1,"(2ES25.15e3)") t(l),savevar0D(l)

           close(1)

        end do
     end if

     return

  end if


! ******************************
! ***   STANDARD 1D ARRAYS   ***
! ******************************

  tiny = 1.d-50

! Loop over all grid levels. I do it from fine
! to coarse grid. The reason for this is that
! I want to catch any possible NaN's on the fine
! grids first.

  do l=lmax,lmin,-1

     if (l<10) then
        write(filen,'(i1)') l
     else
        write(filen,'(i2)') l
     end if


!    ****************************
!    ***   FIND LOCAL NORMS   ***
!    ****************************

     imin = 1

     if (rank==size-1) then
        imax = Nrl(rank)
     else
        imax = Nrl(rank)-ghost
     end if

!    Set NaN catching flag to false.

     foundnan = 0

!    Find minimum, maximum and norm 1.

     vmax(l) = -1.d10
     vmin(l) = +1.d10
     nm1(l)  =  0.d0

     do i=imin,imax

        aux = savevar(l,i)

        if (isnan(aux)) then
           print *
           write(*,"(A,A,A,3I6,E18.8)") ' Found NaN in value of ',varname,' at (proc,level,step,time):',rank,l,step,time
           write(*,"(A,I6,E18.8)") ' at point (i,r):',i,r(l,i)
           print *, 'Aborting! (subroutine save0Ddata.f90)'
           print *
           foundnan = 1.d0
           goto 100
        end if

        if (aux>vmax(l)) vmax(l) = aux
        if (aux<vmin(l)) vmin(l) = aux
        if (abs(aux)>nm1(l)) nm1(l) = abs(aux)

     end do

     100 continue

!    If we found a NaN die.

     aux = 0.0d0
     call MPI_Allreduce(foundnan,aux,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

     if (aux>0.d0) then
        call die
     end if

!    Find norm 2 (L2 norm).

     NN(l)  = 0.d0
     nm2(l) = 0.d0

     do i=imin,imax
        NN(l)  = NN(l) + 1.d0
        nm2(l) = nm2(l) + savevar(l,i)**2
     end do

!    Find total variation.

     tvar(l) = 0.d0

     do i=imin,imax
        tvar(l) = tvar(l) + abs(savevar(l,i)-savevar(l,i-1))
     end do


!    ***************************************************
!    ***   REDUCE QUANTITIES ACROSS ALL PROCESSORS   ***
!    ***************************************************

!    If we are running on one processor just copy.

     if (size==1) then

        gmax(l) = vmax(l)
        gmin(l) = vmin(l)

        gnm1(l) = nm1(l)

        NNtot(l) = NN(l)
        gnm2(l)  = sqrt(nm2(l)/NN(l))

        gtvar(l) = tvar(l)

!       If we are running on several processors
!       reduce quantities across processors.

     else

!       Global maximum.

        call MPI_Allreduce(vmax(l),gmax(l),1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

!       Global minimum.

        call MPI_Allreduce(vmin(l),gmin(l),1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)

!       Global nm1.

        call MPI_Allreduce(nm1(l),gnm1(l),1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

!       Global nm2.

        call MPI_Allreduce(NN(l),NNtot(l),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(nm2(l),gnm2(l),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

        gnm2(l) = sqrt(gnm2(l)/NNtot(l))

!       Global variation.

        call MPI_Allreduce(tvar(l),gtvar(l),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

     end if


!    *****************************
!    ***   SAVE GLOBAL NORMS   ***
!    *****************************

!    Only processor 0 does output.

     if (rank==0) then

        if (l<10) then
           write(filen,'(i1)') l
        else
           write(filen,'(i2)') l
        end if

!       Save maximum value.

        unit = 1000 + l + Nl*(nvar-1)

        if (filestatus == 'replace') then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_max.tl',form='formatted', &
           status='replace')
           write(unit,*) comment//varname//'_max.tl'
        else if (closefiles) then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_max.tl',form='formatted', &
           status='old',position='append')
        end if

        write(unit,"(2ES25.15e3)") t(l),gmax(l)+tiny

        if ((closefiles).or.(step==Nt)) then
           close(unit)
        end if

!       Save minimum value.

        unit = 2000 + l + Nl*(nvar-1)

        if (filestatus == 'replace') then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_min.tl',form='formatted', &
           status='replace')
           write(unit,*) comment//varname//'_min.tl'
        else if (closefiles) then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_min.tl',form='formatted', &
           status='old',position='append')
        end if

        write(unit,"(2ES25.15e3)") t(l),gmin(l)+tiny

        if ((closefiles).or.(step==Nt)) then
           close(unit)
        end if

!       Save norm 1.

        unit = 3000 + l + Nl*(nvar-1)

        if (filestatus == 'replace') then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_nm1.tl',form='formatted', &
           status='replace')
           write(unit,*) comment//varname//'_nm1.tl'
        else if (closefiles) then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_nm1.tl',form='formatted', &
           status='old',position='append')
        end if

        write(unit,"(2ES25.15e3)") t(l),gnm1(l)+tiny

        if ((closefiles).or.(step==Nt)) then
           close(unit)
        end if

!       Save norm 2.

        unit = 4000 + l + Nl*(nvar-1)

        if (filestatus == 'replace') then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_nm2.tl',form='formatted', &
           status='replace')
           write(unit,*) comment//varname//'_nm2.tl'
        else if (closefiles) then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_nm2.tl',form='formatted', &
           status='old',position='append')
        end if

        write(unit,"(2ES25.15e3)") t(l),gnm2(l)+tiny

        if ((closefiles).or.(step==Nt)) then
           close(unit)
        end if

!       Save total variation.

        unit = 5000 + l + Nl*(nvar-1)

        if (filestatus == 'replace') then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_var.tl',form='formatted', &
           status='replace')
           write(unit,*) comment//varname//'_var.tl'
        else if (closefiles) then
           open(unit,file=trim(outdir)//'/'//varname//trim(filen)//'_var.tl',form='formatted', &
           status='old',position='append')
        end if

        write(unit,"(2ES25.15e3)") t(l),gtvar(l)+tiny

        if ((closefiles).or.(step==Nt)) then
           close(unit)
        end if

     end if

  end do


! ********************************
! ***   SAVE VALUE AT ORIGIN   ***
! ********************************

! Notice that processor 0 does output and also owns the origin.
! Also, we only output the value for the finest grid assuming it
! is the most accurate.

  if ((rank==0).and.((ll<0).or.(ll==Nl-1))) then

!    Save central value.

     unit = 6000 + (nvar-1)

     if (filestatus == 'replace') then
        open(unit,file=trim(outdir)//'/'//varname//'_origin.tl',form='formatted', &
        status='replace')
        write(unit,*) comment//varname//'_origin.tl'
     else if (closefiles) then
        open(unit,file=trim(outdir)//'/'//varname//'_origin.tl',form='formatted', &
        status='old',position='append')
     end if

     if (ghost==1) then        ! Second order interpolation.
        aux = 0.5d0*(savevar(Nl-1,0)+savevar(Nl-1,1))
     else if (ghost>1) then    ! Fourth order interpolation.
        aux = (9.d0*(savevar(Nl-1,0)+savevar(Nl-1,1)) - (savevar(Nl-1,-1)+savevar(Nl-1,2)))/16.d0
     end if

     write(unit,"(2ES25.15e3)") t(Nl-1),aux+tiny

     if ((closefiles).or.(step==Nt)) then
        close(unit)
     end if

  end if


! **************************************
! ***   SAVE VALUE AT EDGE OF GRID   ***
! **************************************

! Notice that processor 0 does output, but the last processor
! on the coarsest grid owns the boundary.
!
! Also, the boundary value lives on the coarsest grid, so
! we only output it if ll<=0.

  if (ll<=0) then

     p = size-1

!    Find asymptotic value at last point and send it to processor 0.

     if (rank==p) then
        aux = savevar(0,Nrl(p))
        if (size>1) then
           call MPI_SEND(aux,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end if
     end if

!    Receive asymptotic value and output it.

     if (rank==0) then

        if (size>1) then
           call MPI_RECV(aux,1,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
        end if

        unit = 7000 + (nvar-1)

        if (filestatus == 'replace') then
           open(unit,file=trim(outdir)//'/'//varname//'_bound.tl',form='formatted', &
           status='replace')
           write(unit,*) comment//varname//'_bound.tl'
        else if (closefiles) then
           open(unit,file=trim(outdir)//'/'//varname//'_bound.tl',form='formatted', &
           status='old',position='append')
        end if

        write(unit,"(2ES25.15e3)") t(0),aux+tiny
 
        if ((closefiles).or.(step==Nt)) then
           close(unit)
        end if

     end if

  end if


! ****************************************
! ***   SAVE VALUE AT SPECIAL POINTS   ***
! ****************************************

! Value at special points at the moment only works
! if we are doing output for all time levels, since
! we need to figure out on which level to do output.

  if (ll<0) then

!    At the moment only for single processor runs.

     if ((output_r1/=0.d0).and.(size==1)) then

!       Find out best grid to interpolate.

        do l=0,Nl-1
           if ((dble(Nrtotal)-0.5d0)*dr(l)>output_r1) then
              bestl = l
           end if
        end do

!       call interpolation routine.

        interpvar => savevar
        aux = interp(bestl,output_r1,.false.)

!       Save value.

        unit = 8000 + (nvar-1)

        if (filestatus == 'replace') then
           open(unit,file=trim(outdir)//'/'//varname//'_r1.tl',form='formatted', &
           status='replace')
           write(unit,*) comment//varname//'_r1.tl'
        else if (closefiles) then
           open(unit,file=trim(outdir)//'/'//varname//'_r1.tl',form='formatted', &
           status='old',position='append')
        end if

        write(unit,"(2ES25.15e3)") t(0),aux+tiny

        if ((closefiles).or.(step==Nt)) then
           close(unit)
        end if

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine save0Dvariable
