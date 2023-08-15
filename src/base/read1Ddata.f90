!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/read1Ddata.f90,v 1.4 2021/02/18 20:31:10 malcubi Exp $

  subroutine read1Ddata(varname,outdir)

! *************************
! ***   READ 2D FILES   ***
! *************************

! Ths subroutine reads 2D files. Notice that it
! only reads files on a single grid level "l"
! that is received as input.
!
! The routine receives as parameters the name
! of the file to grid function to be read,
! the output directory, and the grid level l.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical contains

  integer i,l,p           ! Counters.

  real(8) rdata           ! Position data.
  real(8) tdata           ! Time data.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: var

  character(1) first
  character(2) filen
  character(50) varname,outdir,line


! **********************************
! ***   PROCESSOR 0 READS DATA   ***
! **********************************

! Only processor 0 reads the data and later
! distributes it to the other processors.

  if (rank==0) then

!    *********************************
!    ***   LOOP OVER GRID LEVELS   ***
!    *********************************

     do l=0,Nl-1

!       Append grid level to file name.

        if (l<10) then
           write(filen,'(i1)') l
        else
           write(filen,'(i2)') l
        end if


!       *********************
!       ***   READ FILE   ***
!       *********************

!       Open file.

        open(1,file=trim(outdir)//'/'//trim(varname)//trim(filen)//'.rl',form='formatted',status='old')

!       Read first line and check if it is a comment.
!       If it is, then read the current time from it.

        read(1,'(A)') line

        line = adjustl(line)
        first = line(1:1)

        if (first=="#") then
           if (contains(line,"#Time")) then
              i = index(line,'=')
              line = trim(adjustl(line(i+1:len_trim(line))))
              read(line,*) tdata
              t = tdata
           end if
        else
           rewind(1)
        end if

!       Read data.

        do i=1-ghost,Nrtotal
           read(1,*) rdata,var(l,i)
        end do

!       Close file.

        close(1)

     end do

  end if


! ********************************************
! ***   DISTRIBUTE DATA AMONG PROCESSORS   ***
! ********************************************

! For parallel runs we need to distribute the
! data among the different processors.

  if (size==1) then
     savevar = var
  else
     call distribute(0,Nl-1,savevar,var)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine read1Ddata
