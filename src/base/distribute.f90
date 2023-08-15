!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/distribute.f90,v 1.4 2019/05/08 15:39:42 malcubi Exp $

  subroutine distribute(lmin,lmax,varl,varg)

! *******************************************************
! ***   DISTRIBUTE ONE LARGE ARRAY AMONG PROCESSORS   ***
! *******************************************************

! This is a small routine that only makes sense for multi-processor
! runs.  It distributes a large (global) array "varg" of dimension Nrtotal
! owned by processor 0 among local arrays "varl" on all processors.
!
! This is useful for example for initial data that was calculated
! or read only by processor 0.
!
! Notice that we might not want to distribute all grid levels.
! That is the purpose of the parameters "lmin" and "lmax".

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer i,i0,imax,l,p,Naux
  integer lmin,lmax
  integer status(MPI_STATUS_SIZE)

  real(8), dimension (1-ghost:Nrmax)   :: vara
  real(8), dimension (0:Nl-1,1-ghost:Nrmax),target   :: varl   ! Local array.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal)        :: varg   ! Global array.


! ***********************
! ***   PROCESSOR 0   ***
! ***********************

  vara = 0.d0

  Naux = Nrmax + ghost

! Processor 0 is assumed to own the original large array.

  if (rank==0) then

!    Set up local boundaries.

     if (size==1) then
        imax = Nrl(0)
     else
        imax = Nrl(0) - ghost
     end if

!    Copy local data for processor 0.

     do l=lmin,lmax
        do i=1-ghost,imax
           varl(l,i) = varg(l,i)
        end do
     end do

!    Send corresponding data to other processors.

     i0 = imax

     do p=1,size-1

!       Set up local boundaries.

        if (p==size-1) then
           imax = Nrl(p)
        else
           imax = Nrl(p) - ghost
        end if

!       Loop over grid levels.

        do l=lmin,lmax

!          Copy local data for processor p and send it.

           do i=1,imax
              vara(i) = varg(l,i+i0)
           end do

           call MPI_SEND(vara,Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

        end do

!       Increment i0.

        i0 = i0 + imax

     end do


! ****************************
! ***   OTHER PROCESSORS   ***
! ****************************

! Other processor just receive data from processor 0.

  else

     do l=lmin,lmax

        call MPI_RECV(vara,Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,status,ierr)

        do i=1,Nr
           varl(l,i) = vara(i)
        end do

     end do

  end if


! ***********************
! ***   SYNCHRONIZE   ***
! ***********************

  do l=lmin,lmax
     syncvar => varl(l,:)
     call sync
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine distribute

