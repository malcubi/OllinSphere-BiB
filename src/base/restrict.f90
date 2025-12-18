
  subroutine restrict(l,all)

! ********************
! ***   RESTRICT   ***
! ********************

! Restrict data from fine grid into coarse grid.
!
! If the logical variable "all" is set to .true.,
! we restrict all evolving arrays.  Otherwise
! we only restrict the array "restrictvar".

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Declare variables.
  
  implicit none

  logical contains
  logical all

  integer l              ! Grid level.
  integer i,j,k,p        ! Counters.
  integer Ndata          ! Number of data points to send/receive.
  integer imin,jmin,jmax ! Auxiliary quantities.

  integer status(MPI_STATUS_SIZE)

  real(8) r0,interp      ! For interpolation.
  real(8) delta,aux      ! Auxiliary quantities.

  real(8), dimension(1:1) :: auxarr1,auxarr2   ! Auxiliary arrays of dimenson 1 for MPI calls.
  real(8), dimension(0:Nrmax/2) :: w           ! Auxiliary array for efficient communications.


! ************************
! ***   SANITY CHECK   ***
! ************************

! This subroutine should never be called for l=0.

  if (l==0) then
     print *
     print *, 'The subroutine "restrict" should never be called with l=0!'
     print *, 'Aborting! (subroutine restrict)'
     print *
     call die
  end if

! write(*,"(A,I2,A,I2)") ' Restricting level ',l,' to level ',l-1


! ****************************************
! ***   RESTRICT DATA TO COARSE GRID   ***
! ****************************************

! Initialize array w to zero.

  w = 0.d0

! Set Ndata to size of array w.

  Ndata = Nrmax/2 + 1

! Displacement from origin.

  delta = 0.5d0*dr(l)

! I) Send data.

  if (rank==0) then

!    Processor 0 only copies data to itself,
!    so there is no need to use MPI.
!
!    Since here I don't how many variables are evolving,
!    I generate the necessary code automatically at compile
!    time and include it here.  The basic block has the form:
!
!    interpvar => alpha
!    do i=1,Nr-ghost,2
!       r0 = r(l,i) + delta
!       alpha(l-1,i/2+1) = interp(l,r0,.true.)
!    end do
!
!    and is repeated for all evolving variables with the
!    necessary conditionals.

     If (all) then
        include '../auto/restrict_copy.inc'
     else
        interpvar => restrictvar
        do i=1,Nr-(ghost+1),2
           r0 = r(l,i) + delta
           restrictvar(l-1,i/2+1) = interp(l,r0,.true.)
        end do
     end if

  else

!    Processor "rank" on level l must send data
!    to processor int(rank/2) on level (l-1).

     p = rank/2
     ! write(*,'(I3,A,I3,A,I3)') rank,'  sending  ',Ndata,'  points to      ',p

!    Figure out where to start sending data. We send data
!    starting on point i=1 or i=2 depending on how the
!    fine grid is arranged with respect to the coarse grid.
!    This in turn depends on whether we have an even or
!    odd number of points from the origin.

     aux = (r(l,1)-delta)/dr(l)

     if (mod(nint(aux),2)==0) then
        imin = 1
     else
        imin = 2
     end if

!    Send imin and r0. Notice that even though "imax" is an integer
!    I send a REAL8, as otherwise the compiler complains of an
!    argument mismatch. Also, I convert both r0 and imin to
!    1 dimensional arrays first for the same reason.

     r0 = r(l,imin)+delta

     auxarr1(1) = dble(imin)
     auxarr2(1) = r0

     call MPI_SEND(auxarr1,1,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
     call MPI_SEND(auxarr2,1,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)

!    Send data for dynamical arrays.
!
!    Since here I don't how many variables are evolving,
!    I generate the necessary code automatically at compile
!    time and include it here.  The basic block has the form:
!
!    interpvar => alpha
!    do i=imin,Nr-ghost,2
!       r0 = r(l,i) + delta
!       w(i/2) = interp(l,r0,.true.)
!    end do
!    call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
!
!    and is repeated for all evolving variables with the
!    necessary conditionals.

     if (all) then
        include '../auto/restrict_send.inc'
     else
        interpvar => restrictvar
        do i=imin,Nr-(ghost+1),2
           r0 = r(l,i) + delta
           w(i/2) = interp(l,r0,.true.)
        end do
        call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)
     end if

  end if

! II) Receive data.

  if (size>1) then

!    Not all processors receive data, only those with rank
!    such that: rank<size/2.  This is because processors
!    with larger rank don't have any daughter grids.  
!    Processors such that size>2*rank+1 in fact always
!    have two daughters, so must receive data from both.

     if (size>2*rank) then

        if (rank==0) then
!          Processor 0 only receives from the second daugther (j=1).
           jmin = 1
           jmax = 1
        else if (size>2*rank+1) then
!          These processors receive from both daughters (j=0,1).
           jmin = 0
           jmax = 1
        else if (size>2*rank) then
!          These processors receive only from the first daughter (j=0).
           jmin = 0
           jmax = 0
        end if

!       Loop over two possible daughter grids.

        do j=jmin,jmax

!          Processor "rank" on level l-1 receives data
!          from processor 2*rank+j on level l.  Notice
!          that j can take two values (0,1) since a
!          given grid on level l-1 can have one or two 
!          daugthers on level l.

           p = 2*rank+j
           ! write(*,'(I3,A,I3,A,I3)') rank,'  receiving',Ndata,'  points from    ',p

!          Receive imin an r0. Notice that even though "imax" is an integer
!          I receive a REAL8, as otherwise the compiler complains of an
!          argument mismatch.

           call MPI_RECV(auxarr1,1,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(auxarr2,1,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

           imin = nint(auxarr1(1))
           r0 = auxarr2(1)

!          First point displacement.

           k = (1 - imin/2) + nint((r0-r(l-1,1))/dr(l-1))

!          Receive data for dynamical arrays.
!
!          Since here I don't how many variables are evolving,
!          I generate the necessary code automatically at compile
!          time and include it here.  The basic block has the form:
!
!          call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
!          do i=imin,Nrl(p)-ghost,2
!            alpha(l-1,i/2+k) = w(i/2)
!          end do
!
!          and is repeated for all evolving variables with the
!          necessary conditionals.

           if (all) then
              include '../auto/restrict_recv.inc'
           else
              call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)
              do i=imin,Nrl(p)-(ghost+1),2
                 restrictvar(l-1,i/2+k) = w(i/2)
              end do
           end if

        end do

     end if

  end if


! **************************
! ***   FIX SYMMETRIES   ***
! **************************

! Apply symmetries at origin.

  if (rank==0) then
     call symmetries(l)
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine restrict
