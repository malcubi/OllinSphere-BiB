!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/restrictintegral.f90,v 1.2 2026/01/05 19:24:10 malcubi Exp $

  subroutine restrictintegral

! For the case of several grid levels we need to restrict integrated
! quantities from the fine to the coarse grids.
!
! But first we need to substract constant difference at edge of
! the fine grid. Since the integrals are done independently at each
! grid level, we find that there will be small jumps due to accumulated
! mnumerical error when we pass from the end of a fine grid to the
! next coarse grid, Here we correct for that.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer i,l

  real(8) r0,m0,m1,delta
  real(8) interp
  real(8) aux


! ************************************************
! ***   FIND JUMPS FROM FINE TO COARSE GRIDS   ***
! ************************************************

  restrictvar => intvar
  interpvar   => intvar

! Loop over all levels, from fine to coarse.

  do l=Nl-1,1,-1

!    Find edge of fine grid, and evaluate the
!    function "intvar" there.
!
!    For parallel runs broadcast the value of (r0,m0)
!    from the last processor (which has the true edge)
!    to all others. 

     r0 = r(l,Nr) - 0.5d0*dr(l)
     m0 = 0.5d0*(intvar(l,Nr) + intvar(l,Nr-1))

     if (size>1) then
        call MPI_BCAST(r0,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(m0,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)
     end if

!    Interpolate the value of "intvar" in the coarse
!    grid at r0 and call it m1.

     m1 = interp(l-1,r0,.false.)

!    For parallel runs we add contributions from all
!    processors.  Since only 1 processor should own
!    the point r0, this ensures that now all have the
!    correct value of m1.

     if (size>1) then
       call MPI_Allreduce(m1,aux,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
       m1 = aux
     end if

!    Find difference bewteen m1 and m0.

     delta = m1 - m0

!    Subtract delta.

     intvar(l-1,:) = intvar(l-1,:) - delta

!    Restrict.

     call restrict(l,.false.)

!    Fix symmetries.

     if (rank==0) then
        do i=1,ghost
           intvar(l-1,1-i) = intvar(l-1,i)
        end do
     end if

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine restrictintegral
