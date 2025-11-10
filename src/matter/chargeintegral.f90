!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/matter/chargeintegral.f90,v 1.6 2025/10/01 17:51:29 malcubi Exp $

  subroutine chargeintegral

! *****************************************************
! ***   CALCULATION OF INTEGRATED ELECTRIC CHARGE   ***
! *****************************************************

! This is the integrated charge:
!
!       /                     /            1/2      6   2
! Q  =  | echarge dV  =  4 pi | echarge [ A   B  psi ] r dr
!       /                     /
!
! Notice that this assumes that the spacetime is REGULAR at the
! origin, so it will not work for an eternal black holes such
! as Schwarzschild or Reissner-Nordstrom.
!
! NOTE: The contribution from the complex scalar field is just Q = q Nb,
! with q the fundamental charge. So if the only electromagnetic source
! is a complex scalar field it is redundant to make both calculations.
! (the two routines are identical)

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo
  use integrals

! Extra variables.

  implicit none

  integer i,l,i0

  real(8) r0,interp
  real(8) delta,m0,m1
  real(8) R99,eQTOT
  real(8) aux,smallpi


! *******************
! ***   NUMBERS   ***
! *******************

  smallpi = acos(-1.d0)


! **********************
! ***   INITIALIZE   ***
! **********************

  eQ_int = 0.d0


! *********************
! ***   INTEGRATE   ***
! *********************

  auxarray = 4.d0*smallpi*echarge*dsqrt(A)*B*psi**6*r**2

  intvar => auxarray

  do l=0,Nl-1
     eQ_int(l,:) = integral(l)
  end do

! Restrict integral.

  if (Nl>1) then
     intvar => eQ_int
     call restrictintegral
  end if


! ************************************************
! ***   OUTPUT TOTAL ELECTRIC CHARGE AND R99   ***
! ************************************************

! Calculate the radius for which we have 99% of the
! total charge. We do this only on the coarse grid,
! and use linear interpolation.
!
! Notice that at the moment we only do this calculation at
! t=0 (for the initial data).  I might change this later.

  if (t(0)==0.d0) then

     R99 = 0.d0
     eQTOT = eQ_int(0,Nr)

!    Single processor run.

     if (size==1) then

       do i=1,Nr
           if ((abs(eQ_int(0,i-1))<0.99d0*abs(eQTOT)).and.(abs(eQ_int(0,i))>=0.99d0*abs(eQTOT))) then
              R99 = r(0,i-1) + dr(0)*(0.99d0*eQTOT-eQ_int(0,i-1))/(eQ_int(0,i)-eQ_int(0,i-1))
           end if
        end do

!    Multiple processor run. Broadcast the charge number from
!    processor size-1 (the boundary) to all other processors.

     else

        call MPI_BCAST(eQTOT,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)

!       Now each processor tries to find R99.

        if (rank==0) then
           i0 = 1
        else
           i0 = 1+ghost
        end if

        do i=i0,Nr
           if ((eQ_int(0,i-1)<0.99d0*eQTOT).and.(eQ_int(0,i)>=0.99d0*eQTOT)) then
              R99 = r(0,i-1) + dr(0)*(0.99d0*eQTOT-eQ_int(0,i-1))/(eQ_int(0,i)-eQ_int(0,i-1))
           end if
        end do

!       We assume that only one processor found R99 (not necessarily the same one).
!       We then find the maximum values of R99 across all processors.
!
!       Notice that if the integrated charge grows monotonically, then only one processor must
!       have found the a value for R99, for the others that value should remain zero.
!       But, if for some reason the integrated charge is not monotone (numerical error, there
!       is a black hole and the integrated charge makes no sense, we have a non-minimally coupled
!       scalar field, etc.), then we keep the largest value.

        call MPI_Allreduce(R99,aux,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
        R99 = aux

     end if

!    Processor 0 writes result to screen.

     if (rank==0) then
        write(*,'(A,ES23.16)') ' Total electric charge Q        = ',eQTOT
        write(*,'(A,ES23.16)') ' Effective radius R99 from Q(r) = ',R99
        print *
     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine chargeintegral
