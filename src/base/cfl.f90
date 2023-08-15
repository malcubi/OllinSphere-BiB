!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/cfl.f90,v 1.15 2023/03/09 23:29:22 malcubi Exp $

  subroutine cfl(step)

! ****************************
! ***   ADJUST TIME STEP   ***
! ****************************

! This subroutine ajusts the time step in order
! to guarantee that the CFL condition is satisfied.

! Include modules.

  use mpi
  use param
  use arrays
  use procinfo

! Extra variables.

  implicit none

  logical firstcall   ! Is this the first call?

  integer i,l         ! Counters.
  integer step        ! Time step.

  real(8) v_max       ! Maximum speed.
  real(8) aux

  real(8), dimension(0:Nl-1,1-ghost:Nrmax) :: vl,va,vs   ! Propagation speeds.
  real(8), dimension(0:Nl-1) :: vlb,vab                  ! Background propagation speeds (cosmological runs).

  character(20) filestatus

  data firstcall / .true. /


! ***********************
! ***   FIND SPEEDS   ***
! ***********************

! Initialize.

  vl = 0.d0
  va = 0.d0
  vs = 0.d0

  if (cosmic_run) then
     vlb = 0.d0
     vab = 0.d0
  end if

! Loop over levels.

  do l=0,Nl-1

!    Speed of light:  vl = alpha/sqrt(A)/psi^2

     vl(l,:) = abs(alpha(l,:))/sqrt(abs(A(l,:)))/psi2(l,:)

!    Maximal slicing.

     if (slicing=='maximal') then

        va(l,:) = vl(l,:)

!    For Bona-Masso type slicings it is:
!
!    va  =  alpha sqrt[ f(alpha) / A ] / psi^2
!
!    But remember that the code defines:
!
!    falpha  :=  alpha^2 f(alpha)

     else if ((index(slicing,"harmonic")/=0).or. &
              (index(slicing,"1+log")/=0).or. &
              (index(slicing,"shockavoid")/=0).or. &
              (index(slicing,"alphaminus2")/=0)) then

        va(l,:) = sqrt(abs(falpha(l,:)/A(l,:)))/psi2(l,:)

     end if

!    Shift speeds.  For Gammadriver type shifts we have:
!
!    vs  =  sqrt[ (4/3) csi ]

     if ((shift(1:11)=="Gammadriver").and.(drivercsi/=0.d0)) then

        if (.not.cosmic_run) then
           vs(l,:) = sqrt(4.d0/3.d0*drivercsi)
        else
           vs(l,:) = sqrt(4.d0/3.d0*drivercsi)/cosmobg_a(l)
        end if

     end if

!    Add shift contribution.

     if (shift/="none") then
        vl(l,:) = vl(l,:) + abs(beta(l,:))
        va(l,:) = va(l,:) + abs(beta(l,:))
        vs(l,:) = vs(l,:) + abs(beta(l,:))
     end if

!    Cosmological background.

     if (cosmic_run) then

!       Speed of light.

        vlb(l) = cosmobg_alpha(l)/cosmobg_a(l)

!       Slicing gauge speed.

        if (slicing=='maximal') then
           vab(l) = vlb(l)
        else
           vab(l) = sqrt(cosmobg_falpha(l))/cosmobg_a(l)
        end if

     end if

  end do


! ******************************
! ***   FIND MAXIMUM SPEED   ***
! ******************************

  v_max = 0.d0

  do l=0,Nl-1

!    Find maximum speed across grid.

     do i=1-ghost,Nrl(rank)
        if (vl(l,i)>v_max) v_max=vl(l,i)
        if (va(l,i)>v_max) v_max=va(l,i)
        if (vs(l,i)>v_max) v_max=vs(l,i)
     end do

!    Cosmological runs.

     if (cosmic_run) then
        if (vlb(l)>v_max) v_max=vlb(l)
        if (vab(l)>v_max) v_max=vab(l)
     end if

!    For parallel runs find the global maximum.

     if (size>1) then
        call MPI_ALLREDUCE(v_max,aux,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
        v_max = aux
     end if

  end do


! ****************************
! ***   MODIFY TIME STEP   ***
! ****************************

! Calculate the next dt0 using the maximum
! propagation speed.

  if (adjuststepmax) then

!    Here we allow dt0 to increase beyond dtfac*dr0
!    if the speed of ligh is small.

     dt0 = dtfac*dr0/v_max

  else

!    Here we do not allow dt0 to increase beyond dtfac*dr0.
!    This is important in some cases when we need good
!    resolutiion in time even if the speed of light becomes
!    small.

     dt0 = dtfac*dr0*min(1.d0,1.d0/v_max)

  end if

! Adjust time step on different grid levels.

  do l=0,Nl-1
     dt(l) = dt0/2**l
  end do


! *********************************
! ***   SAVE TIME STEP TO FILE  ***
! *********************************

! On first call, replace file. Otherwise just append to it.

  if (firstcall) then
     firstcall = .false.
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Save value of dt0 to file.

  if (rank==0) then

     if (filestatus == 'replace') then
        open(1,file=trim(directory)//'/dt0.tl',form='formatted', &
             status=filestatus)
        write(1,*) '"dt0.tl'
     else
        open(1,file=trim(directory)//'/dt0.tl',form='formatted', &
                status=filestatus,position='append')
     end if
  
     write(1,"(2ES25.15e3)") t(0),dt0
     close(1)

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine cfl

