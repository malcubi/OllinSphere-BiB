!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/geometry/horizon_finder.f90,v 1.21 2023/05/17 22:55:19 malcubi Exp $

  subroutine horizon_finder

! **************************
! ***   HORIZON FINDER   ***
! **************************

! This routine finds the apparent horizon.
!
! The expansion of outgoing null lines is given by:
!
! theta = (2/r + d B / B + 4 d phi) / (psi**2 sqrt(A)) - 2 K_theta^theta
!                 r           r

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  logical contains

  logical firstcall        ! Is this the first call?

  logical, dimension (0:Nl-1) :: flag    ! Have we found a horizon on level l?
 
  integer i,imax           ! Counters.
  integer l                ! Refinement level.
  integer p                ! Processor counter.
  integer Naux
  integer dieflag

  integer status(MPI_STATUS_SIZE)

  integer, dimension (0:Nl-1) :: ah_i   ! Grid point to left of horizon.

  real(8) theta,theta_old  ! Expansion of null geodesics.
  real(8) r0               ! Interpolation.
  real(8) smallpi

  real(8), dimension (0:Nl-1) :: ah_r                     ! Horizon position.
  real(8), dimension (0:Nl-1) :: ah_A,ah_B                ! Metric coefficients A and B at horizon.
  real(8), dimension (0:Nl-1) :: ah_psi,ah_alpha          ! Conformal factor and lapse at horizon.
  real(8), dimension (0:Nl-1) :: ah_eQ                    ! Electirc charge at horizon.
  real(8), dimension (0:Nl-1) :: ah_schw,ah_area,ah_mass  ! Schwarzschild radius, area and mass of horizon.

  real(8), dimension (0:Nl-1,1:9) :: hdata,ndata          ! For data transfer between processors.

  character(5)  filen
  character(20) filestatus

  data firstcall / .true. /


! *******************
! ***   NUMBERS   ***
! *******************

  smallpi = acos(-1.d0)


! **********************
! ***   INITIALIZE   ***
! **********************

  ah_i = 0
  ah_r = 0.d0

  ah_psi   = 0.d0
  ah_A     = 0.d0
  ah_B     = 0.d0
  ah_alpha = 0.d0

  ah_schw = 0.d0
  ah_area = 0.d0
  ah_mass = 0.d0

  if (contains(mattertype,"electric")) then
     ah_eQ = 0.0
  end if

  hdata = 0.d0
  ndata = 0.d0


! ***************************************
! ***   FIND HORIZON ON COARSE GRID   ***
! ***************************************

! Initialize flag for horizon.

  flag = .false.

! Initialize trial horizon position at outermost point.
! Notice that all processors look for a horizon.

  if (rank==size-1) then
     imax = Nr
  else
     imax = Nr-ghost+1
  end if

! Loop over grid levels.

  do l=0,Nl-1

!    Start moving inward from boundary.

     i = imax

     do while ((.not.flag(l)).and.(i>=1))

!       Find expansion.

        theta = (2.d0/r(l,i) + D1_B(l,i)/B(l,i) &
              + 4.d0*D1_phi(l,i))/sqrt(A(l,i))/psi(l,i)**2 - 2.d0*KBPHYS(l,i)

!       If the expansion is already negative at the outermost point
!       then the horizon is outside so we just jump out of here.

        if ((theta<0.d0).and.(i==imax)) goto 100

!       If the expansion changes sign then we have bracketed the horizon.

        if (theta<0.d0) then

!          Set horizon flag and grid point to left of horizon.

           flag(l) = .true.
           ah_i(l) = i

!          Interpolate position of the horizon (linearly).

           ah_r(l) = r(l,i+1) - (r(l,i+1)-r(l,i))/(theta_old-theta)*theta_old
           hdata(l,1) = ah_r(l)

        end if

!       Save value of theta.

        theta_old = theta

!       Move one grid point in.

        i = i - 1

     end do

     if (flag(Nl-1)) ahfound = .true.

     100 continue

  end do


! ***********************************
! ***   FIND HORIZON PROPERTIES   ***
! ***********************************

! Loop over grid levels.

  do l=0,Nl-1

!    If we found a horizon calculate its properties.
!    Here we use linear interpolation.

     if (flag(l)) then

!       Interpolating point.

        i  = ah_i(l)
        r0 = ah_r(l)

!       Sanity check.

        if ((r0<r(l,i)).or.(r0>r(l,i+1))) then
           print *
           print *, 'horizon_finder.f90:  This should never happen'
           print *, 'Aborting! (subroutine horizon_finder)'
           print *
           call die
        end if

!       Interpolate conformal factor at the horizon (required
!       to calculate the radius in Schwarzschild coordinates).

        ah_psi(l)  = psi(l,i) + (r0 - r(l,i))*(psi(l,i+1)-psi(l,i))/dr(l)
        hdata(l,2) = ah_psi(l)

!       Interpolate metric function A at the horizon.

        ah_A(l)    = A(l,i) + (r0 - r(l,i))*(A(l,i+1)-A(l,i))/dr(l)
        hdata(l,3) = ah_A(l)

!       Interpolate metric function B at the horizon
!       (required to calculate the radius in Schwarzschild coordinates).

        ah_B(l)    = B(l,i) + (r0 - r(l,i))*(B(l,i+1)-B(l,i))/dr(l)
        hdata(l,4) = ah_B(l)

!       Interpolate lapse at the horizon.

        ah_alpha(l) = alpha(l,i) + (r0 - r(l,i))*(alpha(l,i+1)-alpha(l,i))/dr(l)
        hdata(l,5)  = ah_alpha(l)

!       Calculate the radius of the horizon in Schwarzschild coordinates:
!
!       r_schw  =  r  psi**2 sqrt(B)

        ah_schw(l) = ah_r(l)*ah_psi(l)**2*sqrt(ah_B(l))
        hdata(l,6) = ah_schw(l)

!       Calculate the area of the horizon:
!
!       area  =  4 pi r_schw**2  = 4 pi r**2 psi**4 B

        ah_area(l) = 4.d0*smallpi*(ah_r(l)*ah_psi(l)**2)**2*ah_B(l)
        hdata(l,7) = ah_area(l)

!       Calculate the irreducible mass of the horizon.
!
!       mass  =  sqrt( area / (16 pi) )

        ah_mass(l) = sqrt(ah_area(l)/16.d0/smallpi)
        hdata(l,8) = ah_mass(l)

!       Interpolate electric charge at horizon if required.

        if (contains(mattertype,"electric")) then
           ah_eQ(l) = eQ_surf(l,i) + (r0 - r(l,i))*(eQ_surf(l,i+1)-eQ_surf(l,i))/dr(l)
           hdata(l,9) = ah_eQ(l)
        end if

     end if

  end do


! ********************************************
! ***   COMMUNICATION BETWEEN PROCESSORS   ***
! ********************************************

! Loop over grid levels.

  do l=0,Nl-1

!    Now figure out which processor found a horizon
!    and send the information to processor 0.

     if (size>1) then

        flag(l) = .false.

        Naux = 8

!       Processor 0.

        if (rank==0) then

           if (hdata(l,1)/=0.d0) flag(l) = .true.

!          Receive data from other processors.

           do p=1,size-1

              call MPI_RECV(ndata(l,:),Naux,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)

!             If a horizon was found by processor p copy the data.

              dieflag = 0

              if (ndata(l,1)/=0.d0) then
                 if (.not.flag(l)) then
                    flag(l) = .true.
                    ah_r(l)     = ndata(l,1)
                    ah_psi(l)   = ndata(l,2)
                    ah_A(l)     = ndata(l,3)
                    ah_B(l)     = ndata(l,4)
                    ah_alpha(l) = ndata(l,5)
                    ah_schw(l)  = ndata(l,6)
                    ah_area(l)  = ndata(l,7)
                    ah_mass(l)  = ndata(l,8)
                    ah_eQ(l)    = ndata(l,9)
                 else
                    print *
                    print *, 'Two different processors found a horizon on level ',l
                    print *, 'This should never happen.'
                    print *, 'Aborting! (subroutine horizon_finder)'
                    dieflag = 1
                 end if
              end if

           end do

!       Other processors send data to processor 0.

        else

           call MPI_SEND(hdata(l,:),Naux,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)

        end if

!       Check if we want to abort the run.

        call MPI_BCAST(dieflag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if (dieflag/=0) call die

     end if

  end do

! Store the ah mass.

  ahmass = ah_mass(Nl-1)


! *****************************
! ***   SAVE HORIZON DATA   ***
! *****************************

! At the moment I output separately the horizon properties for
! each grid level.  I should probably fix this later to output
! only the horizon from the finer grid where it was found.

! On first call, replace file. Otherwise just append to it.

  if (firstcall) then
     firstcall = .false.
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Save horizon data for all grid levels.

  do l=0,Nl-1

     if (l<10) then
        write(filen,'(i1)') l
     else
        write(filen,'(i2)') l
     end if

     if (rank==0) then

!       Save apparent horizon position.

        if (filestatus == 'replace') then
           open(1,file=trim(directory)//'/ah_rad'//trim(filen)//'.tl',form='formatted', &
                status=filestatus)
           write(1,*) '"ah_rad.tl'
        else
           open(1,file=trim(directory)//'/ah_rad'//trim(filen)//'.tl',form='formatted', &
                status=filestatus,position='append')
        end if
  
        if (isnan(ah_r(l))) then
           print *
           print *, 'Found NaN in horizon finder.  Aborting!'
           call die
        end if

        write(1,"(2ES14.6)") t(0),ah_r(l)
        close(1)

!       Save apparent horizon position in Schwarzschild coordinates.

        if (filestatus == 'replace') then
           open(1,file=trim(directory)//'/ah_schw'//trim(filen)//'.tl',form='formatted', &
                status=filestatus)
           write(1,*) '"ah_schw.tl'
        else
           open(1,file=trim(directory)//'/ah_schw'//trim(filen)//'.tl',form='formatted', &
                status=filestatus,position='append')
        end if
  
        write(1,"(2ES14.6)") t(0),ah_schw(l)
        close(1)

!       Save lapse at the apparent horizon position.

        if (filestatus == 'replace') then
           open(1,file=trim(directory)//'/ah_alpha'//trim(filen)//'.tl',form='formatted', &
                status=filestatus)
           write(1,*) '"ah_alpha.tl'
        else
           open(1,file=trim(directory)//'/ah_alpha'//trim(filen)//'.tl',form='formatted', &
                status=filestatus,position='append')
        end if
  
        write(1,"(2ES14.6)") t(0),ah_alpha(l)
        close(1)

!       Save area of the apparent horizon.

        if (filestatus == 'replace') then
           open(1,file=trim(directory)//'/ah_area'//trim(filen)//'.tl',form='formatted', &
                status=filestatus)
           write(1,*) '"ah_area.tl'
        else
           open(1,file=trim(directory)//'/ah_area'//trim(filen)//'.tl',form='formatted', &
                status=filestatus,position='append')
        end if
  
        write(1,"(2ES14.6)") t(0),ah_area(l)
        close(1)

!       Save mass of the apparent horizon.

        if (filestatus == 'replace') then
           open(1,file=trim(directory)//'/ah_mass'//trim(filen)//'.tl',form='formatted', &
                status=filestatus)
           write(1,*) '"ah_mass.tl'
        else
           open(1,file=trim(directory)//'/ah_mass'//trim(filen)//'.tl',form='formatted', &
                status=filestatus,position='append')
        end if
  
        write(1,"(2ES14.6)") t(0),ah_mass(l)
        close(1)

!       Save electric charge at the apparent horizon.

        if (contains(mattertype,"electric")) then

           if (filestatus == 'replace') then
              open(1,file=trim(directory)//'/ah_eQ'//trim(filen)//'.tl',form='formatted', &
                   status=filestatus)
              write(1,*) '"ah_eQ.tl'
           else
              open(1,file=trim(directory)//'/ah_eQ'//trim(filen)//'.tl',form='formatted', &
                   status=filestatus,position='append')
           end if
  
           write(1,"(2ES14.6)") t(0),ah_eQ(l)
           close(1)

        end if

     end if

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine horizon_finder
