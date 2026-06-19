
  subroutine dissipation(l,sym,diss)

! This routine adds Kreiss-Oliger dissipation to the sources.

! **************************************
! ***   ADD DISSIPATION TO SOURCES   ***
! **************************************

! Include modules.

  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer i,l,sym,imax

  real(8) idt,diss


! *******************
! ***   NUMBERS   ***
! *******************

  idt = 1.d0/dt(l)


! ************************
! ***   SANITY CHECK   ***
! ************************

  if (abs(sym)/=1) then
     print *
     print *, 'dissipation: Symmetry parameter must be +-1'
     print *
     call die
  end if

! For the coarsest grid we apply dissipation as far as possible,
! but for finer grids we stop further away since otherwise we
! can mess up the interpolation from the coarse grid.

  if (l==0) then
     imax = Nr - ghost
  else
     imax = Nr - 2*ghost
  end if


! ************************
! ***   SECOND ORDER   ***
! ************************

  if (order=="two") then

!    Interior points:  Second order evolution
!    requires fourth order dissipation.

     do i=1,imax
        sourcevar(l,i) = sourcevar(l,i) - diss*idt*(6.d0*dissipvar(l,i) &
                       - 4.d0*(dissipvar(l,i+1) + dissipvar(l,i-1)) &
                       +      (dissipvar(l,i+2) + dissipvar(l,i-2)))/6.d0
     end do


! ************************
! ***   FOURTH ORDER   ***
! ************************

  else if (order=="four") then

!    Interior points:  Fourth order evolution
!    requires sixth order dissipation.

     do i=1,imax
        sourcevar(l,i) = sourcevar(l,i) - diss*idt*(20.d0*dissipvar(l,i) &
                       - 15.d0*(dissipvar(l,i+1) + dissipvar(l,i-1)) &
                       +  6.d0*(dissipvar(l,i+2) + dissipvar(l,i-2)) &
                       -       (dissipvar(l,i+3) + dissipvar(l,i-3)))/20.d0
     end do


! ***********************
! ***   SIXTH ORDER   ***
! ***********************

  else if (order=="six") then

!    Interior points:  Sixth order evolution
!    requires eighth order dissipation.

     do i=1,imax
        sourcevar(l,i) = sourcevar(l,i) - diss*idt*(70.d0*dissipvar(l,i) &
                       - 56.d0*(dissipvar(l,i+1) + dissipvar(l,i-1)) &
                       + 28.d0*(dissipvar(l,i+2) + dissipvar(l,i-2)) &
                       -  8.d0*(dissipvar(l,i+3) + dissipvar(l,i-3)) &
                       +       (dissipvar(l,i+4) + dissipvar(l,i-4)))/70.d0
     end do


! ************************
! ***   EIGHTH ORDER   ***
! ************************

  else if (order=="eight") then

!    Interior points:  Eighth order evolution
!    requires tenth order dissipation.

     do i=1,imax
        sourcevar(l,i) = sourcevar(l,i) - diss*idt*(252.d0*dissipvar(l,i) &
                       - 210.d0*(dissipvar(l,i+1) + dissipvar(l,i-1)) &
                       + 120.d0*(dissipvar(l,i+2) + dissipvar(l,i-2)) &
                       -  45.d0*(dissipvar(l,i+3) + dissipvar(l,i-3)) &
                       +  10.d0*(dissipvar(l,i+4) + dissipvar(l,i-4)) &
                       -        (dissipvar(l,i+5) + dissipvar(l,i-5)))/252.d0
     end do

  end if


! ********************************
! ***   SYMMETRIES AT ORIGIN   ***
! ********************************

  if (rank==0) then
     do i=1,ghost
        sourcevar(l,1-i) = dble(sym)*sourcevar(l,i)
     end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine dissipation


