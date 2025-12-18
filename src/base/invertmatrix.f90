!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/invertmatrix.f90,v 1.25 2025/10/13 18:45:14 malcubi Exp $

  subroutine invertmatrix(lmin,lmax,u0,u,C0,C1,C2,sym,bound)

! **********************************************
! ***   SOLVE LINEAR POISSON-TYPE EQUATION   ***
! **********************************************

! This subroutine inverts a band-diagonal matrix that
! solves a linear Poisson-type equation of the form:
!
!  2
! d u  +  C0(r) d u  + C1(r) u  =  C2(r)
!  r             r
!
! Notice that we assume that the coefficient of the second
! derivative has already been normalized to 1.
!
! The outer boundary condition depends on the value of the
! parameter "bound". The possible values are:
!
!       robin      (u0 is the value at infinity)
!       dirichlet  (u0 is the boundary value)
!       newmann    (u0 is the derivative at the boundary)
!
! Notice that we might not want to solve the equation on all
! grid levels. That is the purpose of the parameters "lmin"
! and "lmax".
!
! For the solutions we use the routines  "bandec" and "banbks"
! from NUMERICAL RECIPES.  They can be found at the end of this
! file, with very minor changes to suit the situation here.
! They are used one after the other to invert the band-diagonal
! matrix.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  integer i,l,lmin,lmax,NP
  integer sym
  integer i0,iaux

  real(8) u0,ub
  real(8) r0,aux,delta
  real(8) interp           ! Interpolation function.

  real(8), dimension (0:Nl-1,1-ghost:Nrtotal), target :: u
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: C0,C1,C2

  character bound*(*),bd*10

! Arrays for band diagonal matrix inversion.

  integer, dimension (0:Nrtotal) :: indx
  real(8), dimension (0:Nrtotal,-2:2) :: M     ! Matrix elements defining the system.
  real(8), dimension (0:Nrtotal,-2:2) :: ML    ! Lower diagonal matrix returned from bandec.
  real(8), dimension (0:Nrtotal) :: S0         ! Source vector.


! ************************
! ***   SANITY CHECK   ***
! ************************

  if (abs(sym)/=1) then
     print *
     print *, 'Symmetry parameter must be +-1'
     print *, 'Aborting! (subroutine invertmatrix)'
     print *
     call die
  end if

  if (lmax<lmin) then
     print *
     print *, 'We must have lmax >= lmin'
     print *, 'Aborting! (subroutine invertmatrix)'
     print *
     call die
  end if


! *********************************
! ***   LOOP OVER GRID LEVELS   ***
! *********************************

  do l=lmin,lmax


!    **********************
!    ***   INITIALIZE   ***
!    **********************

!    Initialize to zero.

     M  = 0.d0
     ML = 0.d0
     S0 = 0.d0


!    ***********************
!    ***   SOURCE TERM   ***
!    ***********************

!    Remember that we multiply everything with dr**2,
!    so the source term for the matrix inversion is
!    actually dr**2 times the source of the differential
!    equation.

     do i=1,Nrtotal-1
        S0(i) = dr(l)**2*C2(l,i)
     end do


!    ******************
!    ***   ORIGIN   ***
!    ******************

!    At left boundary we impose symmetry conditions.

     M(0,0) = 1.d0
     M(0,1) = - dble(sym)


!    ************************
!    ***   SECOND ORDER   ***
!    ************************

     if (order=="two") then

!       Interior points.

        do i=1,Nrtotal-1
           M(i,-1) = 1.d0 - 0.5d0*C0(l,i)*dr(l)
           M(i, 0) = - 2.d0 + C1(l,i)*dr(l)**2
           M(i,+1) = 1.d0 + 0.5d0*C0(l,i)*dr(l)
        end do


!    ************************
!    ***   FOURTH ORDER   ***
!    ************************

     else

!       Point i=1.  Here we use the symmetries to avoid making
!       reference to the point at i=-1.

        i = 1

        M(i,-1) = + 16.d0 - 8.d0*C0(l,i)*dr(l)
        M(i, 0) = - 30.d0 + 12.d0*C1(l,i)*dr(l)**2
        M(i,+1) = + 16.d0 + 8.d0*C0(l,i)*dr(l) - dble(sym)*(1.d0 - C0(l,i)*dr(l))
        M(i,+2) = - 1.d0 - C0(l,i)*dr(l)

        S0(i) = 12.d0*S0(i)

!       Interior points.

        do i=2,Nrtotal-2

           M(i,-2) = - 1.d0 + C0(l,i)*dr(l)
           M(i,-1) = + 16.d0 - 8.d0*C0(l,i)*dr(l)
           M(i, 0) = - 30.d0 + 12.d0*C1(l,i)*dr(l)**2
           M(i,+1) = + 16.d0 + 8.d0*C0(l,i)*dr(l)
           M(i,+2) = - 1.d0 - C0(l,i)*dr(l)

           S0(i) = 12.d0*S0(i)

        end do

!       Point i = Nrtotal-1. This is only second order!

        i = Nrtotal - 1

        M(i,-1) = 1.d0 - 0.5d0*C0(l,i)*dr(l)
        M(i, 0) = - 2.d0 + C1(l,i)*dr(l)**2
        M(i,+1) = 1.d0 + 0.5d0*C0(l,i)*dr(l)

     end if


!    ******************************
!    ***   BOUNDARY CONDITION   ***
!    ******************************

!    On the coarsest grid the boundary conditions
!    are those provided when the routine was called.

     if (l==0) then

        bd = bound
        ub = u0

!    On finer grids we always use Dirichlet boundary
!    conditions, with the boundary value interpolated
!    from the coarser grid.

     else

        bd = "dirichlet"

        r0 = r(l,Nrtotal)
        interpvar => u
        ub = interp(l-1,r0,.false.)

     end if

!    Apply boundary condition.

     if (bd=="robin") then

!       Robin boundary condition:
!
!       d u  =  (u0 - u) / r
!        r
!
!       Notice that on multiprocessor runs this routine is
!       only called by processor 0, and processor 0 does
!       not own the outer boundary, so we must calculate
!       its position here.
!
!       This is only second order at the moment. The value
!       of aux below corresponds to:
!
!       aux  =  dr/(r(Nrtotal) + r(Nrtotal-1))
!
!       For a standard radial coordinate this turns out to
!       be given by:
!
!       aux = 1/(2*(Nrtotal-1))
!
!       For a transformed radial coordinate, we need to use
!       the boundary position saved in "rtmax", and the
!       asymptotic derivative "drinf".  But this assumes
!       that the transformation is linear far away.

        if (newr) then
           aux = dr(0)/(2.d0*rtmax - dr(0)*drinf)
        else
           aux = 0.5d0/dble(Nrtotal-1)
        end if

        M(Nrtotal,-1) = - 1.d0 + aux
        M(Nrtotal, 0) = + 1.d0 + aux

        S0(Nrtotal) = 2.d0*aux*ub

     else if (bd=="dirichlet") then

!       Dirichlet boundary condition:
!
!       u = u0

        M(Nrtotal,0) = 1.d0

        S0(Nrtotal) = ub

     else if (bd=="newmann") then

!       Newmann boundary condition:
!
!       d u  =  u0
!        r
!
!       This is only second order at the moment!

        M(Nrtotal,-1) = - 1.d0
        M(Nrtotal, 0) = + 1.d0

        S0(Nrtotal) = ub*dr(l)

     else

!       Unknown boundary condition.

        print *
        print *,'Unknown boundary condition.'
        print *,'Aborting! (subroutine invertmatrix)'
        print *
        call die

     end if


!    *********************************
!    ***   CALL MATRIX INVERSION   ***
!    *********************************

!    Call matrix inversion. On output the
!    solution is contained in S0.

     NP = Nrtotal + 1

     call bandec(M,NP,2,2,ML,indx)
     call banbks(M,NP,2,2,ML,indx,S0)

!    Copy solution.

     do i=0,Nrtotal
        u(l,i) = S0(i)
     end do

!    Symmetries at origin.

     do i=1,ghost
        u(l,1-i) = sym*u(l,i)
     end do

  end do


! ********************
! ***   RESTRICT   ***
! ********************

! Restrict solution from fine to coarse grid.
! We don't call the subroutine "restrict"
! since here we are running only on processor 0.
! We use cubic interpolation.

  do l=lmax,lmin+1,-1

     interpvar => u

     do i=1,Nr-ghost,2
        r0 = r(l,i) + 0.5d0*dr(l)
        u(l-1,i/2+1) = interp(l,r0,.true.)
     end do

!    Fix symmetries.

     do i=1,ghost
        u(l-1,1-i) = sym*u(l-1,i)
     end do

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine invertmatrix








  subroutine bandec(a,n,m1,m2,al,indx)

! This routine performs an LU decomposition of a
! band diagonal matrix.
!
! On input, a(n,m1+m2+1) contains the matrix elements
! defining the system.  Here "n" is the number of
! equations, and (m1,m2) indicate how many subdiagonals
! and superdiagonals are different from zero (e.g. for
! a tridiagonal matrix m1=m2=1, and for a pentadiagonal
! (m1=m2=2).
!
! On output, "a" now contains the upper triangular matrix
! and "al" the lower tringular matrix.  The integer vector "indx"
! records the row permutation effected by the partial pivoting
! and is required by subroutine "banbks".

  implicit none

  integer n,m1,m2
  integer i,j,k,l,mm

  integer indx(n)

  real(8) a(n,m1+m2+1),al(n,m1+m2+1)
  real(8) d,dum,TINY

  parameter (TINY=1.d-20)

  mm=m1+m2+1
  l=m1

  do i=1,m1
     do j=m1+2-i,mm
        a(i,j-l)=a(i,j)
     end do
     l=l-1
     do j=mm-l,mm
        a(i,j)=0.d0
     end do
  end do
 
  d=1.d0
  l=m1

  do k=1,n

    dum=a(k,1)
    i=k
 
    if (l.lt.n) l=l+1

    do j=k+1,l
       if (abs(a(j,1)).gt.abs(dum)) then
          dum=a(j,1)
          i=j
       end if
    end do

    indx(k)=i

    if (dum.eq.0.d0) a(k,1)=TINY

    if (i.ne.k) then
       d=-d
       do j=1,mm
          dum=a(k,j)
          a(k,j)=a(i,j)
          a(i,j)=dum
       end do
    end if

    do i=k+1,l
       dum=a(i,1)/a(k,1)
       al(k,i-k)=dum
       do j=2,mm
          a(i,j-1)=a(i,j)-dum*a(k,j)
       end do
       a(i,mm)=0.d0
    end do
 
  end do

  end subroutine bandec








  subroutine banbks(a,n,m1,m2,al,indx,b)

! This routine solves a band-diagonal system of equations.
! The subroutine must be called immediately after "bandec".
!
! On input, "a" contains the upper triangular matrix,
! and "al" the lower triangular.  Here "n" is the number of
! equations, and (m1,m2) indicate how many subdiagonals
! and superdiagonals are different from zero (e.g. for
! a tridiagonal matrix m1=m2=1, and for a pentadiagonal
! (m1=m2=2).
!
! Also on input "b" contains the right hand side source vector,
! and "indx" records the row permutation effected by the partial
! pivoting as obtained in "bandec".
!
! On output, the vector "b" contains the solution.

  implicit none

  integer n,m1,m2
  integer i,k,l,mm
  
  integer indx(n)

  real(8) b(n),a(n,m1+m2+1),al(n,m1+m2+1)
  real(8) dum  

  mm=m1+m2+1
  l=m1

  do k=1,n
     i=indx(k)
     if (i.ne.k) then
        dum=b(k)
        b(k)=b(i)
        b(i)=dum 
     end if
     if (l.lt.n) l=l+1
     do i=k+1,l
        b(i)=b(i)-al(k,i-k)*b(k)
     end do
  end do

  l=1

  do i=n,1,-1
     dum=b(i)
     do k=2,l
        dum=dum-a(i,k)*b(k+i-1)
     end do
     b(i)=dum/a(i,1)
     if (l.lt.mm) l=l+1
  end do

  end subroutine banbks

