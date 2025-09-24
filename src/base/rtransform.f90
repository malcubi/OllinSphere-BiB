!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/rtransform.f90,v 1.20 2025/09/24 17:17:42 malcubi Exp $

  subroutine rtransform

! *****************************************
! ***   TRANSFORMED RADIAL COORDINATE   ***
! *****************************************

! This subroutine calculates the transformed coordinate
! which is defined via the differential relation:
!
! d(r_trans)/dr = f(r)
!
! If the integral of f(r) is not analitic, we need to do
! a numerical integration. 
!
! Once we transform the radial coordinate, the metric
! coefficients (A,B) must be transformed as:
!
! A_new = A_old (dr_trans/dr)**2
!
! B_new = B_old (r_trans/r)**2
!
! For the moment here we always assume that A_old=B_old=1,
! so this won't work for all types of initial data.
!
! Careful with the notation:  Here "r_trans" is in fact
! the original radial coordinate for which the metric
! coefficients are given by A=B=1, and "r" is actually
! the new radial coordinate that the code will use
! for the evolution.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays
  use derivatives
  use integrals

! Extra variables.

  implicit none

  integer i,j,l,nraux
  real(8) draux,ri,rf,r_new0,raux
  real(8) dRGaussian_dr,dRHat_dr
  real(8) zero,one,two,four
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: u

! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  one  = 1.d0
  two  = 2.d0
  four = 4.d0


! **********************************
! ***   SIGMOID TRANSFORMATION   ***
! **********************************

! This transformation is defined via the differential
! relation:
!
! d(r_trans)/dr  =  dinf / [ 1 + exp( beta r**2 + delta ) ]
!
! We intregrate this equation using either the standard finite
! difference integration routine, or a quadrature based on
! Chebyshev polynomials.

  if (rtype=="sigmoid") then

!    We don't need all this code, there is already a routine that
!    calculates integrals in parallel. We could comment it or even
!    delete it.

     goto 100

!    Initialize.

     r_new0 = 0.d0

     nraux = 1000

!    As we are using a staggered grid, we need to add the 
!    contribution from the interval (0,dr/2) to the integral
!    of the new radial coordinate.

     draux = 0.5d0*dr(0)/dble(nraux)

     do j=1,nraux-1
        raux = dble(j*draux)
        if (mod(j,2)==0) then
           r_new0 = r_new0 +  two*dRGaussian_dr(raux,drinf,beta_transf,delta_transf)
        else
           r_new0 = r_new0 + four*dRGaussian_dr(raux,drinf,beta_transf,delta_transf)
        end if

     end do

     r_new0 = r_new0 + (dRGaussian_dr(0.d0,drinf,beta_transf,delta_transf) &
            + dRGaussian_dr(0.5d0*dr(0),drinf,beta_transf,delta_transf))

     r_new0 = draux*r_new0/3.0d0

!    Integral of f(r) along the grid.

     r_trans(0,1) = r_new0

     draux = dr(0)/dble(nraux)

     do i=2,Nrmax

        do j=1,nraux-1
           raux = r(0,i-1) + j*draux
           if (mod(j,2)==0) then
              r_trans(0,i) = r_trans(0,i) + draux/3.d0*two*dRGaussian_dr(raux,drinf,beta_transf,delta_transf)
           else
              r_trans(0,i) = r_trans(0,i) + draux/3.d0*four*dRGaussian_dr(raux,drinf,beta_transf,delta_transf)
           end if
        end do

        ri = r(0,i-1)
        rf = r(0,i)

        r_trans(0,i) = r_trans(0,i) + draux/3.d0*(dRGaussian_dr(ri,drinf,beta_transf,delta_transf) &
                                    +  dRGaussian_dr(rf,drinf,beta_transf,delta_transf))

        r_trans(0,i) = r_trans(0,i) + r_trans(0,i-1)

     end do

!    Integral using the integration routine.

     100 continue

     auxarray = drinf/(one + dexp(beta_transf*r**2 + delta_transf))

     if (integralmethod == "quadrature") then

!       WARNING: Only works in one level.

        if ((Nl>1)) then
           if (rank==0) then
              print *, 'Multilevel refinement not yet implemented for Chebyshev quadrature.'
              print *, 'Aborting!  (subroutine rtransform)'
           end if
           call die
        end if

        call Chebyshev_quadrature(u)

!       Distribute solution among processors

        if (size==1) then
          r_trans = u
        else
          call distribute(0,Nl-1,r_trans,u)
        end if

     else if (integralmethod == "fourth") then

        intvar => auxarray

        do l=Nl-1,0,-1
           r_trans(l,:) = integral(l)
        end do

!       Restrict integral.

        if (Nl>1) then
           intvar => r_trans
           call restrictintegral
        end if

     end if

!    Ghost zones (only for processor 0).

     if (rank==0) then
        do l=0,Nl-1
           do i=1,ghost
              r_trans(l,1-i) = - r_trans(l,i)
           end do
        end do
     end if

!    If we have refinement levels we need to restrict
!    from fine to coarse grids.

     if (Nl>1) then
        call restrict_rtrans
     end if

!    Metric coefficients (A,B) and lambda.

     A = auxarray**2
     B = (r_trans/r)**2

     lambda = (one - A/B)/r**2


! *********************************************
! ***     SMOOTH HEAVISIDE TRANSFORMATION   ***
! *********************************************

! This transformation is defined via the differential
! relation:
!
! d(r_trans)/dr  =  dinf / [ 1 + Exp(delta)*hat(r) ]
!
! hat(r) = 0.5d0*(tanh(beta*(r+gamma))-dtanh(beta*(r-gamma)))/dtanh(beta*gamma)

! We intregate this equation using a fourth order Simpson
! rule with 1000 points on each interval.

  else if (rtype=="smoothstep") then

!    We don't need all this code, there is already a routine that
!    calculates integrals in parallel. We could comment it or even
!    delete it.

     goto 200

!    Initialize.

     r_new0 = 0.d0

     nraux = 1000

!    As we are using a staggered grid, we need to add the 
!    contribution from the interval (0,dr/2) to the integral
!    of the new radial coordinate.

     draux = 0.5d0*dr(0)/dble(nraux)

     do j=1,nraux-1
        raux = dble(j*draux)
        if (mod(j,2)==0) then
           r_new0 = r_new0 +  two*dRHat_dr(raux,drinf,gamma_transf,beta_transf,delta_transf)
        else
           r_new0 = r_new0 + four*dRHat_dr(raux,drinf,gamma_transf,beta_transf,delta_transf)
        end if

     end do

     r_new0 = r_new0 + (dRHat_dr(0.d0,drinf,gamma_transf,beta_transf,delta_transf) &
            + dRHat_dr(0.5d0*dr(0),drinf,gamma_transf,beta_transf,delta_transf))

     r_new0 = draux*r_new0/3.0d0

!    Integral of f(r) along the grid.

     r_trans(0,1) = r_new0

     draux = dr(0)/dble(nraux)

     do i=2,Nrmax

        do j=1,nraux-1
           raux = r(0,i-1) + j*draux
           if (mod(j,2)==0) then
              r_trans(0,i) = r_trans(0,i) + draux/3.d0*two*dRHat_dr(raux,drinf,gamma_transf,beta_transf,delta_transf)
           else
              r_trans(0,i) = r_trans(0,i) + draux/3.d0*four*dRHat_dr(raux,drinf,gamma_transf,beta_transf,delta_transf)
           end if
        end do

        ri = r(0,i-1)
        rf = r(0,i)

        r_trans(0,i) = r_trans(0,i) + draux/3.d0*(dRHat_dr(ri,drinf,gamma_transf,beta_transf,delta_transf) &
                                    +  dRHat_dr(rf,drinf,gamma_transf,beta_transf,delta_transf))

        r_trans(0,i) = r_trans(0,i) + r_trans(0,i-1)

     end do

!    Integral using the integration routine.

     200 continue

     auxarray = drinf/(one + dexp(delta_transf)*0.5d0 &
              *(dtanh(beta_transf*(r+gamma_transf)) &
              - dtanh(beta_transf*(r-gamma_transf))) &
              / dtanh(beta_transf*gamma_transf))

     intvar => auxarray

     do l=Nl-1,0,-1
        r_trans(l,:) = integral(l)
     end do

!    Restrict integral.

     if (Nl>1) then
        intvar => r_trans
        call restrictintegral
     end if

!    Ghost zones (only for processor 0).

     if (rank==0) then
        do l=0,Nl-1
           do i=1,ghost
              r_trans(l,1-i) = - r_trans(l,i)
           end do
        end do
     end if

!    If we have refinement levels we need to restrict
!    from fine to coarse grids.

     if (Nl>1) then
        call restrict_rtrans
     end if

!    Metric coefficients (A,B) and lambda.

     A = auxarray**2
     B = (r_trans/r)**2

     lambda = (one - A/B)/r**2

  else if (rtype=="sinh") then

     r_trans = beta_transf*dsinh(r/delta_transf)/dsinh(one/delta_transf)
     A = (beta_transf*dcosh(r/delta_transf)/dsinh(one/delta_transf)/delta_transf)**2
     B = (r_trans/r)**2

     lambda = (one - A/B)/r**2


! ****************************************
! ***   CHOPTUIK TYPE TRANSFORMATION   ***
! ****************************************

! This transformation is analytical and has the form:
!
! r_trans  =  exp(r + delta)  -  exp(delta)
!
!          +  beta/(beta - delta - r)  -  beta/(beta - delta)

  else if (rtype=="choptuik") then

!    Radial transformation.

     r_trans = exp(abs(r) + delta_transf) - exp(delta_transf) &
             + beta_transf/(beta_transf - delta_transf - abs(r)) &
             - beta_transf/(beta_transf - delta_transf)

!    Ghost zones.

     do i=1,ghost
        r_trans(:,1-i) = - r_trans(:,i)
     end do

!    Metric coefficients (A,B) and lambda.

     A = (exp(abs(r) + delta_transf) + beta_transf/(beta_transf - delta_transf - abs(r))**2)**2
     B = (r_trans/r)**2

     lambda  = (one - A/B)/r**2

  end if


! **********************************************
! ***   DERIVATIVES OF METRIC COEFFICIENTS   ***
! **********************************************

  if (rtype=="sigmoid") then

!    Derivatives of metric coefficients (A,B).

     do l=0,Nl-1
          D1_A(l,:) = -4.d0*beta_transf*r(l,:)*A(l,:)*dexp(beta_transf*r(l,:)**2+delta_transf)/ &
                      (one+dexp(beta_transf*r(l,:)**2+delta_transf))
          D2_A(l,:) = 4.d0*beta_transf*drinf**2*dexp(beta_transf*r(l,:)**2+ delta_transf)*&
                      (-one-two*beta_transf*r(l,:)**2+dexp(beta_transf*r(l,:)**2+delta_transf)*&
                      (-one+4.d0*beta_transf*r(l,:)**2))/(one+dexp(beta_transf*r(l,:)**2+delta_transf))**4
     end do

     do l=0,Nl-1
          D1_B(l,:) = two*r_trans(l,:)/r(l,:)**3*(r(l,:)*drinf/&
                      (one+dexp(beta_transf*r(l,:)**2 + delta_transf))-r_trans(l,:))
          D2_B(l,:) = two/r(l,:)**4*(3.d0*r_trans(l,:)**2+&
                      drinf*r(l,:)**2*(drinf+two*beta_transf*r(l,:)*r_trans(l,:))/&
                      (one+dexp(beta_transf*r(l,:)**2+delta_transf))**2-&
                      two*drinf*r(l,:)*(two+beta_transf*r(l,:)**2)*r_trans(l,:)/&
                      (one+dexp(beta_transf*r(l,:)**2+delta_transf)))
     end do

  else

    diffvar => A

    do l=0,Nl-1
       D1_A(l,:) = diff1(l,+1)
       D2_A(l,:) = diff2(l,+1)
    end do

    diffvar => B

    do l=0,Nl-1
       D1_B(l,:) = diff1(l,+1)
       D2_B(l,:) = diff2(l,+1)
    end do

  end if


! *****************************
! ***   MESSAGE TO SCREEN   ***
! *****************************

  rtmin = r_trans(Nl-1,1)
  rtmax = r_trans(0,Nrl(rank))

  if (size>1) then
     call MPI_BCAST(rtmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(rtmax,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)
  end if

  if (rank==0) then
     print *, "We are using a transformed radial coordinate. Boundary points are:"
     print *, "rmin =",rtmin
     print *, "rmax =",rtmax
     print *
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine rtransform








  function dRGaussian_dr(r,drinf,k,a) 

! *******************************************************
! ***   DERIVATIVE dR/dr FOR SIGMOID TRANSFORMATION   ***
! *******************************************************

  implicit none

  real(8) dRGaussian_dr
  real(8) r,drinf,k,a

  dRGaussian_dr = drinf/(1.d0 + dexp(k*r**2 + a))

  end function dRGaussian_dr









  function dRHat_dr(r,drinf,r0,s,a) 

! ****************************************************************
! ***   DERIVATIVE dR/dr FOR SMOOTH HEAVISIDE TRANSFORMATION   ***
! ****************************************************************

  implicit none

  real(8) dRHat_dr
  real(8) r,drinf,hat,r0,s,a

  hat = 0.5d0*(dtanh(s*(r+r0))-dtanh(s*(r-r0)))/dtanh(s*r0)

  dRHat_dr = drinf/(1.d0 + dexp(a)*hat)


! ***************
! ***   END   ***
! ***************

  end function dRHat_dr






  subroutine restrict_rtrans

! ***************************
! ***   RESTRICT RTRANS   ***
! ***************************

! For the case of several grid levels we need to restrict r_trans
! from the fine to the coarse grids.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer i,i0,l

  real(8) delta,m0,m1
  real(8) r0,interp


! ***********************
! ***   RESTRICTION   ***
! ***********************

  restrictvar => r_trans

  do l=Nl-1,1,-1

!    Substract constant difference at edge of fine grid. Since the
!    integrals are done indendently at each grid level, we find that
!    after restriction there will be small jumps due to accumulated
!    numerical error when we pass from the end on a fine grid to the
!    next coarse grid. Here we correct for that.

     i0 = Nr/2

!    Find delta.

     if (size==1) then

!       Find difference between the value at the edge of the fine grid, 
!       and the value in the coarse grid.

        delta = r_trans(l-1,i0) - 0.5d0*(r_trans(l,2*i0) + r_trans(l,2*i0-1))

     else

!       Find values of radius (r0) and r_trans (m0) at edge of fine grid,
!       and send them to all processors.  Notice that the edge of the grid
!       belongs to the processor with rank=size-1.

        r0 = 0.5d0*(r(l,2*i0)+r(l,2*i0-1))
        m0 = 0.5d0*(r_trans(l,2*i0)+r_trans(l,2*i0-1))

        call MPI_BCAST(r0,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(m0,1,MPI_DOUBLE_PRECISION,size-1,MPI_COMM_WORLD,ierr)

!       Now interpolate value of r_trans at r0 from coarse grid,
!       and send it to all processors.  Notice that the point r0
!       is in the middle of the coarse grid, and belongs to the
!       processor with rank=(size-1)/2.

        interpvar => r_trans
        m1 = interp(l-1,r0,.false.)
        call MPI_BCAST(m1,1,MPI_DOUBLE_PRECISION,(size-1)/2,MPI_COMM_WORLD,ierr)

!       Find difference bewteen m1 and m0.

        delta = m1 - m0

     end if

!    Subtract delta.

     r_trans(l-1,:) = r_trans(l-1,:) - delta

!    Restrict.

     call restrict(l,.false.)

!    Fix symmetries.

     if (rank==0) then
        do i=1,ghost
           r_trans(l-1,1-i) = - r_trans(l-1,i)
        end do
     end if

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine restrict_rtrans









  subroutine Chebyshev_quadrature(u)

! ********************************
! ***   CHEBYSHEV QUADRATURE   ***
! ********************************

! The problem of calculating the integral of f can be restated 
! as the solution of the ODE u'(x)=f(x).
!
! In order to use Chebyshev quadrature, we change the interval of
! integration [a,b] to [-1,1], and set the boundary condition u(-1)=0.
!
! This ODE can be representated by the linear problem:
!
! D1 u = f
!
! where D1 is the actual differential operator, u is the solution,
!  and f the original integrand.
!
! With this restatement the integral of f will be the first value
!  of the column vector u.
!
! u = D1^-1 f
!
! where D1^-1 is the inverse operator of D1.

! Include modules.

  use mpi
  use procinfo
  use param
  use arrays

! Extra variables.

  implicit none

  integer i,j,k,l,nraux

  real(8) zero,one,two,four,smallpi
  real(8) c
  real(8) r1,r2                                         ! Integration limits.

  real(8), dimension(0:NCheb):: xc,caux                 ! Chebyshev points cos(i*pi/N), i=0,1,...,N.
  real(8), dimension(0:NCheb-1):: xcaux                 
  real(8), dimension(0:NCheb-1):: f                     ! Integrand evaluated at collocation points.
  real(8), dimension(0:NCheb):: iD1f                    ! Integrand multiplied by the inverse differential operator.
  real(8), dimension(0:NCheb,0:NCheb) :: D1             ! Derivative matrix.
  real(8), dimension(0:NCheb-1,0:NCheb-1) :: D,D_inv    ! Final differential operator and inverse.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: u      ! Global solution of ODE.
!  LAPLACK AUXILIARIES
!  integer info                                 !LAPACK error flag.
!  integer, dimension(0:NCheb-1) :: ipiv        !LAPACK pivoting array.
!  real(8), dimension(0:NCheb-1) :: work        !work array for LAPACK.
!  external DGETRF
!  external DGETRI


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.d0
  one  = 1.d0
  two  = 2.d0
  four = 4.d0

  smallpi = acos(-1.d0)


! *****************************
! ***   INITIALIZE ARRAYS   ***
! *****************************

  xc   = zero
  caux = zero
  D1   = zero
  f    = zero 
  iD1f = zero

  D = zero
  D_inv = zero

! ipiv = zero
! work = zero
  

! *******************************
! ***   PSEUDOSPECTRAL GRID   ***
! *******************************

! Fill the Chebyshev points and his weights.

  do i=0,NCheb
     xc(i)=cos(dble(i*smallpi/NCheb))
  end do

  do i=0,NCheb
     if ((i==0).or.(i==NCheb)) then
        caux(i) = two
     else
        caux(i) = one
     end if
  end do


! ******************************
! *** DIFFERENTIAL OPERATORS ***
! ******************************

! Fill the derivative operators in D1.

! Upper left and lower right corners.

  D1(0,0)= dble((two*NCheb**2+one)/6.d0)
  D1(NCheb,NCheb) = -dble((two*NCheb**2+one)/6.d0)

! Diagonal.

  do i=1,NCheb-1
     D1(i,i) = -0.5d0*xc(i)/(one-xc(i)**2)
  end do

! Off-diagonal.

  do i=0,NCheb
     do j=0,NCheb
        if (i/=j) D1(i,j) = caux(i)/caux(j) * (-one)**(i+j)/(xc(i)-xc(j))
     end do
  end do

! Before calculating D1 inverse, we need to aply the boundary condition
! to this operator, otherwise D1 is singular. The easiest way to do this
! is to strip the last row and column of D1.

! Strip last row and column of D1.

  do i=0,NCheb-1
    do j=0,NCheb-1
      D(i,j) = D1(i,j)
    end do
  end do

  do i=0,NCheb-1
    xcaux(i) = xc(i)
  end do


! ********************************************
! ***   INVERT THE DIFFERENTIAL OPERATOR   ***
! ********************************************

! This part depends on LAPACK. I comment this and use a simple solver 
! to find the inverse operator.
!
! D_inv = D
! nraux = Ncheb
!
! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.

! call DGETRF(nraux, nraux, D_inv, nraux, ipiv, info)

! if (info /= 0) then
!    stop 'Operator is numerically singular!'
!    call die
! end if

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.

! call DGETRI(nraux, D_inv, nraux, ipiv, work, nraux, info)

! if (info /= 0) then
!    stop 'Matrix inversion failed!'
!    call die
! end if  

  call invrsmtx(D,D_inv,NCheb)


! *****************************************
! ***   CHEBYSHEV QUADRATURE INTEGRAL   ***
! *****************************************

! Evaluate the integrand at the collocation points,
! but we must be careful because Chebyshev quadrature only
! works in the interval [-1,1], so we need to change
! the integrand.

  !do i=1,Nrmax
  if (rank == 0) then
    do i=1,Nrtotal
      if (i==1) then
         r1 = 0.d0
         r2 = (dble(Nmin(rank) + i) - 0.5d0)*dr(0)
      else
         r1 = (dble(Nmin(rank) + i-1) - 0.5d0)*dr(0)
         r2 = (dble(Nmin(rank) + i) - 0.5d0)*dr(0)
      end if

      f = 0.5d0*(r2-r1)*drinf/(one + dexp(beta_transf*(0.5d0*(r2-r1)*xcaux &
        + 0.5d0*(r1+r2))**2 + delta_transf))

      iD1f = matmul(D_inv,f)
      !r_trans(0,i) = r_trans(0,i-1)+iD1f(0)
      u(0,i) = u(0,i-1)+iD1f(0)
    end do
  end if


! ***************
! ***   END   ***
! ***************

  end subroutine Chebyshev_quadrature









  subroutine invrsmtx(A,C,n)

! *************************
! ***   INVERT MATRIX   ***
! *************************

! Calculate the inverse of A using a LU decomposition.
! The inverse matrix is stored in C.
!
! WARNING: Original A will be destroyed.

  implicit none

  real(8) A(n,n), C(n,n)
  real(8) L(n,n), U(n,n), b(n), d(n), x(n)
  real(8) coeff
  integer i, j, k , n

! Step 0: initialization for matrices L and U and b.

  L = 0.0d0
  U = 0.0d0
  b = 0.d0

! Step 1: forward elimination.

  do k=1, n-1
    do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
        a(i,j) = a(i,j)-coeff*a(k,j)
      end do
    end do
  end do

! Step 2: prepare L and U matrices.
! L is the matrix of the elimination coefficient
! + the diagonal elements are 1.0.

  do i=1,n
    L(i,i) = 1.0
  end do

! U matrix is the upper triangular part of A.

  do j=1,n
    do i=1,j
      U(i,j) = a(i,j)
    end do
  end do

! Step 3: compute columns of the inverse matrix C.

  do k=1,n

     b(k)=1.0
     d(1)=b(1)

!    Step 3a: Solve Ld=b using the forward substitution.

     do i=2,n
        d(i)=b(i)
        do j=1,i-1
           d(i) = d(i) - L(i,j)*d(j)
        end do
     end do

!    Step 3b: Solve Ux=d using the back substitution.

     x(n)=d(n)/U(n,n)

     do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
           x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
     end do

!    Step 3c: fill the solutions x(n) into column k of C.

     do i=1,n
        c(i,k) = x(i)
     end do
     b(k)=0.0

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine invrsmtx
