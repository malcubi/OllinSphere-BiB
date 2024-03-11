!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/difeq.f90,v 1.3 2024/03/01 19:38:19 malcubi Exp $

  subroutine difeq(kk,kk1,kk2,jsf,is1,isf,indexv,ss,rr,yy,ll)

! *******************************************************
! ***   MATRIX OF DERIVATIVES FOR RELAXATION METHODS   **
! *******************************************************
!
! This routine generalizes the original Numerical Recipes
! routine for calculating the matrix of derivatives for
! relaxation methods. It is called from "solvde".
! Notice that it must be adapted to each problem!
!
! On entry we have:

! kk   =  current mesh point
! kk1  =  first point
! kk2  =  last point
! ll   =  current grid level
!
! The array yy contains the current guess for the solution.
! The routine sets the value of the matrix ss, which corresponds
! to the derivatives of the system of equation, and the
! equations themselves.
!
! The matrix SS has dimensions (ne,2*ne+1), where ne is
! the total number of equations to solve.
!
! The equations are assumed to be already in finite difference
! form using only 2 points (second order).  If the original
! equations have the form:
!
! du(j)/dr  -  g[r,u(j)] = 0
!
! with j indicating a given equation, then they must be
! represented as:
!
! E(j) := u(j,k) - u(j,k-1)  -  g[r(k+1/2),u(j,k),u(j,k-1)]  = 0
!
! where k numbers the grid points.  Given the equations,
! we define the matrix S of derivatives at each point
! in the following way:
!
! S(j,i)     :=  dE(j)/du(i,k-1)
! S(j,ne+i)  :=  dE(j)/du(i,k)
!
! with i=1,...,ne.  Finally, the component S(j,2*ne+1)
! corresponds to the equation itself:
!
! S(j,2*ne+1)  :=  E(j)
!
! Some care must be taken at the boundaries.  At the first
! point there are nb<ne boundary conditions that must
! have some dependence on the FIRST nb functions evaluated
! only at the first grid point.  The boundary conditions and
! derivatives must be stored in the LAST rows of S(j,i)
! (that is S(j,1) is non zero only for j>(ne-nb)),
! and must involve only the values of the variables
! corresponding to point k (and not k-1), so that i>ne.
!
! At the last grid point the remaining (ne-nb) boundary
! conditions must correspond to the FIRST columns of S.
! But in this case they can involve both points k and k-1.
!
! For example, for 1 boundary condition (nb=1) the only
! non-zero components of S for the first grid point should
! be S(ne,i) with i>ne, and for the last grid the
! only non-zero components should have j<ne.

! Include modules.

  use param
  use arrays

! Extra variables.

  implicit none

  integer, intent(in) :: is1,isf,jsf,kk,kk1,kk2,ll
  integer, dimension(:), intent(in) :: indexv

  real(8) r0,V0
  real(8) smallpi,aux0,aux1,aux2
  real(8) s1,s2,s3,s4,s5,s6

  real(8), dimension(:), intent(in)  :: rr
  real(8), dimension(:,:), intent(out) :: ss
  real(8), dimension(:,:), intent(in)  :: yy


! *******************
! ***   NUMBERS   ***
! *******************

  smallpi = acos(-1.d0)


! ************************
! ***   SCALAR PULSE   ***
! ************************

! In this case we have 2 equations, one for psi and one for dpsi/dr.
! We have chosen y1=dpsi/dr, y2=psi.
!
! The matrix ss is a 2x5 matrix.  The equations are in column 5
! for each variable, while the derivatives with respect to y(j,kk-1)
! are on columns 1 and 2, and with respect to y(j,kk) on columns 4 and 5.

  if (idata=="scalarpulse") then

!    Initialize matrix S.

     ss = 0.d0

!    Boundary condition at first point.  Here we use
!    the fact that dpsi/dr=0 at the origin. If we now
!    expand the solution in Taylor close to the origin
!    to second order we find that we must have:
!
!    y1  + (pi/3) r [ scalar_xi**2 y2 + 2 scalar_V y2^5 ]  ~  0

     if (kk==kk1) then

!       Equation at left boundary.

        aux0 = smallpi/3.d0
        aux1 = aux0*rr(1)*scalar_xi(ll,1)**2
        aux2 = aux0*rr(1)*scalar_V(ll,1)

        ss(2,5) = yy(1,1) + aux1*yy(2,1) + 2.d0*aux2*yy(2,1)**5

!       Derivatives.

        ss(2,3) = 1.d0
        ss(2,4) = aux1 + 10.d0*aux2*yy(2,1)**4
 
!    Boundary condition at last point.  Here we impose
!    the condition that far away we have psi ~ 1 + c/r
!    with c a constant.  This condition reduces to:
!
!    dpsi/dr + (psi-1)/r = 0
!
!    which upon finite differencing becomes:
!
!    y1  +  (y2-1) / r  = 0
!
!    We apply this directly on the last point kk2=kk-1.
!
!    Also, in this case the matrix S only has non-zero
!    components for j=1, as the boundary condition does
!    not involve the point kk (which would correspond to
!    taking j=2).
!
!    The above condition is only imposed for the coarse
!    grid ll=0.  For finer grid we impose a Dirichlet
!    boundary condition by interpolating psi from the
!    coarser grid level (linear interpolation).

     else if (kk>kk2) then

        if (ll==0) then

!          Equation at right boundary.

           aux1 = 1.d0/rr(kk2)
           ss(1,5) = yy(1,kk2) + aux1*(yy(2,kk2) - 1.d0)

!          Derivatives.

           ss(1,3) = 1.d0
           ss(1,4) = aux1

        else

!          Equation at right boundary.

           r0 = (dble(Nrtotal)-0.5d0)*dr(ll)
           aux0 = (r0 - (dble(Nrtotal/2)-0.5d0)*dr(ll-1))/dr(ll-1)
           aux1 = psi(ll-1,Nrtotal/2) + aux0*(psi(ll-1,Nrtotal/2+1)-psi(ll-1,Nrtotal/2))

           ss(1,5) = yy(2,kk2) - aux1

!          Derivatives.

           ss(1,4) = 1.d0

        end if

!    Interior points.

     else

!       Equation 1.  This is just the definition of y1:
!
!       y1  :=  dy2/dr
!
!       so that we find:
!
!       0  =  y2(k) - y2(k-1) - dr/2 (y1(k) + y1(k-1))

        aux1 = 0.5d0*dr(ll)

        ss(1,5) = yy(2,kk) - yy(2,kk-1) - aux1*(yy(1,kk) + yy(1,kk-1))

!       Derivatives.

        ss(1,1) = - aux1
        ss(1,2) = - 1.d0
        ss(1,3) = - aux1
        ss(1,4) = + 1.d0

!       Equation 2.  This is the Hamiltonian constraint, which
!       we write as:
!
!       0  =  y1(k) - y1(k-1)  +  2 dr (y1(k) + y1(k-1)) / (r(k) + r(k-1))
!
!          +  pi dr [ (scalar_xi(k)^2 + scalar_xi(k-1)**2) (y2(k) + y2(k-1)) / 4
!
!          +  (scalar_V(k) + scalar_V(k-1)) (y2(k) + y2(k-1))^5 / 32 ]

        aux0 = smallpi*dr(ll)
        aux1 = aux0*(scalar_xi(ll,kk)**2 + scalar_xi(ll,kk-1)**2)/4.d0
        aux2 = aux0*(scalar_V(ll,kk) + scalar_V(ll,kk-1))

        aux0 = 2.d0*dr(ll)/(rr(kk)+rr(kk-1))
        ss(2,5) = yy(1,kk) - yy(1,kk-1) + aux0*(yy(1,kk) + yy(1,kk-1)) &
                + aux1*(yy(2,kk) + yy(2,kk-1)) + aux2*(yy(2,kk) + yy(2,kk-1))**5/32.d0

!       Derivatives.

        ss(2,1) = - 1.d0 + aux0
        ss(2,2) = aux1 + 5.d0*aux2*(yy(2,kk) + yy(2,kk-1))**4/32.d0
        ss(2,3) = + 1.d0 + aux0
        ss(2,4) = aux1 + 5.d0*aux2*(yy(2,kk) + yy(2,kk-1))**4/32.d0

     end if


! ********************************************
! ***   BOSON STARS IN POLAR AREAL GAUGE   ***
! ********************************************

! In this case we have 5 equations for 5 variables:
!
! y1  =  A
! y2  =  alpha
! y3  =  phi
! y4  =  dphi/dr  =  xi
! y5  =  omega  (eigenvalue)
!
! The matrix SS is a 5x11 matrix.  The equations are in column 11
! for each variable with all terms on the same size (that is EQ=0),
! while the derivatives with respect to y(j,kk-1) are on columns 1-5,
! and with respect to y(j,kk) on columns 6-10. The expressions below
! are obtained with MAPLE, so they look nasty.
!
! Also, notice that at the moment I always take the potential to
! be just a mass term (might change this later).

  else if ((idata=="bosonstar").and.(boson_gauge=="PA")) then

!    Initialize matrix S.

     ss = 0.d0

!    Boundary conditions at first point (4 conditions).
!    These boundary conditions correspond to:
!
!    A(r<<1)     -  (1  +  A2 r^2)  =  0 
!    alpha(r<<1) -  1    =  0
!    phi(r=0)    -  phi0 =  0
!    xi(r=0)             =  0
!
!    Notice that the value of (alpha,phi) are arbitrary
!    at the first point.  For A we do a Taylor expansion
!    to second order from the Hamiltonian constraint to
!    find the value of A2:
!
!    A2 = (4 pi / 3) ( m^2 + omega^2 ) phi0^2
!
!    Remember that the 4 boundary conditions at the left
!    must involve the LAST 4 rows of SS.

     if (kk==kk1) then

!       Equations at left boundary.

        r0 = rr(kk1)

        ss(2,11) = yy(1,1) - 1.d0 - 4.d0*smallpi*r0**2/3.d0*(complex_mass**2 + yy(5,1)**2)*boson_phi0**2
        ss(3,11) = yy(2,1) - 1.d0
        ss(4,11) = yy(3,1) - boson_phi0
        ss(5,11) = yy(4,1) - r0/3.d0*(complex_mass**2 - yy(5,1)**2)*boson_phi0

!       Derivatives. Remember that the values of
!       the variables at the first point correspond
!       to the columns (ne+1,...) since there is
!       no dependence on point kk-1.

        ss(2,6 ) = 1.d0
        ss(2,10) = - 8.d0*r0**2*smallpi/3.d0*yy(5,1)*boson_phi0**2

        ss(3,7 ) = 1.d0
        ss(4,8 ) = 1.d0

        ss(5,9 ) = 1.d0
        ss(5,10) = 2.d0*r0/3.d0*yy(5,1)*boson_phi0

!    Boundary condition at last point (1 condition).
!    At the moment I assume that phi decays as 1/r,
!    so that we have:
!
!    yy4 + yy3/r = 0
!
!    Of course, phi should really decay exponentially, but
!    far away it makes very little difference since phi
!    is very small anyway, and this choice guarantees that
!    if the boundary is not too far phi still decays.
!
!    Since we only have 1 boundary condition here, it must
!    correspond to the FIRST row of SS.

     else if (kk>kk2) then

!       Equation at right boundary.

        r0 = rr(kk2)
        ss(1,11) = yy(3,kk2)/r0 + yy(4,kk2)

!       Derivatives.

        ss(1,8) = 1.d0/r0
        ss(1,9) = 1.d0

!    Interior points.

     else

!       Calculate average radius.

        r0 = 0.5d0*(rr(kk) + rr(kk-1))

!       Equation 1.  This is the Hamiltonian constraint:
!
!       dA/dr  =  A [ (1 - A) / r + 4 pi r A ( 2 V + (omega/alpha)**2 phi**2 + xi**2 / A ) ]

        ss(1,11) = yy(1,kk)-yy(1,kk-1)-4*dr(ll)*(yy(1,kk)/2.D0+yy(1,kk-1)/2.D0) &
                 *(smallpi*((yy(3,kk)/2.D0+yy(3,kk-1)/2.D0)**2*(complex_mass**2 &
                 +(yy(5,kk)+yy(5,kk-1))**2/(yy(2,kk)+yy(2,kk-1))**2) &
                 *(yy(1,kk)/2.D0+yy(1,kk-1)/2.D0)+(yy(4,kk)/2.D0+yy(4,kk-1)/2.D0)**2)*r0**2 &
                 -yy(1,kk)/8.D0-yy(1,kk-1)/8.D0+1.D0/4.D0)/r0

!       Derivatives.

        ss(1,1) = -1.D0-2.D0*dr(ll)*(smallpi*((yy(3,kk)/2+yy(3,kk-1)/2)**2 &
                *(complex_mass**2+(yy(5,kk)+yy(5,kk-1))**2/(yy(2,kk)+yy(2,kk-1))**2) &
                *(yy(1,kk)/2+yy(1,kk-1)/2)+(yy(4,kk)/2+yy(4,kk-1)/2)**2)*r0**2 &
                -yy(1,kk)/8-yy(1,kk-1)/8+1.D0/4.D0)/r0-4*dr(ll)*(yy(1,kk)/2+yy(1,kk-1)/2) &
                *(smallpi*(yy(3,kk)/2+yy(3,kk-1)/2)**2*(complex_mass**2 &
                +(yy(5,kk)+yy(5,kk-1))**2/(yy(2,kk)+yy(2,kk-1))**2)*r0**2/2-1.D0/8.D0)/r0

        ss(1,2) = 8*dr(ll)*(yy(1,kk)/2+yy(1,kk-1)/2)**2*smallpi*(yy(3,kk)/2+yy(3,kk-1)/2)**2 &
                *(yy(5,kk)+yy(5,kk-1))**2/(yy(2,kk)+yy(2,kk-1))**3*r0

        ss(1,3) = -4*dr(ll)*(yy(1,kk)/2+yy(1,kk-1)/2)**2*smallpi*(yy(3,kk)/2+yy(3,kk-1)/2) &
                *(complex_mass**2+(yy(5,kk)+yy(5,kk-1))**2/(yy(2,kk)+yy(2,kk-1))**2)*r0

        ss(1,4) = -4*dr(ll)*(yy(1,kk)/2+yy(1,kk-1)/2)*smallpi*(yy(4,kk)/2+yy(4,kk-1)/2)*r0

        ss(1,5) = -8*dr(ll)*(yy(1,kk)/2+yy(1,kk-1)/2)**2*smallpi*(yy(3,kk)/2+yy(3,kk-1)/2)**2 &
                *(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))**2*r0

        ss(1,6) = 1-2*dr(ll)*(smallpi*((yy(3,kk)/2+yy(3,kk-1)/2)**2*(complex_mass**2 &
                +(yy(5,kk)+yy(5,kk-1))**2/(yy(2,kk)+yy(2,kk-1))**2)*(yy(1,kk)/2+yy(1,kk-1)/2) &
                +(yy(4,kk)/2+yy(4,kk-1)/2)**2)*r0**2-yy(1,kk)/8-yy(1,kk-1)/8+1.D0/4.D0)/r0 &
                -4*dr(ll)*(yy(1,kk)/2+yy(1,kk-1)/2)*(smallpi*(yy(3,kk)/2+yy(3,kk-1)/2)**2 &
                *(complex_mass**2+(yy(5,kk)+yy(5,kk-1))**2/(yy(2,kk)+yy(2,kk-1))**2)*r0**2/2-1.D0/8.D0)/r0

        ss(1,7) = 8*dr(ll)*(yy(1,kk)/2+yy(1,kk-1)/2)**2*smallpi*(yy(3,kk)/2+yy(3,kk-1)/2)**2 &
                *(yy(5,kk)+yy(5,kk-1))**2/(yy(2,kk)+yy(2,kk-1))**3*r0

        ss(1,8) = -4*dr(ll)*(yy(1,kk)/2+yy(1,kk-1)/2)**2*smallpi*(yy(3,kk)/2+yy(3,kk-1)/2) &
                *(complex_mass**2+(yy(5,kk)+yy(5,kk-1))**2/(yy(2,kk)+yy(2,kk-1))**2)*r0

        ss(1,9) = -4*dr(ll)*(yy(1,kk)/2+yy(1,kk-1)/2)*smallpi*(yy(4,kk)/2+yy(4,kk-1)/2)*r0

        ss(1,10) = -8*dr(ll)*(yy(1,kk)/2+yy(1,kk-1)/2)**2*smallpi*(yy(3,kk)/2+yy(3,kk-1)/2)**2 &
                 *(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))**2*r0

!       Equation 2. This is the lapse condition:
!
!       dalpha/dr  =  alpha [ (A - 1) / r + (dA/dr) / (2A) - 8 pi r A V ]
!
!       which upon substituting dA/dr becomes:
!
!       dalpha/dr  =  alpha [ (A - 1) / 2r  -  2 pi r A ( 2 V - (omega/alpha)**2 phi**2 - xi**2 / A ) ]

        ss(2,11) = yy(2,kk)-yy(2,kk-1)+2.D0*(yy(2,kk)/2.D0+yy(2,kk-1)/2.D0)*((-1.D0/4.D0 &
                 +(yy(3,kk)/2.D0+yy(3,kk-1)/2.D0)**2*smallpi*(complex_mass-(yy(5,kk)+yy(5,kk-1)) &
                 /(yy(2,kk)+yy(2,kk-1)))*(complex_mass+(yy(5,kk)+yy(5,kk-1)) &
                 /(yy(2,kk)+yy(2,kk-1)))*r0**2)*(yy(1,kk)/2.D0+yy(1,kk-1)/2.D0) &
                 -r0**2*smallpi*(yy(4,kk)/2.D0+yy(4,kk-1)/2.D0)**2+1.D0/4.D0)*dr(ll)/r0

!       Derivatives.

        ss(2,1) = 2*(yy(2,kk)/2+yy(2,kk-1)/2)*(-1.D0/8.D0+(yy(3,kk)/2+yy(3,kk-1)/2)**2 &
                *smallpi*(complex_mass-(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))) &
                *(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0**2/2)*dr(ll)/r0

        ss(2,2) = -1.d0 + ((-1.D0/4.D0+(yy(3,kk)/2+yy(3,kk-1)/2)**2*smallpi*(complex_mass &
                -(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*(complex_mass+(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1)))*r0**2)*(yy(1,kk)/2+yy(1,kk-1)/2) &
                -r0**2*smallpi*(yy(4,kk)/2+yy(4,kk-1)/2)**2+1.D0/4.D0)*dr(ll)/r0 &
                + 2*(yy(2,kk)/2+yy(2,kk-1)/2)*((yy(3,kk)/2+yy(3,kk-1)/2)**2*smallpi*(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1))**2*(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0**2 &
                -(yy(3,kk)/2+yy(3,kk-1)/2)**2*smallpi*(complex_mass-(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1)))*(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))**2*r0**2) &
                *(yy(1,kk)/2+yy(1,kk-1)/2)*dr(ll)/r0

        ss(2,3) = 2*(yy(2,kk)/2+yy(2,kk-1)/2)*(yy(3,kk)/2+yy(3,kk-1)/2)*smallpi &
                *(complex_mass-(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))) &
                *(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0*(yy(1,kk)/2+yy(1,kk-1)/2)*dr(ll)

        ss(2,4) = -2*(yy(2,kk)/2+yy(2,kk-1)/2)*smallpi*(yy(4,kk)/2+yy(4,kk-1)/2)*r0*dr(ll)

        ss(2,5) = 2*(yy(2,kk)/2+yy(2,kk-1)/2)*(-(yy(3,kk)/2+yy(3,kk-1)/2)**2*smallpi/(yy(2,kk)+yy(2,kk-1)) &
                *(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0**2 &
                +(yy(3,kk)/2+yy(3,kk-1)/2)**2*smallpi*(complex_mass-(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1)))/(yy(2,kk)+yy(2,kk-1))*r0**2)*(yy(1,kk)/2+yy(1,kk-1)/2)*dr(ll)/r0

        ss(2,6) = 2*(yy(2,kk)/2+yy(2,kk-1)/2)*(-1.D0/8.D0+(yy(3,kk)/2+yy(3,kk-1)/2)**2 &
                *smallpi*(complex_mass-(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))) &
                *(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0**2/2)*dr(ll)/r0

        ss(2,7) = 1.d0+((-1.D0/4.D0+(yy(3,kk)/2+yy(3,kk-1)/2)**2*smallpi*(complex_mass-(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1)))*(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0**2) &
                *(yy(1,kk)/2+yy(1,kk-1)/2)-r0**2*smallpi*(yy(4,kk)/2+yy(4,kk-1)/2)**2+1.D0/4.D0)*dr(ll)/r0 &
                +2*(yy(2,kk)/2+yy(2,kk-1)/2)*((yy(3,kk)/2+yy(3,kk-1)/2)**2*smallpi*(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1))**2*(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0**2 &
                -(yy(3,kk)/2+yy(3,kk-1)/2)**2*smallpi*(complex_mass-(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1)))*(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))**2*r0**2)*(yy(1,kk)/2+yy(1,kk-1)/2)*dr(ll)/r0

        ss(2,8) = 2*(yy(2,kk)/2+yy(2,kk-1)/2)*(yy(3,kk)/2+yy(3,kk-1)/2)*smallpi*(complex_mass-(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1)))*(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0 &
                *(yy(1,kk)/2+yy(1,kk-1)/2)*dr(ll)

        ss(2,9) = -2*(yy(2,kk)/2+yy(2,kk-1)/2)*smallpi*(yy(4,kk)/2+yy(4,kk-1)/2)*r0*dr(ll)

        ss(2,10) = 2*(yy(2,kk)/2+yy(2,kk-1)/2)*(-(yy(3,kk)/2+yy(3,kk-1)/2)**2*smallpi &
                 /(yy(2,kk)+yy(2,kk-1))*(complex_mass+(yy(5,kk)+yy(5,kk-1)) &
                 /(yy(2,kk)+yy(2,kk-1)))*r0**2+(yy(3,kk)/2+yy(3,kk-1)/2)**2*smallpi &
                 *(complex_mass-(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))) &
                 /(yy(2,kk)+yy(2,kk-1))*r0**2)*(yy(1,kk)/2+yy(1,kk-1)/2)*dr(ll)/r0

!       Equation 3. This is just the definition of xi:
!
!       y4  :=  dy3/dr

        aux1 = 0.5d0*dr(ll)

        ss(3,11) = yy(3,kk) - yy(3,kk-1) - aux1*(yy(4,kk) + yy(4,kk-1))

!       Derivatives.

        ss(3,3) = - 1.d0
        ss(3,4) = - aux1

        ss(3,8) = + 1.d0
        ss(3,9) = - aux1

!       Equation 4. This is the Klein-Gordon equation:
!
!       dxi/dr  =  - xi [ 2/r + (dalpha/dr)/alpha - (dA/dr)/(2A) ] + A [ V' - (omega/alpha)**2 phi ]
!
!       which upon substituting dA/dr and daplha/dr becomes:
!
!       dxi/dr  =  - xi [ (A + 1)/r - 8 pi r A V ]  +  A  ( V' - (omega/alpha)**2 phi )

        ss(4,11) = ((-4*(yy(1,kk)/2+yy(1,kk-1)/2)*complex_mass**2 &
                 *(yy(3,kk)/2+yy(3,kk-1)/2)**2*r0**2*smallpi*(yy(4,kk)/2+yy(4,kk-1)/2) &
                 -(yy(1,kk)/2+yy(1,kk-1)/2)*(yy(3,kk)/2+yy(3,kk-1)/2) &
                 *(complex_mass-(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))) &
                 *(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0 &
                 +(yy(4,kk)/2+yy(4,kk-1)/2)*(yy(1,kk)/2+yy(1,kk-1)/2+1))*dr(ll)+(yy(4,kk)-yy(4,kk-1))*r0)/r0

!       Derivatives.

        ss(4,1) = (-2*complex_mass**2*(yy(3,kk)/2+yy(3,kk-1)/2)**2*r0**2*smallpi &
                *(yy(4,kk)/2+yy(4,kk-1)/2)-(yy(3,kk)/2+yy(3,kk-1)/2)*(complex_mass &
                -(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*(complex_mass &
                +(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0/2+yy(4,kk)/4+yy(4,kk-1)/4)*dr(ll)/r0

        ss(4,2) = (-(yy(1,kk)/2+yy(1,kk-1)/2)*(yy(3,kk)/2+yy(3,kk-1)/2)*(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1))**2*(complex_mass+(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1)))*r0+(yy(1,kk)/2+yy(1,kk-1)/2)*(yy(3,kk)/2+yy(3,kk-1)/2) &
                *(complex_mass-(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1))**2*r0)*dr(ll)/r0

        ss(4,3) = (-4*(yy(1,kk)/2+yy(1,kk-1)/2)*complex_mass**2*(yy(3,kk)/2+yy(3,kk-1)/2)*r0**2*smallpi &
                *(yy(4,kk)/2+yy(4,kk-1)/2)-(yy(1,kk)/2+yy(1,kk-1)/2)*(complex_mass-(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1)))*(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0/2)*dr(ll)/r0

        ss(4,4) = ((-2*(yy(1,kk)/2+yy(1,kk-1)/2)*complex_mass**2*(yy(3,kk)/2+yy(3,kk-1)/2)**2*r0**2*smallpi &
                +yy(1,kk)/4+yy(1,kk-1)/4+1.D0/2.D0)*dr(ll)-r0)/r0

        ss(4,5) = ((yy(1,kk)/2+yy(1,kk-1)/2)*(yy(3,kk)/2+yy(3,kk-1)/2)/(yy(2,kk)+yy(2,kk-1)) &
               *(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0-(yy(1,kk)/2+yy(1,kk-1)/2) &
               *(yy(3,kk)/2+yy(3,kk-1)/2)*(complex_mass-(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))) &
               /(yy(2,kk)+yy(2,kk-1))*r0)*dr(ll)/r0

        ss(4,6) = (-2*complex_mass**2*(yy(3,kk)/2+yy(3,kk-1)/2)**2*r0**2*smallpi*(yy(4,kk)/2+yy(4,kk-1)/2) &
                -(yy(3,kk)/2+yy(3,kk-1)/2)*(complex_mass-(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))) &
                *(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0/2+yy(4,kk)/4+yy(4,kk-1)/4)*dr(ll)/r0

        ss(4,7) = (-(yy(1,kk)/2+yy(1,kk-1)/2)*(yy(3,kk)/2+yy(3,kk-1)/2)*(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1))**2*(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0 &
                +(yy(1,kk)/2+yy(1,kk-1)/2)*(yy(3,kk)/2+yy(3,kk-1)/2)*(complex_mass-(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1)))*(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))**2*r0)*dr(ll)/r0

        ss(4,8) = (-4*(yy(1,kk)/2+yy(1,kk-1)/2)*complex_mass**2*(yy(3,kk)/2+yy(3,kk-1)/2)*r0**2*smallpi &
                *(yy(4,kk)/2+yy(4,kk-1)/2)-(yy(1,kk)/2+yy(1,kk-1)/2)*(complex_mass-(yy(5,kk)+yy(5,kk-1)) &
                /(yy(2,kk)+yy(2,kk-1)))*(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0/2)*dr(ll)/r0

        ss(4,9) = ((-2*(yy(1,kk)/2+yy(1,kk-1)/2)*complex_mass**2*(yy(3,kk)/2+yy(3,kk-1)/2)**2*r0**2*smallpi &
                +yy(1,kk)/4+yy(1,kk-1)/4+1.D0/2.D0)*dr(ll)+r0)/r0

        ss(4,10) = ((yy(1,kk)/2+yy(1,kk-1)/2)*(yy(3,kk)/2+yy(3,kk-1)/2)/(yy(2,kk)+yy(2,kk-1)) &
                 *(complex_mass+(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1)))*r0-(yy(1,kk)/2+yy(1,kk-1)/2) &
                 *(yy(3,kk)/2+yy(3,kk-1)/2)*(complex_mass-(yy(5,kk)+yy(5,kk-1))/(yy(2,kk)+yy(2,kk-1))) &
                 /(yy(2,kk)+yy(2,kk-1))*r0)*dr(ll)/r0

!       Equation 5. This is the equation for the eigenvalue,
!       which is just:
!
!       dy5/dr  =  0

        ss(5,11) = yy(5,kk) - yy(5,kk-1)

!       Derivatives.

        ss(5,5 ) = - 1.d0
        ss(5,10) = + 1.d0

     end if


! *******************************************
! ***   L-PROCA STARS POLAR AREAL GAUGE   ***
! *******************************************

! In this case we have 7 equations for 7 variables:
!
! y1  =  A
! y2  =  alpha
! y3  =  procaF
! y4  =  procaA
! y5  =  procaB
! y6  =  procaG
! y7  =  omega  (eigenvalue)
!
! The matrix SS is a 7x15 matrix.  The equations are in column 15
! for each variable with all terms on the same size (that is EQ=0),
! while the derivatives with respect to y(j,kk-1) are on columns 1-7,
! and with respect to y(j,kk) on columns 8-14. The expressions below
! are obtained with MAPLE, so they look nasty.

  else if (idata=="l-procastar") then

!    Initialize matrix S.

     ss = 0.d0

!    Case cproca_l=0.

     if (cproca_l==0) then

!       Boundary conditions at first point (6 conditions).
!       For cproca_l=0 the boundary conditions correspond to:
!
!       A(r<<1)     -  ( 1 + A2 r^2)  =  0 
!       alpha(r=0)  -  1     =  0
!
!       procaF(r=0) -  phi0  =  0
!       procaA(r=0)          =  0
!
!       procaB(r=0)          =  0
!       procaG(r=0)          =  0
!
!       The last two are there to guarantee that procaB and
!       procaG remain 0 for l=0.  Notice that the value of
!       (alpha,phi) are arbitrary at the first point.  For A
!       we do a Taylor expansion to second order from the
!       Hamiltonian constraint to find the value of A2:
!
!       A2  =  (1/3) m^2 rho0^2
!
!       At the moment the above conditions are only first
!       order, will fix this later.
!
!       Remember that the boundary conditions at the left
!       must involve the LAST rows of SS.

        if (kk==kk1) then

!          Equations at left boundary.

           r0 = rr(kk1)

           ss(2,15) = yy(1,1) - (1.d0 + (proca_mass*proca_phi0*r0)**2/3.d0)
           ss(3,15) = yy(2,1) - 1.d0
           ss(4,15) = yy(3,1) - proca_phi0
           ss(5,15) = yy(4,1)
           ss(6,15) = yy(5,1)
           ss(7,15) = yy(6,1)

!          Derivatives. Remember that the values of
!          the variables at the first point correspond
!          to the columns (ne+1,...) since there is
!          no dependence on point kk-1.

           ss(2,8)  = 1.d0
           ss(3,9)  = 1.d0
           ss(4,10) = 1.d0
           ss(5,11) = 1.d0
           ss(6,12) = 1.d0
           ss(7,13) = 1.d0
 
!       Boundary conditions at last point (1 condition).
!
!       The boundary conditions here must correspond to the
!       FIRST rows of SS.

        else if (kk>kk2) then

!          Equations at right boundary.

           !ss(1,15) = yy(3,kk2)
           ss(7,15) = yy(7,kk2) - 0.5d0*(omega_left + omega_right)

!          Derivatives.

           !ss(1,10) = 1.d0
           ss(1,14) = 1.d0

!       Interior points.

        else

!          Calculate average radius.

           r0 = 0.5d0*(rr(kk) + rr(kk-1))

!          Equation 1. This is the Hamiltonian constraint:
!
!          dA/dr  =  A [ (1-A)/r  +  8 pi r A rho ]
!
!          with rho the energy density given by:
!
!          rho  =  1/(8 pi) [ A E^2 + m^2 ( (procaF/alpha)^2 + procaA^2/A) ]
!
!          and E the "electric field":
!
!          E  =  - alpha / (omega A) [ m^2 procaA^2 ]


!          Derivatives.


!          Equation 2. This is the lapse condition:
!
!          dalpha/dr  =  alpha [ (1-A)/2r + 4 pi r A SA ]
!
!          with SA given by:
!
!          SA  =  1/(8 pi) [ - A E^2 + m^2 ( (procaF/alpha)^2 + procaA^2/A) ]


           ss(2,15) = yy(2,kk) - yy(2,kk-1)

!          Derivatives.

           ss(2,2) = - 1.d0
           ss(2,9) = + 1.d0

!          Equation 3.  This is the equation for procaF:
!
!          dprocaF/dr  =   omega procaA [ (m alpha/omega)**2 - 1 ]

           aux1 = proca_s0/proca_phi0
           ss(3,15) = yy(3,kk) - proca_phi0*exp(-rr(kk)**2/aux1**2)

!          Derivatives.

           ss(3,10) = 1.d0

!          Equation 4.  This is the equation for procaA:
!
!          dprocaA/dr  =  omega*F*A / alpha**2  -  a [ (A+1)/r + 4 pi r A (SA - rho) ]

           aux1 = proca_s0/proca_phi0
           ss(4,15) = yy(4,kk) + (2.d0*rr(kk)/aux1**2)*proca_phi0*exp(-rr(kk)**2/aux1**2) &
                    /(proca_mass**2/proca_omega - proca_omega)

!          Derivatives.

           ss(4,11) = 1.d0

!          Equation 5. This is the equation for procaB, which
!          is just the definition for procaG:
!
!          dprocaB/dr  =  procaG

           aux1 = 0.5d0*dr(ll)

           ss(5,15) = yy(5,kk) - yy(5,kk-1) - aux1*(yy(6,kk) + yy(6,kk-1))

!          Derivatives.

           ss(5,5)  = - 1.d0
           ss(5,6)  = - aux1

           ss(5,12) = + 1.d0
           ss(5,13) = - aux1

!          Equation 6.  This is the equation for procaG, which
!          for cproca_l=0 we just take to be:
!
!          dprocaG/dr  =  0
 
           ss(6,15) = yy(6,kk) - yy(6,kk-1)

!          Derivatives.

           ss(6,6 ) = - 1.d0
           ss(6,13) = + 1.d0

!          Equation 7. This is the equation for the eigenvalue,
!          which is just:
!
!          domega/dr  =  0

           ss(7,15) = yy(7,kk) - yy(7,kk-1)

!          Derivatives.

           ss(7,7 ) = - 1.d0
           ss(7,14) = + 1.d0

        end if

!    Case cproca_l>0.

     else

     end if

  end if


! ***************
! ***   END   ***
! ***************

  end subroutine difeq
