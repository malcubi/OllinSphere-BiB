!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/tools/FFT.f90,v 1.11 2024/10/02 18:20:10 malcubi Exp $

  program FFT

! This program calculates the discete Fourier
! transform of a set of reat data, using
! the FFT (Fast Fourier Transform) routines
! from Numerical Recipes.
!
! It adjusts the number of points to the largest
! power of 2 smaller than the total.


! **************************************
! ***   DISCRETE FOURIER TRANSFORM   ***
! **************************************

! Declare variables.

  implicit none

  integer i                   ! Counters.
  integer N                   ! Total number of data points.
  integer NMAX                ! Maximum number of data points.
  integer units               ! Units used for output.

  real(8) x                   ! Auxiliary variables.
  real(8) smallpi,twopi       ! Pi and 2*Pi.
  real(8) fac                 ! Factor used for rescaling the Power.

  real(8) G,c,h,hbar          ! Newtons´s constant, speed of light and Planck's constant in SI.
  real(8) Msol                ! Solar mass in SI.

  real(8), allocatable, dimension (:) :: t,data  ! Data arrays for time and function f(t).
  real(8), allocatable, dimension (:) :: cr,ci   ! Real and Imaginary Fourier coefficients.
  real(8), allocatable, dimension (:) :: p       ! Power spectrum.
  real(8), allocatable, dimension (:) :: rdata   ! Reconstructed data.

  character(100) datafile     ! Name of data file.
  character(100) line         ! Name of data file.


! *******************
! ***   NUMBERS   ***
! *******************

  NMAX = 1048576  ! 2**20

  smallpi = acos(-1.d0)
  twopi = 2.d0*smallpi

! Physical constants (useful for rescaling to SI units).

  G = 6.6743d-11
  c = 2.99792458d8
  h = 6.62607015d-34

  hbar = h/twopi
  Msol = 1.988475d30


! *********************
! ***   DATA FILE   ***
! *********************

  print *
  print *, 'Reading data file ...'

! Get name of data file.

  call getarg(1,datafile)

  if (datafile==" ") then
     print *
     print *, 'Missing datafile name.'
     print *, 'Aborting!'
     print *
     stop
  end if

! Normalization.

  print *
  print *, 'Choose normalization for power spectrum output:'
  print *, '1. Direct (cycles per second).'
  print *, '2. Angular (radians per second).'
  print *, '3. Convert from Planck units to physical units (cycles per second).'
  print *, '4. Convert from units c=G=Msol=1 to physical units (cycles per sencond).'

  read(*,*) units

  if (units==1) then
     fac = 1.d0
  else if (units==2) then
     fac = twopi
  else if (units==3) then
     fac = sqrt(c**5/(hbar*G))
  else if (units==4) then
     fac = c**3/(G*Msol) 
  else
     print *
     print *, 'Choice of normalization out of range, aborting.'
     print *
     stop
  end if


! *********************
! ***   READ DATA   ***
! *********************

! We assume the data has two columns, one for
! the time and one for the data.

! Count number of data points.

  open(1,file=datafile)

  N = 0

  do while (N<NMAX)

!    Read data as text and ignore comments.

     read(1,'(a)',END=10) line
     line = adjustl(line)

     if (line(1:1)/="#") then
        N=N+1
     end if

  end do

  10 continue

  close(1)

  print *
  print *,'Number of data points: ',N

! Adjust N to largest power of 2 smaller
! than the size of the data set.

  x = log(dble(N))/log(2.d0)
  N = 2**int(x)

  print *,'Adjusted to power of 2:',N

! Allocate arrays.

  allocate(t(1:N),data(1:N),rdata(1:N))
  allocate(cr(0:N/2),ci(0:N/2),p(0:N/2))

! Read data from file.

  open(1,file=datafile)

  i = 0

  do while(i<N)

!    Read data first as text and ignore
!    comments, then convert to real numbers.

     read(1,'(a)',END=10) line
     line = adjustl(line)

     if (line(1:1)/="#") then
        i = i + 1
        read(line,*) t(i),data(i)
        !print *,i,t(i),data(i)
     end if

  end do

  close(1)


! ******************************************
! ***   CALL FOURIER TRANSFORM ROUTINE   ***
! ******************************************

! Call routine realft. This routine the calls
! four1 (the FFT routine).
!
! On exit the array "data" now has the Fourier
! coefficients.

  call realft(data,N,+1)


! ********************************
! ***   EXTRACT COEFFICIENTS   ***
! ********************************

! On exit data(1) has the purely real k=0 coeffcient,
! and data(2) the purely real k=N/2 coefficient.

  cr(0) = data(1)
  ci(0) = 0.d0

  cr(N/2) = data(2)
  ci(N/2) = 0.d0

! The rest of the coefficients come in pairs, real part
! on even elements of "data", and imaginary part on odd.

  do i=1,N/2-1
     cr(i) = data(2*i+1)
     ci(i) = data(2*i+2)
  end do

! Find power spectrum (square root of norm).

  do i=0,N/2
     p(i) = dsqrt(cr(i)**2 + ci(i)**2)
  end do

! Average. The average of the data is given by the
! coefficient cr(0)/N.

  print *
  print *, 'Function average cr(0)/N: ', cr(0)/dble(N)


! **********************************
! ***   RECONSTRUCTED FUNCTION   ***
! **********************************

! Reconstruct original data (for testing).
! We call again realft for the inverse
! transform.

  call realft(data,N,-1)
  rdata = 2.d0*data/dble(N)


! ******************
! ***   OUTPUT   ***
! ******************

! The Fourier coefficients and power spectrum data are
! saved with the k=0 coefficient set to 0 to eliminate
! the average (which can be quite large).
!
! Notice that ci(0) is always zero, and cr(0)/N
! is the average of the data which we already gave
! as output above.
!
! Also, the first column instead of being directly k,
! is transformed to the 'physical' frequency omega
! given by:
!
! omega  =  2 pi k / (T_final - T_initial)

  cr(0) = 0.d0
  p(0)  = 0.d0

  print *
  print *, 'Saving Fourier transform data ... '

! Save real coefficients.

  open(1,file='CoeffRe.dat',form='formatted',status='replace')

  do i=0,N/2
     write(1,"(2ES20.10)") twopi*dble(i)/(t(N)-t(1)),cr(i)
     !write(1,"(2ES20.10)") dble(i),cr(i)
  end do

  close(1)

! Save imaginary coefficients.

  open(1,file='CoeffIm.dat',form='formatted',status='replace')

  do i=0,N/2
     write(1,"(2ES20.10)") twopi*dble(i)/(t(N)-t(1)),ci(i)
     !write(1,"(2ES20.10)") dble(i),ci(i)
  end do

  close(1)

! Save power spectrum.

  open(1,file='Power.dat',form='formatted',status='replace')

  do i=0,N/2
     write(1,"(2ES25.15)") fac*dble(i)/(t(N)-t(1)),p(i)
     !write(1,"(2ES20.10)") dble(i),p(i)
  end do

  close(1)

! Save reconstructed data.

  open(1,file='rdata.dat',form='formatted',status='replace')

  do i=1,N
     write(1,"(2ES20.10)") t(i),rdata(i)
  end do

  close(1)


! ***************
! ***   END   ***
! ***************

  print *
  print *, 'DONE'
  print *

  end program FFT









  SUBROUTINE realft(data,n,isign)

! Numerical Recipes routine realtf.

! Calculates the Fourier transform of a set of n real-valued data points.
! Replaces this data (which is stored in array data(1:n)) by the positive
! frequency half of its complex Fourier transform. The real-valued first
! and last components of the complex transform are returned as elements
! data(1) and data(2), respectively.
!
! n must be a power of 2.
!
! This routine also calculates the inverse transform of a complex data array if it is
! the transform of real data. (Result in this case must be multiplied by 2/n.)

  INTEGER isign,n
  INTEGER i,i1,i2,i3,i4,n2p3

  REAL(8) data(n)
  REAL(8) c1,c2,h1i,h1r,h2i,h2r,wis,wrs
  REAL(8) theta,wi,wpi,wpr,wr,wtemp

! Numbers.

  theta = 3.141592653589793d0/dble(n/2)
  c1 = 0.5d0

! Forward transform.

  if (isign.eq.1) then

     c2 = - 0.5d0
     call four1(data,n/2,+1)

! Set up for inverse transform.

  else

     c2 = + 0.5d0
     theta = - theta

  end if

! Loop over array to rearange data.
! Case i=1 done separately below.

  wpr = - 2.0d0*sin(0.5d0*theta)**2
  wpi = sin(theta)

  wr = 1.0d0 + wpr
  wi = wpi

  n2p3 = n + 3

  do i=2,n/4

     i1 = 2*i - 1
     i2 = i1 + 1
     i3 = n2p3 - i2
     i4 = i3 + 1

     !wrs = sngl(wr)
     !wis = sngl(wi)

     wrs = wr
     wis = wi

!    The two separate transforms are separated out.

     h1r =  c1*(data(i1) + data(i3))
     h1i =  c1*(data(i2) - data(i4))
     h2r = -c2*(data(i2) + data(i4))
     h2i =  c2*(data(i1) - data(i3))

!    Here they are recombined to form the true
!    transform of the original data.

     data(i1) =  h1r + wrs*h2r - wis*h2i
     data(i2) =  h1i + wrs*h2i + wis*h2r
     data(i3) =  h1r - wrs*h2r + wis*h2i
     data(i4) = -h1i + wrs*h2i + wis*h2r

!    The recurrence.

     wtemp = wr
     wr = wr*wpr - wi*wpi + wr
     wi = wi*wpr + wtemp*wpi + wi

  end do

! Squeeze the first and last data together
! to get them all within the original array.

  if (isign.eq.1) then

     h1r = data(1)

     data(1) = h1r + data(2)
     data(2) = h1r - data(2)

! Inverse transform in case isign=-1.

  else

    h1r = data(1)

    data(1) = c1*(h1r + data(2))
    data(2) = c1*(h1r - data(2))

    call four1(data,n/2,-1)

  end if

  return
  END









  SUBROUTINE four1(data,nn,isign)

! Numerical Recipes routine four1 (fast Fourier transform).

! Replaces data(1:2*nn) by its discrete Fourier transform if isign is input as 1,
! or replaces data(1:2*nn) by nn times its inverse discrete Fourier transform if
! isign is input as −1. data is a complex array of length nn or, equivalently,
! a real array of length 2*nn.
!
! nn MUST be an integer power of 2 (this is not checked!).

  INTEGER isign,nn
  INTEGER i,istep,j,m,mmax,n

  REAL(8) data(2*nn)
  REAL(8) tempi,tempr
  REAL(8) theta,wi,wpi,wpr,wr,wtemp 

  n = 2*nn
  j = 1

! This is the bit-reversal section of the routine.

  do i=1,n,2

!    This is the bit-reversal section of the routine.

     if (j.gt.i) then

!       Exchange two complex numbers.

        tempr = data(j)
        tempi = data(j+1)

        data(j)   = data(i)
        data(j+1) = data(i+1)
        data(i)   = tempr
        data(i+1) = tempi

     end if

     m = n/2

     1 continue

     if ((m.ge.2).and.(j.gt.m)) then
        j = j - m
        m = m/2
        goto 1
     end if

     j = j + m

  end do

! Here begins the Danielson-Lanczos section of the routine.
! Outer loop executed log2 nn times.

  mmax = 2

  2 continue

  if (n.gt.mmax) then

     istep = 2*mmax
     theta = 6.28318530717959d0/dble(isign*mmax)

     wpr = - 2.d0*dsin(0.5d0*theta)**2
     wpi = dsin(theta)

     wr = 1.d0
     wi = 0.d0

!    Here are the two nested inner loops.

     do m=1,mmax,2
        do i=m,n,istep

           j = i + mmax

!          This is the Danielson-Lanczos formula.

           tempr = wr*data(j)   - wi*data(j+1)
           tempi = wr*data(j+1) + wi*data(j)

           data(j)   = data(i)   - tempr
           data(j+1) = data(i+1) - tempi
           data(i)   = data(i)   + tempr
           data(i+1) = data(i+1) + tempi

        end do

!       Trigonometric recurrence.

        wtemp = wr

        wr = wr*wpr - wi*wpi + wr
        wi = wi*wpr + wtemp*wpi + wi

     end do

     mmax = istep

!    Not yet done.

     goto 2

  end if

  return
  END


