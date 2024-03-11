
!       This code is taken directly from Numerical Recipes. It is a driver
!       routine for solving two point boundary value problems by relaxation.
!
!       itmax   =  maximum number of iterations
!       conv    =  convergence criterion
!       slowc   =  fraction of corrections used after each iteration
!       scalv   =  typical sizes for each dependent variable (to weight errors)
!       indexv  =  list the column ordering of variables used to construct the matrix S of derivatives
!
!       The problem involves ne equations for the ne ajustable dependent variables
!       at each point.  At the first point there are nb boundary conditions.
!       There are a total of m mesh points.  y(1:nyj,1:nyk) is the two dimensional
!       array that contains the initial guess for all dependent variables at each point.
!       On each iteration, y is updated by the calculated correction.
!
!       This routine calls a subroutine "difeq" that must be supplied by the user,
!
!       I added an extra parameter "l" that tells the routine on which grid level
!       we are in order to apply the correct boundary conditions. This parameter
!       must be passed to "difeq2, and the corresponding interfaces must be 
!       modified accordingly in the file "nr.f90".


	SUBROUTINE solvde(itmax,conv,slowc,scalv,indexv,nb,x,y,l)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror
	USE nr, ONLY : difeq
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: itmax,nb,l
	REAL(DP), INTENT(IN) :: conv,slowc
	REAL(DP), DIMENSION(:), INTENT(IN) :: scalv
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indexv
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: y
	INTEGER(I4B) :: ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8,&
		j9,jc1,jcf,jv,k,k1,k2,km,kp,m,ne,nvars
	INTEGER(I4B), DIMENSION(size(scalv)) :: kmax
	REAL(DP) :: err,fac
	REAL(DP), DIMENSION(size(scalv)) :: ermax
	REAL(DP), DIMENSION(size(scalv),2*size(scalv)+1) :: s
	REAL(DP), DIMENSION(size(scalv),size(scalv)-nb+1,size(y,2)+1) :: c
	ne=assert_eq(size(scalv),size(indexv),size(y,1),'solvde: ne')
	m=size(y,2)
	k1=1
	k2=m
	nvars=ne*m
	j1=1
	j2=nb
	j3=nb+1
	j4=ne
	j5=j4+j1
	j6=j4+j2
	j7=j4+j3
	j8=j4+j4
	j9=j8+j1
	ic1=1
	ic2=ne-nb
	ic3=ic2+1
	ic4=ne
	jc1=1
	jcf=ic3
	do it=1,itmax
		k=k1
		call difeq(k,k1,k2,j9,ic3,ic4,indexv,s,x,y,l)
		call pinvs(ic3,ic4,j5,j9,jc1,k1,c,s)
		do k=k1+1,k2
			kp=k-1
			call difeq(k,k1,k2,j9,ic1,ic4,indexv,s,x,y,l)
			call red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s)
			call pinvs(ic1,ic4,j3,j9,jc1,k,c,s)
		end do
		k=k2+1
		call difeq(k,k1,k2,j9,ic1,ic2,indexv,s,x,y,l)
		call red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s)
		call pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,s)
		call bksub(ne,nb,jcf,k1,k2,c)
		do j=1,ne
			jv=indexv(j)
			km=imaxloc(abs(c(jv,1,k1:k2)))+k1-1
			ermax(j)=c(jv,1,km)
			kmax(j)=km
		end do
		ermax(:)=ermax(:)/scalv(:)
		err=sum(sum(abs(c(indexv(:),1,k1:k2)),dim=2)/scalv(:))/nvars
		fac=slowc/max(slowc,err)
		y(:,k1:k2)=y(:,k1:k2)-fac*c(indexv(:),1,k1:k2)
		!print *, it,err
                write(*,"(A,I0,A,ES8.2)") ' Iteration: ',it,'    Residual: ',err
		if (err < conv) RETURN
	end do
	call nrerror('itmax exceeded in solvde')
	CONTAINS
!BL
	SUBROUTINE bksub(ne,nb,jf,k1,k2,c)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ne,nb,jf,k1,k2
	REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: c
	INTEGER(I4B) :: im,k,nbf
	nbf=ne-nb
	im=1
	do k=k2,k1,-1
		if (k == k1) im=nbf+1
		c(im:ne,jf,k)=c(im:ne,jf,k)-matmul(c(im:ne,1:nbf,k),c(1:nbf,jf,k+1))
	end do
	c(1:nb,1,k1:k2)=c(1+nbf:nb+nbf,jf,k1:k2)
	c(1+nb:nbf+nb,1,k1:k2)=c(1:nbf,jf,k1+1:k2+1)
	END SUBROUTINE bksub
!BL
	SUBROUTINE pinvs(ie1,ie2,je1,jsf,jc1,k,c,s)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ie1,ie2,je1,jsf,jc1,k
	REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: c
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: s
	INTEGER(I4B) :: i,icoff,id,ipiv,jcoff,je2,jp,jpiv,js1
	INTEGER(I4B), DIMENSION(ie2) :: indxr
	REAL(DP) :: big,piv,pivinv
	REAL(DP), DIMENSION(ie2) :: pscl
	je2=je1+ie2-ie1
	js1=je2+1
	pscl(ie1:ie2)=maxval(abs(s(ie1:ie2,je1:je2)),dim=2)
	if (any(pscl(ie1:ie2) == 0.0)) &
		call nrerror('singular matrix, row all 0 in pinvs')
	pscl(ie1:ie2)=1.0_sp/pscl(ie1:ie2)
	indxr(ie1:ie2)=0
	do id=ie1,ie2
		piv=0.0
		do i=ie1,ie2
			if (indxr(i) == 0) then
				jp=imaxloc(abs(s(i,je1:je2)))+je1-1
				big=abs(s(i,jp))
				if (big*pscl(i) > piv) then
					ipiv=i
					jpiv=jp
					piv=big*pscl(i)
				end if
			end if
		end do
		if (s(ipiv,jpiv) == 0.0) call nrerror('singular matrix in pinvs')
		indxr(ipiv)=jpiv
		pivinv=1.0_sp/s(ipiv,jpiv)
		s(ipiv,je1:jsf)=s(ipiv,je1:jsf)*pivinv
		s(ipiv,jpiv)=1.0
		do i=ie1,ie2
			if (indxr(i) /= jpiv .and. s(i,jpiv) /= 0.0) then
				s(i,je1:jsf)=s(i,je1:jsf)-s(i,jpiv)*s(ipiv,je1:jsf)
				s(i,jpiv)=0.0
			end if
		end do
	end do
	jcoff=jc1-js1
	icoff=ie1-je1
	c(indxr(ie1:ie2)+icoff,js1+jcoff:jsf+jcoff,k)=s(ie1:ie2,js1:jsf)
	END SUBROUTINE pinvs
!BL
	SUBROUTINE red(iz1,iz2,jz1,jz2,jm1,jm2,jmf,ic1,jc1,jcf,kc,c,s)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: iz1,iz2,jz1,jz2,jm1,jm2,jmf,ic1,jc1,jcf,kc
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: s
	REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: c
	INTEGER(I4B) :: ic,l,loff
	loff=jc1-jm1
	ic=ic1
	do j=jz1,jz2
		do l=jm1,jm2
			s(iz1:iz2,l)=s(iz1:iz2,l)-s(iz1:iz2,j)*c(ic,l+loff,kc)
		end do
		s(iz1:iz2,jmf)=s(iz1:iz2,jmf)-s(iz1:iz2,j)*c(ic,jcf,kc)
		ic=ic+1
	end do
	END SUBROUTINE red
	END SUBROUTINE solvde
