!taken from CAMB:
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function rombint2(f,a,b,tol, maxit, minsteps)
        use precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.

! Modified by AL to specify max iterations and minimum number of steps
! (min steps useful to stop wrong results on periodic or sharp functions)
        implicit none
        integer, parameter :: MAXITER=20,MAXJ=5
        dimension g(MAXJ+1)
        real(dl) f
        external f
        real(dl) :: rombint2
        real(dl), intent(in) :: a,b,tol
        integer, intent(in):: maxit,minsteps
     
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
      
        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
        do
          i=i+1
          if (i > maxit.or.(i > 5.and.abs(error) < tol) .and. nint > minsteps) exit
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._dl
          do k=1,nint
            g0=g0+f(a+(k+k-1)*h)
          end do
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._dl
          do j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
          end do  
          if (abs(g0).gt.tol) then
            error=1._dl-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        end do
  
        rombint2=g0
        if (i > maxit .and. abs(error) > tol)  then
          write(*,*) 'Warning: Rombint2 failed to converge; '
          write (*,*)'integral, error, tol:', rombint2,error, tol
        end if
        
end function rombint2

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function rombint(f,a,b,tol)
        use Precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, parameter :: MAXITER=20
        integer, parameter :: MAXJ=5
        dimension g(MAXJ+1)
        real(dl) f
        external f
        real(dl) :: rombint
        real(dl), intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
!

        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
            go to 40
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._dl
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._dl
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1._dl-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)  then
          write(*,*) 'Warning: Rombint failed to converge; '
          write (*,*)'integral, error, tol:', rombint,error, tol
        end if
        
end function rombint

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function rombint_obj(obj,f,a,b,tol, maxit)
        use Precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, intent(in), optional :: maxit
        integer :: MAXITER=20
        integer, parameter :: MAXJ=5
        dimension g(MAXJ+1)
        real obj !dummy
        real(dl) f
        external f
        real(dl) :: rombint_obj
        real(dl), intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
!

        if (present(maxit)) then
            MaxIter = maxit
        end if
        h=0.5d0*(b-a)
        gmax=h*(f(obj,a)+f(obj,b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
            go to 40
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._dl
            do 20 k=1,nint
            g0=g0+f(obj,a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._dl
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1._dl-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint_obj=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)  then
          write(*,*) 'Warning: Rombint failed to converge; '
          write (*,*)'integral, error, tol:', rombint_obj,error, tol
        end if
        
end function rombint_obj

subroutine spline(x,y,n,yp1,ypn,y2)
	use precision
	implicit none
	integer :: n,i,k
	real(dl) :: yp1,ypn,x(n),y(n),y2(n)
	integer, parameter :: NMAX=500
	real(dl) :: p,qn,sig,un,u(NMAX)

	if (yp1.gt..99e30) then
	   y2(1)=0.
	   u(1)=0.
	else
	   y2(1)=-0.5
	   u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	endif
	do i=2,n-1
	   sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
	   p=sig*y2(i-1)+2.
	   y2(i)=(sig-1.)/p
	   u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
	 end do
	if (ypn.gt..99e30) then
	  qn=0.
	  un=0.
	else
	  qn=0.5
	  un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
	endif
	y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
	do k=n-1,1,-1
	  y2(k)=y2(k)*y2(k+1)+u(k)
	end do
end subroutine spline

subroutine splint(xa,ya,y2a,n,x,y)
	use precision
	use ndErrors
	implicit none
	integer :: n,k,khi,klo
	real(dl) :: x,y,xa(n),y2a(n),ya(n),a,b,h

	klo=1
	khi=n
	do while (khi-klo.gt.1) 
	   k=(khi+klo)/2
	   if(xa(k).gt.x)then
		  khi=k
	   else
		  klo=k
	   endif
	end do   
	h=xa(khi)-xa(klo)
	if (h.eq.0.) call criticalError('bad xa input in sbl_splint')
	a=(xa(khi)-x)/h
	b=(x-xa(klo))/h
	y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
end subroutine splint
