module fpInterpolate
	use precision

	implicit none

	type spline_class
		integer :: n
		real(dl), dimension(:), allocatable :: xa, ya, y2a
		contains
		procedure :: initialize => spline_initialize
		procedure :: replace => spline_replace
		procedure :: evaluate => spline_interp
		procedure :: unset => spline_unset
	end type spline_class

	contains

	pure subroutine spline_initialize(obj, n, xa, ya)
		class(spline_class), intent(inout) :: obj
		integer, intent(in) :: n
		real(dl), dimension(n), intent(in) :: xa, ya
		obj%n = n
		allocate(obj%xa(n), obj%ya(n), obj%y2a(n))
		obj%xa = xa
		obj%ya = ya
		call spline(obj%xa, obj%ya, obj%n, 1.d30, 1.d30, obj%y2a) !x, y(x), N, ?, ?, derivatives)
	end subroutine spline_initialize

	pure subroutine spline_replace(obj, n, xa, ya)
		class(spline_class), intent(inout) :: obj
		integer, intent(in) :: n
		real(dl), dimension(n), intent(in) :: xa, ya
		obj%n = n
		obj%xa = xa
		obj%ya = ya
		call spline(obj%xa, obj%ya, obj%n, 1.d30, 1.d30, obj%y2a) !x, y(x), N, ?, ?, derivatives)
	end subroutine spline_replace

	pure subroutine spline_unset(obj)
		class(spline_class), intent(inout) :: obj
		deallocate(obj%xa, obj%ya, obj%y2a)
	end subroutine spline_unset

	pure function spline_interp(obj, x)
		class(spline_class), intent(in) :: obj
		real(dl), intent(in) :: x
		real(dl) :: spline_interp
		call splint(obj%xa, obj%ya, obj%y2a, obj%n, x, spline_interp) !x, y(x), derivatives, N, xout, yout)
	end function spline_interp

	pure subroutine spline(x,y,n,yp1,ypn,y2)
		use precision
		implicit none
		integer, intent(in) :: n
		real(dl), intent(in) :: x(n), y(n), yp1,ypn
		real(dl), intent(out) :: y2(n)
		integer :: i,k
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

	pure subroutine splint(xa,ya,y2a,n,x,y)
		use precision
		implicit none
		integer, intent(in) :: n
		real(dl), intent(in) :: x, xa(n),y2a(n),ya(n)
		real(dl), intent(out) :: y
		integer :: k,khi,klo
		real(dl) :: a,b,h

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
!		if (h.eq.0.) call criticalError('bad xa input in splint')
		a=(xa(khi)-x)/h
		b=(x-xa(klo))/h
		y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
	end subroutine splint
end module fpInterpolate

module utilities
	use precision
	use fpErrors
	use variables
	use fpInterpolate
	use fpInterfaces1

    integer, parameter :: GI = DL
    integer, parameter :: sp_acc = DL

	contains

	subroutine tic(t1)
		implicit none
		real(8), intent(out) :: t1
		call cpu_time(t1)
	end subroutine tic

	subroutine toc(t1, string)
		implicit none
		real(8), intent(in) :: t1
		character(len=*), intent(in) :: string
		real(8) :: t2
		call cpu_time(t2)
		write (*,*) string, " Time Taken -->", real(t2-t1)
	end subroutine toc

	pure function loglinspace(minv, cenv, maxv, numb, numlog)
		real(dl), intent(in) :: minv, maxv, cenv
		integer,  intent(in) :: numb, numlog
		real(dl), dimension(numb) :: loglinspace
		real(dl) :: dx
		integer :: i

		dx = (log10(cenv)-log10(minv)) / (numlog-1)
		do concurrent (I = 1:numlog-1)
			loglinspace(i) = 10.d0**((i-1)*dx + log10(minv))
		end do
		dx = (maxv-cenv) / (numb-numlog)
		do concurrent (I = 1:numb-numlog+1)
			loglinspace(i+numlog-1) = (i-1)*dx + cenv
		end do
		return
	end function loglinspace

	pure function linspace(minv, maxv, numb)
		real(dl), intent(in) :: minv, maxv
		integer,  intent(in) :: numb
		real(dl), dimension(numb) :: linspace

		real(dl) :: dx
		integer :: i

		dx = (maxv-minv) / (numb-1)
		do concurrent (I = 1:numb)
			linspace(i) = (i-1)*dx +minv
		end do
		return
	end function linspace

	pure function logspace(minv, maxv, numb)
		real(dl), intent(in) :: minv, maxv
		integer,  intent(in) :: numb
		real(dl), dimension(numb) :: logspace

		real(dl) :: dx
		integer :: i

		dx = (maxv-minv) / (numb-1)
		do concurrent (I = 1:numb)
			logspace(i) = 10.d0**((i-1)*dx +minv)
		end do
		return
	end function logspace

	pure function geomspace(minv, maxv, numb)
		real(dl), intent(in) :: minv, maxv
		integer,  intent(in) :: numb
		real(dl), dimension(numb) :: geomspace

		geomspace = logspace(log10(minv), log10(maxv), numb)
	end function geomspace

	elemental function E_k_m(k, m)
		real(dl), intent(in) :: k, m
		real(dl) :: E_k_m
		E_k_m = sqrt(k*k+m*m)
	end function E_k_m

	elemental function Ebare_i_dme(x, y, dme2)!for electrons
		real(dl) :: Ebare_i_dme
		real(dl), intent(in) :: x, y, dme2

		Ebare_i_dme = sqrt(x*x+y*y+dme2)
	end function Ebare_i_dme

	elemental function fermiDirac(x)
		real(dl) :: fermiDirac
		real(dl), intent(in) :: x
		fermiDirac = 1.d0/(exp(x) + 1.d0)
	end function fermiDirac

	subroutine openFile(u, fname, overwrite)
		integer :: u
		character(len=*), intent(in) :: fname
		logical :: overwrite
		character(len=5) :: tmpchar

		if (overwrite) then
			write(tmpchar,'(I5)') u
			call addToLog("[output] Opening new file "//fname//" in unit "//tmpchar)
			open(unit=u, file=fname, status="unknown", action="write")
		else
			open(unit=u, file=fname, status="old", position="append", action="write")
		end if
	end subroutine openFile

	subroutine writeCheckpoints(n, x, vars)
		integer :: n
		real(dl) :: x
		real(dl), dimension(n), intent(in) :: vars

		if (.not. checkpoint) return

		if (firstPoint) then
			firstPoint = .false.
			return
		end if

		if (verbose.gt.1) call addToLog("[ckpt] saving checkpoint")

		!save a tmp file and then move to the real one
		open(file=trim(outputFolder)//"/checkpoint.chk_tmp", unit=8643, status="unknown", form="unformatted")
		write(8643) n
		write(8643) x
		write(8643) vars(:)
		close(8643)
		call rename(trim(outputFolder)//"/checkpoint.chk_tmp",trim(outputFolder)//"/checkpoint.chk")
	end subroutine writeCheckpoints

	subroutine readCheckpoints(n, x, vars, p)
		integer, intent(out) :: n
		real(dl), intent(out) :: x
		real(dl), dimension(:), allocatable, intent(out) :: vars
		logical, intent(out) :: p

		if (.not. checkpoint) return

		inquire(file=trim(outputFolder)//"/checkpoint.chk", exist=p)

		if (.not. p) return

		call addToLog("[ckpt] reading checkpoint...")

		!save a tmp file and then move to the real one
		open(file=trim(outputFolder)//"/checkpoint.chk", unit=8643, form="unformatted")
		read(8643) n
		read(8643) x
		allocate(vars(n))
		read(8643) vars(:)
		close(8643)

		call addToLog("[ckpt] ...done!")
	end subroutine readCheckpoints

	subroutine deleteCheckpoints
		integer :: fstat
		open(file=trim(outputFolder)//"/checkpoint.chk", unit=8643, iostat=fstat, status='old')
		if (fstat == 0) close(8643, status='delete')
	end subroutine deleteCheckpoints

	pure function rombint_ri(obj1, obj2, f, a, b, tol, maxit)
		use Precision
		!  Rombint returns the integral from a to b of using Romberg integration.
		!  The method converges provided that f(x) is continuous in (a,b).
		!  f must be real(dl) and must be declared external in the calling
		!  routine.  tol indicates the desired relative accuracy in the integral.
		implicit none

		interface
			pure real(KIND(1.d0)) function f(a, b, c)
				real(KIND(1.d0)), intent(in) :: a, b
				integer, intent(in) :: c
			end function
		end interface

		integer, intent(in), optional :: maxit
		integer :: MAXITER
		integer, parameter :: MAXJ=5
		dimension g(MAXJ+1)
		real(dl), intent(in) :: obj1
		integer, intent(in) :: obj2
		real(dl) :: rombint_ri
		real(dl), intent(in) :: a,b,tol
		integer :: nint, i, k, jmax, j
		real(dl) :: h, gmax, error, g, g0, g1, fourj

		if (present(maxit)) then
			MaxIter = maxit
		else
			MAXITER=20
		end if
		h=0.5d0*(b-a)
		gmax=h*(f(obj1, a, obj2)+f(obj1, b, obj2))
		g(1)=gmax
		nint=1
		error=1.0d20
		i=0
10      i=i+1
		if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
			go to 40
		g0=0._dl
		do 20 k=1,nint
			g0=g0+f(obj1, a+(k+k-1)*h, obj2)
20      continue
		g0=0.5d0*g(1)+h*g0
		h=0.5d0*h
		nint=nint+nint
		jmax=min(i,MAXJ)
		fourj=1._dl
		do 30 j=1,jmax
			fourj=4._dl*fourj
			g1=g0+(g0-g(j))/(fourj-1._dl)
			g(j)=g0
			g0=g1
30      continue
		if (abs(g0).gt.tol) then
			error=1._dl-gmax/g0
		else
			error=gmax
		end if
		gmax=g0
		g(jmax+1)=g0
		go to 10
40      rombint_ri=g0
	end function rombint_ri

	!Newton-Cotes method:
	!approximate 1D/2D integrals using a linear interpolation inside bins
	!and perform analytical integration of the line/plane.
	!save both interpolation and some numerical integration time
	pure function integral_NC_1d(N, dx, f)
		real(dl), dimension(:), allocatable, intent(in) :: dx, f
		integer, intent(in) :: N
		real(dl) :: integral_NC_1d
		integer :: ix

		integral_NC_1d = 0.d0
		do ix=1, N-1
			integral_NC_1d = integral_NC_1d + dx(ix) * (f(ix) + f(ix+1))
		end do
		integral_NC_1d = integral_NC_1d * 0.5d0
	end function integral_NC_1d

	pure function integral_NC_2d(N1, N2, dx, dy, f)
		real(dl), dimension(:), allocatable, intent(in) :: dx, dy
		real(dl), dimension(:,:), allocatable, intent(in) :: f
		integer, intent(in) :: N1, N2
		real(dl) :: integral_NC_2d
		integer :: ix, iy

		integral_NC_2d = 0.d0
		do ix=1, N1-1
			do iy=1, N2-1
				integral_NC_2d = integral_NC_2d &
					+ dx(ix) * dy(iy) * (&
						f(ix, iy) + f(ix+1, iy) + f(ix, iy+1) + f(ix+1, iy+1)&
					)
			end do
		end do
		integral_NC_2d = integral_NC_2d * 0.25d0
	end function integral_NC_2d

	!utilities for Gauss-Laguerre quadrature
	elemental function gammln(xx)
		real(dl) :: gammln
		real(dl), intent(in) :: xx
		!Returns the value ln[Γ(xx)] for xx > 0.
		integer j
		real(dl) :: ser,tmp,x,y
		!Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
		!accuracy is good enough.
		real(dl), parameter :: stp = 2.5066282746310005d0
		real(dl), dimension(6), parameter :: cof = (/76.18009172947146d0,&
		-86.50532032941677d0,24.01409824083091d0,&
		-1.231739572450155d0,0.1208650973866179d-2,&
		-0.5395239384953d-5/)

		x=xx
		y=x
		tmp=x+5.5d0
		tmp=(x+0.5d0)*log(tmp)-tmp
		ser=1.000000000190015d0
		do j=1,6
			y=y+1.d0
			ser=ser+cof(j)/y
		enddo
		gammln=tmp+log(stp*ser/x)
		return
	end function gammln

	subroutine gaulag(x,w,n,alf)
		integer, intent(in) :: n
		real(dl), intent(in) :: alf
		real(dl), intent(out), dimension(:), allocatable :: w, x
		integer, parameter :: MAXIT=20
		real(qp), parameter :: eps=1.d-15
		! Increase EPS if you don't have this precision.
		!Given alf , the parameter α of the Laguerre polynomials, this routine returns arrays x(1:n)
		!and w(1:n) containing the abscissas and weights of the n-point Gauss-Laguerre quadrature
		!formula. The smallest abscissa is returned in x(1), the largest in x(n) .
		integer :: i,its,j
		real(dl) :: nd
		real(qp) :: ai
		real(qp) :: p1,p2,p3,pp,z,z1
		nd = n
		allocate(x(n), w(n))
		do i=1,n	! Loop over the desired roots.
			if(i.eq.1)then	! Initial guess for the smallest root.
				z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
			else if(i.eq.2)then	! Initial guess for the second root.
				z=z+(15.+6.25*alf)/(1.+.9*alf+2.5*n)
			else	! Initial guess for the other roots.
				ai=i-2
				z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))*(z-x(i-2))/(1.+.3*alf)
			endif
			do its=1,MAXIT	! Refinement by Newton's method.
				p1=1.d0
				p2=0.d0
				do j=1,n	! Loop up the recurrence relation to get the Laguerre polynomial evaluated at z.
					p3=p2
					p2=p1
					p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
				enddo
				!p1 is now the desired Laguerre polynomial. We next compute pp, its derivative, by
				!a standard relation involving also p2, the polynomial of one lower order.
				pp=(n*p1-(n+alf)*p2)/z
				z1=z
				z=z1-p1/pp ! Newton's formula.
				if(abs(z-z1).le.EPS)goto 1
			enddo
			call criticalError("too many iterations in gaulag")
1			x(i)=z ! Store the root and the weight.
			w(i)=-exp(gammln(alf+n)-gammln(nd))/(pp*n*p2)
		enddo
		return
	end subroutine

	subroutine get_GLq_vectors(nreal, yv, wv, wv2, verb, alpha, ycut)
		real(dl), dimension(:), allocatable, intent(out) :: yv, wv, wv2
		integer, intent(in) :: nreal, alpha
		logical, intent(in) :: verb
		real(dl), dimension(:), allocatable :: tyv, twv
		real(dl), intent(in) :: ycut
		integer :: ix, iy, effective_Ny
		character(len=300) :: tmpstr

		do ix=1, 1500
			effective_Ny = 0
			call gaulag(tyv, twv, ix, 1.d0*alpha)
			do iy=1, ix
				if (tyv(iy).gt.ycut)then
					effective_Ny = iy-1
					exit
				end if
			end do
			if (nreal .le. effective_Ny) then
				if (verb) then
					write(tmpstr, "('[config] use Gauss-Laguerre, n=',I4,' and selecting the first ',I3,' roots')") ix, effective_Ny
					call addToLog(trim(tmpstr))
				end if
				if (.not. allocated(yv)) &
					allocate(yv(effective_Ny))
				if (.not. allocated(wv)) &
					allocate(wv(effective_Ny))
				yv = tyv(1:effective_Ny)
				wv = twv(1:effective_Ny)
				do iy=1, nreal
					wv(iy) = wv(iy)*exp(yv(iy))
				end do
				if (.not.allocated(wv2)) &
					allocate(wv2(Ny))
				do iy=1, Ny
					wv2(iy) = wv(iy) / yv(iy)**alpha
				end do
				deallocate(tyv, twv)
				return
			end if
			deallocate(tyv, twv)
		end do
		write(tmpstr, "('Cannot find Ny=',I3,' Gauss-Laguerre nodes below y_max=',E11.4,'. Reached N=',I4)") nreal, ycut, ix-1
		call criticalError(trim(tmpstr))
	end subroutine get_GLq_vectors

	!Gauss-Laguerre quadrature: very effective for exponentially-suppressed functions
	pure function integral_GL_1d(wx, f)
		real(dl), dimension(:), allocatable, intent(in) :: wx, f
		real(dl) :: integral_GL_1d

		integral_GL_1d = sum(wx(:) * f(:))
	end function integral_GL_1d

	pure function integral_GL_2d(Nx, wx, wy, f)
		real(dl), dimension(:), allocatable, intent(in) :: wx, wy
		real(dl), dimension(:,:), allocatable, intent(in) :: f
		integer, intent(in) :: Nx
		real(dl) :: integral_GL_2d
		integer :: ix

		integral_GL_2d = 0.d0
		do ix=1, Nx
			integral_GL_2d = integral_GL_2d &
				+ sum(wx(ix) * wy(:) * f(ix, :))
		end do
	end function integral_GL_2d

	pure function integrate_ftqed_ln(func, x, z)
		real(dl) :: integrate_ftqed_ln
		real(dl), intent(in) :: x, z
		procedure (ftqed_ln_integrand) :: func
		real(dl) :: y, k
		integer :: ix, iy
		real(dl), dimension(:,:), allocatable :: fval

		allocate(fval(ln_2dint_Npts, ln_2dint_Npts))
		do ix=1, ln_2dint_Npts
			do iy=1, ln_2dint_Npts
				y = (ln_2dint_y(ix) + ln_2dint_y(iy)) / 2.
				k = (ln_2dint_y(ix) - ln_2dint_y(iy)) / 2.
				if (k.ge.0 .and. y.ge.0) then
					fval(ix, iy) = func(x, z, y, k)
				else
					fval(ix, iy) = 0.d0
				end if
			end do
		end do
		integrate_ftqed_ln = integral_NC_2d(&
			ln_2dint_Npts, &
			ln_2dint_Npts, &
			ln_2dint_dy, &
			ln_2dint_dy, &
			fval &
			)
		do ix=1, ln_2dint_Npts
			do iy=1, ln_2dint_Npts
				y = (ln_2dint_y(ix) - ln_2dint_y(iy)) / 2.
				k = (ln_2dint_y(ix) + ln_2dint_y(iy)) / 2.
				if (k.ge.0 .and. y.ge.0) then
					fval(ix, iy) = func(x, z, y, k)
				else
					fval(ix, iy) = 0.d0
				end if
			end do
		end do
		integrate_ftqed_ln = integrate_ftqed_ln + integral_NC_2d(&
			ln_2dint_Npts, &
			ln_2dint_Npts, &
			ln_2dint_dy, &
			ln_2dint_dy, &
			fval &
			)
		deallocate(fval)
		integrate_ftqed_ln = integrate_ftqed_ln / 2.d0
	end function integrate_ftqed_ln

	!Interpolation of nu density matrix
	pure function get_interpolated_nudens(ndmv, y, nf, ny) result(newmat)
		type(cmplxMatNN) :: newmat
		type(cmplxMatNN), dimension(:), allocatable, intent(in) :: ndmv
		real(dl), intent(in) :: y
		integer, intent(in) :: nf, ny
		integer :: i, j, iy
		real(dl) :: fd0, fd1, fdc

		call allocateCmplxMat(newmat)
		newmat%y = y
		newmat%re = 0.d0
		newmat%im = 0.d0

		iy = 0
		do i=1,ny
			if (y.gt.ndmv(i)%y) &
				iy=i
		end do
		if(iy.lt.1 .or. iy.ge.ny) &
			return
		fd0 = fermiDirac(ndmv(iy)%y)
		fd1 = fermiDirac(ndmv(iy+1)%y)
		fdc = fermiDirac(y)
		do i=1, nf
			newmat%re(i,i) = ( &
				ndmv(iy)%re(i,i)/fd0 &
				+ (y-ndmv(iy)%y) &
					* (ndmv(iy+1)%re(i,i)/fd1 - ndmv(iy)%re(i,i)/fd0) &
					/ (ndmv(iy+1)%y-ndmv(iy)%y) &
				) * fdc
			do j=i+1, nf
#ifdef RHO_OFFDIAG_INTERP_DIV_FD
				newmat%re(i,j) = ( &
					ndmv(iy)%re(i,j)/fd0 &
					+ (y-ndmv(iy)%y) &
						* (ndmv(iy+1)%re(i,j)/fd1 - ndmv(iy)%re(i,j)/fd0) &
						/ (ndmv(iy+1)%y-ndmv(iy)%y) &
				) * fdc
				newmat%im(i,j) = ( &
					ndmv(iy)%im(i,j) &
					+ (y-ndmv(iy)%y) &
						* (ndmv(iy+1)%im(i,j)/fd1 - ndmv(iy)%im(i,j)/fd0) &
						/ (ndmv(iy+1)%y-ndmv(iy)%y) &
				) * fdc
#else
				newmat%re(i,j) = &
					ndmv(iy)%re(i,j) &
					+ (y-ndmv(iy)%y) * (ndmv(iy+1)%re(i,j)-ndmv(iy)%re(i,j))/(ndmv(iy+1)%y-ndmv(iy)%y)
				newmat%im(i,j) = &
					ndmv(iy)%im(i,j) &
					+ (y-ndmv(iy)%y) * (ndmv(iy+1)%im(i,j)-ndmv(iy)%im(i,j))/(ndmv(iy+1)%y-ndmv(iy)%y)
#endif
				newmat%re(j,i) = newmat%re(i,j)
				newmat%im(j,i) = -newmat%im(i,j)
			end do
		end do
	end function
end module utilities
