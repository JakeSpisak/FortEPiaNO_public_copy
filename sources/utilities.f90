module utilities
	use precision
	use ndErrors
	use variables
	
    integer, parameter :: GI = DL
    integer, parameter :: sp_acc = DL
    
	abstract interface
	  function funcX (x)
		use precision
		 real(dl) :: funcX
		 real(dl), intent (in) :: x
	  end function funcX
	end interface
    
	abstract interface
	  function func2_X (x)
		use precision
		 real(dl), dimension(2) :: func2_X
		 real(dl), intent (in) :: x
	  end function func2_X
	end interface
    
	abstract interface
	  function funcXY (x,y)
		use precision
		 real(dl) :: funcXY
		 real(dl), intent (in) :: x,y
	  end function funcXY
	end interface
    
	abstract interface
	  function funcXYZ (x,y,z)
		use precision
		 real(dl) :: funcXYZ
		 real(dl), intent (in) :: x,y,z
	  end function funcXYZ
	end interface
	
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

	function linspace(minv, maxv, numb)
		real(dl), intent(in) :: minv, maxv
		integer,  intent(in) :: numb
		real(dl), dimension(numb) :: linspace
		
		real(dl) :: dx
		integer :: i
		
		dx = (maxv-minv) / (numb-1)
		do i=1, numb
			linspace(i) = (i-1)*dx +minv
		end do
		return
	end function linspace
	
	function logspace(minv, maxv, numb)
		real(dl), intent(in) :: minv, maxv
		integer,  intent(in) :: numb
		real(dl), dimension(numb) :: logspace
		
		real(dl) :: dx
		integer :: i
		
		dx = (maxv-minv) / (numb-1)
		do i=1, numb
			logspace(i) = 10.d0**((i-1)*dx +minv)
		end do
		return
	end function logspace
	
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
	
	subroutine writeCheckpoints(n,x,vars)
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
		open(file=trim(outputFolder)//"/checkpoint.chk_tmp", unit=8643, status="unknown", form="binary")
		write(8643) n
		write(8643) x
		write(8643) vars(:)
		close(8643)
		call rename(trim(outputFolder)//"/checkpoint.chk_tmp",trim(outputFolder)//"/checkpoint.chk")
	end subroutine writeCheckpoints
	
	subroutine readCheckpoints(n,x,vars,p)
		integer, intent(out) :: n
		real(dl), intent(out) :: x
		real(dl), dimension(:), allocatable, intent(out) :: vars
		logical, intent(out) :: p
		
		if (.not. checkpoint) return
		
		inquire(file=trim(outputFolder)//"/checkpoint.chk", exist=p)
		
		if (.not. p) return
		
		call addToLog("[ckpt] reading checkpoint...")
		
		!save a tmp file and then move to the real one
		open(file=trim(outputFolder)//"/checkpoint.chk", unit=8643, form="binary")
		read(8643) n
		read(8643) x
		allocate(vars(n))
		read(8643) vars(:)
		close(8643)
		
		call addToLog("[ckpt] ...done!")
	end subroutine readCheckpoints
	
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

SUBROUTINE D01GCF(N,F,REGION,NPTS,VK,NRAND,ITRANS,RES,ERRr,IFAIL,obj)
!use omp_lib!SG
!C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
!C     MARK 11.5(F77) REVISED. (SEPT 1985.)
!C
!C     MULTIPLE INTEGRATION
!C
!C     THIS SUBROUTINE CALCULATES AN APPROXIMATION TO
!C            D1      DN
!C            I ......I F(X1,...,XN) DX1 .. DXN
!C            C1      CN
!C
!C     USING THE KOROBOV-CONROY NUMBER-THEORETIC METHOD.
!C     (N.M.KOROBOV,NUMBER THEORETIC METHODS IN APPROXIMATE ANALYSIS,
!C     FIZMATGIZ,MOSCOW,1963, H.CONROY,J.CHEM.PHYS.47,(1967),
!C     5307-5813)
!C
!C     INPUT ARGUMENTS
!C
!C     N        -INTEGER SPECIFYING THE NUMBER OF DIMENSIONS.
!C               1 .LE. N .LE. 20.
!C
!C     F        -EXTERNALLY DECLARED REAL USER FUNCTION INTEGRAND
!C               HAVING ARGUMENTS (N,X) WHERE X IS A REAL ARRAY
!C               OF DIMENSION N WHOSE ELEMENTS ARE THE VARIABLE
!C               VALUES.
!C
!C     REGION   -EXTERNALLY DECLARED USER SUBROUTINE WITH ARGUMENTS
!C               (N,X,J,C,D) WHICH CALCULATES THE LOWER LIMIT C
!C               AND THE UPPER LIMIT D CORRESPONDING TO THE ARRAY
!C               VARIABLE X(J). C AND D MAY DEPEND ON X(1)..X(J-1).
!C
!C     NPTS     -INTEGER VARIABLE WHICH SPECIFIES THE KOROBOV RULE
!C               TO BE USED. THERE ARE TWO ALTERNATIVES DEPENDING
!C               ON THE VALUE OF NPTS
!C               (A) 1 .LE. NPTS .LE. 6
!C                   IN THIS CASE ONE OF THE SIX PRESET RULES IS
!C                   CHOSEN USING 2129, 5003, 10007, 20011, 40009 OR
!C                   80021 POINTS ACCORDING TO THE RESPECTIVE VALUE
!C                   OF NPTS.
!C               (B) NPTS .GT. 6
!C                   NPTS IS THE NUMBER OF ACTUAL POINTS TO BE USED
!C                   WITH CORRESPONDING OPTIMAL COEFFICIENTS SUPPLIED
!C                   BY THE USER IN ARRAY VK.
!C
!C     VK       -REAL ARRAY CONTAINING N OPTIMAL COEFFICIENTS
!C               FOR THE CASE NPTS .GT. 6.
!C
!C     NRAND    -INTEGER NUMBER SPECIFYING THE NUMBER OF RANDOM
!C               SAMPLES TO BE GENERATED IN THE ERROR ESTIMATION.
!C               (GENERALLY A SMALL VALUE ,SAY 3 TO 5, IS SUFFICIENT)
!C               THE TOTAL NUMBER OF INTEGRAND EVALUATIONS WILL BE
!C               NRAND*NPTS.
!C
!C     ITRANS   -THIS SHOULD BE SET TO 0 OR 1. THE NORMAL SETTING
!C               IS 0. WHEN SET TO 1 THE PERIODIZING TRANSFORMATION
!C               IS SUPPRESSED (TO COVER CASES WHERE THE INTEGRAND IS
!C               ALREADY PERIODIC OR WHERE THE USER DESIRES TO
!C               SPECIFY A PARTICULAR TRANSFORMATION IN THE
!C               DEFINITION OF F).
!C
!C     IFAIL    -INTEGER VARIABLE SPECIFYING THE ERROR REPORTING
!C               OPTION. IFAIL IS SET TO 0 FOR HARD FAIL AND TO 1 FOR
!C               SOFT FAIL REPORT.
!C
!C
!C     OUTPUT ARGUMENTS
!C
!C     VK       -WHEN THE SUBROUTINE IS CALLED WITH 1 .LE. NPTS
!C               .LE. 6 THEN N ELEMENTS OF ARRAY VK WILL BE SET
!C               TO THE OPTIMAL COEFFICIENTS USED BY THE PRESET RULE.
!C
!C     RES      -APPROXIMATION TO THE VALUE OF THE INTEGRAL.
!C
!C     ERRr     -STANDARD ERROR AS COMPUTED FROM NRAND SAMPLE VALUES.
!C               IF NRAND IS SET TO 1 THEN ERR WILL BE RETURNED WITH
!C               VALUE 0.
!C
!C     IFAIL    -THIS REPORTS THE ERROR CONDITIONS AS FOLLOWS
!C               IFAIL=0  NO ERROR DETECTED.
!C                    =1  N FAILS TO SATISFY 1 .LE. N .LE. 20
!C                    =2  NPTS .LT. 1.
!C                    =3  NRAND .LT. 1.
!C
!C     ************* IMPLEMENTATION-DEPENDENT CODING ****************
!C     THIS ROUTINE REQUIRES AT ONE POINT TO PERFORM EXACT INTEGER
!C     ARITHMETIC WITH AT LEAST 32 BITS AVAILABLE FOR INTERMEDIATE
!C     RESULTS. THE INTEGERS ARE STORED AS FLOATING-POINT VARIABLES
!C     (IN ORDER TO AVOID INTEGER OVERFLOW ON SOME MACHINES). IF THE
!C     MACHINE DOES NOT PROVIDE AT LEAST 32 BITS OF BASIC FLOATING-
!C     POINT PRECISION, THEN ONE STATEMENT IN THE CODE MUST BE
!C     CONVERTED TO ADDITIONAL PRECISION (SEE COMMENTS IN CODE).
!C     **************************************************************
!C
!C     REGION
!C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01GCF')
!C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERRr, RES
      INTEGER           IFAIL, ITRANS, N, NPTS, NRAND
!C     .. Array Arguments ..
      DOUBLE PRECISION  VK(N)
!C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
!C     .. Subroutine Arguments ..
      EXTERNAL          REGION
!C     .. Local Scalars ..
      DOUBLE PRECISION  C, D, DPK, GRAD, PK, PTS, SUM1, SUM2, WT, XA, XX
		type(coll_args) :: obj
!      Double precision obj !-SG-dummy arg
      INTEGER           I, ICHECK, J, K, M, MAXDIM, MAXRUL, NDIM, NP
!C     .. Local Arrays ..
      DOUBLE PRECISION  ALPHA(20), KPTS(6), KR(6,20), X(20)
      CHARACTER*1       P01REC(1)
!C     .. External Functions ..
      DOUBLE PRECISION  G05CAF
      INTEGER           P01ABF
!!!!      EXTERNAL          G05CAF, P01ABF
!C     .. Intrinsic Functions ..
      INTRINSIC         ABS, AINT, MOD, DBLE, SQRT
!C     .. Data statements ..
!C
!C     KOROBOV COEFFICIENTS FOR 1 TO 20 DIMENSIONS.
!C
!C     (ACKNOWLEDGEMENT-CALCULATED  ON THE CRAY I COMPUTER AT SANDIA
!C     NATIONAL LABORATORY, LIVERMORE, CALIFORNIA, USA, THROUGH
!C     THE KIND ASSISTANCE OF DR. R.E. HUDDLESTON)
!C
!C     KR(I,J) CONTAINS KOROBOV COEFFICIENT FOR RULE I DIMENSION J.
!C     ASSOCIATED NO. OF POINTS IS GIVEN BY KPTS(I)
!C
!C
!C     KPTS(I) = NO. OF POINTS ASSOCIATED WITH RULE I.
!C
      DATA              MAXRUL, MAXDIM/6, 20/
      DATA              KR(1,1), KR(2,1), KR(3,1), KR(4,1), KR(5,1), &
                        KR(6,1)/1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0/
      DATA              KR(1,2), KR(2,2), KR(3,2), KR(4,2), KR(5,2), &
                        KR(6,2)/780.D0, 1850.D0, 3822.D0, 6103.D0, &
                        15152.D0, 30954.D0/
      DATA              KR(1,3), KR(2,3), KR(3,3), KR(4,3), KR(5,3), &
                        KR(6,3)/359.D0, 1476.D0, 544.D0, 4104.D0, &
                        16592.D0, 19394.D0/
      DATA              KR(1,4), KR(2,4), KR(3,4), KR(4,4), KR(5,4), &
                        KR(6,4)/766.D0, 792.D0, 1206.D0, 6016.D0, &
                        9023.D0, 15710.D0/
      DATA              KR(1,5), KR(2,5), KR(3,5), KR(4,5), KR(5,5), &
                        KR(6,5)/618.D0, 840.D0, 198.D0, 6019.D0, &
                        12216.D0, 2302.D0/
      DATA              KR(1,6), KR(2,6), KR(3,6), KR(4,6), KR(5,6), &
                        KR(6,6)/41.D0, 2037.D0, 2240.D0, 4167.D0, &
                        4902.D0, 9227.D0/
      DATA              KR(1,7), KR(2,7), KR(3,7), KR(4,7), KR(5,7), &
                        KR(6,7)/596.D0, 229.D0, 2304.D0, 3851.D0, &
                        12506.D0, 3420.D0/
      DATA              KR(1,8), KR(2,8), KR(3,8), KR(4,8), KR(5,8), &
                        KR(6,8)/86.D0, 1578.D0, 436.D0, 4138.D0, &
                        7824.D0, 3824.D0/
      DATA              KR(1,9), KR(2,9), KR(3,9), KR(4,9), KR(5,9), &
                        KR(6,9)/636.D0, 526.D0, 470.D0, 259.D0, 6093.D0, &
                        22300.D0/
      DATA              KR(1,10), KR(2,10), KR(3,10), KR(4,10), &
                        KR(5,10), KR(6,10)/287.D0, 431.D0, 1554.D0, &
                        1117.D0, 12088.D0, 5130.D0/
      DATA              KR(1,11), KR(2,11), KR(3,11), KR(4,11), &
                        KR(5,11), KR(6,11)/707.D0, 1485.D0, 480.D0, &
                        1188.D0, 2399.D0, 11222.D0/
      DATA              KR(1,12), KR(2,12), KR(3,12), KR(4,12), &
                        KR(5,12), KR(6,12)/707.D0, 1450.D0, 1004.D0, &
                        173.D0, 8764.D0, 17698.D0/
      DATA              KR(1,13), KR(2,13), KR(3,13), KR(4,13), &
                        KR(5,13), KR(6,13)/96.D0, 1001.D0, 684.D0, &
                        2919.D0, 5491.D0, 7057.D0/
      DATA              KR(1,14), KR(2,14), KR(3,14), KR(4,14), &
                        KR(5,14), KR(6,14)/49.D0, 1001.D0, 684.D0, &
                        235.D0, 9274.D0, 28739.D0/
      DATA              KR(1,15), KR(2,15), KR(3,15), KR(4,15), &
                        KR(5,15), KR(6,15)/373.D0, 1001.D0, 1447.D0, &
                        3043.D0, 3054.D0, 33207.D0/
      DATA              KR(1,16), KR(2,16), KR(3,16), KR(4,16), &
                        KR(5,16), KR(6,16)/613.D0, 2.D0, 857.D0, &
                        1249.D0, 2648.D0, 27717.D0/
      DATA              KR(1,17), KR(2,17), KR(3,17), KR(4,17), &
                        KR(5,17), KR(6,17)/373.D0, 2.D0, 2.D0, 1249.D0, &
                        2648.D0, 33207.D0/
      DATA              KR(1,18), KR(2,18), KR(3,18), KR(4,18), &
                        KR(5,18), KR(6,18)/2.D0, 2.D0, 2.D0, 2.D0, &
                        2648.D0, 1420.D0/
      DATA              KR(1,19), KR(2,19), KR(3,19), KR(4,19), &
                        KR(5,19), KR(6,19)/2.D0, 2.D0, 2.D0, 2.D0, 2.D0, &
                        1420.D0/
      DATA              KR(1,20), KR(2,20), KR(3,20), KR(4,20), &
                        KR(5,20), KR(6,20)/2.D0, 2.D0, 2.D0, 2.D0, 2.D0, &
                        2.D0/
      DATA              KPTS(1), KPTS(2), KPTS(3), KPTS(4), KPTS(5), &
                        KPTS(6)/2129.D0, 5003.D0, 10007.D0, 20011.D0, &
                        40009.D0, 80021.D0/
!C     .. Executable Statements ..
      NDIM = N
      PTS = DBLE(NPTS)
!C     VALIDITY CHECK
      ICHECK = 1
!!!!      IF (NDIM.LT.1 .OR. NDIM.GT.MAXDIM) GO TO 160
      ICHECK = 2
!!!!      IF (NPTS.LT.1) GO TO 160
      ICHECK = 3
!!!!      IF (NRAND.LT.1) GO TO 160
      ICHECK = 0
      IF (NPTS.GT.MAXRUL) GO TO 40
!C     SELECT KOROBOV VECTOR
      NP = NPTS
      PTS = KPTS(NPTS)
      VK(1) = 1.0D0
      IF (NDIM.EQ.1) GO TO 40
      VK(2) = KR(NP,NDIM)
      IF (NDIM.EQ.2) GO TO 40
      XA = VK(2)
      DO 20 I = 3, NDIM
         VK(I) = AINT(MOD(VK(I-1)*XA,PTS)+0.5D0)
!C        ************* IMPLEMENTATION-DEPENDENT CODING ****************
!C        IF THE MACHINE DOES NOT PROVIDE AT LEAST 32 BITS OF BASIC
!C        PRECISION, THEN INTERMEDIATE WORKING IN THE PRECEDING
!C        STATEMENT MUST BE CHANGED TO ADDITIONAL PRECISION, E.G. IN A
!C        SINGLE PRECISION IMPLEMENTATION CHANGED TO
!C        VK(I) = AINT(SNGL(DMOD(DBLE(VK(I-1))*DBLE(XA),DBLE(PTS)))+
!C        *    0.5)
!C        THE FUNCTION DMOD MUST BE SUPPLIED IF NECESSARY.
!C        **************************************************************
   20 CONTINUE
!C     BEGIN INTEGRATION
   40 SUM1 = 0.0D0
      SUM2 = 0.0D0
      DO 140 M = 1, NRAND
         RES = 0.0D0
!C        CALCULATE RANDOM SHIFT
         DO 60 K = 1, NDIM
            ALPHA(K) = rand(0)!!!!G05CAF(ALPHA(K))
   60    CONTINUE
!C        CALCULATE TRANSFORMED INTEGRAND
         PK = 1.0D0
         DPK = 1.0D0/PTS
 !SG170328
   !$omp parallel do &
   default(shared) &
   private(wt,xx,c,d,grad,j,x) &
   reduction(+:res)
do k1=1,int(pts)                                        !SG-PF
 !SG170328 if(omp_get_thread_num().eq.0) &
 !SG170328		print *,OMP_GET_NUM_THREADS()
! print *,omp_get_thread_num()
   80    WT = DPK
         DO 120 J = 1, NDIM
!            XX = MOD(ALPHA(J)+VK(J)*PK*DPK,1.0D0)
            XX = MOD(ALPHA(J)+VK(J)*dble(k1)*DPK,1.0D0) !SG-PF
            CALL REGION(NDIM,X,J,C,D)
            GRAD = D - C
            WT = WT*GRAD
            IF (ITRANS.NE.0) GO TO 100
!C           PERIODIZING TRANSFORMATION
            WT = WT*6.0D0*XX*(1.0D0-XX)
            XX = XX*XX*(3.0D0-2.0D0*XX)
  100       X(J) = C + GRAD*XX
  120    CONTINUE
         RES = F(NDIM,X,obj)*WT + RES
!!!!         PK = AINT(PK+1.5D0)
enddo                                                   !SG-PF
  !$omp end parallel do
  !SG170328
!!!!         IF (PK.LE.PTS) GO TO 80
         SUM1 = SUM1 + RES
         SUM2 = SUM2 + RES*RES
  140 CONTINUE
!C     MEAN
      RES = SUM1/DBLE(NRAND)
      ERRr = 0.0D0
!!!!      IF (NRAND.EQ.1) GO TO 160
!C     STANDARD ERROR
      ERRr = SQRT(ABS((SUM2-DBLE(NRAND)*RES*RES)/DBLE(NRAND*(NRAND-1))))
!!!!  160 IFAIL = P01ABF(IFAIL,ICHECK,SRNAME,0,P01REC)
      RETURN
      END subroutine D01GCF
      
      
	!taken from CAMB and adapted:
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function rombint_re(obj,f,a,b,tol, maxit)
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
			real(dl) obj
			real(dl) f
			external f
			real(dl) :: rombint_re
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
	40      rombint_re=g0
			if (i.gt.MAXITER.and.abs(error).gt.tol)  then
			  write(*,*) 'Warning: Rombint failed to converge; '
			  write (*,*)'integral, error, tol:', rombint_re,error, tol
			end if
	end function rombint_re

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function rombint_vec(obj,f,a,b,tol, maxit)
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
			real(dl), dimension(:), intent(in) :: obj
			real(dl) f
			external f
			real(dl) :: rombint_vec
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
	40      rombint_vec=g0
			if (i.gt.MAXITER.and.abs(error).gt.tol)  then
			  write(*,*) 'Warning: Rombint failed to converge; '
			  write (*,*)'integral, error, tol:', rombint_vec,error, tol
			end if
			
	end function rombint_vec

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function rombint_nD(obj,f,a,b,tol, maxit)
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
			type(nuDensArgs) :: obj
			real(dl) f
			external f
			real(dl) :: rombint_nD
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
	40      rombint_nD=g0
			if (i.gt.MAXITER.and.abs(error).gt.tol)  then
			  write(*,*) 'Warning: Rombint failed to converge; '
			  write (*,*)'integral, error, tol:', rombint_nD,error, tol
			end if
			
	end function rombint_nD

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function rombint_dRN(obj,f,a,b,tol, maxit)
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
			type(integrRhoNuPar) :: obj
			real(dl) f
			external f
			real(dl) :: rombint_dRN
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
	40      rombint_dRN=g0
			if (i.gt.MAXITER.and.abs(error).gt.tol)  then
			  write(*,*) 'Warning: Rombint failed to converge; '
			  write (*,*)'integral, error, tol:', rombint_dRN,error, tol
			end if
			
	end function rombint_dRN

end module utilities
