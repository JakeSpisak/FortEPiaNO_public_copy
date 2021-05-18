module fpStuff
	use precision
	use constants
	use utilities
	use fpErrors
	use fpInteractions
	use fpEquations
	use fpCosmology
	use linear_interpolation_module

	type nuDensArgs
		real(dl) :: x,z
		integer iFl
	end type nuDensArgs

	type interpNuDens_obj
		type(spline_class), dimension(:,:), allocatable :: re, im
	end type interpNuDens_obj
	type(interpNuDens_obj) :: interpNuDens

	type(spline_class) :: FD_interp
	type(linear_interp_4d) :: D2_interp, D3_interp, PI1_12_interp, PI1_13_interp

	contains

	subroutine deallocateStuff()
		integer :: nf
		nf = flavorNumber
		deallocate(nuMasses, nuFactor, sterile, Gs)
		deallocate(mixMat, mixMatInv)
		deallocate(nuMassesMat, leptonDensities)
		deallocate(dampTermMatrixCoeffNue, dampTermMatrixCoeffNunu)
		deallocate(dampingSetZero)
		deallocate(GL_mat, GR_mat, GLR_vec)
		deallocate(mixingAngles, massSplittings)
	end subroutine deallocateStuff

	pure subroutine deallocateCmplxMat(m)
		type(cmplxMatNN), intent(inout) :: m

		if (m%a) then
			m%a=.false.
			deallocate(m%re, m%im)
		end if
	end subroutine deallocateCmplxMat

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
   !$omp parallel do default(shared) private(wt,xx,c,d,grad,j,x) reduction(+:res)
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

	subroutine allocate_interpNuDens
		real(dl), dimension(:), allocatable :: ndmv_re, ndmv_im
		integer :: i, j, iy
		allocate(interpNuDens%re(flavorNumber, flavorNumber), &
			interpNuDens%im(flavorNumber, flavorNumber))
		allocate(ndmv_re(Ny), ndmv_im(Ny))
		do i=1, flavorNumber
			do iy=1, Ny
				ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, i)
			end do
			call interpNuDens%re(i, i)%initialize(Ny, y_arr, ndmv_re)
			do j=i+1, flavorNumber
				do iy=1, Ny
					ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, j)
					ndmv_im(iy) = nuDensMatVecFD(iy)%im(i, j)
				end do
				call interpNuDens%re(i, j)%initialize(Ny, y_arr, ndmv_re)
				call interpNuDens%im(i, j)%initialize(Ny, y_arr, ndmv_im)
			end do
		end do
	end subroutine allocate_interpNuDens

	subroutine init_interpNuDens
		real(dl), dimension(:), allocatable :: ndmv_re, ndmv_im
		integer :: i, j, iy
		allocate(ndmv_re(Ny), ndmv_im(Ny))
		do i=1, flavorNumber
			do iy=1, Ny
				ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, i)
			end do
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
			do j=i+1, flavorNumber
				do iy=1, Ny
					ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, j)
					ndmv_im(iy) = nuDensMatVecFD(iy)%im(i, j)
				end do
				call interpNuDens%re(i, j)%replace(Ny, y_arr, ndmv_re)
				call interpNuDens%im(i, j)%replace(Ny, y_arr, ndmv_im)
			end do
		end do
	end subroutine init_interpNuDens

	pure function fermiDirac_i(x)
		real(dl) :: fermiDirac_i
		real(dl), intent(in) :: x

		fermiDirac_i = FD_interp%evaluate(x)
	end function fermiDirac_i

	subroutine init_interp_FD
		real(dl), dimension(:), allocatable :: fd_vec, fd_x
		integer :: ix, iflag, interp_nfd
		real(dl) :: x, t1,t2
		real(8) :: timer1

		call addToLog("[interactions] Initializing interpolation for FermiDirac...")
		interp_nfd=100
		allocate(fd_x(interp_nfd), fd_vec(interp_nfd))
		fd_x=logspace(-3.d0, 2.d0, interp_nfd)
		do ix=1, interp_nfd
			fd_vec(ix) = fermiDirac(fd_x(ix))
		end do
		call FD_interp%initialize(interp_nfd, fd_x, fd_vec)

		call random_seed()
		call random_number(x)
		x=10.d0**(x*(2.d0-(-3.d0)) - 3.d0)
		write(*,"(' [interactions] test fermiDirac in ',*(E12.5))") x
		t1 = fermiDirac(x)
		t2 = fermiDirac_i(x)
		write(*,"(' [interactions] comparison (true vs interp): ',*(E17.10))") t1,t2

		if (timing_tests) then
			write (*,*) "[interactions] now doing some timing..."
			call tic(timer1)
			do ix=1, 100000000
				call random_number(x)
				x=10.d0**(x*(2.d0-(-3.d0)) - 3.d0)
				t1 = fermiDirac_i(x)
			end do
			call toc(timer1, "<reset>")

			call tic(timer1)
			do ix=1, 100000000
				call random_number(x)
				x=10.d0**(x*(2.d0-(-3.d0)) - 3.d0)
				t1 = fermiDirac_i(x)
			end do
			call toc(timer1, "<interpolated>")

			call tic(timer1)
			do ix=1, 100000000
				call random_number(x)
				x=10.d0**(x*(2.d0-(-3.d0)) - 3.d0)
				t1 = fermiDirac(x)
			end do
			call toc(timer1, "<full>")
		end if

		call addToLog("[interactions] ...done!")
	end subroutine init_interp_FD

	subroutine time_electron_energyDensity
		real(dl), dimension(:,:), allocatable :: ed_vec, md_vec
		integer :: ix, iz, iflag
		real(dl) :: x,z, t1,t2
		real(8) :: timer1

		call random_seed()

		if (timing_tests) then
			call tic(timer1)
			write (*,*) "[cosmo] now doing some timing for electron energy density..."
			call tic(timer1)
			do ix=1, 1000000
				call random_number(x)
				call random_number(z)
				x=(x_fin-x_in)*x + x_in
				z=0.4d0*z + z_in
				t1 = electrons%energyDensity(x,z, .false.)
			end do
			call toc(timer1, "<interpolated>")

			call tic(timer1)
			do ix=1, 1000000
				call random_number(x)
				call random_number(z)
				x=(x_fin-x_in)*x + x_in
				z=0.4d0*z + z_in
				t1 = electrons%energyDensityFull(x,z, .false.)
			end do
			call toc(timer1, "<full>")
		end if

		call addToLog("[cosmo] ...done!")
	end subroutine time_electron_energyDensity

	subroutine test_speed_coll_int
		real(dl) :: x,z, y1
		type(coll_args) :: collArgs
		integer :: i, ix, Npt!,j, k
		real(dl) :: errr1,errr2, res1,res2,res3,res4, cf
		INTEGER :: IFAIL, ITRANS, N, NPTS, NRAND
		real(dl) ::  VK(2)
		real(dl), dimension(:), allocatable :: ndmv_re
		real(8) :: timer1

		x=0.05d0
		y1=1.2d0
		z=1.06d0
		npts=1
		nrand=1
		n=2
		Npt=1000!number of calculations for each comparison

		allocate(ndmv_re(Ny))

		do i=1, flavorNumber
			do ix=1, Ny
				nuDensMatVecFD(ix)%re(i, i) = 1.d0*i * fermiDirac(y_arr(ix)/z)
				ndmv_re(ix) = nuDensMatVecFD(ix)%re(i, i)
			end do
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
		end do

		collArgs%x = x
		collArgs%z = z
		collArgs%y1 = y1
		collArgs%dme2 = 0.d0!dme2_electron(x, 0.d0, z)
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		collArgs%iy = 12

		write (*,*) "[interactions] timing 2D integrals..."
		call tic(timer1)
		do ix=1, Npt
			call random_number(x)
			call random_number(z)
			collArgs%x = (x_fin-x_in)*x + x_in
			collArgs%z = 0.4d0*z + z_in
			ifail=0
			itrans=0
			res2 = integrate_collint_nue_NC(coll_nue_sc_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
		end do
		call toc(timer1, "<sc re semilin>")
		call tic(timer1)
		do ix=1, Npt
			call random_number(x)
			call random_number(z)
			collArgs%x = (x_fin-x_in)*x + x_in
			collArgs%z = 0.4d0*z + z_in
			ifail=0
			itrans=0
			call D01GCF(n,coll_nue_4_sc_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
		end do
		call toc(timer1, "<sc re D01GCF>")

		call tic(timer1)
		do ix=1, Npt
			call random_number(x)
			call random_number(z)
			collArgs%x = (x_fin-x_in)*x + x_in
			collArgs%z = 0.4d0*z + z_in
			ifail=0
			itrans=0
			res2 = integrate_collint_nue_NC(coll_nue_ann_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
		end do
		call toc(timer1, "<ann re semilin>")
		call tic(timer1)
		do ix=1, Npt
			call random_number(x)
			call random_number(z)
			collArgs%x = (x_fin-x_in)*x + x_in
			collArgs%z = 0.4d0*z + z_in
			ifail=0
			itrans=0
			call D01GCF(n,coll_nue_4_ann_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
		end do
		call toc(timer1, "<ann re D01GCF>")
		
		call tic(timer1)
		do ix=1, 10*Npt
			call random_number(x)
			call random_number(z)
			collArgs%x = (x_fin-x_in)*x + x_in
			collArgs%z = 0.4d0*z + z_in
			ifail=0
			itrans=0
			res1 = integrate_collint_nue_NC(coll_nue_ann_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
			res2 = integrate_collint_nue_NC(coll_nue_sc_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
			res2 = integrate_collint_nue_NC(coll_nue_int, collArgs, F_ab_ann_re, F_ab_sc_re)
		end do
		call toc(timer1, "<reset>")
		call tic(timer1)
		do ix=1, 10*Npt
			call random_number(x)
			call random_number(z)
			collArgs%x = (x_fin-x_in)*x + x_in
			collArgs%z = 0.4d0*z + z_in
			ifail=0
			itrans=0
			res1 = integrate_collint_nue_NC(coll_nue_ann_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
			res2 = integrate_collint_nue_NC(coll_nue_sc_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
		end do
		call toc(timer1, "<sep ann, sc>")
		call tic(timer1)
		do ix=1, 10*Npt
			call random_number(x)
			call random_number(z)
			collArgs%x = (x_fin-x_in)*x + x_in
			collArgs%z = 0.4d0*z + z_in
			ifail=0
			itrans=0
			res2 = integrate_collint_nue_NC(coll_nue_int, collArgs, F_ab_ann_re, F_ab_sc_re)
		end do
		call toc(timer1, "<sum ann+sc>")
		
		deallocate(ndmv_re)
		call addToLog("[interactions] ...done!")
	end subroutine test_speed_coll_int

	function integr_rho_nu(vec,y)
		real(dl) :: integr_rho_nu
		real(dl), intent(in) :: y
		type(nuDensArgs), intent(in) :: vec
		integr_rho_nu = y*y*y * interpNuDens%re(vec%iFl, vec%iFl)%evaluate(y) * fermiDirac(y)
	end function integr_rho_nu

	function nuDensity(z, iFl)
		real(dl) :: nuDensity
		real(dl), intent(in) :: z
		integer, intent(in) :: iFl
		type(nuDensArgs) :: vec

		vec%x=0.d0
		vec%z=z
		vec%iFl=iFl
		nuDensity = rombint_nD(vec, integr_rho_nu, y_min, y_max, 1d-3, maxiter) / PISQ
	end function nuDensity

	pure function integrate_dRhoNu(params, y)
		real(dl) :: integrate_dRhoNu
		real(dl), intent(in) :: y
		type(spline_class), intent(in) :: params

		integrate_dRhoNu = params%evaluate(y)
	end function integrate_dRhoNu

#ifndef NO_INTERPOLATION
	subroutine dz_o_dx_old(x,z, ydot, n)!eq 17 from doi:10.1016/S0370-2693(02)01622-2
		real(dl), intent(in) :: x,z
		integer, intent(in) :: n
		real(dl), dimension(n), intent(inout) :: ydot
		real(dl) :: dzodx, tmp
		real(dl), dimension(2) :: coeffs
		type(spline_class) :: params
		real(dl), dimension(:), allocatable :: der
		integer :: ix, m

		coeffs = dzodxcoef_interp_func(x/z)

		allocate(der(Ny))

		der = 0.d0
		do m=1, Ny
			do ix=1, flavorNumber
				der(m) = der(m) + ydot((m-1)*flavNumSqu + ix) * nuFactor(ix)
			end do
			der(m) = y_arr(m)**3 * der(m) * fermiDirac(y_arr(m))
		end do
		call params%initialize(Ny, y_arr, der)
		tmp = rombint_spli(params, integrate_dRhoNu, y_min, y_max, 1.d-3, maxiter)

		dzodx = coeffs(1) - coeffs(2)/ z**3 * tmp

		ydot(n)=dzodx
	end subroutine dz_o_dx_old

	subroutine test_dzodx_speed
		real(dl), dimension(:), allocatable :: fd_vec, fd_x
		integer :: ix, iflag, interp_nfd
		real(dl) :: w,x,z
		real(8) :: timer1
		integer, parameter :: n=901
		real(dl), dimension(n) :: ydot
		integer :: m

		ydot = 0.d0
		do m=1, Ny
			ydot((m-1)*flavNumSqu + 1) = cos(0.02d0*y_arr(m))
			ydot((m-1)*flavNumSqu + 2) = y_arr(m)/20.d0
			ydot((m-1)*flavNumSqu + 3) = 1.d0
		end do

		call addToLog("[equations] Testing dz_o_dx speed...")
		call random_seed()
		call random_number(x)
		call random_number(z)
		call random_number(w)
		call tic(timer1)
		do ix=1, 1000000
			call random_number(x)
			call random_number(z)
			call random_number(w)
			x=(x_fin-x_in)*x + x_in
			z=0.4d0*z + z_in
			w=0.4d0*w + z_in
			call dz_o_dx_old(x, z, ydot, n)
			call dz_o_dx(x, w, z, ydot, n)
		end do
		call toc(timer1, "<reset>")

		call tic(timer1)
		do ix=1, 1000000
			call random_number(x)
			call random_number(z)
			x=(x_fin-x_in)*x + x_in
			z=0.4d0*z + z_in
			call dz_o_dx_old(x, z, ydot, n)
		end do
		call toc(timer1, "<rombint>")

		call tic(timer1)
		do ix=1, 1000000
			call random_number(x)
			call random_number(z)
			call random_number(w)
			x=(x_fin-x_in)*x + x_in
			z=0.4d0*z + z_in
			w=0.4d0*w + z_in
			call dz_o_dx(x, w, z, ydot, n)
		end do
		call toc(timer1, "<linearized>")

		call addToLog("[equations] ...done!")
	end subroutine test_dzodx_speed
#endif

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function rombint_vec(obj,f,a,b,tol, maxit)
		use Precision
		!  Rombint returns the integral from a to b of using Romberg integration.
		!  The method converges provided that f(x) is continuous in (a,b).
		!  f must be real(dl) and must be declared external in the calling
		!  routine.  tol indicates the desired relative accuracy in the integral.
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

		if (present(maxit)) then
			MaxIter = maxit
		end if
		h=0.5d0*(b-a)
		gmax=h*(f(obj,a)+f(obj,b))
		g(1)=gmax
		nint=1
		error=1.0d20
		i=0
10      i=i+1
		if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
			go to 40
		!  Calculate next trapezoidal rule approximation to integral.
		g0=0._dl
		do 20 k=1,nint
			g0=g0+f(obj,a+(k+k-1)*h)
20      continue
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
30      continue
		if (abs(g0).gt.tol) then
			error=1._dl-gmax/g0
		else
			error=gmax
		end if
		gmax=g0
		g(jmax+1)=g0
		go to 10
40      rombint_vec=g0
	end function rombint_vec

	function rombint_nD(obj,f,a,b,tol, maxit)
			use Precision
	!  Rombint returns the integral from a to b of using Romberg integration.
	!  The method converges provided that f(x) is continuous in (a,b).
	!  f must be real(dl) and must be declared external in the calling
	!  routine.  tol indicates the desired relative accuracy in the integral.

			implicit none
			integer, intent(in), optional :: maxit
			integer :: MAXITER
			integer, parameter :: MAXJ=5
			dimension g(MAXJ+1)
			real(dl), intent(in) :: a,b,tol
			type(nuDensArgs), intent(in) :: obj
			real(dl) f
			external f
			real(dl) :: rombint_nD
			integer :: nint, i, k, jmax, j
			real(dl) :: h, gmax, error, g, g0, g1, fourj

			if (present(maxit)) then
				MaxIter = maxit
			else
				MaxIter = 20
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
	pure function rombint_spli(obj,f,a,b,tol, maxit)
			use Precision
	!  Rombint returns the integral from a to b of using Romberg integration.
	!  The method converges provided that f(x) is continuous in (a,b).
	!  f must be real(dl) and must be declared external in the calling
	!  routine.  tol indicates the desired relative accuracy in the integral.
	!
		implicit none

		interface
			pure real(dl) function f(a, b)
				use precision
				use fpInterpolate
				type(spline_class), intent(in) :: a
				real(KIND(1.d0)), intent(in) :: b
			end function
		end interface

		integer, intent(in), optional :: maxit
		type(spline_class), intent(in) :: obj
		real(dl), intent(in) :: a,b,tol
		real(dl) :: rombint_spli
		integer :: MAXITER
		integer, parameter :: MAXJ=5
		dimension g(MAXJ+1)
		integer :: nint, i, k, jmax, j
		real(dl) :: h, gmax, error, g, g0, g1, fourj

			if (present(maxit)) then
				MaxIter = maxit
			else
				MAXITER=20
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
	40      rombint_spli=g0
!			if (i.gt.MAXITER.and.abs(error).gt.tol)  then
!			  write(*,*) 'Warning: Rombint failed to converge; '
!			  write (*,*)'integral, error, tol:', rombint_spli,error, tol
!			end if
	end function rombint_spli

	subroutine test_nuDens_speed
		real(dl), dimension(:), allocatable :: fd_vec, fd_x
		integer :: ix, iflag, interp_nfd
		real(dl) :: x, t1,t2
		real(8) :: timer1

		call addToLog("[cosmology] Testing nuDensity speed...")
		call random_seed()
		call random_number(x)
		call tic(timer1)
		do ix=1, 1000000
			call random_number(x)
			x=10.d0**(x/6.d0)
			t1 = nuDensity(x, 1)
			t1 = nuDensityNC(1, 1)
		end do
		call toc(timer1, "<reset>")

		call tic(timer1)
		do ix=1, 1000000
			call random_number(x)
			x=10.d0**(x/6.d0)
			t1 = nuDensity(x, 1)
		end do
		call toc(timer1, "<rombint>")

		call tic(timer1)
		do ix=1, 1000000
			call random_number(x)
			x=10.d0**(x/6.d0)
			t1 = nuDensityNC(1, 1)
		end do
		call toc(timer1, "<linearized>")

		call addToLog("[cosmology] ...done!")
	end subroutine test_nuDens_speed

	!kappa function for damping terms (Bennett:2020zkv)
	elemental function kappa_damp(a, b, c)
		real(dl) :: kappa_damp
		real(dl), intent(in) :: a, b, c
		kappa_damp = a*a*(c-b) - 2.d0/3.d0*a* (c**3-b**3) + (c**5-b**5)/5.d0
	end function kappa_damp

	!function that returns the integrand of the nunu damping coefficient in Bennett:2020zkv terms using K kernels
	elemental function nunu_damp_integrand(p, p1, q)
		real(dl) :: nunu_damp_integrand
		real(dl), intent(in) :: p, q, p1
		real(dl) :: q1
		real(dl) :: as, bs, cs, ap, bp, cp
		real(dl) :: fq, fp1, fq1
		q1 = p + q - p1

		cs = p + q
		as = cs**2
		bs = max(abs(p-q), abs(p1-q1))

		ap = (q - p1)**2
		bp = abs(p1 - q)
		cp = min(p+q1, q+p1)

		fq = fermiDirac(q)
		fp1 = fermiDirac(p1)
		fq1 = fermiDirac(q1)

		nunu_damp_integrand = &
			(kappa_damp(as, bs, cs) + 2.d0 * kappa_damp(ap, bp, cp)) &
			* ((1.d0-fq)*fp1*fq1 + fq*(1.d0-fp1)*(1.d0-fq1))
	end function nunu_damp_integrand

	!integral of the above function
	pure function dy_damping(y)
		real(dl) :: dy_damping
		real(dl), intent(in) :: y
		integer :: ia, ib
		integer, parameter :: nydamp = 2**10
		real(dl), dimension(:), allocatable :: ya, dya
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(ya(nydamp), dya(nydamp))
		ya = linspace(0.d0, 100.d0, nydamp)
		do ia=1, nydamp-1
			dya(ia) = ya(ia+1) - ya(ia)
		end do
		allocate(fy2_arr(nydamp, nydamp))
		fy2_arr = 0.d0
		do ia=1, nydamp
			do ib=1, nydamp
				if (ya(ib) .gt. ya(ia)-y) &
					fy2_arr(ia, ib) = nunu_damp_integrand(y, ya(ia), ya(ib))
			end do
		end do
		dy_damping = integral_NC_2d(nydamp, nydamp, dya, dya, fy2_arr) / y**3
		deallocate(ya, dya, fy2_arr)
	end function dy_damping

	!function that returns the integrand of the nunu damping coefficient in Bennett:2020zkv terms using PI kernels
	elemental function nunu_damp_integrand_pi(p, p1, q)
		real(dl) :: nunu_damp_integrand_pi
		real(dl), intent(in) :: p, q, p1
		real(dl), dimension(2) :: pi2
		real(dl) :: q1, fq, fp1, fq1

		nunu_damp_integrand_pi = 0.d0
		q1 = p + q - p1
		if ( &
			.not.(q1.lt.0.d0 &
			.or. p .gt.p1+q +q1 &
			.or. q .gt.p +p1+q1 &
			.or. p1.gt.p +q +q1 &
			.or. q1.gt.p +q +p1 &
			) &
		) then
			pi2 = PI2_ne_f(p, q, p1, q1, q, q1)
			fq = fermiDirac(q)
			fp1 = fermiDirac(p1)
			fq1 = fermiDirac(q1)
			nunu_damp_integrand_pi = &
				(pi2(2) + 2.d0 * pi2(1)) &
				* ((1.d0-fq)*fp1*fq1 + fq*(1.d0-fp1)*(1.d0-fq1))
		end if
	end function nunu_damp_integrand_pi

	!integral of the above function
	pure function dy_damping_pi(y)
		real(dl) :: dy_damping_pi
		real(dl), intent(in) :: y
		integer :: ia, ib
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(fy2_arr(Ny, Ny))
		fy2_arr = 0.d0
		do ia=1, Ny
			do ib=1, Ny
				fy2_arr(ia, ib) = nunu_damp_integrand_pi(y, y_arr(ia), y_arr(ib))
			end do
		end do
		dy_damping_pi = integral_GL_2d(Ny, w_gl_arr2, w_gl_arr2, fy2_arr) / y**3
		deallocate(fy2_arr)
	end function dy_damping_pi

	pure function F_nu_sc_da(n1, n2, n3, n4, i, j)
		real(dl) :: F_nu_sc_da
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		integer :: k
		real(dl), dimension(:), allocatable :: t2r, t4r
		real(dl) :: s24r

		F_nu_sc_da = 0.d0

		if (i.ne.j) &
			return
		allocate(t2r(flavorNumber), t4r(flavorNumber))
		do k=1, flavorNumber
			s24r = n2%re(k, k)*n4%re(k, k)
			t2r(k) = n4%re(k, k) - s24r
			t4r(k) = n2%re(k, k) - s24r
		end do
		F_nu_sc_da = &
			n3%re(i,i) * (t2r(i) + sum(t2r)) &
			+ (1.d0 - n3%re(i,i)) * (t4r(i) + sum(t4r))
		deallocate(t2r, t4r)
		F_nu_sc_da = 2.d0 * F_nu_sc_da ! +h.c.
	end function F_nu_sc_da

	pure function F_nu_pa_da(n1, n2, n3, n4, i, j)
		real(dl) :: F_nu_pa_da
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		integer :: k
		real(dl), dimension(:), allocatable :: t2r, t4r
		real(dl) :: sBr

		F_nu_pa_da = 0.d0

		if (i.ne.j) &
			return

		allocate(t2r(flavorNumber), t4r(flavorNumber))
		!first line: sAr -> 12, sBr->34
		do k=1, flavorNumber
			sBr = n4%re(k, k)*n3%re(k, k)
			t2r(k) = sBr
			t4r(k) = sBr + 1.d0 - n4%re(k, k) - n3%re(k, k)
		end do
		F_nu_pa_da = F_nu_pa_da &
			+ (1.d0 - n2%re(i,i)) * (t2r(i) + sum(t2r)) &
			+ n2%re(i,i) * (t4r(i) + sum(t4r))
		!second line: sAr -> 13, sBr->24
		do k=1, flavorNumber
			sBr = n4%re(k, k)*n2%re(k, k)
			t2r(k) = n4%re(k, k) - sBr
			t4r(k) = n2%re(k, k) - sBr
		end do
		F_nu_pa_da = F_nu_pa_da &
			+ n3%re(i,i) * (t2r(i) + sum(t2r)) &
			+ (1.d0 - n3%re(i,i)) * (t4r(i) + sum(t4r))
		deallocate(t2r, t4r)
		F_nu_pa_da = 2.d0 * F_nu_pa_da ! +h.c.
	end function F_nu_pa_da


	subroutine init_interp_d123
		real(dl), dimension(:,:,:,:), allocatable :: d2, d3, pi1_12, pi1_13
		real(dl) :: y1,  y2,  y3,  y4, d1a, d1b, d2a, d2b, d3a, d3b
		integer  :: ix1, ix2, ix3, ix4, iflag
		real(dl), dimension(2) :: PI1
		real(dl), dimension(:), allocatable :: yarr
		real(8) :: timer1

		call addToLog("[interactions] Initializing interpolation for D1/2/3 functions...")
		allocate(yarr(interp_ny))
		yarr = logspace(log10(y_min), log10(y_max), interp_ny)
		allocate(d2(interp_ny,interp_ny,interp_ny,interp_ny), &
			d3(interp_ny,interp_ny,interp_ny,interp_ny), &
			pi1_12(interp_ny,interp_ny,interp_ny,interp_ny), &
			pi1_13(interp_ny,interp_ny,interp_ny,interp_ny))
		do ix1=1, interp_ny
			do ix2=1, interp_ny
				do ix3=1, interp_ny
					do ix4=1, interp_ny
						d2(ix1,ix2,ix3,ix4) = d2_full(yarr(ix1), yarr(ix2), yarr(ix3), yarr(ix4))
						d3(ix1,ix2,ix3,ix4) = d3_full(yarr(ix1), yarr(ix2), yarr(ix3), yarr(ix4))

						pi1_12(ix1,ix2,ix3,ix4) = PI1_12_full(yarr(ix1), yarr(ix2), yarr(ix3), yarr(ix4))
						pi1_13(ix1,ix2,ix3,ix4) = PI1_13_full(yarr(ix1), yarr(ix2), yarr(ix3), yarr(ix4))
					end do
				end do
			end do
		end do
		call PI1_12_interp%initialize(yarr, yarr, yarr, yarr, pi1_12, iflag)
		deallocate(pi1_12)
		call PI1_13_interp%initialize(yarr, yarr, yarr, yarr, pi1_13, iflag)
		deallocate(pi1_13)

		call D2_interp%initialize(yarr, yarr, yarr, yarr, d2, iflag)
		deallocate(d2)
		call D3_interp%initialize(yarr, yarr, yarr, yarr, d3, iflag)
		deallocate(d3)

		call random_seed()
		call random_number(y1)
		call random_number(y2)
		call random_number(y3)
		call random_number(y4)
		y1=(y_max-y_min)*y1 + y_min
		y2=(y_max-y_min)*y2 + y_min
		y3=(y_max-y_min)*y3 + y_min
		y4=(y_max-y_min)*y4 + y_min

		write(*,"(' [interactions] test PI_1, D2_f, D3_f in ',*(E12.5))") y1,y2,y3,y4
		d2a = D2_f(y1,y2,y3,y4)
		d3a = D3_f(y1,y2,y3,y4)
		d2b = D2_full(y1,y2,y3,y4)
		d3b = D3_full(y1,y2,y3,y4)
		pi1 = PI1_full(y1,y2,y3,y4)
		write(*,"(' [interactions] comparison pi1_12_i (true vs interp): ',*(E17.10))") PI1_12_full(y1,y2,y3,y4), pi1_12_i(y1,y2,y3,y4)
		write(*,"(' [interactions] comparison pi1_13_i (true vs interp): ',*(E17.10))") PI1_13_full(y1,y2,y3,y4), pi1_13_i(y1,y2,y3,y4)
		write(*,"(' [interactions] comparison D2_f (true vs interp): ',*(E17.10))") d2a,d2b
		write(*,"(' [interactions] comparison D3_f (true vs interp): ',*(E17.10))") d3a,d3b

		if (timing_tests) then
			call tic(timer1)
			write (*,*) "[interactions] now doing some timing..."
			call tic(timer1)
			do ix1=1, 10000000
				call random_number(y1)
				call random_number(y2)
				call random_number(y3)
				call random_number(y4)
				y1=(y_max-y_min)*y1 + y_min
				y2=(y_max-y_min)*y2 + y_min
				y3=(y_max-y_min)*y3 + y_min
				y4=(y_max-y_min)*y4 + y_min
!				d1a = d1_cases(y1,y2,y3,y4)
				d1a = d2_cases(y1,y2,y3,y4)
				d1b = d3_cases(y1,y2,y3,y4)
!				d2a = PI1_12_i(y1,y2,y3,y4)
!				d3a = PI1_13_i(y1,y2,y3,y4)
			end do
			call toc(timer1, "<reset>")

			call tic(timer1)
			do ix1=1, 100000000
				call random_number(y1)
				call random_number(y2)
				call random_number(y3)
				call random_number(y4)
				y1=(y_max-y_min)*y1 + y_min
				y2=(y_max-y_min)*y2 + y_min
				y3=(y_max-y_min)*y3 + y_min
				y4=(y_max-y_min)*y4 + y_min
!				d1a = d1_cases(y1,y2,y3,y4)
				d1a = d2_cases(y1,y2,y3,y4)
				d1b = d3_cases(y1,y2,y3,y4)
!				d2a = PI1_12_i(y1,y2,y3,y4)
!				d3a = PI1_13_i(y1,y2,y3,y4)
			end do
			call toc(timer1, "<cases>")

			call tic(timer1)
			do ix1=1, 100000000
				call random_number(y1)
				call random_number(y2)
				call random_number(y3)
				call random_number(y4)
				y1=(y_max-y_min)*y1 + y_min
				y2=(y_max-y_min)*y2 + y_min
				y3=(y_max-y_min)*y3 + y_min
				y4=(y_max-y_min)*y4 + y_min
!				d1a = d1_full(y1,y2,y3,y4)
				d1a = d2_full(y1,y2,y3,y4)
				d1b = d3_full(y1,y2,y3,y4)
!				d2a = PI1_12_full(y1,y2,y3,y4)
!				d3a = PI1_13_full(y1,y2,y3,y4)
			end do
			call toc(timer1, "<full>")

			call tic(timer1)
			do ix1=1, 100000000
				call random_number(y1)
				call random_number(y2)
				call random_number(y3)
				call random_number(y4)
				y1=(y_max-y_min)*y1 + y_min
				y2=(y_max-y_min)*y2 + y_min
				y3=(y_max-y_min)*y3 + y_min
				y4=(y_max-y_min)*y4 + y_min
!				d1a = D1_bis(y1,y2,y3,y4)
				d1a = D2_bis(y1,y2,y3,y4)
				d1b = D3_bis(y1,y2,y3,y4)
!				d2a = PI1_12_full(y1,y2,y3,y4)
!				d3a = PI1_13_full(y1,y2,y3,y4)
			end do
			call toc(timer1, "<pablo>")
		end if

		call addToLog("[interactions] ...done!")
	end subroutine init_interp_d123

	elemental function D1_bis(y1, y2, y3, y4)
		implicit none
		real(dl) :: D1_bis
		real(dl), intent(in) :: y1, y2, y3, y4
		real(dl) a12,a34
		a12=y1+y2
		a34=y3+y4
		D1_bis=a12+a34-abs(a12-a34)-abs(y1-y2+y3-y4)-abs(y1-y2-y3+y4)
	end function D1_bis

	function D2_f(y1, y2, y3, y4)
		implicit none
		real(dl) :: D2_f
		real(dl), intent(in) :: y1, y2, y3, y4
		integer :: iflag

		D2_f=0.d0
		call D2_interp%evaluate(y1,y2,y3,y4, D2_f)
	end function D2_f

	elemental function D2_cases(q1, q2, q3, q4)
	!10.1103/PhysRevD.94.033009 eq.D2
		implicit none
		real(dl) :: D2_cases
		real(dl), intent(in) :: q1, q2, q3, q4
		real(dl) :: y1, y2, y3, y4
		real(dl) :: a

		if (q3.le.q4) then
			y1=q4
			y2=q3
		else
			y1=q3
			y2=q4
		end if
		if (q1.le.q2) then
			y3=q2
			y4=q1
		else
			y3=q1
			y4=q2
		end if
		if (y1+y2 .ge. y3+y4) then
			if (y1+y4 .ge. y2+y3) then
				if (y1 .le. y2+y3+y4) then
					a = (y1-y2)
					D2_cases = - ( a**3 &
						+ 2.d0 * (y3**3 + y4**3) &
						- 3.d0 * a * (y3**2 + y4**2) &
						) / 3.d0
				else
					D2_cases = 0.d0
				end if
			else
				D2_cases = - y4**3 / 0.75d0
			end if
		else
			if (y1+y4 .le. y2+y3) then
				if (y3 .le. y2+y1+y4) then
					a = (y1+y2)
					D2_cases = - ( - a**3 &
						- 2.d0 * (y3**3 - y4**3) &
						+ 3.d0 * a * (y3**2 + y4**2) &
						) / 3.d0
				else
					D2_cases = 0.d0
				end if
			else
				D2_cases = - y2 / 1.5d0 * ( &
					3.d0 * (y3**2 + y4**2 - y1**2) &
					- y2**2 )
			end if
		end if
	end function D2_cases

    elemental real(dl) function D2_bis(y1,y2,y3,y4)
		implicit none
		real(dl),intent(in)::y1,y2,y3,y4
		real(dl) a12,a13,a14,a34,a123,a124,a134,a234,a1234,s12,s13,s14
		real(dl) a12_2,sum1,sum2

		a12_2=y1+y2
		a34=y3+y4

		a123=a12_2+y3-y4
		a124=a12_2-y3+y4
		a134=y1-y2+a34
		a234=-y1+y2+a34
		a1234=a12_2+a34
		a12=a12_2-a34
		a13=y1-y2+y3-y4
		a14=y1-y2-y3+y4
		s12=a12*a12*sign(1.d0,a12)
		s13=a13*a13*sign(1.d0,a13)
		s14=a14*a14*sign(1.d0,a14)

		sum1=-a123*a123 + a1234*a1234 - a124*a124 + s12
		sum2=- a134*a134 + a234*a234 + s13 + s14

		D2_bis=(a123**3 - a1234**3 + a124**3 + a134**3 + a234**3 + &
			3.d0*(sum1 + sum2)*y1 + &
			3.d0*(sum1 - sum2)*y2 - abs(a12)**3 - abs(a13)**3 - &
			abs(a14)**3 + 6.d0*y1*y2*(a12_2 - 3.d0*a34 - abs(a12) + &
			abs(a13) + abs(a14)))/6.d0
	end function D2_bis

	function D3_f(y1, y2, y3, y4)
		implicit none
		real(dl) :: D3_f
		real(dl), intent(in) :: y1, y2, y3, y4
		integer :: iflag

		D3_f=0.d0
		call D3_interp%evaluate(y1,y2,y3,y4, D3_f)
	end function D3_f

	elemental function D3_cases(q1, q2, q3, q4)
	!10.1103/PhysRevD.94.033009
		implicit none
		real(dl) :: D3_cases
		real(dl), intent(in) :: q1, q2, q3, q4
		real(dl) :: y1, y2, y3, y4
		real(dl) :: a

		if (q3.le.q4) then
			y1=q4
			y2=q3
		else
			y1=q3
			y2=q4
		end if
		if (q1.le.q2) then
			y3=q2
			y4=q1
		else
			y3=q1
			y4=q2
		end if
		if (y1+y2 .ge. y3+y4) then
			if (y1+y4 .ge. y2+y3) then
				if (y1 .le. y2+y3+y4) then
					D3_cases = 1.d0 / 15.d0 * ( &
						y1**5 - y2**5 - y3**5 - y4**5 + 5.d0 * &
						( y1**2 * (y2**3 + y3**3 + y4**3) &
						- y1**3 * (y2**2 + y3**2 + y4**2) &
						+ y2**2 * y3**2 * (y2 + y3) &
						+ y2**2 * y4**2 * (y2 + y4) &
						+ y4**2 * y3**2 * (y4 + y3) ) )
				else
					D3_cases = 0.d0
				end if
			else
				D3_cases = y4**3 / 7.5d0 * (5.d0 * (y1**2 + y2**2 + y3**2) - y4**2 )
			end if
		else
			if (y1+y4 .le. y2+y3) then
				if (y3 .le. y2+y1+y4) then
					D3_cases = 1.d0 / 15.d0 * ( &
						y3**5 - y4**5 - y1**5 - y2**5 + 5.d0 * &
						( y3**2 * (y4**3 + y1**3 + y2**3) &
						- y3**3 * (y4**2 + y1**2 + y2**2) &
						+ y4**2 * y1**2 * (y4 + y1) &
						+ y4**2 * y2**2 * (y4 + y2) &
						+ y2**2 * y1**2 * (y2 + y1) ) )
				else
					D3_cases = 0.d0
				end if
			else
				D3_cases = y2**3 / 7.5d0 * (5.d0 * (y1**2 + y3**2 + y4**2) - y2**2 )
			end if
		end if

	end function D3_cases

    elemental real(dl) function D3_bis(y1,y2,y3,y4)
		implicit none
		real(dl),intent(in)::y1,y2,y3,y4
		real(dl) a12,a13,a14,a123,a124,a134,a234,a1234,P1234,s12,s13,s14
		real(dl) absa12,absa13,absa14,absa12e3,absa13e3,absa14e3
		real(dl) ampp,appm,apmp,a123e2,a124e2,a134e2,a234e2,a1234e2
		real(dl) a123e3,a124e3,a134e3,a234e3,a1234e3
		real(dl) a123e5,a124e5,a134e5,a234e5,a1234e5

		a12=y1+y2-y3-y4
		a13=y1-y2+y3-y4
		a14=y1-y2-y3+y4
		a123=y1+y2+y3-y4
		a124=y1+y2-y3+y4
		a134=y1-y2+y3+y4
		a234=-y1+y2+y3+y4
		a1234=y1+y2+y3+y4
		s12=a12**2*sign(1.0d0,a12)
		s13=a13**2*sign(1.0d0,a13)
		s14=a14**2*sign(1.0d0,a14)
		P1234=-120.d0*y1*y2*y3*y4

		absa12=abs(a12)
		absa13=abs(a13)
		absa14=abs(a14)
		absa12e3=absa12**3
		absa13e3=absa13**3
		absa14e3=absa14**3
		ampp=-absa12e3 + absa13e3 + absa14e3
		appm=absa12e3 + absa13e3 - absa14e3
		apmp=absa12e3 - absa13e3 + absa14e3
		a123e2=a123*a123
		a124e2=a124*a124
		a134e2=a134*a134
		a234e2=a234*a234
		a1234e2=a1234*a1234
		a123e3=a123e2*a123
		a124e3=a124e2*a124
		a134e3=a134e2*a134
		a234e3=a234e2*a234
		a1234e3=a1234e2*a1234

		a123e5=a123e3*a123e2
		a124e5=a124e3*a124e2
		a134e5=a134e3*a134e2
		a234e5=a234e3*a234e2
		a1234e5=a1234e3*a1234e2

		D3_bis = (&
			4.d0*(a1234e5 - a123e5 - a124e5 - a134e5 - a234e5) - &
			absa12**5 - absa13**5 - absa14**5 + &
			(a123 + a1234 + a124 + a134 + a234 + absa12 + absa13 + absa14)*P1234 + &
			5.d0*(a12**3*s12+a13**3*s13 +a14**3*s14) + &
			20.d0*((-a1234e3 + a123e3 + a124e3 - a134e3 - a234e3 + ampp)*y1*y2 + &
			(-a1234e3 + a123e3 - a124e3 + a134e3 - a234e3 + apmp)*y1*y3 + &
			(-a1234e3 + a123e3 - a124e3 - a134e3 + a234e3 + appm)*y2*y3 + &
			(-a1234e3 - a123e3 + a124e3 + a134e3 - a234e3 + appm)*y1*y4 + &
			(-a1234e3 - a123e3 + a124e3 - a134e3 + a234e3 + apmp)*y2*y4 + &
			(-a1234e3 - a123e3 - a124e3 + a134e3 + a234e3 + ampp)*y3*y4) + &
			60.d0*((a1234e2 + a123e2 - a124e2 + a134e2 + a234e2 - s12 + s13 - s14)*y1*y2*y4 + &
			(a1234e2 - a123e2 + a124e2 + a134e2 + a234e2 - s12 - s13 + s14)*y1*y2*y3 + &
			(a1234e2 + a123e2 + a124e2 - a134e2 + a234e2 + s12 - s13 - s14)*y1*y3*y4 + &
			(a1234e2 + a123e2 + a124e2 + a134e2 - a234e2 + s12 + s13 + s14)*y2*y3*y4) &
			)/120.d0
	end function D3_bis

	function PI1_12_i(y1, y2, y3, y4)
		implicit none
		real(dl) :: PI1_12_i
		real(dl), intent(in) :: y1, y2, y3, y4
		integer :: iflag

		PI1_12_i=0.d0
		call PI1_12_interp%evaluate(y1,y2,y3,y4, PI1_12_i)
	end function PI1_12_i

	function PI1_13_i(y1, y2, y3, y4)
		implicit none
		real(dl) :: PI1_13_i
		real(dl), intent(in) :: y1, y2, y3, y4
		integer :: iflag

		PI1_13_i=0.d0
		call PI1_13_interp%evaluate(y1,y2,y3,y4, PI1_13_i)
	end function PI1_13_i

	pure function PI1_full (y1, y2, y3, y4) !1: (y1,y3), 	2: (y1, y2)
		real(dl), dimension(2) :: PI1_full
		real(dl), intent(in) :: y1, y2, y3, y4
		real(dl) :: d1

		d1 = D1_full(y1, y2, y3, y4)

		!pi_1(y1, y3): nu+e -> nu+e
		PI1_full(1) = y1 * y3 * d1 + D2_full(y1, y3, y2, y4)

		!pi_1(y1, y2): nu+nu -> e+e
		PI1_full(2) = y1 * y2 * d1 - D2_full(y1, y2, y3, y4)
	end function PI1_full

	function coll_nue_4_ann_int_re(ndim, ve, obj)
		!annihilation
		integer, intent(in) :: ndim
		real(dl), dimension(2), intent(in) :: ve
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_4_ann_int_re
		integer :: ix, iy
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y2,y3,y4, f3,f4, E3, E4, dme2
		type(cmplxMatNN) :: nB

		coll_nue_4_ann_int_re = 0.d0

		dme2 = obj%dme2
		y3 = ve(1)
		E3 = Ebare_i_dme(obj%x,y3, dme2)
		f3 = fermiDirac(E3 / obj%z)
		y4 = ve(2)
		E4 = Ebare_i_dme(obj%x,y4, dme2)
		f4 = fermiDirac(E4 / obj%z)

		y2 = E3 + E4 - obj%y1
		if (.not. (y2.lt.0.d0 &
				.or. obj%y1.gt.y3+y2+y4 &
				.or. y2.gt.obj%y1+y3+y4 &
				.or. y3.gt.obj%y1+y2+y4 &
				.or. y4.gt.obj%y1+y2+y3)) then
			call allocateCmplxMat(nB)
			do ix=1, flavorNumber
				nB%re(ix,ix) = interpNuDens%re(ix,ix)%evaluate(y2)
				nB%im(ix,ix) = 0.d0
				do iy=ix+1, flavorNumber
					nB%re(ix,iy) = interpNuDens%re(ix,iy)%evaluate(y2)
					nB%im(ix,iy) = interpNuDens%im(ix,iy)%evaluate(y2)
					nB%re(iy,ix) = nB%re(ix,iy)
					nB%im(iy,ix) = -nB%im(ix,iy)
				end do
			end do
			pi2_vec = PI2_nn_f (obj%y1, y2, y3, y4, E3, E4)
			coll_nue_4_ann_int_re = coll_nue_4_ann_int_re + &
				y3/E3 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab_ann_re(nuDensMatVecFD(obj%iy),nB,f3,f4, 1, 1, obj%ix1,obj%ix2) + &
					pi2_vec(2) * F_ab_ann_re(nuDensMatVecFD(obj%iy),nB,f3,f4, 2, 2, obj%ix1,obj%ix2) + &
					(obj%x*obj%x + dme2) * PI1_12_full(obj%y1, y2, y3, y4) * ( &
						F_ab_ann_re(nuDensMatVecFD(obj%iy),nB,f3,f4, 2, 1, obj%ix1,obj%ix2) + &
						F_ab_ann_re(nuDensMatVecFD(obj%iy),nB,f3,f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_4_ann_int_re

	function coll_nue_4_sc_int_re(ndim, ve, obj)
		!scattering, summing positron and electrons
		integer, intent(in) :: ndim
		real(dl), dimension(2), intent(in) :: ve
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_4_sc_int_re
		integer :: ix, iy
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y2,y3,y4, f2,f4, E2, E4, dme2
		type(cmplxMatNN) :: nB

		coll_nue_4_sc_int_re = 0.d0

		dme2 = obj%dme2
		y2 = ve(1)
		E2 = Ebare_i_dme(obj%x,y2, dme2)
		f2 = fermiDirac(E2 / obj%z)
		y4 = ve(2)
		E4 = Ebare_i_dme(obj%x,y4, dme2)
		f4 = fermiDirac(E4 / obj%z)

		y3 = obj%y1 + E2 - E4
		if (.not. (y3.lt.0.d0 &
				.or. obj%y1.gt.y2+y3+y4 &
				.or. y2.gt.obj%y1+y3+y4 &
				.or. y3.gt.obj%y1+y2+y4 &
				.or. y4.gt.obj%y1+y2+y3)) then
			call allocateCmplxMat(nB)
			do ix=1, flavorNumber
				nB%re(ix,ix) = interpNuDens%re(ix,ix)%evaluate(y3)
				nB%im(ix,ix) = 0.d0
				do iy=ix+1, flavorNumber
					nB%re(ix,iy) = interpNuDens%re(ix,iy)%evaluate(y3)
					nB%im(ix,iy) = interpNuDens%im(ix,iy)%evaluate(y3)
					nB%re(iy,ix) = nB%re(ix,iy)
					nB%im(iy,ix) = -nB%im(ix,iy)
				end do
			end do
			pi2_vec = PI2_ne_f (obj%y1, y2, y3, y4, E2, E4)
			coll_nue_4_sc_int_re = coll_nue_4_sc_int_re + &
				y2/E2 * &
				y4/E4 * &
				( &
					( pi2_vec(1) + pi2_vec(2) ) * ( & !F_sc^LL + F_sc^RR
						F_ab_sc_re(nuDensMatVecFD(obj%iy),nB,f2,f4, 1, 1, obj%ix1,obj%ix2) + &
						F_ab_sc_re(nuDensMatVecFD(obj%iy),nB,f2,f4, 2, 2, obj%ix1,obj%ix2) &
					) - &
					2.d0 * (obj%x*obj%x + dme2) * PI1_13_full(obj%y1, y2, y3, y4) * ( & !F_sc^RL and F_sc^LR
						F_ab_sc_re(nuDensMatVecFD(obj%iy),nB,f2,f4, 2, 1, obj%ix1,obj%ix2) + &
						F_ab_sc_re(nuDensMatVecFD(obj%iy),nB,f2,f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_4_sc_int_re

	function coll_nue_4_int_re(ndim, ve, obj)
		integer, intent(in) :: ndim
		real(dl), dimension(2), intent(in) :: ve
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_4_int_re
		coll_nue_4_int_re = &
			coll_nue_4_sc_int_re(ndim, ve, obj) &
			+ coll_nue_4_ann_int_re(ndim, ve, obj)
	end function coll_nue_4_int_re

	function coll_nue_int_im(ndim, ve, obj)
		integer, intent(in) :: ndim
		real(dl), dimension(2), intent(in) :: ve
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_int_im
		integer :: ix, iy
		real(dl), dimension(2) :: pi2_vec
		type(cmplxMatNN) :: nX
		real(dl) :: yA, y2,y3,y4, f2,f3,f4, E2, E3, E4, dme2

		coll_nue_int_im = 0.0
		call allocateCmplxMat(nX)

		yA = ve(1)
		y4 = ve(2)
		dme2 = obj%dme2
		E4 = Ebare_i_dme(obj%x,y4, dme2)

		!scattering, summing positrons and electrons
		y2 = yA
		E2 = Ebare_i_dme(obj%x,y2, dme2)
		y3 = obj%y1 + E2 - E4
		f2 = fermiDirac(E2 / obj%z)
		f4 = fermiDirac(E4 / obj%z)
		if (y3.lt.0.d0 &
			.or. obj%y1.gt.y2+y3+y4 &
			.or. y2.gt.obj%y1+y3+y4 &
			.or. y3.gt.obj%y1+y2+y4 &
			.or. y4.gt.obj%y1+y2+y3) then
			coll_nue_int_im = coll_nue_int_im + 0.d0
		else
!			fd = fermiDirac(y3 / obj%z)
			do ix=1, flavorNumber
				nX%re(ix,ix) = interpNuDens%re(ix,ix)%evaluate(y3)! * fd
				nX%im(ix,ix) = 0.d0
				do iy=ix+1, flavorNumber
					nX%re(ix,iy) = interpNuDens%re(ix,iy)%evaluate(y3)
					nX%im(ix,iy) = interpNuDens%im(ix,iy)%evaluate(y3)
					nX%re(iy,ix) = nX%re(ix,iy)
					nX%im(iy,ix) = -nX%im(ix,iy)
				end do
			end do

			pi2_vec = PI2_ne_f (obj%y1, y2, y3, y4, E2, E4)

			coll_nue_int_im = coll_nue_int_im + &
				y2/E2 * &
				y4/E4 * &
				( &
					( pi2_vec(1) * pi2_vec(2) ) * ( & !F_sc^LL + F_sc^RR
						F_ab_sc_im(nuDensMatVecFD(obj%iy),nX,f2,f4, 1, 1, obj%ix1,obj%ix2) + &
						F_ab_sc_im(nuDensMatVecFD(obj%iy),nX,f2,f4, 2, 2, obj%ix1,obj%ix2) &
					) - &
					2.d0 * (obj%x*obj%x + dme2) * PI1_13_full(obj%y1, y2, y3, y4) * ( &!F_sc^RL and F_sc^LR
						F_ab_sc_im(nuDensMatVecFD(obj%iy),nX,f2,f4, 2,1, obj%ix1,obj%ix2) + &
						F_ab_sc_im(nuDensMatVecFD(obj%iy),nX,f2,f4, 1,2, obj%ix1,obj%ix2) ) &
				)
		end if

		!annihilation
		y3=yA
		E3 = Ebare_i_dme(obj%x,y3, dme2)
		y2 = E3 + E4 - obj%y1
		f3 = fermiDirac(E3 / obj%z)
		f4 = fermiDirac(E4 / obj%z)
		if (y2.lt.0.d0 &
			.or. obj%y1.gt.y3+y2+y4 &
			.or. y2.gt.obj%y1+y3+y4 &
			.or. y3.gt.obj%y1+y2+y4 &
			.or. y4.gt.obj%y1+y2+y3) then
			coll_nue_int_im = coll_nue_int_im + 0.d0
		else
!			fd = fermiDirac(y2 / obj%z)
			do ix=1, flavorNumber
				nX%re(ix,ix) = interpNuDens%re(ix,ix)%evaluate(y2)! * fd
				nX%im(ix,ix) = 0.d0
				do iy=ix+1, flavorNumber
					nX%re(ix,iy) = interpNuDens%re(ix,iy)%evaluate(y2)
					nX%im(ix,iy) = interpNuDens%im(ix,iy)%evaluate(y2)
					nX%re(iy,ix) = nX%re(ix,iy)
					nX%im(iy,ix) = -nX%im(ix,iy)
				end do
			end do

			pi2_vec = PI2_nn_f (obj%y1, y2, y3, y4, E3, E4)

			coll_nue_int_im = coll_nue_int_im + &
				y3/E3 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab_ann_im(nuDensMatVecFD(obj%iy),nX,f3,f4, 1, 1, obj%ix1,obj%ix2) + &
					pi2_vec(2) * F_ab_ann_im(nuDensMatVecFD(obj%iy),nX,f3,f4, 2, 2, obj%ix1,obj%ix2) + &
					(obj%x*obj%x + dme2) * PI1_12_full(obj%y1, y2, y3, y4) * ( &
						F_ab_ann_im(nuDensMatVecFD(obj%iy),nX,f3,f4, 2, 1, obj%ix1,obj%ix2) + &
						F_ab_ann_im(nuDensMatVecFD(obj%iy),nX,f3,f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_int_im

	pure function coll_nue_ann_int_w(iy, yx, obj, F_ab_ann, F_ab_sc)
		use fpInterfaces1
		implicit None
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		integer, intent(in) :: iy
		real(dl), intent(in) :: yx
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_ann_int_w
		coll_nue_ann_int_w = coll_nue_ann_int(iy, yx, obj, F_ab_ann)
	end function coll_nue_ann_int_w

	pure function coll_nue_sc_int_w(iy, yx, obj, F_ab_ann, F_ab_sc)
		use fpInterfaces1
		implicit None
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		integer, intent(in) :: iy
		real(dl), intent(in) :: yx
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_sc_int_w
		coll_nue_sc_int_w = coll_nue_sc_int(iy, yx, obj, F_ab_sc)
	end function coll_nue_sc_int_w

	pure SUBROUTINE region(ndim,x,j,c,d)
		use precision
		IMPLICIT NONE
		REAL (dl), INTENT (OUT) :: c, d
		INTEGER, INTENT (IN)    :: j, ndim
		REAL (dl), INTENT (IN)  :: x(ndim)
		c = y_min
		d = y_max
		RETURN
	END SUBROUTINE region
end module fpStuff
