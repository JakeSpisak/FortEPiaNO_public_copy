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
	
	!interpolation 2d parts are copied from CosmoMC/Interpolation.f90
    Type, abstract :: TSaveLoadStateObject
    contains
    procedure :: SaveState
    procedure :: LoadState
    end Type TSaveLoadStateObject
    
    Type, extends(TSaveLoadStateObject) :: TInterpolator
        logical :: initialized =.false.
        character(len=50) :: nameString
    contains
    procedure :: FirstUse => TInterpolator_FirstUse
    procedure :: Error => TInterpolator_Error
    end Type TInterpolator
	
    Type, extends(TInterpolator) :: TSpline1D
        REAL(GI), private, allocatable :: y2(:)
        REAL(GI), allocatable :: x(:), y(:)
        integer :: nx=0
    contains
    procedure :: Init => TSpline1D_Init
    procedure :: Value => TSpline1D_Value !one point
    procedure :: Clear => TSpline1D_Clear
    FINAL :: TSpline1D_Free
    end Type TSpline1D

    Type, extends(TInterpolator) :: TInterpGrid2D
        !      ALGORITHM 760, COLLECTED ALGORITHMS FROM ACM.
        !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
        !      VOL. 22, NO. 3, September, 1996, P.  357--361.
        REAL(GI), private, allocatable :: wk(:,:,:)
        REAL(GI), allocatable :: x(:), y(:)
        REAL(GI), allocatable :: z(:,:)
        integer :: nx=0, ny=0
    contains
    procedure :: Init => TInterpGrid2D_Init
    procedure :: Value => TInterpGrid2D_Value !one point
    procedure :: Values => TInterpGrid2D_Values !array of points
    procedure :: Clear => TInterpGrid2D_Clear
    procedure, private :: InitInterp => TInterpGrid2D_InitInterp
    FINAL :: TInterpGrid2D_Free
    end Type TInterpGrid2D
	
	contains
	
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
	
    subroutine LoadState(this)
		class(TSaveLoadStateObject) :: this
	!    class(TFileStream) :: F
    end subroutine LoadState

    subroutine SaveState(this)
		class(TSaveLoadStateObject) :: this
	!    class(TFileStream) :: F
    end subroutine SaveState
	
    subroutine TInterpolator_FirstUse(this)
		class(TInterpolator) this

		this%Initialized = .true.
		call this%error('TInterpolator not initialized: '//trim(this%nameString))

    end subroutine TInterpolator_FirstUse

    subroutine TInterpolator_error(this,S,v1,v2)
		class(TInterpolator):: this
		character(LEN=*), intent(in) :: S
		class(*), intent(in), optional :: v1, v2

		call criticalError('Interpolation error! '//trim(S))!,v1,v2

    end subroutine TInterpolator_error
    
    subroutine TSpline1D_init(this, x, y, nome)
		class(TSpline1D):: this
		real(dl), intent(in)      :: x(:)
		real(dl), intent(in)      :: y(:)
		character(len=*), intent(in)      :: nome
		
		this%nameString = nome
		this%nx = size(x)
		allocate(this%x, source = x)
		allocate(this%y, source = y)
		allocate(this%y2(this%nx))
		this%Initialized = .true.
        call spline(this%x, this%y, this%nx, 1.d30, 1.d30, this%y2)
    end subroutine TSpline1D_init
    
    function TSpline1D_Value(this, x)
		class(TSpline1D):: this
		REAL(dl), INTENT(IN)      :: x
		REAL(dl) :: TSpline1D_Value
		
		if (.not. this%Initialized)  call this%FirstUse
		call splint(this%x, this%y, this%y2, this%nx, x, TSpline1D_Value)
    end function TSpline1D_Value
    
    subroutine TSpline1D_clear (this)
		class(TSpline1D):: this

		if (allocated(this%y2)) then
			deallocate(this%x)
			deallocate(this%y)
			deallocate(this%y2)
		end if
		this%Initialized = .false.
    
    end subroutine TSpline1D_clear
    
    subroutine TSpline1D_Free(this)
		Type(TSpline1D):: this

		call this%Clear()
    end subroutine TSpline1D_Free
	
    subroutine TInterpGrid2D_InitInterp(this)
    class(TInterpGrid2D):: this

    allocate(this%Wk(3,this%NX,this%NY))
    CALL rgpd3p(this%nx, this%ny, this%x, this%y, this%z, this%wk)
    this%Initialized = .true.

    end subroutine TInterpGrid2D_InitInterp

    !F2003 wrappers by AL, 2013
    subroutine TInterpGrid2D_Init(this, x, y, z, nome)
    class(TInterpGrid2D):: this
    REAL(GI), INTENT(IN)      :: x(:)
    REAL(GI), INTENT(IN)      :: y(:)
    REAL(GI), INTENT(IN)      :: z(:,:)
	character(len=*), intent(in)      :: nome
	
	this%nameString = nome

    call this%Clear()
    this%nx = size(x)
    this%ny = size(y)
    allocate(this%x, source = x)
    allocate(this%y, source = y)
    allocate(this%z, source = z)

    call this%InitInterp()

    end subroutine TInterpGrid2D_Init

    subroutine TInterpGrid2D_Clear(this)
    class(TInterpGrid2D):: this

    if (allocated(this%Wk)) then
        deallocate(this%x)
        deallocate(this%y)
        deallocate(this%Wk)
        deallocate(this%z)
    end if
    this%Initialized = .false.

    end subroutine TInterpGrid2D_Clear

    subroutine TInterpGrid2D_Free(this)
    Type(TInterpGrid2D):: this

    call this%Clear()

    end subroutine TInterpGrid2D_Free

    function TInterpGrid2D_Value(this,x,y,error) result (res)
    !Z matrix not stored internally to save mem, so must pass again
    class(TInterpGrid2D) :: this
    real(GI) res, z(1), xx(1),yy(1)
    real(GI), intent(in) :: x,y
    integer, intent(inout), optional :: error

    xx(1)=x
    yy(1)=y
    call this%Values(1,xx,yy,z,error)
    res = z(1)

    end function TInterpGrid2D_Value

    subroutine TInterpGrid2D_Values(this, nip, x,y,z, error)
    !Z matrix not stored internally to save mem, so must pass again
    class(TInterpGrid2D) :: this
    integer, intent(in) :: nip
    real(GI), intent(out):: z(*)
    real(GI), intent(in) :: x(*),y(*)
    integer, intent(inout), optional :: error
    integer md,ier

    md=2
    if (.not. this%Initialized)  call this%FirstUse

    call rgbi3p(this%Wk,md, this%nx, this%ny, this%x, this%y, this%z, nip, x, y, z, ier)
    if (present(error)) then
        error=ier
    elseif (ier/=0) then
        call this%Error('error interpolating value')
    end if

    end subroutine TInterpGrid2D_Values


    SUBROUTINE rgbi3p(Wk,md, nxd, nyd, xd, yd, zd, nip, xi, yi, zi, ier)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2003-06-11  Time: 10:11:03

    ! Rectangular-grid bivariate interpolation
    ! (a master subroutine of the RGBI3P/RGSF3P subroutine package)

    ! Hiroshi Akima
    ! U.S. Department of Commerce, NTIA/ITS
    ! Version of 1995/08

    ! This subroutine performs interpolation of a bivariate function, z(x,y), on a
    ! rectangular grid in the x-y plane.  It is based on the revised Akima method.

    ! In this subroutine, the interpolating function is a piecewise function
    ! composed of a set of bicubic (bivariate third-degree) polynomials, each
    ! applicable to a rectangle of the input grid in the x-y plane.
    ! Each polynomial is determined locally.

    ! This subroutine has the accuracy of a bicubic polynomial, i.e., it
    ! interpolates accurately when all data points lie on a surface of a
    ! bicubic polynomial.

    ! The grid lines can be unevenly spaced.

    ! The input arguments are
    !   MD  = mode of computation
    !       = 1 for new XD, YD, or ZD data (default)
    !       = 2 for old XD, YD, and ZD data,
    !   NXD = number of the input-grid data points in the x coordinate
    !         (must be 2 or greater),
    !   NYD = number of the input-grid data points in the y coordinate
    !         (must be 2 or greater),
    !   XD  = array of dimension NXD containing the x coordinates of the
    !         input-grid data points (must be in a monotonic increasing order),
    !   YD  = array of dimension NYD containing the y coordinates of the
    !         input-grid data points (must be in a monotonic increasing order),
    !   ZD  = two-dimensional array of dimension NXD*NYD
    !         containing the z(x,y) values at the input-grid data points,
    !   NIP = number of the output points at which interpolation
    !         of the z value is desired (must be 1 or greater),
    !   XI  = array of dimension NIP containing the x coordinates
    !         of the output points,
    !   YI  = array of dimension NIP containing the y coordinates
    !         of the output points.

    ! The output arguments are
    !   ZI  = array of dimension NIP where the interpolated z
    !         values at the output points are to be stored,
    !   IER = error flag
    !       = 0 for no errors
    !       = 1 for NXD = 1 or less
    !       = 2 for NYD = 1 or less
    !       = 3 for identical XD values or XD values out of sequence
    !       = 4 for identical YD values or YD values out of sequence
    !       = 5 for NIP = 0 or less.

    ! N.B. The workspace has been removed from the argument list.
    !   WK  = three dimensional array of dimension 3*NXD*NYD used internally
    !         as a work area.

    ! The very fisrt call to this subroutine and the call with a new XD, YD, and
    ! ZD array must be made with MD=1.  The call with MD=2 must be preceded by
    ! another call with the same XD, YD, and ZD arrays.  Between the call with
    ! MD=2 and its preceding call, the WK array must not be disturbed.

    ! The constant in the PARAMETER statement below is
    !   NIPIMX = maximum number of output points to be processed at a time.
    ! The constant value has been selected empirically.

    ! This subroutine calls the RGPD3P, RGLCTN, and RGPLNL subroutines.


    ! Specification statements
    !     .. Parameters ..

    INTEGER, INTENT(IN)   :: md
    INTEGER, INTENT(IN)   :: nxd
    INTEGER, INTENT(IN)   :: nyd
    REAL(GI), INTENT(IN)      :: xd(nxd)
    REAL(GI), INTENT(IN)      :: yd(nyd)
    REAL(GI), INTENT(IN)  :: zd(nxd,nyd)
    INTEGER, INTENT(IN)   :: nip
    REAL(GI), INTENT(IN)  :: xi(nip)
    REAL(GI), INTENT(IN)  :: yi(nip)
    REAL(GI), INTENT(OUT)  :: zi(nip)
    INTEGER, INTENT(OUT)  :: ier
    REAL(GI), INTENT(INOUT)  :: wk(3,nxd,nyd)

    !     ..
    !     .. Local Scalars ..
    INTEGER, PARAMETER  :: nipimx=51

    INTEGER  :: iip, ix, iy, nipi
    !     ..
    !     .. Local Arrays ..
    INTEGER  :: inxi(nipimx), inyi(nipimx)

    !     ..
    !     .. External Subroutines ..
    ! EXTERNAL         rglctn, rgpd3p, rgplnl
    !     ..
    !     .. Intrinsic Functions ..
    ! INTRINSIC        MIN
    !     ..

    ! Preliminary processing
    ! Error check
    IF (nxd <= 1) GO TO 40
    IF (nyd <= 1) GO TO 50
    DO  ix = 2,nxd
        IF (xd(ix) <= xd(ix-1)) GO TO 60
    END DO
    DO  iy = 2,nyd
        IF (yd(iy) <= yd(iy-1)) GO TO 70
    END DO
    IF (nip <= 0) GO TO 80
    ier = 0

    ! Calculation
    ! Estimates partial derivatives at all input-grid data points (for MD=1).
    IF (md /= 2) THEN
        CALL rgpd3p(nxd, nyd, xd, yd, zd, wk)
    END IF

    ! DO-loop with respect to the output point
    ! Processes NIPIMX output points, at most, at a time.
    DO  iip = 1,nip,nipimx
        nipi = MIN(nip- (iip-1),nipimx)
        ! Locates the output points.
        CALL rglctn(nxd, nyd, xd, yd, nipi, xi(iip), yi(iip), inxi, inyi)

        ! Calculates the z values at the output points.
        CALL rgplnl(nxd, nyd, xd, yd, zd, wk, nipi, xi(iip), yi(iip), inxi, inyi, &
            zi(iip))
    END DO
    RETURN

    ! Error exit
40  WRITE (*,FMT=9000)
    ier = 1
    GO TO 90
50  WRITE (*,FMT=9010)
    ier = 2
    GO TO 90
60  WRITE (*,FMT=9020) ix,xd(ix)
    ier = 3
    GO TO 90
70  WRITE (*,FMT=9030) iy,yd(iy)
    ier = 4
    GO TO 90
80  WRITE (*,FMT=9040)
    ier = 5
90  WRITE (*,FMT=9050) nxd,nyd,nip
    RETURN

    ! Format statements for error messages
9000 FORMAT (/' *** RGBI3P Error 1: NXD = 1 or less')
9010 FORMAT (/' *** RGBI3P Error 2: NYD = 1 or less')
9020 FORMAT (/' *** RGBI3P Error 3: Identical XD values or',  &
        ' XD values out of sequence'/ '    IX =', i6, ',  XD(IX) =', e11.3)
9030 FORMAT (/' *** RGBI3P Error 4: Identical YD values or',  &
        ' YD values out of sequence',/,'    IY =',i6,',  YD(IY) =', e11.3)
9040 FORMAT (/' *** RGBI3P Error 5: NIP = 0 or less')
9050 FORMAT ('    NXD =', i5,',  NYD =', i5,',  NIP =', i5/)
    END SUBROUTINE rgbi3p



    SUBROUTINE rgsf3p(wk, md, nxd, nyd, xd, yd, zd, nxi, xi, nyi, yi, zi, ier)

    ! Rectangular-grid surface fitting
    ! (a master subroutine of the RGBI3P/RGSF3P subroutine package)

    ! Hiroshi Akima
    ! U.S. Department of Commerce, NTIA/ITS
    ! Version of 1995/08

    ! This subroutine performs surface fitting by interpolating values of a
    ! bivariate function, z(x,y), on a rectangular grid in the x-y plane.
    ! It is based on the revised Akima method.

    ! In this subroutine, the interpolating function is a piecewise function
    ! composed of a set of bicubic (bivariate third-degree) polynomials, each
    ! applicable to a rectangle of the input grid in the x-y plane.
    ! Each polynomial is determined locally.

    ! This subroutine has the accuracy of a bicubic polynomial, i.e., it fits the
    ! surface accurately when all data points lie on a surface of a bicubic
    ! polynomial.

    ! The grid lines of the input and output data can be unevenly spaced.

    ! The input arguments are
    !   MD  = mode of computation
    !       = 1 for new XD, YD, or ZD data (default)
    !       = 2 for old XD, YD, and ZD data,
    !   NXD = number of the input-grid data points in the x
    !         coordinate (must be 2 or greater),
    !   NYD = number of the input-grid data points in the y
    !         coordinate (must be 2 or greater),
    !   XD  = array of dimension NXD containing the x coordinates
    !         of the input-grid data points (must be in a
    !         monotonic increasing order),
    !   YD  = array of dimension NYD containing the y coordinates
    !         of the input-grid data points (must be in a
    !         monotonic increasing order),
    !   ZD  = two-dimensional array of dimension NXD*NYD
    !         containing the z(x,y) values at the input-grid data points,
    !   NXI = number of output grid points in the x coordinate
    !         (must be 1 or greater),
    !   XI  = array of dimension NXI containing the x coordinates
    !         of the output grid points,
    !   NYI = number of output grid points in the y coordinate
    !         (must be 1 or greater),
    !   YI  = array of dimension NYI containing the y coordinates
    !         of the output grid points.

    ! The output arguments are
    !   ZI  = two-dimensional array of dimension NXI*NYI, where the interpolated
    !         z values at the output grid points are to be stored,
    !   IER = error flag
    !       = 0 for no error
    !       = 1 for NXD = 1 or less
    !       = 2 for NYD = 1 or less
    !       = 3 for identical XD values or XD values out of sequence
    !       = 4 for identical YD values or YD values out of sequence
    !       = 5 for NXI = 0 or less
    !       = 6 for NYI = 0 or less.

    ! N.B. The workspace has been removed from the argument list.
    !   WK  = three-dimensional array of dimension 3*NXD*NYD used internally
    !         as a work area.

    ! The very first call to this subroutine and the call with a new XD, YD, or
    ! ZD array must be made with MD=1.  The call with MD=2 must be preceded by
    ! another call with the same XD, YD, and ZD arrays.  Between the call with
    ! MD=2 and its preceding call, the WK array must not be disturbed.

    ! The constant in the PARAMETER statement below is
    !   NIPIMX = maximum number of output points to be processed at a time.
    ! The constant value has been selected empirically.

    ! This subroutine calls the RGPD3P, RGLCTN, and RGPLNL subroutines.


    ! Specification statements
    !     .. Parameters ..

    INTEGER, INTENT(IN)   :: md
    INTEGER, INTENT(IN)   :: nxd
    INTEGER, INTENT(IN)   :: nyd
    REAL(GI), INTENT(IN)      :: xd(nxd)
    REAL(GI), INTENT(IN)      :: yd(nyd)
    REAL(GI), INTENT(IN OUT)  :: zd(nxd,nyd)
    INTEGER, INTENT(IN)   :: nxi
    REAL(GI), INTENT(IN OUT)  :: xi(nxi)
    INTEGER, INTENT(IN)   :: nyi
    REAL(GI), INTENT(IN)      :: yi(nyi)
    REAL(GI), INTENT(IN OUT)  :: zi(nxi,nyi)
    INTEGER, INTENT(OUT)  :: ier
    REAL(GI), INTENT(INOUT)  :: wk(3,nxd,nyd)

    !     ..
    !     .. Local Scalars ..
    INTEGER, PARAMETER  :: nipimx=51

    INTEGER  :: ix, ixi, iy, iyi, nipi
    !     ..
    !     .. Local Arrays ..
    REAL(GI)     :: yii(nipimx)
    INTEGER  :: inxi(nipimx), inyi(nipimx)

    !     ..
    !     .. External Subroutines ..
    ! EXTERNAL         rglctn,rgpd3p,rgplnl
    !     ..
    !     .. Intrinsic Functions ..
    ! INTRINSIC        MIN
    !     ..

    ! Preliminary processing
    ! Error check
    IF (nxd <= 1) GO TO 60
    IF (nyd <= 1) GO TO 70
    DO  ix = 2,nxd
        IF (xd(ix) <= xd(ix-1)) GO TO 80
    END DO
    DO  iy = 2,nyd
        IF (yd(iy) <= yd(iy-1)) GO TO 90
    END DO
    IF (nxi <= 0) GO TO 100
    IF (nyi <= 0) GO TO 110
    ier = 0

    ! Calculation
    ! Estimates partial derivatives at all input-grid data points
    ! (for MD=1).
    IF (md /= 2) THEN
        CALL rgpd3p(nxd, nyd, xd, yd, zd, wk)
    END IF

    ! Outermost DO-loop with respect to the y coordinate of the output grid points
    DO  iyi = 1,nyi
        DO  ixi = 1,nipimx
            yii(ixi) = yi(iyi)
        END DO

        ! Second DO-loop with respect to the x coordinate of the output grid points
        ! Processes NIPIMX output-grid points, at most, at a time.
        DO  ixi = 1,nxi,nipimx
            nipi = MIN(nxi- (ixi-1), nipimx)
            ! Locates the output-grid points.
            CALL rglctn(nxd, nyd, xd, yd, nipi, xi(ixi), yii, inxi, inyi)

            ! Calculates the z values at the output-grid points.
            CALL rgplnl(nxd, nyd, xd, yd, zd, wk, nipi, xi(ixi), yii, inxi, inyi,  &
                zi(ixi,iyi))
        END DO
    END DO
    RETURN

    ! Error exit
60  WRITE (*,FMT=9000)
    ier = 1
    GO TO 120
70  WRITE (*,FMT=9010)
    ier = 2
    GO TO 120
80  WRITE (*,FMT=9020) ix,xd(ix)
    ier = 3
    GO TO 120
90  WRITE (*,FMT=9030) iy,yd(iy)
    ier = 4
    GO TO 120
100 WRITE (*,FMT=9040)
    ier = 5
    GO TO 120
110 WRITE (*,FMT=9050)
    ier = 6
120 WRITE (*,FMT=9060) nxd,nyd,nxi,nyi
    RETURN

    ! Format statements for error messages
9000 FORMAT (/' *** RGSF3P Error 1: NXD = 1 or less')
9010 FORMAT (/' *** RGSF3P Error 2: NYD = 1 or less')
9020 FORMAT (/' *** RGSF3P Error 3: Identical XD values or',  &
        ' XD values out of sequence',/,'    IX =',i6,',  XD(IX) =', e11.3)
9030 FORMAT (/' *** RGSF3P Error 4: Identical YD values or',  &
        ' YD values out of sequence',/,'    IY =',i6,',  YD(IY) =', e11.3)
9040 FORMAT (/' *** RGSF3P Error 5: NXI = 0 or less')
9050 FORMAT (/' *** RGSF3P Error 6: NYI = 0 or less')
9060 FORMAT ('    NXD =', i5, ',  NYD =', i5, ',  NXI =', i5,',  NYI =', i5 /)
    END SUBROUTINE rgsf3p



    !     ..
    ! Statement Function definitions
    ! z2f(xx1,xx2,zz0,zz1) = (zz1-zz0)*xx2/xx1 + zz0
    ! z3f(xx1,xx2,xx3,zz0,zz1,zz2) = ((zz2-zz0)* (xx3-xx1)/xx2 -  &
    !    (zz1-zz0)* (xx3-xx2)/xx1)* (xx3/ (xx2-xx1)) + zz0

    FUNCTION z2f(xx1, xx2, zz0, zz1) RESULT(fn_val)

    REAL(GI), INTENT(IN)  :: xx1, xx2, zz0, zz1
    REAL(GI)              :: fn_val

    fn_val = (zz1-zz0)*xx2/xx1 + zz0
    RETURN
    END FUNCTION z2f



    FUNCTION z3f(xx1, xx2, xx3, zz0, zz1, zz2) RESULT(fn_val)

    REAL(GI), INTENT(IN)  :: xx1, xx2, xx3, zz0, zz1, zz2
    REAL(GI)              :: fn_val

    fn_val = ((zz2-zz0)*(xx3-xx1)/xx2 - (zz1-zz0)*(xx3-xx2)/xx1) *  &
        (xx3/(xx2-xx1)) + zz0
    RETURN
    END FUNCTION z3f



    SUBROUTINE rgpd3p(nxd, nyd, xd, yd, zd, pdd)

    ! Partial derivatives of a bivariate function on a rectangular grid
    ! (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)

    ! Hiroshi Akima
    ! U.S. Department of Commerce, NTIA/ITS
    ! Version of 1995/08

    ! This subroutine estimates three partial derivatives, zx, zy, and
    ! zxy, of a bivariate function, z(x,y), on a rectangular grid in
    ! the x-y plane.  It is based on the revised Akima method that has
    ! the accuracy of a bicubic polynomial.

    ! The input arguments are
    !   NXD = number of the input-grid data points in the x
    !         coordinate (must be 2 or greater),
    !   NYD = number of the input-grid data points in the y
    !         coordinate (must be 2 or greater),
    !   XD  = array of dimension NXD containing the x coordinates of the
    !         input-grid data points (must be in a monotonic increasing order),
    !   YD  = array of dimension NYD containing the y coordinates of the
    !         input-grid data points (must be in a monotonic increasing order),
    !   ZD  = two-dimensional array of dimension NXD*NYD
    !         containing the z(x,y) values at the input-grid data points.

    ! The output argument is
    !   PDD = three-dimensional array of dimension 3*NXD*NYD,
    !         where the estimated zx, zy, and zxy values at the
    !         input-grid data points are to be stored.


    ! Specification statements
    !     .. Scalar Arguments ..

    INTEGER, INTENT(IN)  :: nxd
    INTEGER, INTENT(IN)  :: nyd
    REAL(GI), INTENT(IN)     :: xd(nxd)
    REAL(GI), INTENT(IN)     :: yd(nyd)
    REAL(GI), INTENT(IN)     :: zd(nxd,nyd)
    REAL(GI), INTENT(OUT)    :: pdd(3,nxd,nyd)

    !     ..
    !     .. Local Scalars ..
    REAL(GI) :: b00, b00x, b00y, b01, b10, b11, cx1, cx2, cx3, cy1, cy2,  &
        cy3, disf, dnm, dz00, dz01, dz02, dz03, dz10, dz11, dz12,  &
        dz13, dz20, dz21, dz22, dz23, dz30, dz31, dz32, dz33,  &
        dzx10, dzx20, dzx30, dzxy11, dzxy12, dzxy13, dzxy21,  &
        dzxy22, dzxy23, dzxy31, dzxy32, dzxy33, dzy01, dzy02,  &
        dzy03, epsln, pezx, pezxy, pezy, smpef, smpei, smwtf,  &
        smwti, sx, sxx, sxxy, sxxyy, sxy, sxyy, sxyz, sxz, sy, syy,  &
        syz, sz, volf, wt, x0, x1, x2, x3, y0, y1, y2,  &
        y3, z00, z01, z02, z03, z10, z11, z12, z13, z20, z21, z22,  &
        z23, z30, z31, z32, z33, zxdi, zxydi, zydi
    INTEGER :: ipex, ipey, ix0, ix1, ix2, ix3, iy0, iy1, iy2, iy3, nx0, ny0
    !     ..
    !     .. Local Arrays ..
    REAL(GI)    :: b00xa(4), b00ya(4), b01a(4), b10a(4), cxa(3,4), cya(3,4),   &
        sxa(4), sxxa(4), sya(4), syya(4), xa(3,4), ya(3,4),   &
        z0ia(3,4), zi0a(3,4)

    INTEGER, parameter :: idlt(3,4) = RESHAPE([-3, -2, -1, -2, -1, 1,-1, 1,2, 1 ,2, 3 ], [ 3, 4 ] )
    !AL Jun 14: Fixed F90 translation

    ! Calculation
    ! Initial setting of some local variables
    nx0 = MAX(4,nxd)
    ny0 = MAX(4,nyd)

    ! Double DO-loop with respect to the input grid points
    DO  iy0 = 1,nyd
        DO  ix0 = 1,nxd
            x0 = xd(ix0)
            y0 = yd(iy0)
            z00 = zd(ix0,iy0)

            ! Part 1.  Estimation of ZXDI
            ! Initial setting
            smpef = 0.0
            smwtf = 0.0
            smpei = 0.0
            smwti = 0.0
            ! DO-loop with respect to the primary estimate
            DO  ipex = 1,4
                ! Selects necessary grid points in the x direction.
                ix1 = ix0 + idlt(1,ipex)
                ix2 = ix0 + idlt(2,ipex)
                ix3 = ix0 + idlt(3,ipex)
                IF ((ix1 < 1) .OR. (ix2 < 1) .OR. (ix3 < 1) .OR.  &
                    (ix1 > nx0) .OR. (ix2 > nx0) .OR. (ix3 > nx0)) CYCLE
                ! Selects and/or supplements the x and z values.
                x1 = xd(ix1) - x0
                z10 = zd(ix1,iy0)
                IF (nxd >= 4) THEN
                    x2 = xd(ix2) - x0
                    x3 = xd(ix3) - x0
                    z20 = zd(ix2,iy0)
                    z30 = zd(ix3,iy0)
                ELSE IF (nxd == 3) THEN
                    x2 = xd(ix2) - x0
                    z20 = zd(ix2,iy0)
                    x3 = 2*xd(3) - xd(2) - x0
                    z30 = z3f(x1,x2,x3,z00,z10,z20)
                ELSE IF (nxd == 2) THEN
                    x2 = 2*xd(2) - xd(1) - x0
                    z20 = z2f(x1,x2,z00,z10)
                    x3 = 2*xd(1) - xd(2) - x0
                    z30 = z2f(x1,x3,z00,z10)
                END IF
                dzx10 = (z10-z00)/x1
                dzx20 = (z20-z00)/x2
                dzx30 = (z30-z00)/x3
                ! Calculates the primary estimate of partial derivative zx as
                ! the coefficient of the bicubic polynomial.
                cx1 = x2*x3/ ((x1-x2)* (x1-x3))
                cx2 = x3*x1/ ((x2-x3)* (x2-x1))
                cx3 = x1*x2/ ((x3-x1)* (x3-x2))
                pezx = cx1*dzx10 + cx2*dzx20 + cx3*dzx30
                ! Calculates the volatility factor and distance factor in the x
                ! direction for the primary estimate of zx.
                sx = x1 + x2 + x3
                sz = z00 + z10 + z20 + z30
                sxx = x1*x1 + x2*x2 + x3*x3
                sxz = x1*z10 + x2*z20 + x3*z30
                dnm = 4.0*sxx - sx*sx
                b00 = (sxx*sz-sx*sxz)/dnm
                b10 = (4.0*sxz-sx*sz)/dnm
                dz00 = z00 - b00
                dz10 = z10 - (b00+b10*x1)
                dz20 = z20 - (b00+b10*x2)
                dz30 = z30 - (b00+b10*x3)
                volf = dz00**2 + dz10**2 + dz20**2 + dz30**2
                disf = sxx

                ! Calculates the EPSLN value, which is used to decide whether or
                ! not the volatility factor is essentially zero.
                epsln = (z00**2+z10**2+z20**2+z30**2)*1.0E-12
                ! Accumulates the weighted primary estimates of zx and their weights.
                IF (volf > epsln) THEN
                    ! - For a finite weight.
                    wt = 1.0/ (volf*disf)
                    smpef = smpef + wt*pezx
                    smwtf = smwtf + wt
                ELSE
                    ! - For an infinite weight.
                    smpei = smpei + pezx
                    smwti = smwti + 1.0
                END IF

                ! Saves the necessary values for estimating zxy
                xa(1,ipex) = x1
                xa(2,ipex) = x2
                xa(3,ipex) = x3
                zi0a(1,ipex) = z10
                zi0a(2,ipex) = z20
                zi0a(3,ipex) = z30
                cxa(1,ipex) = cx1
                cxa(2,ipex) = cx2
                cxa(3,ipex) = cx3
                sxa(ipex) = sx
                sxxa(ipex) = sxx
                b00xa(ipex) = b00
                b10a(ipex) = b10
            END DO

            ! Calculates the final estimate of zx.
            IF (smwti < 0.5) THEN
                ! - When no infinite weights exist.
                zxdi = smpef/smwtf
            ELSE
                ! - When infinite weights exist.
                zxdi = smpei/smwti
            END IF
            ! End of Part 1.

            ! Part 2.  Estimation of ZYDI
            ! Initial setting
            smpef = 0.0
            smwtf = 0.0
            smpei = 0.0
            smwti = 0.0
            ! DO-loop with respect to the primary estimate
            DO  ipey = 1,4
                ! Selects necessary grid points in the y direction.
                iy1 = iy0 + idlt(1,ipey)
                iy2 = iy0 + idlt(2,ipey)
                iy3 = iy0 + idlt(3,ipey)
                IF ((iy1 < 1) .OR. (iy2 < 1) .OR. (iy3 < 1) .OR.  &
                    (iy1 > ny0) .OR. (iy2 > ny0) .OR. (iy3 > ny0)) CYCLE
                ! Selects and/or supplements the y and z values.
                y1 = yd(iy1) - y0
                z01 = zd(ix0,iy1)
                IF (nyd >= 4) THEN
                    y2 = yd(iy2) - y0
                    y3 = yd(iy3) - y0
                    z02 = zd(ix0,iy2)
                    z03 = zd(ix0,iy3)
                ELSE IF (nyd == 3) THEN
                    y2 = yd(iy2) - y0
                    z02 = zd(ix0,iy2)
                    y3 = 2*yd(3) - yd(2) - y0
                    z03 = z3f(y1,y2,y3,z00,z01,z02)
                ELSE IF (nyd == 2) THEN
                    y2 = 2*yd(2) - yd(1) - y0
                    z02 = z2f(y1,y2,z00,z01)
                    y3 = 2*yd(1) - yd(2) - y0
                    z03 = z2f(y1,y3,z00,z01)
                END IF
                dzy01 = (z01-z00)/y1
                dzy02 = (z02-z00)/y2
                dzy03 = (z03-z00)/y3
                ! Calculates the primary estimate of partial derivative zy as
                ! the coefficient of the bicubic polynomial.
                cy1 = y2*y3/ ((y1-y2)* (y1-y3))
                cy2 = y3*y1/ ((y2-y3)* (y2-y1))
                cy3 = y1*y2/ ((y3-y1)* (y3-y2))
                pezy = cy1*dzy01 + cy2*dzy02 + cy3*dzy03
                ! Calculates the volatility factor and distance factor in the y
                ! direction for the primary estimate of zy.
                sy = y1 + y2 + y3
                sz = z00 + z01 + z02 + z03
                syy = y1*y1 + y2*y2 + y3*y3
                syz = y1*z01 + y2*z02 + y3*z03
                dnm = 4.0*syy - sy*sy
                b00 = (syy*sz-sy*syz)/dnm
                b01 = (4.0*syz-sy*sz)/dnm
                dz00 = z00 - b00
                dz01 = z01 - (b00+b01*y1)
                dz02 = z02 - (b00+b01*y2)
                dz03 = z03 - (b00+b01*y3)
                volf = dz00**2 + dz01**2 + dz02**2 + dz03**2
                disf = syy

                ! Calculates the EPSLN value, which is used to decide whether or
                ! not the volatility factor is essentially zero.
                epsln = (z00**2+z01**2+z02**2+z03**2)*1.0E-12
                ! Accumulates the weighted primary estimates of zy and their weights.
                IF (volf > epsln) THEN
                    ! - For a finite weight.
                    wt = 1.0/ (volf*disf)
                    smpef = smpef + wt*pezy
                    smwtf = smwtf + wt
                ELSE
                    ! - For an infinite weight.
                    smpei = smpei + pezy
                    smwti = smwti + 1.0
                END IF
                ! Saves the necessary values for estimating zxy
                ya(1,ipey) = y1
                ya(2,ipey) = y2
                ya(3,ipey) = y3
                z0ia(1,ipey) = z01
                z0ia(2,ipey) = z02
                z0ia(3,ipey) = z03
                cya(1,ipey) = cy1
                cya(2,ipey) = cy2
                cya(3,ipey) = cy3
                sya(ipey) = sy
                syya(ipey) = syy
                b00ya(ipey) = b00
                b01a(ipey) = b01
            END DO

            ! Calculates the final estimate of zy.
            IF (smwti < 0.5) THEN
                ! - When no infinite weights exist.
                zydi = smpef/smwtf
            ELSE
                ! - When infinite weights exist.
                zydi = smpei/smwti
            END IF
            ! End of Part 2.

            ! Part 3.  Estimation of ZXYDI
            ! Initial setting
            smpef = 0.0
            smwtf = 0.0
            smpei = 0.0
            smwti = 0.0
            ! Outer DO-loops with respect to the primary estimates in the x direction
            DO  ipex = 1,4
                ix1 = ix0 + idlt(1,ipex)
                ix2 = ix0 + idlt(2,ipex)
                ix3 = ix0 + idlt(3,ipex)
                IF ((ix1 < 1) .OR. (ix2 < 1) .OR. (ix3 < 1) .OR.  &
                    (ix1 > nx0) .OR. (ix2 > nx0) .OR. (ix3 > nx0)) CYCLE
                ! Retrieves the necessary values for estimating zxy in the x direction.
                x1 = xa(1,ipex)
                x2 = xa(2,ipex)
                x3 = xa(3,ipex)
                z10 = zi0a(1,ipex)
                z20 = zi0a(2,ipex)
                z30 = zi0a(3,ipex)
                cx1 = cxa(1,ipex)
                cx2 = cxa(2,ipex)
                cx3 = cxa(3,ipex)
                sx = sxa(ipex)
                sxx = sxxa(ipex)
                b00x = b00xa(ipex)
                b10 = b10a(ipex)

                ! Inner DO-loops with respect to the primary estimates in the y direction
                DO  ipey = 1,4
                    iy1 = iy0 + idlt(1,ipey)
                    iy2 = iy0 + idlt(2,ipey)
                    iy3 = iy0 + idlt(3,ipey)
                    IF ((iy1 < 1) .OR. (iy2 < 1) .OR. (iy3 < 1) .OR. (iy1 > ny0) .OR.  &
                        (iy2 > ny0) .OR. (iy3 > ny0)) CYCLE
                    ! Retrieves the necessary values for estimating zxy in the y direction.
                    y1 = ya(1,ipey)
                    y2 = ya(2,ipey)
                    y3 = ya(3,ipey)
                    z01 = z0ia(1,ipey)
                    z02 = z0ia(2,ipey)
                    z03 = z0ia(3,ipey)
                    cy1 = cya(1,ipey)
                    cy2 = cya(2,ipey)
                    cy3 = cya(3,ipey)
                    sy = sya(ipey)
                    syy = syya(ipey)
                    b00y = b00ya(ipey)
                    b01 = b01a(ipey)
                    ! Selects and/or supplements the z values.
                    IF (nyd >= 4) THEN
                        z11 = zd(ix1,iy1)
                        z12 = zd(ix1,iy2)
                        z13 = zd(ix1,iy3)
                        IF (nxd >= 4) THEN
                            z21 = zd(ix2,iy1)
                            z22 = zd(ix2,iy2)
                            z23 = zd(ix2,iy3)
                            z31 = zd(ix3,iy1)
                            z32 = zd(ix3,iy2)
                            z33 = zd(ix3,iy3)
                        ELSE IF (nxd == 3) THEN
                            z21 = zd(ix2,iy1)
                            z22 = zd(ix2,iy2)
                            z23 = zd(ix2,iy3)
                            z31 = z3f(x1,x2,x3,z01,z11,z21)
                            z32 = z3f(x1,x2,x3,z02,z12,z22)
                            z33 = z3f(x1,x2,x3,z03,z13,z23)
                        ELSE IF (nxd == 2) THEN
                            z21 = z2f(x1,x2,z01,z11)
                            z22 = z2f(x1,x2,z02,z12)
                            z23 = z2f(x1,x2,z03,z13)
                            z31 = z2f(x1,x3,z01,z11)
                            z32 = z2f(x1,x3,z02,z12)
                            z33 = z2f(x1,x3,z03,z13)
                        END IF
                    ELSE IF (nyd == 3) THEN
                        z11 = zd(ix1,iy1)
                        z12 = zd(ix1,iy2)
                        z13 = z3f(y1,y2,y3,z10,z11,z12)
                        IF (nxd >= 4) THEN
                            z21 = zd(ix2,iy1)
                            z22 = zd(ix2,iy2)
                            z31 = zd(ix3,iy1)
                            z32 = zd(ix3,iy2)
                        ELSE IF (nxd == 3) THEN
                            z21 = zd(ix2,iy1)
                            z22 = zd(ix2,iy2)
                            z31 = z3f(x1,x2,x3,z01,z11,z21)
                            z32 = z3f(x1,x2,x3,z02,z12,z22)
                        ELSE IF (nxd == 2) THEN
                            z21 = z2f(x1,x2,z01,z11)
                            z22 = z2f(x1,x2,z02,z12)
                            z31 = z2f(x1,x3,z01,z11)
                            z32 = z2f(x1,x3,z02,z12)
                        END IF
                        z23 = z3f(y1,y2,y3,z20,z21,z22)
                        z33 = z3f(y1,y2,y3,z30,z31,z32)
                    ELSE IF (nyd == 2) THEN
                        z11 = zd(ix1,iy1)
                        z12 = z2f(y1,y2,z10,z11)
                        z13 = z2f(y1,y3,z10,z11)
                        IF (nxd >= 4) THEN
                            z21 = zd(ix2,iy1)
                            z31 = zd(ix3,iy1)
                        ELSE IF (nxd == 3) THEN
                            z21 = zd(ix2,iy1)
                            z31 = z3f(x1,x2,x3,z01,z11,z21)
                        ELSE IF (nxd == 2) THEN
                            z21 = z2f(x1,x2,z01,z11)
                            z31 = z2f(x1,x3,z01,z11)
                        END IF
                        z22 = z2f(y1,y2,z20,z21)
                        z23 = z2f(y1,y3,z20,z21)
                        z32 = z2f(y1,y2,z30,z31)
                        z33 = z2f(y1,y3,z30,z31)
                    END IF
                    ! Calculates the primary estimate of partial derivative zxy as
                    ! the coefficient of the bicubic polynomial.
                    dzxy11 = (z11-z10-z01+z00)/ (x1*y1)
                    dzxy12 = (z12-z10-z02+z00)/ (x1*y2)
                    dzxy13 = (z13-z10-z03+z00)/ (x1*y3)
                    dzxy21 = (z21-z20-z01+z00)/ (x2*y1)
                    dzxy22 = (z22-z20-z02+z00)/ (x2*y2)
                    dzxy23 = (z23-z20-z03+z00)/ (x2*y3)
                    dzxy31 = (z31-z30-z01+z00)/ (x3*y1)
                    dzxy32 = (z32-z30-z02+z00)/ (x3*y2)
                    dzxy33 = (z33-z30-z03+z00)/ (x3*y3)
                    pezxy = cx1* (cy1*dzxy11+cy2*dzxy12+cy3*dzxy13) +  &
                        cx2* (cy1*dzxy21+cy2*dzxy22+cy3*dzxy23) +  &
                        cx3* (cy1*dzxy31+cy2*dzxy32+cy3*dzxy33)
                    ! Calculates the volatility factor and distance factor in the x
                    ! and y directions for the primary estimate of zxy.
                    b00 = (b00x+b00y)/2.0
                    sxy = sx*sy
                    sxxy = sxx*sy
                    sxyy = sx*syy
                    sxxyy = sxx*syy
                    sxyz = x1* (y1*z11+y2*z12+y3*z13) + x2* (y1*z21+y2*z22+y3*z23) +  &
                        x3* (y1*z31+y2*z32+y3*z33)
                    b11 = (sxyz-b00*sxy-b10*sxxy-b01*sxyy)/sxxyy
                    dz00 = z00 - b00
                    dz01 = z01 - (b00+b01*y1)
                    dz02 = z02 - (b00+b01*y2)
                    dz03 = z03 - (b00+b01*y3)
                    dz10 = z10 - (b00+b10*x1)
                    dz11 = z11 - (b00+b01*y1+x1* (b10+b11*y1))
                    dz12 = z12 - (b00+b01*y2+x1* (b10+b11*y2))
                    dz13 = z13 - (b00+b01*y3+x1* (b10+b11*y3))
                    dz20 = z20 - (b00+b10*x2)
                    dz21 = z21 - (b00+b01*y1+x2* (b10+b11*y1))
                    dz22 = z22 - (b00+b01*y2+x2* (b10+b11*y2))
                    dz23 = z23 - (b00+b01*y3+x2* (b10+b11*y3))
                    dz30 = z30 - (b00+b10*x3)
                    dz31 = z31 - (b00+b01*y1+x3* (b10+b11*y1))
                    dz32 = z32 - (b00+b01*y2+x3* (b10+b11*y2))
                    dz33 = z33 - (b00+b01*y3+x3* (b10+b11*y3))
                    volf = dz00**2 + dz01**2 + dz02**2 + dz03**2 +  &
                        dz10**2 + dz11**2 + dz12**2 + dz13**2 +  &
                        dz20**2 + dz21**2 + dz22**2 + dz23**2 +  &
                        dz30**2 + dz31**2 + dz32**2 + dz33**2
                    disf = sxx*syy
                    ! Calculates EPSLN.
                    epsln = (z00**2 + z01**2 + z02**2 + z03**2 + z10**2 +   &
                        z11**2 + z12**2 + z13**2 + z20**2 + z21**2 + z22**2 +   &
                        z23**2 + z30**2 + z31**2 + z32**2 + z33**2)* 1.0E-12
                    ! Accumulates the weighted primary estimates of zxy and their weights.
                    IF (volf > epsln) THEN
                        ! - For a finite weight.
                        wt = 1.0/ (volf*disf)
                        smpef = smpef + wt*pezxy
                        smwtf = smwtf + wt
                    ELSE
                        ! - For an infinite weight.
                        smpei = smpei + pezxy
                        smwti = smwti + 1.0
                    END IF
                END DO
            END DO

            ! Calculates the final estimate of zxy.
            IF (smwti < 0.5) THEN
                ! - When no infinite weights exist.
                zxydi = smpef/smwtf
            ELSE
                ! - When infinite weights exist.
                zxydi = smpei/smwti
            END IF
            ! End of Part 3

            pdd(1,ix0,iy0) = zxdi
            pdd(2,ix0,iy0) = zydi
            pdd(3,ix0,iy0) = zxydi
        END DO
    END DO
    RETURN
    END SUBROUTINE rgpd3p



    SUBROUTINE rglctn(nxd, nyd, xd, yd, nip, xi, yi, inxi, inyi)

    ! Location of the desired points in a rectangular grid
    ! (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)

    ! Hiroshi Akima
    ! U.S. Department of Commerce, NTIA/ITS
    ! Version of 1995/08

    ! This subroutine locates the desired points in a rectangular grid
    ! in the x-y plane.

    ! The grid lines can be unevenly spaced.

    ! The input arguments are
    !   NXD  = number of the input-grid data points in the x
    !          coordinate (must be 2 or greater),
    !   NYD  = number of the input-grid data points in the y
    !          coordinate (must be 2 or greater),
    !   XD   = array of dimension NXD containing the x coordinates of the
    !          input-grid data points (must be in a monotonic increasing order),
    !   YD   = array of dimension NYD containing the y coordinates of the
    !          input-grid data points (must be in a monotonic increasing order),
    !   NIP  = number of the output points to be located (must be 1 or greater),
    !   XI   = array of dimension NIP containing the x coordinates
    !          of the output points to be located,
    !   YI   = array of dimension NIP containing the y coordinates
    !          of the output points to be located.

    ! The output arguments are
    !   INXI = integer array of dimension NIP where the interval
    !          numbers of the XI array elements are to be stored,
    !   INYI = integer array of dimension NIP where the interval
    !          numbers of the YI array elements are to be stored.
    ! The interval numbers are between 0 and NXD and between 0 and NYD,
    ! respectively.


    ! Specification statements
    !     .. Scalar Arguments ..

    INTEGER, INTENT(IN)   :: nxd
    INTEGER, INTENT(IN)   :: nyd
    REAL(GI), INTENT(IN)      :: xd(nxd)
    REAL(GI), INTENT(IN)      :: yd(nyd)
    INTEGER, INTENT(IN)   :: nip
    REAL(GI), INTENT(IN)      :: xi(nip)
    REAL(GI), INTENT(IN)      :: yi(nip)
    INTEGER, INTENT(OUT)  :: inxi(nip)
    INTEGER, INTENT(OUT)  :: inyi(nip)

    !     ..
    !     .. Local Scalars ..
    REAL(GI)     :: xii, yii
    INTEGER  :: iip, imd, imn, imx, ixd, iyd, nintx, ninty

    !     ..
    ! DO-loop with respect to IIP, which is the point number of the output point
    DO  iip = 1,nip
        xii = xi(iip)
        yii = yi(iip)
        ! Checks if the x coordinate of the IIPth output point, XII, is
        ! in a new interval.  (NINTX is the new-interval flag.)
        IF (iip == 1) THEN
            nintx = 1
        ELSE
            nintx = 0
            IF (ixd == 0) THEN
                IF (xii > xd(1)) nintx = 1
            ELSE IF (ixd < nxd) THEN
                IF ((xii < xd(ixd)) .OR. (xii > xd(ixd+1))) nintx = 1
            ELSE
                IF (xii < xd(nxd)) nintx = 1
            END IF
        END IF

        ! Locates the output point by binary search if XII is in a new interval.
        ! Determines IXD for which XII lies between XD(IXD) and XD(IXD+1).
        IF (nintx == 1) THEN
            IF (xii <= xd(1)) THEN
                ixd = 0
            ELSE IF (xii < xd(nxd)) THEN
                imn = 1
                imx = nxd
                imd = (imn+imx)/2
10              IF (xii >= xd(imd)) THEN
                    imn = imd
                ELSE
                    imx = imd
                END IF
                imd = (imn+imx)/2
                IF (imd > imn) GO TO 10
                ixd = imd
            ELSE
                ixd = nxd
            END IF
        END IF
        inxi(iip) = ixd

        ! Checks if the y coordinate of the IIPth output point, YII, is
        ! in a new interval.  (NINTY is the new-interval flag.)
        IF (iip == 1) THEN
            ninty = 1
        ELSE
            ninty = 0
            IF (iyd == 0) THEN
                IF (yii > yd(1)) ninty = 1
            ELSE IF (iyd < nyd) THEN
                IF ((yii < yd(iyd)) .OR. (yii > yd(iyd+1))) ninty = 1
            ELSE
                IF (yii < yd(nyd)) ninty = 1
            END IF
        END IF

        ! Locates the output point by binary search if YII is in a new interval.
        ! Determines IYD for which YII lies between YD(IYD) and YD(IYD+1).
        IF (ninty == 1) THEN
            IF (yii <= yd(1)) THEN
                iyd = 0
            ELSE IF (yii < yd(nyd)) THEN
                imn = 1
                imx = nyd
                imd = (imn+imx)/2
20              IF (yii >= yd(imd)) THEN
                    imn = imd
                ELSE
                    imx = imd
                END IF
                imd = (imn+imx)/2
                IF (imd > imn) GO TO 20
                iyd = imd
            ELSE
                iyd = nyd
            END IF
        END IF
        inyi(iip) = iyd
    END DO
    RETURN
    END SUBROUTINE rglctn



    SUBROUTINE rgplnl(nxd, nyd, xd, yd, zd, pdd, nip, xi, yi, inxi, inyi, zi)

    ! Polynomials for rectangular-grid bivariate interpolation and surface fitting
    ! (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)

    ! Hiroshi Akima
    ! U.S. Department of Commerce, NTIA/ITS
    ! Version of 1995/08

    ! This subroutine determines a polynomial in x and y for a rectangle of the
    ! input grid in the x-y plane and calculates the z value for the desired
    ! points by evaluating the polynomial for rectangular-grid bivariate
    ! interpolation and surface fitting.

    ! The input arguments are
    !   NXD  = number of the input-grid data points in the x
    !          coordinate (must be 2 or greater),
    !   NYD  = number of the input-grid data points in the y
    !          coordinate (must be 2 or greater),
    !   XD   = array of dimension NXD containing the x coordinates of the
    !          input-grid data points (must be in a monotonic increasing order),
    !   YD   = array of dimension NYD containing the y coordinates of the
    !          input-grid data points (must be in a monotonic increasing order),
    !   ZD   = two-dimensional array of dimension NXD*NYD
    !          containing the z(x,y) values at the input-grid data points,
    !   PDD  = three-dimensional array of dimension 3*NXD*NYD
    !          containing the estimated zx, zy, and zxy values
    !          at the input-grid data points,
    !   NIP  = number of the output points at which interpolation
    !          is to be performed,
    !   XI   = array of dimension NIP containing the x coordinates
    !          of the output points,
    !   YI   = array of dimension NIP containing the y coordinates
    !          of the output points,
    !   INXI = integer array of dimension NIP containing the
    !          interval numbers of the input grid intervals in the
    !          x direction where the x coordinates of the output points lie,
    !   INYI = integer array of dimension NIP containing the
    !          interval numbers of the input grid intervals in the
    !          y direction where the y coordinates of the output points lie.

    ! The output argument is
    !   ZI   = array of dimension NIP, where the interpolated z
    !          values at the output points are to be stored.


    ! Specification statements
    !     .. Scalar Arguments ..

    INTEGER, INTENT(IN)  :: nxd
    INTEGER, INTENT(IN)  :: nyd
    REAL(GI), INTENT(IN)     :: xd(nxd)
    REAL(GI), INTENT(IN)     :: yd(nyd)
    REAL(GI), INTENT(IN)     :: zd(nxd,nyd)
    REAL(GI), INTENT(IN)     :: pdd(3,nxd,nyd)
    INTEGER, INTENT(IN)  :: nip
    REAL(GI), INTENT(IN)     :: xi(nip)
    REAL(GI), INTENT(IN)     :: yi(nip)
    INTEGER, INTENT(IN)  :: inxi(nip)
    INTEGER, INTENT(IN)  :: inyi(nip)
    REAL(GI), INTENT(OUT)    :: zi(nip)

    !     ..
    !     .. Local Scalars ..
    REAL(GI) :: a, b, c, this, dx, dxsq, dy, dysq, p00, p01, p02, p03, p10, p11,  &
        p12, p13, p20, p21, p22, p23, p30, p31, p32, p33, q0, q1, q2,  &
        q3, u, v, x0, xii, y0, yii, z00, z01, z0dx, z0dy, z10, z11,  &
        z1dx, z1dy, zdxdy, zii, zx00, zx01, zx0dy, zx10, zx11,  &
        zx1dy, zxy00, zxy01, zxy10, zxy11, zy00, zy01, zy0dx, zy10, zy11, zy1dx
    INTEGER :: iip, ixd0, ixd1, ixdi, ixdipv, iyd0, iyd1, iydi, iydipv
    !     ..
    !     .. Intrinsic Functions ..
    ! INTRINSIC        MAX
    !     ..

    ! Calculation
    ! Outermost DO-loop with respect to the output point
    DO  iip = 1,nip
        xii = xi(iip)
        yii = yi(iip)
        IF (iip == 1) THEN
            ixdipv = -1
            iydipv = -1
        ELSE
            ixdipv = ixdi
            iydipv = iydi
        END IF
        ixdi = inxi(iip)
        iydi = inyi(iip)

        ! Retrieves the z and partial derivative values at the origin of
        ! the coordinate for the rectangle.
        IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
            ixd0 = MAX(1,ixdi)
            iyd0 = MAX(1,iydi)
            x0 = xd(ixd0)
            y0 = yd(iyd0)
            z00 = zd(ixd0,iyd0)
            zx00 = pdd(1,ixd0,iyd0)
            zy00 = pdd(2,ixd0,iyd0)
            zxy00 = pdd(3,ixd0,iyd0)
        END IF

        ! Case 1.  When the rectangle is inside the data area in both the
        ! x and y directions.
        IF ((ixdi > 0 .AND. ixdi < nxd) .AND. (iydi > 0 .AND. iydi < nyd)) THEN
            ! Retrieves the z and partial derivative values at the other three
            ! vertices of the rectangle.
            IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
                ixd1 = ixd0 + 1
                dx = xd(ixd1) - x0
                dxsq = dx*dx
                iyd1 = iyd0 + 1
                dy = yd(iyd1) - y0
                dysq = dy*dy
                z10 = zd(ixd1,iyd0)
                z01 = zd(ixd0,iyd1)
                z11 = zd(ixd1,iyd1)
                zx10 = pdd(1,ixd1,iyd0)
                zx01 = pdd(1,ixd0,iyd1)
                zx11 = pdd(1,ixd1,iyd1)
                zy10 = pdd(2,ixd1,iyd0)
                zy01 = pdd(2,ixd0,iyd1)
                zy11 = pdd(2,ixd1,iyd1)
                zxy10 = pdd(3,ixd1,iyd0)
                zxy01 = pdd(3,ixd0,iyd1)
                zxy11 = pdd(3,ixd1,iyd1)
                ! Calculates the polynomial coefficients.
                z0dx = (z10-z00)/dx
                z1dx = (z11-z01)/dx
                z0dy = (z01-z00)/dy
                z1dy = (z11-z10)/dy
                zx0dy = (zx01-zx00)/dy
                zx1dy = (zx11-zx10)/dy
                zy0dx = (zy10-zy00)/dx
                zy1dx = (zy11-zy01)/dx
                zdxdy = (z1dy-z0dy)/dx
                a = zdxdy - zx0dy - zy0dx + zxy00
                b = zx1dy - zx0dy - zxy10 + zxy00
                c = zy1dx - zy0dx - zxy01 + zxy00
                this = zxy11 - zxy10 - zxy01 + zxy00
                p00 = z00
                p01 = zy00
                p02 = (2.0* (z0dy-zy00)+z0dy-zy01)/dy
                p03 = (-2.0*z0dy+zy01+zy00)/dysq
                p10 = zx00
                p11 = zxy00
                p12 = (2.0* (zx0dy-zxy00)+zx0dy-zxy01)/dy
                p13 = (-2.0*zx0dy+zxy01+zxy00)/dysq
                p20 = (2.0* (z0dx-zx00)+z0dx-zx10)/dx
                p21 = (2.0* (zy0dx-zxy00)+zy0dx-zxy10)/dx
                p22 = (3.0* (3.0*a-b-c)+this)/ (dx*dy)
                p23 = (-6.0*a+2.0*b+3.0*c-this)/ (dx*dysq)
                p30 = (-2.0*z0dx+zx10+zx00)/dxsq
                p31 = (-2.0*zy0dx+zxy10+zxy00)/dxsq
                p32 = (-6.0*a+3.0*b+2.0*c-this)/ (dxsq*dy)
                p33 = (2.0* (2.0*a-b-c)+this)/ (dxsq*dysq)
            END IF

            ! Evaluates the polynomial.
            u = xii - x0
            v = yii - y0
            q0 = p00 + v* (p01+v* (p02+v*p03))
            q1 = p10 + v* (p11+v* (p12+v*p13))
            q2 = p20 + v* (p21+v* (p22+v*p23))
            q3 = p30 + v* (p31+v* (p32+v*p33))
            zii = q0 + u* (q1+u* (q2+u*q3))
            ! End of Case 1

            ! Case 2.  When the rectangle is inside the data area in the x
            ! direction but outside in the y direction.
        ELSE IF ((ixdi > 0.AND.ixdi < nxd) .AND.  &
            (iydi <= 0.OR.iydi >= nyd)) THEN
        ! Retrieves the z and partial derivative values at the other
        ! vertex of the semi-infinite rectangle.
        IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
            ixd1 = ixd0 + 1
            dx = xd(ixd1) - x0
            dxsq = dx*dx
            z10 = zd(ixd1,iyd0)
            zx10 = pdd(1,ixd1,iyd0)
            zy10 = pdd(2,ixd1,iyd0)
            zxy10 = pdd(3,ixd1,iyd0)
            ! Calculates the polynomial coefficients.
            z0dx = (z10-z00)/dx
            zy0dx = (zy10-zy00)/dx
            p00 = z00
            p01 = zy00
            p10 = zx00
            p11 = zxy00
            p20 = (2.0* (z0dx-zx00)+z0dx-zx10)/dx
            p21 = (2.0* (zy0dx-zxy00)+zy0dx-zxy10)/dx
            p30 = (-2.0*z0dx+zx10+zx00)/dxsq
            p31 = (-2.0*zy0dx+zxy10+zxy00)/dxsq
        END IF
        ! Evaluates the polynomial.
        u = xii - x0
        v = yii - y0
        q0 = p00 + v*p01
        q1 = p10 + v*p11
        q2 = p20 + v*p21
        q3 = p30 + v*p31
        zii = q0 + u* (q1+u* (q2+u*q3))
        ! End of Case 2

        ! Case 3.  When the rectangle is outside the data area in the x
        ! direction but inside in the y direction.
        ELSE IF ((ixdi <= 0.OR.ixdi >= nxd) .AND.  &
            (iydi > 0 .AND. iydi < nyd)) THEN
        ! Retrieves the z and partial derivative values at the other
        ! vertex of the semi-infinite rectangle.
        IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
            iyd1 = iyd0 + 1
            dy = yd(iyd1) - y0
            dysq = dy*dy
            z01 = zd(ixd0,iyd1)
            zx01 = pdd(1,ixd0,iyd1)
            zy01 = pdd(2,ixd0,iyd1)
            zxy01 = pdd(3,ixd0,iyd1)
            ! Calculates the polynomial coefficients.
            z0dy = (z01-z00)/dy
            zx0dy = (zx01-zx00)/dy
            p00 = z00
            p01 = zy00
            p02 = (2.0*(z0dy-zy00)+z0dy-zy01)/dy
            p03 = (-2.0*z0dy+zy01+zy00)/dysq
            p10 = zx00
            p11 = zxy00
            p12 = (2.0*(zx0dy-zxy00) + zx0dy - zxy01)/dy
            p13 = (-2.0*zx0dy + zxy01 + zxy00)/dysq
        END IF

        ! Evaluates the polynomial.
        u = xii - x0
        v = yii - y0
        q0 = p00 + v* (p01 + v*(p02+v*p03))
        q1 = p10 + v* (p11 + v*(p12+v*p13))
        zii = q0 + u*q1
        ! End of Case 3

        ! Case 4.  When the rectangle is outside the data area in both the
        ! x and y direction.
        ELSE IF ((ixdi <= 0 .OR. ixdi >= nxd) .AND.  &
            (iydi <= 0 .OR. iydi >= nyd)) THEN
        ! Calculates the polynomial coefficients.
        IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
            p00 = z00
            p01 = zy00
            p10 = zx00
            p11 = zxy00
        END IF
        ! Evaluates the polynomial.
        u = xii - x0
        v = yii - y0
        q0 = p00 + v*p01
        q1 = p10 + v*p11
        zii = q0 + u*q1
        END IF
        ! End of Case 4
        zi(iip) = zii
    END DO

    RETURN
    END SUBROUTINE rgplnl

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
          !!$omp parallel &
          !!$omp default(shared) &
          !!$omp private(k1,wt,xx,x) &
          !!$omp reduction(+:RES)
!do k1=1,int(pts)                                        ! linea modificada por mi
   80    WT = DPK
         DO 120 J = 1, NDIM
            XX = MOD(ALPHA(J)+VK(J)*PK*DPK,1.0D0)
!            XX = MOD(ALPHA(J)+VK(J)*dble(k1)*DPK,1.0D0) !linea anterior modificada por mi: PK --> dble(k1)
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
!enddo                                                   ! linea modificada por mi
          !!$omp end parallel do
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
end module utilities

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
