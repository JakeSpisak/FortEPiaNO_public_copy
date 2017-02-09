module ndConfig
	use precision
	use FileIni
	implicit none

	integer :: num_args, ixa
	character(len=300), dimension(:), allocatable :: args
	
!	real(dl) :: minim_prec = 1d-4
!	integer :: itmax_root=5000
	integer :: verbose = 1
	
	character(LEN=*), parameter :: mainPath="/home/gariazzo/software/nuDensity/"

!	integer, parameter :: numNuisance = 6
!	type nuisanceParam
!		character(LEN=15) :: paramName
!		real(dl) :: defaultVal
!		logical :: fixed, setPrior
!		real(dl) :: mean, devst
!	end type nuisanceParam
!	type (nuisanceParam), dimension(numNuisance) :: nuisanceParams
!	integer :: varyingNuisance
!	integer :: atmFluxNormIX, pi2kaRatioIX, nu2anuRatioIX, CRSpectralIndexIX, deltaIX, medianEnergyIX
	
!	integer, parameter :: EMuRecoBins=10 !10 for the 2D analysis, 28 for the 1D (not used)
!	integer, parameter :: cosThetaBins=21
!	real(dl), dimension(2) :: EMuRecoInterval  = [400.,20000.]
!	real(dl), dimension(2) :: cosThetaInterval = [-1., 0.2]
!	real(dl), dimension(EMuRecoBins) :: EMuRecoEdges
!	real(dl), dimension(cosThetaBins) :: cosThetaEdges
	
	contains
	
!	subroutine setParam(param, pname, defaultval, fixed, prior, mean, devst)
!		integer, intent(in) :: param
!		character(LEN=*), intent(in) :: pname
!		logical, intent(in) :: fixed, prior
!		real(dl), intent(in) :: mean, devst, defaultval
		
!		nuisanceParams(param)%paramName = pname
!		nuisanceParams(param)%fixed = fixed
!		if (.not. fixed) varyingNuisance = varyingNuisance + 1
!		nuisanceParams(param)%setPrior = prior
!		nuisanceParams(param)%defaultVal = defaultval
!		nuisanceParams(param)%mean = mean
!		nuisanceParams(param)%devst = devst
!	end subroutine setParam
	
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
	
	subroutine initConfig()
		character(len=300) :: tmparg
		
		if (verbose>0) print *, '[config] init configuration'
		num_args = command_argument_count()
		
!		atmFluxNormIX    = 1
!		pi2kaRatioIX     = 2
!		nu2anuRatioIX    = 3
!		CRSpectralIndexIX= 4
!		deltaIX          = 5
!		medianEnergyIX   = 6
!		call setParam(atmFluxNormIX,     "atmFluxNorm",     1.0d0, .false., .true., 1.d0, 0.4d0)
!		call setParam(pi2kaRatioIX,      "pi2kaRatio",      1.0d0, .false., .true., 1.d0, 0.1d0)
!		call setParam(nu2anuRatioIX,     "nu2anuRatio",     1.0d0, .false., .true., 1.d0, 2.5d-2)
!		call setParam(CRSpectralIndexIX, "CRSpectralIndex", 0.0d0, .false., .true., 0.d0, 5.d-2)
!		call setParam(deltaIX,           "delta",           0.0d0, .false., .true., 0.d0, 5.d-2)
!		call setParam(medianEnergyIX,    "medianEnergy",    2.0d3, .true., .false., 0.d0, 0.d0) !in GeV
		
!		EMuRecoEdges = logspace(log10(EMuRecoInterval(1)), log10(EMuRecoInterval(2)), EMuRecoBins)
!		cosThetaEdges= linspace(cosThetaInterval(1),       cosThetaInterval(2),       cosThetaBins)
		
		allocate(args(num_args))
		do ixa = 1, num_args
			call get_command_argument(ixa,args(ixa))
		end do
		if(num_args.gt.0) then
			print *, "[config] reading additional configuration from "//mainPath//trim(args(1))
!			call ini_file_open(mainPath//trim(args(1)), mainPath//"log/"//trim(args(1)))
!			tmparg=trim(read_ini_char('pionFluxFile'))
!			if (trim(tmparg)/="") pionFluxFile=tmparg
!			tmparg=trim(read_ini_char('kaonFluxFile'))
!			if (trim(tmparg)/="") kaonFluxFile=tmparg
!			tmparg=trim(read_ini_char('outMinInfo'))
!			if (trim(tmparg)/="") outMinInfo=tmparg
!			tmparg=trim(read_ini_char('outPlotNames'))
!			if (trim(tmparg)/="") outPlotNames=tmparg
			
			verbose = read_ini_int('verbose',verbose)
		end if
		call ini_file_close()
		
	end subroutine initconfig
end module ndConfig
