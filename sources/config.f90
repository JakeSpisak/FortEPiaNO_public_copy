module ndConfig
	use precision
	use FileIni
	use constants
	use ndMatrices
	implicit none

	integer :: num_args, ixa
	character(len=300), dimension(:), allocatable :: args
	
!	real(dl) :: minim_prec = 1d-4
!	integer :: itmax_root=5000
	integer :: verbose = 1
	
	character(LEN=*), parameter :: mainPath="/home/gariazzo/software/nuDensity/"

	logical :: massOrdering
	integer :: flavorNumber
	real(dl) :: m_lightest
	real(dl) :: theta12, dm12
	real(dl) :: theta13, theta23, dm23, deltaCP13
	real(dl) :: theta14, theta24, theta34, dm14
	real(dl), dimension(:), allocatable :: nuMasses
	real(dl), dimension(:,:), allocatable :: mixMat, mixMatInv, nuMassesMat, leptonDensities
	
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
	
	subroutine allocateStuff()
		integer :: nf
		nf = flavorNumber
		allocate(nuMasses(nf))
		allocate(mixMat(nf,nf), mixMatInv(nf,nf))
		allocate(nuMassesMat(nf,nf), leptonDensities(nf,nf))
	end subroutine allocateStuff
	
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
	
	subroutine setMassMatrix()
		integer :: nf
		integer :: ix
		real(dl), dimension(:), allocatable :: mv
		real(dl), dimension(:,:), allocatable :: m1, m2
		
		nf = flavorNumber
		allocate(mv(nf), m1(nf,nf), m2(nf,nf))
		
		if (nf .eq. 2) then
			mv(1) = m_lightest
			mv(2) = sqrt(m_lightest*m_lightest + dm12)
		end if
		
		if (massOrdering) then
			mv(1) = m_lightest
			mv(2) = sqrt(m_lightest*m_lightest + dm12)
			mv(3) = sqrt(mv(2)*mv(2) + dm23)
		else
			mv(3) = m_lightest
			mv(2) = sqrt(m_lightest*m_lightest + dm23)
			mv(1) = sqrt(mv(2)*mv(2) - dm12)
		end if
		if (nf.gt.3) then
			mv(4) = sqrt(mv(1)*mv(1) + dm14)
		end if
		
		call createDiagMat(m1, nf, mv)
		call multMat(mixMat, m1, m2)
		call multMat(m2, mixMatInv, nuMassesMat)
		call printMat(nuMassesMat)
	end subroutine setMassMatrix
	
	subroutine setMixingMatrix()
		integer :: nf
		real(dl), dimension(:,:), allocatable :: m1, m2, m3
		integer :: ix
		real(dl) :: tmp
		
		nf = flavorNumber
		
		if (verbose.gt.1) write (*,"(' [config] creating mixing matrix, using ',I2,' flavors')") nf
		allocate(m1(nf,nf), m2(nf,nf), m3(nf,nf))
		
		if (nf.gt.3) then
			call createRotMat(m1, nf, 3, 4, theta34)
			call createRotMat(m2, nf, 2, 4, theta24)
			call multMat(m1, m2, m3)
			m1=m3
			call createRotMat(m2, nf, 1, 4, theta14)
			call multMat(m1, m2, m3)
			m1=m3
		else
			call createRotMat(m1, nf, 1, 2, 0.0d0)
		end if
		if (nf.gt.2) then
			call createRotMat(m2, nf, 2, 3, theta23)
			call multMat(m1, m2, m3)
			m1=m3
			call createRotMat(m2, nf, 1, 3, theta13)
			call multMat(m1, m2, m3)
			m1=m3
		else
			call createRotMat(m1, nf, 1, 2, 0.0d0)
		end if
		
		call createRotMat(m2, nf, 1, 2, theta12)
		call multMat(m1, m2, m3)
		
		mixMat=m3
		call inverseMat(mixMat, mixMatInv)
	end subroutine setMixingMatrix
	
	subroutine initConfig()
		character(len=300) :: tmparg
		
		if (verbose>0) write(*,*) '[config] init configuration'
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
			write (*,*) "[config] reading additional configuration from "//mainPath//trim(args(1))
			call ini_file_open(mainPath//trim(args(1)), mainPath//trim(args(1))//".log")
			verbose = read_ini_int('verbose',verbose)
			
			m_lightest   = read_ini_real('m_lightest', 0.0d0)
			massOrdering = read_ini_logical('massOrdering', .true.)
			
			flavorNumber = read_ini_int('flavorNumber', i_flavorNumber)
			call allocateStuff
			
			theta12      = read_ini_real('theta12', i_theta12)
			dm12         = read_ini_real('dm12', i_dm12)
			if (flavorNumber .gt. 2) then
				theta13      = read_ini_real('theta13', i_theta13)
				theta23      = read_ini_real('theta23', i_theta23)
				dm23         = read_ini_real('dm23', i_dm23)
				deltaCP13    = read_ini_real('deltaCP13', i_deltaCP13)
			end if
			if (flavorNumber .gt. 3) then
				theta14      = read_ini_real('theta14', zero)
				theta24      = read_ini_real('theta24', zero)
				theta24      = read_ini_real('theta34', zero)
				dm14         = read_ini_real('dm14', zero)
			end if
			if (flavorNumber .gt. 4) then
				write(*,*) "[config] WARNING: only up to 4 neutrino flavors are supported. Using N=4"
				flavorNumber = 4
			end if
			
			call setMixingMatrix()
			call setMassMatrix()
			
!			tmparg=trim(read_ini_char('pionFluxFile'))
!			if (trim(tmparg)/="") pionFluxFile=tmparg
!			tmparg=trim(read_ini_char('kaonFluxFile'))
!			if (trim(tmparg)/="") kaonFluxFile=tmparg
!			tmparg=trim(read_ini_char('outMinInfo'))
!			if (trim(tmparg)/="") outMinInfo=tmparg
!			tmparg=trim(read_ini_char('outPlotNames'))
!			if (trim(tmparg)/="") outPlotNames=tmparg
			
		end if
		call ini_file_close()
		
	end subroutine initconfig
end module ndConfig
