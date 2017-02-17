module ndConfig
	use precision
	use constants
	use variables
	use ndErrors
	use FileIni
	use ndMatrices
	use ndInteractions
	implicit none

	integer :: num_args, ixa
	character(len=300), dimension(:), allocatable :: args
	
!	real(dl) :: minim_prec = 1d-4
!	integer :: itmax_root=5000
	
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
	
	subroutine allocateStuff()
		integer :: nf
		nf = flavorNumber
		allocate(nuMasses(nf))
		allocate(mixMat(nf,nf), mixMatInv(nf,nf))
		allocate(nuMassesMat(nf,nf), leptonDensities(nf,nf))
		allocate(GL_mat(nf,nf), GR_mat(nf,nf), GLR_vec(2, nf,nf))
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
		real(dl), dimension(:,:), allocatable :: m1
		
		nf = flavorNumber
		allocate(mv(nf), m1(nf,nf))
		
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
		call tripleProdMat(mixMat, m1, mixMatInv, nuMassesMat)
		call printMat(nuMassesMat)
	end subroutine setMassMatrix
	
	subroutine setMixingMatrix()
		integer :: nf
		real(dl), dimension(:,:), allocatable :: m1, m2, m3, m4
		character(len=1) :: nfstr
		integer :: ix
		real(dl) :: tmp
		
		nf = flavorNumber
		write(nfstr,"(I1)") nf
		
		if (verbose.gt.1) &
			call addToLog("[config] creating mixing matrix, using "//nfstr//" flavors")
		allocate(m1(nf,nf), m2(nf,nf), m3(nf,nf), m4(nf,nf))
		
		if (nf.gt.3) then
			call createRotMat(m1, nf, 3, 4, theta34)
			call createRotMat(m2, nf, 2, 4, theta24)
			call createRotMat(m3, nf, 1, 4, theta14)
			call tripleProdMat(m1, m2, m3, m4)
			m1=m4
		else
			call createIdentityMat(m1, nf)
		end if
		if (nf.gt.2) then
			call createRotMat(m2, nf, 2, 3, theta23)
			call createRotMat(m3, nf, 1, 3, theta13)
			call tripleProdMat(m1, m2, m3, m4)
			m1=m4
		else
			call createRotMat(m1, nf, 1, 2, 0.0d0)
		end if
		
		call createRotMat(m2, nf, 1, 2, theta12)
		call multiplyMat(m1, m2, m3)
		
		mixMat=m3
		call inverseMat(mixMat, mixMatInv)
	end subroutine setMixingMatrix
	
	subroutine init_matrices
		real(dl), dimension(:), allocatable :: diag_el
		integer :: ix
		real(dl) :: fdm
		
		allocate(diag_el(flavorNumber))
		
		!GL
		diag_el = sin2thW - 0.5d0
		diag_el(1) = sin2thW + 0.5d0
		call createDiagMat(GL_mat, flavorNumber, diag_el)
		GLR_vec(1,:,:) = GL_mat
		
		!GR
		diag_el = sin2thW
		call createDiagMat(GR_mat, flavorNumber, diag_el)
		GLR_vec(2,:,:) = GR_mat
		
		!lepton densities
		leptonDensities = 0.d0
!		leptonDensities(1,1) = ...

		!identity (Matrix complex)
		allocate(idMat(flavorNumber,flavorNumber))
		call createIdentityMat(idMat, flavorNumber)
		
		!nu density matrix
		allocate(nuDensMatVec(Ny))
		open(unit=3154, file='test_fd.dat') !save the initial Fermi-Dirac distribution for nu
		do ix=1, Ny
			allocate(nuDensMatVec(ix)%re(flavorNumber,flavorNumber), nuDensMatVec(ix)%im(flavorNumber,flavorNumber))
			nuDensMatVec(ix)%y = y_arr(ix)
			nuDensMatVec(ix)%re(:,:) = 0.d0
			nuDensMatVec(ix)%im(:,:) = 0.d0
			fdm = fermiDirac_massless(y_arr(ix),z_in)
			write(3154,"(2F14.7)") y_arr(ix), fdm
			nuDensMatVec(ix)%re(1,1) = fdm
			if (flavorNumber.ne.2 .or. (.not. only_1a_1s)) then
				nuDensMatVec(ix)%re(2,2) = fdm
				if (flavorNumber.gt.2) &
					nuDensMatVec(ix)%re(3,3) = fdm
				if (flavorNumber.gt.3) &
					nuDensMatVec(ix)%re(4,4) = 0.d0
			end if
		end do
		close(3154)
		
		deallocate(diag_el)
	end subroutine init_matrices
	
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
			call addToLog("[config] reading additional configuration from "//mainPath//trim(args(1)))
			call ini_file_open(mainPath//trim(args(1)), mainPath//trim(args(1))//".log")
			
			verbose = read_ini_int('verbose',verbose)
			maxiter = read_ini_int('maxiter',100)
			toler   = read_ini_real('tolerance', 1.d-3)
			
			Nx = read_ini_int('Nx',100000)
			Ny = read_ini_int('Ny',100)
			allocate(x_arr(Nx), y_arr(Ny))
			
			x_in    = read_ini_real('x_in', 0.01d0)
			x_fin   = read_ini_real('x_fin', 40.d0)
			x_arr = logspace(log10(x_in), log10(x_fin), Nx)

			y_min   = read_ini_real('y_min', 0.0d0)
			y_max   = read_ini_real('y_max', 20.0d0)
			y_arr = logspace(log10(y_min), log10(y_max), Ny)

			z_in    = read_ini_real('z_in', 1.00003d0)
			
			!read mixing parameters and create matrices
			m_lightest   = read_ini_real('m_lightest', 0.0d0)
			massOrdering = read_ini_logical('massOrdering', .true.)
			
			only_1a_1s = read_ini_logical('only_1a_1s', .false.)
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
				call error("[config] WARNING: only up to 4 neutrino flavors are supported. Using N=4")
				flavorNumber = 4
			end if
			
			call setMixingMatrix()
			call setMassMatrix()
			
			!read cosmo parameters
			hubbleParam = read_ini_real('hubbleParam', i_HubbleParam)
			photonTemperatureToday = read_ini_real('photonTemperatureToday', i_photonTempToday)
			
			!create other matrices
			call init_matrices
			
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
