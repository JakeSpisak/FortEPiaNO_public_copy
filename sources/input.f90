module fpinput
	use precision
	use variables
	use fpErrors
	use FileIni
	implicit none

	contains

	subroutine readBasicParams
		character(len=300) :: tmparg

		flavorNumber = read_ini_int('flavorNumber', i_flavorNumber)
		force_replace = read_ini_logical('force_replace', .false.)
		tmparg=trim(read_ini_char('outputFolder'))
		if (trim(tmparg)/="") outputFolder=tmparg
		num_threads = read_ini_int('num_threads', 0)
	end subroutine readBasicParams

	subroutine readNuDef
		integer :: ix
		character(len=300) :: tmparg

		do ix=1, flavorNumber
			write(tmparg,"('nuFactor',I1)") ix
			nuFactor(ix) = read_ini_real(trim(tmparg), 1.d0)
			write(tmparg,"('sterile',I1)") ix
			if (ix.le.3) then
				sterile(ix) = read_ini_logical(trim(tmparg), .false.)
			else
				sterile(ix) = read_ini_logical(trim(tmparg), .true.)
			end if
		end do
	end subroutine readNuDef

	subroutine readNuMixing
		integer :: ix, iy
		character(len=300) :: tmpstr

		giveSinSq = read_ini_logical('givesinsq', .true.)
		massSplittings = 0.d0
		massSplittings(2) = read_ini_real('dm21', i_dm21)
		if (flavorNumber .gt. 2) then
			massSplittings(3) = read_ini_real('dm31', i_dm31)
		end if
		do ix = 4, flavorNumber
			write(tmpstr,"('dm',I1,'1')") ix
			massSplittings(ix) = read_ini_real(trim(tmpstr), zero)
		end do
		if (giveSinSq) then
			mixingAngles(1,2) = asin(sqrt(read_ini_real('theta12', i_theta12)))
			if (flavorNumber .gt. 2) then
				mixingAngles(1,3) = asin(sqrt(read_ini_real('theta13', i_theta13)))
				mixingAngles(2,3) = asin(sqrt(read_ini_real('theta23', i_theta23)))
			end if
			do ix=4, flavorNumber
				do iy=1, ix-1
					write(tmpstr,"('theta',2I1)") iy, ix
					mixingAngles(iy, ix) = asin(sqrt(read_ini_real(trim(tmpstr), zero)))
				end do
			end do
		else
			mixingAngles(1,2) = read_ini_real('theta12', asin(sqrt(i_theta12)))
			if (flavorNumber .gt. 2) then
				mixingAngles(1,3) = read_ini_real('theta13', asin(sqrt(i_theta13)))
				mixingAngles(2,3) = read_ini_real('theta23', asin(sqrt(i_theta23)))
			end if
			do ix=4, flavorNumber
				do iy=1, ix-1
					write(tmpstr,"('theta',2I1)") iy, ix
					mixingAngles(iy, ix) = read_ini_real(trim(tmpstr), zero)
				end do
			end do
		end if
	end subroutine readNuMixing

#ifdef LOW_REHEATING
	subroutine readLowReheating
		Trh = read_ini_real('Trh', 0.d0)
	end subroutine readLowReheating
#endif

	subroutine readCollint
		integer :: ix, iy
		character(len=300) :: tmpstr

		collint_damping_type = read_ini_int('collint_damping_type', 1)
		collint_diagonal_zero = read_ini_logical('collint_diagonal_zero', .false.)
		collint_offdiag_damping = read_ini_logical('collint_offdiag_damping', .true.)
		collint_d_no_nue = read_ini_logical('collint_d_no_nue', .false.)
		collint_d_no_nunu = read_ini_logical('collint_d_no_nunu', .false.)
		collint_od_no_nue = read_ini_logical('collint_od_no_nue', .false.)
		collint_od_no_nunu = read_ini_logical('collint_od_no_nunu', .false.)
		damping_read_zero = .true.
		do ix=1, flavorNumber
			do iy=ix+1, flavorNumber
				write(tmpstr,"('damping_',2I1,'_zero')") ix, iy
				dampingSetZero(ix,iy) = read_ini_logical(trim(tmpstr), .false.)
			end do
		end do
	end subroutine readCollint

	subroutine readFTQED
		ftqed_temperature_corr = read_ini_logical('ftqed_temperature_corr',.true.)
		ftqed_log_term = read_ini_logical('ftqed_log_term',.false.)
		ftqed_ord3 = read_ini_logical('ftqed_ord3',.true.)
		ftqed_e_mth_leptondens = read_ini_logical('ftqed_e_mth_leptondens',.true.)
	end subroutine readFTQED

	subroutine readXY
		x_in  = read_ini_real('x_in', 0.01d0)
		x_fin = read_ini_real('x_fin', 40.d0)
		Nx = read_ini_int('Nx',100)
		use_gauss_laguerre = read_ini_logical('use_gauss_laguerre', .true.)
		Ny = read_ini_int('Ny',30)
		y_min = read_ini_real('y_min', 0.01d0)
		y_max = read_ini_real('y_max', 20.0d0)
	end subroutine readXY

	subroutine readOutputConfig
		checkpoint = read_ini_logical('checkpoint', .true.)
		save_BBN = read_ini_logical('save_BBN', .true.)
		save_fd = read_ini_logical('save_fd', .true.)
		save_energy_entropy_evolution = read_ini_logical('save_energy_entropy_evolution', .true.)
		save_Neff = read_ini_logical('save_Neff', .true.)
		save_nuDens_evolution = read_ini_logical('save_nuDens_evolution', .true.)
		save_number_evolution = read_ini_logical('save_number_evolution', .true.)
		save_w_evolution = read_ini_logical('save_w_evolution', .true.)
		save_z_evolution = read_ini_logical('save_z_evolution', .true.)
		intermediateSteps%output = read_ini_logical('save_intermediate_steps', .false.)
	end subroutine readOutputConfig

	subroutine readPrecision
		maxiter = read_ini_int('maxiter', 100)
		toler_jkyg = read_ini_real('tolerance_jkyg', 1.d-7)
		dlsoda_atol_z = read_ini_real('dlsoda_atol_z', 1.d-6)
		dlsoda_atol_d = read_ini_real('dlsoda_atol_d', 1.d-6)
		dlsoda_atol_o = read_ini_real('dlsoda_atol_o', 1.d-6)
		dlsoda_rtol = read_ini_real('dlsoda_rtol', 1.d-6)

		interp_nx = read_ini_int('interp_nx', interp_nx0)
		interp_nz = read_ini_int('interp_nz', interp_nz0)
		interp_nxz = read_ini_int('interp_nxz', interp_nxz0)
		interp_zmin = read_ini_real('interp_zmin', interp_zmin0)
		interp_zmax = read_ini_real('interp_zmax', interp_zmax0)
	end subroutine readPrecision

	subroutine readVerbosity
		verbose = read_ini_int('verbose', verbose)
		Nprintderivs = read_ini_real('Nprintderivs', Nprintderivs)
	end subroutine readVerbosity
end module fpinput
