module fpOutput
	use precision
	use variables
	use fpInterfaces1
	use fpErrors
	use utilities
	use fpCosmology
	use fpMatter
	implicit none

	real(dl) :: deriv_counter

	integer :: Nsave
	logical :: first_store = .true.
	real(dl), dimension(:), allocatable :: storeX, storeNorm, storeTmp
	real(dl), dimension(:,:), allocatable :: storeY, storeYdot, storeHeff, storeComm, storeCTNue, storeCTNunu
	character(len=14), parameter :: intermfmt = "(1P, *(E14.6))"
	integer, parameter :: iu = 8972

	contains

	subroutine nuDens_to_file(u, ix, iy, x, mat, reim, fname)
		integer, intent(in) :: u, ix, iy
		real(dl), intent(in) :: x
		logical, intent(in) :: reim!true for real, false for imaginary part
		type(cmplxMatNN), dimension(:), allocatable, intent(in) :: mat
		character(len=*), intent(in) :: fname
		integer :: m
		real(dl), dimension(:), allocatable :: tmpvec

		allocate(tmpvec(Ny))
		call openFile(u, trim(fname), firstWrite)
		if (reim) then
			do m=1, nY
				tmpvec(m)=mat(m)%re(ix, iy)
			end do
		else
			do m=1, nY
				tmpvec(m)=mat(m)%im(ix, iy)
			end do
		end if
		write(u, multidblfmt) x, tmpvec
		deallocate(tmpvec)
		close(u)
	end subroutine nuDens_to_file

	subroutine saveRelevantInfo(x, vec)
		real(dl), intent(in) :: x
		real(dl), dimension(:), intent(in) :: vec
		type(cmplxMatNN), dimension(:), allocatable :: rho_mass
		complex(dl), dimension(:,:), allocatable :: tmpComplMat, transfMat
		real(dl), dimension(maxFlavorNumber) :: nuEnDens
		integer :: k, i, j, iy
		real(dl) :: neff, z, w
		character(len=200) :: fname
		procedure (nuDensity_integrator), pointer :: nuDensityInt, nuNumDensInt

		if (use_gauss_laguerre) then
			nuDensityInt => nuDensityGL
			nuNumDensInt => nuNumberDensityGL
		else
			nuDensityInt => nuDensityNC
			nuNumDensInt => nuNumberDensityNC
		end if

		write(fname, '(A,'//dblfmt//')') '[output] Saving info at x=', x
		call addToLog(trim(fname))!not a filename but the above string

		w = vec(ntot-1)+1.d0
		z = vec(ntot)+1.d0
		if (save_nuDens_evolution) then
			!density matrix in flavor space
			do k=1, flavorNumber
				write(fname, '(A,I1,A)') trim(outputFolder)//'/nuDens_diag', k, '.dat'
				call nuDens_to_file(iu, k, k, x, nuDensMatVecFD, .true., trim(fname))
			end do
			if (has_offdiagonal()) then
				do i=1, flavorNumber-1
					do j=i+1, flavorNumber
						write(fname, '(A,I1,I1,A)') trim(outputFolder)//'/nuDens_nd_', i, j, '_re.dat'
						call nuDens_to_file(iu, i, j, x, nuDensMatVecFD, .true., trim(fname))

						write(fname, '(A,I1,I1,A)') trim(outputFolder)//'/nuDens_nd_', i, j, '_im.dat'
						call nuDens_to_file(iu, i, j, x, nuDensMatVecFD, .false., trim(fname))
					end do
				end do
			end if
			!density matrix in mass space
			allocate(rho_mass(Ny))
			call updateMatterDensities(x, z)
			!$omp parallel do shared(rho_mass) private(iy) schedule(static)
			do iy=1, Ny
				rho_mass(iy) = rho_diag_mass(iy)
			end do
			!$omp end parallel do
			do k=1, flavorNumber
				write(fname, '(A,I1,A)') trim(outputFolder)//'/nuDens_mass', k, '.dat'
				call nuDens_to_file(iu, k, k, x, rho_mass, .true., trim(fname))
			end do
			deallocate(rho_mass)
		end if
		if (save_energy_entropy_evolution) then
			do k=1, flavorNumber
				nuEnDens(k) = nuDensityInt(k, k)*nuFactor(k)
			end do
			call openFile(iu, trim(outputFolder)//'/energyDensity.dat', firstWrite)
			write(iu, multidblfmt) x, z, &
				photonDensity(z), &
				electrons%energyDensity(x, z, .false.), &
#ifdef DO_MUONS
				muons%energyDensity(x, z, .false.), &
#endif
				nuEnDens(1:flavorNumber)
			close(iu)
			call openFile(iu, trim(outputFolder)//'/entropy.dat', firstWrite)
			write(iu, multidblfmt) x, z, &
				photonEntropy(z), &
				electrons%entropy(x, z), &
#ifdef DO_MUONS
				muons%entropy(x, z), &
#endif
				nuEnDens(1:flavorNumber)*four_thirds/w
			close(iu)
		end if
		if (save_number_evolution) then
			do k=1, flavorNumber
				nuEnDens(k) = nuNumDensInt(k, k)*nuFactor(k)
			end do
			call openFile(iu, trim(outputFolder)//'/numberDensity.dat', firstWrite)
			write(iu, multidblfmt) x, z, &
				photonNumberDensity(z), &
				electrons%numberDensity(x, z, .false.), &
#ifdef DO_MUONS
				muons%numberDensity(x, z, .false.), &
#endif
				nuEnDens(1:flavorNumber)
			close(iu)
		end if
		if (save_z_evolution) then
			call openFile(iu, trim(outputFolder)//'/z.dat', firstWrite)
			if (save_w_evolution) then
				write(iu, multidblfmt) x, z, vec(ntot-1)+1.d0
			else
				write(iu, multidblfmt) x, z
			end if
			close(iu)
		end if
		if (save_Neff) then
			neff = Neff_from_rho_z(vec(ntot)+1.d0)
			call openFile(iu, trim(outputFolder)//'/Neff.dat', firstWrite)
			write(iu, multidblfmt) &
				x, neff/zid**4, neff
			close(iu)
		end if
		if (firstWrite) &
			call ens_header
		firstWrite=.false.
	end subroutine saveRelevantInfo

	subroutine ens_header
		integer :: k
		character(len=200) :: tmpstr

		if (save_energy_entropy_evolution .or. save_number_evolution) then
			call openFile(iu, trim(outputFolder)//'/ens_header.dat', firstWrite)
			write(tmpstr, "('nu (1 to ', I1, ')')") flavorNumber
			write(iu, "(*(A17))") "x", "z", &
				"photon", &
				"electron", &
#ifdef DO_MUONS
				"muon", &
#endif
				trim(tmpstr)
			close(iu)
		end if
	end subroutine

	subroutine allocateStoreVar1D(vec, N)
		real(dl), dimension(:), allocatable :: vec
		integer :: N
		if (.not.allocated(vec)) &
			allocate(vec(N))
	end subroutine allocateStoreVar1D

	subroutine allocateStoreVar2D(mat, Nx, Nv)
		real(dl), dimension(:,:), allocatable :: mat
		integer :: Nx, Nv
		if (.not.allocated(mat)) &
			allocate(mat(Nx, Nv))
	end subroutine allocateStoreVar2D

	subroutine allocateStoreVars
		if (Nprintderivs.ge.100) then
			Nsave = Nprintderivs
		else
			Nsave = 100
		end if
		call allocateStoreVar1D(storeX, Nsave)
		call allocateStoreVar1D(storeNorm, Nsave)
		call allocateStoreVar1D(storeTmp, ntotrho)
		call allocateStoreVar2D(storeY, Nsave, ntot)
		call allocateStoreVar2D(storeYdot, Nsave, ntot)
		call allocateStoreVar2D(storeHeff, Nsave, ntotrho)
		call allocateStoreVar2D(storeComm, Nsave, ntotrho)
		call allocateStoreVar2D(storeCTNue, Nsave, ntotrho)
		call allocateStoreVar2D(storeCTNunu, Nsave, ntotrho)
	end subroutine allocateStoreVars

	subroutine deallocateStoreVars
		if (intermediateSteps%output) &
			call intermediateToFiles(int(mod(deriv_counter, 1.d0*Nsave)))

		if(allocated(storeX))&
			deallocate(storeX)
		if(allocated(storeNorm))&
			deallocate(storeNorm)
		if(allocated(storeTmp))&
			deallocate(storeTmp)
		if(allocated(storeY))&
			deallocate(storeY)
		if(allocated(storeYdot))&
			deallocate(storeYdot)
		if(allocated(storeHeff))&
			deallocate(storeHeff)
		if(allocated(storeComm))&
			deallocate(storeComm)
		if(allocated(storeCTNue))&
			deallocate(storeCTNue)
		if(allocated(storeCTNunu))&
			deallocate(storeCTNunu)
	end subroutine deallocateStoreVars

	subroutine intermediateToFiles(N)
		integer :: i, N
		call openFile(iu, trim(outputFolder)//'/intermXF.dat', first_store)
		do i=1, N
			write(iu, intermfmt) storeX(i), storeNorm(i)
		end do
		close(iu)
		call openFile(iu, trim(outputFolder)//'/intermY.dat', first_store)
		do i=1, N
			write(iu, intermfmt) storeY(i, :)
		end do
		close(iu)
		call openFile(iu, trim(outputFolder)//'/intermYdot.dat', first_store)
		do i=1, N
			write(iu, intermfmt) storeYdot(i, :)
		end do
		close(iu)
		call openFile(iu, trim(outputFolder)//'/intermHeff.dat', first_store)
		do i=1, N
			write(iu, intermfmt) storeHeff(i, :)
		end do
		close(iu)
		call openFile(iu, trim(outputFolder)//'/intermComm.dat', first_store)
		do i=1, N
			write(iu, intermfmt) storeComm(i, :)
		end do
		close(iu)
		call openFile(iu, trim(outputFolder)//'/intermCTNue.dat', first_store)
		do i=1, N
			write(iu, intermfmt) storeCTNue(i, :)
		end do
		close(iu)
		call openFile(iu, trim(outputFolder)//'/intermCTNunu.dat', first_store)
		do i=1, N
			write(iu, intermfmt) storeCTNunu(i, :)
		end do
		close(iu)
		first_store=.false.
	end subroutine intermediateToFiles

	subroutine saveIntermediateSteps
		integer :: cix

		if (.not.intermediateSteps%output) &
			return

		cix=mod(deriv_counter-1, 1.d0*Nsave)+1

		!store things
		storeX(cix) = intermediateSteps%x
		storeNorm(cix) = intermediateSteps%norm
		storeY(cix, :) = intermediateSteps%yvec
		storeYdot(cix, :) = intermediateSteps%ydot
		call mat_2_vec(intermediateSteps%Heff, Ny, storeTmp)
		storeHeff(cix,:) = storeTmp
		call mat_2_vec(intermediateSteps%commutator, Ny, storeTmp)
		storeComm(cix,:) = storetmp
		call mat_2_vec(intermediateSteps%colltermsNue(:), Ny, storeTmp)
		storeCTNue(cix,:) = storeTmp
		call mat_2_vec(intermediateSteps%colltermsNunu(:), Ny, storeTmp)
		storeCTNunu(cix,:) = storeTmp

		if (cix.eq.Nsave) then
			call intermediateToFiles(Nsave)
		end if
	end subroutine saveIntermediateSteps

	subroutine finalresults
		procedure (nuDensity_integrator), pointer :: nuDensityInt
		real(dl) :: ndeq, tmp, w, z
		real(dl) :: totrhonu, Neff
		real(dl), dimension(:), allocatable :: tmpvec
		integer :: ix, iy
		type(cmplxMatNN), dimension(:), allocatable :: rho_mass

		if (use_gauss_laguerre) then
			nuDensityInt => nuDensityGL
		else
			nuDensityInt => nuDensityNC
		end if

		w = nuDensVec(ntot-1) + 1.d0
		z = nuDensVec(ntot) + 1.d0

		!save final diagonal elements of the neutrino density matrix, in flavor basis
		call openFile(9876, trim(outputFolder)//'/rho_final.dat', .true.)
		allocate(tmpvec(flavorNumber))
		do iy=1, nY
			do ix=1, flavorNumber
				tmpvec(ix)=nuDensMatVecFD(iy)%re(ix, ix)
			end do
			write(9876, multidblfmt) nuDensMatVecFD(iy)%y, tmpvec
		end do
		close(9876)
		!save final diagonal elements of the neutrino density matrix, in mass basis
		allocate(rho_mass(Ny))
		call updateMatterDensities(x_arr(Nx), z)
		!$omp parallel do shared(rho_mass) private(iy) schedule(static)
		do iy=1, Ny
			rho_mass(iy) = rho_diag_mass(iy)
		end do
		!$omp end parallel do
		call openFile(9876, trim(outputFolder)//'/rho_final_mass.dat', .true.)
		do iy=1, nY
			do ix=1, flavorNumber
				tmpvec(ix)=rho_mass(iy)%re(ix, ix)
			end do
			write(9876, multidblfmt) nuDensMatVecFD(iy)%y, tmpvec
		end do
		close(9876)
		deallocate(rho_mass)
		deallocate(tmpvec)

		call openFile(9876, trim(outputFolder)//'/resume.dat', .true.)
		if (save_w_evolution) then
			write(*,"('final w = ',F11.8)") w
			write(9876,"('final w = ',F11.8)") w
		end if
		write(*,"('final z = ',F11.8)") z
		write(9876,"('final z = ',F11.8)") z
		totrhonu = allNuDensity()
		Neff = Neff_from_rho_z(z)
		do ix=1, flavorNumber
			tmp = nuDensityInt(ix, ix)*nuFactor(ix)/totrhonu * Neff
			write(*,"('deltaNeff_',I1,'  = ',F9.6)") ix, tmp
			write(9876,"('deltaNeff_',I1,'  = ',F9.6)") ix, tmp
		end do
		write(*,"('Neff    = ',F9.6)") Neff
		write(9876,"('Neff    = ',F9.6)") Neff

		close(9876)
	end subroutine finalresults

end module fpOutput
