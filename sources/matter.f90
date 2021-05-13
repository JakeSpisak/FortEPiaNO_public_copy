module fpMatter
	use precision
	use constants
	use variables
	use fpInterfaces1
	use fpCosmology
	use diagonalize
	implicit none

	contains

	subroutine updateMatterDensities(x, z)
		real(dl), intent(in) :: x, z
		real(dl) :: ldf
		integer :: ix, iy
		procedure (nuDensity_integrator), pointer :: nuDensityInt

		if (use_gauss_laguerre) then
			nuDensityInt => nuDensityGL
		else
			nuDensityInt => nuDensityNC
		end if

		leptonDensities = 0.d0
		ldf = leptDensFactor / x**6
		leptonDensities(1,1) = ldf * ( &
			electrons%energyDensity(x, z, ftqed_e_mth_leptondens) &
			+ electrons%pressure(x, z, ftqed_e_mth_leptondens) &
		)
#ifndef NO_MUONS
		if (flavorNumber.gt.2) &
			leptonDensities(2,2) = ldf * ( &
				muons%energyDensity(x, z, .false.) &
				+ muons%pressure(x, z, .false.) &
			)
#endif

		nuDensities%re = 0.d0
		nuDensities%im = 0.d0
		do ix=1, flavorNumber
			if (.not.sterile(ix)) then !this implements the Gs matrix, without computing the integral
				nuDensities%re(ix, ix) = nuDensities%re(ix, ix) + nuDensityInt(ix, ix)
				do iy=ix+1, flavorNumber
					if (.not.sterile(iy)) then !this implements the Gs matrix, without computing the integral
						nuDensities%re(ix, iy) = nuDensities%re(ix, iy) + nuDensityInt(ix, iy)
						nuDensities%im(ix, iy) = nuDensities%im(ix, iy) + nuDensityInt(ix, iy, .false.)
					end if
				end do
			end if
			do iy=ix+1, flavorNumber
				nuDensities%re(iy, ix) = nuDensities%re(ix, iy)
				nuDensities%im(iy, ix) = - nuDensities%im(ix, iy)
			end do
		end do
		ldf = ldf*4.d0/3.d0
		nuDensities%re(:,:) = nuDensities%re(:,:) * ldf * (cos2thW_Z)
		nuDensities%im(:,:) = nuDensities%im(:,:) * ldf * (cos2thW_Z)
	end subroutine updateMatterDensities

	pure function H_eff(y)
		real(dl), intent(in) :: y
		type(cmplxMatNN) :: H_eff

		call allocateCmplxMat(H_eff)

		!missing: term for NC!
		H_eff%re = 0.d0 &
			+ nuMassesMat(:,:)/(2.d0*y) &
			+ leptonDensities(:,:) * y &
			+ nuDensities%re(:,:) * y
		H_eff%im = 0.d0 &
			+ nuDensities%im(:,:) * y
	end function H_eff

	pure function H_eff_cmplx(y)
		complex(dl), dimension(maxFlavorNumber, maxFlavorNumber) :: H_eff_cmplx
		real(dl), intent(in) :: y
		type(cmplxMatNN) :: H
		integer :: i, j

		H = H_eff(y)
		H_eff_cmplx(:,:) = cmplx(0.d0, 0.d0)
		do i=1, flavorNumber
			do j=1, flavorNumber
				H_eff_cmplx(i, j) = cmplx(H%re(i,j), H%im(i,j))
			end do
		end do
	end function H_eff_cmplx

	pure function rho_diag_mass(iy)
		type(cmplxMatNN) :: rho_diag_mass
		integer, intent(in) :: iy
		complex(dl), dimension(maxFlavorNumber, maxFlavorNumber) :: tmpComplMat, transfMat
		real(dl), dimension(maxFlavorNumber) :: tmpvec
		integer :: i, k

		call allocateCmplxMat(rho_diag_mass)
		rho_diag_mass%re(:,:) = 0.d0
		rho_diag_mass%im(:,:) = 0.d0
		tmpvec = 0.d0

		transfMat(:,:) = cmplx(0.d0, 0.d0)
		tmpComplMat = H_eff_cmplx(y_arr(iy))
		call HEigensystem(flavorNumber, tmpComplMat, flavorNumber, tmpvec, transfMat, flavorNumber, 0)
		do k=1, flavorNumber
			do i=1, flavorNumber
				rho_diag_mass%re(k, k) = rho_diag_mass%re(k, k) &
					+ dble(conjg(transfMat(i, k))*transfMat(i, k)) &
						* nuDensMatVecFD(iy)%re(i, i)
			end do
		end do
	end function rho_diag_mass

end module fpMatter
