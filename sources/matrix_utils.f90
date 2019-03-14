module ndMatrices
	use precision
	use ndErrors
	implicit none

	contains

	pure function numRows(mat)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat
		integer :: numRows
		integer, dimension(2) :: s

		s = shape(mat)
		numRows = s(1)
	end function numRows

	pure function numColumns(mat)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat
		integer :: numColumns
		integer, dimension(2) :: s

		s = shape(mat)
		numColumns = s(2)
	end function numColumns

	pure subroutine createEmptyMat(mat, nc)
		integer, intent(in) :: nc
		real(dl), dimension(:,:), allocatable, intent(out) :: mat

		if (.not. allocated(mat)) &
			allocate(mat(nc, nc))

		mat = 0.d0
	end subroutine createEmptyMat

	pure subroutine createIdentityMat(mat, nc)
		integer, intent(in) :: nc
		real(dl), dimension(:,:), allocatable, intent(out) :: mat
		integer :: ix

		call createEmptyMat(mat, nc)
		do ix = 1, nc
			mat(ix, ix) = 1.d0
		end do
	end subroutine createIdentityMat

	pure subroutine createDiagMat(mat, nc, d)
		integer, intent(in) :: nc
		real(dl), dimension(:), intent(in) :: d
		real(dl), dimension(:,:), allocatable, intent(out) :: mat
		integer :: ix

		call createEmptyMat(mat, nc)
		do ix=1, nc
			mat(ix,ix) = d(ix)
		end do
	end subroutine createDiagMat

	pure subroutine createRotMat(mat, nc, i1, i2, theta)
		integer, intent(in) :: nc, i1, i2
		real(dl), intent(in) :: theta
		real(dl), dimension(:,:), allocatable, intent(out) :: mat
		integer :: ix1, ix2

		if (i1.lt.i2) then
			ix1=i1
			ix2=i2
		else
			ix2=i1
			ix1=i2
		end if

		call createIdentityMat(mat, nc)

		mat(ix1, ix1) = cos(theta)
		mat(ix2, ix2) = cos(theta)
		mat(ix1, ix2) = sin(theta)
		mat(ix2, ix1) = -sin(theta)
	end subroutine createRotMat

	pure function traceMat(inM)
		real(dl), dimension(:,:), allocatable, intent(in) :: inM
		real(dl) :: traceMat
		integer :: i

		traceMat = 0.d0
		do i=1,minval(shape(inM))
			traceMat = traceMat + inM(i,i)
		end do
	end function traceMat

	pure subroutine transposeMat(inM, outM)
		real(dl), dimension(:,:), allocatable, intent(in) :: inM
		real(dl), dimension(:,:), allocatable, intent(out) :: outM
		integer, dimension(2) :: sizeM

		sizeM=shape(inM)

		if (.not. allocated(outM)) &
			allocate(outM(sizeM(2), sizeM(1)))

		outM = transpose(inM)
	end subroutine transposeMat

	pure subroutine minoreMat(inM, ix1, ix2, outM)
		integer, intent(in) :: ix1, ix2
		real(dl), dimension(:,:), allocatable, intent(in) :: inM
		real(dl), dimension(:,:), allocatable, intent(out) :: outM
		integer :: i,j,r,c
		integer, dimension(2) :: sizeM

		sizeM=shape(inM)

		if (.not. allocated(outM)) allocate(outM(sizeM(1)-1, sizeM(2)-1))

		outM = 0.d0
		r=0
		do i=1, sizeM(1)
			c=0
			if (i.ne.ix1) then
				r=r+1
				do j=1, sizeM(2)
					if (j .ne. ix2) then
						c=c+1
						outM(r,c) = inM(i,j)
					end if
				end do
			end if
		end do
	end subroutine minoreMat

	recursive function determMat(mat) result(determ)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat
		real(dl), dimension(:,:), allocatable :: mi
		real(dl) :: determ
		integer :: ix
		integer, dimension(2) :: s

		s=shape(mat)
		determ = 0.d0

		if (s(1) .ne. s(2)) &
			call criticalError("non-squared matrices have no det!")

		if (s(1) .eq. 1) then
			determ = mat(1,1)
		else if (s(1) .eq. 2) then
			determ = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
		else
			do ix=1, s(2)
				call minoreMat(mat, 1, ix, mi)
				determ = determ + determMat(mi) * mat(1,ix) * (1 - 2 * modulo(1+ix, 2) )
			end do
		end if
	end function determMat

	function cofactorMat(mat, i, j)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat
		real(dl), dimension(:,:), allocatable :: mi
		integer, intent(in) :: i,j
		real(dl) :: cofactorMat

		call minoreMat(mat, i, j, mi)
		cofactorMat = determMat(mi) * (1 - 2 * modulo(i+j, 2) )! for (-1)**(i+j)
	end function cofactorMat

	subroutine inverseMat(inM, outM)
		real(dl), dimension(:,:), allocatable, intent(in) :: inM
		real(dl), dimension(:,:), allocatable, intent(out) :: outM
		integer :: i,j,n
		integer, dimension(2) :: s
		real(dl) :: det

		s=shape(inM)
		det = determMat(inM)
		if (det.eq.0) then
			call criticalError("the matrix is degenerate!")
		end if
		if (s(1) .ne. s(2)) then
			call criticalError("non-squared matrices have no inverse!")
		else
			n=s(1)
		end if

		if (.not. allocated(outM)) allocate(outM(n,n))
		outM=0
		do i=1, n
			do j=1, n
				outM(i,j) = cofactorMat(inM, j, i)/det
			end do
		end do
	end subroutine inverseMat

	pure subroutine tripleProdMat(mat1, mat2, mat3, outMat)
		real(dl), dimension(:,:), intent(in) :: mat1, mat2, mat3
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat
		integer, dimension(2) :: size1, size3

		size1=shape(mat1)
		size3=shape(mat3)

		if (.not. allocated(outMat)) &
			allocate(outMat(size1(1), size3(2)))

		outMat = matmul(matmul(mat1, mat2), mat3)
	end subroutine tripleProdMat

	pure subroutine quadrupleProdMat(mat1, mat2, mat3, mat4, outMat)
		real(dl), dimension(:,:), intent(in) :: mat1, mat2, mat3, mat4
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat
		integer, dimension(2) :: size1, size4

		size1=shape(mat1)
		size4=shape(mat4)

		if (.not. allocated(outMat)) &
			allocate(outMat(size1(1), size4(2)))

		outMat = matmul(matmul(matmul(mat1, mat2), mat3), mat4)
	end subroutine quadrupleProdMat

	pure subroutine Commutator(mat1, mat2, outmat)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat1, mat2
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat

		call CommAntiComm(mat1, mat2, outmat, -1)
	end subroutine Commutator

	pure subroutine AntiCommutator(mat1, mat2, outmat)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat1, mat2
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat

		call CommAntiComm(mat1, mat2, outmat, 1)
	end subroutine antiCommutator

	pure subroutine CommAntiComm(mat1, mat2, outmat, s)
		!use the sign s to switch between commutator and anticommutator
		real(dl), dimension(:,:), allocatable, intent(in) :: mat1, mat2
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat
		integer, intent(in) :: s
		integer, dimension(2) :: size1
		integer :: d

		size1=shape(mat1)
		d = size1(1)

		allocate(outMat(d,d))
		outMat = matmul(mat1, mat2) + s * matmul(mat2, mat1)
	end subroutine

	subroutine printMat(m)
		real(dl), dimension(:,:), allocatable, intent(in) :: m
		integer, dimension(2) :: s
		integer :: i
		character(len=20) :: frmt

		s = shape(m)
		write (frmt,"('(',I4,'(E14.6))')") s(2)
		do i=1, s(1)
			write (*,frmt) m(i,:)
		end do
	end subroutine printMat

	subroutine printVec(vec,n)
		real(dl), dimension(:), intent(in) :: vec
		integer :: n
		character(len=20) :: frmt

		write (frmt,"('(',I4,'(E14.6))')") n
		write(*,trim(frmt)) vec(:)
	end subroutine printVec
end module ndMatrices
