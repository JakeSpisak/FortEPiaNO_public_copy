module sgTestUtils
	use precision
	implicit None

	integer :: totalTests = 0
	integer :: successTests = 0
	integer :: failedTests = 0
	logical :: blocking = .true.

	contains

	subroutine resetTestCounter
		totalTests = 0
		successTests = 0
		failedTests = 0
	end subroutine resetTestCounter

	subroutine printTestBlockName(testname)
		character(len=*), intent(in) :: testname
		character(len=39) :: lname
		write(lname,*) testname
		write(*,*) ""
		write(*,"('###########################################################')")
		write(*,"('# starting tests: ',A,' #')") adjustl(lname)
		write(*,"('###########################################################')")
	end subroutine printTestBlockName

	subroutine printTotalTests()
		write(*,*) ""
		write(*,"('###########################################################')")
		write(*,"('#',I4,' tests performed: ',I4,' successfully and ',I4,' failed. #')") &
			totalTests, successTests, failedTests
		write(*,"('###########################################################')")
		write(*,*) ""
	end subroutine printTotalTests

	subroutine assert_logical(testname, num1, num2)
		character(len=*), intent(in) :: testname
		logical, intent(in) :: num1, num2

		totalTests = totalTests + 1
		if (num1 .eqv. num2) then
			write(*, fmt="(a)", advance="no") "."
			successTests = successTests + 1
		else
			print *, testname, " failed"
			write (*,*) num1, num2
			failedTests = failedTests + 1
			if (blocking) &
				call exit(1)
		end if
	end subroutine assert_logical

	subroutine assert_double(testname, num1, num2, tol)
		character(len=*), intent(in) :: testname
		real(dl), intent(in) :: num1, num2, tol

		totalTests = totalTests + 1
		if (abs(num1 - num2) .lt. tol) then
			write(*, fmt="(a)", advance="no") "."
			successTests = successTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, num1 - num2, tol
			failedTests = failedTests + 1
			if (blocking) &
				call exit(1)
		end if
	end subroutine assert_double

	subroutine assert_double_verb(testname, num1, num2, tol)
		character(len=*), intent(in) :: testname
		real(dl), intent(in) :: num1, num2, tol

		totalTests = totalTests + 1
		if (abs(num1 - num2) .lt. tol) then
			write(*, fmt="(a,'  ',E14.7)") testname, num1 - num2
			successTests = successTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, num1 - num2, tol
			failedTests = failedTests + 1
			if (blocking) &
				call exit(1)
		end if
	end subroutine assert_double_verb

	subroutine assert_double_rel(testname, num1, num2, tol)
		character(len=*), intent(in) :: testname
		real(dl), intent(in) :: num1, num2, tol

		totalTests = totalTests + 1
		if (abs((num1 - num2)/num1) .lt. tol) then
			write(*, fmt="(a)", advance="no") "."
			successTests = successTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, (num1 - num2)/num1, tol
			failedTests = failedTests + 1
			if (blocking) &
				call exit(1)
		end if
	end subroutine assert_double_rel

	subroutine assert_double_rel_safe(testname, num1, num2, tola, tolr)
		character(len=*), intent(in) :: testname
		real(dl), intent(in) :: num1, num2, tola, tolr

		totalTests = totalTests + 1
		if(abs(num2).lt.tola) then
			if (abs(num1 - num2) .lt. tola) then
				write(*, fmt="(a)", advance="no") "."
				successTests = successTests + 1
			else
				print *, testname, " failed (safe-abs)"
				write (*,"(*(E17.9))") num1, num2, num1 - num2, tola
				failedTests = failedTests + 1
				if (blocking) &
					call exit(1)
			end if
		else
			if (abs((num1 - num2)/num1) .lt. tolr) then
				write(*, fmt="(a)", advance="no") "."
				successTests = successTests + 1
			else
				print *, testname, " failed (safe-rel)"
				write (*,"(*(E17.9))") num1, num2, (num1 - num2)/num1, tolr
				failedTests = failedTests + 1
				if (blocking) &
					call exit(1)
			end if
		end if
	end subroutine assert_double_rel_safe

	subroutine assert_double_rel_verb(testname, num1, num2, tol)
		character(len=*), intent(in) :: testname
		real(dl), intent(in) :: num1, num2, tol

		totalTests = totalTests + 1
		if (abs((num1 - num2)/num1) .lt. tol) then
			write(*, fmt="(a,'  ',E14.7)") testname, (num1 - num2)/num1
			successTests = successTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, (num1 - num2)/num1, tol
			failedTests = failedTests + 1
			if (blocking) &
				call exit(1)
		end if
	end subroutine assert_double_rel_verb
end module sgTestUtils
