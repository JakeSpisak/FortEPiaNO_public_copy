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

	subroutine assert_double(testname, num1, num2, tol)
		character(len=*) :: testname
		real(dl) :: num1, num2, tol

		totalTests = totalTests + 1
		if (abs(num1 - num2) .lt. tol) then
			write(*, fmt="(a)", advance="no") "."
			successTests = successTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, num1 - num2, tol
			failedTests = failedTests + 1
			if (blocking) &
				call exit()
		end if
	end subroutine assert_double

	subroutine assert_double_verb(testname, num1, num2, tol)
		character(len=*) :: testname
		real(dl) :: num1, num2, tol

		totalTests = totalTests + 1
		if (abs(num1 - num2) .lt. tol) then
			write(*, fmt="(a,'  ',E14.7)") testname, num1 - num2
			successTests = successTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, num1 - num2, tol
			failedTests = failedTests + 1
			if (blocking) &
				call exit()
		end if
	end subroutine assert_double_verb

	subroutine assert_double_rel(testname, num1, num2, tol)
		character(len=*) :: testname
		real(dl) :: num1, num2, tol

		totalTests = totalTests + 1
		if (abs((num1 - num2)/num1) .lt. tol) then
			write(*, fmt="(a)", advance="no") "."
			successTests = successTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, (num1 - num2)/num1, tol
			failedTests = failedTests + 1
			if (blocking) &
				call exit()
		end if
	end subroutine assert_double_rel

	subroutine assert_double_rel_verb(testname, num1, num2, tol)
		character(len=*) :: testname
		real(dl) :: num1, num2, tol

		totalTests = totalTests + 1
		if (abs((num1 - num2)/num1) .lt. tol) then
			write(*, fmt="(a,'  ',E14.7)") testname, (num1 - num2)/num1
			successTests = successTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, (num1 - num2)/num1, tol
			failedTests = failedTests + 1
			if (blocking) &
				call exit()
		end if
	end subroutine assert_double_rel_verb
end module sgTestUtils
