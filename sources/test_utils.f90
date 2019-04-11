module sgTestUtils
	use precision
	implicit None

	integer :: totalTests = 0

	contains

	subroutine resetTestCounter
		totalTests = 0
	end subroutine resetTestCounter

	subroutine printTotalTests()
		write(*,*) ""
		write(*,"(I7,' tests have been performed successfully.')") totalTests
		write(*,*) ""
	end subroutine printTotalTests

	subroutine assert_double(testname, num1, num2, tol)
		character(len=*) :: testname
		real(dl) :: num1, num2, tol

		if (abs(num1 - num2) .lt. tol) then
			write(*, fmt="(a)", advance="no") "."
			totalTests = totalTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, num1 - num2, tol
			call exit()
		end if
	end subroutine assert_double

	subroutine assert_double_verb(testname, num1, num2, tol)
		character(len=*) :: testname
		real(dl) :: num1, num2, tol

		if (abs(num1 - num2) .lt. tol) then
			write(*, fmt="(a,'  ',E14.7)") testname, num1 - num2
			totalTests = totalTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, num1 - num2, tol
			call exit()
		end if
	end subroutine assert_double_verb

	subroutine assert_double_rel(testname, num1, num2, tol)
		character(len=*) :: testname
		real(dl) :: num1, num2, tol

		if (abs((num1 - num2)/num1) .lt. tol) then
			write(*, fmt="(a)", advance="no") "."
			totalTests = totalTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, (num1 - num2)/num1, tol
			call exit()
		end if
	end subroutine assert_double_rel

	subroutine assert_double_rel_verb(testname, num1, num2, tol)
		character(len=*) :: testname
		real(dl) :: num1, num2, tol

		if (abs((num1 - num2)/num1) .lt. tol) then
			write(*, fmt="(a,'  ',E14.7)") testname, (num1 - num2)/num1
			totalTests = totalTests + 1
		else
			print *, testname, " failed"
			write (*,"(*(E17.9))") num1, num2, (num1 - num2)/num1, tol
			call exit()
		end if
	end subroutine assert_double_rel_verb
end module sgTestUtils
