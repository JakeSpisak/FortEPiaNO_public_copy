module FileIni
	use precision

	integer, parameter :: fileunitin=1234
	integer, parameter :: fileunitout=1235
	
	character(LEN=*), parameter :: floatformat = '(*(E16.7))'

	contains
	
	subroutine ini_file_open(filename, logname)
		character(len=*), intent(in) :: filename, logname
		
		open(unit=fileunitin, file=trim(filename))
		open(unit=fileunitout, file=trim(logname))
	end subroutine ini_file_open
	
	subroutine ini_file_close
		close(fileunitin)
		close(fileunitout)
	end subroutine ini_file_close
	
	function read_ini_int(paramname, default_val)
		integer read_ini_int
		character(len=*), intent(in) :: paramname
		integer, intent(in) :: default_val
		character(len=300) :: tmpstr
		character(len=16)  :: tmpval
		integer temp
		
		write(tmpval,*) default_val
		tmpstr=read_ini_par_def(paramname,tmpval)
		read(tmpstr,*) temp
		read_ini_int=temp
		
		call writelog_int(paramname, read_ini_int)
	end function read_ini_int
	
	function read_ini_real(paramname, default_val)
		real(dl) read_ini_real
		character(len=*), intent(in) :: paramname
		real(dl),intent(in) :: default_val
		character(len=300) :: tmpstr
		character(len=16)  :: tmpval
		real(dl) temp
		
		write(tmpval, floatformat) default_val
		tmpstr=read_ini_par_def(paramname,tmpval)
		read(tmpstr,*) temp
		read_ini_real=temp
		
		call writelog_real(paramname, read_ini_real)
	end function read_ini_real

	function read_ini_logical(paramname, default_val)
		logical :: read_ini_logical
		character(len=*), intent(in) :: paramname
		logical, intent(in) :: default_val
		character(len=300) :: tmpstr
		character(len=16) :: tmpval
		logical temp

		write(tmpval, *) default_val
		tmpstr=read_ini_par_def(paramname,tmpval)
		read(tmpstr,*) temp
		read_ini_logical=temp
		
		call writelog_logical(paramname, read_ini_logical)
	end function read_ini_logical

	function read_ini_char(paramname)
		character(len=300) read_ini_char
		character(len=*) paramname
		
		read_ini_char=trim(read_ini_par(paramname))
		
		call writelog_char(paramname, read_ini_char)
	end function read_ini_char
	
	function read_ini_par_def(paramname, default_val)
		character(len=300) read_ini_par_def
		character(len=*) paramname, default_val
		
		read_ini_par_def=trim(read_ini_par(paramname))//default_val
	end function read_ini_par_def

	function read_ini_par(paramname)
		character(len=300) read_ini_par
		integer ii,jj
		character(len=300) currline
		character(len=*) paramname
		
		read_ini_par=""
		
		rewind(unit=fileunitin)
	10	read (fileunitin,'(A)',end=30) currline
		jj = index (currline, trim(paramname))
		if (jj .eq. 1) then
			goto 20
		endif
		goto 10
	20	ii=len(paramname)
		jj=index(currline,'#')
		if (jj.gt.0) currline=currline(:jj-1)
		currline=adjustl(currline(ii+1:))
		read_ini_par=adjustl(currline(index(currline,'=')+1:))
	30	continue
	end function read_ini_par
	
	subroutine writelog_char(paramname, parval)
		character(len=*), intent(in) :: paramname, parval
	
		write (fileunitout,*) trim(paramname)//' = '//trim(parval)
	end subroutine writelog_char
	
	subroutine writelog_real(paramname, parval)
		character(len=*), intent(in) :: paramname
		real(dl), intent(in) :: parval
		character(len=16) :: tmp
	
		write (tmp, floatformat) parval
		write (fileunitout,*) trim(paramname)//' = '//trim(tmp)
	end subroutine writelog_real
	
	subroutine writelog_int(paramname, parval)
		character(len=*), intent(in) :: paramname
		integer, intent(in) :: parval
		character(len=16) :: tmp
	
		write (tmp, *) parval
		write (fileunitout,*) trim(paramname)//' = '//trim(tmp)
	end subroutine writelog_int
	
	subroutine writelog_logical(paramname, parval)
		character(len=*), intent(in) :: paramname
		logical, intent(in) :: parval
		character(len=3) :: tmp
	
		write (tmp, *) parval
		write (fileunitout,*) trim(paramname)//' = '//trim(tmp)
	end subroutine writelog_logical

end module FileIni
