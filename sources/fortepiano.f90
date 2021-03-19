program fortepiano
	use fpversion
	use precision
	use fpConfig
	use fpErrors
	use fpEquations
	implicit none

	write(*,*) "This is FortEPiaNO version "//version

	call initConfig

	call solver
	call finalresults

	call addToLog("Finished. closing log file.")
	call closeLogFile
	call renameLogFile(trim(outputFolder)//'/messages.log')
end program fortepiano
