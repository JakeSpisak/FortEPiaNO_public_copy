program fortepiano
	use precision
	use fpConfig
	use fpErrors
	use fpEquations
	implicit none

	call initConfig

	call solver
	call finalresults

	call addToLog("Finished. closing log file.")
	call closeLogFile
	call renameLogFile(trim(outputFolder)//'/messages.log')
end program fortepiano
