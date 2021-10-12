/*
	Generate simulated environment
*/
clear
set obs 1000				// number of cities

// set parameters
gen F 	= 1
gen b0 	= 1
gen b1 	= 0
gen a0 	= 5
gen a1 	= 1

// initialize endogenous variables 
foreach var in N Lerner Herf lnHerf lnLernerObs Q{
    g `var' = .
}


set seed 155133
g unif = uniform() - 0.5

// loop through scenarios for (a) and (b)
forval i = 1/2{
    
	if (`i' == 1){
		set seed 24531
		gen nu	= runiform(-1, 1)
		gen eta = 0
	}
	else{
	    set seed 781385
	    replace nu 	= 0
		replace eta = runiform(-1, 1)
	}

	/*
		Calculate relevant endogenous variables
	*/
	
	replace N 			= ((a0 - b0 + nu - eta)/sqrt(F*a1)) - 1
	replace Q 			= (1/a1)*(a0 - b0 + nu - eta)*(N/(N+1))
	replace Lerner		= (a0 - b0 + nu - eta)/(a0 + nu + (b0 + eta)*N) 
	replace Herf		= 1 / N 
	replace lnHerf 		= ln(Herf)
	replace lnLernerObs = ln(Lerner) + .1*unif

	/*
		Conduct SCP regression and save results
	*/

	reg lnLernerObs lnHerf
	
	// save regression output
	loc N`i' = e(N)
	loc b`i' : di %4.3f _b[lnHerf]
	loc seb`i' : di %5.4f _se[lnHerf]
	loc a`i' : di %4.3f _b[_cons]
	loc sea`i' : di %5.4f _se[_cons]
	loc R`i' : di %4.3f e(r2)

}

/*
	Print output to table2.tex
*/ 

file write table2	///
	_tab "$\alpha$ 	& `a1' 		& `a2' 		\\"	_newline 	///
	_tab "			& (`sea1')	& (`sea2')	\\"	_newline 	///
	_tab "&&								\\"	_newline 	///
	_tab "$\beta$ 	& `b1' 		& `b2' 		\\"	_newline 	///
	_tab "			& (`seb1')	& (`seb2')	\\"	_newline 	///
	_tab "&& 								\\"	_newline 	///
	_tab "$ R^2$	& `R1'		& `R2'		\\" _newline	///
	_tab "N 		& `N1'		& `N2'		\\" _newline	///
	_tab "&& \\\hline"							_newline 	///
	"\end{tabular}"

file close table2 