/*
	This file conducts all of the analysis necessary for question 2 of problem
	set 2 for Econ 761
	
	Date created:  04 Oct 2021
	Last modified: 08 Oct 2021
	Author: Danny Edgel
*/

/*
	Housekeeping
*/

clear
capture file close table1
capture file close table2


// initialize output tables
file open table1 using table1.tex, write replace
file write table1	///
	"\begin{tabular}{r|ccccc}" _newline											///
	_tab "& & & \multicolumn{2}{c}{Test: $\beta_{\loge{H}}=1$} & \\" _newline	///
	_tab " & $\beta_{\loge{H}}$ & se$\left(\beta_{\loge{H}}\right)$" _newline	///
			"& \textit{F}-score  & \textit{p}-score & N \\\hline"	 _newline
			
file open table2 using table2.tex, write replace
file write table2	///
	"\begin{tabular}{r|cc}"	_newline											///
	_tab "& $\nu\sim U[-1,1]$ & $\eta\sim U[-1,1]$	\\\hline && \\" _newline


//______________________________________________________________________________

/*
	Question 2
*/
do ps2_q2.do

//______________________________________________________________________________

/*
	Question 3
*/

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
	
	replace N 			= sqrt((1/4) + (((a0 - b0 + nu - eta)^2)/(F*a1))) 
	replace Q 			= (1/a1)*(a0 - b0 + nu - eta)*N*(N + 1)
	replace Lerner		= -(a1*Q)/((a0 + nu - a1*Q)*N)
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

//______________________________________________________________________________