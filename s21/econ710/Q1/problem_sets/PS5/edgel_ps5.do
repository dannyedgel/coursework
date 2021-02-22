/*
This file is used to conduct all empirical exercises from Problem Set 5 of
Econ710q1.

Date created:  22 Feb 2021
Last modified: 22 Feb 2021
Author: Danny Edgel
*/
set more off
capture file close q2i
capture file close q2ii


// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ710\Q1\problem_sets\PS5"

// set number of simulations
loc N = 1	//10000


// save values for coefficients
loc a0 = 1
loc d0 = 1
loc b0 = 1/100


// set seed for reproducability
set seed 546153


// save possible values of T and rho in local macros for looping
loc T_vals		50 150 250
loc rho_vals	0.7 0.9 0.95

// initialize output file for 2(i)
file open q2i using 2i.tex, write replace
file write q2i																	///
	"\begin{tabular}{cc|ccc}" _newline "\hline\hline"				_newline	///
		_tab "\multicolumn{5}{c}{ } \\"								_newline 	///
		_tab "T & $\rho_0$ & $\alpha_0$ & $\beta_0$ & $\rho_1$"					///
				"\\\hline" 											_newline 	///
		_tab "&&&& \\"												_newline 	



// loop through possible simulations for each question

loc c = 0		// column counter
//set trace on
foreach T in `T_vals'{
	forval i = 1/`: word count `rho_vals''{
		
		loc rho = `: word `i' of `rho_vals''
		loc c = `c' + 1	
		
		// initialize current setup's matrix of coefficient estimates
		mat dat_`T'_`i' = J(`N', 3, .)
		mat coln dat_`T'_`i' =  "alpha0" "beta0" "rho1"
		
		forval sim = 1/`N'{
			
			// every 1,000 simulations, display progress
			loc prog "T = `T'; rho = `rho'; sim `sim' of `N'..."
			if ( mod(`sim',1000) == 0 ) di "`prog'"
			
			qui{
				// clear currently-loaded data, set observations according to T 
				clear
				set obs `T'
				
				// generate simulated errors and data for Y and X
				g 		V = rnormal()
				g 		X = rnormal() if _n == 1
				
				forval t = 2/`T'{
					replace X = 0.3*X[_n-1] + V if _n == `t'
				}
				
				
				g 		U = rnormal()
				g		Y = rnormal() if _n == 1
				g 		t = _n
				
				forval t = 2/`T'{
					replace	Y = `a0' + t*`b0' + `d0'*X + `rho'*Y[_n-1] + U 	///
						if _n == `t'
				}
				
				
				// establish data as a time series
				tsset t
				
				// estimate OLS coefficients 
				reg Y t X l.Y
				
				// save coefficients in current model's matrix
				mat dat_`T'_`i'[`sim', 1] = _b[_cons]
				mat dat_`T'_`i'[`sim', 1] = _b[t]
				mat dat_`T'_`i'[`sim', 1] = _b[l.Y]
				
				// only output first simulation's estimates
				if ( `sim' == 1 ){
					// save coefficients and their SEs
					loc ba0`c' = _b[_cons]
					loc sa0`c' = _se[_cons]
					loc bb0`c' = _b[t]
					loc sb0`c' = _se[t]
					loc br1`c' = _b[l.Y]
					loc sr1`c' = _se[l.Y]
				}
			
			} // end quietly statement
		
		} // end simulation loop
	} // end rho loop
}  // end T loop

/*
			(i) Report OLS coefficients and C.I.'s
*/

loc c = 0
foreach T in `T_vals'{
	foreach rho in `rho_vals'{
		
		loc c = `c' + 1	
		
		// calculate c.i.'s and format how each value is displayed
		loc a_lb : di %4.3f `=`ba0`c''-1.96*`sa0`c'''
		loc a_ub : di %4.3f `=`ba0`c''+ 1.96*`sa0`c'''
		loc b_lb : di %4.3f `=`bb0`c''-1.96*`sb0`c'''
		loc b_ub : di %4.3f `=`bb0`c''+ 1.96*`sb0`c'''
		loc r_lb : di %4.3f `=`br1`c''-1.96*`sr1`c'''
		loc r_ub : di %4.3f `=`br1`c''+ 1.96*`sr1`c'''
		
		loc a0 : di %4.3f `ba0`c''
		loc b0 : di %4.3f `bb0`c''
		loc r1 : di %4.3f `br1`c''
		
		file write q2i															///
			_tab "`T' & `rho' & `a0' & `b0' & `r1' \\"				_newline 	///
			_tab "& & [`a_lb',`a_ub'] & [`b_lb',`b_ub'] & [`r_lb',`r_ub']"		///
				 "\\ &&&& \\"										_newline	
	}
}

file write q2i _tab "\multicolumn{5}{c}{ } \\\hline\hline" _newline 			///
					"\end{tabular}"
file close q2i

/*
		(ii) 
*/