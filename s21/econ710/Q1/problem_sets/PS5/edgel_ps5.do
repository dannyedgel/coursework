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
loc N = 10000


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
		_tab "T & $\rho_1$ & $\alpha_0$ & $\delta_0$ & $\rho_1$"					///
				"\\\hline" 											_newline 	///
		_tab "&&&& \\"												_newline 	
		
		
// initialize output file for 2(ii)
file open q2ii using 2ii.tex, write replace
file write q2ii																	///
	"\begin{tabular}{cc|cccccc}" _newline "\hline\hline"			_newline	///
		_tab "\multicolumn{8}{c}{ } \\"								_newline 	///
		_tab "&&\multicolumn{2}{c}{$\alpha_0$}"									///
				"&\multicolumn{2}{c}{$\delta_0$}"								///
				"&\multicolumn{2}{c}{$\rho_1$} \\"					_newline 	///
		_tab "T &  $\rho_1$ & Mean & Coverage & Mean & Coverage &"				///
				"Mean & Coverage \\\hline" 							_newline 	///
		_tab "&&&&&&& \\"											_newline 	



// loop through possible simulations for each question

loc c = 0		// column counter
//set trace on
foreach T in `T_vals'{
	forval i = 1/`: word count `rho_vals''{
		
		loc rho = `: word `i' of `rho_vals''
		loc c = `c' + 1	
		
		// initialize current setup's matrix of coefficient estimates
		mat dat_`T'_`i' = J(`N', 9, .)
		mat coln dat_`T'_`i' =  a0 a0_lb a0_ub d0 d0_lb d0_ub r1 r1_lb r1_ub
		
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
				
				// save coefficients and 95% c.i.'s in current model's matrix
				mat dat_`T'_`i'[`sim', 1] = _b[_cons]
				mat dat_`T'_`i'[`sim', 2] = _b[_cons] - 1.96*_se[_cons]
				mat dat_`T'_`i'[`sim', 3] = _b[_cons] + 1.96*_se[_cons]
				mat dat_`T'_`i'[`sim', 4] = _b[X]
				mat dat_`T'_`i'[`sim', 5] = _b[X] - 1.96*_se[X]
				mat dat_`T'_`i'[`sim', 6] = _b[X] + 1.96*_se[X]
				mat dat_`T'_`i'[`sim', 7] = _b[l.Y]
				mat dat_`T'_`i'[`sim', 8] = _b[l.Y] - 1.96*_se[l.Y]
				mat dat_`T'_`i'[`sim', 9] = _b[l.Y] + 1.96*_se[l.Y]
				
				// only output first simulation's estimates
				if ( `sim' == 1 ){
					// save coefficients and their SEs
					loc ba0`c' = _b[_cons]
					loc sa0`c' = _se[_cons]
					loc bd0`c' = _b[t]
					loc sd0`c' = _se[t]
					loc br1`c' = _b[l.Y]
					loc sr1`c' = _se[l.Y]
				}
			
			} // end quietly statement
		
		} // end simulation loop
		
		di "" // print blank line when model changes
		
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
		loc b_lb : di %4.3f `=`bd0`c''-1.96*`sd0`c'''
		loc b_ub : di %4.3f `=`bd0`c''+ 1.96*`sd0`c'''
		loc r_lb : di %4.3f `=`br1`c''-1.96*`sr1`c'''
		loc r_ub : di %4.3f `=`br1`c''+ 1.96*`sr1`c'''
		
		loc a0 : di %4.3f `ba0`c''
		loc d0 : di %4.3f `bd0`c''
		loc r1 : di %4.3f `br1`c''
		
		// write current loop's results to table
		file write q2i															///
			_tab "`T' & `rho' & `a0' & `d0' & `r1' \\"				_newline 	///
			_tab "& & [`a_lb',`a_ub'] & [`b_lb',`b_ub'] & [`r_lb',`r_ub']"		///
				 "\\ &&&& \\"										_newline	
	}
}

// close LaTeX table and save (i) file
file write q2i _tab "\\\hline\hline" _newline "\end{tabular}"
file close q2i



/*
		(ii) Report simulated means and coverage ratios
*/

foreach T in `T_vals'{
	forval i = 1/`: word count `rho_vals''{
		
		loc rho = `: word `i' of `rho_vals''
		
		// load data from current loop's matrix
		clear
		qui svmat dat_`T'_`i', n(col)
		
		// generate means and coverage rates for each coefficient
		foreach x in a0 d0 r1{
			qui sum `x'
			loc `x'_m = r(mean)
			
			qui count if `x'_lb <= ``x'_m' & `x'_ub >= ``x'_m'
			loc `x'_c = r(N) / `N'
			
			loc `x'_m : di %4.3f ``x'_m'
			loc `x'_c : di %4.3f ``x'_c'
		}
		
		// write current loop's results to table
		file write q2ii															///
			_tab "`T' & `rho' & `a0_m' & `a0_c' & `d0_m' & `d0_c' &"			///
					"`r1_m' & `r1_c' \\"							_newline 	///
				 " &&&&&&& \\"										_newline	
		
	}
}

// close LaTeX table and save (ii) file
file write q2ii _tab "\\\hline\hline" _newline "\end{tabular}"
file close q2ii