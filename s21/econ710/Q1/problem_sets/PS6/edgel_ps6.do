/*
This file is used to conduct all empirical exercises from Problem Set 6 of
Econ710q1.

Date created:  01 Mar 2021
Last modified: 01 Mar 2021
Author: Danny Edgel
*/
set more off
capture file close q3


// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ710\Q1\problem_sets\PS6"

// set number of simulations and within-panel time observations
loc N = 10000
loc T = 4


// save values for coefficients
loc b0 = 1
loc d1 = 0
loc d2 = 1
loc d3 = 1
loc d4 = 1


// set seed for reproducability
set seed 846152


// save possible values of n and rho in local macros for looping
loc n_vals		40 70 100
loc phi_vals	0 0.8


// initialize output file for q3
file open q3 using q3.tex, write replace
file write q3																	///
	"\begin{tabular}{cc|ccccccc}" _newline "\hline\hline \\"		_newline	///
		_tab "&&&\multicolumn{2}{c}{OLS} &\multicolumn{2}{c}{FE, Robust}"		///
				"&\multicolumn{2}{c}{FE, Cluster} \\"				_newline 	///
		_tab "n & $\phi$ && Mean & Coverage & Mean & Coverage & "				///
				"Mean & Coverage \\\hline" 							_newline 	///
		_tab "&&&&&&&& \\"											_newline 	
		


// loop through possible simulations for each question

loc c = 0		// column counter
//set trace on
foreach n in `n_vals'{
	forval i = 1/`: word count `phi_vals''{
		
		loc phi = `: word `i' of `phi_vals''
		loc c = `c' + 1	
		
		// initialize current setup's matrix of coefficient estimates
		foreach m in ols rb cl{
			mat `m'_dat_`n'_`i' = J(`N', 12, .)
			mat coln `m'_dat_`n'_`i' =	b0 b0_lb b0_ub d2 d2_lb d2_ub d3 d3_lb 	///
										d3_ub d4 d4_lb d4_ub
		}
		
		
		forval sim = 1/`N'{
			
			// every 1,000 simulations, display progress
			loc prog "n = `n'; phi = `phi'; sim `sim' of `N'..."
			if ( mod(`sim',1000) == 0 ) di "`prog'"
			
			qui{
			
				// clear currently-loaded data, set observations according to T 
				clear
				set obs `n'
				
				// generate time and panel variables
				g i = floor((_n-1)/4) + 1
				g t = mod(_n-1, 4) + 1
				
				
				// declare data as a panel
				xtset t i
				sort i t
				
				// generate simulated errors
				g	u 	= rnormal()
				g	eps	= `phi'*rnormal() + u if t == 1
				bys i (t): replace eps = `phi'*eps[_n-1] + u if t > 1
				
				// generate model variables
				g 		X = 0 if inlist(t, 1, 2) == 1
				replace X = 1 if inlist(t, 3, 4) == 1
				
				g ai = rnormal()
				bys i (t): replace ai = ai[1] if _n > 1
				g ai_ind = ( ai > 0.6 )
				
				replace X = X*ai_ind 
				drop ai_ind
				
				g dt = 0
				forval t = 1/`T'{
					replace dt = `d`t'' if t == `t'
				}
				
				g Y = X*`b0' + dt + ai + eps
				
				
				// estimate OLS coefficients, alternatively using robust and 
				// cluster-robust standard errors
				foreach m in ols rb cl{
					
					if ( "`m'" == "rb" ) 	loc spec "robust nocons"
					else					loc spec "cl(i) nocons"
					
					
					if ( "`m'" == "ols" )	reg Y X 2.t 3.t 4.t
					else 					reg Y X 2.t 3.t 4.t i.i, `spec' 
					
					// save coefficients and 95% c.i.'s in current model's matrix
					mat `m'_dat_`n'_`i'[`sim',  1] = _b[X]
					mat `m'_dat_`n'_`i'[`sim',  2] = _b[X] - 1.96*_se[X]
					mat `m'_dat_`n'_`i'[`sim',  3] = _b[X] + 1.96*_se[X]
					mat `m'_dat_`n'_`i'[`sim',  4] = _b[2.t]
					mat `m'_dat_`n'_`i'[`sim',  5] = _b[2.t] - 1.96*_se[2.t]
					mat `m'_dat_`n'_`i'[`sim',  6] = _b[2.t] + 1.96*_se[2.t]
					mat `m'_dat_`n'_`i'[`sim',  7] = _b[3.t]
					mat `m'_dat_`n'_`i'[`sim',  8] = _b[3.t] - 1.96*_se[3.t]
					mat `m'_dat_`n'_`i'[`sim',  9] = _b[3.t] + 1.96*_se[3.t]
					mat `m'_dat_`n'_`i'[`sim', 10] = _b[4.t]
					mat `m'_dat_`n'_`i'[`sim', 11] = _b[4.t] - 1.96*_se[4.t]
					mat `m'_dat_`n'_`i'[`sim', 12] = _b[4.t] + 1.96*_se[4.t]

				} // end m loop
				
			} // end quietly statement
			
		} // end simulation loop
		
		di "" // print blank line when model changes
		
	} // end rho loop
}  // end T loop



/*
		Report simulated means and coverage ratios
*/

foreach n in `n_vals'{
	forval i = 1/`: word count `phi_vals''{
						
		loc phi = `: word `i' of `phi_vals''
		
		foreach m in ols rb cl{
			
			// load data from current loop's matrix
			clear
			qui svmat `m'_dat_`n'_`i', n(col)
			
			// generate mean and coverage rate for beta0
			qui sum b0
			loc b0_m = r(mean)
			
			qui count if b0_lb <= `b0_m' & b0_ub >= `b0_m'
			loc b0_c = r(N) / `N'
			
			loc `m'_b0_m : di %4.3f `b0_m'
			loc `m'_b0_c : di %4.3f `b0_c'
		}
		
		// write current loop's results to table
		file write q3															///
			_tab "`n' & `phi' && `ols_b0_m' & `ols_b0_c' "						///
					"& `rb_b0_m' & `rb_b0_c' & `cl_b0_m' & `cl_b0_c' \\"		///
					_newline " &&&&&&&& \\"	_newline	
		
	}
}

// close LaTeX table and save (ii) file
file write q3 _tab "\\\hline\hline" _newline "\end{tabular}"
file close q3