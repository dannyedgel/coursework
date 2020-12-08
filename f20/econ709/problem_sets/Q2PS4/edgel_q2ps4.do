/*
This file is used to conduct all empirical exercises from Problem Set 4 of
Econ709b.

Date created:  06 Dec 2020
Last modified: 07 Dec 2020
Author: Danny Edgel
*/
set more off

// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\f20\econ709\problem_sets\Q2PS4"

// open CPS data
use cps09mar, clear


/*
	3.24: Use the data set from Section 3.22 and the sub-sample used for 
		equation (3.50)
*/

//	Housekeeping: generate subsample from Hansen, 3.25
//	(note: code following * comments directly copied from Hansen)


* Clear memory and load the data
clear
use cps09mar.dta

* Generate transformations
gen wage 		= ln(earnings/(hours*week))
gen experience 	= age - education - 6
gen exp2 		= (experience^2)/100

* Create indicator for subsamples
gen mbf 	= (race == 2) 	& (marital <= 2) 		& (female == 1)
gen mbf12 	= (mbf == 1) 	& (experience == 12)
gen sam 	= (race == 4) 	& (marital == 7) 		& (female == 0)

// label variables
lab var education 	"Education"
lab var experience 	"Experience"
lab var exp2		"$\text{Experience}^2$"

// subset data to subsample
preserve
keep if sam == 1 & experience < 45

//*** (a) Estimate the regression with robust SE's; output results

reg wage education experience exp2, robust
outreg2 using table1a.tex, replace ct("log(wage)") tex(fragment) lab 	///
	adds(Sum-of-squared Errors, e(rss))

	
//*** (b) compute and output thetahat

loc thetahat = _b[education] / (_b[experience] + (1/5)*_b[exp2])

// output Tex file
file open resultsfile using "1b.tex", write replace
file write resultsfile														///
	"\[" 														_newline	///
		_tab "\hat{\theta} \approx `=round(`thetahat',.001)'"	_newline	///
	"\]"
file close resultsfile


//*** (d) calculate asymptotic SE of thetahat and output

// save variance-covariance matrix from regression
mat V = e(V)

// define theta's gradient vector
loc gb1 = 1/(_b[experience] + (1/5)*_b[exp2])
loc gb2 = -_b[education] / (_b[experience] + (1/5)*_b[exp2])^2
loc gb3 = -_b[education] / (5*(_b[experience] + (1/5)*_b[exp2])^2)

mat gb = ( `gb1' \ `gb2' \ `gb3' \ 0 )

// calculate SE
mat thetavar = gb' * V * gb
loc thetase = round(sqrt(el(thetavar,1,1)),.001)
di `thetahat'
di `thetase'

// calculate 90% c.i.
loc ci_low 	= round(`thetahat' - 1.645*`thetase',.001)
loc ci_high	= round(`thetahat' + 1.645*`thetase',.001)
loc ci_high	 `=substr("`ci_high'",1,strpos("`ci_high'",".")+3)'

// output Tex file
file open resultsfile using "1d.tex", write replace
file write resultsfile														///
	"\begin{align*}" 											_newline	///
		_tab "s(\hat{\theta})  &\approx `thetase' \\"			_newline	///
		_tab "\text{90\% c.i.} &= [\hat{\theta}-1.645s(\hat{\theta})," 		///
				"\hat{\theta}+1.645s(\hat{\theta})]" 						///
		_tab "\approx [`ci_low',`ci_high']"						_newline	///
	"\end{align*}" 												
file close resultsfile

// restore full data set
restore
