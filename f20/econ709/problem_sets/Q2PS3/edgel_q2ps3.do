/*
This file is used to conduct all empirical exercises from Problem Set 3 of
Econ709b.

Date created:  23 Nov 2020
Last modified: 23 Nov 2020
Author: Danny Edgel
*/
set more off

// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\f20\econ709\problem_sets\Q2PS3"

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

// (a) Estimate equation (3.50) and compute the equation R2 and sum of 
// squared errors.

// run regression to obtain R^2 and SSE
reg wage education experience exp2 
outreg2 using table1b.tex, replace ct("log(wage)") tex(fragment) lab 	///
	adds(Sum-of-squared Errors, e(rss))

/*
	(b) Re-estimate the slope on education using the residual regression 
		approach. Regress log(wage) on experience and its square, regress 
		education on experience and its square, and the residuals on the
		residuals. Report the estimates from this final regression, along 
		with the equation
*/


// log(wage) on experience
reg wage experience exp2 
predict e_wage, residual 	// save residuals from wage-experience regression 
outreg2 using table1b.tex, append ct("log(wage)") tex(fragment) lab 	///
	adds(Sum-of-squared Errors, e(rss))

// education on experience
reg education experience exp2
predict e_educ, residual 	// save residuals from education-experience regression 
outreg2 using table1b.tex, append ct("education") tex(fragment) lab	///
	adds(Sum-of-squared Errors, e(rss))

lab var e_educ "$\hat{\varepsilon}_{educ}$"

// wage residuals on education residuals 
reg e_wage e_educ
outreg2 using table1b.tex, append ct("$\hat{\varepsilon}_{wage}$") ///
	tex(fragment) lab adds(Sum-of-squared Errors, e(rss))



/*
	Exercise 3.25 Estimate equation (3.50) as in part (a) of the previous 
	question.
*/

// re-run regression for (3.50) and obtain yhat and ehat
reg wage education experience exp2 
predict yhat, xb
predict ehat, residual 

// for each of the following, simply save the answers to each sub-question 
// in a local macro, then report them all at once at the end of the code

// (a)
qui sum ehat		// generates summary table
loc a = round(r(sum),.001)

// (b)
g x1e = education * ehat
qui sum x1e
loc b = round(r(sum),.001)

// (c)
g x2e = experience * ehat
qui sum x2e
loc c = round(r(sum),.001)

// (d)
g x1sq_e = (education^2) * ehat
qui sum x1sq_e
loc d = round(r(sum),.001)

// (e)
g x2sq_e = (experience^2) * ehat
qui sum x2sq_e
loc e = round(r(sum),.001)

// (f)
g yhat_ehat = yhat * ehat
qui sum yhat_ehat 
loc f = round(r(sum),.001)

// (g)
g ehat_sq = ehat^2
qui sum ehat_sq 
loc g = round(r(sum),.001)



// display results:
di 	"(a) `a'" _newline	///
	"(b) `b'" _newline	///
	"(c) `c'" _newline	///
	"(d) `d'" _newline	///
	"(e) `e'" _newline	///
	"(f) `f'" _newline	///
	"(g) `g'" 
	
// output Tex file
file open resultsfile using "3.25_results.tex", write replace
file write resultsfile															///
	"\begin{itemize}" 												_newline	///
		_tab "\item[(a)] $\sum_{i=1}^n\hat{e}_i = `a'$"				_newline	///
		_tab "\item[(b)] $\sum_{i=1}^nX_{1i}\hat{e}_i = `b'$"		_newline	///
		_tab "\item[(c)] $\sum_{i=1}^nX_{2i}\hat{e}_i = `c'$"		_newline	///
		_tab "\item[(d)] $\sum_{i=1}^nX_{1i}^2\hat{e}_i = `d'$"		_newline	///
		_tab "\item[(e)] $\sum_{i=1}^nX_{2i}^2\hat{e}_i = `e'$"		_newline	///
		_tab "\item[(f)] $\sum_{i=1}^n\hat{Y}_i\hat{e}_i = `f'$"	_newline	///
		_tab "\item[(g)] $\sum_{i=1}^n\hat{e}^2_i = `g'$"			_newline	///
	"\end{itemize}"
file close resultsfile

// restore full data set
restore
