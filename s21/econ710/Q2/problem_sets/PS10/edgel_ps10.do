/*
This file is used to conduct all empirical exercises from Problem Set 10 of
Econ710q2.

Date created:  08 Apr 2021
Last modified: 08 Apr 2021
Author: Danny Edgel
*/
set more off
capture file close table23_8


// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ710\Q2\problem_sets\PS10"


/*
	Exercise 23.8
*/

// open the Papageorgiou, Saam, and Schulte (2017) data set
use PSS2017, clear



// save model estimates in local macros
loc Y 	EG_total
loc X1	EC_c_alt
loc X2	EC_d_alt

// generate log of Y for CES estimation
g logY = log(`Y')

// narrow to a subset of fully populated observations
keep if !missing(logY) & !missing(`X1') & !missing(`X2')


// estimate CES production function with NLLS
nl (logY = {beta} + ({v}/{rho})*log({a}*`X1'^{rho}+(1-{a})*`X2'^{rho})),	///
	initial(beta 1.2 v 1 rho .5 a .5)
	
// save coefficient vector and variance-covariance matrix
mat b = e(b)
mat V = e(V)

// calculate and save sigma's moments
nlcom 1/(1-[rho]_cons)
mat sig_b = r(b)
mat sig_V = r(V)

loc b_sigma :  di %4.3f sig_b[1,1]
loc se_sigma : di %4.3f sqrt(sig_V[1,1])


	
// save formatted estimates and SEs for each parameter
forval i =1/4{
	
	// extract parameter name from matrix column equations
	loc param `: word `i' of `: coleq b''
	
	// save parameter estimate and SE
	loc b_`param'  : di %4.3f b[1,`i']
	loc se_`param' : di %4.3f sqrt(V[`i',`i'])
}



// output results in table
file open table23_8 using table23_8.tex, write replace
file write table23_8															///
	"\begin{tabular}{ccc}" _newline "\hline\hline"					_newline	///
		_tab " Parameter & Estimate & Standard Error \\\hline"		_newline	///
		_tab	"$\rho$ & `b_rho' & `se_rho' \\"					_newline	///
				"$\nu$ & `b_v' & `se_v' \\"							_newline	///
				"$\alpha$ & `b_a' & `se_a' \\"						_newline	///
				"$\beta$ & `b_beta' & `se_beta' \\"					_newline	///
				"$\sigma$ & `b_sigma' & `se_sigma' \\\hline"		_newline	///
	"\end{tabular}"					
		
		
// close and save tex file
file close table23_8


/*
	Exercise 23.8
*/

// open the March 2009 CPS data set
use cps09mar, clear


// subset to Hispanic women
keep if female == 1 & hisp == 1


// generate secondary variables
g wage 	= earnings / (hours*week)
g lwage	= log(wage)


// perform a quantile regression of log(wages) on education for a few different
// quantiles
forval q = .1(.2).9{
	loc i = `q'*10
	
	qui qreg lwage education, q(`q')
	loc b0_`i' = _b[_cons]
	loc b1_`i' = _b[education]
}

// plot each quantile's curve
mylabels 0(10)50, prefix("$") local(ylabs)
tw	///
	function y = exp(`b0_1'+`b1_1'*x), ra(education)	||	///
	function y = exp(`b0_3'+`b1_3'*x), ra(education)	||	///
	function y = exp(`b0_5'+`b1_5'*x), ra(education)	||	///
	function y = exp(`b0_7'+`b1_7'*x), ra(education)	||	///
	function y = exp(`b0_9'+`b1_9'*x), ra(education)		///
		ylab(`ylabs', angle(horizontal) grid nogextend) ysc(r(0 52))			///
		plotregion(margin(tiny) lcolor(none)) ytitle("")						///
		title("Quantile Regression Results") xtitle("Education (years)") 		///
		subtitle("log(wage) on education for Hispanic women")					///
		leg( region(lstyle(none)) r(1)											///
			 lab(1 "10%") lab(2 "30%") lab(3 "50%") lab(4 "70%") lab(5 "90%"))
graph export fig24_14.png, replace
