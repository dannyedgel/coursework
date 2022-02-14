/*
This file is used to conduct all empirical exercises from Problem Set 12 of
Econ710q2.

Date created:  21 Apr 2021
Last modified: 22 Apr 2021
Author: Danny Edgel
*/
set more off
capture file close ex27_9b
capture file close ex28_12

// to run 27.9, set the below = 1
loc run_27_9 = 0


// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ710\Q2\problem_sets\PS12"


if (`run_27_9' == 1){
	/*
		Exercise 27.9
	*/

	// open the CHJ (2004) CPS data set
	use CHJ2004, clear

	// transform variables/generate secondary variables
	g		tinkind = transfers / 1000
	replace income 	= income / 1000
	g		dincome = (income-1)*(income > 1)

	// label variables of interest
	lab var tinkind "In-kind Transfers"
	lab var income "Income ('000)"
	lab var dincome "(Income-1)$\one{\text{Income}>1}$"

	// (a) run and report OLS
	reg tinkind income dincome, robust
	outreg2 using table27_9.tex, tex(fragment) lab replace ctitle("(a)")


	// (b) calculate and report the share of censored observations
	file open ex27_9b using 27_9b.tex, write replace
	qui count if tinkind == 0
	loc censored = r(N)
	qui count
	loc fwrite: di %3.2f (`censored'/`=r(N)')*100
	file write ex27_9b "`fwrite'\%"
	file close ex27_9b

	// (c) estimate with censored observations omitted and report
	preserve
	keep if tinkind > 0
	reg tinkind income dincome, robust
	outreg2 using table27_9.tex, tex(fragment) lab append 	///
		ctitle("(c)")

	restore

	// (d) estimate with tobit and report
	tobit tinkind income dincome
	outreg2 using table27_9.tex, tex(fragment) lab append ctitle("(d)")


	// (e) estimate with CLAD
	clad tinkind income dincome
	outreg2 using table27_9.tex, tex(fragment) lab append ctitle("(e)")

}

/*
	Exercise 28.12
*/

// open the March 2020 CPS data set
use cps09mar, clear

// write simple function for calculating and saving formatted return
// coefficients and standard errors to simplify rest of the notation
capture program drop expreturns
program expreturns, rclass
	
	syntax anything, LOCal(str)
	
	qui reg `anything', robust
	
	if ( "`: word 6 of `anything''" == "exp5" ){
		qui nlcom _b[experience]*30 + _b[exp2]*(30^2) + _b[exp3]*(30^3) + 	///
			_b[exp4]*(30^4) + _b[exp5]*(30^5) + _b[exp6]*(30^6)
	}
	else if ( "`: word 4 of `anything''" == "exp3" ){
		qui nlcom _b[experience]*30 + _b[exp2]*(30^2) + _b[exp3]*(30^3) + 	///
			_b[exp4]*(30^4)
	}
	else if ( "`: word 3 of `anything''" == "exp2" ){
		qui nlcom _b[experience]*30 + _b[exp2]*(30^2)
	}
	else di as err "model misspecified"
	mat b = r(b)
	mat V = r(V)	
	
	qui reg `anything', robust
	qui estat ic

	mat S = r(S)
	
	loc b 	= b[1,1]
	loc se	= sqrt(V[1,1])
	loc AIC = S[1,5]
	loc BIC = S[1,6]
	
	loc `local'_b 	: di %2.0f `b'*100
	loc `local'_se	: di %2.0f `se'*100
	
	return local `local'_b  	= ``local'_b'
	return local `local'_se 	= ``local'_se'
	return local `local'_AIC	= round(`AIC')
	return local `local'_BIC	= round(`BIC')
	
	end


// subset to Hispanic women
keep if female == 1 & hisp == 1 // & race == 4

// generate secondary variables
g lwage = log(earnings/(hours*week))
g experience = age - education - 6
forval i = 2/6{
	g exp`i' = experience^(`i')
}
g married 	= (inlist(marital,1,2,3) == 1)
g college 	= (education >= 16)
g colsp9	= (education - 9)*(education > 9)

foreach i in 12 13 14 16 18 20{
	g educ`i' = (education == `i')
	loc edu_dummies `edu_dummies' educ`i'
}

// save common covariates in local macro
loc X married i.region 
loc i = 0

// run regressions and generate returns
expreturns lwage experience exp2 college `X', loc(reg`++i')

loc reg`i'_b 	= r(reg`i'_b)
loc reg`i'_se 	= r(reg`i'_se)
loc reg`i'_AIC 	= r(reg`i'_AIC)
loc reg`i'_BIC 	= r(reg`i'_BIC)


expreturns lwage experience exp2 education colsp9 `X', loc(reg`++i')

loc reg`i'_b 	= r(reg`i'_b)
loc reg`i'_se 	= r(reg`i'_se)
loc reg`i'_AIC 	= r(reg`i'_AIC)
loc reg`i'_BIC 	= r(reg`i'_BIC)



expreturns lwage experience exp2 `edu_dummies' `X', loc(reg`++i')

loc reg`i'_b 	= r(reg`i'_b)
loc reg`i'_se 	= r(reg`i'_se)
loc reg`i'_AIC 	= r(reg`i'_AIC)
loc reg`i'_BIC 	= r(reg`i'_BIC)


expreturns lwage experience exp2 exp3 exp4 college `X', loc(reg`++i')

loc reg`i'_b 	= r(reg`i'_b)
loc reg`i'_se 	= r(reg`i'_se)
loc reg`i'_AIC 	= r(reg`i'_AIC)
loc reg`i'_BIC 	= r(reg`i'_BIC)


expreturns lwage experience exp2 exp3 exp4 education colsp9 `X', loc(reg`++i')

loc reg`i'_b 	= r(reg`i'_b)
loc reg`i'_se 	= r(reg`i'_se)
loc reg`i'_AIC 	= r(reg`i'_AIC)
loc reg`i'_BIC 	= r(reg`i'_BIC)



expreturns lwage experience exp2 exp3 exp4 `edu_dummies' `X', loc(reg`++i')

loc reg`i'_b 	= r(reg`i'_b)
loc reg`i'_se 	= r(reg`i'_se)
loc reg`i'_AIC 	= r(reg`i'_AIC)
loc reg`i'_BIC 	= r(reg`i'_BIC)

expreturns lwage experience exp2 exp3 exp4 exp5 exp6 college `X', loc(reg`++i')

loc reg`i'_b 	= r(reg`i'_b)
loc reg`i'_se 	= r(reg`i'_se)
loc reg`i'_AIC 	= r(reg`i'_AIC)
loc reg`i'_BIC 	= r(reg`i'_BIC)


expreturns lwage experience exp2 exp3 exp4 exp5 exp6 education colsp9 `X', 	///
	loc(reg`++i')

loc reg`i'_b 	= r(reg`i'_b)
loc reg`i'_se 	= r(reg`i'_se)
loc reg`i'_AIC 	= r(reg`i'_AIC)
loc reg`i'_BIC 	= r(reg`i'_BIC)



expreturns lwage experience exp2 exp3 exp4 exp5 exp6 `edu_dummies' `X', 	///
	loc(reg`++i')

loc reg`i'_b 	= r(reg`i'_b)
loc reg`i'_se 	= r(reg`i'_se)
loc reg`i'_AIC 	= r(reg`i'_AIC)
loc reg`i'_BIC 	= r(reg`i'_BIC)


// output results in a latex table
loc i = 1
file open ex28_12 using table28.12.tex, write replace
file write ex28_12	///
	"\begin{tabular}{lccccccccc}" _newline	"\hline\hline" ///
	_tab "\multicolumn{10}{c}{}	\\"									_newline	///
	_tab "&Model 1 &Model 2 &Model 3 &Model 4 &Model 5 &Model 6"	_newline 	/// 
			"&Model 7 &Model 8 &Model 9 \\\cline{2-10}"				_newline 	///
	_tab "Return & `reg`i++'_b'\% & `reg`i++'_b'\% & `reg`i++'_b'\% &"			///
			"`reg`i++'_b'\% & `reg`i++'_b'\% & `reg`i++'_b'\% & `reg`i++'_b'\%"	///
			"& `reg`i++'_b'\% & `reg`i++'_b'\% \\"	_newline
			
loc i = 1
file write ex28_12	///
	_tab "s.e. & `reg`i++'_se' & `reg`i++'_se' & `reg`i++'_se' & `reg`i++'_se'"	///
			"& `reg`i++'_se' & `reg`i++'_se' & `reg`i++'_se' & `reg`i++'_se' &"	///
			"`reg`i++'_se' \\"	_newline 
			
loc i = 1
file write ex28_12	///
	_tab "BIC & `reg`i++'_BIC' & \textbf{`reg`i++'_BIC'} & `reg`i++'_BIC' &"	///
			"`reg`i++'_BIC' & `reg`i++'_BIC' & `reg`i++'_BIC' & `reg`i++'_BIC'"	///
			"& `reg`i++'_BIC' & `reg`i++'_BIC' \\" _newline
			
loc i = 1
file write ex28_12	///
	_tab "AIC & `reg`i++'_AIC' & `reg`i++'_AIC' & `reg`i++'_AIC' &"				///
			"`reg`i++'_AIC' & `reg`i++'_AIC' & `reg`i++'_AIC' &"				///
			"`reg`i++'_AIC' & \textbf{`reg`i++'_AIC'} & `reg`i++'_AIC' \\"		///
															_newline			///
	_tab "Education & College & Spline & Dummy & College & Spline & Dummy &"	///
			"College & Spline & Dummy \\"	_newline							///
	_tab "Experience & 2 & 2 & 2 & 4 & 4 & 4 & 6 & 6 & 6"						///
	"\end{tabular}"
	

file close ex28_12