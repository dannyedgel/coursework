/*
This file is used to conduct all empirical exercises from Problem Set 12 of
Econ710q2.

Date created:  21 Apr 2021
Last modified: 21 Apr 2021
Author: Danny Edgel
*/
set more off
capture file close ex27_9b

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


// subset to Hispanic women
keep if female == 1 & hisp == 1

// generate secondary variables
g lwage = log(earnings/(hours*week))
g experience = age - education - 6
forval i = 2/6{
	g exp`i' = experience^(`i')
}
g married 	= (inlist(marital,1,2,3) == 1)
g college 	= (education >= 16)
g colsp9	= (education - 9)*(education >= 9)

foreach i in 12 13 14 15 18 20{
	g educ`i' = (education == `i')
	loc edu_dummies `edu_dummies' educ`i'
}

// save common covariates in local macro
loc X married i.region 


// run regressions and generate returns
reg lwage experience exp2 college `X', robust 

lincom _b[experience]*30 + _b[exp2]*(30^2)
