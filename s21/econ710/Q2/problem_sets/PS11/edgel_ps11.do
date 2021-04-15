/*
This file is used to conduct all empirical exercises from Problem Set 11 of
Econ710q2.

Date created:  15 Apr 2021
Last modified: 15 Apr 2021
Author: Danny Edgel
*/
set more off
//capture file close table23_8


// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ710\Q2\problem_sets\PS11"


/*
	Exercise 25.15
*/

// open the March 2009 CPS data set
use cps09mar, clear

// generate dummy variable for Black individuals
g black = (race == 2)

// label variables for the table
lab var union "$\one{union}$"
lab var black "$\one{Black}$"
lab var hisp "$\one{Hispanic}$"


// estimate probit for union status, conditioning on age, race, and education
probit union age education black hisp

// export results to tex file
outreg2 using table_25_15.tex, tex(fragment) lab replace ct("Probit")

// calculate average marginal effects and output
margins, dydx(*) post
outreg2 using table_25_15.tex, tex(fragment) lab append ct("AME")
