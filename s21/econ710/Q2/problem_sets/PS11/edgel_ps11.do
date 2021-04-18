/*
This file is used to conduct all empirical exercises from Problem Set 11 of
Econ710q2.

Date created:  15 Apr 2021
Last modified: 18 Apr 2021
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

// subset to men
keep if female == 0

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

/*
	Exercise 25.17
*/

// re-open the March 2020 CPS data set
use cps09mar, clear

// subset to women with a college degree
keep if female == 1 & education >= 16

// generate dependent variable
g Y = (inlist(marital,1,2,3) == 1)

// generate spline variables
g age25 = (age >= 25)*(age - 25)
g age30 = (age >= 30)*(age - 30)
g age35 = (age >= 35)*(age - 35)
g age40 = (age >= 40)*(age - 40)
g age50 = (age >= 50)*(age - 50)
g age60 = (age >= 60)*(age - 60)


// estimate a probit model against age
probit Y age age25 age30 age35 age40 age50 age60

// plot probabilities -- generate empirical averages for each age group
collapse (count) n_obs = female (sum) n_married = Y, by(age)

g Y_avg = n_married / n_obs
g age25 = (age >= 25)*(age - 25)
g age30 = (age >= 30)*(age - 30)
g age35 = (age >= 35)*(age - 35)
g age40 = (age >= 40)*(age - 40)
g age50 = (age >= 50)*(age - 50)
g age60 = (age >= 60)*(age - 60)

	
	
predict Y_pred, pr
tw	///
	line Y_pred age || scatter Y_avg age, msize(tiny) mcol(black) 	///
	xlab(20(5)80) ylab(0(.1)1, angle(horizontal) grid nogextend) xtitle("Age") 	///
	plotregion(margin(tiny) lcolor(none)) ytitle("")						///
	title("Probability of Being Married: Women with College Degrees")		///
	subtitle("Fitted values, probit with linear splines")					///
	note("Splines at age 25, 30, 35, 40, 50, and 60")						///
	leg(lab(1 "Fitted values") lab(2 "Empirical proportion") 	///
			region(lstyle(none)))
graph export fig_25_17.png, replace
