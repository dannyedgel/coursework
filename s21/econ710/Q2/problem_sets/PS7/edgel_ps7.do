/*
This file is used to conduct all empirical exercises from Problem Set 7 of
Econ710q2.

Date created:  21 Mar 2021
Last modified: 21 Mar 2021
Author: Danny Edgel
*/
set more off


// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ710\Q2\problem_sets\PS7"


/*
	Exercise 13.28
*/

// open the Card (1995) data set
use card1995, clear


// generate secondary variables, rename and/or label relevant variables
ren lwage76 lwage
ren ed76 edu
g exp = ag - edu - 6
g exp2 = (exp^2)/100
ren reg76r south
ren smsa76r urban
ren nearc4a public
ren nearc4b private

lab var edu "education"
lab var exp "experience"
lab var exp2 "experience$^2/100$"

// generate model estimates for (a)
ivregress 2sls lwage exp exp2 south black urban (edu = public private), r 
outreg2 using q13_28.tex, replace ct("(a)-2SLS") tex(fragment) lab nocons
ivregress gmm lwage exp exp2 south black urban (edu = public private), r 
outreg2 using q13_28.tex, append ct("(a)-GMM") tex(fragment) lab nocons 		///
	adds(J, e(J))


// generate additional instruments
g pub_age 		= public*age76
g pub_age_sq 	= (public*(age76^2))/100

// save expanded instrument set in local macro
loc Z public private pub_age pub_age_sq

// generate model estimates for (b)
ivregress 2sls lwage exp exp2 south black urban (edu = `Z'), r 
outreg2 using q13_28.tex, append ct("(b)-2SLS") tex(fragment) lab nocons
ivregress gmm lwage exp exp2 south black urban (edu = `Z'), r 
outreg2 using q13_28.tex, append ct("(b)-GMM") tex(fragment) lab nocons 		///
	adds(J, e(J))



/*
	Exercise 17.15
*/

// open Arellano & Bond (1991) data set
use ab1991, clear


// estimate model using Arellano-Bond one-step GMM
xtabond k, lags(1) vce(robust)
outreg2 using q17_15.tex, replace ct("(a) Arellano-Bond") tex(fragment) nocons


// estimate model using Blundell-Bond one-step GMM
xtdpdsys k, lags(1) vce(robust)
outreg2 using q17_15.tex, append ct("(b) Blundell-Bond") tex(fragment) nocons


