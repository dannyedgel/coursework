/*
This file is used to conduct all empirical exercises from Problem Set 8 of
Econ710q2.

Date created:  25 Mar 2021
Last modified: 25 Mar 2021
Author: Danny Edgel
*/
set more off


// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ710\Q2\problem_sets\PS8"


/*
	Exercise 19.9
*/

// open the Card (1995) data set
use invest1993, clear


//rename relevant variables and subset data
ren vala Q
ren inva I

lab var Q "Q"
lab var I "I"

keep if Q <= 5

// estimate with Nadaraya-Watson
lpoly I Q, ci title("Nadaraya-Watson Estimator") ylab(, angle(horizontal)) 		///
	noscatter leg(region(lwidth(none)) lab(2 "m̂(x)"))
graph export "fig19_9a.png", replace
	
	
// estimate with local-linear
lpoly I Q, deg(1) ci ylab(, angle(horizontal)) noscatter	///
	title("Local-Linear Estimator") leg(region(lwidth(none)) lab(2 "m̂(x)"))
graph export "fig19_9b.png", replace

	

/*
	Exercise 19.11
*/

// FRED quarterly data set
use FRED-QD, clear

// establish data as a time series
tsset time

// generate growth rate variable, Yt, and its lag
g Y 	= 100*((gdpc1/l.gdpc1)^4-1)
g Ylag 	= l.Y


// estimate with Nadaraya-Watson
lpoly Y Ylag, ci title("Nadaraya-Watson Estimator") ylab(, angle(horizontal)) 		///
	noscatter leg(region(lwidth(none)) lab(2 "m̂(x)"))
graph export "fig19_11b.png", replace


// estimate with local-linear
lpoly Y Ylag, deg(1) ci ylab(, angle(horizontal)) noscatter	///
	title("Local-Linear Estimator") leg(region(lwidth(none)) lab(2 "m̂(x)"))
graph export "fig19_11c.png", replace
