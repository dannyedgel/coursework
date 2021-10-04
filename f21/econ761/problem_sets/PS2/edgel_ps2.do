/*
	This file conducts all of the analysis necessary for question 2 of problem
	set 2 for Econ 761
	
	Date created:  04 Oct 2021
	Last modified: 04 Oct 2021
	Author: Danny Edgel
*/

/*
	Generate simulated environment
*/
set obs 1000				// number of cities


// assign city characteristics at random

set seed 543186
gen unif = uniform()


gen N = int(unif*10+1)	// number of firms in the city

// in 500 cities (the first 500, say), firms can collude perfectly when N<=8
gen collude = (_n <= 500 & N <= 8)


// set parameters
gen F 	= 1
gen c0 	= 1
gen c1 	= .9
gen b0 	= 1
gen b1 	= 0
gen z 	= 0
gen h 	= 0

/*
	Calculate relevant indexes
*/

// Eqm Lerner index for Cournot firms
gen LernerCN = 

// Eqm Lerner index for Monopoly
gen LernerM = {insert formula as function of c0, b0, etc.}

// Herfindahl index for symmetric firms and fixed N
gen Herf = 1/N
 
// apply Lerner index based on the conduct of firms in the city
gen Lerner = collude * LernerM + (1-collude)*LernerCN
 
set seed 155133
replace unif = uniform() - 0.5

gen lnLernerObs = ln(Lerner) + .1*unif