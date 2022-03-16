/*
	This file completes all coding tasks for problem set 3 of the first
	quarter of Econ 717
	
	Date created:  14 Mar 2022
	Last modified: 15 Mar 2022
	Author: Danny Edgel (edgel@wisc.edu)
*/
capture log c

/*
	Housekeeping
*/

// open log file
log using edgel_ps3_log.log, replace

// add user-written functions
adopath + "C:\Users\edgel\Google Drive\Code\Stata\functions"

// open data set
use "Economics 717 Miron and Tetelbaum Data", clear

// save the set of independent variables in a local macro
// also save a set that omits "married" and "muslim"
loc X1	age2 educ black hisp married nodegree
loc X2	`X1' re74 re75

// save outreg options
loc opts "tex(frag) nor noobs noas"


// save list of files in a local macro; open all files in write mode 
loc files 
foreach f in `files'{
	capture file close `f'
	file open `f' using `f'.tex, write replace
}


		
/*
	Problems
*/

// 1) establish the data set as a panel

	xtset state year

	
// 2) create a binary treatment indicator that equals one for state-years in 
//    which the MLDA equals 21 and 0 otherwise

	g mlda21 = (mlda == 21)

	
// 3) obtain a naive estimate of the ATT with a simple reg on the treatment

	reg rate18_20ht mlda21, robust
	// outreg2 using table1.tex, replace `opts' addtext(Year FE, No, State FE, No)

	
// 4) repeat (3) by separately adding state and year FE's

	reg rate18_20ht mlda21 i.year, robust
	// outreg2 using table1.tex, append `opts' addtext(Year FE, Yes, State FE, No)	///
	//	drop(i.year)

	reg rate18_20ht mlda21 i.state, robust
	// outreg2 using table1.tex, append `opts' addtext(Year FE, No, State FE, Yes)	///
	//	drop(i.state)

	
// 5) repeat (3) including *both* state and year FE

	reg rate18_20ht mlda21 i.state i.year, robust cl(state)
	// outreg2 using table1.tex, append `opts' addtext(Year FE, Yes, State FE, Yes)	///
	//	drop(i.state i.year)


// 6) repeat (5) without clustering by state

	reg rate18_20ht mlda21 i.state i.year, robust
	// outreg2 using table1.tex, append `opts' addtext(Year FE, Yes, State FE, Yes)	///
	//	drop(i.state i.year)


// 7) repeat (5) but omit post-1990 data

	reg rate18_20ht mlda21 i.state i.year if year <= 1990, robust cl(state)
	// outreg2 using table1.tex, append `opts' addtext(Year FE, Yes, State FE, Yes)	///
	//	drop(i.state i.year)
	/*
		NOTE: MLDA = 21 for all states after 1987, so the change in results
			is due to Goodman-Bacon's example
	*/

// 8) perform a pre-program test as specified in the problem set

	// generate placebo 
	g placebo82 = (mldayr == 1987 & year >= 1982) 

	// identify states that always have mlda = 21
	egen min_mlda = min(mlda)
	
	// run regression with placebo on states that switch in 1987 and the 
	// always 21 states
	reg rate18_20ht placebo82 i.state i.year 	///
		if min_mlda == 21 | mldayr == 1987, robust cl(state)
	// outreg2 using table1.tex, append `opts' 	///
	//	addtext(Year FE, Yes, State FE, Yes) drop(i.state i.year)
	
	
// 9) run DiD on voluntary switch states with always-takers as a control

	reg rate18_20ht mlda21 i.state i.year 	///
		if inlist(state, 21, 23) == 1 | min_mlda == 21, robust cl(state)




/*
	Finishing up
*/

	
// close the log and save tables
foreach f in `files'{
	file close `f'
}

log c
