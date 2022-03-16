/*
	This file completes all coding tasks for problem set 3 of the first
	quarter of Econ 717
	
	Date created:  14 Mar 2022
	Last modified: 16 Mar 2022
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
loc opts "tex(frag) nor noobs noas nocon"


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
	outreg2 using table1.tex, replace `opts' ctitle(Q3)	///
		addtext(Year FE, No, State FE, No, Cluster, No) keep(mlda21)

	
// 4) repeat (3) by separately adding state and year FE's

	reg rate18_20ht mlda21 i.year, robust
	outreg2 using table1.tex, append `opts' keep(mlda21) ctitle(Q4a)	///
		addtext(Year FE, Yes, State FE, No, Cluster, No)
		

	reg rate18_20ht mlda21 i.state, robust
	outreg2 using table1.tex, append `opts' keep(mlda21) ctitle(Q4b)	///
		addtext(Year FE, No, State FE, Yes, Cluster, No)
		

	
// 5) repeat (3) including *both* state and year FE

	reg rate18_20ht mlda21 i.state i.year, robust cl(state)
	outreg2 using table1.tex, append `opts' keep(mlda21) ctitle(Q5)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes)
		


// 6) repeat (5) without clustering by state

	reg rate18_20ht mlda21 i.state i.year, robust
	outreg2 using table1.tex, append `opts' keep(mlda21) ctitle(Q6)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, No)
		


// 7) repeat (5) but omit post-1990 data

	reg rate18_20ht mlda21 i.state i.year if year <= 1990, robust cl(state)
	outreg2 using table1.tex, append `opts' keep(mlda21) ctitle(Q7)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes)
		
	/*
		NOTE: MLDA = 21 for all states after 1987, so the change in results
			is due to Goodman-Bacon's example
	*/

// 8) perform a pre-program test as specified in the problem set

	// generate placebo 
	g placebo82 = (mldayr == 1987 & year >= 1982) 

	// identify states that always have mlda = 21
	egen min_mlda = min(mlda), by(state)
	
	// run regression with placebo on states that switch in 1987 and the 
	// always 21 states
	reg rate18_20ht placebo82 i.state i.year 	///
		if min_mlda == 21 | mldayr == 1987, robust cl(state)
	outreg2 using table1.tex, append `opts' ctitle(Q8)	sortvar(mlda21)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes) keep(placebo82)
	
	
// 9) run DiD on voluntary switch states with always-takers as a control

	reg rate18_20ht mlda21 i.state i.year if state == 21 | min_mlda == 21,	///
		robust cl(state)
		
	outreg2 using table2.tex, replace `opts' ctitle(MD Only) sortvar(mlda21)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes) keep(mlda21)

	reg rate18_20ht mlda21 i.state i.year if state == 23 | min_mlda == 21,	///
		robust cl(state)
		
	outreg2 using table2.tex, append `opts' ctitle(MI Only) sortvar(mlda21)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes) keep(mlda21)

// 10) repeat (9) with separate dummies for two post-treatment periods

	// generate separate treatment indicators
	bys state mlda (year): 	///
		g mlda21_14 = (year - year[1] <= 3 & mlda == 21 & min_mlda == 18)
	bys state mlda (year): 	///
		g mlda_later = (year - year[1] > 3 & mlda == 21 & min_mlda == 18)

	
	// run analysis
	reg rate18_20ht mlda21_14 mlda_later i.state i.year 	///
		if state == 21 | min_mlda == 21, robust cl(state)
		
	outreg2 using table2.tex, append `opts' ctitle(MD Only) sortvar(mlda21)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes) 	///
		keep(mlda21_14 mlda_later)

	reg rate18_20ht mlda21_14 mlda_later i.state i.year 	///
		if state == 23 | min_mlda == 21, robust cl(state)
		
	outreg2 using table2.tex, append `opts' ctitle(MI Only) sortvar(mlda21)	///
		addtext(Year FE, Yes, State FE, Yes, Cluster, Yes) 	///
		keep(mlda21_14 mlda_later)


/*
	Finishing up
*/

	
// close the log and save tables
foreach f in `files'{
	file close `f'
}

log c
