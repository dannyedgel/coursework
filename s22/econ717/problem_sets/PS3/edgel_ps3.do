/*
	This file completes all coding tasks for problem set 3 of the first
	quarter of Econ 717
	
	Date created:  14 Mar 2022
	Last modified: 17 Mar 2022
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

// save outreg options
loc opts "tex(frag) nor noobs noas nocon"


// save list of files in a local macro; open all files in write mode 
loc files table3
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


// Bonus question 1

	// perform Goodman-Bacon decomposition
	bacondecomp rate18_20ht mlda21, ddetail nograph cl(state) robust
				
	// save results from decomposition and load them
	foreach x in b V dd wt sumdd{
		mat `x' = e(`x')
	}
	mat dat = J(`=colsof(dd)', 2, .)
	mat rown dat = `: coln dd'
	mat coln dat = dd wt
	forval r = 1/`=colsof(dd)'{
		mat dat[`r', 1] = dd[1, `r']
		mat dat[`r', 2] = wt[1, `r']
	}
	
	clear
	svmat dat, n(col)
	g category = ""
	forval i = 1/`=colsof(dd)'{
		replace category = "`: word `i' of `: coln dd''" if _n == `i'
	}
	g b = `=b[1,1]'
	g se = sqrt(`=V[1,1]')
	
	// generate comparison groups from category variable
	g first = substr(category, 1, strpos(category, "_") - 1)
	g last  = substr(category, strpos(category, "_") + 1, .)
	
	g 		group = "Treated vs. Always 21" 	if last == "Always"
	//replace group = "Timing"		if last != "Always"
	
	destring first, g(fyear) force
	destring last, g(lyear) force
	
	replace group = "Pre vs. Post" if fyear < lyear & !missing(lyear)
	replace group = "Post vs. Pre" if fyear > lyear & !missing(lyear)
	
	// plot the estimates
	loc note "Dashed lines indicate the 95% confidence interval"
	loc note "`note' for the TWFEDD estimate"
	tw	///
		scatter dd wt if group == "Treated vs. Always 21", m(t)	||	///
		scatter dd wt if group == "Pre vs. Post", m(oh)			||	///
		scatter dd wt if group == "Post vs. Pre", m(x)				///
			graphregion(color(white)) ylab(-15(5)15, angle(horizontal))		///
			xlab(0(0.1)0.6) title("2x2 DiD Estimates vs. TWFEDD Weight")	///
			xtitle("Weight") ytitle("") note("`note'")				///
			yline(`=b[1]', lc(red)) yline(0, lw(thin) lc(black))	///
			yline(`=b[1] + 1.96*se[1]', lp(dash) lc(red) lw(vthin))	///
			yline(`=b[1] - 1.96*se[1]', lp(dash) lc(red) lw(vthin))	///
			text(`=b[1] + se[1]' 0.3 								///
				"TWFEDD Estimate = `=round(b[1], .01)'", 			///
				color(red) size(small))								///
			leg(lab(1 "Treated vs. Always 21") lab(2 "Pre vs. Post") 	///
				nobox r(1) lab(3 "Post vs. Pre") region(color(white)))
				
	graph export figure1.png, replace
	
	// write a table of each DD estimate by group with its weight
	egen grpwt = sum(wt), by(group)
	g dd_weighted = dd*(wt/grpwt)
	collapse (mean) dd (sum) dd_weighted wt, by(group)
	
	file write table3	///
		"\begin{tabular}{r|ccc}"									_newline ///
			_tab " &\multicolumn{2}{c}{$\hat{\beta}^{2x2}$} & \\"	_newline ///
			_tab " & Simple Avg & Weighted Avg & TWFEDD Weight"				 /// 
					"\\\hline &&\\" _newline
	
	qui count
	forval i = 1/`=r(N)'{
		file write table3 _tab "`=group[`i']' & $`: di %3.2f dd[`i']'$"	///
								"& $`: di %3.2f dd_weighted[`i']'$" 	///
								"& $`: di %4.3f wt[`i']'$ \\" _newline
	}
	
	file write table3 "\end{tabular}"
	
		
		
/*
	Finishing up
*/

	
// close the log and save tables
foreach f in `files'{
	file close `f'
}

log c
