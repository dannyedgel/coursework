/*
	This file completes all coding tasks for problem set 2 of the first
	quarter of Econ 717
	
	Date created:  25 Feb 2022
	Last modified: 03 Mar 2022
	Author: Danny Edgel (edgel@wisc.edu)
*/
capture log c

/*
	Housekeeping
*/

// open log file
log using edgel_ps2_log.log, replace

// add user-written functions
adopath + "C:\Users\edgel\Google Drive\Code\Stata\functions"

// open data set
use "Economics 717 Spring 2022 NSW Data", clear

// save the set of independent variables in a local macro
// also save a set that omits "married" and "muslim"
loc X1	age2 educ black hisp married nodegree
loc X2	`X1' re74 re75

// save outreg options
loc opts "tex(frag) nor noobs noas"


// save list of files in a local macro; open all files in write mode 
loc files table2 q6a q6b q7a q7b q8 table3
foreach f in `files'{
	capture file close `f'
	file open `f' using `f'.tex, write replace
}

// declare temporary variables
tempvar meanvar meanvar2

// generate age^2 variable
g age2 = age^2
		
/*
	Problems
*/

// 0: drop the obs from the PSID comparison group
drop if sample == 3


// 1: Generate an experimental impact estimate by running a regression of 
//    earnings in 1978 on the treatment dummy and broader set of covariates
reg re78 treated `X2', robust
outreg2 using table1.tex, replace `opts'

// save the experimental impact estimate
loc att = _b[treated]


// 2: drop the experimental treatment group
drop if treated == 1 & sample == 1


// 3: Estimate propensity scores using a probit model, with each set of 
// 	  covariates
g 		d = 1 if sample == 1
replace d = 0 if sample == 2
probit d `X1'
predict pscorea, pr
probit d `X2'
predict pscoreb, pr
	
	
// 4: compare the distributions of the propensity scores across subsamples

foreach s in a b{
	forval i = 0/2{
		
		// generate summary statistics
		if (`i' < 2) 	qui sum pscore`s' if sample == 2 & d == `i', d
		else			qui sum pscore`s' if sample == 1, d
		
		loc j = `i' + 1
		foreach x in mean min p5 p25 p50 p75 p95 max{
			
			// save formatted statistic in relevant macro for printing to .tex
			// file in the "finishing up" section at the end of the .do file
			loc `x'`s'`j' : di %4.3f r(`x')
			
		} // x loop
	} // i loop
} // s loop


// 5: Construct a histogram of the estimated propensity scores for the combined
// 	  experimental control group and for the CPS comparison group

tw	///
	hist pscorea if d == 0, width(0.05) color(red%30)		||	///
	hist pscorea if d == 1, width(0.05) color(blue%30)			///
		ylab(, angle(horizontal)) graphregion(color(white))		///
		xlab(0(.1)0.8) leg(off) name(pscorea_hist, replace)		///
		xtitle("Coarse Score") nodraw ytitle("Count")
tw	///
	hist pscoreb if d == 0, width(0.05) color(red%30)		||		///
	hist pscoreb if d == 1, width(0.05) color(blue%30)				///
		ylab(, angle(horizontal)) graphregion(color(white))			///
		leg(lab(1 "d=0") lab(2 "d=1") region(lstyle(none)) nobox) 	///
		xlab(0(.1)0.8) xtitle("Rich Score") ytitle("Count")			///
		name(pscoreb_hist, replace) nodraw

graph combine pscorea_hist pscoreb_hist, r(2) xcom ysize(8) iscale(0.6)	///
	graphregion(color(white))
graph export q5.png, as(png) replace

	
// 6: Construct non-experimental bias estimates for both sets of estimated 
//	  propensity scores using single nearest neighbor matching without 
//    replacement

// generate all results for each pscore
foreach s in a b{
	psmatch2 d, out(re78) p(pscore`s') n(1) com norepl
	loc att6 : di %7.1fc r(att)
	
	// count number of obs dropped
	qui count if _support == 0
	
	// output dropped obs and ATT
	if (r(N) == 1) 	loc x "is"
	else			loc x "are"
	
	file write q6`s' "Using pscore`s', `=r(N)' observation, or "			///
			"$`: di %3.1f `=100*r(N)/_N''$\%, `x' dropped, and the average"	///
			" treatment on the treated is $`att6'$."

	drop _*
}

// 7: repeat problem 6 with replacement
foreach s in a b{
	psmatch2 d, out(re78) p(pscore`s') n(1) com
	loc att7 : di %7.1fc r(att)
	
	// for rich score, save att and its SE for comparison in q9
	if ("`s'" == "b"){
	    loc q7_att `att7'
		loc q7_b  : di %7.1fc 100*(`q7_att' - `att')/`att'
	} 
	
	// count number of obs dropped
	qui count if _support == 0
	
	// output dropped obs and ATT
	if (r(N) == 1) 	loc x "is"
	else			loc x "are"
	
	file write q7`s' "Using pscore`s', `=r(N)' observation, or "			///
			"$`: di %3.1f `=100*r(N)/_N''$\%, `x' dropped, and the average"	///
			" treatment on the treated is $`att7'$."

	drop _pscore _treated _support _re78 _id _n1 _nn _pdif
	ren _weight _weight`s'
}

// 8: estimate the standardized difference in real earnings for 1974 and 1975

loc y = 74
foreach y in 74 75{
	qui psmatch2 d, out(re78) p(pscoreb) n(1) com
	pstest re`y', both t(d)
	loc red`y' : di %3.1f `=100*(r(meanbiasbef)-r(meanbiasaft))/r(meanbiasbef)'	
}
file write q8 "The mean standardized bias reduction from matching is "	///
		"$`red74'$\% and $`red75'$\% for real earnings in 1974 and "	///
		"1975, respectively."
		
		
// 9: Create propensity score matching estimates using the rich propensity 
// 	  scores and kernel matching with a Gaussian (normal) kernel and bandwidths 
//	  of 0.02, 0.2 and 2.0. Impose common support condition

loc bws 0.02 0.2 2.0

forval i = 1/`: word count `bws''{
    loc bw : word `i' of `bws'
	
    psmatch2 d, out(re78) p(pscoreb) com k(normal) bw(`bw')
	
	loc q9_att`i'	: di %7.1fc r(att)
	loc q9_b`i'		: di %7.1fc 100*(`q9_att`i'' - `att')/`att'
}


// 10: Repeat Q9 using local linear, instead of kernel, matching 

forval i = 1/`: word count `bws''{
    loc bw : word `i' of `bws'
	
    psmatch2 d, out(re78) p(pscoreb) com llr bw(`bw')
	
	loc q10_att`i'	: di %7.1fc r(att)
	loc q10_b`i'	: di %7.1fc 100*(`q10_att`i'' - `att')/`att'
}

// 11: Obtain an estimate of the bias by estimating a linear regression of real 
// earnings in 1978 on the variables in the rich propensity scores and a 
// treatment dummy using all of the observations in the control group and the 
// CPS comparison group

reg re78 d `X2', robust

loc q11_att : di %7.1fc _b[d]
loc q11_b   : di %7.1fc 100*(`q11_att`i'' - `att')/`att'


// 12: Obtain an estimate of the bias by estimating a linear regression of real 
//	   earnings in 1978 on the variables in the rich propensity scores using 
//	   only the untreated observations. Use the predicted values from this 
//	   regression, evaluated at the covariate values associated with each 
//	   treated unit, as the estimated expected counterfactual outcomes for the 
//	   treated units

reg re78 `X2' if d == 0, robust
predict re78_cf, xb
g q12_att = re78 - re78_cf 
reg q12_att if d == 1

loc q12_att : di %7.1fc _b[_cons]
loc q12_b   : di %7.1fc 100*(`q12_att`i'' - `att')/`att'


// 13: Use inverse probability weighting, both by having the weights sum to one
// 	   and not

qui count if d == 1
loc n1 = r(N)
qui count if d == 0
loc n0 = r(N)
qui sum d 
loc P = r(mean)

g eff = (1/`n1')*re78*d - 	///
	(1/`n0')*((1-`P')/`P')*((pscoreb*re78*(1-d))/(1-pscoreb))
qui sum eff
loc q13_att2 : di %7.1fc r(sum)
loc q13_b2   : di %7.1fc 100*(`q13_att2' - `att')/`att'

g den = (pscoreb*(1-d))/((`n0')*(1-pscoreb))
qui sum den
loc den = r(sum)

replace eff = (1/`n1')*re78*d - 	///
	(1/`n0')*((`den')^(-1))*((pscoreb*re78*(1-d))/(1-pscoreb))
qui sum eff
loc q13_att1 : di %7.1fc r(sum)
loc q13_b1   : di %7.1fc 100*(`q13_att1' - `att')/`att'

/*
	Finishing up
*/


file write table2	///
	"\begin{tabular}{rcccccc}"										_newline ///
	_tab "&\multicolumn{2}{c}{CPS Untreated} "								 ///
			"&&\multicolumn{2}{c}{Experimental Control} \\"			_newline ///
	_tab "& $ pscorea$ & $ pscoreb$ && $ pscorea$"							 ///
			" & $ pscoreb$ \\\cline{2-3}\cline{5-6}"				_newline
	
foreach x in Mean Min p5 p25 p50 p75 p95 Max{
	loc tex "`x'"
	loc x `=lower("`x'")'
	
	forval j = 1(2)3{
		foreach s in a b{
			loc tex "`tex' & ``x'`s'`j''"
		}
		if (`j' < 3) loc tex "`tex' &"
	}
	
	file write table2 "`tex' \\" _newline
} // x loop	
file write table2 "\end{tabular}"


file write table3	///
	"\begin{tabular}{rcc}"										_newline ///
	_tab " 					& ATT 		& Bias 		\\\hline"	_newline ///
	_tab "					&			&			\\"			_newline ///	
	_tab "NN w/ replacement & `q7_att'	& `q7_b'\%	\\ &&\\"	_newline ///
	_tab "Gaussian Kernel 	&			&			\\"
	
forval i = 1/`: word count `bws''{
	file write table3 _newline _tab 	///
		"\textit{bw = `: word `i' of `bws''} & `q9_att`i'' & `q9_b`i''\% \\"
}

file write table3 _newline _tab " &&\\ Local Linear &&\\"
forval i = 1/`: word count `bws''{
	file write table3 _newline _tab 	///
		"\textit{bw = `: word `i' of `bws''} & `q10_att`i'' & `q10_b`i''\% \\"	
}

file write table3 	///
	_newline _tab "&&\\ Linear Regression 	& `q11_att' 	& `q11_b'\%	\\"	///
	_newline _tab "&&\\ Q12 Estimate 		& `q12_att' 	& `q12_b'\%	\\"	///
	_newline _tab "&&\\ IPW					&				&			\\"	///
	_newline _tab "		\textit{Scaled}		& `q13_att1'	& `q13_b1'\%\\" ///
	_newline _tab "		\textit{Unscaled}	& `q13_att2'	& `q13_b2'\%"

file write table3 _newline "\end{tabular}"
	
// close the log and save tables
foreach f in `files'{
	file close `f'
}

log c
