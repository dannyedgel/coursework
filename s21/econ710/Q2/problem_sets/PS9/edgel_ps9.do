/*
This file is used to conduct all empirical exercises from Problem Set 9 of
Econ710q2.

Date created:  01 Apr 2021
Last modified: 02 Apr 2021
Author: Danny Edgel
*/
set more off
capture file close table21_6
capture file close table21_8


// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ710\Q2\problem_sets\PS9"

// choose number of grid points for plotting regression function
loc n = 10001


/*
	Exercise 20.11
*/

// open the March 2009 CPS data set
use cps09mar, clear

// generate secondary variables
g wage 	= earnings / (hours*week)
g lwage	= log(wage)


// scale experience and generate polynomials
sum education
loc x_scale = r(max)			// for re-scaling later
g edu_scaled = education / r(max)

g x 	= edu_scaled
forval i = 2/6{
	g x`i' = x^(`i')
}

// (a) estimate polynomial regression
reg lwage x x2 x3 x4 x5 x6, robust

// (b) plot regression function: create new frame with x grid and estimated y
capture frame drop plotframe
frame create plotframe
frame plotframe{
	set obs `n'
	
	// generate regressor grid
	g x = _n / `n'
	forval i = 2/6{
		g x`i' = x^(`i')
	}
	
	// estimate y and se for grid
	predict yhat, xb
	predict y_se, stdp
	
	// generate c.i. bands, re-scaling y
	g y_low 	= exp(yhat - 1.96*y_se)
	g y_high	= exp(yhat + 1.96*y_se)
	replace yhat = exp(yhat)
	
	// re-scale x
	g education = x*`x_scale'
	
	// plot function
	loc note "6th-order polynomial; log of wage used as dependent variable,"	///
		"with education scaled from 0 to 1"
	mylabels 5(5)40, prefix("$") loc(ylabs)
	tw	///
		rarea y_high y_low education,  col(gs10) lc(gs8)	||	///
		line yhat education, lc(blue)							///
			title("Polynomial Regression Function") note("`note'") 				///
			subtitle("Hourly Wage on Experience") xlab(0(2)20) leg(off) 		///
			xtitle("Education") ylab(`ylabs', angle(horizontal)) 	
	graph export fig20_11b.png, replace
}



/*
	Exercise 20.15
*/

// open Reinhard and Rogoff (2010) data set
use RR2010, clear

// generate lagged variables and splines
g gdpm1 	= l.gdp
g debtm1	= l.debt
g debt_sp60	= (debtm1 - 60)*(debtm1>=60)
g debt_sp40	= (debtm1 - 40)*(debtm1>=40)
g debt_sp80	= (debtm1 - 80)*(debtm1>=80)

// label variables according to model
lab var gdp "$ Y_t$"
lab var gdpm1 "$ Y_{t-1}$"
lab var debtm1 "$ D_{t-1}$"
lab var debt_sp40 "$(D_{t-1}-40)\one{D_{t-1}\geq 40}$"
lab var debt_sp60 "$(D_{t-1}-60)\one{D_{t-1}\geq 60}$"
lab var debt_sp80 "$(D_{t-1}-80)\one{D_{t-1}\geq 80}$"

// *** (a) estimate three models

// estimate linear model, generating AIC and outputting it to regression table: 
reg gdp gdpm1 debtm1
estat ic
mat es_ic = r(S)
local AIC: display %4.1f es_ic[1,5]
outreg2 using tbl20_15.tex, replace lab tex(fragment) addstat(AIC, `AIC')

// estimate with splines
reg gdp gdpm1 debtm1 debt_sp60
estat ic
mat es_ic = r(S)
local AIC: display %4.1f es_ic[1,5]
outreg2 using tbl20_15.tex, append lab tex(fragment) addstat(AIC, `AIC')

reg gdp gdpm1 debtm1 debt_sp40 debt_sp80
estat ic
mat es_ic = r(S)
local AIC: display %4.1f es_ic[1,5]
outreg2 using tbl20_15.tex, append lab tex(fragment) addstat(AIC, `AIC')


// *** (b) plot the model with one spline


// re-estimate model
reg gdp gdpm1 debtm1 debt_sp60

// create a matrix with each debt coefficient and its ci bands, along with range
// of debt-to-gdp 
sum debtm1
loc Dmax = round(r(max),10)
mat dat = J(`Dmax'/10 + 1, 4, . )

// label columns of matrix
mat coln dat = debtm1 yhat y_low y_high 


forval D = 0(10)`Dmax'{
	di "D=`D'"
	// generate estimates of the effect of the lagged debt ratio on GDP growth for
	// D < 60
	if (`D' < 60){
		lincom _b[_cons] + _b[debtm1]*`D'
	}
	else{
		lincom _b[_cons] + _b[debtm1]*`D' + _b[debt_sp60]*(`D' - 60)
	}

	
	// fill in debt observation
	mat dat[`D'/10+1,1] = `D'
	
	// fill in debt effect and c.i.
	mat dat[`D'/10+1,2] = r(estimate)
	mat dat[`D'/10+1,3] = r(estimate) - 1.96*r(se)
	mat dat[`D'/10+1,4] = r(estimate) + 1.96*r(se)
	
}


// preserve current data, load matrix as a data set
preserve
clear
svmat dat, n(col)

// sort data on debt
sort debtm1

// plot function
loc note "Lagged growth included as a regressor"
mylabels -10(5)15, suffix("%") loc(ylabs)
tw	///
	rarea y_high y_low debtm1,  col(gs10) lc(gs10)	||	///
	line yhat debtm1, lc(blue)							///
		title("Estimated Effect of Debt-to-GDP Ratio on Growth") note("`note'") ///
		subtitle("Linear Spline at 60%") xlab(0(20)120) leg(off) 				///
		xtitle("Debt, % of GDP (last period)") ylab(`ylabs', angle(horizontal)) ///
		xline(60, lp(dash) lc(red) lw(thin)) yline(0, lc(black))
graph export fig20_15b.png, replace 

restore
	

/*
	Exercise 21.6
*/

// open Ludwig & Miller (2007) data set
use LM2007, clear
	
	
// initialize output table
file open table21_6 using table21_6.tex, write replace
file write table21_6															///
	"\begin{tabular}{lccccc}" _newline "\hline\hline"				_newline	///
		_tab " & \multicolumn{2}{c}{$ h=4$} && \multicolumn{2}{c}{ $ h=12$} \\"	///
				"\cline{2-3}\cline{5-6}"										///
		_tab " & Baseline & Covariates & &Baseline & Covariates \\"				///
				"\cline{2-6}"										_newline	 	
		
// generate RDD cutoff variables
g D		= (povrate60 >= 59.2)
g Dmx 	= (povrate60 - 59.2)*D

// store regression variables in local macros
loc Y 	mort_age59_related_postHS
loc RDD povrate60 D Dmx
loc X 	census1960_pctblack census1960_pcturban

// run baseline model and estimate theta for h = 4
reg `Y' `RDD' 		if abs(povrate60 - 59.2) <= 4, robust
loc b2 	= _b[D]
loc se2	= _se[D]

// run covariate model and estimate theta for h = 4
reg `Y' `RDD' `X' 	if abs(povrate60 - 59.2) <= 4, robust
loc b3 	 = _b[D]
loc se3	 = _se[D]
loc b13	 = _b[census1960_pctblack] 
loc se13 = _se[census1960_pctblack]
loc b23	 = _b[census1960_pcturban] 
loc se23 = _se[census1960_pcturban]


// run baseline model and estimate theta for h = 12
reg `Y' `RDD' 		if abs(povrate60 - 59.2) <= 12, robust
loc b4 	= _b[D]
loc se4	= _se[D]

// run covariate model and estimate theta for h = 12
reg `Y' `RDD' `X' 	if abs(povrate60 - 59.2) <= 12, robust
loc b5 	 = _b[D]
loc se5	 = _se[D]
loc b15	 = _b[census1960_pctblack] 
loc se15 = _se[census1960_pctblack]
loc b25	 = _b[census1960_pcturban] 
loc se25 = _se[census1960_pcturban]

// loop through estimates, formatting them
foreach i in 2 3 4 5 13 15 23 25{
	foreach x in b se{
		loc `x'`i': di %4.3f ``x'`i''
	}
}

// output results in table
file write table21_6															///
		_tab "$\hat{\theta}$ & `b2' & `b3' & &`b4' & `b5' \\"		_newline	///
		_tab "$ s(\hat{\theta})$ &(`se2') &(`se3') & & (`se4') & (`se5') \\"	///
																	_newline	///
		_tab "\% Black &  & `b13' & & & `b15' \\"					_newline	///
		_tab "$ s(\hat{\beta}_1)$ & &(`se13') & &  & (`se15') \\"	_newline	///
		_tab "\% Urban &  & `b23' & & & `b25' \\"					_newline	///
		_tab "$ s(\hat{\beta}_2)$ &&(`se23') &&& (`se25') \\\hline"	_newline	///
	"\end{tabular}"					
		
		
// close and save tex file
file close table21_6		


/*
	Exercise 21.8 -- repeat 21.6, but with different dependent variable
*/

// initialize output table
file open table21_8 using table21_8.tex, write replace
file write table21_8															///
	"\begin{tabular}{lccccc}" _newline "\hline\hline"				_newline	///
		_tab " & \multicolumn{2}{c}{$ h=4$} && \multicolumn{2}{c}{ $ h=12$} \\"	///
				"\cline{2-3}\cline{5-6}"										///
		_tab " & Baseline & Covariates & &Baseline & Covariates \\"				///
				"\cline{2-6}"		
	
loc Y mort_age25plus_related_postHS

// run baseline model and estimate theta for h = 4
reg `Y' `RDD' 		if abs(povrate60 - 59.2) <= 4, robust
loc b2 	= _b[D]
loc se2	= _se[D]

// run covariate model and estimate theta for h = 4
reg `Y' `RDD' `X' 	if abs(povrate60 - 59.2) <= 4, robust
loc b3 	 = _b[D]
loc se3	 = _se[D]
loc b13	 = _b[census1960_pctblack] 
loc se13 = _se[census1960_pctblack]
loc b23	 = _b[census1960_pcturban] 
loc se23 = _se[census1960_pcturban]


// run baseline model and estimate theta for h = 12
reg `Y' `RDD' 		if abs(povrate60 - 59.2) <= 12, robust
loc b4 	= _b[D]
loc se4	= _se[D]

// run covariate model and estimate theta for h = 12
reg `Y' `RDD' `X' 	if abs(povrate60 - 59.2) <= 12, robust
loc b5 	 = _b[D]
loc se5	 = _se[D]
loc b15	 = _b[census1960_pctblack] 
loc se15 = _se[census1960_pctblack]
loc b25	 = _b[census1960_pcturban] 
loc se25 = _se[census1960_pcturban]

// loop through estimates, formatting them
foreach i in 2 3 4 5 13 15 23 25{
	foreach x in b se{
		loc `x'`i': di %4.3f ``x'`i''
	}
}

// output results in table
file write table21_8															///
		_tab "$\hat{\theta}$ & `b2' & `b3' & &`b4' & `b5' \\"		_newline	///
		_tab "$ s(\hat{\theta})$ &(`se2') &(`se3') & & (`se4') & (`se5') \\"	///
																	_newline	///
		_tab "\% Black &  & `b13' & & & `b15' \\"					_newline	///
		_tab "$ s(\hat{\beta}_1)$ & &(`se13') & &  & (`se15') \\"	_newline	///
		_tab "\% Urban &  & `b23' & & & `b25' \\"					_newline	///
		_tab "$ s(\hat{\beta}_2)$ &&(`se23') &&& (`se25') \\\hline"	_newline	///
	"\end{tabular}"					
		
		
// close and save tex file
file close table21_8		