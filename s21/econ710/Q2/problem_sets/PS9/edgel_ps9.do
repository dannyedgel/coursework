/*
This file is used to conduct all empirical exercises from Problem Set 9 of
Econ710q2.

Date created:  01 Apr 2021
Last modified: 01 Apr 2021
Author: Danny Edgel
*/
set more off


// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ710\Q2\problem_sets\PS9"

// choose number of grid points for plotting regression function
loc n = 10001


/*
	Exercise 20.11
*/
/*
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

*/

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
	rarea y_high y_low debtm1,  col(gs10) lc(gs8)	||	///
	line yhat debtm1, lc(blue)							///
		title("Estimated Effect of Debt-to-GDP Ratio on Growth") note("`note'") ///
		subtitle("Linear Spline at 60%") xlab(0(20)120) leg(off) 				///
		xtitle("Debt, % of GDP (last period)") ylab(`ylabs', angle(horizontal)) ///
		xline(60, lp(dash) lc(red) lw(thin)) yline(0, lc(black))
graph export fig20_15b.png, replace 

restore
	
	
	

	