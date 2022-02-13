/*
	This file completes all coding tasks for problem set 1 of the first
	quarter of Econ 717
	
	Date created:  10 Feb 2022
	Last modified: 12 Feb 2022
	Author: Danny Edgel (edgel@wisc.edu)
*/
capture log c
capture file close table2
capture file close table3

/*
	Housekeeping
*/

// establish directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s22\econ717\problem_sets\PS1"

// open log file
log using edgel_ps1_log.log, replace

// add user-written functions
adopath + "C:\Users\edgel\Google Drive\Code\Stata\functions"

// open data set
use Field+et+al+%282010%29+Analysis+Sample, clear

// save the set of independent variables in a local macro
loc X	Treated Client_Age Client_Married Client_Education HH_Income muslim 	///
		Hindu_SC_Kat

// save outreg options
loc opts "tex(frag) nor noobs noas"

// save epsilon (for use in 7c and 8)
loc eps = 1e-5
		
/*
	Problems
*/

// 1) Drop observations with missing values of client age, marital status, 
//		education, or household income
drop if 	///
	missing(Client_Age) 		| ///
	missing(Client_Married) 	| ///
	missing(Client_Education)	| ///
	missing(HH_Income)
	
// 	Estimate a linear probability model with loan take-up as the dependent 
//	variable and client age, client marital status, client years of 
//	schooling, client household size, client household income, indicators 
//	for Muslim and Hindu scheduled caste, and experimental treatment status

// 2) Estimate the LPM with OLS standard errors
reg taken_new `X'
outreg2 using table1.tex, replace 	/// output estimates of coefficients and SEs
				ctitle("Q2") `opts'	//  to table 1
				
qui sigdig _b[Client_Age]
loc q7_lpm = r(value) // save for table 3

// 3) Repeat the LPM estimation with robust SEs
reg taken_new `X', robust
outreg2 using table1.tex, append ctitle("Q3") `opts'

// 4) Determine the extent to which predicted probabilities fall outside the 
// range of the dependent variable
// NOTE: predictions later summarized in more detail at the end of the .do
// file, and in table 2
predict pred_lpm, xb


// 5) Estimate the LPM with variance-weighted least squares
predict y_stdp, stdp
vwls taken_new `X', sd(y_stdp)
outreg2 using table1.tex, append ctitle("Q5") `opts'


// 6) Estimate the model using logit and probit -- output estimates to table 1
// and save both predicted probabilities and their summary statistics
// Also save the coefficient for Client_Age and the predicted depvar
foreach m in Logit Probit{
	`=lower("`m'")' taken_new `X'
	outreg2 using table1.tex, append ctitle("`m'") `opts'
	predict pred_`=lower("`m'")', pr
	predict xb_`=lower("`m'")', xb
	loc beta_`=lower("`m'")' = _b[Client_Age]
}


// 7) Calculate mean partial derivatives for the logit and probit models
// (save results, at 3 sigdigs, for table 3)

// LPM: the result is simply the coefficient for Client_Age; saved above

// logit: simply calculate the derivative of the distribution function for each
// obs, multiplied by the coefficient for Client_Age, then take the mean
	tempvar meanvar
	g `meanvar' = `beta_logit'*(exp(xb_logit)/((1 + exp(xb_logit))^2))
	sum `meanvar'
	qui sigdig r(mean)
	loc q7_logit = r(value)

// Probit: use 4 different methods

	// 7a) use dprobit
	dprobit taken_new `X'

	mat x = e(dfdx)
	qui sigdig `=x[1, `=colnumb(x, "Client_Age")']'
	loc q7a = r(value)

	// 7b) calculate "by hand"
	replace `meanvar' = `beta_probit'*normal(xb_probit)
	sum `meanvar'
	qui sigdig r(mean)
	loc q7b = r(value)
	
	// 7c) change Client_Age by a little bit, then re-calculate predicted
	//     probabilities
	qui probit taken_new `X'
	replace Client_Age = Client_Age + `eps'
	predict pred2_probit, pr 
	replace `meanvar' = (pred2_probit - pred_probit)/`eps'
	sum `meanvar'
	qui sigdig r(mean)
	loc q7c = r(value)
	
	// 7d) use the margins command
	margins, dydx(Client_Age)
	mat x = r(table)
	qui sigdig `=x[1, 1]'
	loc q7d = r(value)
	

// 8) Re-estimate LPM using a quartic for Client_Age, then calculate average
//    derivative as in 7c
replace Client_Age = Client_Age - `eps'
forval i = 2/4{
	g ca`i' = Client_Age^`i'
}

reg taken_new `X' ca*
predict q8a, xb

replace Client_Age = Client_Age + `eps'
forval i = 2/4{
	replace ca`i' = Client_Age^`i'
}
predict q8b, xb

replace `meanvar' = (q8b - q8a)/`eps'
sum `meanvar'
qui sigdig r(mean)
loc q8 = r(value)
	
/*
	Finishing up
*/

// format the predicted probabilities and write them to table2.tex
loc models lpm logit probit
forval i = 1/3{
	di _newline "Predicted probabilities for `: word `i' of `models'':"
	sum pred_`: word `i' of `models'', d
	foreach x in mean sd min p5 p25 p50 p75 p95 max{
		loc `x'`i' : di %4.3f r(`x')
	}
}


file open table2 using table2.tex, write replace
file write table2	///
	"\begin{tabular}{rccc}"							_newline ///
	_tab "& LPM & Logit & Probit \\\hline"			_newline ///
	_tab "Mean & `mean1' & `mean2' & `mean3' \\"	_newline ///
	_tab "Std. Dev. & `sd1' & `sd2' & `sd3' \\"		_newline ///
	_tab "Min & `min1' & `min2' & `min3' \\"		_newline ///
	_tab "p5 & `p51' & `p52' & `p53' \\"			_newline ///
	_tab "p25 & `p251' & `p252' & `p253' \\"		_newline /// 
	_tab "Median & `p501' & `p502' & `p503' \\"		_newline ///
	_tab "p75 & `p751' & `p752' & `p753' \\"		_newline ///
	_tab "p95 & `p951' & `p952' & `p953' \\"		_newline ///
	_tab "Max & `max1' & `max2' & `max3'"			_newline ///
	"\end{tabular}"
	
file open table3 using table3.tex, write replace
file write table3	///
	"\begin{tabular}{rcccc}"										_newline ///
	_tab "& LPM & Logit & Probit & Quartic LPM \\\hline"			_newline ///
	_tab "Mean Partial Effect & `q7_lpm' & `q7_logit' & - & -\\"	_newline ///
	_tab "7a) & - & - & `q7a' & - \\"								_newline ///
	_tab "7b) & - & - & `q7b' & - \\"								_newline ///
	_tab "7c) & - & - & `q7c' & `q8' \\"							_newline ///
	_tab "7d) & - & - & `q7d' & - \\"								_newline /// 
	"\end{tabular}"
	
// close the log and save tables 2 and 3
log c
file close table2
file close table3