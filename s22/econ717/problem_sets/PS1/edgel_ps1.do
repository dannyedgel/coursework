/*
	This file completes all coding tasks for problem set 1 of the first
	quarter of Econ 717
	
	Date created:  10 Feb 2022
	Last modified: 15 Feb 2022
	Author: Danny Edgel (edgel@wisc.edu)
*/
capture log c

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
// also save a set that omits "married" and "muslim"
loc X	Treated Client_Age Client_Married Client_Education HH_Income muslim 	///
		Hindu_SC_Kat
loc X2	Treated Client_Age Client_Education HH_Income Hindu_SC_Kat

// save outreg options
loc opts "tex(frag) nor noobs noas"

// save epsilon (for use in 7c and 8)
loc eps = 1e-5

// save list of files in a local macro; open all files in write mode 
loc files table2 table3 table4 q9 q10 q13 q14 q14a table6
foreach f in `files'{
	capture file close `f'
	file open `f' using `f'.tex, write replace
}

// declare temporary variables
tempvar meanvar meanvar2
		
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
predict e_lpm, resid  // save for q15

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
	if ("`m'" == "Probit") loc q9a = e(ll) // save log-likelihood for use in q9
}


// 7) Calculate mean partial derivatives for the logit and probit models
// (save results, at 3 sigdigs, for table 3)

// LPM: the result is simply the coefficient for Client_Age; saved above

// logit: simply calculate the derivative of the distribution function for each
// obs, multiplied by the coefficient for Client_Age, then take the mean
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
	replace `meanvar' = `beta_probit'*normalden(xb_probit)
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
predict pred_lpmq, xb

// write the R^2 to q9a.tex for interpreting the result from (9)
qui sigdig e(r2)

replace Client_Age = Client_Age + `eps'
forval i = 2/4{
	replace ca`i' = Client_Age^`i'
}
predict q8b, xb

replace `meanvar' = (q8b - pred_lpmq)/`eps'
sum `meanvar'
qui sigdig r(mean)
loc q8 = r(value)
	
replace Client_Age = Client_Age - `eps'
forval i = 2/4{
	replace ca`i' = Client_Age^`i'
}
	
// 9) Calculate the LRI for the new LPM and output it to a tex file
qui probit taken_new
sigdig `=1 - (`q9a'/`=e(ll)')'

file write q9 "$`=r(value)'$"

// 10) calculate correct prediction rates for the each model, using both
// 		a 50% threshold and the population rate

qui sum taken_new
loc pop_mean = r(mean)
file write q10 "`: di %4.3f `pop_mean''"
foreach m in lpm lpmq logit probit{
    qui count if 	(pred_`m' >= 0.5 & taken_new == 1) | 		///
					(pred_`m' < 0.5  & taken_new == 0)
					
    loc q10_50_`m' : di %4.3f r(N)/_N
	
    qui count if 	(pred_`m' >= `pop_mean' & taken_new == 1) | ///
					(pred_`m' < `pop_mean'  & taken_new == 0)
					
    loc q10_sm_`m' : di %4.3f r(N)/_N
}

// 11) Repeat (10) using a model that only includes some of the sample, but
//     tests on the rest of the sample
qui count if imidlineid >= 1400
foreach m in lpm lpmq logit probit{
    
	qui{
		if 		("`m'" == "lpm")  	reg taken_new `X'		if imidlineid < 1400
		else if ("`m'" == "lpmq") 	reg taken_new `X' ca*	if imidlineid < 1400
		else						`m' taken_new `X'		if imidlineid < 1400
		
		if (strpos("`m'", "lpm") > 0) 	predict pred, xb 
		else							predict pred, pr
	}
	
    qui count if   ((pred >= 0.5 & taken_new == 1) | 		///
					(pred < 0.5  & taken_new == 0))
					
    loc q11_50_`m' : di %4.3f r(N)/_N
	
    qui count if   ((pred >= `pop_mean' & taken_new == 1) | ///
					(pred < `pop_mean'  & taken_new == 0))
					
    loc q11_sm_`m' : di %4.3f r(N)/_N
	
	drop pred
}
	

// 12) Re-estimate probit model with an interaction between married and Muslim
// (also re-estimate the baseline probit, for comparison)
qui probit taken_new `X2' i.Client_Married i.muslim
outreg2 using table5.tex, replace `opts' addstat("LRI", `e(r2_p)') 	///
	ctitle("Baseline")

probit taken_new `X2' 1.Client_Married##1.muslim
outreg2 using table5.tex, append `opts' addstat("LRI", `e(r2_p)') 	///
	ctitle("Interaction Effect")

// 13) Calculate the mean finite difference for Client_Married and muslim, for 
//     both the model with the interaction term and the one without

// first calculate each variable's XB for all but the vars of interest
replace `meanvar' = _b[_cons]
foreach ivar in `X2'{
    replace `meanvar' = `meanvar' + _b[`ivar']*`ivar'
}
g `meanvar2' = `meanvar'

// now calculate the interaction effect
replace `meanvar' = 	///
	normal(`meanvar' + _b[1.Client_Married] + _b[1.muslim] + 					///
	_b[1.Client_Married#1.muslim]) - normal(`meanvar' + _b[1.muslim]) - 		///
	normal(`meanvar' + _b[1.Client_Married]) + normal(`meanvar')
	
qui sum `meanvar' 
loc interaction = r(mean)
file write q13 "$`: di %4.3f `interaction''$"
file write q14 "$`: di %5.4f r(sd)'$, or "	///
				"$`: di %3.1f abs(100*(r(sd)/r(mean)))'$\% of the mean"

// now calculate the effects for each case
foreach v in Client_Married muslim{
	replace `meanvar' = normal(`meanvar2' + _b[1.`v']) - normal(`meanvar2')
	qui sum `meanvar'    
	loc `v'_eff = r(mean)
}

loc q13i_10 : di %4.3f `Client_Married_eff'
loc q13i_01 : di %4.3f `muslim_eff'
loc q13i_11 : di %4.3f `Client_Married_eff' + `muslim_eff' + `interaction'

// repeat for the model without interaction effects
qui probit taken_new `X2' i.Client_Married i.muslim
replace `meanvar2' = _b[_cons]
foreach ivar in `X2'{
    replace `meanvar2' = `meanvar2' + _b[`ivar']*`ivar'
}
foreach v in Client_Married muslim{
	replace `meanvar' = normal(`meanvar2' + _b[1.`v']) - normal(`meanvar2')
	qui sum `meanvar'    
	loc `v'_eff = r(mean)
}

loc q13n_10 : di %4.3f `Client_Married_eff'
loc q13n_01 : di %4.3f `muslim_eff'
loc q13n_11 : di %4.3f `Client_Married_eff' + `muslim_eff'


// 14) count number of married Muslims in the sample to illustrate why the SD
//     is so low 
qui count if Client_Married == 1 & muslim == 1
file write q14a "`=r(N)', or $`: di %3.1f 100*(r(N)/_N)'$\% of the sample."

// 15) regress the squared residuals from the LPM on its covariance to determine
//  if there is any evidence of heteroskedasticity
g e2_lpm = e_lpm^2
reg e2_lpm `X'
outreg2 using table7.tex, replace `opts' ctitle(" ")


// 16) run a probit regression, allowing for heteroskedasticity in age and educ
	
hetprob taken_new `X2' i.Client_Married i.muslim, 	///
	het(Client_Age Client_Education)
outreg2 using table5.tex, append `opts' ctitle("HetProbit";"lnsigma")
	
/*
	Finishing up
*/

// format the predicted probabilities and write them to table2.tex
loc models lpm lpmq logit probit
forval i = 1/4{
	di _newline "Predicted probabilities for `: word `i' of `models'':"
	sum pred_`: word `i' of `models'', d
	foreach x in mean sd min p5 p25 p50 p75 p95 max{
		loc `x'`i' : di %4.3f r(`x')
	}
}


file write table2	///
	"\begin{tabular}{rcccc}"								_newline ///
	_tab "& LPM & Quartic LPM & Logit & Probit \\\hline"	_newline ///
	_tab "Mean & `mean1' & `mean2' & `mean3' & `mean4' \\"	_newline ///
	_tab "Std. Dev. & `sd1' & `sd2' & `sd3' & `sd4' \\"		_newline ///
	_tab "Min & `min1' & `min2' & `min3' & `min4' \\"		_newline ///
	_tab "p5 & `p51' & `p52' & `p53' & `p54' \\"			_newline ///
	_tab "p25 & `p251' & `p252' & `p253' & `p254' \\"		_newline /// 
	_tab "Median & `p501' & `p502' & `p503' & `p504' \\"	_newline ///
	_tab "p75 & `p751' & `p752' & `p753' & `p754' \\"		_newline ///
	_tab "p95 & `p951' & `p952' & `p953' & `p954' \\"		_newline ///
	_tab "Max & `max1' & `max2' & `max3' & `max4'"			_newline ///
	"\end{tabular}"
	
file write table3	///
	"\begin{tabular}{rcccc}"										_newline ///
	_tab "& LPM & Logit & Probit & Quartic LPM \\\hline"			_newline ///
	_tab "Mean Partial Effect & `q7_lpm' & `q7_logit' & - & -\\"	_newline ///
	_tab "7a) & - & - & `q7a' & - \\"								_newline ///
	_tab "7b) & - & - & `q7b' & - \\"								_newline ///
	_tab "7c) & - & - & `q7c' & `q8' \\"							_newline ///
	_tab "7d) & - & - & `q7d' & - \\"								_newline /// 
	"\end{tabular}"
	
file write table4	///
	"\begin{tabular}{rcccc}"										_newline ///
	_tab 			" & \multicolumn{2}{c}{In-Sample} 		"				 ///
					" & \multicolumn{2}{c}{Out-of-Sample} \\"		_newline ///
	_tab 			" & $\geq0.5$ 		& $\geq\hat{p}$     "				 ///
					" & $\geq0.5$ 		& $\geq\hat{p}$\\\hline"	_newline ///
	_tab "LPM 		  & `q10_50_lpm'    & `q10_sm_lpm'"						 ///
					" & `q11_50_lpm'    & `q11_sm_lpm' \\"			_newline ///
	_tab "Quartic LPM & `q10_50_lpmq'   & `q10_sm_lpmq'"		 			 ///
					" & `q11_50_lpmq'   & `q11_sm_lpmq'   \\"		_newline ///
	_tab "Logit 	  & `q10_50_logit'  & `q10_sm_logit'"			 		 ///
					" & `q11_50_logit'  & `q11_sm_logit'  \\"		_newline ///
	_tab "Probit 	  & `q10_50_probit' & `q10_sm_probit'"					 ///
					" & `q11_50_probit' & `q11_sm_probit' \\"		_newline ///
	"\end{tabular}"

file write table6 	///
	"\begin{tabular}{cc|cc}"										_newline ///
	_tab "&&\multicolumn{2}{c}{\small{Mean Finite Difference}} \\"	_newline ///
	_tab "$\one{Client\_Married}$ & $\one{Muslim}$"							 ///
		 "& w/o interaction & w/ interaction \\\hline"				_newline ///
	_tab "1 & 0 & `q13n_10' & `q13i_10' \\"							_newline ///
	_tab "0 & 1 & `q13n_01' & `q13i_01' \\"							_newline ///
	_tab "1 & 1 & `q13n_11' & `q13i_11' \\"							_newline ///
	"\end{tabular}"
	
	
// close the log and save tables
log c

foreach f in `files'{
	file close `f'
}