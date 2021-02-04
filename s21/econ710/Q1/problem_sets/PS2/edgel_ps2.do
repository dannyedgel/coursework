/*
This file is used to conduct all empirical exercises from Problem Set 2 of
Econ710q1.

Date created:  03 Feb 2021
Last modified: 03 Feb 2021
Author: Danny Edgel
*/
set more off
capture file close resultsfile
// establish problem set directory
cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ710\Q1\problem_sets\PS2"

// open CPS data
use AE80, clear

// save variables in local macros using pset names
loc x1 morekids
loc z1 samesex
loc x2 agefstm agem1 boy1st boy2nd blackd hispd othraced

/*
	3(v): Estimate reduced form regression of X_1 on Z_1 and X_2
	
	NOTE:	- X_1 = morekids
			- Z_1 = samesex
			- X_2 = agefstm agem1 boy1st boy2nd blackd hispd othraced
	
	A&E note: "Other covariates in the models are Age, Age at first birth, plus 
	indicators for Boy 1st, Boy 2nd, Black, Hispanic, and Other race"
*/

reg `x1' `z1' `x2'
outreg2 using table3v.tex, replace ct("$\one{3\text{+kids}}$") tex(fragment) lab

/*
	3(vi): Replicate Table 7, columns (1), (2), (7), and (8) in AE98
*/

// initialize LaTeX table
file open resultsfile using table3vi.tex, write replace
file write resultsfile															///
	"\begin{tabular}{lcccc}" _newline "\hline\hline"				_newline	///
		_tab "\multicolumn{5}{c}{ } \\"								  _newline 	///
		_tab "& \multicolumn{2}{c}{All women} &"								///
			 "\multicolumn{2}{c}{Husbands of married women} \\" 	  _newline 	///
		_tab "& (1) & (2) & (7) & (8) \\\hline"						  _newline 	///
		_tab "\multicolumn{5}{c}{ } \\"								  _newline 	///
		_tab "Esimation method & OLS & 2SLS & OLS & 2SLS \\"		  _newline	///
		_tab "\multicolumn{5}{c}{ } \\"								  _newline 	///
		_tab "Instrument for \textit{More than}& \textemdash"					///
			 "&\textit{Same sex} &\textemdash &\textit{Same sex} \\"   _newline	///
		_tab "\textit{ 2 children} &\multicolumn{4}{c}{ } \\"		   _newline ///
		_tab "\multicolumn{5}{c}{ } \\"								  _newline 	///
		_tab "Dependent variable: &\multicolumn{4}{c}{ } \\"

// run regressions, saving coefficients and SEs of interest as local macros

ren weeksm1 weeksm 
ren weeksd1 weeksd

loc ylist	workedm workedd weeksm weeksd hourswm hourswd incomem incomed 	///
			famincl	// dependent variables, in order
			
loc row = 1
foreach y of varlist `ylist'{
    
	// determine which parent's variable is being used
	if 		(substr("`y'",-1,1) == "d") loc p "d"
	else if (substr("`y'",-1,1) == "l") loc p "l"
	else								loc p "m"
	
	// OLS regression
	reg `y' `x2' `x1'
	
	// extract X_1 coefficient and SE
	loc b_ols`row'`p' 	= round(_b[`x1'], 0.001)
	loc se_ols`row'`p' 	= round(_se[`x1'], 0.001)
	
	// 2SLS regression
	ivregress 2sls `y' `x2' (`x1' = `z1')
	
	// extract X_1 coefficient and SE
	loc b_2sls`row'`p' 	= round(_b[`x1'], 0.001)
	loc se_2sls`row'`p' = round(_se[`x1'], 0.001)
	
	// fix issue with SE not rounding correctly
	if (length("se_2sls`row'`p'") > 8){
		loc se_2sls`row'`p' `=substr("`se_2sls`row'`p''",1,5)'
	}
	
	// if dad variable, update row
	if ("`p'" == "d") loc row = `row' + 1
	
}							

// finish LaTeX table
loc r = 1
file write resultsfile															///
		_tab "\textit{ Worked for pay} & `b_ols`r'm' & `b_2sls`r'm'"			///
			 "& `b_ols`r'd' & `b_2sls`r'd' \\"						   _newline ///
		_tab "& \small{(`se_ols`r'm')} & \small{(`se_2sls`r'm')}"				///
			 "& \small{(`se_ols`r'd')} & \small{(`se_2sls`r++'d')} \\" _newline ///
		_tab "\multicolumn{5}{c}{ } \\"								  _newline 	///
		_tab "\textit{ Weeks worked} & `b_ols`r'm' & `b_2sls`r'm'"				///
			 "& `b_ols`r'd' & `b_2sls`r'd' \\"						   _newline ///
		_tab "& \small{(`se_ols`r'm')} & \small{(`se_2sls`r'm')}"				///
			 "& \small{(`se_ols`r'd')} & \small{(`se_2sls`r++'d')} \\" _newline ///
		_tab "\multicolumn{5}{c}{ } \\"								  _newline 	///
		_tab "\textit{ Hours/Week} & `b_ols`r'm' & `b_2sls`r'm'"				///
			 "& `b_ols`r'd' & `b_2sls`r'd' \\"						   _newline ///
		_tab "& \small{(`se_ols`r'm')} & \small{(`se_2sls`r'm')}"				///
			 "& \small{(`se_ols`r'd')} & \small{(`se_2sls`r++'d')} \\" _newline ///
		_tab "\multicolumn{5}{c}{ } \\"								  _newline 	///
		_tab "\textit{ Labor income} & `b_ols`r'm' & `b_2sls`r'm'"				///
			 "& `b_ols`r'd' & `b_2sls`r'd' \\"						   _newline ///
		_tab "& \small{(`se_ols`r'm')} & \small{(`se_2sls`r'm')}"				///
			 "& \small{(`se_ols`r'd')} & \small{(`se_2sls`r++'d')} \\" _newline ///
		_tab "\multicolumn{5}{c}{ } \\"								  _newline 	///
		_tab " ln(\textit{Family income}) & `b_ols`r'l' & `b_2sls`r'l'"			///
			 "& \textemdash & \textemdash \\"						   _newline ///
		_tab "& \small{(`se_ols`r'l')} & \small{(`se_2sls`r'l')}"				///
			 "& \textemdash & \textemdash \\"			   			   _newline ///
		_tab "\multicolumn{5}{c}{ } \\\hline\hline"					   _newline ///	
	"\end{tabular}"
	
// close and save LaTex table
file close resultsfile

