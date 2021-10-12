/*
	This file conducts all of the analysis necessary for question 2 of problem
	set 2 for Econ 761
	
	Date created:  04 Oct 2021
	Last modified: 09 Oct 2021
	Author: Danny Edgel
*/

/*
	Housekeeping
*/

clear
capture file close table1
capture file close table2


// initialize output tables
file open table1 using table1.tex, write replace
file write table1	///
	"\begin{tabular}{r|ccccc}" _newline											///
	_tab "& & & \multicolumn{2}{c}{Test: $\beta_{\loge{H}}=1$} & \\" _newline	///
	_tab " & $\beta_{\loge{H}}$ & se$\left(\beta_{\loge{H}}\right)$" _newline	///
			"& \textit{F}-score  & \textit{p}-score & N \\\hline"	 _newline
			
file open table2 using table2.tex, write replace
file write table2	///
	"\begin{tabular}{r|cc}"	_newline											///
	_tab "& $\nu\sim U[-1,1]$ & $\eta\sim U[-1,1]$	\\\hline && \\" _newline


//______________________________________________________________________________

/*
	Question 2
*/
do ps2_q2.do

//______________________________________________________________________________

/*
	Question 3
*/
do ps2_q3.do


//______________________________________________________________________________