/*
	This file conducts all of the analysis necessary for question 2 of problem
	set 2 for Econ 761
	
	Date created:  04 Oct 2021
	Last modified: 06 Oct 2021
	Author: Danny Edgel
*/

/*
	Housekeeping
*/

clear
capture file close table1


// initialize output table
file open table1 using table1.tex, write replace
file write table1	///
	"\begin{tabular}{r|ccccc}" _newline											///
	_tab " & $\beta_{\loge{H}}$ & se$\left(\beta_{\loge{H}}\right)$"			///
			"& F-score, $\beta_{\loge{H}}=1$ & "								///
			"Pr$\left(\beta_{\loge{H}} = 1\right)$ & N \\\hline"


/*
	Generate simulated environment
*/
set obs 1000				// number of cities


// assign city characteristics at random

set seed 543186
gen unif = uniform()


gen N = int(unif*10+1)	// number of firms in the city

// in 500 cities (the first 500, say), firms can collude perfectly when N<=8
gen collude = (_n <= 500 & N <= 8)


// set parameters
gen F 	= 1
gen c0 	= 1
gen c1 	= .9
gen b0 	= 1
gen b1 	= 0
gen z 	= 0
gen xi 	= 0
gen a0 	= 3
gen a1 	= 1
gen nu	= 0
gen eta = 0



// run twice; first with equation 3, then with equation 1
foreach eqn in 3 1{
	
	// print sample description
	file write table1 	///
		_newline _tab "&&&&& \\ \textit{Equation (`eqn')} & & & & & \\"
	
	/*
		Calculate relevant indexes
	*/
	
	if (`eqn' == 3){
		// Eqm Lerner index for Cournot firms
		gen LernerCN`eqn' = c1*exp(-c0-xi)

		// Eqm Lerner index for Monopoly
		gen LernerM`eqn' = c1*exp(-c0-xi)
	}
	if (`eqn' == 1){
		// Eqm Lerner index for Cournot firms
		gen LernerCN`eqn' = ((a0-b0+nu-eta)*N)/(2*(a0+nu) - (a0-b0+nu-eta)*N)

		// Eqm Lerner index for Monopoly
		gen LernerM`eqn' = (a0 - b0) / (a0 + b0 + 2*nu)
	}
	
	// Herfindahl index for symmetric firms and fixed N
	gen Herf`eqn' 	= 1/N
	gen lnHerf`eqn' = ln(Herf`eqn')
	 
	// apply Lerner index based on the conduct of firms in the city
	gen Lerner`eqn' = collude * LernerM`eqn' + (1-collude)*LernerCN`eqn'
	 
	set seed 155133
	replace unif = uniform() - 0.5

	gen lnLernerObs`eqn' = ln(Lerner`eqn') + .1*unif
	
	// generate list of sample descriptions
	loc samp_descs Cournot Collusion Pooled
	
	// run analysis both pooled and separate
	forval pooled = 0/2{
		
		// save sample description
		loc samp `: word `=`pooled'+1' of `samp_descs''

		//set trace on
		/*
			Conduct regression analysis 
		*/

		// SCP regression
		if (`pooled' == 2){
			reg lnLernerObs`eqn' lnHerf`eqn'
			loc N = e(N)
		}
		else{
			reg lnLernerObs`eqn' lnHerf`eqn' if collude == `pooled'
			loc N = e(N)
		}
		
		
		
		loc beta : 	di %4.3f _b[lnHerf`eqn']
		loc se : 	di %4.3f _se[lnHerf`eqn']
		
		test _b[lnHerf`eqn'] = 1
		
		loc F : di %7.0fc r(F)
		loc p : di %4.3f r(p)
		
		file write table1 	///
			_newline _tab "`samp' & `beta' & `se' & `F' & `p' & `N' \\"

	} // pool loop
	
}	// eqn loop


// finish and close LaTeX table
file write table1 _newline "\end{tabular}"
file close table1