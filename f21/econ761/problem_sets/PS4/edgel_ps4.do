/*
	This file conducts all of the analysis necessary for question 2 of problem
	set 4 for Econ 761
	
	Date created:  09 Nov 2021
	Last modified: 16 Nov 2021
	Author: Danny Edgel
*/

capture file close table1
capture file close table2


cd "C:\Users\edgel\Google Drive\UW-Madison\f21\econ761\problem_sets\PS4\"

/*
	import data and assign variable names
*/


// rename variables
loc fnames data78 XMat 
				
loc varnames1	statefips countyfips lpop78 linc urban_pop_sh					///
				small_stores_count lsales_pop irrelevant1 irrelevant2 urban_pop ///
				inc sales pop78 pop80

loc varnames2 	cid lpop lsales_pop urban_pop_sh d_midwest ldist d_south 		///
				d_kmart d_walmart small_stores_count kmart_count_w 				///
				walmart_count_w
				
forval i = 1/2{
	qui import delimited Data/1997/`: word `i' of `fnames''.out, clear
	forval j = 1/`: word count `varnames`i'''{
		ren v`j' `: word `j' of `varnames`i'''
	}
	qui save Data/`: word `i' of `fnames'', replace
}



/*
	Question 1: run probit regressions of entry by Wal-Mart
*/

// choose covariates
loc X urban_pop_sh lpop d_midwest d_south 

// run probit for Walmart and save output
probit d_walmart 	d_kmart `X'

loc b1_wm : di %4.3f _b[d_kmart]
loc s1_wm : di %4.3f _se[d_kmart]
forval i = 1/`: word count `X''{
	loc b`=`i'+1'_wm : di %4.3f _b[`: word `i' of `X'']
	loc s`=`i'+1'_wm : di %4.3f _se[`: word `i' of `X'']
}
loc r2_wm : di %4.3f e(r2_p)
loc ll_wm : di %7.0f e(ll)
loc N_wm : di %7.0fc e(N)


// run probit for K-Mart and save output
probit d_kmart		d_walmart `X'

loc b1_km : di %4.3f _b[d_walmart]
loc s1_km : di %4.3f _se[d_walmart]
forval i = 1/`: word count `X''{
	loc b`=`i'+1'_km : di %4.3f _b[`: word `i' of `X'']
	loc s`=`i'+1'_km : di %4.3f _se[`: word `i' of `X'']
}
loc r2_km : di %4.3f e(r2_p)
loc ll_km : di %7.0f e(ll)
loc N_km : di %7.0fc e(N)


// output table 1
loc i = 1
file open table1 using table1.tex, write replace
file write table1	///
	"\begin{tabular}{rcc}"										_newline	///
	_tab "& (1) &  (2) \\"										_newline 	///
	_tab "& $\one{WalMart}$ & $\one{KMart}$ \\\hline && \\"		_newline 	///
	_tab "$\one{KMart}$ & `b`i'_wm' &             \\"			_newline 	///
	_tab "              & (`s`i'_wm') &             \\"			_newline 	///
	_tab "$\one{WalMart}$ &         & `b`i'_km'   \\"			_newline 	///
	_tab "   			  &         & (`s`i++'_km') \\"			_newline 	///
	_tab "\% Urban & `b`i'_wm' & `b`i'_km'   \\"				_newline 	///
	_tab " 		   & (`s`i'_wm') & (`s`i++'_km') \\"			_newline 	///
	_tab "$\loge{Population}$ & `b`i'_wm' & `b`i'_km' \\"		_newline 	///
	_tab " 					  & (`s`i'_wm') & (`s`i++'_km') \\"	_newline 	///
	_tab "$\one{Midwest}$ & `b`i'_wm' & `b`i'_km' \\"			_newline 	///
	_tab " 				  & (`s`i'_wm') & (`s`i++'_km') \\"		_newline 	///
	_tab "$\one{South}$ & `b`i'_wm' & `b`i'_km' \\"				_newline 	///
	_tab " 			    & (`s`i'_wm') & (`s`i++'_km') \\ && \\"	_newline 	///
	_tab "Pseudo $ R^2$ & `r2_wm' & `r2_km' \\"					_newline	///
	_tab "Log-likelihood & `ll_wm' & `ll_km' \\"				_newline	///
	_tab "Obs.				& `N_wm' & `N_km' \\\hline"			_newline	///
	"\end{tabular}"
file close table1 




/*
	Question 2: IV Probit using distance to Benton county
*/


// probit for K-Mart entry, including instrumented Wal-Mart entry
ivprobit d_kmart `X' (d_walmart = ldist), mle

loc b1_km : di %4.3f _b[d_walmart]
loc s1_km : di %4.3f _se[d_walmart]
forval i = 1/`: word count `X''{
	loc b`=`i'+1'_km : di %4.3f _b[`: word `i' of `X'']
	loc s`=`i'+1'_km : di %4.3f _se[`: word `i' of `X'']
}
loc r2_km : di %4.3f e(r2_p)
loc ll_km : di %7.0f e(ll)
loc N_km : di %7.0fc e(N)

// Also include WalMart regression? Same as Q1, but include distance to Benton
// county?
probit d_walmart 	d_kmart ldist  `X'

loc b1_wm : di %4.3f _b[d_kmart]
loc s1_wm : di %4.3f _se[d_kmart]
loc b2_wm : di %4.3f _b[ldist]
loc s2_wm : di %4.3f _se[ldist]
forval i = 1/`: word count `X''{
	loc b`=`i'+2'_wm : di %4.3f _b[`: word `i' of `X'']
	loc s`=`i'+2'_wm : di %4.3f _se[`: word `i' of `X'']
}
loc ll_wm : di %7.0f e(ll)
loc N_wm : di %7.0fc e(N)

// output table 2 (not used)
loc i = 1
file open table2 using table2.tex, write replace
file write table2	///
	"\begin{tabular}{rcc}"										_newline	///
	_tab "& (1) &  (2) \\"										_newline 	///
	_tab "& $\one{WalMart}$ & $\one{KMart}$ \\\hline && \\"		_newline 	///
	_tab "$\one{KMart}$ & `b`i'_wm' &             \\"			_newline 	///
	_tab "              & (`s`i'_wm') &             \\"			_newline 	///
	_tab "$\loge{Dist}$ & `b2_wm' &             	\\"			_newline 	///
	_tab "              & (`s2_wm') &             \\"			_newline 	///
	_tab "$\one{WalMart}$ &         & `b`i'_km'   \\"			_newline 	///
	_tab "   			  &         & (`s`i++'_km') \\"			_newline 	///
	_tab "\% Urban & `b`i'_wm' & `b`i'_km'   \\"				_newline 	///
	_tab " 		   & (`s`i'_wm') & (`s`i++'_km') \\"			_newline 	///
	_tab "$\loge{Population}$ & `b`i'_wm' & `b`i'_km' \\"		_newline 	///
	_tab " 					  & (`s`i'_wm') & (`s`i++'_km') \\"	_newline 	///
	_tab "$\one{Midwest}$ & `b`i'_wm' & `b`i'_km' \\"			_newline 	///
	_tab " 				  & (`s`i'_wm') & (`s`i++'_km') \\"		_newline 	///
	_tab "$\one{South}$ & `b`i'_wm' & `b`i'_km' \\"				_newline 	///
	_tab " 			    & (`s`i'_wm') & (`s`i++'_km') \\ && \\"	_newline 	///
	_tab "Pseudo $ R^2$ & `r2_wm' & - \\"						_newline	///
	_tab "Log-likelihood & `ll_wm' & `ll_km' \\"				_newline	///
	_tab "Obs.				& `N_wm' & `N_km' \\\hline"			_newline	///
	"\end{tabular}"
file close table2



/*
	Question 3: Bresnahan and Reiss analysis 
*/

// ordered probit on number of large players
g nlarge = d_walmart + d_kmart

oprobit nlarge `X'

forval i = 1/`: word count `X''{
		loc b`i'l : di %4.3f _b[`: word `i' of `X'']
		loc s`i'l : di %4.3f _se[`: word `i' of `X'']
}
loc lll : di %7.0f e(ll)
loc r2l : di %4.3f e(r2_p)
loc Nl  : di %7.0fc e(N)


// ordered probit on total number of players
g ntotal = nlarge + small_stores_count

oprobit ntotal `X'

forval i = 1/`: word count `X''{
		loc b`i'a : di %4.3f _b[`: word `i' of `X'']
		loc s`i'a : di %4.3f _se[`: word `i' of `X'']
}
loc lla : di %7.0f e(ll)
loc r2a : di %4.3f e(r2_p)
loc Na  : di %7.0fc e(N)



// output table 2
loc i = 1
file open table2 using table2.tex, write replace
file write table2	///
	"\begin{tabular}{rcc}"									_newline	///
	_tab "& \# Large $ & \# Total \\\hline && \\"			_newline 	///
	_tab "\% Urban & `b`i'l' & `b`i'a'   \\"				_newline 	///
	_tab " 		   & (`s`i'l') & (`s`i++'a') \\"			_newline 	///
	_tab "$\loge{Population}$ & `b`i'l' & `b`i'a' \\"		_newline 	///
	_tab " 					  & (`s`i'l') & (`s`i++'a') \\"	_newline 	///
	_tab "$\one{Midwest}$ & `b`i'l' & `b`i'a' \\"			_newline 	///
	_tab " 				  & (`s`i'l') & (`s`i++'a') \\"		_newline 	///
	_tab "$\one{South}$ & `b`i'l' & `b`i'a' \\"				_newline 	///
	_tab " 			    & (`s`i'l') & (`s`i++'a') \\ && \\"	_newline 	///
	_tab "Pseudo $ R^2$ & `r2l' & `r2a' \\"					_newline	///
	_tab "Log-likelihood & `lll' & `lla' \\"				_newline	///
	_tab "Obs.				& `Nl' & `Na' \\\hline"			_newline	///
	"\end{tabular}"
file close table2


/*
	Question 4: two-stage estimation a la Bajari et al. (2012)
*/


// K-Mart first stage
ivprobit d_kmart `X' (d_walmart = ldist), mle
predict pr_kmart, pr

loc b1_km : di %4.3f _b[d_walmart]
loc s1_km : di %4.3f _se[d_walmart]
forval i = 1/`: word count `X''{
	loc b`=`i'+1'_km : di %4.3f _b[`: word `i' of `X'']
	loc s`=`i'+1'_km : di %4.3f _se[`: word `i' of `X'']
}
loc r2_km : di %4.3f e(r2_p)
loc ll_km : di %7.0f e(ll)
loc N_km : di %7.0fc e(N)

// Wal-Mart first stage
probit d_walmart 	d_kmart ldist  `X'
predict pr_walmart, pr

loc b1_wm : di %4.3f _b[d_kmart]
loc s1_wm : di %4.3f _se[d_kmart]
loc b2_wm : di %4.3f _b[ldist]
loc s2_wm : di %4.3f _se[ldist]
forval i = 1/`: word count `X''{
	loc b`=`i'+2'_wm : di %4.3f _b[`: word `i' of `X'']
	loc s`=`i'+2'_wm : di %4.3f _se[`: word `i' of `X'']
}
loc ll_wm : di %7.0f e(ll)
loc N_wm : di %7.0fc e(N)

// K-Mart second stage
probit d_kmart pr_walmart `X',

loc b1_km2 : di %4.3f _b[pr_walmart]
loc s1_km2 : di %4.3f _se[pr_walmart]
forval i = 1/`: word count `X''{
	loc b`=`i'+1'_km2 : di %4.3f _b[`: word `i' of `X'']
	loc s`=`i'+1'_km2 : di %4.3f _se[`: word `i' of `X'']
}
loc r2_km2 : di %4.3f e(r2_p)
loc ll_km2 : di %7.0f e(ll)
loc N_km2 : di %7.0fc e(N)

// Wal-Mart second stage
probit d_walmart pr_kmart ldist  `X'

loc b1_wm2 : di %4.3f _b[pr_kmart]
loc s1_wm2 : di %4.3f _se[pr_kmart]
loc b2_wm2 : di %4.3f _b[ldist]
loc s2_wm2 : di %4.3f _se[ldist]
forval i = 1/`: word count `X''{
	loc b`=`i'+2'_wm2 : di %4.3f _b[`: word `i' of `X'']
	loc s`=`i'+2'_wm2 : di %4.3f _se[`: word `i' of `X'']
}
loc ll_wm2 : di %7.0f e(ll)
loc N_wm2 : di %7.0fc e(N)

// output table 3
loc i = 1
file open table3 using table3.tex, write replace
file write table3	///
	"\begin{tabular}{rcccc}"										_newline	///
	_tab "& \multicolumn{2}{c}{1st Stage} &" 								///
							"\multicolumn{2}{c}{2nd Stage} \\"	_newline 	///
	_tab "& (1) &  (2) & (3) &(4) \\"							_newline 	///
	_tab "& $\one{WalMart}$ & $\one{KMart}$ " 								///
			"& $\one{WalMart}$ & $\one{KMart}$\\\hline && \\"	_newline 	///
	_tab "$\Delta_i$ & `b`i'_wm' & `b`i'_km' &" 							///
				"`b`i'_wm2' & `b`i++'_km2'           \\"		_newline 	///
	_tab "$\loge{Dist}$ & `b2_wm' &  & `b2_wm2' &    	\\"		_newline 	///
	_tab "\% Urban & `b`i'_wm' & `b`i'_km' & `b`i'_wm2' & `b`i++'_km2'  \\"	///		
																_newline 	///
	_tab "$\loge{Population}$ & `b`i'_wm' & `b`i'_km' &" 					///
								"`b`i'_wm2' & `b`i++'_km2' \\"	_newline 	///
	_tab "$\one{Midwest}$ &`b`i'_wm'&`b`i'_km'&`b`i'_wm2' &`b`i++'_km2' \\"	///
																_newline 	///
	_tab "$\one{South}$ & `b`i'_wm'&`b`i'_km' &`b`i'_wm2' &`b`i++'_km2' \\"	///
																_newline 	///
	_tab "Log-likelihood & `ll_wm' & `ll_km' & `ll_wm2' & `ll_km2' \\"		///
																_newline	///
	_tab "Obs.	& `N_wm' & `N_km' & `N_wm2' & `N_km2' \\\hline"	_newline	///
	"\end{tabular}"
file close table3





