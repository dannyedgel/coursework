

capture file close table1


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
loc ll_wm : di %7.0f e(ll_0)
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
loc ll_km : di %7.0f e(ll_0)
loc N_km : di %7.0fc e(N)


// output table 1
loc i = 1
file open table1 using table1.tex, write replace
file write table1	///
	"\begin{tabular}{rcc}"										_newline	///
	_tab "& $\one{WalMart}$ & $\one{KMart}$ \\\hline && \\"		_newline 	///
	_tab "$\one{KMart}$ & `b`i'_wm' &             \\"			_newline 	///
	_tab "              & (`s`i'_wm') &             \\"			_newline 	///
	_tab "$\one{WalMart}$ &         & `s`i'_km'   \\"			_newline 	///
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
	_tab "Obs.				& `N_wm' & `N_km' \\"				_newline	///
	_tab "&& \\\hline"											_newline	///
	"\end{tabular}"
file close table1 




/*
	Question 2: IV Probit using distance to Benton county
*/


// probit for K-Mart entry, including instrumented Wal-Mart entry
ivprobit d_kmart `X' (d_walmart = ldist)

// Also include WalMart regression? Same as Q1, but include distance to Benton
// county?

