


cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ713\Q1\problem_sets\PS2"


tw	///
	function y = 9 - (1/5)*x, range(0 10) lc(blue)	||	///
	function y = 9 - (2/5)*x, range(0 10) lc(red) 		///
		ylab(0(2)10, angle(horizontal))					///
		xtitle(Quantity) ytitle(Price)					///
		legend( lab(1 Demand) lab(2 Marginal revenue) region(lwidth(none)) )
		
graph export 4a.png, replace




tw	///
	function y = (1/4)*x^2, range(0 6) lc(red) 			///
		ylab(0(2)10, angle(horizontal)) title(Supply)	///
		xtitle(Quantity) ytitle(Price)					///
		legend( lab(1 Supply) region(lwidth(none)) )
		
graph export 4b.png, replace


tw	///
	function y = 9 - (1/2)*x, range(0 4) lc(red)		||	///
	function y = 9 - (1/4)*4, range(4 5) lc(red)		||	///
	function y = 9 + 1/4 - (1/2)*x, range(5 10) lc(red)		///
		ylab(0(2)10, angle(horizontal)) xtitle(Quantity)	///
		ytitle(Price) title(Marginal Revenue) legend(off)
		
graph export 4c.png, replace