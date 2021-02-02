

cd "C:\Users\edgel\Google Drive\UW-Madison\s21\econ713\Q1\problem_sets\PS1"

clear
set obs 30

g i = _n

g demand_rate = 3 + .01*i if mod(i,2) == 1
g supply_rate = 3 + .01*i if mod(i,2) == 0

g rate = 3 + .01*i

gsort -rate
g demand = 0
replace demand = demand + 1 if demand_rate[_n-1] > rate & !missing(demand_rate[_n-1])
replace demand = demand + demand[_n-1] if !missing(demand[_n-1])


sort rate
g supply = 0
replace supply = supply + 1 if supply_rate[_n-1] < rate
replace supply = supply + supply[_n-1] if !missing(supply[_n-1])
edit

g transactions = floor(i/2)

tw	line demand supply rate, 													///
		title("Classroom Credit Market")										///
		xlab(3.02(0.04)3.30, format(%3.2fc)) xtitle("Interest Rate")  			///
		ylab(0(2)15,angle(horizontal) format(%7.0fc)) ytitle("")				///
		legend(region(lwidth(none)) r(1) lab(1 "Demand") lab(2 "Supply"))
		
graph export figure1.png
edit