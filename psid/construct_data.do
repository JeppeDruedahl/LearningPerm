clear
cap log close
set mem 36m
set matsize 800
set more off

*TABLE 6 *
/*BASELINE:	010010*/	/*Whole sample*/
/*COL. 2:	010011*/ 	/*College*/
/*COL. 3:	010012*/ 	/*No college*/
/*COL. 4:	010013*/ 	/*Born 1940s*/
/*COL. 5:	010014*/	/*Born 1930s*/
*TABLE 7 *
/*COL. 2:	200010*/ 	/*y=earnings*/
/*COL. 3:	400010*/ 	/*y=male earnings*/
*TABLE 8 *
/*COL. 2:	810010*/  	/*y excludes help*/
/*COL. 6:	010110*/ 	/*Add SEO*/

version 6.0
				scalar incdef  =0		/* 0=earnings+transfers, 1=earnings+transfer+financial income		*/
								/* 2=earnings, 3=earnings+financial income, 4=male earnings			*/
								/* 5=male+female earnings, 6=male earnings net of own taxes			*/
								/* 7=avhy, 8=income net of intrerg. transfers					*/
				scalar aftertax=1		/* 0=income before tax, 1=income after tax 					*/
				scalar totexp  =0		/* 0=non dur cons, 1=tot exp, 2=food, 3=non dur cons+ services from dur */
								/* 4=tot exp+services from dur							*/
				scalar SEO     =0		/* 0=representative sample, 1=representative sample+SEO, 2=Just SEO 	*/
				scalar agesel  =1		/* 0=20-65, 1=30-65, 2=25-65 3=30-60 [OURS]								*/
				scalar group   =0		/* 0=whole sample, 1=no college, 2= college, 3=born 1940s, 4=born 1930s, 5=years before 1986, 6= years after 1985  */

u data3,clear

******Does not use osbervations with topcoded income or financial income or federal taxes paid******
replace asset=. if trunca==1
replace y    =. if truncy==1
replace ftax =. if ftax==99999 & year>=1982 & year<=1985

gen ratio_ya=(y-asset)/y
replace ratio_ya=0 if ratio_ya<0			/*Ratio of non-financial income to income*/
								/*We have total federal taxes, ftax*/
								/*We assume that taxes paid on non-financial income*/
								/*are equal to (ratio_ya*ftax)*/

if SEO==0		{keep if seo==0}      
else if SEO==1	{keep if seo==0|seo==1}
else			{keep if seo==1}

if agesel==0	{drop if age<20|age>65}      
else if agesel==1	{drop if age<30|age>65}
else if agesel==2	{drop if age<25|age>65}
else if agesel==2	{drop if age<25|age>60}
else if agesel==3	{drop if age<30|age>60}

ren ncomp fsize

*Add a CPI price index and transform variables in real term
*This series is drawn from http://woodrow.mpls.frb.fed.us/economy/calc/hist1913.html
*CPI Base year is chained; 1982-1984 = 100
*2000 U.S. Department Of Labor, Bureau of Labor Statistics, Washington, D.C. 20212
*Consumer Price Index All Urban Consumers - (CPI-U) U.S. city average All items 1982-84=100

gen     price=1.403 if year==1992
replace price=1.362 if year==1991
replace price=1.307 if year==1990
replace price=1.240 if year==1989
replace price=1.183 if year==1988
replace price=1.136 if year==1987
replace price=1.096 if year==1986
replace price=1.076 if year==1985
replace price=1.039 if year==1984
replace price=0.996 if year==1983
replace price=0.965 if year==1982
replace price=0.909 if year==1981
replace price=0.824 if year==1980
replace price=0.726 if year==1979
replace price=0.652 if year==1978
replace price=0.606 if year==1977
replace price=0.569 if year==1976
replace price=0.538 if year==1975
replace price=0.493 if year==1974
replace price=0.444 if year==1973
replace price=0.418 if year==1972
replace price=0.405 if year==1971
replace price=0.388 if year==1970
replace price=0.367 if year==1969
replace price=0.348 if year==1968
replace price=0.334 if year==1967
replace price=1.445 if year==1993
replace price=1.482 if year==1994
replace price=1.524 if year==1995
replace price=1.569 if year==1996

gen ratio=ly/y						/*Ratio of male earnings to total income*/
replace ratio=1 if ratio>1 & ratio!=.

if       incdef==0	{replace     y=(y-asset)/price}
else	if incdef==1	{replace     y=(y)/price}
else	if incdef==2	{replace     y=(y-trhw-troth-asset)/price}
else	if incdef==3	{replace     y=(y-trhw-troth)/price}
else	if incdef==4	{replace     y=ly/price}
else	if incdef==5	{replace     y=(ly+wly)/price}
else	if incdef==6	{replace     y=(ly-ratio*ftax)/price}
else	if incdef==7	{replace     y=log(avhy)/price}
else	if incdef==8	{replace     y=(y-asset-help)/price}

replace    ly=ly   /price
replace   wly=wly  /price
replace  ftax=ftax /price

replace lc    =lc-ln(price)
replace lC    =lC-ln(price)
replace lcs   =lcs-ln(price)
replace lCs   =lCs-ln(price)


if totexp==0    		{drop lC lcs lCs}      
else if totexp==1       {drop lc lcs lCs
                 		 ren lC lc}
else if totexp==2		{drop lc lC lcs lCs
				 gen totf=food+fout+fstmp
				 gen lc=ln(totf/pf)}
else if totexp==3		{drop lc lC lCs
				 ren lcs lc}
else if totexp==4		{drop lc lC lcs
				 ren lCs lc}

if      aftertax==0  			{gen logy=ln(y)}      
else if aftertax==1 & incdef==0  	{gen logy=ln((y-ratio_ya*ftax))}      
else if aftertax==1 & incdef==1  	{gen logy=ln((y-ftax))}      
else if aftertax==1 & incdef==6  	{gen logy=ln(y)}      
else if aftertax==1 & incdef==8  	{gen logy=ln((y-ratio_ya*ftax))}      

**In the experiments in Table 7, col.2/3 the "artificial measures of income" we consider are "before taxes" **

tab year,gen(yrd)

gen black =race==2
gen other =race>=3
gen selfe =self==3
gen empl  =empst==1
gen unempl=empst==2
gen retir =empst==3
tab state,gen(stated)
gen veteran   =vet==1
gen disability=disab==1
gen kidsout   =outkid==1
gen bigcity=smsa==1|smsa==2
tab yb,gen(ybd)
tab fsize,gen(fd)
tab kids,gen(chd)
replace smsa=6 if smsa>6
tab smsa,gen(cityd)

gen extra=(tyoth)>0
gen ybw=year-agew
tab ybw,gen(ybwd)
tab age,gen(aged)

qui sum age
scalar age_min = r(min)
scalar age_max = r(max)

/* Our slighty changed version of this cleaning to include again age effects */

#delimit;
xi	i.educ*i.year i.white*i.year i.black*i.year i.other*i.year i.kidsout*i.year i.region*i.year i.extra*i.year
	i.year*i.bigcity i.year*i.kids i.fsize*i.year i.empl*i.year i.unempl*i.year i.retir*i.year;
#delimit cr

/*Controls present:	year dummies, AGE dummies, education dummies, race dummies,# of kids dummies, #of family member dummies,
				empl. status dummies, dummy for income recipient other than h/w, region dummies, dummy for big city, dummy for kids not in FU*/
/*Interactions present: 1)educ*year 2)race dummies*year 3)empl. status dummies*year 4)region*year 5)big city*year*/	

local age_min = age_min
local age_max = age_max

cap drop logc
g logc  = lc
global VARS logy logc
set more off
foreach var in $VARS {
	cap drop  `var'_res
		
		* regression
		qui reg `var' aged* yrd* edd*  white black other fd* chd* empl unempl retir extra regd* bigcity kidsout Iey_* Iwy_* Iby_* Ioy_* Ieyb* Iuy_*  Irya*  Iry_* Iyb_*
		
		predict `var'_res if e(sample), residual
		* add age-dummies
		cap drop agei
		g agei = 0
		forvalues a=`age_min'/`age_max' {
		local agedum = `a'-`age_min'+1
			replace agei = agei + _b[aged`agedum'] if age==`a' 
		}
		replace `var'_res = `var'_res + _b[_cons] + agei 
		
}


g income 	   = y
g income_res   = exp(logy_res)	
g cons 	   = exp(logc)			 	
g cons_res = exp(logc_res)

drop I*

if      group==1	{keep if educ<3}
else if group==2	{keep if educ==3}
else if group==3	{keep if yb>=1940 & yb<=1949}
else if group==4	{keep if yb>=1930 & yb<=1939}
else if group==5	{keep if year< 1986}
else if group==6	{keep if year>=1986}

sort person year
xtset person year

* save dataset for further analysis
keep person year age income* cons*
save BPPdata.dta,replace

*******************************
* save MATLAB .txt data files *
*******************************

sort person year
format income* cons* %12.5f

local group = group
quietly {
	set more off
	
		foreach var in age income_res cons_res {
			preserve
			keep person year `var'
			reshape wide `var', i(person) j(year)

			outfile `var'* using "`var'`group'.txt", comma wide replace 
			restore
		}

}