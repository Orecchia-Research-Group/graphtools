/////////////////////////////////////
/////////////////////////////////////
/////				/////
/////	B. SEARCH FOR EDGES	/////
/////				/////
/////////////////////////////////////
/////////////////////////////////////


************* PROCEDURES *****************************************************************************
* 
* This code implements PART B of our record linkage algorithm:
* (A) Preprocess the data (Merge together, clean, and remove invalid results)
* (B) Search for potential edges (Multiple blocking strategies, including one executed on GPUs)
* (C) Score potential edges (Weights optimized using training data are taken as an input here)
* (D) Scored edges are passed to fuzzyGraph.R for graph-based entity resolution in R.
* (E) Resolved entities are then assessed for performance in fuzzyGraph_validation.do in Stata.
* 	(Parts D and E are currently implemented using a grid search over a range of input parameters for
*	the graph-based record linkage (part D), with the choice of "best" result guided by part E.
*
******************************************************************************************************



cd "/restricted/projectnb/salabs/RecordLinkage/HIV_cohort/alg_output_files"

log using "B_search_26apr2017.log", append


timer on 2
set more off


********************************************************************************************
* B.1.  Pass "HIV_demog_EM_ID.csv" to fuzzyMatch probabilistic linkage algorithm
********************************************************************************************

* Call program manually for now.

* Katia implements this on GPUs






********************************************************************************************
* B.2.  Conduct deterministic linkage on first and last names.
********************************************************************************************

* B.2.1. MULTIPLE FIRST NAMES
* Check for >1 words in first, and link to other names with just one of each. 

*just do this for exact matches on DOB/last

use "HIV_demog_EM_ID.dta", clear

* split first name
split first, p("-" " ") gen(_f)

* reshape long
reshape long _f, i(EM_ID)
drop if missing(_f)

bys dob_year dob_month dob_day last _f EM_ID: keep if _n==1
bysort dob_year dob_month dob_day last _f: gen N_copies = _N
drop if N_copies==1

* Expand to all pairs. Thank you Nick Cox! http://www.stata.com/support/faqs/data-management/expanding-datasets-to-all-pairs/
expand N_copies
sort dob_year dob_month dob_day last _f EM_ID
by dob_year dob_month dob_day last _f EM_ID: gen j=_n
by dob_year dob_month dob_day last _f: gen double EM_ID2 = EM_ID[N_copies * j]
drop if EM_ID == EM_ID2
drop if EM_ID > EM_ID2
bys EM_ID EM_ID2: keep if _n==1
keep EM_ID EM_ID2
gen DOB_sim = 1
gen last_sim = 1
gen first_sim = .95
gen match_type = "Mult first"
save "HIV_EMID_edges_multfirst.dta", replace




* B.2.2 MULTIPLE LAST NAMES
* Check for >1 words in last, and link to other names with just one of each. 
* just do this for exact matches on DOB/first

use "HIV_demog_EM_ID.dta", clear

* split last name
split last, p("-" " ") gen(_l)

* reshape long
reshape long _l, i(EM_ID)
drop if missing(_l)

* do not match on name components that are common prefixes
drop if _l == "MC" | _l == "MAC" | _l == "LA" | _l == "LE" | _l == "VAN" | _l == "VEN" | _l == "VON" | _l == "DE" | _l == "DEN" | _l == "DER" | _l == "DA" | _l == "DU" | _l == "O" 

bys dob_year dob_month dob_day first _l EM_ID: keep if _n==1
bys dob_year dob_month dob_day first _l: gen N_copies = _N
drop if N_copies==1

* Expand to all pairs. Thank you Nick Cox! http://www.stata.com/support/faqs/data-management/expanding-datasets-to-all-pairs/
expand N_copies
sort dob_year dob_month dob_day first _l EM_ID
by dob_year dob_month dob_day first _l EM_ID: gen j=_n
by dob_year dob_month dob_day first _l: gen double EM_ID2 = EM_ID[N_copies * j]
drop if EM_ID == EM_ID2
drop if EM_ID > EM_ID2
bys EM_ID EM_ID2: keep if _n==1
keep EM_ID EM_ID2
gen DOB_sim = 1
gen last_sim = .95
gen first_sim = 1
gen match_type = "Mult last"
save "HIV_EMID_edges_multlast.dta", replace




* B.2.3 CHECK FOR FIRST / LAST INVERSION

use "HIV_demog_EM_ID.dta", clear
expand 2, gen(dup)
sort EM_ID dup
by EM_ID: replace first = last[1] if dup==1
sort EM_ID dup
by EM_ID: replace last = first[1] if dup==1
bys dob_year dob_month dob_day last first: gen inversion= _N - 1

keep if inversion==1
drop dup
bys dob_year dob_month dob_day last first: gen j = _n
reshape wide EM_ID, i(dob_year dob_month dob_day last first) j(j)
drop if EM_ID1 == EM_ID2
drop if EM_ID1 > EM_ID2
rename EM_ID1 EM_ID
keep EM_ID EM_ID2
gen DOB_sim = 1
gen last_sim = .975
gen first_sim = .975
gen match_type = "Inversion"
save "HIV_EMID_edges_inversions.dta", replace





* B.2.4 MATCHING WITH NICKNAMES
* exact-matching on last, DOB

** NICKNAMES

	/* *create stata nicknames file
	insheet using "../Nicknames_HE2RO.csv", clear
	drop match
	gen match = busimatch if busimatch == patmatch
	replace match = finalcall if busimatch != patmatch
	replace match = 0 if match == .
	keep if match ==1
	keep i_first first
	drop if i_first == "DANIEL" & first == "DAVID"
	save "../Nicknames_HE2RO.dta", replace
	bys i_first: gen J = _N
	bys i_first: gen j = _n
	reshape wide first, i(i_first) j(j)
	rename i_first first
	save "Nicknames_HE2RO_wide.dta", replace
	*/

use "HIV_demog_EM_ID.dta", clear
joinby first using "../alg_input_files/Nicknames_HE2RO.dta", unmatched(none)
rename first orig_first
rename i_first first
rename EM_ID EM_ID2
joinby dob_year dob_month dob_day last first using "HIV_demog_EM_ID.dta"

keep EM_ID EM_ID2
egen double max_EM_ID = rmax(EM_ID EM_ID2)
egen double min_EM_ID = rmin(EM_ID EM_ID2)
replace EM_ID = min_EM_ID
replace EM_ID2 = max_EM_ID
keep EM_ID EM_ID2
gen DOB_sim = 1
gen last_sim = 1
gen first_sim = .95
gen match_type = "Nickname"
save "HIV_EMID_edges_nicknames.dta", replace




*****************************************************************************
** APPEND EM_ID EDGE LISTS, EXPAND W/GENDER + FACILITY, CREATE OTHER EDGES
*****************************************************************************

**************************************
* B.3 IMPORT GPU FUZZY_MATCH OUTPUT **
**************************************

set more off

insheet using "/restricted/projectnb/salabs/demog_output/Full_HIV_cohort_June_4_2017_part.csv", clear
rename v1 EM_ID
rename v2 EM_ID2
rename v3 DOB_sim
rename v4 first_sim
rename v5 last_sim
gen match_type = "Fuzzy"
compress
save "HIV_edges_fuzzy_d20.dta", replace

insheet using "/restricted/projectnb/salabs/demog_output/Full_HIV_cohort_June_4_2017_max_matches_results.csv", clear
rename v1 EM_ID
rename v2 EM_ID2
rename v3 DOB_sim
rename v4 first_sim
rename v5 last_sim
gen match_type = "Fuzzy"
compress
save "HIV_edges_fuzzy_max.dta", replace






*************************************************
* B.4 Append all edge lists at the EM_ID level
*************************************************

use "HIV_edges_fuzzy_d20.dta", clear

append using "HIV_edges_fuzzy_max.dta"

append using "HIV_EMID_edges_nicknames.dta"

append using "HIV_EMID_edges_inversions.dta"

append using "HIV_EMID_edges_multfirst.dta"

append using "HIV_EMID_edges_multlast.dta"

compress

encode match_type, gen(matchtype)
drop match_type

drop DOB_sim

save "HIV_EMID_edges.dta", replace /*formerly called Full275_all_edges.dta*/




*******************************************************************************************
* B.5 Expand EM_ID edges by gender and facility to get edges at the EM_ID_plus level
*******************************************************************************************

* first create smaller EM_ID_plus dataset, excluding first and last name (because we will just be searching within EM_ID's, so first_sim = last_sim = 1)
use "HIV_EM_ID_plus.dta", clear
drop last first test_year test_month test_day
compress
save "HIV_EM_ID_plus_small.dta", replace

	
* expand within EM_ID's
use "HIV_EMID_edges.dta", clear
count
*188,953,879
joinby EM_ID using "HIV_EM_ID_plus_small.dta"

local vars "EM_ID  gender dob_year dob_month dob_day age_flag facility_code province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'1
	}
rename EM_ID2 EM_ID
count
*244,455,015
joinby EM_ID using "HIV_EM_ID_plus_small.dta"
local vars "EM_ID  gender dob_year dob_month dob_day age_flag facility_code province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'2
	}
count
*317,057,168

drop EM_ID1 EM_ID2

compress

save "HIV_EMIDplus_edges.dta", replace /*formerly Full275_all_edges_plus.dta, renamed because we are now adding more edges*/




***************************************************************
* B.6 Exact match on first/last/DOB; fuzzy on sex/facility
***************************************************************

use "HIV_EM_ID_plus_small.dta", clear
local vars "gender dob_year dob_month dob_day age_flag facility_code province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'1
	}
joinby EM_ID using "HIV_EM_ID_plus_small.dta"
local vars "gender dob_year dob_month dob_day age_flag facility_code province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'2
	}
drop EM_ID
compress
drop if EM_ID_plus1 == EM_ID_plus2
count
*32,413,676
gen matchtype = 90
gen last_sim = 1 
gen first_sim = 1
save "HIV_edges_sex_facility.dta", replace




***************************************************************
* B.7 Exact match on first/last/gender/facility; fuzzy on DOB
***************************************************************

*** Add edges for links with EM_ID, but between EM_ID_plus (newly created when we expand by gender and facility
use "HIV_EM_ID_plus.dta", clear 
drop test_year test_month test_day EM_ID
sort last first gender facility_code
save "HIV_EM_ID_plus_sortLFGFac.dta", replace

local vars "dob_year dob_month dob_day age_flag province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'1
	}

joinby last first gender facility_code using "HIV_EM_ID_plus_sortLFGFac.dta"
local vars "gender dob_year dob_month dob_day age_flag facility_code province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'2
	}

drop if EM_ID_plus1 == EM_ID_plus2

gen gender1 = gender2
gen facility_code1 = facility_code2
drop first last
compress
count
*62,579,774
gen matchtype = 91
gen last_sim = 1 
gen first_sim = 1
save "HIV_edges_DOB.dta", replace



***************************************************************
* B.8 Exact match on DOB/last/gender/facility; fuzzy on first
***************************************************************

* fuzzy on first
use "HIV_EM_ID_plus.dta", clear 
drop test_year test_month test_day EM_ID
sort dob_year dob_month dob_day last gender facility_code
save "HIV_EM_ID_plus_sortDLGFac.dta", replace

local vars "first age_flag province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'1
	}

joinby dob_year dob_month dob_day last gender facility_code using "HIV_EM_ID_plus_sortDLGFac.dta"
local vars "gender first dob_year dob_month dob_day age_flag facility_code province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'2
	}

drop if EM_ID_plus1 == EM_ID_plus2	
gen gender1 = gender2
gen dob_year1 = dob_year2
gen dob_month1 = dob_month2
gen dob_day1 = dob_day2
gen facility_code1 = facility_code2
drop last
compress
count
*10,586,682
gen matchtype = 92
*JW score first names
jarowinkler first1 first2, gen(first_sim)

******* AUGMENT JARO-WINKLER WITH THESE ADDITIONAL AD HOC RULES FOR SHORT NAMES (OFTEN INITIALS)
	
	*clean out Afrikaans titles (**let's move this up above to pre-procesing in next version)
	local title "MEV MRS MR MISS MNR AUPA OUPA"
	foreach t of local title {
		replace first1 = "" if first1 == "`t'"
		replace first2 = "" if first2 == "`t'"
		}

	*recode first_sim for short names (values are ad hoc)
	replace first1 = subinstr(first1," ","",.) if length(first1) <=3
	replace first2 = subinstr(first2," ","",.) if length(first2) <=3
	replace first_sim = 0.85 if length(first1) == 1 & strpos(first2,first1)==1
	replace first_sim = 0.85 if length(first2) == 1 & strpos(first1,first2)==1
	replace first_sim = 0.85 if length(first1) == 2 & strpos(first2,substr(first1,1,1))==1
	replace first_sim = 0.85 if length(first2) == 2 & strpos(first1,substr(first2,1,1))==1
	replace first_sim = 0.9 if length(first1) == 2 & strpos(first2,substr(first1,1,1))==1 & strpos(word(first2,2),substr(first1,1,1))==1 
	replace first_sim = 0.9 if length(first2) == 2 & strpos(first1,substr(first2,1,1))==1 & strpos(word(first1,2),substr(first2,1,1))==1 
	replace first_sim = 0.85 if length(first2)==2 & length(first1)==2 & substr(first1,1,1)==substr(first2,2,1) & substr(first1,2,1)==substr(first2,1,1)
	replace first_sim = 0.85 if strpos(first1,first2) & length(first2)>=2 & first_sim==.
	replace first_sim = 0.85 if strpos(first2,first1) & length(first1)>=2 & first_sim==.

	*code the other "missings" as zeroes
	replace first_sim = 0 if first_sim==.
	
	*code as missing if firstNA
	replace first_sim = . if first1 == "" | first2 == ""

*drop identifiers
drop first1 first2
gen last_sim = 1 
save "HIV_edges_first.dta", replace





***************************************************************
* B.9 Exact match on DOB/first/gender/facility; fuzzy on last
***************************************************************

* fuzzy on last
use "HIV_EM_ID_plus.dta", clear 
drop test_year test_month test_day EM_ID
sort dob_year dob_month dob_day first gender facility_code
save "HIV_EM_ID_plus_sortDFGFac.dta", replace

local vars "age_flag last province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'1
	}

joinby dob_year dob_month dob_day first gender facility_code using "HIV_EM_ID_plus_sortDFGFac.dta"
local vars "last gender dob_year dob_month dob_day age_flag facility_code province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'2
	}

drop if EM_ID_plus1 == EM_ID_plus2
gen gender1 = gender2
gen dob_year1 = dob_year2
gen dob_month1 = dob_month2
gen dob_day1 = dob_day2
gen facility_code1 = facility_code2
drop first
compress
count
*6,840,290
gen matchtype = 93
*JW score first names
jarowinkler last1 last2, gen(last_sim)

******* AUGMENT JARO-WINKLER WITH ADDITIONAL AD HOC RULES FOR LAST NAMES WHERE ONE INCLUDES THE OTHER
	replace last_sim = .85 if strpos(last1,last2) & length(last2)>=3 & last_sim==.
	replace last_sim = .85 if strpos(last2,last1) & length(last1)>=3 & last_sim==.
	replace last_sim = 0 if last_sim == .

drop last1 last2
gen first_sim = 1 
save "HIV_edges_last.dta", replace




***************************************************************
* B.10 Exact match on last/DOB/sex/facility; first missing
***************************************************************

*get count of EM_ID_plus in the non-missing data
use EM_ID_plus using "HIV_EM_ID_plus.dta", clear
count
*59,124,566
scalar define N_EM_ID_plus = r(N)

*define EM_ID_plus so they keep counting upwards
use "HIV_demog_all_firstNA.dta", clear
drop if num == 1
drop num episode_no
cap drop EM_ID_plus
bys dob_year dob_month dob_day last first gender facility_code: keep if _n ==1 
gen double EM_ID_plus = _n + N_EM_ID_plus
count
*511,528 
scalar define N_EM_ID_plus_firstNA = r(N)
save "HIV_EM_ID_plus_firstNA.dta", replace

	preserve
	merge 1:m dob_year dob_month dob_day last first gender facility_code using "HIV_demog_all_firstNA.dta" 
	/*  Result                           # of obs.
	    -----------------------------------------
	    not matched                             0
	    matched                           996,900  (_merge==3)
	    -----------------------------------------
	*/
	drop _merge
	save "HIV_demog_all_firstNA.dta", replace
	restore

*rename vars
local vars "age_flag province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'1
	}

joinby dob_year dob_month dob_day last gender facility_code using "HIV_EM_ID_plus_sortDLGFac.dta"
local vars "gender dob_year dob_month dob_day age_flag facility_code province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'2
	}
drop if EM_ID_plus1 == EM_ID_plus2
gen dob_year1 = dob_year2
gen dob_month1 = dob_month2
gen dob_day1 = dob_day2
gen gender1 = gender2
gen facility_code1 = facility_code2
gen p_YOB1 = p_YOB2
gen p_first1 = .
gen p_last1  = p_last2
gen p_gender1 = p_gender2
gen p_facility1 = p_facility2
gen p_province1 = p_province2

drop last
compress
count
*95,306
gen matchtype = 94
gen last_sim = 1
gen first_sim = . 
save "HIV_edges_firstNA.dta", replace




***************************************************************
* B.11 Exact match on last/first/sex/facility; DOB missing
***************************************************************

*define EM_ID_plus so they keep counting upwards
use "HIV_demog_all_dob1800.dta", clear
drop if num == 1
drop num episode_no
cap drop EM_ID_plus
bys dob_year dob_month dob_day last first gender facility_code: keep if _n ==1 
gen double EM_ID_plus = _n + N_EM_ID_plus + N_EM_ID_plus_firstNA
save "HIV_EM_ID_plus_dob1800.dta", replace

	preserve
	merge 1:m dob_year dob_month dob_day last first gender facility using "HIV_demog_all_dob1800.dta" 
	/*  Result                           # of obs.
	    -----------------------------------------
	    not matched                             0
	    matched                         4,192,130  (_merge==3)
	    -----------------------------------------
	*/
	drop _merge
	save "HIV_demog_all_dob1800.dta", replace
	restore
	
*rename vars
local vars "dob_year dob_month dob_day age_flag province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'1
	}

joinby last first gender facility_code using "HIV_EM_ID_plus_sortLFGFac.dta"
local vars "gender dob_year dob_month dob_day age_flag facility_code province p_YOB p_first p_last p_gender p_facility p_province EM_ID_plus"
foreach v of varlist `vars' {
	rename `v' `v'2
	}
drop if EM_ID_plus1 == EM_ID_plus2
drop last first
gen gender1 = gender2
gen facility_code1 = facility_code2
gen p_YOB1 = .
gen p_first1 = p_first2
gen p_last1  = p_last2
gen p_gender1 = p_gender2
gen p_facility1 = p_facility2
gen p_province1 = p_province2
compress
count
*3,349,716
gen matchtype = 95
gen last_sim = 1
gen first_sim = 1 
save "HIV_edges_dob1800.dta", replace


*clean up
rm "HIV_EM_ID_plus_sortLFGFac.dta"
rm "HIV_EM_ID_plus_sortDLGFac.dta"
rm "HIV_EM_ID_plus_sortDFGFac.dta"





*******************************************************************
* B.12 Append all EM_ID_plus edges
*******************************************************************

* Merge back in with EMIDplus_edges

use "HIV_EMIDplus_edges.dta", clear
append using "HIV_edges_sex_facility.dta"
append using "HIV_edges_DOB.dta"
append using "HIV_edges_first.dta"
append using "HIV_edges_last.dta"
append using "HIV_edges_firstNA.dta"
append using "HIV_edges_dob1800.dta"

* count all edges
count
*432,922,612 edges

compress

label define matchtype 90 "Fuzzy_sex_fac", add
label define matchtype 91 "Fuzzy_DOB", add
label define matchtype 92 "Fuzzy_first", add
label define matchtype 93 "Fuzzy_last", add
label define matchtype 94 "Fuzzy_firstNA", add
label define matchtype 95 "Fuzzy_DOB1800", add
label values matchtype matchtype

tab matchtype
/*    matchtype |      Freq.     Percent        Cum.
--------------+-----------------------------------
        Fuzzy |310,469,082       71.71       71.71
    Inversion |    471,299        0.11       71.82
   Mult first |  2,100,920        0.49       72.31
    Mult last |     55,924        0.01       72.32
     Nickname |  3,959,943        0.91       73.24
Fuzzy_sex_fac | 32,413,676        7.49       80.72
    Fuzzy_DOB | 62,579,774       14.46       95.18
  Fuzzy_first | 10,586,682        2.45       97.62
   Fuzzy_last |  6,840,290        1.58       99.20
Fuzzy_firstNA |     95,306        0.02       99.23
Fuzzy_DOB1800 |  3,349,716        0.77      100.00
--------------+-----------------------------------
        Total |432,922,612      100.00
*/

save "HIV_all_edges_plus.dta", replace




*******************************************************************
* B.13 Create complete list of vertices
*******************************************************************

*** Create EM_ID_plus list including missings
use "HIV_EM_ID_plus.dta", clear
append using "HIV_EM_ID_plus_firstNA.dta"
append using "HIV_EM_ID_plus_dob1800.dta"
count 
*62,808,665
isid EM_ID_plus
save "HIV_EM_ID_plus_ALL.dta", replace

*** VERTICES (Just EM_ID_plus)
keep EM_ID_plus
compress
save "HIV_EMIDplus_vertices.dta", replace

*** Check that all edges are between vertices represented in vertex list
rename EM_ID_plus EM_ID_plus1
merge 1:m EM_ID_plus1 using "HIV_all_edges_plus.dta", keepusing(EM_ID_plus2)
rename EM_ID_plus2 EM_ID_plus
*merge2 should be zero
/*    Result                           # of obs.
    -----------------------------------------
    not matched                    13,517,907
        from master                13,517,907  (_merge==1)
        from using                          0  (_merge==2)

    matched                       432,922,612  (_merge==3)
    -----------------------------------------*/
drop _merge
merge m:1 EM_ID_plus using "HIV_EMIDplus_vertices.dta"
drop _merge
*merge1 should be equal to the prior merge1
/*
    Result                           # of obs.
    -----------------------------------------
    not matched                    27,951,108
        from master                13,517,907  (_merge==1) EM_ID_plus that do not appear as EM_ID_plus1
        from using                 14,433,201  (_merge==2) EM_ID_plus that do not appear as EM_ID_plus2

    matched                       432,922,612  (_merge==3) EM_ID_plus that do appear as either EM_ID_plus1 or  EM_ID_plus that do not appear as EM_ID_plus2
    -----------------------------------------
*/


*** Create file linking EM_ID_plus to Episode_no's
use episode_no EM_ID_plus using "HIV_demog_all_clean.dta", clear
append using "HIV_demog_all_dob1800.dta", keep(episode_no EM_ID_plus)
append using "HIV_demog_all_firstNA.dta", keep(episode_no EM_ID_plus)
merge m:1 EM_ID_plus using "HIV_EMIDplus_vertices.dta" /*confirm that all match*/
drop _merge
compress
save "HIV_EMIDplus_episode.dta", replace

count /*116,569,286*/
bys EM_ID_plus: keep if _n==1
count /*62,808,665*/  

timer off 2
timer list 2

log close
