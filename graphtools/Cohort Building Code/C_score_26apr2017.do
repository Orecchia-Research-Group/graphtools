/////////////////////////////////////
/////////////////////////////////////
/////				/////
/////	C. SCORE EDGES		/////
/////				/////
/////////////////////////////////////
/////////////////////////////////////


************* PROCEDURES *****************************************************************************
* 
* This code implements PART C of our record linkage algorithm:
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

log using "C_score_26apr2017.log", append


timer on 3
set more off

*
***********************************************************************************************
* C.1 Generate similarity score for each characteristic using Fellegi-Sunter formula
***********************************************************************************************

* Get m-probabilities from Training Data
use "/restricted/projectnb/salabs/RecordLinkage/man.match.data/for.graphlink/Pooled_RA_matched_Training_MMredux2_forBaichuan.dta", clear
rename facility_id1 facility_code1
rename facility_id2 facility_code2

keep if r_manmatch==1
count
gen e_first = first_sim <.9
gen e_last = last_sim <.9
gen e_dob_Yg10 = dob_year1 != dob_year2 & abs(dob_year1 - dob_year2) >=10
gen e_dob_Yl10 = dob_year1 != dob_year2 & abs(dob_year1 - dob_year2) <10
gen e_dob_M = dob_month1 != dob_month2 if age_flag1 != 1 & age_flag2 != 1
gen e_dob_D = dob_day1 != dob_day2 if age_flag1 != 1 & age_flag2 != 1
gen e_gender = gender1 != gender2 if gender1 != "U" & gender2 != "U"
gen e_facility = facility_code1 != facility_code2 if facility_code2 <. & facility_code1 <.
gen e_province = province1 != province2

/*
su e_*
    Variable |       Obs        Mean    Std. Dev.       Min        Max
-------------+--------------------------------------------------------
     e_first |      3926    .1164035     .320749          0          1
      e_last |      3926    .0438105    .2046993          0          1
  e_dob_Yg10 |      3926     .005349      .07295          0          1
  e_dob_Yl10 |      3926    .0382068    .1917197          0          1
     e_dob_M |      3742    .0916622    .2885871          0          1
     e_dob_D |      3742    .1087654    .3113862          0          1
-------------+--------------------------------------------------------
    e_gender |      3866    .0271599    .1625701          0          1
  e_facility |      3917     .457493    .4982535          0          1
  e_province |      3926    .0509424    .2199082          0          1
*/

su e_first
	scalar define m_first = 1 - r(mean)
su e_last
	scalar define m_last = 1 - r(mean)
su e_dob_Yg10
	scalar define m_Y = 1 - r(mean)
su e_dob_M
	scalar define m_M = 1 - r(mean)
su e_dob_D
	scalar define m_D = 1 - r(mean)
su e_gender
	scalar define m_gender = 1 - r(mean)
su e_facility
	scalar define m_facility = 1 - r(mean)
su e_province
	scalar define m_province = 1 - r(mean)




use "HIV_all_edges_plus.dta", clear

** Define max(p1,p2)
gen p_gender = max(p_gender1,p_gender2)
gen p_YOB = max(p_YOB1,p_YOB2)
gen p_last = max(p_last1,p_last2)
gen p_first = max(p_first1,p_first2)
gen p_province = max(p_province1,p_province2)
gen p_facility = max(p_facility1,p_facility2)
drop p_gender1 p_gender2 p_YOB1 p_YOB2 p_last1 p_last2 p_first1 p_first2 p_province1 p_province2 p_facility1 p_facility2


**Gender
gen gender_sim_score = log(m_gender/p_gender)/log(2) if gender1==gender2
	replace gender_sim_score = log((1-m_gender)/(1-p_gender))/log(2) if gender1 != gender2
	replace gender_sim_score = 0 if gender1=="U" | gender2=="U"

drop gender1 gender2 p_gender


**Last
*last_sim_score is J-W weighted average, with 4 decay (see p133 in Herzog 2007)
gen w_last_agree = log(m_last/p_last)/log(2)
gen w_last_disagree = log( (1-m_last) / (1-p_last) )  /  log(2)
gen last_sim_score = max(w_last_disagree, w_last_agree - 4*(w_last_agree - w_last_disagree)*(1-last_sim))
replace last_sim_score = 0 if last_sim_score == .
drop p_last w_last_agree w_last_disagree


**First
*first_sim_score is J-W weighted average, with 4 decay
gen w_first_agree = log(m_first/p_first)/log(2)
gen w_first_disagree = log((1-m_first)/(1-p_first))/log(2)
gen first_sim_score = max(w_first_disagree, w_first_agree - 4*(w_first_agree - w_first_disagree)*(1-first_sim))
replace first_sim_score = 0 if first_sim_score == .
drop w_first_agree w_first_disagree p_first


**Province
gen prov_sim_score = log(m_province/p_province)/log(2) if province1==province2
	replace prov_sim_score = log((1-m_province)/(1-p_province))/log(2) if province1 != province2 & province1 != "" & province2 != ""
	replace prov_sim_score = 0 if province1 == "" | province2 == ""
drop p_province province1 province2


**Facility
gen facility_sim_score = log(m_facility/p_facility)/log(2) if facility_code1==facility_code2
	replace facility_sim_score = log((1-m_facility)/(1-p_facility))/log(2) if facility_code1!=facility_code2 & facility_code1 != "" & facility_code2 != ""
	replace facility_sim_score = 0 if facility_code1 == "" | facility_code2 == ""

drop facility_code1 facility_code2 p_facility



**DOB
// Generate variable to store the DOB_cases
gen DOB_case = .
gen DOB_sim_score = .

scalar define p_M = 1/12
scalar define p_D = 1/30

// Recode Jan 1 M/D as missing
replace dob_month1 = . if (dob_month1 == 1 & dob_day1 == 1)
replace dob_day1 = . if (dob_month1 == . & dob_day1 == 1)
replace dob_month2 = . if (dob_month2 == 1 & dob_day2 == 1)
replace dob_day2 = . if (dob_month2 == . & dob_day2 == 1)


//// Y ==
// DOB_case 1: Y ==, M ==, D ==
replace DOB_case = 1 if (dob_year1==dob_year2 & dob_year1 != .) & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1==dob_day2 & dob_day1 != .)
replace DOB_sim_score = log( (m_Y/p_YOB)*(m_M/p_M)*(m_D/p_D)) / log(2) if DOB_case == 1

// DOB_case 2: Y ==, M ==, D !=
replace DOB_case = 2 if (dob_year1==dob_year2 & dob_year1 != .) & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1 != dob_day2 & dob_day1 != . & dob_day2 != .)
replace DOB_sim_score = log( (m_Y/p_YOB)*(m_M/p_M)*((1-m_D)/(1-p_D))) / log(2) if DOB_case == 2

// DOB_case 3: Y ==, M ==, D Flag
replace DOB_case = 3 if (dob_year1==dob_year2 & dob_year1 != .) & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1 == . | dob_day2 == .)
replace DOB_sim_score = log( (m_Y/p_YOB)*(m_M/p_M)) / log(2)  if DOB_case == 3

// DOB_case 4: Y ==, M !=, D ==
replace DOB_case = 4 if (dob_year1==dob_year2 & dob_year1 != .) & (dob_month1 != dob_month2 & dob_month1 != . & dob_month2 != .) & (dob_day1==dob_day2 & dob_day1 != .)
replace DOB_sim_score = log( (m_Y/p_YOB)*((1-m_M)/(1-p_M))*(m_D/p_D)) / log(2) if DOB_case == 4

// DOB_case 5: Y ==, M !=, D !=
replace DOB_case = 5 if (dob_year1==dob_year2 & dob_year1 != .) & (dob_month1 != dob_month2 & dob_month1 != . & dob_month2 != .) & (dob_day1 != dob_day2 & dob_day1 != . & dob_day2 != .)
replace DOB_sim_score = log( (m_Y/p_YOB)*((1-m_M)/(1-p_M))*((1-m_D)/(1-p_D))) / log(2) if DOB_case == 5

// DOB_case 6: Y ==, M !=, D Flag
replace DOB_case = 6 if (dob_year1==dob_year2 & dob_year1 != .) & (dob_month1 != dob_month2 & dob_month1 != . & dob_month2 != .) & (dob_day1 == . | dob_day2 == .)
replace DOB_sim_score = log( (m_Y/p_YOB)*((1-m_M)/(1-p_M))) / log(2) if DOB_case == 6

// DOB_case 7: Y ==, M Flag, D ==
replace DOB_case = 7 if (dob_year1==dob_year2 & dob_year1 != .) & (dob_month1 == . | dob_month2 == .) & (dob_day1==dob_day2 & dob_day1 != .)
replace DOB_sim_score = log( (m_Y/p_YOB)*((m_D)/(p_D))) / log(2) if DOB_case == 7

// DOB_case 8: Y ==, M Flag, D !=
replace DOB_case = 8 if (dob_year1==dob_year2 & dob_year1 != .) & (dob_month1 == . | dob_month2 == .) & (dob_day1 != dob_day2 & dob_day1 != . & dob_day2 != .)
replace DOB_sim_score = log( (m_Y/p_YOB)*((1-m_D)/(1-p_D))) / log(2) if DOB_case == 8

// DOB_case 9: Y ==, M Flag, D Flag
replace DOB_case = 9 if (dob_year1==dob_year2 & dob_year1 != .) & (dob_month1 == . | dob_month2 == .) & (dob_day1 == . | dob_day2 == .)
replace DOB_sim_score = log((m_Y/p_YOB))/log(2) if DOB_case == 9

//// Y !=
// DOB_case 10: Y !=, M ==, D ==
replace DOB_case = 10 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1==dob_day2 & dob_day1 != .)
replace DOB_sim_score = log( ((1-m_Y)/(1-p_YOB))*(m_M/p_M)*(m_D/p_D)) / log(2) if DOB_case == 10

// DOB_case 11: Y !=, M ==, D !=
replace DOB_case = 11 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1 != dob_day2 & dob_day1 != . & dob_day2 != .)
replace DOB_sim_score = log(((1-m_Y)/(1-p_YOB))*(m_M/p_M)*((1-m_D)/(1-p_D))) / log(2) if DOB_case == 11

// DOB_case 12: Y !=, M ==, D Flag
replace DOB_case = 12 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1 == . | dob_day2 == .)
replace DOB_sim_score = log(((1-m_Y)/(1-p_YOB))*(m_M/p_M)) / log(2) if DOB_case == 12

// DOB_case 13: Y !=, M !=, D ==
replace DOB_case = 13 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (dob_month1 != dob_month2 & dob_month1 != . & dob_month2 != .) & (dob_day1==dob_day2 & dob_day1 != .)
replace DOB_sim_score = log(((1-m_Y)/(1-p_YOB))*((1-m_M)/(1-p_M))*(m_D/p_D)) / log(2) if DOB_case == 13

// DOB_case 14: Y !=, M !=, D !=
replace DOB_case = 14 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (dob_month1 != dob_month2 & dob_month1 != . & dob_month2 != .) & (dob_day1 != dob_day2 & dob_day1 != . & dob_day2 != .)
replace DOB_sim_score = log(((1-m_Y)/(1-p_YOB))*((1-m_M)/(1-p_M))*((1-m_D)/(1-p_D))) / log(2) if DOB_case == 14

// DOB_case 15: Y !=, M !=, D Flag
replace DOB_case = 15 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (dob_month1 != dob_month2 & dob_month1 != . & dob_month2 != .) & (dob_day1 == . | dob_day2 == .)
replace DOB_sim_score = log(((1-m_Y)/(1-p_YOB))*((1-m_M)/(1-p_M))) / log(2) if DOB_case == 15

// DOB_case 16: Y !=, M Flag, D ==
replace DOB_case = 16 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (dob_month1 == . | dob_month2 == .) & (dob_day1==dob_day2 & dob_day1 != .)
replace DOB_sim_score = log(((1-m_Y)/(1-p_YOB))*(m_D/p_D)) / log(2)  if DOB_case == 16

// DOB_case 17: Y !=, M Flag, D !=
replace DOB_case = 17 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (dob_month1 == . | dob_month2 == .) & (dob_day1 != dob_day2 & dob_day1 != . & dob_day2 != .)
replace DOB_sim_score = log(((1-m_Y)/(1-p_YOB))*((1-m_D)/(1-p_D))) / log(2)  if DOB_case == 17

// DOB_case 18: Y !=, M Flag, D Flag   // Note, many of these are Jan 1st M/D of birth, which we've replaced as missing. //
replace DOB_case = 18 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (dob_month1 == . | dob_month2 == .) & (dob_day1 == . | dob_day2 == .)
replace DOB_sim_score = log(((1-m_Y)/(1-p_YOB))) / log(2)  if DOB_case == 18

//// Y Flag (i.e. missing)
// DOB_case 19: Y Flag, M ==, D ==
replace DOB_case = 19 if (dob_year1 == . | dob_year2 == .) & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1==dob_day2 & dob_day1 != .)
replace DOB_sim_score = log((m_M/p_M)*(m_D/p_D)) / log(2)  if DOB_case == 19

// DOB_case 20: Y Flag, M ==, D !=
replace DOB_case = 20 if (dob_year1 == . | dob_year2 == .) & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1 != dob_day2 & dob_day1 != . & dob_day2 != .)
replace DOB_sim_score = log((m_M/p_M)*((1-m_D)/(1-p_D))) / log(2) if DOB_case == 20

// DOB_case 21: Y Flag, M ==, D Flag
replace DOB_case = 21 if (dob_year1 == . | dob_year2 == .) & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1 == . | dob_day2 == .)
replace DOB_sim_score = log((m_M/p_M)) / log(2) if DOB_case == 21

// DOB_case 22: Y Flag, M !=, D ==
replace DOB_case = 22 if (dob_year1 == . | dob_year2 == .) & (dob_month1 != dob_month2 & dob_month1 != . & dob_month2 != .) & (dob_day1==dob_day2 & dob_day1 != .)
replace DOB_sim_score = log(((1-m_M)/(1-p_M))*(m_D/p_D)) / log(2) if DOB_case == 22

// DOB_case 23: Y Flag, M !=, D !=
replace DOB_case = 23 if (dob_year1 == . | dob_year2 == .) & (dob_month1 != dob_month2 & dob_month1 != . & dob_month2 != .) & (dob_day1 != dob_day2 & dob_day1 != . & dob_day2 != .)
replace DOB_sim_score = log(((1-m_M)/(1-p_M))*((1-m_D)/(1-p_D))) / log(2) if DOB_case == 23

// DOB_case 24: Y Flag, M !=, D Flag
replace DOB_case = 24 if (dob_year1 == . | dob_year2 == .) & (dob_month1 != dob_month2 & dob_month1 != . & dob_month2 != .) & (dob_day1 == . | dob_day2 == .)
replace DOB_sim_score = log(((1-m_M)/(1-p_M))) / log(2) if DOB_case == 24

// DOB_case 25: Y Flag, M Flag, D ==
replace DOB_case = 25 if (dob_year1 == . | dob_year2 == .) & (dob_month1 == . | dob_month2 == .) & (dob_day1==dob_day2 & dob_day1 != .)
replace DOB_sim_score = log((m_D/p_D)) / log(2) if DOB_case == 25

// DOB_case 26: Y Flag, M Flag, D !=
replace DOB_case = 26 if (dob_year1 == . | dob_year2 == .) & (dob_month1 == . | dob_month2 == .) & (dob_day1 != dob_day2 & dob_day1 != . & dob_day2 != .)
replace DOB_sim_score = log(((1-m_D)/(1-p_D))) / log(2) if DOB_case == 26

// DOB_case 27: Y Flag, M Flag, D Flag
replace DOB_case = 27 if (dob_year1 == . | dob_year2 == .) & (dob_month1 == . | dob_month2 == .) & (dob_day1 == . | dob_day2 == .)
replace DOB_sim_score = 0 if DOB_case == 27


////// Corrections for special DOB_cases

//// Month/Day Inversions
// Give 50% credit to match, 50% to non-match

// DOB_case 28: Y ==, M/D inversion (correction to DOB_case 5, Y==, M!=, D!=)
replace DOB_case = 28 if DOB_case==5 & (dob_month1 == dob_day2) & (dob_month2 == dob_day1)
replace DOB_sim_score = .5*log( (m_Y/p_YOB)*(m_M/p_M)*(m_D/p_D)) / log(2) + (1-.5)*log( (m_Y/p_YOB)*((1-m_M)/(1-p_M))*((1-m_D)/(1-p_D))) / log(2) ///
	if DOB_case == 28

// DOB_case 29: Y !=, M/D inversion (correction to DOB_case 14, Y!=, M!=, D!=)
replace DOB_case = 29 if DOB_case==14 & (dob_month1 == dob_day2) & (dob_month2 == dob_day1)
replace DOB_sim_score = .5*log( ((1-m_Y)/(1-p_YOB))*(m_M/p_M)*(m_D/p_D)) / log(2) + (1-.5)*log( ((1-m_Y)/(1-p_YOB))*((1-m_M)/(1-p_M))*((1-m_D)/(1-p_D))) / log(2) ///
	if DOB_case == 29

// DOB_case 30: Y Flag, M/D inversion (correction to DOB_case 23, Y Flag, M!=, D!=)
replace DOB_case = 30 if DOB_case==32 & (dob_month1 == dob_day2) & (dob_month2 == dob_day1)
replace DOB_sim_score = .5*log((m_M/p_M)*(m_D/p_D)) / log(2) + (1-.5)*log(((1-m_M)/(1-p_M))*((1-m_D)/(1-p_D))) / log(2) ///
	if DOB_case == 30

//// Similar (but not equal) year, day
// DOB_case 31: Y !=, M ==, D ==, but Y-diff <10y (correction to DOB_case 10)
// Triangular kernel function, with penalty for non-equal
// When M== and D==, give credit for Y closeness, but with a 50% penalty for being non-equal 
// 50% penalty, 10% decrement per year diff
replace DOB_case = 31 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (abs(dob_year2 - dob_year1)<10) & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1==dob_day2 & dob_day1 != .)
	gen A = log( ((1-m_Y)/(1-p_YOB))*(m_M/p_M)*(m_D/p_D)) / log(2)
	gen B = log( ((m_Y)/(p_YOB))*(m_M/p_M)*(m_D/p_D)) / log(2)
replace DOB_sim_score = max(A, A + 0.5*(B-A)*(1 - 0.1*abs(dob_year2-dob_year1))) if DOB_case==31
drop A B

// DOB_case 32: Y !=, M Flag, D Flag, but Y-diff <10y (correction to DOB_case 18)
// NOTE: many of these were Jan 1 M/D of birth
// NOTE: some scores increase because they were previously considered non-matches; some decrease because they were previously considered matches (incorrectly)
// Triangular kernel function, with penalty for non-equal
// When M== and D==, give credit for Y closeness, but with a 50% penalty for being non-equal 
// 50% penalty, 10% decrement per year diff
replace DOB_case = 32 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (abs(dob_year2 - dob_year1)<10) & ((dob_month1== . & dob_day1 == .) | (dob_month2 == . & dob_day2 == .))
	gen A = log( ((1-m_Y)/(1-p_YOB))) / log(2)
	gen B = log( ((m_Y)/(p_YOB))) / log(2)
replace DOB_sim_score = max(A, A + 0.5*(B-A)*(1 - 0.1*abs(dob_year2-dob_year1))) if DOB_case==32
drop A B


// DOB_case 33: Y==, M==, D!= but D-diff <5 days (correction to DOB_case 2)
// 50% penalty, 20% decrement per day diff
replace DOB_case = 33 if (dob_year1 == dob_year2 & dob_year1 != .)  & (dob_month1==dob_month2 & dob_month1 != .) & (dob_day1!=dob_day2) & abs(dob_day2 - dob_day1) <5
	gen A = log( (m_Y/p_YOB)*(m_M/p_M)*((1-m_D)/(1-p_D))) / log(2)
	gen B = log( ((m_Y)/(p_YOB))*(m_M/p_M)*(m_D/p_D)) / log(2)
replace DOB_sim_score = max(A, A + 0.5*(B-A)*(1 - 0.2*abs(dob_day2-dob_day1))) if DOB_case==33
drop A B



////// AGE FLAG - DAY AND MONTH SHOULD BE CONSIDERED MISSING, YEAR SHOULD BE TREATED AS FUZZY
// AGE-FLAG indicates that M/D of birth date == M/D of test_date, which occurs in CDW when a patient gives an 
// age but no DOB; Note: AGE-FLAG over rules all matches on M / D

// DOB_case 34: Y ==, M?, D?, Age Flag (correction to DOB_cases 1-6)
// Just treat as M = missing, D = missing
replace DOB_case = 34 if (dob_year1==dob_year2 & dob_year1 != .) & (age_flag1 == 1 | age_flag2 == 1)
replace DOB_sim_score = log((m_Y/p_YOB))/log(2) if DOB_case == 34

// DOB_case 35: Y !=, M?, D?, Age Flag, but Y-diff <10y (correction to DOB_case 18)
// When Age-Flag, we might expect year to be fuzzy. Here we reduce the penalty for not matching to 25%, giving 75% credit to cases where the year is different
// NOTE: some scores increase because they were previously considered non-matches; some decrease because they were previously considered matches (incorrectly)
// 25% penalty, 10% decrement per year diff
replace DOB_case = 35 if (dob_year1 != dob_year2 & dob_year1 != . & dob_year2 != .)  & (abs(dob_year2 - dob_year1)<10) & (age_flag1 == 1 | age_flag2 == 1)
	gen A = log( ((1-m_Y)/(1-p_YOB))) / log(2)
	gen B = log( ((m_Y)/(p_YOB))) / log(2)
replace DOB_sim_score = max(A, A + 0.75*(B-A)*(1 - 0.1*abs(dob_year2-dob_year1))) if DOB_case==35
drop A B

tab DOB_case

/*   DOB_case |      Freq.     Percent        Cum.
------------+-----------------------------------
          1 |131,294,803       30.33       30.33
          2 | 42,216,181        9.75       40.08
          4 | 23,028,512        5.32       45.40
          5 | 56,434,369       13.04       58.43
          9 | 25,343,303        5.85       64.29
         10 |  1,059,599        0.24       64.53
         11 |  2,441,396        0.56       65.10
         13 |  2,137,722        0.49       65.59
         14 | 20,960,882        4.84       70.43
         18 |  4,337,560        1.00       71.43
         28 |  3,384,231        0.78       72.22
         29 |    107,940        0.02       72.24
         31 |  5,878,930        1.36       73.60
         32 |    788,364        0.18       73.78
         33 | 16,354,560        3.78       77.56
         34 | 81,726,937       18.88       96.44
         35 | 15,427,323        3.56      100.00
------------+-----------------------------------
      Total |432,922,612      100.00*/
      
drop DOB_case




***********************************************************************************************
* C.2 Construct total similarity score as weighted average of component similarity scores
***********************************************************************************************

*tot_sim_score is the vector of ones
*r_jimmy is optimized to max PPV across all values of Sen, using bootstrap aggregation (BAGging) to avoid over-fitting training data 

* NOTE: see scored_training_data.do for description of methods used to get weights

gen tot_sim_score = DOB_sim_score + first_sim_score + last_sim_score + gender_sim_score + prov_sim_score + facility_sim_score
*gen r_jimmy_old = 1.433*first_sim_score + 1.194*last_sim_score + 1.131*gender_sim_score + 1.264*DOB_sim_score + 0.821*prov_sim_score + 0.597*facility_sim_score
gen r_jimmy = 1.0547*first_sim_score + 1.0969*last_sim_score + 1.1856*gender_sim_score + 1.2794*DOB_sim_score + 0.8955*prov_sim_score + 0.7152*facility_sim_score

** Translate score into predicted probability of a match using training data (see "scored_training_data.do")
*gen p = exp(-11.13 + 0.366*r_jimmy) / (1 + exp(-11.13 + 0.366*r_jimmy))

compress

count
*432,922,612

save "HIV_all_edges_tot_sim.dta", replace







*******************************************************************
* C.3 Identify and remove duplicates
*******************************************************************

use "HIV_all_edges_tot_sim.dta", clear

* summary stats, including duplicates
tab matchtype


/* 
    matchtype |      Freq.     Percent        Cum.
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

tabstat r_jimmy, by(matchtype) stat(n mean p25 p50 p75)
/*
    matchtype |         N      mean       p25       p50       p75
--------------+--------------------------------------------------
        Fuzzy |  3.10e+08  21.37542  13.73311  19.21334  26.04547
    Inversion |    471299  48.68397  43.07122  48.10329  53.99965
   Mult first |   2100920  47.48475  42.03206  46.79132   52.4221
    Mult last |     55924  48.12758   42.2393  47.41458  53.71954
     Nickname |   3959943  44.25827  39.33115  44.41992   49.2976
Fuzzy_sex_fac |  3.24e+07  47.10983  42.41461  46.86071  51.35619
    Fuzzy_DOB |  6.26e+07   28.0078   14.1103  28.24055  40.16695
  Fuzzy_first |  1.06e+07  46.10568   39.5642  46.14291  52.78141
   Fuzzy_last |   6840290  45.61346  37.09605  47.69495  54.04376
Fuzzy_firstNA |     95306  39.34991   35.4568  39.96412  43.85266
Fuzzy_DOB1800 |   3349716  24.41828  19.53535  23.11914  28.26265
--------------+--------------------------------------------------
        Total |  4.33e+08  25.64534  14.73766  21.48799  36.57946
-----------------------------------------------------------------
*/

* identify and remove duplicates (for some duplicates, there may be differences in first_sim and/or last_sim if there was both a nickname and a JW score. Keep the higher valued edge.)
* order EM_ID_plus1 EM_ID_plus2 pairs so that EM_ID_plus1 < EM_ID_plus2
gen double temp = EM_ID_plus2
replace EM_ID_plus2 = EM_ID_plus1 if EM_ID_plus1 > EM_ID_plus2
replace EM_ID_plus1 = temp if EM_ID_plus1 > temp
drop temp
* sort and keep highest scored edge (using r_jimmy)
gsort EM_ID_plus1 EM_ID_plus2 -r_jimmy
by EM_ID_plus1 EM_ID_plus2: keep if _n==1
* count all unique edges
count
*305,576,121
save "HIV_all_edges_tot_sim.dta", replace




**************************************
* C.4 Quick summary stats of edges
**************************************

* Quick summary stats of edges; duplicates excluded

tab matchtype


/*    matchtype |      Freq.     Percent        Cum.
--------------+-----------------------------------
        Fuzzy |248,936,431       81.46       81.46
    Inversion |    350,374        0.11       81.58
   Mult first |  2,028,678        0.66       82.24
    Mult last |     50,142        0.02       82.26
     Nickname |  1,973,170        0.65       82.91
Fuzzy_sex_fac | 16,206,838        5.30       88.21
    Fuzzy_DOB | 27,706,098        9.07       97.28
  Fuzzy_first |  2,822,228        0.92       98.20
   Fuzzy_last |  2,057,140        0.67       98.87
Fuzzy_firstNA |     95,306        0.03       98.90
Fuzzy_DOB1800 |  3,349,716        1.10      100.00
--------------+-----------------------------------
        Total |305,576,121      100.00
*/

tabstat r_jimmy, by(matchtype) stat(n mean p25 p50 p75)
/*     matchtype |         N      mean       p25       p50       p75
--------------+--------------------------------------------------
        Fuzzy |  2.49e+08  20.98757  13.83295  19.19637  25.58463
    Inversion |    350374  48.68246  43.08085  48.10092  53.99532
   Mult first |   2028678  47.11602  41.85851  46.52211   51.8796
    Mult last |     50142  47.59323  41.92123  46.92505   53.0013
     Nickname |   1973170  43.31832  38.35236  43.68353  48.75554
Fuzzy_sex_fac |  1.62e+07  47.10983  42.41461  46.86071  51.35619
    Fuzzy_DOB |  2.77e+07  26.05972  12.92156  24.61351  37.84451
  Fuzzy_first |   2822228  43.37237  36.84603  42.92631  49.50465
   Fuzzy_last |   2057140  42.39894   32.6714  44.04477  52.15284
Fuzzy_firstNA |     95306  39.34991   35.4568  39.96412  43.85266
Fuzzy_DOB1800 |   3349716  24.41828  19.53535  23.11914  28.26265
--------------+--------------------------------------------------
        Total |  3.06e+08   23.5809  14.37151  20.36167  30.24183
-----------------------------------------------------------------*/



********************************************************
* C.5 Create edgelist to export to PART D, fuzzyGraph.R
********************************************************

* Drop excess variables
*use "HIV_all_edges_tot_sim.dta", clear
keep EM_ID_plus1 EM_ID_plus2 r_jimmy
order EM_ID_plus1 EM_ID_plus2 r_jimmy
save "HIV_all_edges_R.dta", replace







********************************************************
* C.5 Create alternate edgelist and vertices list,
*	restricting to CD4 and VL
********************************************************

* Keep only those episode_no's associated with CD4 or VL
use episode_no using "../raw_data/pre2015tests/tests_demo_may2015_2.dta", clear
append using  "../raw_data/post2015tests/Bill_viral_loads.dta", keep(episode_no)
append using  "../raw_data/post2015tests/Bill_cd4.dta", keep(episode_no)
bys episode_no: keep if _n==1
	
merge 1:1 episode_no using "HIV_EMIDplus_episode.dta"
   
    /*     Result                           # of obs.
    -----------------------------------------
    not matched                    71,128,901
        from master                   564,645  (_merge==1)
        from using                 70,564,256  (_merge==2)

    matched                        46,005,030  (_merge==3)
    ----------------------------------------- */
    
keep if _merge==3
drop _merge
save "HIV_EMIDplus_episode_CD4VL.dta", replace

* Create extract of HIV pre-processed data, limited to CD4/VL results
merge 1:1 episode_no using "HIV_demog_all_names.dta"
keep if _merge==3
drop _merge
compress
save "HIV_demog_all_names_CD4VL.dta", replace

* Collapse back to EM_ID_plus 
keep EM_ID_plus
bys EM_ID_plus: keep if _n==1
save "HIV_EMIDplus_CD4VL.dta", replace   /*CD4/VL vertices*/

* Double-merge with edges and keep edges where both sides are CD4/VL
use "HIV_EMIDplus_CD4VL.dta", clear
rename EM_ID_plus EM_ID_plus1
merge 1:m EM_ID_plus1 using "HIV_all_edges_R.dta"
keep if _merge==3
drop _merge
rename EM_ID_plus2 EM_ID_plus
merge m:1 EM_ID_plus using "HIV_EMIDplus_CD4VL.dta"
keep if _merge==3
drop _merge
rename EM_ID_plus EM_ID_plus2
keep EM_ID_plus1 EM_ID_plus2 r_jimmy
order EM_ID_plus1 EM_ID_plus2 r_jimmy
count
*104,400,240
save "HIV_CD4VL_edges_R.dta", replace


timer off 3
timer list 3


log close


exit
