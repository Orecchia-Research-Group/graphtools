

****************************************************************************************************************
****************************************************************************************************************
**** 
****   PART C: LINK BU_UNIQ_ID TO EPISODE_NO AND EXPORT. VALIDATE AGAINST MANUALLY-MATCHED DATASET.
**** 
****************************************************************************************************************
****************************************************************************************************************

* the script will call stata using:
* do matching_algorithm_19feb2016_validation_loop.do diameter t_low t_high weight
* Note: diameter = 2, 3
* t_low = 18,20,22,23,25
* t_high = 25,27,28,30,40
* weight = 3 (ones), 4(r), 5(r1)



cd "/restricted/projectnb/salabs/RecordLinkage/HIV_cohort/alg_output_files"

local d = `1'
local w = `2'
local m = `3'


local d100 = `d'*100

/*

*** C.1 Import fuzzyGraph output

	/*  ** Older version of wd_fuzzyGraph.R exported data as labeled values; but now the values are the labels.
forvalues i = 1/`m' {
	use "wd_G_out_final_dta/G_out_final_A_wd0`d'_`w'_`i'.dta", clear
	decode ref_id, gen(temp)
	destring temp, gen(EM_ID_plus)
	rename new_id BU_uniq_ID
	drop temp ref_id
	save "wd_G_out_final_dta/G_out_final_A_wd0`d'_`w'_`i'.dta", replace
	}
	*/


use "wd_G_out_final_dta/G_out_final_A_wd0`d'_`w'_1.dta", clear
gen old = 1

forvalues i = 2/`m' {
	append using "wd_G_out_final_dta/G_out_final_A_wd0`d'_`w'_`i'.dta"
	su new_id if old==1
	replace new_id = new_id + r(max) if old == .
	replace old = 1
	}
drop old

rename new_id BU_uniq_ID
rename ref_id EM_ID_plus


*** Expand to Episode_No
merge 1:m EM_ID_plus using "HIV_EMIDplus_episode.dta"
*merge 1:m EM_ID_plus using "/restricted/projectnb/salabs/RecordLinkage/Feb2016/Full275M/Full275_EMIDplus_episode_CD4VLworkup.dta"
*merge 1:m EM_ID_plus using "/restricted/projectnb/salabs/RecordLinkage/Feb2016/Full275M/Full275_EMIDplus_episode.dta"
keep if _merge ==3
drop _merge
merge 1:1 episode_no using "HIV_demog_all_names.dta"
*merge 1:1 episode_no using "/restricted/projectnb/salabs/RecordLinkage/Feb2016/Full275M/demog_oct2015r_2_names_CD4VLworkup.dta"
*merge 1:1 episode_no using "/restricted/projectnb/salabs/RecordLinkage/Feb2016/Full275M/demog_oct2015r_2_names.dta"
keep if _merge ==3
drop _merge first last
compress

merge 1:1 episode_no using "HIV_EMIDplus_episode_CD4VL.dta" 
gen m3 = _merge==3
drop _merge
bys BU_uniq_ID: egen anyCD4VL = max(m3)

count
*116,569,286
keep if anyCD4VL == 1
count
*72,683,441
drop anyCD4VL m3

save "BUuniqID_wd_`d100'_`w'.dta", replace

*/

*local w = 3
*local d100 = "50m90"
*local d100 = 60m85
*local d100 = 60m80


** QUICK DESCRIPTIVE ANALYSIS
scalar drop _all
use "BUuniqID_wd_`d100'_`w'.dta", clear
count
scalar define est_N_episodes_wd = r(N)

count if episode_no ==""
scalar define est_null_episodes_wd = r(N)

sort BU_uniq_ID test_year test_month test_day
by BU_uniq_ID: gen first_test_date = _n==1
count if first_test_date==1
scalar define est_N_BU_uniq_ID_wd = r(N)

* BU_uniq_ID's with >X EMID_plus
bys EM_ID_plus: gen first_EM_ID_plus = _n==1
bys BU_uniq_ID: egen EMIDplus_per_BU_uniq_ID = total(first_EM_ID_plus)
count if EMIDplus_per_BU_uniq_ID >=10 & first_test_date==1
scalar define est_EMIDplus_10_wd = r(N)
count if EMIDplus_per_BU_uniq_ID >=25 & first_test_date==1
scalar define est_EMIDplus_25_wd = r(N)

* BU_uniq_ID's with >X episodes
by BU_uniq_ID: gen N_per_BU_uniq_ID = _N
count if N_per_BU_uniq_ID >=25 & first_test_date==1
scalar define est_N25_wd = r(N)
count if N_per_BU_uniq_ID >=50 & first_test_date==1
scalar define est_N50_wd = r(N)
count if N_per_BU_uniq_ID >=100 & first_test_date==1
scalar define est_N100_wd = r(N)
count if N_per_BU_uniq_ID >=1000 & first_test_date==1
scalar define est_N1000_wd = r(N)
count if N_per_BU_uniq_ID >=10000 & first_test_date==1
scalar define est_N10000_wd = r(N)
count if N_per_BU_uniq_ID >=100000 & first_test_date==1
scalar define est_N100000_wd = r(N)

* prov_jump, facility_jump, mmi
preserve
	drop if facility==""
	gen jump_prov=0
		by BU_uniq_ID: replace jump_prov = 1 if province != province[_n-1] & _n>1
		by BU_uniq_ID: egen n_jump_prov = sum(jump_prov)
		gen n_jump_prov_5 = n_jump_prov >=5
	gen jump_fac=0
		by BU_uniq_ID: replace jump_fac = 1 if facility != facility[_n-1] & _n>1
		by BU_uniq_ID: egen n_jump_fac = sum(jump_fac)
		gen n_jump_fac_5 = n_jump_fac >=5
	gen female = 0 if gender=="M"
		replace female = 1 if gender =="F"
		by BU_uniq_ID: egen mean_fem = mean(female)
		gen mixed_sex = 0
		replace mixed_sex = 1 if mean_fem > 0 & mean_fem <1
	gen jump_fac_sex=0
		by BU_uniq_ID: replace jump_fac_sex = 1 if facility != facility[_n-1] & gender != gender[_n-1] & _n>1
		by BU_uniq_ID: egen n_jump_fac_sex = sum(jump_fac_sex)
		gen n_jump_fac_sex1 = n_jump_fac_sex >=1

	su n_jump_prov_5 if first_test_date == 1
	scalar define est_JumpProv5_wd = r(mean)

	su n_jump_fac_5 if first_test_date == 1
	scalar define est_JumpFac5_wd = r(mean)

	su mixed_sex if first_test_date == 1
	scalar define est_MixedSex_wd = r(mean)

	su n_jump_fac_sex1 if first_test_date == 1
	scalar define est_JumpFacSex_wd = r(mean)
restore






*** Calculate Sen and PPV in manual match training data 

* identify ref_episodes in training set
*use "/restricted/projectnb/salabs/RecordLinkage/man.match.data/MMredux_2_UNBLINDED_final.dta", clear
use "../alg_input_files/Training_Data.dta", clear
keep ref_episode
bys ref_episode: keep if _n==1
*count
gen episode_no = ref_episode 

* merge ref_episodes with the big BU_ID file
merge 1:1 episode_no using "BUuniqID_wd_`d100'_`w'.dta", keepusing(BU_uniq_ID)
gen to_include = _merge==3
drop _merge

* drop obeservations not in the original dataset used for manual matching
merge m:1 episode_no using "/restricted/projectnb/salabs/RecordLinkage/Dec2014/Bill_demo_clean_orig2014data_episodes.dta"
keep if _merge==3
drop _merge

* and keep only the CD4/VL
merge m:1 episode_no using "HIV_EMIDplus_episode_CD4VL.dta", keepusing(episode_no)
keep if _merge==3
drop _merge

* identify other episodes linked to ref_episodes through BU_ID (or EM_ID or UP_ID)
bys BU_uniq_ID: egen compmatch = max(to_include)
* drop all other episodes
keep if compmatch == 1 

* apply ref_episode_no to all episodes
bys BU_uniq_ID: egen BU_ref_episode_no = mode(ref_episode_no)
replace ref_episode_no = BU_ref_episode_no if ref_episode_no == ""
drop BU_ref 

* merge ref_episode with training dataset as a normal episode
*merge 1:m ref_episode_no episode_no using "/restricted/projectnb/salabs/RecordLinkage/man.match.data/MMredux_2_UNBLINDED_final.dta", keepusing(manmatch r_manmatch is_match __merge)
merge 1:m ref_episode_no episode_no using "../alg_input_files/Training_Data.dta", keepusing(manmatch r_manmatch)
sort ref_episode_no episode_no

* recode matching concepts
*computer match
replace compmatch = 0 if compmatch==.
*original manual match, excl. MMredux1a and MMredux2
*replace is_match = . if __merge==1
*Man Match, including MMredux1a
replace manmatch = 0 if manmatch==.
*replace manmatch = . if __merge==1
*Man Match, including MMredux1a and MMredux2
replace r_manmatch = 0 if r_manmatch ==.
*replace r_manmatch = . if __merge==1

* drop episodes that are neither a manual match nor a compmatch
drop if r_manmatch==0 & manmatch==0 & compmatch ==0

* drop ref episodes
drop if ref_episode==episode_no

*** Calculate Sensitivity and PPV

*su is_match if compmatch ==1 /*PPV*/
*scalar define est_PPV1_wd = r(mean)
*su compmatch if is_match ==1 /*sensitivity*/
*scalar define est_Sen1_wd = r(mean)

su manmatch if compmatch ==1 /*PPV*/
scalar define est_PPV2_wd = r(mean)
su compmatch if manmatch ==1 /*sensitivity*/
scalar define est_Sen2_wd = r(mean)

su r_manmatch if compmatch ==1 /*PPV*/
scalar define est_PPV3_wd = r(mean)
su compmatch if r_manmatch ==1 /*sensitivity*/
scalar define est_Sen3_wd = r(mean)



/* Manual inspection of results
drop _merge
merge 1:1 episode_no ref_episode_no using "../alg_input_files/Training_Data_names.dta"
keep if _merge == 3

sort ref_episode_no compmatch r_manmatch

browse

save "man_inspection_wd_85.dta", replace

rename compmatch compmatch85

drop _merge
merge 1:1 episode_no ref_episode_no using "man_inspection_wd_60.dta", keepusing(compmatch)

rename compmatch compmatch60

move compmatch60 compmatch85
sort ref_episode compmatch60 compmatch85

*/




/* 
**Bootstrap? This is helpful for getting the variance correct (though the asymptotic CI is fine), but has no bearing on the point estiamte.
set seed 6819
bootstrap PPV3 = r(mean), reps(1000) cluster(ref_episode) saving(bs_ppv3.dta, replace): su r_manmatch if compmatch ==1
estat bootstrap, v all
preserve
use bs_ppv3.dta, clear
su PPV3
hist PPV3
restore

set seed 6819
bootstrap Sen3 = r(mean), reps(1000) cluster(ref_episode) saving(bs_sen3.dta, replace): su compmatch if r_manmatch ==1 
estat bootstrap, v all
preserve
use bs_sen3.dta, clear
su Sen3
hist Sen3
restore
*/



/* Run this once before starting start_tasks.sh
clear
set obs 1000
gen diameter = .
gen weight = .
gen t_low = .
gen t_high = .
gen long N_episodes = .
gen long N_null_episodes = .
gen long N_BU_uniq_ID = .
gen Sen1 = .
label var Sen1 "Sen vs. original man-match data"
gen PPV1 = .
label var PPV1 "PPV vs. original man-match data"
gen Sen2 = .
label var Sen2 "Sen vs. man-match data, incl. MMredux1a"
gen PPV2 = .
label var PPV2 "PPV vs. man-match data, incl. MMredux1a"
gen Sen3 = .
label var Sen3 "Sen vs. man-match data, incl. MMredux1a & MMredux2"
gen PPV3 = .
label var PPV3 "PPV vs. man-match data, incl. MMredux1a & MMredux2"
gen N_EMIDplus_10 = .
gen N_EMIDplus_25 = .
gen N25 = .
gen N50 = .
gen N100 = .
gen N1000 = .
gen N10000 = .
gen N100000 = .
gen JumpProv5 = .
gen JumpFac5 = .
gen MixedSex = .
gen JumpFacSex = .
save "wd_ART_Validation_Loop_Out.dta", replace
*/


use "wd_ART_Validation_Loop_Out.dta", clear
count if diameter < .
local nextobs = r(N)+1
dis `nextobs'
replace diameter = `d' in `nextobs'
replace weight = `w' in `nextobs'
*replace t_low = `t_low' in `nextobs'
*replace t_high = `t_high' in `nextobs'
replace N_episodes = est_N_episodes_wd in `nextobs'
replace N_null_episodes = est_null_episodes_wd in `nextobs'
replace N_BU_uniq_ID = est_N_BU_uniq_ID_wd in `nextobs'
*replace PPV1 = est_PPV1_wd in `nextobs'
*replace Sen1 = est_Sen1_wd in `nextobs'
replace PPV2 = est_PPV2_wd in `nextobs'
replace Sen2 = est_Sen2_wd in `nextobs'
replace PPV3 = est_PPV3_wd in `nextobs'
replace Sen3 = est_Sen3_wd in `nextobs'
replace N_EMIDplus_10 = est_EMIDplus_10_wd in `nextobs'
replace N_EMIDplus_25 = est_EMIDplus_25_wd in `nextobs'
replace N25 = est_N25_wd in `nextobs'
replace N50 = est_N50_wd in `nextobs'
replace N100 = est_N100_wd in `nextobs'
replace N1000 = est_N1000_wd in `nextobs'
replace N10000 = est_N10000_wd in `nextobs'
replace N100000 = est_N100000_wd in `nextobs'
replace JumpProv5 = est_JumpProv5_wd in `nextobs'
replace JumpFac5 = est_JumpFac5_wd in `nextobs'
replace MixedSex = est_MixedSex_wd in `nextobs'
replace JumpFacSex = est_JumpFacSex_wd in `nextobs'
save "wd_ART_Validation_Loop_Out.dta", replace

scalar drop _all

*rm "BUuniqID_wd_`d100'_`w'.dta"




/* ANALYSIS (Comment out if we run again)

cd "/restricted/projectnb/salabs/RecordLinkage/Feb2016/optim_fuzzyGraph/"

use "Validation_Loop_Out.dta", clear

sort diameter weight t_low t_high

bys diameter weight t_low t_high: gen n =_n

gen spec = _n

count if diameter <.
drop if n>1 & spec < r(N) 

*F-stat indifference curves
gen F_Sen = .8 + .0002*_n
gen F95_PPV = .95/(2-.95/F_Sen)
gen F90_PPV = .90/(2-.90/F_Sen)
gen F85_PPV = .85/(2-.85/F_Sen)
gen F80_PPV = .8/(2-.8/F_Sen)
replace F95_PPV = 1 if F95_PPV >1
replace F90_PPV = 1 if F90_PPV >1
replace F85_PPV = 1 if F85_PPV >1
replace F80_PPV = 1 if F80_PPV >1

*f = PPV*2 / (1+PPV/Sen)

twoway 	(area F95_PPV F_Sen if F95_PPV >.8, col(gs15) lw(vvthin)) ///
	(area F90_PPV F_Sen if F90_PPV >.8, col(gs14) lw(vvthin)) ///
	(area F85_PPV F_Sen if F85_PPV >.78, col(gs13) lw(vvthin)) ///
	(scatter PPV3 Sen3 if PPV3>.8, mlabel(spec) mlabpos(0) msymbol(none) mlabcol(navy)) ///
	(scatter PPV3 Sen3 if spec==65, mlabpos(0) msymbol(Oh) mcol(sandb) msize(8)) ///
	, xlabel(.8(.05)1, labsize(med)) ylabel(.8(.05)1, nogrid labsize(med)) ///
	ytitle(Positive Predictive Value, size(medlarge)) xtitle(Sensitivity, size(medlarge)) ///
	legend(off) graphregion(col(white))

/*twoway (scatter PPV1 Sen1, mlabel(spec) mlabpos(0) msymbol(none)) ///
	(line F95_PPV F_Sen if F95_PPV >.8 & F95_PPV<1, lw(vvthin) lc(gs10)) ///
	(line F90_PPV F_Sen if F90_PPV >.8 & F90_PPV<1, lw(vvthin) lc(gs10)) ///
	(line F85_PPV F_Sen if F85_PPV >.8 & F85_PPV<1, lw(vvthin) lc(gs10)) ///
	, xlabel(.8(.05)1) ylabel(.8(.05)1, nogrid) legend(off)*/

twoway 	(area F95_PPV F_Sen if F95_PPV >.8, col(gs15) lw(vvthin)) ///
	(area F90_PPV F_Sen if F90_PPV >.8, col(gs14) lw(vvthin)) ///
	(area F85_PPV F_Sen if F85_PPV >.8, col(gs13) lw(vvthin)) ///
	(scatter PPV2 Sen2 if PPV2>.8, mlabel(spec) mlabpos(0) msymbol(none) mlabcol(navy)) ///
	, xlabel(.8(.05)1, labsize(med)) ylabel(.8(.05)1, nogrid labsize(med)) ///
	ytitle(Positive Predictive Value, size(medlarge)) xtitle(Sensitivity, size(medlarge)) ///
	legend(off) graphregion(col(white))
	
twoway 	(area F95_PPV F_Sen if F95_PPV >.8, col(gs15) lw(vvthin)) ///
	(area F90_PPV F_Sen if F90_PPV >.8, col(gs14) lw(vvthin)) ///
	(area F85_PPV F_Sen if F85_PPV >.8, col(gs13) lw(vvthin)) ///
	(scatter PPV3 Sen3 if PPV3>.8, mlabel(spec) mlabpos(0) msymbol(none) mlabcol(navy)) ///
	, xlabel(.8(.05)1, labsize(med)) ylabel(.8(.05)1, nogrid labsize(med)) ///
	ytitle(Positive Predictive Value, size(medlarge)) xtitle(Sensitivity, size(medlarge)) ///
	legend(off) graphregion(col(white))
	
	
gen N_BU_M = N_BU/1000000	
twoway (scatter Sen3 N_BU_M, mlabel(spec) mlabpos(0) msymbol(none)) ///
	(scatter Sen3 N_BU_M if spec==65, mlabpos(0) msymbol(Oh) mcol(sandb) msize(8)) ///
	, xlabel(, labsize(med)) ylabel(.8(.05)1, nogrid labsize(med)) ///
	ytitle(Sensitivity, size(medlarge)) xtitle(Millions of Unique Patients , size(medlarge)) ///
	legend(off) graphregion(col(white))

twoway (scatter PPV3 N_BU_M if PPV3 >.8, mlabel(spec) mlabpos(0) msymbol(none)) ///
	(scatter PPV3 N_BU_M if spec==65, mlabpos(0) msymbol(Oh) mcol(sandb) msize(8)) ///
	, xlabel(, labsize(med)) ylabel(.8(.05)1, nogrid labsize(med)) ///
	ytitle(Positive Predictive Value, size(medlarge)) xtitle(Millions of Unique Patients , size(medlarge)) ///
	legend(off) graphregion(col(white))

twoway (scatter PPV3 JumpProv5 if PPV3 >.8, mlabel(spec) mlabpos(0) msymbol(none)), ///
	ylabel(.75(.05)1) 
twoway (scatter PPV3 JumpFac5 if PPV3 >.8, mlabel(spec) mlabpos(0) msymbol(none)), ///
	ylabel(.75(.05)1) 	
twoway (scatter PPV3 MixedSex if PPV3 >.8, mlabel(spec) mlabpos(0) msymbol(none)), ///
	ylabel(.75(.05)1) 	
twoway (scatter PPV3 JumpFacSex if PPV3 >.8, mlabel(spec) mlabpos(0) msymbol(none)), ///
	ylabel(.75(.05)1) 	

* Can we predict Sen and PPV?
corr PPV3 JumpProv5 JumpFac5 MixedSex JumpFacSex N_BU
cap drop PPV3_hat f1
pca JumpProv5 JumpFac5 MixedSex JumpFacSex 
predict f1
reg PPV3 f1
predict PPV3_hat

cap drop N_BU_sq Sen3_hat
gen N_BU_sq = (N_BU - 8000000)^2
reg Sen3 N_BU_uniq N_BU_sq
predict Sen3_hat
 
twoway 	(area F95_PPV F_Sen if F95_PPV >.8, col(gs15) lw(vvthin)) ///
	(area F90_PPV F_Sen if F90_PPV >.8, col(gs14) lw(vvthin)) ///
	(area F85_PPV F_Sen if F85_PPV >.8, col(gs13) lw(vvthin)) ///
	(scatter PPV3_hat Sen3_hat) (scatter PPV3 Sen3)	///
	, xlabel(.8(.05)1, labsize(med)) ylabel(.8(.05)1, nogrid labsize(med)) ///
	ytitle(Positive Predictive Value, size(medlarge)) xtitle(Sensitivity, size(medlarge)) ///
	legend(off) graphregion(col(white))
	
*/
