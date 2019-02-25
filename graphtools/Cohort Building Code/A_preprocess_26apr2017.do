**********************************************************************************************
*
*	<<< BU-HE2RO-NHLS RECORD LINKAGE ALGORITHM >>>
*
*	"A South African National HIV Cohort: Graph-Based Record Linkage of Routine Laboratory Data"
*
*	Collaboration between Boston University (BU), the Health Economics and Epidemiology Research Office (HE2RO).
*	and the National Health Laboratory Service (NHLS) of South Africa
*
*	Contact: Jacob Bor, jbor@bu.edu
*	Written by: Jacob Bor, Katia Oleinik, Bill MacLeod, Jimmy Potter
*	In consultation with: Sergio Carmona, Sue Candy, Wendy Stevens, Jaco Grobler, Matt Fox, Mhairi Maskew,
*	and additional HE2RO/BU/NHLS colleagues
*
*	Version: 26 April 2017, Previous versions: 7 Sept 2016, 19 Feb 2016
*
*	STATA-MP v12.1, SCC4.BU.EDU
*
**********************************************************************************************



************* Version notes *****************************************************************
*
* Algorithm is now (as of 26-apr-2017) adapted to be implemented on data updated through December 2016,
* including the following lab results: CD4, VL, HIV_PCR_elisa, ALT, Haemoglobin, Creatinine, and CrAG.
* All tests except CrAG available since 2004. CrAG available since 2015.
* The lab tests to be linked are all HIV-relevant laboratory tests, many may have been conducted for
* patients without HIV. Prior to linking the data, we cannot exclude these. However, a significant
* fraction of the total results linked will be omitted from the resulting National HIV Cohort
*
**********************************************************************************************




/////////////////////////////////////
/////////////////////////////////////
/////				/////
/////	A. PRE-PROCESSING	/////
/////				/////
/////////////////////////////////////
/////////////////////////////////////


************* PROCEDURES *****************************************************************************
* 
* This code implements PART A of our record linkage algorithm:
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


log using "A_preprocess_26apr2017.log", append

timer on 1
set more off


************************************************
* A.1.     Prep datasets
************************************************

 
** A.1.1 Read in full 275M through mid-2015

	* FULL DEMOG DATA - ALL 275 MILLION
	* sed -e  "s/\"//g" demog_oct2015_reduced.csv > demog_oct2015r_2.csv
	* sed -e "s/[\\]|//g" demog_oct2015r_2.csv > demog_oct2015r_3.csv
	* mv ../../Nov2015data/demog_oct2015r_3.csv ../raw_data/pre2015demog
	* Now, in Stata...
	*insheet using "../../Nov2015data/demog_oct2015r_2.csv", delimiter("|") clear names
	
	insheet using "../raw_data/pre2015demog/demog_oct2015r_3.csv", delimiter("|") clear names
	
	compress
		
	keep  episode_no  patient_name_lis patient_surname patient_firstname date_of_birth gender date_specimen unique_patient_id
	replace patient_name_lis = itrim(patient_name_lis)
	replace patient_surname = itrim(patient_surname)
	replace patient_firstname = itrim(patient_firstname)
	keep if episode_no !="" /*209 dropped*/
	gen len = length(episode_no) 
	keep if len==10 /*221 dropped */
	drop len
	*DOB
	gen dob_year = substr(date_of_birth,1,4)
	destring dob_year, replace force
	gen dob_month = substr(date_of_birth,6,2)
	destring dob_month, replace force
	gen dob_day = substr(date_of_birth,9,2)
	destring dob_day, replace force
	*Test Date
	gen test_year = substr(date_specimen,1,4)
	destring test_year, replace force
	gen test_month = substr(date_specimen,6,2)
	destring test_month, replace force
	gen test_day = substr(date_specimen,9,2)
	destring test_day, replace force
	*Clean up
	drop if date_specimen == "" | date_of_birth == ""
	drop date_specimen date_of_birth
	gen age_flag = dob_month==test_month & dob_day==test_day
	*bys episode_no: keep if _n==1 /*already unique*/
	destring unique_patient_id, replace
	compress
	count /*275,622,859*/
	save "../raw_data/pre2015demog/demog_oct2015.dta", replace
	
	* Merge in pre-2015 facility data: this will not be necessary eventually
	insheet using "../raw_data/pre2015tests/FACILITY_EPISODE.csv", delimiter("|") clear
	rename v1 episode_no
	rename v2 facility_code
	drop in 1
	bys episode_no: keep if _n==1
	merge 1:1 episode_no using "../raw_data/pre2015demog/demog_oct2015.dta"
	drop if _merge ==1 /* dropped 3870*/
	drop _merge
	save "../raw_data/pre2015demog/demog_oct2015.dta", replace 

	*get province
	gen location_code = facility_code
	merge m:1 location_code using "/restricted/projectnb/salabs/facilities/TARGET_DIM_LOCATIONS_link_codes.dta", keepusing(location_name_prov_uniq)
	/*   Result                           # of obs.
	    -----------------------------------------
	    not matched                     5,293,566
		from master                 5,292,787  (_merge==1)
		from using                        779  (_merge==2)

	    matched                       270,330,072  (_merge==3)
	    -----------------------------------------
	    */
	gen province = substr(location_name_prov_uniq,1,2)
	drop if _merge == 2 /*dropped 779*/
	count if _merge ==1 /*Note: XXX facilities did not match our list */
	drop location_name_prov_uniq location_code
	drop _merge
	bys episode_no: keep if _n==1
	compress
	count /*275622859*/
	save "../raw_data/pre2015demog/demog_oct2015.dta", replace 
	

** A.1.2 From Full275, keep only the episode_no's that correspond to plausibly HIV-related lab results

	insheet using "../raw_data/pre2015tests/Alanine_Transaminase.csv", delimiter("|") clear
		gen test_name = "ALT"
		drop if alanine_transaminase=="NULL" /*72,771,118*/
		save "../raw_data/pre2015tests/ALT.dta", replace
	insheet using "../raw_data/pre2015tests/Haemoglobin.csv", delimiter("|") clear
		gen test_name = "Hemoglobin"
		drop if haemoglobin=="NULL"	/*44,925,802 */	
		save "../raw_data/pre2015tests/hemoglobin.dta", replace
	insheet using "../raw_data/pre2015tests/CD4_Counts.csv", delimiter("|") clear
		rename episode_no episode_no2
		rename episode_id episode_no
		rename episode_no2 episode_id
		gen test_name = "CD4"
		drop if result_numeric == "NULL" /*163,402*/		
		save "../raw_data/pre2015tests/CD4.dta", replace
	insheet using "../raw_data/pre2015tests/HIV_PCR_Elisa.csv", delimiter("|") clear
		rename episode_no episode_no2
		rename episode_id episode_no
		rename episode_no2 episode_id
		gen test_name = "HIV_Elisa"
		drop if hiv_pcr_result == "NULL" & ///
			hiv_elisa_screening=="NULL" & ///
			hiv_elisa_screening_titre =="NULL" & ///
			hiv_elisa_confimatory=="NULL" & ///
			hiv_elisa_confimatory_titre=="NULL"  /* 0 */
		save "../raw_data/pre2015tests/HIV_Elisa.dta", replace
	insheet using "../raw_data/pre2015tests/HIV_Viral_Loads.csv", delimiter("|") clear
		gen test_name = "Viral Load"
		drop if result_raw=="NULL" & result_numeric=="NULL" & result_log=="NULL" /*233,010*/
		save "../raw_data/pre2015tests/VL.dta", replace
	insheet using "../raw_data/pre2015tests/CRAG.csv", delimiter("|") clear
		gen test_name = "CRAG"
		drop if crag_result=="NULL" /* 1342 */
		save "../raw_data/pre2015tests/CRAG.dta", replace
	insheet using "../raw_data/pre2015tests/Bill_crt_clearence.csv", delimiter("|") clear
		gen test_name = "Creatinine"
		drop if creatinine_clearence ==. /* 0 */
		save "../raw_data/pre2015tests/CRT.dta", replace

	*Append test data together
	use "../raw_data/pre2015tests/CRT.dta", clear
	keep  episode_no taken_date_id test_name
	append using  "../raw_data/pre2015tests/CRAG.dta", keep( episode_no taken_date_id test_name)
	append using  "../raw_data/pre2015tests/VL.dta", keep( episode_no taken_date_id test_name)
	append using  "../raw_data/pre2015tests/HIV_Elisa.dta", keep( episode_no taken_date_id test_name)
	append using  "../raw_data/pre2015tests/CD4.dta", keep( episode_no taken_date_id test_name)
	append using  "../raw_data/pre2015tests/hemoglobin.dta", keep( episode_no taken_date_id test_name)
	append using  "../raw_data/pre2015tests/ALT.dta", keep( episode_no taken_date_id test_name)
	count /*141952880*/
	tab test_name
	/*  test_name |      Freq.     Percent        Cum.
	------------+-----------------------------------
		ALT | 22,223,208       15.66       15.66
		CD4 | 30,192,947       21.27       36.93
	       CRAG |    699,612        0.49       37.42
	 Creatinine | 10,878,283        7.66       45.08
	  HIV_Elisa | 11,024,498        7.77       52.85
	 Hemoglobin | 50,068,524       35.27       88.12
	 Viral Load | 16,865,808       11.88      100.00
	------------+-----------------------------------
	      Total |141,952,880      100.00*/
	drop test_name
	bys episode_no: keep if _n==1
	count /*106,608,072*/
	
	drop taken_date_id
	
	* Merge with the pre2015 demog data
	merge 1:1 episode_no using "../raw_data/pre2015demog/demog_oct2015.dta"
	
	/*    Result                           # of obs.
	    -----------------------------------------
	    not matched                   182,238,731
		from master                 6,611,972  (_merge==1)  THESE ARE EPISODE_NO'S FOR LABS WE DO CARE ABOUT THAT ARE NOT IN THE DEMOG FILE, BUT THESE ARE NEARLY ALL AFTER OCT 2015
		from using                175,626,759  (_merge==2)  THESE ARE EPISODE_NO's ASSOCIATED WITH OTHER LAB RESULTS WE DON"T CARE ABOUT
	    matched                        99,996,100  (_merge==3)
	*/		  
	compress
	keep if _merge == 3
	drop _merge
	bys episode_no: keep if _n==1
	count /*99,996,100*/
	compress
	save "../raw_data/pre2015demog/HIV_demog_oct2015.dta", replace 




** A.1.3 Read in demographic data for HIV-related lab results for 2015-2016
	* sed -e  "s/\"//g" Bill_demo.csv  > Bill_demo_noquote.csv 
	* sed -e "s/[\\]|//g" Bill_demo_noquote.csv > demog_dec2016.csv
	* Now, in Stata...
	insheet using "../raw_data/post2015demog/demog_dec2016.csv", delimiter("|") clear names
	keep  episode_no  patient_name_lis patient_surname patient_firstname date_of_birth gender unique_patient_id national_id_no folder_no facility_code
	replace patient_name_lis = itrim(patient_name_lis)
	replace patient_surname = itrim(patient_surname)
	replace patient_firstname = itrim(patient_firstname)
	keep if episode_no !="" /*0 dropped*/
	gen len = length(episode_no) 
	keep if len==10 /*0 dropped */
	drop len
	*DOB
	gen dob_year = substr(date_of_birth,1,4)
	destring dob_year, replace force
	gen dob_month = substr(date_of_birth,6,2)
	destring dob_month, replace force
	gen dob_day = substr(date_of_birth,9,2)
	destring dob_day, replace force
	*Clean up
	drop if date_of_birth=="" /* 0 dropped */
	compress
	bys episode_no: keep if _n==1
	count /*40,033,684*/
	save "../raw_data/post2015demog/demog_dec2016.dta", replace


** A.1.4 Merge in test data to obtain specimen dates

	*Read in test data
	insheet using "../raw_data/post2015tests/Bill_bloods.csv", delimiter("|") clear
		gen test_name = "ALT"
		drop if alanine_transaminase=="NULL"   /*32,266,300*/
		save "../raw_data/post2015tests/Bill_ALT.dta", replace
	insheet using "../raw_data/post2015tests/Bill_bloods2.csv", delimiter("|") clear
		gen test_name = "Hemoglobin"
		drop if haemoglobin=="NULL"	/*24,043,913*/			
		save "../raw_data/post2015tests/Bill_hemoglobin.dta", replace
	insheet using "../raw_data/post2015tests/Bill_cd4.csv", delimiter("|") clear
		rename episode_no episode_no2
		rename episode_id episode_no
		rename episode_no2 episode_id
		gen test_name = "CD4"
	*	drop if result_numeric == "NULL" /*14*/
		save "../raw_data/post2015tests/Bill_cd4z.dta", replace
	insheet using "../raw_data/post2015tests/Bill_hiv_elisa.csv", delimiter("|") clear
		rename episode_no episode_no2
		rename episode_id episode_no
		rename episode_no2 episode_id
		gen test_name = "HIV_Elisa"
		drop if hiv_pcr_result == "NULL" & ///
			hiv_elisa_screening=="NULL" & ///
			hiv_elisa_screening_titre =="NULL" & ///
			hiv_elisa_confimatory=="NULL" & ///
			hiv_elisa_confimatory_titre=="NULL" /*0*/
	save "../raw_data/post2015tests/Bill_hiv_elisa.dta", replace
	insheet using "../raw_data/post2015tests/Bill_viral_loads.csv", delimiter("|") clear
		gen test_name = "Viral Load"
		drop if result_raw=="NULL" & result_numeric=="NULL" & result_log=="NULL" /*947*/
		save "../raw_data/post2015tests/Bill_viral_loads.dta", replace
	insheet using "../raw_data/post2015tests/CRAG.csv", delimiter("|") clear
		gen test_name = "CRAG"
		drop if crag_result=="NULL" 	/*64*/	
		save "../raw_data/post2015tests/Bill_CRAG.dta", replace
	insheet using "../raw_data/post2015tests/Bill_crt_clearance.csv", delimiter("|") clear
		gen test_name = "Creatinine"
		drop if creatinine_clearence == . 	/*0*/	
		save "../raw_data/post2015tests/Bill_crt_clearance.dta", replace

	*Append test data together
	use "../raw_data/post2015tests/Bill_crt_clearance.dta", clear
	keep  episode_no taken_date_id test_name
	append using  "../raw_data/post2015tests/Bill_CRAG.dta", keep( episode_no taken_date_id test_name)
	append using  "../raw_data/post2015tests/Bill_viral_loads.dta", keep( episode_no taken_date_id test_name)
	append using  "../raw_data/post2015tests/Bill_hiv_elisa.dta", keep( episode_no taken_date_id test_name)
	append using  "../raw_data/post2015tests/Bill_cd4.dta", keep( episode_no taken_date_id test_name)
	append using  "../raw_data/post2015tests/Bill_hemoglobin.dta", keep( episode_no taken_date_id test_name)
	append using  "../raw_data/post2015tests/Bill_ALT.dta", keep( episode_no taken_date_id test_name)
	count /*42,622,147*/
	tab test_name
	/*  test_name |      Freq.     Percent        Cum.
	------------+-----------------------------------
		ALT |  6,659,231       15.62       15.62
		CD4 |  7,041,595       16.52       32.14
	       CRAG |    678,195        1.59       33.74
	 Creatinine |  3,357,548        7.88       41.61
	  HIV_Elisa |  2,093,648        4.91       46.53
	 Hemoglobin | 14,881,618       34.92       81.44
	 Viral Load |  7,910,312       18.56      100.00
	------------+-----------------------------------
	      Total | 42,622,147      100.00 */
	drop test_name
	bys episode_no: keep if _n==1
	count /*29099119*/
	
	*Merge with demog data
	merge 1:1 episode_no using "../raw_data/post2015demog/demog_dec2016.dta"
	keep if _merge == 3 /*ANYTHING LOST HERE? WHY SO MANY DON"T MERGE?*/

	/*     Result                           # of obs.
	    -----------------------------------------
	    not matched                    10,934,575
		from master                         5  (_merge==1)
		from using                 10,934,570  (_merge==2)

	    matched                        29,099,114  (_merge==3)
	    -----------------------------------------
	*/

	drop _merge

	*NOTE: DATE_SPECIMEN IS SAME AS TAKEN DATE; THIS IS THE RIGHT ONE TO USE HERE BC IT IS WHAT"S USED FOR AGE FLAGS
	gen test_year = floor(taken_date_id*0.0001)
	gen test_month = floor((taken_date_id - test_year*10000)*.01) 
	gen test_day = taken_date_id - test_year*10000 - test_month*100
	drop if taken_date_id == . | date_of_birth == "" /* 0 dropped */
	drop taken_date_id date_of_birth
	gen age_flag = dob_month==test_month & dob_day==test_day

	keep episode_no patient_name_lis patient_surname patient_firstname gender dob_year dob_month dob_day test_year test_month test_day age_flag facility_code unique_patient_id national_id_no folder_no
	compress

	*get province, remove invalid facilities
	gen location_code = facility_code
	merge m:1 location_code using "/restricted/projectnb/salabs/facilities/TARGET_DIM_LOCATIONS_link_codes.dta", keepusing(location_name_prov_uniq)

	  /*   Result                           # of obs.
	    -----------------------------------------
	    not matched                     2,047,974
		from master                 2,041,528  (_merge==1)
		from using                      6,446  (_merge==2)

	    matched                        27,057,586  (_merge==3)
	    -----------------------------------------*/

	gen province = substr(location_name_prov_uniq,1,2)
	drop if _merge == 2
	count if _merge ==1 /*Note: 2,041,528 results with facilities did not match our list */
	drop _merge
	drop location_name_prov_uniq location_code
	compress
	count /* 29,099,114 */

	save "../raw_data/post2015demog/HIV_demog_dec2016.dta", replace
	


** A.1.5 Merge pre2015 and post2015 demog data and eliminate duplicates
	use "../raw_data/pre2015demog/HIV_demog_oct2015.dta", clear

	merge 1:1 episode_no using "../raw_data/post2015demog/HIV_demog_dec2016.dta"
	
	/*    Result                           # of obs.
	    -----------------------------------------
	    not matched                   105,900,290
		from master                88,398,638  (_merge==1)
		from using                 17,501,652  (_merge==2)

	    matched                        11,597,462  (_merge==3)
	    -----------------------------------------*/
    
	count
	*117497752
	
	drop _merge
	*replace age_flag = dob_month==test_month & dob_day==test_day
	*note, pwd = /alg_output_files/
		
	preserve
	keep episode_no unique_patient_id national_id_no folder_no
	save "HIV_demog_all_otherIDs.dta", replace
	restore
	
	drop unique_patient_id national_id_no folder_no
	save "HIV_demog_all.dta", replace

	




*************************************************************
* A.2     Pre-process / clean lab first & last names
*************************************************************


*** A.2.0 Extract names for preprocessing

use "HIV_demog_all.dta", clear

rename patient_surname last
cap gen comma_pos = strpos(patient_name_lis,",")
gen first = substr(patient_name_lis,comma_pos+1,.)

* indicator for whether there is a number >=2 digits in the name (use below as a screen for bad results)
gen num = 0
forvalues i = 0/9 {
	forvalues j= 0/9 {
		replace num = 1 if strpos(patient_name_lis,`"`i'`j'"')>0
		}
	}

keep episode_no last first num 
compress
count
*117,497,752
save "HIV_demog_all_names.dta", replace





*** A.2.1 Clean LAST NAME

use  "HIV_demog_all_names.dta", clear
	
* replace special characters used to separate names with spaces
replace last = subinstr(last,"."," ",.)
replace last = subinstr(last,","," ",.)
replace last = subinstr(last,"/"," ",.)
replace last = subinstr(last,"\"," ",.)
replace last = subinstr(last,":"," ",.)
replace last = subinstr(last,";"," ",.)
replace last = subinstr(last,"-"," ",.)

* omit all other numeric and special characters
egen _last = sieve(last), keep(alphabetic space)
replace last = trim(_last)
drop _last

* omit titles
replace last = subinstr(last,"MRS ","",.) if substr(last,1,4) == "MRS "
replace last = subinstr(last," MRS","",.) if substr(last,-4,4) == " MRS"
replace last = subinstr(last,"MR ","",.) if substr(last,1,3) == "MR "
replace last = subinstr(last," MR","",.) if substr(last,-3,3) == " MR"
replace last = subinstr(last,"MS ","",.) if substr(last,1,3) == "MS "
replace last = subinstr(last," MS","",.) if substr(last,-3,3) == " MS"
replace last = subinstr(last,"MISS ","",.) if substr(last,1,5) == "MISS "
replace last = subinstr(last," MISS","",.) if substr(last,-5,5) == " MISS"
replace last = subinstr(last,"DR ","",.) if substr(last,1,3) == "DR "
replace last = subinstr(last," DR","",.) if substr(last,-3,3) == " DR"
replace last = subinstr(last,"DOCTOR ","",.) if substr(last,1,7) == "DOCTOR "
replace last = subinstr(last," DOCTOR","",.) if substr(last,-7,7) == " DOCTOR"
replace last = subinstr(last,"PARTNER","",.)

* standardize common prefixes to last names
replace last = subinstr(last,"VD ","VAN DER ",.) if strpos(last,"VD") ==1
replace last = subinstr(last,"V D ","VAN DER ",.) if strpos(last,"V D") ==1
replace last = subinstr(last,"VD,","VAN DER ",.) if strpos(last,"VD,") ==1
replace last = subinstr(last,"V ","VAN ",.) if strpos(last,"V ") ==1
replace last = subinstr(last,"O  ","O",.) if strpos(last,"O ")==1 /* O'Brien*/
replace last = subinstr(last,"O ","O",.) if strpos(last,"O ")==1 /* O'Brien*/

* omit single initial at end of last name
replace last = trim(last)
replace last = substr(last,1,length(last)-2) if substr(last,-2,1)==" " & substr(last,-1,1) == substr(first,1,1)

* omit double initials at end
replace last = substr(last,1,length(last)-3) if substr(last,-3,1)==" " & substr(last,-2,1) == substr(first,1,1)	

* reclassify "No Name" and single initials as missing (we will exclude those missing both first and last)
replace last = "" if length(last) == 1
replace last = "" if strpos(last,"UNKNOWN")
replace last = "" if strpos(last,"UNKOWN")
replace last = "" if last == "NO NAME"
replace last = "" if last == "NO SURNAME"
replace last = "" if last == "FICTITIOUS"
replace last = "" if last == "ANONYMOUS"
replace last = "" if last == "MALE"
replace last = "" if last == "FEMALE"
replace last = "" if strpos(last,"TWIN")
replace last = "" if strpos(last,"TRIPLET")
replace last = "" if strpos(last,"UKNOWN")
replace last = "" if strpos(last,"UNKOWN")
replace last = "" if strpos(last,"UNKWON")
replace last = "" if strpos(last,"UKNWON")
replace last = "" if strpos(last,"UNKWOWN")
replace last = "" if strpos(last,"MOTHER")

* reclassify non-patient names as "INVALID"
replace last = "INVALID" if length(last) == 2 & num==1 /* defined above */
replace last = "INVALID" if strpos(last,"CODE")
replace last = "INVALID" if strpos(last,"STUDY")
replace last = "INVALID" if strpos(last,"SURVEY")
replace last = "INVALID" if strpos(last,"PROGRAMME")
replace last = "INVALID" if strpos(last,"PROJECT")
replace last = "INVALID" if strpos(last,"PATIENT")
replace last = "INVALID" if strpos(last,"EMER ")>0
replace last = "INVALID" if strpos(last,"EMERG")>0
replace last = "INVALID" if last == "WC"
replace last = "INVALID" if last == "AB"
replace last = "INVALID" if last == "EM"
replace last = "INVALID" if last == "EMT"
replace last = "INVALID" if last == "FS"
replace last = "INVALID" if last == "GP"
replace last = "INVALID" if last == "KZ"
replace last = "INVALID" if last == "LP"
replace last = "INVALID" if last == "MP"
replace last = "INVALID" if last == "EC"
replace last = "INVALID" if last == "E J"
replace last = "INVALID" if last == "NW"
replace last = "INVALID" if last == "NC"
replace last = "INVALID" if last == "HCW"
replace last = "INVALID" if last == "HIV"
replace last = "INVALID" if last == "HN"
replace last = "INVALID" if last == "GWA"
replace last = "INVALID" if last == "KLR"
replace last = "INVALID" if last == "KNS"
replace last = "INVALID" if last == "LIRA"
replace last = "INVALID" if last == "MADA"
replace last = "INVALID" if last == "MED"
replace last = "INVALID" if last == "M W"
replace last = "INVALID" if last == "MGMH"
replace last = "INVALID" if last == "MW"
replace last = "INVALID" if last == "MMMRH"
replace last = "INVALID" if last == "NEW"
replace last = "INVALID" if last == "R K K"
replace last = "INVALID" if strpos(last,"STAFF")>0
replace last = "INVALID" if last == "OT"
replace last = "INVALID" if last == "TCDBS"
replace last = "INVALID" if last == "TCFL"
replace last = "INVALID" if last == "TCODBS"
replace last = "INVALID" if last == "TCOFL"
replace last = "INVALID" if last == "VCT"
replace last = "INVALID" if strpos(last,"CTK")>0
replace last = "INVALID" if last == "LADY"
replace last = "INVALID" if last == "SARI"
replace last = "INVALID" if last == "ABGB"
replace last = "INVALID" if last == "TEST"
replace last = "INVALID" if last == "SECURITY"
replace last = "INVALID" if strpos(last,"SAMPLE")
replace last = "INVALID" if strpos(last,"NHLS")
replace last = "INVALID" if strpos(last,"HAEM")
replace last = "INVALID" if strpos(last,"TRIAL")
replace last = "INVALID" if strpos(last,"CYCLE")
replace last = "INVALID" if strpos(last,"COAG")
replace last = "INVALID" if strpos(last,"EQA")
replace last = "INVALID" if strpos(last,"OLD LOT")
replace last = "INVALID" if strpos(last,"CHEMISTRY")
replace last = "INVALID" if strpos(last,"VALID")
replace last = "INVALID" if strpos(last,"HUMAN")
replace last = "INVALID" if strpos(last,"UNIDENTIFIED")
replace last = "INVALID" if strpos(last,"UNINDENTIFIED")
replace last = "INVALID" if strpos(last,"UNIDENTIFEID")
replace last = "INVALID" if strpos(last,"UNIDENTIFY")
replace last = "INVALID" if strpos(last,"BLOOD GAS")
replace last = "INVALID" if strpos(last,"BLOODGAS")
replace last = "INVALID" if strpos(last,"MACHINE")
replace last = "INVALID" if strpos(last,"E COLI")
replace last = "INVALID" if strpos(last,"ECOLI")
replace last = "INVALID" if strpos(last,"CONTROL")
replace last = "INVALID" if strpos(last,"WATER")
replace last = "INVALID" if strpos(last,"SYSTEM")
replace last = "INVALID" if strpos(last,"BYPASSED")
replace last = "INVALID" if last == "ZAA"
replace last = "INVALID" if last == "ZAB"
replace last = "INVALID" if last == "ZAC"
replace last = "INVALID" if last == "ZAD"
replace last = "INVALID" if last == "ZAE"
replace last = "INVALID" if last == "ZAF"
replace last = "INVALID" if last == "ZAG"
replace last = "INVALID" if last == "ZAH"
replace last = "INVALID" if last == "ZAI"
replace last = "INVALID" if last == "ZAJ"
replace last = "INVALID" if last == "POS QC"
replace last = "INVALID" if last == "BLACK MALE"
replace last = "INVALID" if last == "BLACKMALE"
replace last = "INVALID" if last == "BLACK FEMALE"
replace last = "INVALID" if last == "BLACKFEMALE"
replace last = "INVALID" if last == "WHITE MALE"
replace last = "INVALID" if last == "WHITEMALE"
replace last = "INVALID" if last == "WHITE FEMALE"
replace last = "INVALID" if last == "WHITEFEMALE"
replace last = "INVALID" if last == "INDIAN MALE"
replace last = "INVALID" if last == "INDIANMALE"
replace last = "INVALID" if last == "INDIAN FEMALE"
replace last = "INVALID" if last == "INDIANFEMALE"
replace last = "INVALID" if last == "COLORED MALE"
replace last = "INVALID" if last == "COLOREDMALE"
replace last = "INVALID" if last == "COLORED FEMALE"
replace last = "INVALID" if last == "COLOREDFEMALE"
replace last = "INVALID" if last == "FOODHANDLER"
replace last = "INVALID" if last == "FOOD HANDLER"
replace last = "INVALID" if last == "DOH"
replace last = "INVALID" if last == "DO H"
replace last = "INVALID" if last == "NDOH"
replace last = "INVALID" if last == "PATIENT"
replace last = "INVALID" if last == "UNKNWON"
replace last = "INVALID" if last == "TEST"
replace last = "INVALID" if last == "MEDITECH"
replace last = "INVALID" if last == "MILK"
replace last = "INVALID" if last == "INVALID"
replace last = "INVALID" if last == "URINE"
replace last = "INVALID" if last == "BLOOD"
replace last = "INVALID" if strpos(last,"CARDIAC")
replace last = "INVALID" if last == "NO FORM"

* >2 letter name, no vowels
replace last = "INVALID" if length(last)>=2 & strpos(last,"A")==0 & strpos(last,"E")==0 & strpos(last,"I")==0 & strpos(last,"O")==0 & strpos(last,"U")==0 & strpos(last,"Y")==0 

replace last = trim(last)
replace last = itrim(last)

save  "HIV_demog_all_names.dta", replace






*** A.2.2 CLEAN FIRST NAME

use  "HIV_demog_all_names.dta", clear
* replace numeric and special characters

* replace special characters used to separate names with spaces
replace first = subinstr(first,"."," ",.)
replace first = subinstr(first,","," ",.)
replace first = subinstr(first,"/"," ",.)
replace first = subinstr(first,"\"," ",.)
replace first = subinstr(first,":"," ",.)
replace first = subinstr(first,";"," ",.)
replace first = subinstr(first,"-"," ",.)

* omit all other numeric and special characters
egen _first = sieve(first), keep(alphabetic space)
replace first = trim(_first)
replace first = itrim(first)
drop _first

* Omit initials 
*  1 letter initials at beginning
replace first = substr(first,1,1)+ substr(first,3,.) if substr(first,2,1)==" " & substr(first,4,1)==" "
replace first = substr(first,3,.) if strpos(first," ")==2 & substr(first,1,1)==substr(first,3,1)
replace first = substr(first,3,.) if strpos(first," ")==2 & substr(first,1,1)== substr(word(first,-1),1,1)
replace first = substr(first,1,length(first)-2) if strpos(first," ")==length(first)-2 & substr(first,-1,1)==substr(first,1,1)
*  2 letter initials at beginning
replace first = substr(first,4,.) if strpos(first," ")==3 & substr(first,1,1)==substr(first,4,1)
replace first = substr(first,4,.) if strpos(first," ")==3 & substr(first,2,1)==substr(first,4,1)
replace first = substr(first,4,.) if strpos(first," ")==3 & substr(first,2,1)== substr(word(first,-1),1,1)
*  3 letter initials at beginning
replace first = substr(first,5,.) if wordcount(first)==4 & strpos(first," ")==4 & substr(first,1,1)==substr(first,5,1)
*  1 letter middle initial at end
replace first = substr(first,1,length(first)-2) if strpos(first," ")==length(first)-2 & length(first) >2
* drop 1 letter initial in middle
local alphabet "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z"
foreach a of local alphabet {
	replace first = substr(first,1,strpos(first,`" `a' "')) + substr(first,strpos(first,`" `a' "')+4,.) ///
	if strpos(first,`" `a' "') > 0 & substr(first,strpos(first,`" `a' "')+2,1)==`"`a'"'	
	}	

*drop non-sensical pieces of names 
	replace first = regexr(first,"AAA(.?)P","")
	replace first = regexr(first,"AAG(.?)A","")
	replace first = regexr(first,"AAH(.?)A","")

* drop titles
replace first = subinstr(first,"MRS ","",.) if substr(first,1,4) == "MRS "
replace first = subinstr(first," MRS","",.) if substr(first,-4,4) == " MRS"
replace first = subinstr(first,"MM ","",.) if substr(first,1,3) == "MM "
replace first = subinstr(first," MM","",.) if substr(first,-3,3) == " MM"
replace first = subinstr(first,"MR ","",.) if substr(first,1,3) == "MR "
replace first = subinstr(first," MR","",.) if substr(first,-3,3) == " MR"
replace first = subinstr(first,"MISS ","",.) if substr(first,1,5) == "MISS "
replace first = subinstr(first," MISS","",.) if substr(first,-5,5) == " MISS"
replace first = subinstr(first,"MS ","",.) if substr(first,1,3) == "MS "
replace first = subinstr(first," MS","",.) if substr(first,-3,3) == " MS"
replace first = subinstr(first,"DR ","",.) if substr(first,1,3) == "DR "
replace first = subinstr(first," DR","",.) if substr(first,-3,3) == " DR"
replace first = subinstr(first,"DOCTOR ","",.) if substr(first,1,7) == "DOCTOR "
replace first = subinstr(first," DOCTOR","",.) if substr(first,-7,7) == " DOCTOR"

replace first = trim(first)
replace first = itrim(first)

* reclassify "No Name" and single initials as missing (we will exclude those missing both first and last)
replace first = "" if strpos(first,"NT STATED")
replace first = "" if strpos(first,"NOT STATED")
replace first = "" if strpos(first,"NOTSTATED")
replace first = "" if strpos(first,"NO ANSWER")
replace first = "" if strpos(first,"NOANSWE")
replace first = "" if strpos(first,"NO LAST")
replace first = "" if strpos(first,"NOLAST")
replace first = "" if strpos(first,"NO FIRS")
replace first = "" if strpos(first,"NOFIRS")
replace first = "" if strpos(first,"NO NAME")
replace first = "" if strpos(first,"NONAME")
replace first = "" if strpos(first,"NOT SUPPLIED")
replace first = "" if strpos(first,"UNKNOWN")
replace first = "" if strpos(first,"UNKOWN")
replace first = "" if strpos(first,"FICTITIOUS")
replace first = "" if strpos(first,"MOTHER TO")
replace first = "" if strpos(first,"MOTHER OF")
replace first = "" if strpos(first,"FATHER TO")
replace first = "" if strpos(first,"FATHER OF")
replace first = "" if strpos(first,"BABY")
replace first = "" if strpos(first,"TWIN")
replace first = "" if first == "MALE"
replace first = "" if first == "FEMALE"
replace first = "" if strpos(first,"TWIN")
replace first = "" if strpos(first,"TRIPLET")
replace first = "" if strpos(first,"UKNOWN")
replace first = "" if strpos(first,"UNKOWN")
replace first = "" if strpos(first,"UNKWON")
replace first = "" if strpos(first,"UKNWON")
replace first = "" if strpos(first,"UNKWOWN")
replace first = "" if strpos(first,"MOTHER")

* reclassify non-patient names as INVALID
replace first = "INVALID" if strpos(first,"STUDY")
replace first = "INVALID" if strpos(first,"PROGRAMME")
replace first = "INVALID" if strpos(first,"PROJECT")
replace first = "INVALID" if strpos(first,"SURVEY")
replace first = "INVALID" if strpos(first,"SUVREY")
replace first = "INVALID" if strpos(first,"SRUVEY")
replace first = "INVALID" if first == "HSRC"
replace first = "INVALID" if first == "FOSA" & last == "FOSA"
replace first = "INVALID" if first == "ABAE" & last == "ABAE"
replace first = "INVALID" if first == "ABGB" & last == "ABGB"
replace first = "INVALID" if first == "OHSC"
replace first = "INVALID" if strpos(first,"CTKN") > 0
replace first = "INVALID" if strpos(first,"SAMPLE")
replace first = "INVALID" if strpos(first,"NHLS")
replace first = "INVALID" if strpos(first,"HAEM")
replace first = "INVALID" if strpos(first,"TRIAL")
replace first = "INVALID" if strpos(first,"CYCLE")
replace first = "INVALID" if strpos(first,"COAG")
replace first = "INVALID" if strpos(first,"EQA")
replace first = "INVALID" if strpos(first,"OLD LOT")
replace first = "INVALID" if strpos(first,"CHEMISTRY")
replace first = "INVALID" if strpos(first,"VALID")
replace first = "INVALID" if strpos(first,"HUMAN")
replace first = "INVALID" if strpos(first,"UNIDENTIFIED")
replace first = "INVALID" if strpos(first,"UNINDENTIFIED")
replace first = "INVALID" if strpos(first,"UNIDENTIFEID")
replace first = "INVALID" if strpos(first,"UNIDENTIFY")
replace first = "INVALID" if strpos(first,"BLOOD GAS")
replace first = "INVALID" if strpos(first,"BLOODGAS")
replace first = "INVALID" if strpos(first,"MACHINE")
replace first = "INVALID" if strpos(first,"E COLI")
replace first = "INVALID" if strpos(first,"ECOLI")
replace first = "INVALID" if strpos(first,"CONTROL")
replace first = "INVALID" if strpos(first,"WATER")
replace first = "INVALID" if strpos(first,"SYSTEM")
replace first = "INVALID" if strpos(first,"BYPASSED")
replace first = "INVALID" if first == "ZAA"
replace first = "INVALID" if first == "ZAB"
replace first = "INVALID" if first == "ZAC"
replace first = "INVALID" if first == "ZAD"
replace first = "INVALID" if first == "ZAE"
replace first = "INVALID" if first == "ZAF"
replace first = "INVALID" if first == "ZAG"
replace first = "INVALID" if first == "ZAH"
replace first = "INVALID" if first == "ZAI"
replace first = "INVALID" if first == "ZAJ"
replace first = "INVALID" if first == "POS QC"
replace first = "INVALID" if first == "BLACK MALE"
replace first = "INVALID" if first == "BLACKMALE"
replace first = "INVALID" if first == "BLACK FEMALE"
replace first = "INVALID" if first == "BLACKFEMALE"
replace first = "INVALID" if first == "WHITE MALE"
replace first = "INVALID" if first == "WHITEMALE"
replace first = "INVALID" if first == "WHITE FEMALE"
replace first = "INVALID" if first == "WHITEFEMALE"
replace first = "INVALID" if first == "INDIAN MALE"
replace first = "INVALID" if first == "INDIANMALE"
replace first = "INVALID" if first == "INDIAN FEMALE"
replace first = "INVALID" if first == "INDIANFEMALE"
replace first = "INVALID" if first == "COLORED MALE"
replace first = "INVALID" if first == "COLOREDMALE"
replace first = "INVALID" if first == "COLORED FEMALE"
replace first = "INVALID" if first == "COLOREDFEMALE"
replace first = "INVALID" if first == "FOODHANDLER"
replace first = "INVALID" if first == "FOOD HANDLER"
replace first = "INVALID" if first == "DOH"
replace first = "INVALID" if first == "DO H"
replace first = "INVALID" if first == "NDOH"
replace first = "INVALID" if first == "PATIENT"
replace first = "INVALID" if first == "UNKNWON"
replace first = "INVALID" if first == "TEST"
replace first = "INVALID" if first == "MEDITECH"
replace first = "INVALID" if first == "MILK"
replace first = "INVALID" if first == "INVALID"
replace first = "INVALID" if first == "URINE"
replace first = "INVALID" if first == "BLOOD"
replace first = "INVALID" if strpos(first,"CARDIAC")
replace first = "INVALID" if first == "NO FORM"
replace first = "INVALID" if first == "FORM"

* Screen for invalid names based on absence of other name, presence of double digit number, and length of name
replace first = "INVALID" if length(first) == 1 & num==1 & last == "" /* defined above */
replace first = "INVALID" if length(first) == 2 & num==1 & last == "" 
replace first = "INVALID" if length(first) == 3 & num==1 & last == "" 
replace first = "INVALID" if length(first) == 4 & num==1 & last == "" 
replace last = "INVALID" if first == ""  & length(last)==1 & num==1
replace last = "INVALID" if first == ""  & length(last)==2 & num==1
replace last = "INVALID" if first == ""  & length(last)==3 & num==1
replace last = "INVALID" if first == ""  & length(last)==4 & num==1

*get rid of some non-sensical names
replace first = "INVALID" if first == "ABA" | first == "ACA" | first == "AC A" | first == "AB A" | first == "AD A"
replace last = "INVALID" if last == "ABA" | last == "ACA" | last == "AC A" | last == "AB A" | last == "AD A"
replace last = "INVALID" if last == "ADA" & first == "ADA"
replace last = "INVALID" if last == "A DA" & first == "ADA"
replace last = "INVALID" if last == "AD A" & first == "AD A"	

save  "HIV_demog_all_names.dta", replace



*** A.2.3 Remove one and two letter segments of first names IF MORE THAN ONE segment
split first, p("-" " ") gen(_f)
gen length_fi = .
gen max_length_fi = 0
forvalues i = 1/7 {
	replace max_length_fi = length(_f`i') if length(_f`i') > max_length_fi & length(_f`i') < .
	}
forvalues i = 1/7 {
	replace length_fi = length(_f`i')
	replace _f`i' = "" if length_fi < 3 & _f2 !="" & max_length_fi >=3
	}
replace first = _f1+" "+_f2+" "+_f3+" "+_f4+" "+_f5+" "+_f6+" "+_f7
replace first = trim(itrim(first))
drop _f*
drop length_fi max_length_fi

save "HIV_demog_all_names.dta", replace






***************************************
* A.3	Additional Pre-Processing
***************************************

*** A.3.1 Distribution of characteristics (identify rare names etc)

use "HIV_demog_all_names.dta", clear

merge 1:1 episode_no using "HIV_demog_all.dta", keepusing(episode_no gender dob_year dob_month dob_day test_year test_month test_day age_flag facility_code province)
keep if _merge == 3
drop _merge

count
*117,497,752

* recode a few missing DOBs as 1800-01-01
replace dob_year = 1800 if dob_year==. /*31*/
replace dob_month = 1 if dob_month==. /*31*/
replace dob_day = 1 if dob_day==. /*31*/

* recode a few missing and invalid genders as "U"
replace gender = "U" if gender != "M" & gender != "F" /*83*/

* distribution of characteristics
bys dob_year: gen p_YOB = _N/r(N)  /*number of episodes with same YOB*/
bys first: gen p_first = _N/r(N)  /*number of episodes with same first name*/
bys last: gen p_last = _N/r(N)  /*number of episodes with same last name*/
bys gender: gen p_gender = _N/r(N)  /*number of episodes with same gender*/
bys facility: gen p_facility = _N/r(N)  /*number of episodes with same facility*/
bys province: gen p_province = _N/r(N)  /*number of episodes with same province*/
compress
save "HIV_demog_all_names.dta", replace	



*** A.3.2 Separate out invalid lab episodes or episodes with missing information

* DROP INVALID RESULTS OR RESULTS WITH MISSING DOB AND MISSING LAST;
* THESE RESULTS WILL BE NEARLY IMPOSSIBLE TO MATCH AND MOST ARE INVALID. 
* MOVE RESULTS WITH MISSING DOB, BUT PRESENT FIRST/LAST, TO DATASET FOR EXACT MATCHING LATER.
* MOVE RESULTS WITH MISSING FIRST, BUT PRESENT DOB/LAST, TO DATASET FOR EXACT MATCHING LATER.

** Drop invalid results
use  "HIV_demog_all_names.dta", clear
count
*117497752
drop if first == "INVALID" | last == "INVALID"
*(434024 observations deleted)
drop if first == "" & num == 1
*(73410 observations deleted)
drop if dob_year == 1800 & num == 1
*(31184 observations deleted)
drop if first == "" & dob_year == 1800
*(68228 observations deleted)
drop if last == ""
*(321620 observations deleted)
count
*116569286

save "HIV_demog_all_names.dta", replace   
* THIS DATASET CONTAINS THE IDENTIFYING INFO, PROBABILITIES, AND FACILITY INFO
* FOR ALL VALID LAB RESULTS.



** Partition into "clean", "invalid", "dob1800", and "firstNA"
use "HIV_demog_all_names.dta", clear

* invalid records
merge 1:1 episode_no using "HIV_demog_all.dta", keepusing(dob_year)

preserve
keep if _merge == 2
drop _merge
save "HIV_demog_all_INVALID.dta", replace
restore

drop if _merge == 2
*(928466 observations deleted)
drop _merge

* dob1800 records
preserve
keep if dob_year == 1800
save "HIV_demog_all_dob1800.dta", replace
restore

drop if dob_year == 1800
*(4192130 observations deleted)

* firstNA records
preserve
keep if first == ""
save "HIV_demog_all_firstNA.dta", replace
restore

drop if first == ""
*(996900 observations deleted)

* clean records
count
*111,380,256

save "HIV_demog_all_clean.dta", replace









*****************************************
* A.4 Create Exact Match Identifiers
*****************************************

*** A.4.1 Create Exact Match Identifiers, but retain at level of Episode_No

* THIS DATASET IS AT THE LEVEL OF THE EPISODE NUMBER. I NOW CREATE IDENTIFIERS FOR THE 
* DOB-LAST-FIRST EM_ID AND FOR THE DOB-LAST-FIRST-GENDER-FACILITY EM_ID_PLUS
use "HIV_demog_all_clean.dta", clear
*keep last first dob_year dob_month dob_day gender facility
sort dob_year dob_month dob_day last first gender facility 
*by DOB last first: gen labs_per_patient = _N
preserve
	by dob_year dob_month dob_day last first: keep if _n ==1 
	gen double EM_ID = _n
	label var EM_ID "Exact matched ID - DOB / Last / First"
	save "HIV_EM_ID.dta", replace
restore
merge m:1 dob_year dob_month dob_day last first using "HIV_EM_ID.dta"
drop _merge
save "HIV_demog_all_clean.dta", replace

use "HIV_demog_all_clean.dta", clear
cap drop EM_ID_plus
preserve
	bys dob_year dob_month dob_day last first gender facility: keep if _n ==1 
	gen double EM_ID_plus = _n
	label var EM_ID_plus "Exact matched ID - DOB / Last / First / Gender / Facility" 
	drop num episode_no
	save "HIV_EM_ID_plus.dta", replace
	outsheet using "HIV_EM_ID_plus.csv", replace noquote		
restore
merge m:1 dob_year dob_month dob_day last first gender facility using "HIV_EM_ID_plus.dta" 
drop _merge
count /* n=111,380,256*/
save "HIV_demog_all_clean.dta", replace



*** A.4.2 Create First/Last/DOB extract for fuzzyMatch on the GPUs
use "HIV_demog_all_clean.dta", clear
keep dob_year dob_month dob_day last first EM_ID
bys EM_ID: keep if _n==1
order EM_ID dob_year dob_month dob_day last first
save "HIV_demog_EM_ID.dta", replace


*** Export for Katia/fuzzyMatch.cu as csv file
outsheet using "HIV_demog_EM_ID.csv", replace noquote

timer off 1
timer list 1




