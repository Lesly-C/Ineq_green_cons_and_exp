*******************************************************
* Replication code
* Authors: Lesly Cassin, Paolo Melindi-Ghidi, Fabien Prieur
* Article: "The impact of income inequality on public environmental 
*           expenditure with green consumers"
* Journal: European Journal of Political Economy (forthcoming)
*
* Corresponding Author: lesly.cassin@univ-lorraine.fr
* Please cite as: [à compléter avec référence finale]
*******************************************************

*---------------------------------------------------*
* Prerequisites
*---------------------------------------------------*
* 1. Ensure that raw datasets are inside "initial_datasets/"
* 2. Ensure "stata_files/" folder exists for outputs
* 3. Run this file before "data_analysis.do"
*---------------------------------------------------*

**# Bookmark #1: Complete list of countries 

// demo_pjan_selected.xlsx : Eurostat population data (subset of countries)
// We use population as the base dataset since it is the most complete 
// in terms of country coverage (ensures consistent country list)

// list_countries.xlsx : Reference list of countries with iso2, iso3, country name
// Used to merge and keep consistent country codes across datasets

import excel "initial_datasets/eurostat_data/demo_pjan_selected.xlsx", firstrow case(lower) clear
keep geo
save "stata_files/countries_eurostat.dta", replace

import excel "initial_datasets/eurostat_data/list_countries.xlsx", firstrow case(lower) clear
merge m:m geo using "stata_files/countries_eurostat.dta"
drop if _merge != 3
drop _merge
save "stata_files/countries_eurostat.dta", replace


**# Bookmark #2: Public expenditures

/*---------------------------------------------------
GOV_10A_EXP - General Government Total Expenditure
Source: gov_10a_exp_total_selected.xlsx
Unit: Million Euro
---------------------------------------------------*/

import excel "initial_datasets/eurostat_data/gov_10a_exp_total_selected.xlsx", firstrow clear
destring gov_10_a_tot_meur*, replace force
reshape long gov_10_a_tot_meur, i(geo) j(time)
save "stata_files/gov_10a_exp_GG_tot.dta", replace   // General government total expenditure


/*---------------------------------------------------
GOV_10A_EXP - Environmental Protection Expenditure
Source: gov_10a_exp_GG_selected.xlsx
Unit: Million Euro
---------------------------------------------------*/

import excel "initial_datasets/eurostat_data/gov_10a_exp_GG_selected.xlsx", firstrow clear
destring gen_gov_exp_env_prot*, replace force
reshape long gen_gov_exp_env_prot, i(geo) j(time)

// Replace invalid/missing values (<0) with missing
replace gen_gov_exp_env_prot = . if gen_gov_exp_env_prot < 0   // 2 values replaced
save "stata_files/gov_10a_exp_GG.dta", replace   // General government env. protection expenditure


/*---------------------------------------------------
Merge: Total expenditure + Env. protection expenditure
---------------------------------------------------*/

merge m:m geo time using "stata_files/gov_10a_exp_GG_tot.dta"
drop _merge
save "stata_files/gov_10a_exp.dta", replace


**# Bookmark #3: Main independent variables

/*---------------------------------------------------
GDP (NAMA_10_GDP)
Source: nama_10_gdp_selected_meur.xlsx
Unit: million euro (current)
---------------------------------------------------*/
import excel "initial_datasets/eurostat_data/nama_10_gdp_selected_meur.xlsx", firstrow clear
destring meur_gdp*, replace force
reshape long meur_gdp, i(geo) j(time)
save "stata_files/nama_10_gdp.dta", replace


/*---------------------------------------------------
Population (DEMO_PJAN)
Source: demo_pjan_selected.xlsx
Unit: persons (headcount)
---------------------------------------------------*/
import excel "initial_datasets/eurostat_data/demo_pjan_selected.xlsx", firstrow clear 
destring pop_total*, replace force 
reshape long pop_total, i(geo) j(time)
save "stata_files/demo_pjan.dta", replace


/*---------------------------------------------------
Inequalities: Gini
Source: ilc_di12_selected.xlsx
Unit: 0–100 scale
---------------------------------------------------*/
import excel "initial_datasets/eurostat_data/ilc_di12_selected.xlsx", firstrow clear 
destring gini*, replace force 
reshape long gini, i(geo) j(time)
save "stata_files/ilc_di12.dta", replace


/*---------------------------------------------------
Inequalities: Income by decile groups
Source: ilc_di01_D1_selected.xlsx ... ilc_di01_D10_selected.xlsx
Unit: % of total GDP in each decile
---------------------------------------------------*/

* Step 1: Start with D1
import excel "initial_datasets/eurostat_data/ilc_di01_D1_selected.xlsx", firstrow clear 
destring inc_dec01*, replace force 
reshape long inc_dec01, i(geo) j(time)
save "stata_files/ilc_di01_eur.dta", replace

* Step 2: Loop over D2–D10 and merge
foreach d of numlist 2/10 {
    local D = string(`d', "%02.0f")   // 02, 03, ..., 10
    
    import excel "initial_datasets/eurostat_data/ilc_di01_D`d'_selected.xlsx", firstrow clear
    destring inc_dec`D'*, replace force 
    reshape long inc_dec`D', i(geo) j(time)
    
    merge m:m geo time using "stata_files/ilc_di01_eur.dta"
    drop _merge
    save "stata_files/ilc_di01_eur.dta", replace
}



**# Bookmark #4: Environmental variables 

/*---------------------------------------------------
EPI indicators (REC, SNM, GHP)
Source: initial_datasets/epi_indicators.csv
Unit: index values (n/a converted to missing)
---------------------------------------------------*/

foreach var in REC SNM GHP {
    insheet using "initial_datasets/epi_indicators/`var'_ind_na.csv", clear
    
    drop code country
    reshape long `=lower("`var'")'ind, i(iso) j(j)
    rename j time 
    rename `=lower("`var'")'ind `=lower("`var'")'
    destring `=lower("`var'")', replace force   // n/a → missing
    
    capture confirm file "stata_files/environmental_var.dta"
    if _rc != 0 {
        save "stata_files/environmental_var.dta", replace
    }
    else {
        merge m:m iso time using "stata_files/environmental_var.dta"
        drop _merge
        save "stata_files/environmental_var.dta", replace
    }
}


/*---------------------------------------------------
Organic food consumption (FIBL dataset)
Unit: EUR per capita
---------------------------------------------------*/
import excel "initial_datasets/fibl_dataset_organic_food_consumption_pc_2023.xlsx", firstrow clear 
destring organic_food_consumption_pc, replace force 
merge m:m iso time using "stata_files/environmental_var.dta"
drop _merge
save "stata_files/environmental_var.dta", replace


**# Bookmark #5: Merge all files

/* Merge all constructed files by country and year
   Using m:m since some sources contain multiple rows per (geo,time) 
   (NOTE: check for duplicates with isid geo time at the end) 
*/

use "stata_files/countries_eurostat.dta", replace 

merge m:m geo using "stata_files/gov_10a_exp.dta"
drop if _merge!=3 
drop _merge

merge m:m geo time using "stata_files/nama_10_gdp.dta"
drop _merge

merge m:m geo time using "stata_files/demo_pjan.dta" 
drop _merge

merge m:m geo time using "stata_files/ilc_di12.dta"
drop _merge

merge m:m geo time using "stata_files/ilc_di01_eur.dta"
drop _merge

recast str2 geo, force

merge m:m iso time using "stata_files/environmental_var.dta"
drop _merge

// Drop rows with missing country info (incomplete merges)
drop if country==""

save "stata_files/dataset.dta", replace

// Optional: save in older format for compatibility with Stata 13+
// saveold "stata_files/dataset_V13.dta", version(13) replace

