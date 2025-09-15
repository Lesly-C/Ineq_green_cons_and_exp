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
**# Bookmark #1 - Prerequisites
*---------------------------------------------------*
* 1. Ensure "dataset.dta" inside the folder "stata_files" or optionally run dataset_creation.do
* 2. Create folders "figures" and "tables" for graph and table exports
* 3. Ensure required user-written packages are installed
*---------------------------------------------------*

// Load dataset
use "stata_files/dataset.dta", clear

*---------------------------------------------------*
* Installation of required packages (auto-check)
*---------------------------------------------------*
cap which xttest3
if _rc ssc install xttest3

cap which xtserial
if _rc net install st0039, from("http://www.stata-journal.com/software/sj3-2")

cap which xtscc
if _rc ssc install xtscc

cap which xtabond2
if _rc ssc install xtabond2

cap which xtgcause
if _rc net install st0455, from("http://www.stata-journal.com/software/sj16-3")

cap which ivreg2
if _rc ssc install ivreg2, replace

cap which xtoverid
if _rc ssc install xtoverid

cap which utest
if _rc ssc install utest, replace

cap which xtcsd
if _rc ssc install xtcsd, replace

cap which esttab
if _rc ssc install estout




*---------------------------------------------------*
**# Bookmark #2 - Panel declaration
*---------------------------------------------------*

{ // Panel declaration

	/*---------------------------------------------------
	Country selection
	---------------------------------------------------*/
	// Keep only selected European countries
	keep if country=="Austria" || country=="Belgium" || country=="Bulgaria" || country=="Croatia" || ///
			 country=="Cyprus" || country=="Czechia" || country=="Denmark" || country=="Estonia" || ///
			 country=="Finland" || country=="France" || country=="Germany" || country=="Greece" || ///
			 country=="Hungary" || country=="Iceland" || country=="Ireland" || country=="Italy" || ///
			 country=="Latvia" || country=="Lithuania" || country=="Luxembourg" || country=="Malta" || ///
			 country=="Netherlands (the)" || country=="Norway" || country=="Poland" || country=="Portugal" || ///
			 country=="Romania" || country=="Slovakia" || country=="Slovenia" || country=="Spain" || ///
			 country=="Sweden" || country=="Switzerland" || ///
			 country=="United Kingdom of Great Britain and Northern Ireland (the)"

	/*
	Included in the full analysis:
	Austria, Belgium, Bulgaria, Switzerland, Czechia, Germany, Denmark, Estonia,
	Greece, Spain, Finland, France, Croatia, Hungary, Ireland, Italy, Lithuania,
	Luxembourg, Latvia, Netherlands (the), Norway, Poland, Romania, Sweden, Slovenia

	Excluded from the full analysis (missing data):
	Cyprus, Iceland, Malta, Portugal, Slovakia, United Kingdom
	*/

	// Restrict time period
	drop if time > 2021   // exclude Covid years


	/*---------------------------------------------------
	Panel declaration
	---------------------------------------------------*/
	encode geo, gen(geo1)
	sort geo1 time
	xtset geo1 time


	/*---------------------------------------------------
	Dependent variables
	---------------------------------------------------*/
	gen GEEP_pc     = gen_gov_exp_env_prot*1000000/pop_total
	gen log_GEEP_pc = log(GEEP_pc)

	// GFE grouping: based on mean depvar
	by geo1: egen mean_depvar = mean(log_GEEP_pc)
	xtile gfe_group = mean_depvar, n(4)
	drop mean_depvar

	/*---------------------------------------------------
	Main independent variables
	---------------------------------------------------*/
	// Gini
	gen gini_2 = gini^2

	gen L5_gini   = L5.gini
	gen L5_gini_2 = L5.gini_2

	// Palma ratio
	gen palma   = inc_dec10 / (inc_dec01 + inc_dec02 + inc_dec03 + inc_dec04)
	gen palma_2 = palma^2

	gen L5_palma   = L5.palma
	gen L5_palma_2 = L5.palma_2

	// GDP per capita
	gen gdp_pcap_eur     = meur_gdp*10^6 / pop_total
	gen log_gdp_pcap_eur = log(gdp_pcap_eur)

	// Organic consumption
	gen log_organic = log(organic_food_consumption_pc)


	/*---------------------------------------------------
	Controls
	---------------------------------------------------*/
	// Public expenditures share
	gen share_pub_exp = gov_10_a_tot_meur / meur_gdp

	// Population (level + growth)
	gen log_pop = log(pop_total)
	by geo1: gen pop_growth_1 = pop_total/pop_total[_n-1]  // variation (ratio)


	/*---------------------------------------------------
	Label variables
	---------------------------------------------------*/
	// Inequalities
	label variable palma       "$ Palma $"
	label variable palma_2     "$ Palma^2 $"
	label variable L5_palma    "$ Palma_{t-5} $"
	label variable L5_palma_2  "$ Palma_{t-5}^2 $"

	label variable gini        "$ Gini $"
	label variable gini_2      "$ Gini^2 $"
	label variable L5_gini     "$ Gini_{t-5} $"
	label variable L5_gini_2   "$ Gini_{t-5}^2 $"

	// Economy
	label variable log_gdp_pcap_eur "$ GDP $"
	label variable share_pub_exp    "$ Pub.Exp/GDP $"

	// Demography
	label variable pop_growth_1 "$ \Delta pop $"

	// Environment
	label variable log_organic "$ Organic $" 
	label variable rec   "$ REC $" 
	label variable snm   "$ SNM $" 
	label variable ghp   "$ GHP $"

}


/*---------------------------------------------------
Save final dataset with the new variables
---------------------------------------------------*/
sort geo1 time
save "stata_files/database_for_analysis.dta", replace   


*---------------------------------------------------*
**# Bookmark #3 - Stylized facts
*---------------------------------------------------*

{ // Stylized Facts - Section 3


	/* Table 1: Summary statistics


	Table 1 reports descriptive statistics for all main variables used in the empirical analysis. 
	These statistics are obtained with the summarize command in Stata.

	Variables included:
	  - GEEP_pc : General government EEP per capita
	  - palma : Palma ratio
	  - gini : Gini index
	  - gdp_pcap_eur : GDP per capita (EUR)
	  - share_pub_exp : Public expenditures / GDP
	  - pop_total : Population (millions)
	  - rec : Renewable energy consumption (percentage)
	  - snm : Sustainable nutrient management (score)
	  - ghp : Greenhouse gas performance (score)
	  - organic_food_consumption_pc : Organic food consumption per capita

	The command below produces N, mean, std. dev., min, max which are used to construct Table 1 in the article.
	-------------------------------------------------------------*
	*/

	summarize GEEP_pc palma gini gdp_pcap_eur share_pub_exp pop_total ///
			  rec snm ghp organic_food_consumption_pc


	* Figure 1: Organic consumption over time
	twoway qfitci organic_food_consumption_pc time, ///
		xtitle("Time") ///
		ytitle("Organic consumption per capita (EUR)")
	graph export "figures/organic_time_2025.eps", replace

	* Figure 2: Environmental expenditure vs income indicators
	// NOTE: "GEEP_pc" = Dependent variable created above (gen_gov_exp_env_prot per capita)
	* 2a. vs GDP per capita
	twoway scatter GEEP_pc gdp_pcap_eur, ///
		ytitle("Gen. Gov. EEP per capita (EUR)", size(large)) ///
		xtitle("GDP per capita", size(large)) ///
		ylabel(, labsize(large)) ///
		xlabel(, labsize(large))
	graph export "figures/gov_exp_vs_gdp_2025.eps", replace

	* 2b. vs Palma ratio
	twoway scatter GEEP_pc palma, ///
		ytitle("Gen. Gov. EEP per capita (EUR)", size(large)) ///
		xtitle("Palma ratio", size(large)) ///
		ylabel(, labsize(large)) ///
		xlabel(, labsize(large))
	graph export "figures/gov_exp_vs_palma_2025.eps", replace

	* 2c. vs Gini index
	twoway scatter GEEP_pc gini, ///
		ytitle("Gen. Gov. EEP per capita (EUR)", size(large)) ///
		xtitle("Gini index", size(large))  ///
		ylabel(, labsize(large)) ///
		xlabel(, labsize(large))
	graph export "figures/gov_exp_vs_gini_2025.eps", replace

}


********************************************************************************
********************************** Palma Ratio *********************************
********************************************************************************

*------------------------------------------------------------------*
**# Bookmark #4 - Main analysis with Palma ratio - Section 3
*------------------------------------------------------------------*

{ // PALMA RATIO - Main analysis

    { // Tests and final model: Within Fixed Effects (WFE)

        *---------------------------------------------------*
        * 1. Random vs Fixed effects: Hausman test
        *---------------------------------------------------*
        xtreg log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic i.time, fe
        estimates store fixed

        xtreg log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic i.time, re
        estimates store random

        hausman fixed random, sigmamore

        /* Results interpretation:

        H0: No systematic difference in coefficients between FE and RE
        chi2(17) = 74.13, Prob > chi2 = 0.0000

        → Reject H0 at 1% significance level.
        → There is significant correlation between unobserved individual effects (country-specific) 
          and explanatory variables.
        → Fixed effects model is preferred.
        */


        *---------------------------------------------------*
        * 2. Time fixed effects test
        *---------------------------------------------------*
        xi: xtreg log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur share_pub_exp ///
                 pop_growth_1 rec snm ghp log_organic i.time, fe 
        testparm _Itime*

        /* Results interpretation:

        H0: All year dummies coefficients = 0 → no global year effect
        H1: At least one year effect is significant

        F(20,313) = 2.62, Prob > F = 0.0002
        → Reject H0: year fixed effects are collectively significant
        */


        *---------------------------------------------------*
        * 3. Heteroskedasticity test
        *---------------------------------------------------*
        xtreg log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic i.time, fe
        xttest3

        /* Results interpretation:
		
        Modified Wald test for groupwise heteroskedasticity
        H0: σi^2 = σ^2 for all i → homoskedasticity
        Prob>chi2 = 0.0000 → Reject H0
        → Significant groupwise heteroskedasticity. Standard errors need correction.
        */


        *---------------------------------------------------*
        * 4. Autocorrelation test (Wooldridge)
        *---------------------------------------------------*
        xtserial log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur ///
                 share_pub_exp pop_growth_1 rec snm ghp log_organic

        /* Results interpretation:
        H0: No first-order autocorrelation
        Prob > F = 0.0000 → Reject H0
        → AR(1) autocorrelation detected in residuals.
        */


        *---------------------------------------------------*
        * 5. Overidentification test (FE vs RE with clustered SE)
        *---------------------------------------------------*
        xtreg log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic, vce(cluster geo1) fe
        estimates store fixed

        xtreg log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic, vce(cluster geo1) re
        estimates store random

        xtoverid

        /* Results interpretation:
        Sargan-Hansen statistic  83.095, P-value = 0.0000
        → Reject H0: Fixed effects model preferred
        */


        *---------------------------------------------------*
        * 6. Cross-sectional dependence test (Pesaran)
        *---------------------------------------------------*
        xi: xtreg log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur share_pub_exp ///
                 pop_growth_1 rec snm ghp log_organic i.time, vce(cluster geo1) fe 
        xtcsd, pesaran abs

        /* Results interpretation:
        Pesaran's test: statistic = -2.790, Pr = 0.0053
        → Reject H0: significant cross-sectional dependence
        */


        *---------------------------------------------------*
        * 7. Final model estimation
        *---------------------------------------------------*
        xtscc log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 i.time, fe
        estimates store PALMA_WFEP

        xtscc log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic i.time, fe
        estimates store PALMA_WFEC


        *---------------------------------------------------*
        * 8. U-shape test
        *---------------------------------------------------*
        utest L5_palma L5_palma_2

        /* Results interpretation:
        Extreme point: 1.4336
        Interval around extreme:
            Lower bound: 0.744
            Upper bound: 1.749
        Slope left:  -0.757, p=0.003
        Slope right: 0.346, p=0.029
        → Supports existence of U-shape, but few observations in upper branch
        */

    } // end Tests and final model: WFE

	
	{ // Tests and final model: Grouped (Within-) Fixed Effects (GFE) 

        *---------------------------------------------------*
        * 1. GFE model with controls only
        *---------------------------------------------------*
        xtscc log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur ///
              share_pub_exp pop_growth_1 ///
              i.gfe_group#i.time, fe
        estimates store PALMA_GFEP

        *---------------------------------------------------*
        * 2. GFE model with controls + environmental variables
        *---------------------------------------------------*
        xtscc log_GEEP_pc L5_palma L5_palma_2 log_gdp_pcap_eur ///
              share_pub_exp pop_growth_1 rec snm ghp log_organic ///
              i.gfe_group#i.time, fe
        estimates store PALMA_GFEC

        *---------------------------------------------------*
        * 3. U-shape test
        *---------------------------------------------------*
        utest L5_palma L5_palma_2

        /* Results interpretation:

        Extreme point: 1.5299

        Slope to the left of the extreme: -1.0944, p = 0.0024 → significant decrease of dependent variable as L5_palma increases
        Slope to the right of the extreme: 0.3044, p = 0.0821 → increase not statistically significant

        Overall U-shape test:
        t-value = 1.44, p-value = 0.0821

        → Not significant at 5%, marginally significant at 10%
        → Cannot confidently conclude presence of a U-shape
        */

    } // end Tests and final model: GFE

	
    { // Tests and final model: System Generalized method of moments (GMM)

        *---------------------------------------------------*
        * 1. Unit root tests (Fisher-type) for Palma ratio and controls
        *---------------------------------------------------*
        xtunitroot fisher palma, pperron lags(1) demean
        xtunitroot fisher gdp_pcap, pperron lags(1) demean
        xtunitroot fisher share_pub_exp, pperron lags(1) demean
        xtunitroot fisher pop_growth_1, pperron lags(1) demean
        xtunitroot fisher rec, pperron lags(1) demean
        xtunitroot fisher snm, pperron lags(1) demean
        xtunitroot fisher ghp, pperron lags(1) demean
        xtunitroot fisher log_organic, pperron lags(1) demean

        /* Results interpretation:
        H0: All panels contain unit roots
        Ha: At least one panel is stationary
        → Reject H0 for all variables, series are stationary after first difference or demeaned
        */


        *---------------------------------------------------*
        * 2. GMM estimation
        *---------------------------------------------------*
        xtabond2 log_GEEP_pc L.log_GEEP_pc L5_palma L5_palma_2 ///
                 log_gdp_pcap_eur share_pub_exp pop_growth_1 time, ///
                 gmm(L.log_GEEP_pc, lag(2 3) collapse) ///
                 gmm(L5_palma L5_palma_2, lag(6 6) collapse) ///
                 gmm(log_gdp_pcap_eur, lag(1 1) collapse) ///
                 iv(share_pub_exp pop_growth_1, eq(level)) ///
                 twostep robust small
        estimates store PALMA_GMMP

        xtabond2 log_GEEP_pc L.log_GEEP_pc L5_palma L5_palma_2 ///
                 log_gdp_pcap_eur share_pub_exp pop_growth_1 rec snm ghp ///
                 log_organic time, ///
                 gmm(L.log_GEEP_pc, lag(2 3) collapse) ///
                 gmm(L5_palma L5_palma_2, lag(6 6) collapse) ///
                 gmm(log_gdp_pcap_eur log_organic, lag(1 1) collapse) ///
                 iv(share_pub_exp pop_growth_1 rec snm ghp time, eq(level)) ///
                 twostep robust small
        estimates store PALMA_GMMC


        *---------------------------------------------------*
        * 3. U-shape test
        *---------------------------------------------------*
        utest L5_palma L5_palma_2

        /* Results interpretation:
        Extreme point: 1.2455
        Left slope:  -2.371, p=0.0071
        Right slope: +2.379, p=0.0075
        Overall U-shape test: t=2.62, p=0.0075

        → Strong evidence of a significant U-shaped relationship.
        → Minimum at L5_palma ≈ 1.25
        → Left of minimum: L5_palma increase → dependent variable decreases
        → Right of minimum: L5_palma increase → dependent variable increases
        */

    } // end Tests and final model

	
	{ // Tables included in the article  

	
		/* ---------------Note on Table 3 – Within-FE specification tests:------------------
		
		   This table summarizes the main specification tests for the Within-FE (WFE) models,
		   including:
			 - Hausman test (Random vs Fixed effects),
			 - Time-fixed effects test,
			 - Modified Wald test for groupwise heteroskedasticity,
			 - Wooldridge test for AR(1) autocorrelation,
			 - Sargan-Hansen overidentification test (cluster-robust),
			 - Pesaran's cross-sectional dependence (CD) test.

		   The results in Table 3 are directly based on the outputs of these commands
		   already computed in the WFE section.
		----------------------------------------------------------------------------------*/

		/* ---------------Note on Table 4 – System-GMM specification tests:------------------
		   

		   This table summarizes the main specification tests for the System-GMM models,
		   including:
			 - Arellano–Bond test for AR(1) and AR(2) autocorrelation,
			 - Sargan test of overidentifying restrictions,
			 - Hansen J-test (robust),
			 - Difference-in-Hansen tests for instrument exogeneity.

		   The results in Table 4 are directly based on the outputs of these commands
		   already computed in the GMM section.
		----------------------------------------------------------------------------------*/

		/* -----Note on Table 5: Table 5 General government EEP per capita (LaTeX Output)-----
		
		This command exports the regression results from multiple models
		(Within-FE, Grouped-FE, System-GMM) to a LaTeX table using esttab.
		
		Models included:
		  - PALMA_WFEP : Within-FE, partial specification
		  - PALMA_WFEC : Within-FE, full specification
		  - PALMA_GFEP : Grouped-FE, partial specification
		  - PALMA_GFEC : Grouped-FE, full specification
		  - PALMA_GMMP : System-GMM, partial specification
		  - PALMA_GMMC : System-GMM, full specification
		
		Notes on table formatting:
		  - Stars indicate significance: * p<0.10, ** p<0.05, *** p<0.01
		  - Dropped variables: all time dummies (tim*) are excluded
		  - Dependent variable, GDP per capita, and organic consumption per capita are in log
		  - Driscoll-Kraay robust standard errors reported in parentheses for column (1)
		  - Additional stats displayed: N (observations), r2_w (within R-squared)
		----------------------------------------------------------------------------------*/

		esttab PALMA_WFEP PALMA_WFEC PALMA_GFEP PALMA_GFEC PALMA_GMMP PALMA_GMMC ///
			using "tables/PALMA.tex", ///
			star(* 0.1 ** 0.05 *** 0.01) ///
			title("General government EEP per capita") ///
			drop(tim*) ///
			nodepvars ///
			addnotes("Driscoll-Kraay standard errors are in parentheses in column (1). The dependent variable, GDP per capita, and organic consumption per capita are expressed in log.") ///
			label ///
			stats(N r2_w) ///
			tex replace
			
			
		/* --------------------------------Note on Table 6:--------------------------------
		   
		   Table 6 in the article summarizes the results of the U-shape tests (utest) 
		   that were conducted for the three model specifications: Within-FE (WFE), 
		   Grouped-FE (GFE), and System-GMM. 

		   The table reports:
			 - The estimated extreme point (vertex of the parabola),
			 - The slope, t-value, and p-value on either side of the extreme point,
			 - The number of observations below and above the thresholds, 
			   corresponding to the lower and upper branches of the curve.

		   No additional code is required to replicate Table 6; it is based on the utest 
		   outputs already computed in the preceding sections of the code.
		----------------------------------------------------------------------------------*/

	}

	
	{ // Graphical representation: marginal effects + 90% CI - Figure 3

		/********************************************************************/
		/* 1.  Compute marginal effects + 90% CI for GFE, WFE, GMM         */
		/********************************************************************/

		// ---------- GFE ----------
		preserve
		estimates restore PALMA_GFEC
		clear
		set obs 30
		range palma_seq 0.7 1.8 30

		matrix b = e(b)
		matrix V = e(V)

		scalar b_g  = b[1,"L5_palma"]
		scalar b_g2 = b[1,"L5_palma_2"]
		scalar v_g   = V["L5_palma","L5_palma"]
		scalar v_g2  = V["L5_palma_2","L5_palma_2"]
		scalar cov12 = V["L5_palma","L5_palma_2"]

		gen double effect_gfe = b_g + 2*b_g2*palma_seq
		gen double var_gfe    = v_g + 4*palma_seq^2*v_g2 + 4*palma_seq*cov12
		gen double se_gfe     = sqrt(var_gfe)
		gen double lb_gfe     = effect_gfe - 1.645*se_gfe
		gen double ub_gfe     = effect_gfe + 1.645*se_gfe

		tempfile gfe_marg
		save `gfe_marg', replace
		restore

		// ---------- WFE ----------
		preserve
		estimates restore PALMA_WFEC
		clear
		set obs 30
		range palma_seq 0.7 1.8 30

		matrix b = e(b)
		matrix V = e(V)

		scalar b_g  = b[1,"L5_palma"]
		scalar b_g2 = b[1,"L5_palma_2"]
		scalar v_g   = V["L5_palma","L5_palma"]
		scalar v_g2  = V["L5_palma_2","L5_palma_2"]
		scalar cov12 = V["L5_palma","L5_palma_2"]

		gen double effect_wfe = b_g + 2*b_g2*palma_seq
		gen double var_wfe    = v_g + 4*palma_seq^2*v_g2 + 4*palma_seq*cov12
		gen double se_wfe     = sqrt(var_wfe)
		gen double lb_wfe     = effect_wfe - 1.645*se_wfe
		gen double ub_wfe     = effect_wfe + 1.645*se_wfe

		tempfile wfe_marg
		save `wfe_marg', replace
		restore

		// ---------- GMM ----------
		preserve
		estimates restore PALMA_GMMC
		clear
		set obs 30
		range palma_seq 0.7 1.8 30

		matrix b = e(b)
		matrix V = e(V)

		scalar b_g  = b[1,"L5_palma"]
		scalar b_g2 = b[1,"L5_palma_2"]
		scalar v_g   = V["L5_palma","L5_palma"]
		scalar v_g2  = V["L5_palma_2","L5_palma_2"]
		scalar cov12 = V["L5_palma","L5_palma_2"]

		gen double effect_gmm = b_g + 2*b_g2*palma_seq
		gen double var_gmm    = v_g + 4*palma_seq^2*v_g2 + 4*palma_seq*cov12
		gen double se_gmm     = sqrt(var_gmm)
		gen double lb_gmm     = effect_gmm - 1.645*se_gmm
		gen double ub_gmm     = effect_gmm + 1.645*se_gmm

		tempfile gmm_marg
		save `gmm_marg', replace
		restore


		/********************************************************************/
		/* 2.  Plot marginal effects separately                              */
		/********************************************************************/

		// -------- GMM --------
		use `gmm_marg', clear
		twoway ///
			(rcap lb_gmm ub_gmm palma_seq, lcolor(green) lwidth(medthick)) ///
			(line effect_gmm palma_seq, lcolor(green) lpattern(solid) lwidth(medthick)) ///
			, ///
			yline(0, lpattern(dash)) ///
			yscale(range(-4 4)) ///
			ylabel(-4(1)4, angle(0) labsize(small)) ///
			xscale(range(0.7 1.8)) ///
			xlabel(0.7(0.1)1.8, labsize(small) angle(0)) ///
			ytitle("Marginal effect") ///
			xtitle("Palma Ratio") ///
			title("Sys-GMM - Marginal effects (90% CI)") ///
			legend(off) ///
			xsize(6) ysize(3)
		graph save figures/marginal_effect_PALMA_GMM_IC90bar.gph, replace

		// -------- GFE --------
		use `gfe_marg', clear
		twoway ///
			(rcap lb_gfe ub_gfe palma_seq, lcolor(blue) lwidth(medthick)) ///
			(line effect_gfe palma_seq, lcolor(blue) lpattern(solid) lwidth(medthick)) ///
			, ///
			yline(0, lpattern(dash)) ///
			yscale(range(-4 4)) ///
			ylabel(-4(1)4, angle(0) labsize(small)) ///
			xscale(range(0.7 1.8)) ///
			xlabel(0.7(0.1)1.8, labsize(small) angle(0)) ///
			ytitle("Marginal effect") ///
			xtitle("Palma Ratio") ///
			title("Grouped FE - Marginal effects (90% CI)") ///
			legend(off) ///
			xsize(6) ysize(3)
		graph save figures/marginal_effect_PALMA_GFE_IC90bar.gph, replace

		// -------- WFE --------
		use `wfe_marg', clear
		twoway ///
			(rcap lb_wfe ub_wfe palma_seq, lcolor(red) lwidth(medthick)) ///
			(line effect_wfe palma_seq, lcolor(red) lpattern(solid) lwidth(medthick)) ///
			, ///
			yline(0, lpattern(dash)) ///
			yscale(range(-4 4)) ///
			ylabel(-4(1)4, angle(0) labsize(small)) ///
			xscale(range(.7 1.8)) ///
			xlabel(0.7(0.1)1.8, labsize(small) angle(0)) ///
			ytitle("Marginal effect") ///
			xtitle("Palma Ratio") ///
			title("Within-DK - Marginal effects (90% CI)") ///
			legend(off) ///
			xsize(6) ysize(3)
		graph save figures/marginal_effect_PALMA_WFE_IC90bar.gph, replace


		/********************************************************************/
		/* 3.  Distribution of Palma ratio                                   */
		/********************************************************************/

		preserve
		use "stata_files/database_for_analysis.dta", clear
		keep if !missing(L5_palma, palma)

		kdensity L5_palma, n(100) gen(kx kd)

		twoway area kd kx if kx>=0.5 & kx<=2, ///
			color(gs12%40) ///
			plotregion(margin(l=8 r=2 b=4 t=4)) ///
			ytitle("Density") ///
			xtitle("Palma Ratio") ///
			title("Distribution - Palma ratio") ///
			xscale(range(.7 1.8)) ///
			xlabel(0.7(0.1)1.8, labsize(small) angle(0)) ///
			xsize(6) ysize(2.5)
		graph save figures/density_PALMA_IC90bar.gph, replace
		restore


		/********************************************************************/
		/* 4.  Final combined figure                                        */
		/********************************************************************/

		graph combine ///
			"figures/density_PALMA_IC90bar.gph" ///
			"figures/marginal_effect_PALMA_WFE_IC90bar.gph" ///
			"figures/marginal_effect_PALMA_GFE_IC90bar.gph" ///
			"figures/marginal_effect_PALMA_GMM_IC90bar.gph", ///
			rows(4) xsize(6) ysize(12) imargin(0 0 0 0) ///
			saving("figures/combo_PALMA_IC90bar_separate.gph", replace)

		graph export "figures/combo_PALMA_IC90bar_separate.eps", replace

	}


}

*------------------------------------------------------------------*
**# Bookmark #5 - Appendix B.2 - Panel VAR & Granger Causality
*------------------------------------------------------------------*

{ // PALMA RATIO - Panel VAR & Granger Causality

	* Load dataset if necessary
	use "stata_files/database_for_analysis.dta", clear

	* Define panel structure
	xtset geo1 time 

	* Create shorter variable names to avoid long names
	gen gdppcap = log_gdp_pcap_eur  

	*---------------------------------------------------------------*
	* Panel Vector Autoregression (PVAR) estimation
	* Corresponds to Table 8 in the article
	* Endogenous variables: log_GEEP_pc, palma, palma_2, gdppcap, log_organic
	* Exogenous variables: share_pub_exp, pop_growth_1, rec, snm, ghp
	* GMM instruments: lags 1 to 5 of the endogenous variables
	*---------------------------------------------------------------*
	pvar log_GEEP_pc palma palma_2 gdppcap log_organic, ///
		 lags(1) ///
		 exog(share_pub_exp pop_growth_1 rec snm ghp) ///
		 gmmstyle instlags(1/5)

	*---------------------------------------------------------------*
	* Granger causality Wald test
	* Corresponds to Table 9 in the article
	* Tests whether each excluded variable Granger-causes each equation
	*---------------------------------------------------------------*
	pvargranger

	*---------------------------------------------------------------*
	* Notes on results
	* - Table 8: PVAR coefficients with standard errors in parentheses
	*   Endogenous variables with lag subscript _{t-1} correspond to first lags
	*   Exogenous variables: Pub.Exp/GDP, Δpop, REC, SNM, GHP
	*   log_GEEP_pc, gdppcap, and log_organic are in logarithms
	*   Model estimated with GMM using lags 1-5 as instruments
	*
	* - Table 9: Granger causality Wald test
	*   H0: Excluded variable does not Granger-cause the dependent variable
	*   H1: Excluded variable Granger-causes the dependent variable
	*------------------------------------------------------------------*

}


********************************************************************************
********************************** Gini Index **********************************
********************************************************************************

use "stata_files/database_for_analysis.dta", clear

*------------------------------------------------------------------*
**# Bookmark #6 - Appendix B.3 - Gini
*------------------------------------------------------------------*

{ // Analysis with Gini index - Appendix B.3

	{ // Tests and final model: Within Fixed Effects (WFE)

        *---------------------------------------------------*
        * 1. Random vs Fixed effects: Hausman test
        *---------------------------------------------------*
        xtreg log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic i.time, fe
        estimates store fixed

        xtreg log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic i.time, re
        estimates store random

        hausman fixed random, sigmamore

        /* Results interpretation:

        H0: No systematic difference in coefficients between FE and RE
        chi2(17) = 207.34, Prob > chi2 = 0.0000

        → Reject H0 at 1% significance level.
        → There is significant correlation between unobserved individual effects (country-specific) 
          and explanatory variables.
        → Fixed effects model is preferred.
        */


        *---------------------------------------------------*
        * 2. Time fixed effects test
        *---------------------------------------------------*
        xi: xtreg log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur share_pub_exp ///
                 pop_growth_1 rec snm ghp log_organic i.time, fe 
        testparm _Itime*

        /* Results interpretation:

        H0: All year dummies coefficients = 0 → no global year effect
        H1: At least one year effect is significant

        F(20,313) = 1.97, Prob > F = 0.0072
        → Reject H0: year fixed effects are collectively significant
        */


        *---------------------------------------------------*
        * 3. Heteroskedasticity test
        *---------------------------------------------------*
        xtreg log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic i.time, fe
        xttest3

        /* Results interpretation:
		
        Modified Wald test for groupwise heteroskedasticity
        H0: σi^2 = σ^2 for all i → homoskedasticity
        Prob>chi2 = 0.0000 → Reject H0
        → Significant groupwise heteroskedasticity. Standard errors need correction.
        */


        *---------------------------------------------------*
        * 4. Autocorrelation test (Wooldridge)
        *---------------------------------------------------*
        xtserial log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur ///
                 share_pub_exp pop_growth_1 rec snm ghp log_organic

        /* Results interpretation:
        H0: No first-order autocorrelation
        Prob > F = 0.0000 → Reject H0
        → AR(1) autocorrelation detected in residuals.
        */


        *---------------------------------------------------*
        * 5. Overidentification test (FE vs RE with clustered SE)
        *---------------------------------------------------*
        xtreg log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic, vce(cluster geo1) fe
        estimates store fixed

        xtreg log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic, vce(cluster geo1) re
        estimates store random

        xtoverid

        /* Results interpretation:
        Sargan-Hansen statistic  86.406, P-value = 0.0000
        → Reject H0: Fixed effects model preferred
        */


        *---------------------------------------------------*
        * 6. Cross-sectional dependence test (Pesaran)
        *---------------------------------------------------*
        xi: xtreg log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur share_pub_exp ///
                 pop_growth_1 rec snm ghp log_organic i.time, vce(cluster geo1) fe 
        xtcsd, pesaran abs

        /* Results interpretation:
        Pesaran's test: statistic = -2.428, Pr = 0.0152
        → Reject H0: significant cross-sectional dependence
        */


        *---------------------------------------------------*
        * 7. Final model estimation
        *---------------------------------------------------*
        xtscc log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 i.time, fe
        estimates store GINI_WFEP

        xtscc log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur share_pub_exp ///
              pop_growth_1 rec snm ghp log_organic i.time, fe
        estimates store GINI_WFEC


        *---------------------------------------------------*
        * 8. U-shape test
        *---------------------------------------------------*
        utest L5_gini L5_gini_2

        /* Results interpretation:
        Extreme point: 30.72743
        Interval around extreme:
            Lower bound: 20
            Upper bound: 38.9
        Slope left:  -.0728748, p=0.0000613
        Slope right: .0555188, p=0.002228
        → Supports existence of U-shape, but few observations in upper branch
        */

    } // end Tests and final model: WFE

	
	{ // Tests and final model: Grouped (Within-) Fixed Effects (GFE) 

        *---------------------------------------------------*
        * 1. GFE model with controls only
        *---------------------------------------------------*
        xtscc log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur ///
              share_pub_exp pop_growth_1 ///
              i.gfe_group#i.time, fe
        estimates store GINI_GFEP

        *---------------------------------------------------*
        * 2. GFE model with controls + environmental variables
        *---------------------------------------------------*
        xtscc log_GEEP_pc L5_gini L5_gini_2 log_gdp_pcap_eur ///
              share_pub_exp pop_growth_1 rec snm ghp log_organic ///
              i.gfe_group#i.time, fe
        estimates store GINI_GFEC

        *---------------------------------------------------*
        * 3. U-shape test
        *---------------------------------------------------*
        utest L5_gini L5_gini_2

        /* Results interpretation:

        Extreme point: 32.81234

        Slope to the left of the extreme: -.0841621, p = 0.00005 → significant decrease of dependent variable as L5_gini increases
        Slope to the right of the extreme: .0399888, p = 0.0066848 → significant increase of dependent variable as L5_gini increases

        Overall U-shape test:
        t-value = 2.70, p-value = 0.00668

        → Supports existence of U-shape, but few observations in upper branch
        */

    } // end Tests and final model: GFE

	
    { // Tests and final model: System Generalized method of moments (GMM)

        *---------------------------------------------------*
        * 1. Unit root tests (Fisher-type) for Palma ratio and controls
        *---------------------------------------------------*
        xtunitroot fisher gini, pperron lags(1) demean
        xtunitroot fisher gdp_pcap, pperron lags(1) demean
        xtunitroot fisher share_pub_exp, pperron lags(1) demean
        xtunitroot fisher pop_growth_1, pperron lags(1) demean
        xtunitroot fisher rec, pperron lags(1) demean
        xtunitroot fisher snm, pperron lags(1) demean
        xtunitroot fisher ghp, pperron lags(1) demean
        xtunitroot fisher log_organic, pperron lags(1) demean

        /* Results interpretation:
        H0: All panels contain unit roots
        Ha: At least one panel is stationary
        → Reject H0 for all variables, series are stationary after first difference or demeaned
        */


        *---------------------------------------------------*
        * 2. GMM estimation
        *---------------------------------------------------*
        xtabond2 log_GEEP_pc L.log_GEEP_pc L5_gini L5_gini_2 ///
                 log_gdp_pcap_eur share_pub_exp pop_growth_1 time, ///
                 gmm(L.log_GEEP_pc, lag(2 2) collapse) ///
				 gmm(L5_gini L5_gini_2, lag(6 6) collapse) ///
				 gmm(log_gdp_pcap_eur, lag(1 1) collapse) ///
				 iv(share_pub_exp pop_growth_1 time, eq(level)) ///
				 twostep robust small
        estimates store GINI_GMMP

        xtabond2 log_GEEP_pc L.log_GEEP_pc L5_gini L5_gini_2 ///
                 log_gdp_pcap_eur share_pub_exp pop_growth_1 ///
                 rec snm ghp log_organic time, ///
				 gmm(L.log_GEEP_pc, lag(2 2) collapse) ///
				 gmm(L5_gini L5_gini_2, lag(6 6) collapse) ///
				 gmm(log_gdp_pcap_eur, lag(1 1) collapse) ///
				 iv(log_organic share_pub_exp pop_growth_1 rec snm ghp time, eq(level)) ///
				 twostep robust small
        estimates store GINI_GMMC


        *---------------------------------------------------*
        * 3. U-shape test
        *---------------------------------------------------*
        utest L5_gini L5_gini_2

		/* Results interpretation:
		   Extreme point: 28.937
		   Left slope (below extreme point):  -0.196, p=0.0565
		   Right slope (above extreme point): +0.219, p=0.0190
		   Overall U-shape test: t=1.645, p=0.0565

		   → Evidence of a potential U-shaped relationship, but significance is marginal.
		   → Minimum of the parabola at L5_gini ≈ 28.94
		   → Left of minimum: L5_gini increase → dependent variable tends to decrease (weak significance)
		   → Right of minimum: L5_gini increase → dependent variable increases (significant)
		   → Overall U-shape test is marginally non-significant at 5% level, but close to significance.
		*/

    } // end Tests and final model

	
	{ // Tables included in the article  

	
		/* ---------------Note on Table 10 – Within-FE specification tests:------------------
		
		   This table summarizes the main specification tests for the Within-FE (WFE) models,
		   including:
			 - Hausman test (Random vs Fixed effects),
			 - Time-fixed effects test,
			 - Modified Wald test for groupwise heteroskedasticity,
			 - Wooldridge test for AR(1) autocorrelation,
			 - Sargan-Hansen overidentification test (cluster-robust),
			 - Pesaran's cross-sectional dependence (CD) test.

		   The results in Table 10 are directly based on the outputs of these commands
		   already computed in the WFE section.
		----------------------------------------------------------------------------------*/

		/* ---------------Note on Table 11 – System-GMM specification tests:------------------
		   

		   This table summarizes the main specification tests for the System-GMM models,
		   including:
			 - Arellano–Bond test for AR(1) and AR(2) autocorrelation,
			 - Sargan test of overidentifying restrictions,
			 - Hansen J-test (robust),
			 - Difference-in-Hansen tests for instrument exogeneity.

		   The results in Table 11 are directly based on the outputs of these commands
		   already computed in the GMM section.
		----------------------------------------------------------------------------------*/

		/* ---------Note on Table 12: General government EEP per capita (LaTeX Output)--------
		
		This command exports the regression results from multiple models
		(Within-FE, Grouped-FE, System-GMM) to a LaTeX table using esttab.
		
		Models included:
		  - GINI_WFEP : Within-FE, partial specification
		  - GINI_WFEC : Within-FE, full specification
		  - GINI_GFEP : Grouped-FE, partial specification
		  - GINI_GFEC : Grouped-FE, full specification
		  - GINI_GMMP : System-GMM, partial specification
		  - GINI_GMMC : System-GMM, full specification
		
		Notes on table formatting:
		  - Stars indicate significance: * p<0.10, ** p<0.05, *** p<0.01
		  - Dropped variables: all time dummies (tim*) are excluded
		  - Dependent variable, GDP per capita, and organic consumption per capita are in log
		  - Driscoll-Kraay robust standard errors reported in parentheses for column (1)
		  - Additional stats displayed: N (observations), r2_w (within R-squared)
		----------------------------------------------------------------------------------*/

		esttab GINI_WFEP GINI_WFEC GINI_GFEP GINI_GFEC GINI_GMMP GINI_GMMC ///
			using "tables/GINI.tex", ///
			star(* 0.1 ** 0.05 *** 0.01) ///
			title("General government EEP per capita") ///
			drop(tim*) ///
			nodepvars ///
			addnotes("Driscoll-Kraay standard errors are in parentheses in column (1). The dependent variable, GDP per capita, and organic consumption per capita are expressed in log.") ///
			label ///
			stats(N r2_w) ///
			tex replace
			

	}

	
}



**# Bookmark # FIN
