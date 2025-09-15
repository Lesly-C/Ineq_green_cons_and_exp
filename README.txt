README – Stata Code for Analysis
================================

Overview
--------
This repository contains the Stata code and data files used for the analysis presented in the article. 
The code was generated using Stata 18.5 with datasets downloaded in 2022 and 2023 from Eurostat, FiBL, and EPI. 

File structure
--------------
- ineq_green_cons/
    - dataset_creation.do   : Code to create the dataset from raw files.
    - data_analysis.do      : Main code performing all analyses, generating tables and figures.

- ineq_green_cons/initial_datasets/
    - Raw data files (Excel and CSV) used to build the database.

- ineq_green_cons/stata_files/
    - dataset.dta : Pre-built Stata dataset compatible with Stata 18 (recommended for direct use).

- ineq_green_cons/figures/ and ineq_green_cons/tables/
    - Output folders where figures and LaTeX tables will be saved by the code (initially empty).

Instructions
------------
1. Download all files from the repository, preserving the folder structure.

2. Create output folders if they do not exist:
       figures/
       tables/

3. (Optional) Create the dataset by running dataset_creation.do:
       - Open dataset_creation.do in the Do-file Editor.
       - Run dataset_creation.do to recreate the dataset from the raw files.

   This step is recommended if you use an older version of Stata.  
   Alternatively, you can use the pre-built dataset located in stata_files/ (Stata 18 version).

4. Run the main analysis:
       - Open data_analysis.do in the Do-file Editor.
       - Run data_analysis.do to generate the tables and figures.

5. Output:
   - Figures 1, 2, 3 are saved in ineq_green_cons/figures/
   - Tables 5 and 12 are exported directly as LaTeX files in ineq_green_cons/tables/
   - Other tables are created manually in LaTeX, with results collected from the Stata output:

     Palma ratio analysis
       - Table 1 : lines 213–214
       - Table 3 : lines 264–366 (different tests for the WFE analysis)
       - Table 4 : lines 461–478 (statistics from the Sys-GMM analysis)
       - Table 6 : post-estimation U-test on lines 382 (WFE), 418 (GFE), 484 (GMM)
       - Table 8 : Panel VAR analysis, lines 793–796
       - Table 9 : Granger causality Wald tests, line 803

     Gini analysis
       - Table 10 : lines 835–930 (different tests for the WFE analysis)
       - Table 11 : lines 1031–1048 (statistics from the Sys-GMM analysis)

Support
-------
If you encounter any issues or have questions regarding the code or data, 
please do not hesitate to contact Lesly Cassin: lesly.cassin@univ-lorraine.fr
