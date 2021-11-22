# Data-Science-for-Biomedical-Informatics-Directory
BMIN 503

This project aims to conduct a large-scale proteomic analysis to screen potential plasma proteins that can A) differentiate PD from normal control B) differentiate 
fast cognitive progressors from slow cognitive progressors based on change in MoCA score change over time C) differentiate fast motor progressors from slow motor
progressors based on change in UPDRS score over time 

Background: Datasets used for this project were downloaded from NIH Parkinson's Disease Biomarker Program database. A total of 319 PD + Normal Control Subjects were enrolled in this study through two different cohorts - University of Texas Southwestern School of Medicine and Penn State Hershey School of Medicine. A total of ~1100 plasma proteins were measured in all 319 subjects through a high through-put proteomic platform known as SomaScan. Clinical and demographic information for these patients were also collected at each cohort. 

Analysis Plan:

A. Cross Section Analysis

1) Data Cleaning and Outlier Detection

Since those samples were collected from two different cohorts, it would be helpful to start with a principle component analysis (PCA) to detect any potential batch effect.
Besides, PCA analysis can help us quickly identify if there's any significant outliers so that we can remove outliers before formal analysis 

2) Differential Expression Analysis 

Compare the proteins expression level between normal control and subjects diagnosed with Parkinson's Disease. Benjamini-Hochberg method was used to adjust for 
multiple comparison. Proteins that pass the adjusted p-value of 0.05 were considered as proteins that were significantly associated with the phenotype of Parkinson's
Disease. Fold change was calculated using log2(Control-PD). Vocalno plot was used to visualize the difference in protein expression between PD and control group.

3). Gene Set Enrichment Analysis (GSEA)

After identifying plasma proteins that differentiate PD from control, we created a dictionary to map protein names into their corresponding gene symbols. Using KEGG
pathway files, GESA analysis was conducted to identify any potential underlying mechanisms associated with the proteins identified in the 2) and their associated
genes. 

4). Classifcation Analysis 

Using the top 3 proteins identified through the differential expression analysis, we next trained a logistic-regression based classification model, validated with 5-fold cross validation over 50 iterations, to understand how well the top 3 proteins we identified - BSP, OMD and ApoM - together with age and sex, can help differentiate PD subjects from control subjects.  

B) Longitudinal Analysis (Cognitive Function)

11/21/21


C) Longitudinal Analysis (Motor Function)

11/21/21
