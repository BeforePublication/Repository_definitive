#Title: Monetary policy models: Lessons from the Eurozone crisis

#This folder contains all the data and computer procedures used in the manuscript "Monetary policy models: Lessons from the Eurozone crisis".

#Files in this folder are the following: 1. data_MPMLEC_HSSCOMMS.xlsx; 2. Eviews_procedures.wf1; 3. data_R_procedures.csv; 4. R_procedures.r


#The excel file "data_MPMLEC_HSSCOMMS.xlsx" contains and explains the dataset used in the manuscript. The different sheets provided in the Excel file "data_MPMLEC_HSSCOMMS.xlsx" are the following:

#Sheet "Raw_Baseline": contains the raw data used to build the VECM baseline model.
#Sheet "Source_Baseline": enumerates and describes all the data variables in "Raw_Baseline" and provides their URL links.
#Sheet "Raw_Extension": contains the raw data used to build the VECM extended model.
#Sheet "Source_Extension": enumerates and describes all the data variables in "Raw_Extension" and provides their URL links.
#Sheet "Model_Baseline": contains the data as used in the VECM baseline model.
#Sheet "Model_Extension": contains the data as used in the VECM extended model.
#Sheet "Variables_Notes": explains how series in "Model_Baseline"-"Model_Extended" are defined from series in "Raw_Baseline"-"Raw_Extension".
#Sheet "Figures_monthly": contains the data used for the figures with a monthly periodicity.
#Sheet "Figures_quarterly": contains the data used for the figures with a quarterly periodicity.
#Sheet "Source_Figures": enumerates and describes all the data variables and provides their URL links.

#The Eviews workfile "Eviews_procedures.wf1" contains the data series and procedures used to run the cointegration and VECM analyses. Data series in this workfile are those detailed in "data_MPMLEC_HSSCOMMS.xlsx".

#The csv file "data__R_procedures.csv" contains the dataset used in the R procedures. Data series in this workfile are those detailed in "data_MPMLEC_HSSCOMMS.xlsx". Excel data are converted to csv to improve manageability and sharing.

#The file "R_procedures.r" contains the R script to generate impulse-response figures and analyses.