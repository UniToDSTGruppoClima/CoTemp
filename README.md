# Co.Temp
## Comparing series of Temperatures
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1472761.svg)](https://doi.org/10.5281/zenodo.1472761)

This program aims to offer useful information on temperature series (i.e.: identify the biases in the daily maximum and minimum temperature, compare manual and automatic weather stations, etc.) and works in three steps: statistical analysis and characterization of the daily series of maximum (Tx) and minimum temperature (Tn), computation of monthly-aggregated data for both Tx and Tn and comparison between temperature classes of events (like heat wave, cold wave and normal events).

###### 1] Statistical analysis
In this first step, for each raw series (in example, a set of two Tx series and another set of two Tn series that came from different instruments), a statistical analysis is carried out. Mean, median, first and third quartiles, minimum and maximum values are calculated and missing values are identified. A time series plot and a density plot are created for each pair of series. Always for each pair of series, the t test, the Kolmogorov-Smirnov test and the Wilcoxon's Rank Sum test are carried out. Furthermore, the Root Mean Squared Error (RMSE), the Mean Error (ME), and the correlation coefficient by Spearman's method are calculated for every pair of series.

###### 2] Computation of monthly data
In the second step, any values that are missing in one series are also set to be missing (as NA) in its counterpart. For each pair of series, the daily difference is calculated in addition to the monthly difference series and the monthly Percentage Relative Error (PRE). A values of PRE>0 shows an overestimation of the second monthly Tx|Tn series (called candidate) over the first one (reference), while PRE<0 highlights an underestimation of the candidate over the reference series. On the monthly PRE series is also computed the trend with its slope, to identify if the shift between the pair of series increases, decreases or is constant.

###### 3] Classification of events
In the third step, daily data is divided in 5 classes, that are extremely cold events, cold events, mean events, warm events and extremely warm events, using percentiles. The percentiles are estimated by combining the daily series pair (candidate and reference). 

When you start the program, (if interactive mode is available) you will be asked for an input text file and a folder where results will be stored. If you are running in batch mode, please modify the code to use the correct input file and create the right output folder.

## INPUT
The text file has to be formatted in five TAB-separated columns. The first row of the file has to contain the headers (column names) while the first column has to contain the dates (in DD/MM/YYYY format). Column two is the Tx reference serie and column three is the Tn reference serie. Column four and five should contain Tx and Tn candidate series. Missing values must be marked as NA. It is very important to start the file from the first of January (of any year) and end it at the 31st of December (of any year). See the attached file called example_.txt (where _ is a number).

## OUTPUT
- 01_TN_statistics_raw.csv - Main statistics on input Tn serie.
- 01_TN_summary_raw.csv - Summary of input Tn series and their difference.
- 01_TX_statistics_raw.csv - Main statistics on input Tx serie.
- 01_TX_summary_raw.csv - Summary of input Tx series and their difference.
- 01_TN_plot_raw.pdf - Plots of input Tn series.
- 01_TX_plot_raw.pdf - Plots of input Tx series.
- 02_TN_statistics_boxplot_monthly_diff.csv - Monthly differences of candidate and reference Tn.
- 02_TN_trend_relative error.csv - Informations on trend of Tn.
- 02_TN_boxplot_monthly_diff.pdf - Boxplot of monthly series, differences and PRE of Tn.
- 02_TN_plot_diff.pdf - Plot and distribution of difference in Tn series.
- 02_TN_relative_error.pdf - Plots of Tn PRE.
- 02_TX_statistics_boxplot_monthly_diff.csv - Monthly differences of candidate and reference Tx.
- 02_TX_trend_relative error.csv - Informations on trend of Tx.
- 02_TX_boxplot_monthly_diff.pdf - Boxplot of monthly series, differences and PRE of Tx.
- 02_TX_plot_diff.pdf - Plot and distribution of difference in Tx series.
- 02_TX_relative_error.pdf - Plots of Tx PRE.
- 03_TN_classes_candidate.csv - Statistics on the classes of candidate Tn.
- 03_TN_classes_reference.csv - Statistics on the classes of reference Tn.
- 03_TN_cold.csv - Statistics on Tn events classified as cold.
- 03_TN_extr_cold.csv - Statistics on Tn events classified as extremely cold.
- 03_TN_extr_warm.csv - Statistics on Tn events classified as extremely warm.
- 03_TN_frequencies_in_classes.csv - Frequencies of Tn events in every class.
- 03_TN_mean.csv - Statistics on Tn events classified as mean.
- 03_TN_warm.csv - Statistics on Tn events classified as warm.
- 03_TN_boxplot_diff.pdf - Boxplot of each class for Tn events.
- 03_TX_classes_candidate.csv - Statistics on the classes of candidate Tx.
- 03_TX_classes_reference.csv - Statistics on the classes of reference Tx.
- 03_TX_cold.csv - Statistics on Tx events classified as cold.
- 03_TX_extr_cold.csv - Statistics on Tx events classified as extremely cold.
- 03_TX_extr_warm.csv - Statistics on Tx events classified as extremely warm.
- 03_TX_frequencies_in_classes.csv - Frequencies of Tn events in every class.
- 03_TX_mean.csv - Statistics on Tx events classified as mean.
- 03_TX_warm.csv - Statistics on Tx events classified as warm.
- 03_TX_boxplot_diff.pdf - Boxplot of each class for Tx events.
- TN_Results.xlsx - Collection of previous results for Tn in a single Excel file.
- TX_Results.xlsx - Collection of previous results for Tx in a single Excel file.
