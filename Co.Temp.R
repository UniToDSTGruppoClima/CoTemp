###############################################################################
############################################################################### 
#######################         Co.Temp         ###############################
################### Comparing series of Temperatures ##########################
###############################################################################
###############################################################################
###                                                                         ###
### Copyright (C) 2018-2019                                                 ###
### Fiorella Acquaotta, Diego Guenzi, Diego Garzena and Simona Fratianni    ###
###                                                                         ###
### This program is free software: you can redistribute it and/or modify    ###
### it under the terms of the GNU General Public License as published by    ###
### the Free Software Foundation, either version 3 of the License, or       ###
### (at your option) any later version.                                     ###
###                                                                         ###
### This program is distributed in the hope that it will be useful,         ###
### but WITHOUT ANY WARRANTY; without even the implied warranty of          ###
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ###
### GNU General Public License for more details.                            ###
###                                                                         ###
### You should have received a copy of the GNU General Public License       ###
### along with this program.  If not, see <http://www.gnu.org/licenses/>    ###
###                                                                         ###
###############################################################################
###############################################################################


###############################  Co.Temp  #####################################
# Description of the software
#
# This program aims to offer useful information on temperature series (i.e.:
# identify the biases in the daily maximum and minimum temperature, compare manual
# and automatic weather stations, etc.) and works in three steps: statistical
# analysis and characterization of the daily series of maximum (Tx) and minimum
# temperature (Tn), computation of monthly-aggregated data for both Tx and Tn and
# comparison between temperature classes of events (like heat wave, cold wave and
# normal events).
#
# 1] Statistical analysis
# In this first step, for each raw series (in example, a set of two Tx series and
# another set of two Tn series that came from different instruments), a statistical
# analysis is carried out. Mean, median, first and third quartiles, minimum and
# maximum values are calculated and missing values are identified. A time series
# plot and a density plot are created for each pair of series. Always for each pair
# of series, the t test, the Kolmogorov-Smirnov test and the Wilcoxon's Rank Sum
# test are carried out. Furthermore, the Root Mean Squared Error (RMSE) and the
# correlation coefficient by Spearman's method are calculated for every pair of
# series.
# 2] Computation of monthly data
# In the second step, any values that are missing in one series are also set to be
# missing (as NA) in its counterpart. For each pair of series, the daily difference
# is calculated in addition to the monthly difference series and the monthly
# Relative Ratio (RR). A values of RR>0 shows an overestimation of the
# second monthly Tx|Tn series (called candidate) over the first one (reference),
# while RR<1 highlights an underestimation of the candidate over the reference
# series. On the monthly RR series is also computed the trend with its slope, to
# identify if the shift between the pair of series increases, decreases or is
# constant. 
# 3] Classification of events
# In the third step, daily data is divided in 5 classes, that are extremely cold
# events, cold events, mean events, warm events and extremely warm events, using
# percentiles. The percentiles are estimated by combining the daily series pair
# (candidate and reference). 
#
# When you start the program, (if interactive mode is available) you will be 
# asked for an input text file and a folder where results will be stored. If
# you are running in batch mode, please modify the code to use the correct 
# input file and create the right output folder.
#
# INPUT
# The text file has to be formatted in five TAB-separated columns. The first
# row of the file has to contain the headers (column names) while the first
# column has to contain the dates (in DD/MM/YYYY format). Column two is the
# Tx reference series and column three is the Tn reference series. Column four
# and five should contain Tx and Tn candidate series. Missing values must be
# marked as NA. It is very important to start the file from the first of
# January (of any year) and end it at the 31st of December (of any year). See
# the attached file called example_.txt (where _ is a number).
#
# OUTPUT
# 01_TN_statistics_raw.csv - Main statistics on input Tn series
# 01_TN_summary_raw.csv - Summary of input Tn series and their difference
# 01_TX_statistics_raw.csv - Main statistics on input Tx series
# 01_TX_summary_raw.csv - Summary of input Tx series and their difference
# 01_TN_plot_raw.pdf - Plots of input Tn series
# 01_TX_plot_raw.pdf - Plots of input Tx series
# 02_TN_statistics_boxplot_monthly_diff.csv - Monthly differences of candidate and reference Tn
# 02_TN_trend_relative_ratio.csv - Information on trend of Tn
# 02_TN_boxplot_monthly_diff.pdf - Boxplot of monthly series, differences and RR of Tn
# 02_TN_plot_diff.pdf - Plot and distribution of difference in Tn series
# 02_TN_relative_ratio.pdf - Plots of Tn RR
# 02_TX_statistics_boxplot_monthly_diff.csv - Monthly differences of candidate and reference Tx
# 02_TX_trend_relative_ratio.csv - Information on trend of Tx
# 02_TX_boxplot_monthly_diff.pdf - Boxplot of monthly series, differences and RR of Tx
# 02_TX_plot_diff.pdf - Plot and distribution of difference in Tx series
# 02_TX_relative_ratio.pdf - Plots of Tx RR
# 03_TN_classes_candidate.csv - Statistics on the classes of candidate Tn
# 03_TN_classes_reference.csv - Statistics on the classes of reference Tn
# 03_TN_cold.csv - Statistics on Tn events classified as cold
# 03_TN_extr_cold.csv - Statistics on Tn events classified as extremely cold
# 03_TN_extr_warm.csv - Statistics on Tn events classified as extremely warm
# 03_TN_frequencies_in_classes.csv - Frequencies of Tn events in every class
# 03_TN_mean.csv - Statistics on Tn events classified as mean
# 03_TN_warm.csv - Statistics on Tn events classified as warm
# 03_TN_boxplot_diff.pdf - Boxplot of each class for Tn events
# 03_TX_classes_candidate.csv - Statistics on the classes of candidate Tx
# 03_TX_classes_reference.csv - Statistics on the classes of reference Tx
# 03_TX_cold.csv - Statistics on Tx events classified as cold
# 03_TX_extr_cold.csv - Statistics on Tx events classified as extremely cold
# 03_TX_extr_warm.csv - Statistics on Tx events classified as extremely warm
# 03_TX_frequencies_in_classes.csv - Frequencies of Tn events in every class
# 03_TX_mean.csv - Statistics on Tx events classified as mean
# 03_TX_warm.csv - Statistics on Tx events classified as warm
# 03_TX_boxplot_diff.pdf - Boxplot of each class for Tx events
# TN_Results.xlsx - Collection of previous results for Tn in a single Excel file
# TX_Results.xlsx - Collection of previous results for Tx in a single Excel file
#
# Following R-packages have to be installed:
#    "MASS", "timeDate", "timeSeries", "fBasics", "zoo", "hydroGOF", "xts",
#    "hydroTSM", "zyp", "xlsx"
#
# This code has been written under R-Version 3.3.0 and 3.4.3; for older or
# newer versions problems might occur.
#
###############################################################################
# Versions:
# 
# v1.0 - 20180319: First public release after code cleaning and review
# v1.1 - 20190628: Minor code modification and optimization, according to
#                  the reviewers of the methodological paper
# 
###############################################################################


###############################################################################
####                         INITIALIZATION                                ####
###############################################################################

start.time = Sys.time()
library(MASS)
library(timeDate)
library(timeSeries)
library(fBasics)
library(zoo)
library(hydroGOF)
library(xts)
library(hydroTSM)
library(zyp)
library(xlsx)

###############################################################################
####                       FUNCTION DEFINITION                             ####
###############################################################################

classification <- function(serie) {
  ##########################################################
  # Returns summary and length of the serie 
  #
  # INPUT
  # serie: a single temperature serie to analyze
  #
  # OUTPUT
  # Summary and lenght of the serie
  ##########################################################
  
  sum_serie = summary(serie[,4])
  len_serie = length(serie[,4])
  info_serie = c(sum_serie,len_serie)
  
  return(info_serie)
} # end of function classification


statistics <- function(serie1, serie2) {
  ##########################################################
  # Compute RMSE, cor and summary on common class events on both series
  #
  # INPUT
  # serie1: a single candidate temperature serie to analyze
  # serie2: a single reference temperature serie to analyze
  #
  # OUTPUT
  # Summary, lenght, RMSE and cor of the series
  ##########################################################
  
  rmse_common = rmse(serie1, serie2)
  cor_common = cor(serie1, serie2)
  sum_CAN = summary(serie1)
  sum_REF = summary(serie2)
  l_common = length(serie1)
  info_common = cbind(sum_CAN, sum_REF, l_common, rmse_common, cor_common)
  
  return(info_common)
} # end of function statistics


analyze <- function(name, Tserie) {
  ##########################################################
  # This function is the core of Co.Temp and does all the necessary
  # computations to produce the output previously described in the general
  # program description
  #
  # INPUT
  # name: the name of the serie going to be analyzed
  # Tserie: a single temperature serie to be analyzed
  #
  # OUTPUT
  # No output, except 1 xlsx and 12 csv files produces as side effects
  ##########################################################
  
  # Single series and differences
  can_serie = Tserie[,1:4] 
  ref_serie = Tserie[,c(-4,-6,-7)]
  dif_day = (can_serie[,4] - ref_serie[,4])
  Tserie = cbind(Tserie, dif_day)
  
  # Remove NAs
  Tserie1_na = cbind(can_serie, ref_serie[,4], dif_day)
  Tserie1_na = na.omit(Tserie1_na)
  Tserie_na = na.omit(Tserie)
  can_na = Tserie_na[,4]
  ref_na = Tserie_na[,5]
  
  # Create summary for original series, series without NAs and difference serie
  sum_can = summary(can_serie[,4])
  sum_ref = summary(ref_serie[,4])
  sum_can_na = summary(Tserie1_na[,4])
  sum_can_na[7] = 0
  sum_ref_na = summary(Tserie1_na[,5])
  sum_ref_na[7] = 0
  sum_dif_day = summary(dif_day)
  summ = c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","NA's")
  stat = suppressWarnings(cbind(summ, sum_can, sum_can_na, sum_ref, sum_ref_na, sum_dif_day))
  stat = as.data.frame(stat[,-1])
  for (i in c(1:5)) {
    stat[,i] = as.numeric(as.character(stat[,i]))
  }
  write.csv(stat, file=paste("01_",name,"_summary_raw.csv",sep=""), row.names=T)
  
  # Compute Zooreg
  raw_can = zooreg(can_serie[,4], start=as.Date(start.tab))
  raw_ref = zooreg(ref_serie[,4], start=as.Date(start.tab))
  raw_diff = zooreg(dif_day, start=as.Date(start.tab))
  
  # Save PDF plot
  pdf(paste("01_",name,"_plot_raw.pdf",sep=""))
  par(mfcol=c(2,2))
  plot(raw_can, ylab="Temperature (°C)", xlab="Years", 
       main=paste(name,"candidate station"), col="blue")
  plot(raw_ref, ylab="Temperature (°C)", xlab="Years", 
       main=paste(name,"reference station"), col="red")
  plot(raw_can, ylab="Temperature (°C)", xlab="Years", 
       main="Candidate (blue) and reference (red)", col="blue")
  lines(raw_ref, col="red")
  plot(density(ref_serie[,4],na=T), main="Candidate (blue) and reference (red)", col="red")
  lines(density(can_serie[,4],na=T), col="blue")
  dev.off()
  
  # Compute statistics: RMSE, T, KS, Wilcox and Cor between the two series
  rmse_er = round(rmse(ref_serie[,4],can_serie[,4],na.rm=T), digits=2)
  test_t = t.test(ref_serie[,4], can_serie[,4])
  test_ks = suppressWarnings(ks.test(ref_serie[,4], can_serie[,4]))
  test_wil = wilcox.test(ref_serie[,4], can_serie[,4])
  test_cor = cor.test(ref_serie[,4], can_serie[,4], method="spearman", exact=FALSE)
  test_names = c("RMSE", "T-Test p-value", "Kolmogorov-Smirnov p-value",
                 "Wilcoxon p-value", "Spearman rho", "Spearman p-value")
  test_res = c(rmse_er, test_t$p.value, test_ks$p.value, test_wil$p.value, 
               as.numeric(test_cor$estimate), test_cor$p.value)
  test_res = round(test_res, digits=2)
  test_fin = cbind(test_names, test_res)
  colnames(test_fin) = c("Test name", "Test result")
  test_fin = as.data.frame(test_fin)
  test_fin[,2] = as.numeric(as.character(test_fin[,2]))
  write.csv(test_fin, file=paste("01_",name,"_statistics_raw.csv",sep=""), row.names=F)
  
  # Compute monthly data
  dif_m = apply.monthly(raw_diff, mean, na.rm=T)
  Tserie_can_m = apply.monthly(raw_can, mean, na.rm=T)
  Tserie_ref_m = apply.monthly(raw_ref, mean, na.rm=T)
  dif_mm = matrix(dif_m, ncol=12, byrow=T)
  sum_dif_mon = summary(dif_mm)
  Tserie_can_mm = matrix(Tserie_can_m, ncol=12, byrow=T)
  Tserie_ref_mm = matrix(Tserie_ref_m, ncol=12, byrow=T)
  
  # Monthly relative ratio
  err_month = abs(Tserie_can_mm) / abs(Tserie_ref_mm)

  # Save Monthly boxplots PDF
  pdf(paste("02_",name,"_boxplot_monthly_diff.pdf",sep=""))
  par(mfcol=c(3,1))
  bp = boxplot(dif_mm, main="Monthly difference", ylab="Temperature difference (°C)", xlab="Months")
  boxplot(Tserie_can_mm, main=paste(name,"monthly candidate"), ylab="Temperature (°C)", xlab="Months")
  boxplot(Tserie_ref_mm, main=paste(name,"monthly reference"), ylab="Temperature (°C)", xlab="Months")
  dev.off()
  
  # Compute statistics on boxplots
  box_diff = round(bp$stats, digits=2)
  mm_dif = matrix(box_diff, nrow=5, ncol=12,
                  dimnames=list(c("Max","3rd Qu.","Median","1st Qu.","Min"),
                                c("J","F","M","A","M'","J'","J''","A'","S","O","N","D")))
  mm_dif = as.data.frame(mm_dif)
  write.csv(mm_dif, file=paste("02_",name,"_statistics_boxplot_monthly_diff.csv",sep=""), row.names=T)
  
  # Find values of Inf (where reference station has data equal to zero)
  pdf(paste("02_",name,"_relative_ratio.pdf",sep=""))
  par(mfcol=c(2,1))
  M = matrix(err_month, nrow=nrow(err_month), ncol=ncol(err_month), byrow=F)
  threshold = Inf
  M[M==threshold] = NA
  boxplot(M, main=paste(name,"relative ratio"), ylab="RR", xlab="Months")
  vect = as.vector(t(M))
  raw_e = zooreg(vect, frequency=12, start=(tab[1,6]))
  numb = c(1:length(vect))
  plot(numb, raw_e, type="l", main=paste(name,"relative ratio"), ylab="RR", xlab="Index")
  trend = zyp.trend.vector(vect, conf.intervals=T, preserve.range.for.sig.test=T)
  abline(trend[11], trend[2])
  dev.off()
  
  # Save trends
  trend = as.data.frame(trend)
  write.csv(trend, file=paste("02_",name,"_trend_relative_ratio.csv",sep=""), row.names=T)
  
  # Differences plot
  pdf(paste("02_",name,"_plot_diff.pdf",sep=""))
  par(mfcol=c(2,1))
  plot(raw_diff, type="l", main=paste(name,"difference (Candidate - Reference)"),
       ylab="Temperature difference (°C)", xlab="Years")
  abline(0, 0, col="red")
  hist(raw_diff, freq=T, main="History", xlab="Temperature difference (°C)")
  dev.off()
  
  # Classification
  data = c(ref_serie[,4], can_serie[,4])
  quant = quantile(data, c(.05,.20,.80,.95), na.rm=T)
  q05 = quant[1]
  q20 = quant[2]
  q80 = quant[3]
  q95 = quant[4]
  
  # Candidate Classes definition - Extremely cold
  extr_cold_can = Tserie_na[can_na<=q05,]
  Can_Extr_Cold = classification(extr_cold_can[,1:4])
  
  # Candidate Classes definition - Cold
  cold_can = Tserie_na[can_na>q05 & can_na<=q20,] 
  Can_Cold = classification(cold_can[,1:4])
  
  # Candidate Classes definition - Mean
  mean_can = Tserie_na[can_na>q20 & can_na<=q80,]
  Can_Mean = classification(mean_can[,1:4])
  
  # Candidate Classes definition - Warm
  warm_can = Tserie_na[can_na>q80 & can_na<=q95,]
  Can_Warm = classification(warm_can[,1:4])
  
  # Candidate Classes definition - Extremely Warm
  extr_warm_can = Tserie_na[can_na>q95,]
  Can_Extr_Warm = classification(extr_warm_can[,1:4])
  
  # Final result for candidate serie
  Candidate = c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","Length")
  mat_fin_can = cbind(Candidate, Can_Extr_Cold, Can_Cold, Can_Mean,Can_Warm, Can_Extr_Warm)
  mat_fin_can = as.data.frame(mat_fin_can)
  for (i in c(2:6)) {
    mat_fin_can[,i] = as.numeric(as.character(mat_fin_can[,i]))
  }
  write.csv(mat_fin_can, file=paste("03_",name,"_classes_candidate.csv",sep=""), row.names=F)
  
  # Reference Classes definition - Extremely cold
  extr_cold_ref = Tserie_na[ref_na<=q05,]
  extr_cold_ref = cbind(extr_cold_ref[,1:3], extr_cold_ref[5])
  Ref_Extr_Cold = classification(extr_cold_ref)
  
  # Reference Classes definition - Cold
  cold_ref = Tserie_na[ref_na>q05 & ref_na<=q20,]
  cold_ref = cbind(cold_ref[,1:3], cold_ref[5])
  Ref_Cold = classification(cold_ref)
  
  # Reference Classes definition - Mean
  mean_ref = Tserie_na[ref_na>q20 & ref_na<=q80,]
  mean_ref = cbind(mean_ref[,1:3], mean_ref[5])
  Ref_Mean = classification(mean_ref)
  
  # Reference Classes definition - Warm
  warm_ref = Tserie_na[ref_na>q80 & ref_na<=q95,]
  warm_ref = cbind(warm_ref[,1:3], warm_ref[5])
  Ref_Warm = classification(warm_ref)
  
  # Reference Classes definition - Extremely Warm
  extr_warm_ref = Tserie_na[ref_na>q95,]
  extr_warm_ref = cbind(extr_warm_ref[,1:3], extr_warm_ref[5])
  Ref_Extr_Warm = classification(extr_warm_ref)
  
  # Final result for reference serie
  Reference = c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","Length")
  mat_fin_ref = cbind(Reference, Ref_Extr_Cold, Ref_Cold, Ref_Mean, Ref_Warm, Ref_Extr_Warm)
  mat_fin_ref = as.data.frame(mat_fin_ref)
  for (i in c(2:6)) {
    mat_fin_ref[,i] = as.numeric(as.character(mat_fin_ref[,i]))
  }
  write.csv(mat_fin_ref, file=paste("03_",name,"_classes_reference.csv",sep=""), row.names=F)
  
  # Compute the difference in % in numbers of events between series in the same class
  freq = rbind(as.numeric(mat_fin_can[7,2:6]), as.numeric(mat_fin_ref[7,2:6]))
  perc = (freq[2,] - freq[1,]) / freq[2,]
  freq = rbind(freq, perc)
  colnames(freq) = c("Extr_Cold","Cold","Mean","Warm","Extr_Warm")
  rownames(freq) = c("Candidate events","Reference events","Increase in %")
  freq = round(freq, digits=4)
  write.csv(freq, file=paste("03_",name,"_frequencies_in_classes.csv",sep=""), row.names=T)
  
  # Compute where both stations have events in the same classes (common class events)
  extr_cold_common = Tserie_na[Tserie_na[,4]<=q05 & Tserie_na[,5]<=q05,]
  cold_common = Tserie_na[Tserie_na[,4]>q05 & Tserie_na[,4]<=q20 & Tserie_na[,5]>q05 & Tserie_na[,5]<=q20,]
  mean_common = Tserie_na[Tserie_na[,4]>q20 & Tserie_na[,4]<=q80 & Tserie_na[,5]>q20 & Tserie_na[,5]<=q80,]
  warm_common = Tserie_na[Tserie_na[,4]>q80 & Tserie_na[,4]<=q95 & Tserie_na[,5]>q80 & Tserie_na[,5]<=q95,]
  extr_warm_common = Tserie_na[Tserie_na[,4]>q95 & Tserie_na[,5]>q95,]
  
  # Compute RMSE, ME, cor and summary on common class events
  extr_cold_common_stat = as.data.frame(statistics(extr_cold_common[,4], extr_cold_common[,5]))
  cold_common_stat = as.data.frame(statistics(cold_common[,4],cold_common[,5]))
  mean_common_stat = as.data.frame(statistics(mean_common[,4],mean_common[,5]))
  warm_common_stat = as.data.frame(statistics(warm_common[,4],warm_common[,5]))
  extr_warm_common_stat = as.data.frame(statistics(extr_warm_common[,4],extr_warm_common[,5]))
  write.csv(extr_cold_common_stat, file=paste("03_",name,"_extr_cold.csv",sep=""), row.names=T)
  write.csv(cold_common_stat, file=paste("03_",name,"_cold.csv",sep=""), row.names=T)
  write.csv(mean_common_stat, file=paste("03_",name,"_mean.csv",sep=""), row.names=T)
  write.csv(warm_common_stat, file=paste("03_",name,"_warm.csv",sep=""), row.names=T)
  write.csv(extr_warm_common_stat, file=paste("03_",name,"_extr_warm.csv",sep=""), row.names=T)
  
  # Plot of differences on common class events
  pdf(paste("03_",name,"_boxplot_diff.pdf",sep=""))
  par(mfrow = c(3,2))
  boxplot(extr_cold_common[,7], col="darkblue", xlab="Extremely cold",
          ylab="Temperature (°C)", main=paste("Extremely cold",name))
  boxplot(cold_common[,7],col="blue", xlab="Cold", ylab="Temperature (°C)", main=paste("Cold",name))
  boxplot(mean_common[,7],col="green", xlab="Mean", ylab="Temperature (°C)", main=paste("Mean",name))
  boxplot(warm_common[,7],col="red", xlab="Warm", ylab="Temperature (°C)", main=paste("Warm",name))
  boxplot(extr_warm_common[,7],col="darkred", xlab="Extremely warm",
          ylab="Temperature (°C)", main=paste("Extremely warm",name))
  dev.off()
  
  # Prepare single file output in xlsx
  myexcel = paste(name,"_Results.xlsx",sep="")
  write.xlsx(stat, myexcel, sheetName="01 - Raw data")
  wb = loadWorkbook(myexcel)
  sheets = getSheets(wb)
  sheet = sheets[["01 - Raw data"]]
  addDataFrame(test_fin, sheet, startRow=11, startColumn=1, row.names=FALSE)
  sheet = createSheet(wb, sheetName = "02 - Monthly data")
  addDataFrame(mm_dif, sheet, startRow=1, startColumn=1)
  addDataFrame(trend, sheet, startRow=9, startColumn=1, col.names=FALSE)
  sheet = createSheet(wb, sheetName = "03 - Classes")
  addDataFrame(mat_fin_can, sheet, startRow=1, startColumn=1, row.names=FALSE)
  addDataFrame(mat_fin_ref, sheet, startRow=11, startColumn=1, row.names=FALSE)
  addDataFrame(freq, sheet, startRow=21, startColumn=1)
  addDataFrame(extr_cold_common_stat, sheet, startRow=27, startColumn=1)
  addDataFrame(cold_common_stat, sheet, startRow=36, startColumn=1)
  addDataFrame(mean_common_stat, sheet, startRow=45, startColumn=1)
  addDataFrame(warm_common_stat, sheet, startRow=54, startColumn=1)
  addDataFrame(extr_warm_common_stat, sheet, startRow=63, startColumn=1)
  saveWorkbook(wb, myexcel)
  
  return(NULL)
} # end of function analyze


###############################################################################
####                              MAIN PROGRAM                             ####
###############################################################################

# Read input data and output folder creation
if (interactive()) { # User has to choose file in input and path of results
  tab = read.table(file.choose(), header=T, na.strings="NA")
  setwd(choose.dir(getwd(), "Choose a folder to save your results:"))
} else { # File in input and path of results are hardcoded
  setwd("/data/test") # Path where input file is located
  tab = read.table("example1.txt", header=T, na.strings="NA") # Input file name
  dir.create("./results/", showWarnings = FALSE) # Results folder
  setwd("./results/") # Results folder
}

# Clean input and prepare other variables
tab$date = as.Date(tab[,1], format='%d/%m/%Y', tz="GMT")
tab = tab[complete.cases(tab[,1]), ]
tab$Y = as.numeric(format(tab$date, "%Y"))
tab$M = as.numeric(format(tab$date, "%m"))
tab$D = as.numeric(format(tab$date, "%d"))
start.tab = as.Date(tab$date[1])

# Compute TX and TN on both series
TX = data.frame(tab$Y, tab$M, tab$D, tab[4], tab[2])
colnames(TX) = c("Y", "M", "D", "CAN", "REF")
TX$time = ISOdate(TX[[1]], TX[[2]], TX[[3]], 0)
TN = data.frame(tab$Y, tab$M, tab$D, tab[5], tab[3])
colnames(TN) = c("Y", "M", "D", "CAN", "REF")
TN$time = ISOdate(TN[[1]], TN[[2]], TN[[3]], 0)

# Do analysis on TX and on TN
analyze("TX",TX)
analyze("TN",TN)

# Print runtime stats and clean all the environment
cat("Elaboration completed. You can find your results in", getwd(), "\n")
end.time = Sys.time()
cat(capture.output(end.time - start.time))
rm(list=ls())
