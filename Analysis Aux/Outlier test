library(car)
#Bonferroni Outlier Test
out <- outlierTest(modPH, cutoff = Inf, n.max = Inf)
table_pvalue <-out$p #Save the p values in a vector
filter_p <- table_pvalue[table_pvalue<0.05] # filter to have a TRUE/FALSE list of the p-values < 0.05
count_outliers<- sum(table_pvalue<0.05, na.rm = TRUE) #count outliers to be removed
alldata_noOut[alldata_noOut$PlotMN %in% alldata_noOut$PlotMN[as.numeric(names(filter_p))],] <- NA #remove the plots that have p-values < 0.05
message(count_outliers, " Outliers will be removed")
