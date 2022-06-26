#Jan 31 2020
#creating pyramid plots for deer tick data in CT and recording results of multiple breakup, 
#absolute range, relative range, stabilitiy time, proportion significant, proportion significantly wrong functions
#
#dataset source: Zenodo
# study link:
#   https://zenodo.org/record/1476091#.Xij-T8hKg2w
# dataset citation:
#   Damie Pak, Steven B. Jacobs, & Joyce M. Sakamoto. (2018). Raw Dataset of Pennsylvania Passive Surveillance [Data set]. Zenodo. http://doi.org/10.5281/zenodo.1476091

#--------------------------------------------------------------------------#

#call needed libaries
library("ggplot2")
library("openxlsx")

#call to bad breakup algorithm developed by Dr. Christie Bahlai of Kent State
source("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/bad_breakup_2/R_model/bad_breakup_script.R")

#call to regime shift detector developed by Dr. Christie Bahlai of Kent State and Dr. Elise Zipkin
source("D:/Ixodes_scapularis_research_2019/regime_shift_tick_algorithms/tick_regime_09192019/monarch_regime-master_07162019/regime_shift_detector.R")

#------------------------------------------------------------------------#
# Deer ticks found on people in CT towns
# life stage is not specified
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from CT
tick_data <- read.xlsx("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/ct_data/ct_agricultural_experiment_station_tick_data_01292020.xlsx")

#create vector of towns in CT sampled
towns <- unique(tick_data$Town)
towns <- na.omit(towns)
towns <- as.character(towns)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on adult deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
towns.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
num.phases = c()

sink(paste("ct_tick_count_results_",Sys.Date(),".txt", sep = ""))
print("Deer ticks found on people in CT towns")
for (i in towns) {
  tryCatch({
    #location: i
    tick_loc <- tick_data[tick_data$Town == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "Total.Identified")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #data now cleaned
    
    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(ticks_found) >= 10) {
      
      #i returns location name
      print(i)
      towns.with.10.years.data <- c(towns.with.10.years.data, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, ticks_found$Year[1])
      end.year <- c(end.year, ticks_found$Year[nrow(ticks_found)])
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(tryCatch(multiple_breakups(ticks_found), error=function(err) NA))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- tryCatch(stability_time(ticks_found), error=function(err) NA)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- tryCatch(abs_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- tryCatch(relative_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- tryCatch(proportion_significant(ticks_found, significance = 0.05), error=function(err) NA)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- tryCatch(proportion_wrong(ticks_found, significance = 0.05), error=function(err) NA)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <-  bestfit(ticks_modeled, "AIC")
      num_phases <- tryCatch(bf$Nfits, error=function(err) NA)
      print(num_phases)
      print(bf)
      if(towns.with.10.years.data[towns.with.10.years.data==i] == i) {
        num.phases <- c(num.phases, number_phases)
      }
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Deer ticks found on people in ", i, ", CT ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("Deer ticks found on people in ", i, ", CT", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

ct_tick_count_data_results <- data.frame(towns.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
ct_tick_count_data_results <- na.omit(ct_tick_count_data_results)
write.csv(ct_tick_count_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/ct_tick_count_results_",Sys.Date(),".csv",sep = ""))



#------------------------------------------------------------------------#
# Percent of deer ticks found with Borrelia in CT towns
# life stage is not specified
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from CT
tick_data <- read.xlsx("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/ct_data/ct_agricultural_experiment_station_tick_data_01292020.xlsx")

#create vector of towns in CT sampled
towns <- unique(tick_data$Town)
towns <- na.omit(towns)
towns <- as.character(towns)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on adult deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
towns.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
num.phases = c()

sink(paste("ct_tick_infection_results_",Sys.Date(),".txt", sep = ""))
print("Percent deer ticks positive for Borrelia burgdorferi in CT towns")
for (i in towns[6:7]) {
  tryCatch({
    #location: i
    tick_loc <- tick_data[tick_data$Town == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "Positive.for.Borrelia(%)")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #data now cleaned
    
    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(ticks_found) >= 10) {
      
      #i returns location name
      print(i)
      towns.with.10.years.data <- c(towns.with.10.years.data, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, ticks_found$Year[1])
      end.year <- c(end.year, ticks_found$Year[nrow(ticks_found)])
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(tryCatch(multiple_breakups(ticks_found), error=function(err) NA))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- tryCatch(stability_time(ticks_found), error=function(err) NA)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- tryCatch(abs_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- tryCatch(relative_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- tryCatch(proportion_significant(ticks_found, significance = 0.05), error=function(err) NA)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- tryCatch(proportion_wrong(ticks_found, significance = 0.05), error=function(err) NA)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <-  bestfit(ticks_modeled, "AIC")
      num_phases <- tryCatch(bf$Nfits, error=function(err) NA)
      print(num_phases)
      print(bf)
      if(towns.with.10.years.data[towns.with.10.years.data==i] == i) {
        num.phases <- c(num.phases, number_phases)
      }
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Percent deer ticks positive for Borrelia burgdorferi in ", i, ", CT ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("Percent deer ticks positive for Borrelia burgdorferi in ", i, ", CT", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

ct_tick_infection_data_results <- data.frame(towns.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
ct_tick_infection_data_results <- na.omit(ct_tick_infection_data_results)
write.csv(ct_tick_infection_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/ct_tick_infection_results_",Sys.Date(),".csv",sep = ""))
