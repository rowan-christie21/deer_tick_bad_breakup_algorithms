#Oct 17th 2019
#creating pyramid plots for deer tick data in Monmouth County and recording results of multiple breakup, 
#absolute range, relative range, stabilitiy time, proportion significant, proportion significantly wrong functions
#
#March 14th 2020
#Calculating proportion wrong before stability time for all towns
#
#dataset source: PLoS One
# study link:
#   https://datadryad.org/stash/dataset/doi:10.5061/dryad.d1c8046
# dataset citation:
#   Ostfeld RS, Levi T, Keesing F, Oggenfuss K, Canham CD (2018) Data from: Tick-borne disease risk in a forest food web. Dryad Digital Repository.  https://doi.org/10.5061/dryad.d1c8046

#--------------------------------------------------------------------------#

#call needed libaries
library("ggplot2")

#call to bad breakup algorithm developed by Dr. Christie Bahlai of Kent State
source("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/bad_breakup_2/R_model/bad_breakup_script.R")

#call to regime shift detector developed by Dr. Christie Bahlai of Kent State and Dr. Elise Zipkin
source("D:/Ixodes_scapularis_research_2019/regime_shift_tick_algorithms/tick_regime_09192019/monarch_regime-master_07162019/regime_shift_detector.R")

#------------------------------------------------------------------------#
# Adult deer ticks found on people in Monmouth County towns
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from Monmouth County
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/monmouth_county/Monmouth_County_tick_data_10162019.csv", head = T)

#specify adult deer tick life stage
tick_adults <- tick_data[tick_data$life.stage == "adult",]

#create vector of towns in Monmouth County sampled
towns <- unique(tick_adults$SUB_MUNI)
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
proportion.wrong.before.stability = c()
num.phases = c()

sink(paste("monmouth_county_adult_tick_density_data_results",Sys.Date(),".txt", sep = ""))
print("Adult deer ticks found on people in Monmouth County Submunicipalities")
for (i in towns) {
  tryCatch({
    #location: i
    tick_loc <- tick_adults[tick_adults$SUB_MUNI == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "COUNT")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #aggregate all counts per year in single row
    frame=data.frame(ticks_found)
    ticks_found <- aggregate(frame['COUNT'], by=frame['Year'], sum)
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
      print(multiple_breakups(ticks_found))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- stability_time(ticks_found)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- abs_range(ticks_found, only_significant = FALSE, significance = 0.05)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- relative_range(ticks_found, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- proportion_significant(ticks_found, significance = 0.05)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- proportion_wrong(ticks_found, significance = 0.05)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      #returns proportion wrong before stability time
      print("proportion wrong before stability time")
      pwbst <- proportion_wrong_before_stability(ticks_found, significance = 0.05)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Adult deer tick found on people in ", i, ", Monmouth County ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("Adult deer ticks found on people in ", i, ", Monmouth County", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

monmouth_county_adult_tick_data_results <- data.frame(towns.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, proportion.wrong.before.stability, num.phases)
write.csv(monmouth_county_adult_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/monmouth_county_adult_tick_data_density_results_",Sys.Date(),".csv",sep = ""))

##########################################
#now run the aggregated count in Monmouth County as a whole for deer ticks found on people
print("Monmouth County")

#read in survey data from Monmouth County
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/monmouth_county/Monmouth_County_tick_data_10162019.csv", head = T)

#specify adult deer tick life stage
tick_adults <- tick_data[tick_data$life.stage == "adult",]

#location: Monmouth County
tick_loc <- tick_data
#get only year and count column
ticks_found <- tick_loc[,c("Year", "COUNT")]
#omit nas
ticks_found <- na.omit(ticks_found)
#aggregate all counts per year in single row
frame=data.frame(ticks_found)
ticks_found <- aggregate(frame['COUNT'], by=frame['Year'], sum)
#data now cleaned

#print start and end year
print(ticks_found$Year[1])
print(ticks_found$Year[nrow(ticks_found)])
#Multiple breakups
print("Multiple Breakups")
print(multiple_breakups(ticks_found))
#Stability
print("Stability")
print(stability_time(ticks_found))
#Absolute range
print("Absolute Range")
print(abs_range(ticks_found, only_significant = FALSE, significance = 0.05))
#Relative Range
print("Relative Range")
print(relative_range(ticks_found, only_significant = FALSE, significance = 0.05))
#Proportion Significant
print("Proportion Significant")
print(proportion_significant(ticks_found, significance = 0.05))
#Proportion Significantly Wrong
print("Proportion Significantly Wrong")
print(proportion_wrong(ticks_found, significance = 0.05))
print("Proportion significantly wrong before stability time")
print(proportion_wrong_before_stability(ticks_found, significance = 0.05))
#Number of phases
print("Number of phases")
#add column for adundance next year to determine best model
ticks_modeled <- addNt1(ticks_found)
bf <- bestfit(ticks_modeled, "AIC")
print(bf$Nfits)
print(bf)
sink()

png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Adult deer ticks found on people in Monmouth County ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
pyramid_plot(ticks_found, title="Adult deer ticks found on people in Monmouth County", plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st))
dev.off()

#------------------------------------------------------------------------#
# Nymph deer ticks found on people in Monmouth County towns
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from Monmouth County
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/monmouth_county/Monmouth_County_tick_data_10162019.csv", head = T)

#specify nymph deer tick life stage
tick_nymphs <- tick_data[tick_data$life.stage == "nymph",]

#create vector of towns in Monmouth County sampled
towns <- unique(tick_nymphs$SUB_MUNI)
towns <- na.omit(towns)
towns <- as.character(towns)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on nymph deer tick data for each county with at least 10 years data
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
proportion.wrong.before.stability = c()
num.phases = c()

sink(paste("monmouth_county_nymph_tick_density_data_results",Sys.Date(),".txt", sep = ""))
print("nymph deer ticks found on people in Monmouth County Submunicipalities")
for (i in towns) {
  tryCatch({
    #location: i
    tick_loc <- tick_nymphs[tick_nymphs$SUB_MUNI == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "COUNT")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #aggregate all counts per year in single row
    frame=data.frame(ticks_found)
    ticks_found <- aggregate(frame['COUNT'], by=frame['Year'], sum)
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
      print(multiple_breakups(ticks_found))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- stability_time(ticks_found)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- abs_range(ticks_found, only_significant = FALSE, significance = 0.05)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- relative_range(ticks_found, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- proportion_significant(ticks_found, significance = 0.05)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- proportion_wrong(ticks_found, significance = 0.05)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      #returns proportion wrong before stability time
      print("proportion wrong before stability time")
      pwbst <- proportion_wrong_before_stability(ticks_found, significance = 0.05)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/nymph deer tick found on people in ", i, ", Monmouth County ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("Nymph deer ticks found on people in ", i, ", Monmouth County", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

monmouth_county_nymph_tick_data_results <- data.frame(towns.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, proportion.wrong.before.stability, num.phases)
write.csv(monmouth_county_nymph_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/monmouth_county_nymph_tick_data_density_results_",Sys.Date(),".csv",sep = ""))

##########################################
#now run the aggregated count in Monmouth County as a whole for deer ticks found on people
print("Monmouth County")

#read in survey data from Monmouth County
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/monmouth_county/Monmouth_County_tick_data_10162019.csv", head = T)

#specify nymph deer tick life stage
tick_nymphs <- tick_data[tick_data$life.stage == "nymph",]

#location: Monmouth County
tick_loc <- tick_data
#get only year and count column
ticks_found <- tick_loc[,c("Year", "COUNT")]
#omit nas
ticks_found <- na.omit(ticks_found)
#aggregate all counts per year in single row
frame=data.frame(ticks_found)
ticks_found <- aggregate(frame['COUNT'], by=frame['Year'], sum)
#data now cleaned

#print start and end year
print(ticks_found$Year[1])
print(ticks_found$Year[nrow(ticks_found)])
#Multiple breakups
print("Multiple Breakups")
print(multiple_breakups(ticks_found))
#Stability
print("Stability")
print(stability_time(ticks_found))
#Absolute range
print("Absolute Range")
print(abs_range(ticks_found, only_significant = FALSE, significance = 0.05))
#Relative Range
print("Relative Range")
print(relative_range(ticks_found, only_significant = FALSE, significance = 0.05))
#Proportion Significant
print("Proportion Significant")
print(proportion_significant(ticks_found, significance = 0.05))
#Proportion Significantly Wrong
print("Proportion Significantly Wrong")
print(proportion_wrong(ticks_found, significance = 0.05))
print("Proportion significantly wrong before stability time")
print(proportion_wrong_before_stability(ticks_found, significance = 0.05))
#Number of phases
print("Number of phases")
#add column for adundance next year to determine best model
ticks_modeled <- addNt1(ticks_found)
bf <- bestfit(ticks_modeled, "AIC")
print(bf$Nfits)
print(bf)
sink()

png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/nymph deer ticks found on people in Monmouth County ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
pyramid_plot(ticks_found, title="Nymph deer ticks found on people in Monmouth County", plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st))
dev.off()


#------------------------------------------------------------------------#
# larval deer ticks found on people in Monmouth County towns
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from Monmouth County
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/monmouth_county/Monmouth_County_tick_data_10162019.csv", head = T)

#specify larval deer tick life stage
tick_larvals <- tick_data[tick_data$life.stage == "larva",]

#create vector of towns in Monmouth County sampled
towns <- unique(tick_larvals$SUB_MUNI)
towns <- na.omit(towns)
towns <- as.character(towns)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on larval deer tick data for each county with at least 10 years data
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
proportion.wrong.before.stability = c()
num.phases = c()

sink(paste("monmouth_county_larval_tick_density_data_results",Sys.Date(),".txt", sep = ""))
print("larval deer ticks found on people in Monmouth County Submunicipalities")
for (i in towns) {
  tryCatch({
    #location: i
    tick_loc <- tick_larvals[tick_larvals$SUB_MUNI == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "COUNT")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #aggregate all counts per year in single row
    frame=data.frame(ticks_found)
    ticks_found <- aggregate(frame['COUNT'], by=frame['Year'], sum)
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
      print(multiple_breakups(ticks_found))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- stability_time(ticks_found)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- abs_range(ticks_found, only_significant = FALSE, significance = 0.05)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- relative_range(ticks_found, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- proportion_significant(ticks_found, significance = 0.05)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- proportion_wrong(ticks_found, significance = 0.05)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      #returns proportion wrong before stability time
      print("proportion wrong before stability time")
      pwbst <- proportion_wrong_before_stability(ticks_found, significance = 0.05)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/larval deer tick found on people in ", i, ", Monmouth County ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("Larval deer ticks found on people in ", i, ", Monmouth County", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

monmouth_county_larval_tick_data_results <- data.frame(towns.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, proportion.wrong.before.stability, num.phases)
write.csv(monmouth_county_larval_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/monmouth_county_larval_tick_data_density_results_",Sys.Date(),".csv",sep = ""))

##########################################
#now run the aggregated count in Monmouth County as a whole for deer ticks found on people
print("Monmouth County")

#read in survey data from Monmouth County
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/monmouth_county/Monmouth_County_tick_data_10162019.csv", head = T)

#specify larval deer tick life stage
tick_larvals <- tick_data[tick_data$life.stage == "larva",]

#location: Monmouth County
tick_loc <- tick_data
#get only year and count column
ticks_found <- tick_loc[,c("Year", "COUNT")]
#omit nas
ticks_found <- na.omit(ticks_found)
#aggregate all counts per year in single row
frame=data.frame(ticks_found)
ticks_found <- aggregate(frame['COUNT'], by=frame['Year'], sum)
#data now cleaned

#print start and end year
print(ticks_found$Year[1])
print(ticks_found$Year[nrow(ticks_found)])
#Multiple breakups
print("Multiple Breakups")
print(multiple_breakups(ticks_found))
#Stability
print("Stability")
print(stability_time(ticks_found))
#Absolute range
print("Absolute Range")
print(abs_range(ticks_found, only_significant = FALSE, significance = 0.05))
#Relative Range
print("Relative Range")
print(relative_range(ticks_found, only_significant = FALSE, significance = 0.05))
#Proportion Significant
print("Proportion Significant")
print(proportion_significant(ticks_found, significance = 0.05))
#Proportion Significantly Wrong
print("Proportion Significantly Wrong")
print(proportion_wrong(ticks_found, significance = 0.05))
print("Proportion significantly wrong before stability time")
print(proportion_wrong_before_stability(ticks_found, significance = 0.05))
#Number of phases
print("Number of phases")
#add column for adundance next year to determine best model
ticks_modeled <- addNt1(ticks_found)
bf <- bestfit(ticks_modeled, "AIC")
print(bf$Nfits)
print(bf)
sink()

png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/larval deer ticks found on people in Monmouth County ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
pyramid_plot(ticks_found, title="Larval deer ticks found on people in Monmouth County", plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st))
dev.off()



