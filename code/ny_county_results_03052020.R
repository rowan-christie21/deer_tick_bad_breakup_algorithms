#Oct 3rd 2019
#creating pyramid plots for deer tick data in NY counties and recording results of multiple breakup function
#dataset source: New York Deptartment of Health
#
#March 5th 2020
#Calculating proportion wrong before stability time for all grids
#
# dataset for tick adults used:
#   https://health.data.ny.gov/Health/Deer-Tick-Surveillance-Adults-Oct-to-Dec-excluding/vzbp-i2d4
# citation:
#   New York State Department of Health Office of Public Health. 2019. Deer Tick Surveillance: Adults (Oct to Dec) excluding Powassan virus: Beginning 2008. https://health.data.ny.gov/Health/Deer-Tick-Surveillance-Nymphs-May-to-Sept-excludin/kibp-u2ip
# dataset for tick nymphs used:
#   https://health.data.ny.gov/Health/Deer-Tick-Surveillance-Nymphs-May-to-Sept-excludin/kibp-u2ip
# citation:
#   New York State Department of Health Office of Public Health. 2019. Access Nymph Deer Tick Collection Data by County (Excluding Powassan Virus). https://health.data.ny.gov/Health/Deer-Tick-Surveillance-Nymphs-May-to-Sept-excludin/kibp-u2ip

#--------------------------------------------------------------------------#

#call needed libaries
library("ggplot2")

#call to bad breakup algorithm developed by Dr. Christie Bahlai of Kent State
source("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/bad_breakup_2/R_model/bad_breakup_script.R")

#-------------------------------------------------#
#Tick adult density

#read in tick adult data
tick_adults_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/NY_Dept_health/Deer_Tick_Surveillance__Adults__Oct_to_Dec__excluding_Powassan_virus__Beginning_2008_07222019.csv", head = T)

#create vector of counties in NY state in the dataset
counties <- unique(tick_adults_data$County)
counties <- na.omit(counties)
counties <- as.character(counties)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on adult deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
counties.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
proportion.wrong.before.stability = c()
num.phases = c()

sink(paste("ny_county_adult_tick_data_results",Sys.Date()))
for (i in counties) {
  tryCatch({
    #location: i
    county <- tick_adults_data[tick_adults_data$County == i,]
    #get only year and ticks found column
    tick_adults <- county[,c("Year", "Tick.Population.Density")]
    #omit nas
    tick_adults <- na.omit(tick_adults)
    #data now cleaned

    #if number of years in location exceeds 9 run functions
    if(nrow(tick_adults) >= 10) {
      #i returns county name
      print(i)
      counties.with.10.years.data <- c(counties.with.10.years.data, i)
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, tick_adults$Year[1])
      end.year <- c(end.year, tick_adults$Year[nrow(tick_adults)])
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(multiple_breakups(tick_adults))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- stability_time(tick_adults)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- abs_range(tick_adults, only_significant = FALSE, significance = 0.05)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- relative_range(tick_adults, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- proportion_significant(tick_adults, significance = 0.05)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- proportion_wrong(tick_adults, significance = 0.05)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      #returns proportion wrong before stability time
      print("proportion wrong before stability time")
      pwbst <- proportion_wrong_before_stability(tick_adults, significance = 0.05)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(tick_adults)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot to png
      #png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Adult deer tick density in ", i, ", NY ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      #print(pyramid_plot(tick_adults, title=paste("Adult deer tick density in ", i, ", NY", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      #dev.off()

    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

county_adult_tick_data_results <- data.frame(counties.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, proportion.wrong.before.stability, num.phases)
write.csv(county_adult_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/ny_county_adult_tick_data_results_",Sys.Date(),".csv",sep = ""))

#-------------------------------------------------#
#Tick adult pathogen presence

#read in tick adult data
tick_adults_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/NY_Dept_health/Deer_Tick_Surveillance__Adults__Oct_to_Dec__excluding_Powassan_virus__Beginning_2008_07222019.csv", head = T)

#create vector of counties in NY state in the dataset
counties <- unique(tick_adults_data$County)
counties <- na.omit(counties)
counties <- as.character(counties)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on adult deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
counties.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
proportion.wrong.before.stability = c()
num.phases = c()

sink(paste("ny_county_adult_tick_pathogen_presence_data_results",Sys.Date()))
for (i in counties) {
  tryCatch({
    #location: i
    county <- tick_adults_data[tick_adults_data$County == i,]
    #get only year and ticks found column
    tick_pathogen <- county[,c("Year", "B..burgdorferi....")]
    #omit nas
    tick_pathogen <- na.omit(tick_pathogen)
    #data now cleaned
    
    #if number of years in location exceeds 9 run functions
    if(nrow(tick_pathogen) >= 10) {
      #i returns county name
      print(i)
      counties.with.10.years.data <- c(counties.with.10.years.data, i)
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, tick_pathogen$Year[1])
      end.year <- c(end.year, tick_pathogen$Year[nrow(tick_pathogen)])
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(multiple_breakups(tick_pathogen))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- stability_time(tick_pathogen)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- abs_range(tick_pathogen, only_significant = FALSE, significance = 0.05)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- relative_range(tick_pathogen, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- proportion_significant(tick_pathogen, significance = 0.05)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- proportion_wrong(tick_pathogen, significance = 0.05)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      #returns proportion wrong before stability time
      print("proportion wrong before stability time")
      pwbst <- proportion_wrong_before_stability(tick_pathogen, significance = 0.05)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(tick_pathogen)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot to png
      #png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/B Burgdorferi % in deer tick adults in ", i, ", NY ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      #print(pyramid_plot(tick_pathogen, title=paste("B. Burgdorferi % in deer tick adults in ", i, ", NY", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      #dev.off()
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

county_adult_tick_data_results <- data.frame(counties.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, proportion.wrong.before.stability, num.phases)
write.csv(county_adult_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/ny_county_adult_tick_pathogen_presence_data_results_",Sys.Date(),".csv",sep = ""))

#-------------------------------------------------#
#Tick nymph density

#read in tick nymph data
tick_nymphs_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/NY_Dept_health/Deer_Tick_Surveillance__Nymphs__May_to_Sept__excluding_Powassan_virus__Beginning_2008_07222019.csv", head = T)

#create vector of counties in NY state in the dataset
counties <- unique(tick_nymphs_data$County)
counties <- na.omit(counties)
counties <- as.character(counties)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on nymph deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
counties.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
proportion.wrong.before.stability = c()
num.phases = c()

sink(paste("ny_county_nymph_tick_data_results",Sys.Date()))
for (i in counties) {
  tryCatch({
    #location: i
    county <- tick_nymphs_data[tick_nymphs_data$County == i,]
    #get only year and ticks found column
    tick_nymphs <- county[,c("Year", "Tick.Population.Density")]
    #omit nas
    tick_nymphs <- na.omit(tick_nymphs)
    #data now cleaned
    
    #if number of years in location exceeds 9 run functions
    if(nrow(tick_nymphs) >= 10) {
      #i returns county name
      print(i)
      counties.with.10.years.data <- c(counties.with.10.years.data, i)
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, tick_nymphs$Year[1])
      end.year <- c(end.year, tick_nymphs$Year[nrow(tick_nymphs)])
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(multiple_breakups(tick_nymphs))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- stability_time(tick_nymphs)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- abs_range(tick_nymphs, only_significant = FALSE, significance = 0.05)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- relative_range(tick_nymphs, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- proportion_significant(tick_nymphs, significance = 0.05)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- proportion_wrong(tick_nymphs, significance = 0.05)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      #returns proportion wrong before stability time
      print("proportion wrong before stability time")
      pwbst <- proportion_wrong_before_stability(tick_nymphs, significance = 0.05)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(tick_nymphs)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot to png
      #png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/nymph deer tick density in ", i, ", NY ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      #print(pyramid_plot(tick_nymphs, title=paste("Nymph deer tick density in ", i, ", NY", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      #dev.off()
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

county_nymph_tick_data_results <- data.frame(counties.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, proportion.wrong.before.stability, num.phases)
write.csv(county_nymph_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/ny_county_nymph_tick_data_results_",Sys.Date(),".csv",sep = ""))

#-------------------------------------------------#
#Tick nymph pathogen presence

#read in tick nymph data
tick_nymphs_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/Deer_Tick_Surveillance__Nymphs__May_to_Sept__excluding_Powassan_virus__Beginning_2008.csv", head = T)

#create vector of counties in NY state in the dataset
counties <- unique(tick_nymphs_data$County)
counties <- na.omit(counties)
counties <- as.character(counties)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on nymph deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
counties.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
proportion.wrong.before.stability = c()
num.phases = c()

sink(paste("ny_county_nymph_tick_pathogen_presence_data_results",Sys.Date()))
for (i in counties) {
  tryCatch({
    #location: i
    county <- tick_nymphs_data[tick_nymphs_data$County == i,]
    #get only year and ticks found column
    tick_pathogen <- county[,c("Year", "B..burgdorferi....")]
    #omit nas
    tick_pathogen <- na.omit(tick_pathogen)
    #data now cleaned
    
    #if number of years in location exceeds 9 run functions
    if(nrow(tick_pathogen) >= 10) {
      #i returns county name
      print(i)
      counties.with.10.years.data <- c(counties.with.10.years.data, i)
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, tick_pathogen$Year[1])
      end.year <- c(end.year, tick_pathogen$Year[nrow(tick_pathogen)])
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(multiple_breakups(tick_pathogen))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- stability_time(tick_pathogen)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- abs_range(tick_pathogen, only_significant = FALSE, significance = 0.05)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- relative_range(tick_pathogen, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- proportion_significant(tick_pathogen, significance = 0.05)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- proportion_wrong(tick_pathogen, significance = 0.05)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      #returns proportion wrong before stability time
      print("proportion wrong before stability time")
      pwbst <- proportion_wrong_before_stability(tick_pathogen, significance = 0.05)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(tick_pathogen)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot to png
      #png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/B Burgdorferi % in deer tick nymphs in ", i, ", NY ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      #print(pyramid_plot(tick_pathogen, title=paste("B. Burgdorferi % in deer tick nymphs in ", i, ", NY", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      #dev.off()
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

county_nymph_tick_data_results <- data.frame(counties.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, proportion.wrong.before.stability, num.phases)
write.csv(county_nymph_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/ny_county_nymph_tick_pathogen_presence_data_results_",Sys.Date(),".csv",sep = ""))


