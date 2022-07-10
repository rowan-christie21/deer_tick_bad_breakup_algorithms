#Jan 22 2020
#creating pyramid plots for deer tick data in PA and recording results of multiple breakup, 
#absolute range, relative range, stabilitiy time, proportion significant, proportion significantly wrong functions
#
#Jan 27 2020
#Trying to resolve out of bounds subscript errors, and making sure PA dataset is completely subsetted
#
#dataset source: Zenodo
# study link:
#   https://zenodo.org/record/1476091#.Xij-T8hKg2w
# dataset citation:
#   Damie Pak, Steven B. Jacobs, & Joyce M. Sakamoto. (2018). Raw Dataset of Pennsylvania Passive Surveillance [Data set]. Zenodo. http://doi.org/10.5281/zenodo.1476091

#--------------------------------------------------------------------------#

#call needed libaries
library("ggplot2")
library(openxlsx)

#call to bad breakup algorithm developed by Dr. Christie Bahlai of Kent State
source("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/bad_breakup_2/R_model/bad_breakup_script.R")

#call to regime shift detector developed by Dr. Christie Bahlai of Kent State and Dr. Elise Zipkin
source("D:/Ixodes_scapularis_research_2019/regime_shift_tick_algorithms/tick_regime_09192019/monarch_regime-master_07162019/regime_shift_detector.R")

#sink(paste("pa_adult_tick_density_data_results",Sys.Date(),".txt", sep = ""))
#------------------------------------------------------------------------#
#Tick adult passive survellience count
print("Tick adult passive survellience results")

#read in tick adult data
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/main_dat_tick_FINAL_PA_01222020.csv", head = T)

#create vector of county in NY state in the dataset
county <- unique(tick_data$County.where.tick.found)
county <- na.omit(county)
county <- as.character(county)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on adult deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
county.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
num.phases = c()

pa_adult_passive_surv_dataset <- createWorkbook("pa_adult_passive_surv_dataset")
addWorksheet(pa_adult_active_surv_dataset, "pa_dataset")
writeData(pa_adult_active_surv_dataset, sheet = "pa_dataset", tick_data)

modifications <- c("Modifications",
                   "Subsetted tick species Ixodes scapularis",
                   "Subsetted by county",
                   "Subsetted by Host = Human, which represents passive survellience",
                   "Subsetted by adults",
                   "Removed nas",
                   "Aggregated adults by year")

addWorksheet(pa_adult_passive_surv_dataset, "modifications")
writeData(pa_adult_passive_surv_dataset, sheet = "modifications", modifications)

for (i in county) {
  tryCatch({
    #location: i
    county <- tick_data[tick_data$County.where.tick.found == i,]
    
    #subsetting species to only include Ixodes scapularis
    i_scapularis <- county[county$Tick_Taxon == "I. scapularis",]
    
    #subsetting host between human and not attached to host to define
    #different variables for passive and active sampling method
    passive_surv <- i_scapularis[i_scapularis$Host == "Human",]

    #get only year and ticks found column
    tick_adults <- county[,c("Year", "Adults")]
    #omit nas
    tick_adults <- na.omit(tick_adults)
    #aggregate all counts per year in single row
    frame=data.frame(tick_adults)
    tick_adults <- aggregate(frame['Adults'], by=frame['Year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run functions
    if(nrow(tick_adults) >= 10) {
      addWorksheet(pa_adult_passive_surv_dataset, i)
      writeData(pa_adult_passive_surv_dataset, sheet = i, tick_adults)
      
      #i returns grid name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      
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

      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(tick_adults)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)

      #print pyramid plot to png
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Adult deer tick density in county ", i, ", PA ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(tick_adults, title=paste("Adult deer tick density in county ", i, ", PA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

saveWorkbook(pa_adult_passive_surv_dataset, paste("pa_adult_passive_dataset", Sys.Date(),".xlsx", sept=""), overwrite = TRUE)

county_adult_passive_tick_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
write.csv(county_adult_passive_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/pa_adult_passive_tick_data_density_results_",Sys.Date(),".csv",sep = ""))

#------------------------------------------------------------------------#
#Tick adult active sampling count
print("Tick adult active sampling results")

#read in tick adult data
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/main_dat_tick_FINAL_PA_01222020.csv", head = T)

#create vector of county in NY state in the dataset
county <- unique(tick_data$County.where.tick.found)
county <- na.omit(county)
county <- as.character(county)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on adult deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
county.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
num.phases = c()

pa_adult_active_surv_dataset <- createWorkbook("pa_adult_active_surv_dataset")
addWorksheet(pa_adult_active_surv_dataset, "pa_dataset")
writeData(pa_adult_active_surv_dataset, sheet = "pa_dataset", tick_data)

modifications <- c("Modifications",
                   "Subsetted tick species Ixodes scapularis",
                   "Subsetted by county",
                   "Subsetted by Host = not attached to host, which represents active sampling",
                   "Subsetted by adults",
                   "Removed nas",
                   "Aggregated adults by year")

addWorksheet(pa_adult_active_surv_dataset, "modifications")
writeData(pa_adult_active_surv_dataset, sheet = "modifications", modifications)

for (i in county) {
  tryCatch({
    #location: i
    county <- tick_data[tick_data$County.where.tick.found == i,]
    
    #subsetting species to only include Ixodes scapularis
    i_scapularis <- county[county$Tick_Taxon == "I. scapularis",]
    
    #subsetting host between human and not attached to host to define
    #different variables for active and active sampling method
    active_surv <- i_scapularis[i_scapularis$Host == "not attached to host",]
    
    #get only year and ticks found column
    tick_adults <- county[,c("Year", "Adults")]
    #omit nas
    tick_adults <- na.omit(tick_adults)
    #aggregate all counts per year in single row
    frame=data.frame(tick_adults)
    tick_adults <- aggregate(frame['Adults'], by=frame['Year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run functions
    if(nrow(tick_adults) >= 10) {
      addWorksheet(pa_adult_active_surv_dataset, i)
      writeData(pa_adult_active_surv_dataset, sheet = i, tick_adults)
      
      #i returns grid name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      
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
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(tick_adults)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot to png
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Adult deer tick density in county ", i, ", PA ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(tick_adults, title=paste("Adult deer tick density in county ", i, ", PA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

saveWorkbook(pa_adult_active_surv_dataset, paste("pa_adult_active_dataset", Sys.Date(),".xlsx", sept=""), overwrite = TRUE)

county_adult_active_tick_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
write.csv(county_adult_active_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/pa_adult_active_tick_data_density_results_",Sys.Date(),".csv",sep = ""))

#------------------------------------------------------------------------#
#Tick nymph passive survellience count
print("Tick nymph passive survellience results")

#read in tick nymph data
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/main_dat_tick_FINAL_PA_01222020.csv", head = T)

#create vector of county in NY state in the dataset
county <- unique(tick_data$County.where.tick.found)
county <- na.omit(county)
county <- as.character(county)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on nymph deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
county.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
num.phases = c()

pa_nymph_passive_surv_dataset <- createWorkbook("pa_nymph_passive_surv_dataset")
addWorksheet(pa_nymph_passive_surv_dataset, "pa_dataset")
writeData(pa_nymph_passive_surv_dataset, sheet = "pa_dataset", tick_data)

modifications <- c("Modifications",
                   "Subsetted tick species Ixodes scapularis",
                   "Subsetted by county",
                   "Subsetted by Host = Human, which represents passive survellience",
                   "Subsetted by nymphs",
                   "Removed nas",
                   "Aggregated nymphs by year")

addWorksheet(pa_nymph_passive_surv_dataset, "modifications")
writeData(pa_nymph_passive_surv_dataset, sheet = "modifications", modifications)

for (i in county) {
  tryCatch({
    #location: i
    county <- tick_data[tick_data$County.where.tick.found == i,]
    
    #subsetting species to only include Ixodes scapularis
    i_scapularis <- county[county$Tick_Taxon == "I. scapularis",]
    
    #subsetting host between human and not attached to host to define
    #different variables for passive and passive sampling method
    passive_surv <- i_scapularis[i_scapularis$Host == "Human",]
    
    #get only year and ticks found column
    tick_nymphs <- county[,c("Year", "Nymph")]
    #omit nas
    tick_nymphs <- na.omit(tick_nymphs)
    #aggregate all counts per year in single row
    frame=data.frame(tick_nymphs)
    tick_nymphs <- aggregate(frame['Nymph'], by=frame['Year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run functions
    if(nrow(tick_nymphs) >= 10) {
      addWorksheet(pa_nymph_passive_surv_dataset, i)
      writeData(pa_nymph_passive_surv_dataset, sheet = i, tick_nymphs)
      
      #i returns grid name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      
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
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(tick_nymphs)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot to png
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/nymph deer tick density in county ", i, ", PA ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(tick_nymphs, title=paste("nymph deer tick density in county ", i, ", PA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

saveWorkbook(pa_nymph_passive_surv_dataset, paste("pa_nymph_passive_dataset", Sys.Date(),".xlsx", sept=""), overwrite = TRUE)

county_nymph_passive_tick_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
write.csv(county_nymph_passive_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/pa_nymph_passive_tick_data_density_results_",Sys.Date(),".csv",sep = ""))

#------------------------------------------------------------------------#
#Tick Nymph active sampling count
print("Tick Nymph active sampling results")

#read in tick Nymph data
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/main_dat_tick_FINAL_PA_01222020.csv", head = T)

#create vector of county in NY state in the dataset
county <- unique(tick_data$County.where.tick.found)
county <- na.omit(county)
county <- as.character(county)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on Nymph deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
county.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
num.phases = c()

pa_Nymph_active_surv_dataset <- createWorkbook("pa_Nymph_active_surv_dataset")
addWorksheet(pa_Nymph_active_surv_dataset, "pa_dataset")
writeData(pa_Nymph_active_surv_dataset, sheet = "pa_dataset", tick_data)

modifications <- c("Modifications",
                   "Subsetted tick species Ixodes scapularis",
                   "Subsetted by county",
                   "Subsetted by Host = not attached to host, which represents active sampling",
                   "Subsetted by Nymphs",
                   "Removed nas",
                   "Aggregated Nymphs by year")

addWorksheet(pa_Nymph_active_surv_dataset, "modifications")
writeData(pa_Nymph_active_surv_dataset, sheet = "modifications", modifications)

for (i in county) {
  tryCatch({
    #location: i
    county <- tick_data[tick_data$County.where.tick.found == i,]
    
    #subsetting species to only include Ixodes scapularis
    i_scapularis <- county[county$Tick_Taxon == "I. scapularis",]
    
    #subsetting host between human and not attached to host to define
    #different variables for active and active sampling method
    active_surv <- i_scapularis[i_scapularis$Host == "not attached to host",]
    
    #get only year and ticks found column
    tick_Nymphs <- county[,c("Year", "Nymph")]
    #omit nas
    tick_Nymphs <- na.omit(tick_Nymphs)
    #aggregate all counts per year in single row
    frame=data.frame(tick_Nymphs)
    tick_Nymphs <- aggregate(frame['Nymph'], by=frame['Year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run functions
    if(nrow(tick_Nymphs) >= 10) {
      addWorksheet(pa_Nymph_active_surv_dataset, i)
      writeData(pa_Nymph_active_surv_dataset, sheet = i, tick_Nymphs)
      
      #i returns grid name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, tick_Nymphs$Year[1])
      end.year <- c(end.year, tick_Nymphs$Year[nrow(tick_Nymphs)])
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(multiple_breakups(tick_Nymphs))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- stability_time(tick_Nymphs)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- abs_range(tick_Nymphs, only_significant = FALSE, significance = 0.05)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- relative_range(tick_Nymphs, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- proportion_significant(tick_Nymphs, significance = 0.05)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- proportion_wrong(tick_Nymphs, significance = 0.05)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(tick_Nymphs)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot to png
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Nymph deer tick density in county ", i, ", PA ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(tick_Nymphs, title=paste("Nymph deer tick density in county ", i, ", PA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

saveWorkbook(pa_Nymph_active_surv_dataset, paste("pa_Nymph_active_dataset", Sys.Date(),".xlsx", sept=""), overwrite = TRUE)

county_Nymph_active_tick_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
write.csv(county_Nymph_active_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/pa_Nymph_active_tick_data_density_results_",Sys.Date(),".csv",sep = ""))

#------------------------------------------------------------------------#
#Tick Larva passive survellience count
print("Tick Larva passive survellience results")

#read in tick Larva data
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/main_dat_tick_FINAL_PA_01222020.csv", head = T)

#create vector of county in NY state in the dataset
county <- unique(tick_data$County.where.tick.found)
county <- na.omit(county)
county <- as.character(county)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on Larva deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
county.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
num.phases = c()

pa_Larva_passive_surv_dataset <- createWorkbook("pa_Larva_passive_surv_dataset")
addWorksheet(pa_Larva_passive_surv_dataset, "pa_dataset")
writeData(pa_Larva_passive_surv_dataset, sheet = "pa_dataset", tick_data)

modifications <- c("Modifications",
                   "Subsetted tick species Ixodes scapularis",
                   "Subsetted by county",
                   "Subsetted by Host = Human, which represents passive survellience",
                   "Subsetted by Larvas",
                   "Removed nas",
                   "Aggregated Larvas by year")

addWorksheet(pa_Larva_passive_surv_dataset, "modifications")
writeData(pa_Larva_passive_surv_dataset, sheet = "modifications", modifications)

for (i in county) {
  tryCatch({
    #location: i
    county <- tick_data[tick_data$County.where.tick.found == i,]
    
    #subsetting species to only include Ixodes scapularis
    i_scapularis <- county[county$Tick_Taxon == "I. scapularis",]
    
    #subsetting host between human and not attached to host to define
    #different variables for passive and passive sampling method
    passive_surv <- i_scapularis[i_scapularis$Host == "Human",]
    
    #get only year and ticks found column
    tick_Larvas <- county[,c("Year", "Larva")]
    #omit nas
    tick_Larvas <- na.omit(tick_Larvas)
    #aggregate all counts per year in single row
    frame=data.frame(tick_Larvas)
    tick_Larvas <- aggregate(frame['Larva'], by=frame['Year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run functions
    if(nrow(tick_Larvas) >= 10) {
      addWorksheet(pa_Larva_passive_surv_dataset, i)
      writeData(pa_Larva_passive_surv_dataset, sheet = i, tick_Larvas)
      
      #i returns grid name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, tick_Larvas$Year[1])
      end.year <- c(end.year, tick_Larvas$Year[nrow(tick_Larvas)])
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(multiple_breakups(tick_Larvas))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- stability_time(tick_Larvas)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- abs_range(tick_Larvas, only_significant = FALSE, significance = 0.05)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- relative_range(tick_Larvas, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- proportion_significant(tick_Larvas, significance = 0.05)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- proportion_wrong(tick_Larvas, significance = 0.05)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(tick_Larvas)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot to png
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Larva deer tick density in county ", i, ", PA ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(tick_Larvas, title=paste("Larva deer tick density in county ", i, ", PA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

saveWorkbook(pa_Larva_passive_surv_dataset, paste("pa_Larva_passive_dataset", Sys.Date(),".xlsx", sept=""), overwrite = TRUE)

county_Larva_passive_tick_data_results <- data.frame(county.with.10.years.data, start.year, end.year)#, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
write.csv(county_Larva_passive_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/pa_Larva_passive_tick_data_density_results_",Sys.Date(),".csv",sep = ""))

#------------------------------------------------------------------------#
#Tick Larva active sampling count
print("Tick Larva active sampling results")

#read in tick Larva data
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/main_dat_tick_FINAL_PA_01222020.csv", head = T)

#create vector of county in NY state in the dataset
county <- unique(tick_data$County.where.tick.found)
county <- na.omit(county)
county <- as.character(county)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on Larva deer tick data for each county with at least 10 years data
start.year = c()
end.year = c()
county.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()
proportion.significant = c()
proportion.wrong = c()
num.phases = c()

pa_Larva_active_surv_dataset <- createWorkbook("pa_Larva_active_surv_dataset")
addWorksheet(pa_Larva_active_surv_dataset, "pa_dataset")
writeData(pa_Larva_active_surv_dataset, sheet = "pa_dataset", tick_data)

modifications <- c("Modifications",
                   "Subsetted tick species Ixodes scapularis",
                   "Subsetted by county",
                   "Subsetted by Host = not attached to host, which represents active sampling",
                   "Subsetted by Larvas",
                   "Removed nas",
                   "Aggregated Larvas by year")

addWorksheet(pa_Larva_active_surv_dataset, "modifications")
writeData(pa_Larva_active_surv_dataset, sheet = "modifications", modifications)

for (i in county) {
  tryCatch({
    #location: i
    county <- tick_data[tick_data$County.where.tick.found == i,]
    
    #subsetting species to only include Ixodes scapularis
    i_scapularis <- county[county$Tick_Taxon == "I. scapularis",]
    
    #subsetting host between human and not attached to host to define
    #different variables for active and active sampling method
    active_surv <- i_scapularis[i_scapularis$Host == "not attached to host",]
    
    #get only year and ticks found column
    tick_Larvas <- county[,c("Year", "Larva")]
    #omit nas
    tick_Larvas <- na.omit(tick_Larvas)
    #aggregate all counts per year in single row
    frame=data.frame(tick_Larvas)
    tick_Larvas <- aggregate(frame['Larva'], by=frame['Year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run functions
    if(nrow(tick_Larvas) >= 10) {
      addWorksheet(pa_Larva_active_surv_dataset, i)
      writeData(pa_Larva_active_surv_dataset, sheet = i, tick_Larvas)
      
      #i returns grid name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, tick_Larvas$Year[1])
      end.year <- c(end.year, tick_Larvas$Year[nrow(tick_Larvas)])
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(multiple_breakups(tick_Larvas))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- stability_time(tick_Larvas)
      print(st)
      stability.time <- c(stability.time, st)
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- abs_range(tick_Larvas, only_significant = FALSE, significance = 0.05)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- relative_range(tick_Larvas, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- proportion_significant(tick_Larvas, significance = 0.05)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- proportion_wrong(tick_Larvas, significance = 0.05)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(tick_Larvas)
      bf <- bestfit(ticks_modeled, "AIC")
      print(bf$Nfits)
      print(bf)
      num.phases <- c(num.phases, bf$Nfits)
      
      #print pyramid plot to png
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Larva deer tick density in county ", i, ", PA ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(tick_Larvas, title=paste("Larva deer tick density in county ", i, ", PA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

saveWorkbook(pa_Larva_active_surv_dataset, paste("pa_Larva_active_dataset", Sys.Date(),".xlsx", sept=""), overwrite = TRUE)

county_Larva_active_tick_data_results <- data.frame(county.with.10.years.data, start.year, end.year)#, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
write.csv(county_Larva_active_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/pa_Larva_active_tick_data_density_results_",Sys.Date(),".csv",sep = ""))

