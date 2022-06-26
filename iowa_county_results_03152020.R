#Feb 5 2020
#creating pyramid plots for deer tick data in Iowa and recording results of multiple breakup, 
#absolute range, relative range, stabilitiy time, proportion significant, proportion significantly wrong functions
#
#March 15th 2020
#Calculating proportion wrong before stability time for all grids
#
#dataset source: Zenodo
# study link:
#   https://datadryad.org/stash/dataset/doi:10.5061/dryad.n2k66
# dataset citation:
# Oliver, Jonathan D.; Bennett, Steve W.; Beati, Lorenza; Bartholomay, Lyric C. (2018), Data from: Range expansion and increasing Borrelia burgdorferi infection of the tick Ixodes scapularis (Acari: Ixodidae) in Iowa, 1990-2013, Dryad, Dataset, https://doi.org/10.5061/dryad.n2k66
#
#--------------------------------------------------------------------------#

#call needed libaries
library("ggplot2")
library("openxlsx")

#call to bad breakup algorithm developed by Dr. Christie Bahlai of Kent State
source("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/bad_breakup_2/R_model/bad_breakup_script.R")

#call to regime shift detector developed by Dr. Christie Bahlai of Kent State and Dr. Elise Zipkin
source("D:/Ixodes_scapularis_research_2019/regime_shift_tick_algorithms/tick_regime_09192019/monarch_regime-master_07162019/regime_shift_detector.R")

#------------------------------------------------------------------------#
# Adult deer ticks found on people in Iowa counties
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from IA
tick_data <- read.xlsx("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/iowa_data/IA_LDS1990-2013_Iscap_02032020.xlsx")

#subset by adults
tick_adults <- tick_data[tick_data$Stage == "A",]

#create vector of county in IA sampled
county <- unique(tick_adults$County)
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
proportion.wrong.before.stability = c()
num.phases = c()

adult_tick_count_vect = c()

sink(paste("iowa_tick_adult_count_results_",Sys.Date(),".txt", sep = ""))
print("Adult deer ticks found on people in Iowa counties")
for (i in county) {
  tryCatch({
    #location: i
    tick_loc <- tick_adults[tick_adults$County == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "Count")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #aggregate all counts per year in single row
    frame=data.frame(ticks_found)
    ticks_found <- aggregate(frame['Count'], by=frame['Year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(ticks_found) >= 10) {
      
      #i returns location name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      adult_tick_count_vect <- c(adult_tick_count_vect, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, ticks_found$Year[1])
      end.year <- c(end.year, ticks_found$Year[nrow(ticks_found)])
      adult_tick_count_vect <- c(adult_tick_count_vect, paste("start year for ", i, ":",  ticks_found$Year[1]))
      adult_tick_count_vect <- c(adult_tick_count_vect, paste("end year for ", i, ":", ticks_found$Year[nrow(ticks_found)]))
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(tryCatch(multiple_breakups(ticks_found), error=function(err) NA))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- tryCatch(stability_time(ticks_found), error=function(err) NA)
      print(st)
      stability.time <- c(stability.time, st)
      adult_tick_count_vect <- c(adult_tick_count_vect, paste("stability time for ", i, ":", st))
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- tryCatch(abs_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      adult_tick_count_vect <- c(adult_tick_count_vect, paste("absolute range min for ", i, ":", ar[1]))
      adult_tick_count_vect <- c(adult_tick_count_vect, paste("absolute range max for ", i, ":", ar[2]))
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- tryCatch(relative_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      adult_tick_count_vect <- c(adult_tick_count_vect, paste("relative range min for ", i, ":", rr[1]))
      adult_tick_count_vect <- c(adult_tick_count_vect, paste("relative range max for ", i, ":", rr[2]))
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- tryCatch(proportion_significant(ticks_found, significance = 0.05), error=function(err) NA)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      adult_tick_count_vect <- c(adult_tick_count_vect, paste("proportion significant for ", i, ":", ps))
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- tryCatch(proportion_wrong(ticks_found, significance = 0.05), error=function(err) NA)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      adult_tick_count_vect <- c(adult_tick_count_vect, paste("proportion significantly wrong for ", i, ":", psw))
      
      #returns proportion wrong before stability time
      print("proportion wrong before stability time")
      pwbst <- tryCatch(proportion_wrong_before_stability(ticks_found, significance = 0.05), error=function(err) NA)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      adult_tick_count_vect <- c(adult_tick_count_vect, paste("proportion significantly wrong before stability time for ", i, ":", pwbst))

      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <-  bestfit(ticks_modeled, "AIC")
      num_phases <- tryCatch(bf$Nfits, error=function(err) NA)
      print(num_phases)
      print(bf)
      if(county.with.10.years.data[county.with.10.years.data==i] == i) {
        num.phases <- c(num.phases, num_phases)
        adult_tick_count_vect <- c(adult_tick_count_vect, paste("number of phases for ", i, ":", num_phases))
      }
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Adult deer ticks found on people in ", i, ", IA ", (Sys.Date()),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("Adult deer ticks found on people in ", i, ", IA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

write.csv(adult_tick_count_vect, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_adult_tick_count_text_results_",Sys.Date(),".csv",sep = ""))

iowa_adult_tick_count_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
iowa_adult_tick_count_data_results <- na.omit(iowa_adult_tick_count_data_results)
write.csv(iowa_adult_tick_count_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_adult_tick_count_results_",(Sys.Date()),".csv",sep = ""))


#------------------------------------------------------------------------#
# nymph deer ticks found on people in Iowa counties
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from IA
tick_data <- read.xlsx("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/iowa_data/IA_LDS1990-2013_Iscap_02032020.xlsx")

#subset by nymphs
tick_nymphs <- tick_data[tick_data$Stage == "N",]

#create vector of county in IA sampled
county <- unique(tick_nymphs$County)
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

nymph_tick_count_vect = c()

sink(paste("iowa_tick_nymph_count_results_",(Sys.Date()),".txt", sep = ""))
print("nymph deer ticks found on people in Iowa counties")
for (i in county) {
  tryCatch({
    #location: i
    tick_loc <- tick_nymphs[tick_nymphs$County == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "Count")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #aggregate all counts per year in single row
    frame=data.frame(ticks_found)
    ticks_found <- aggregate(frame['Count'], by=frame['Year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(ticks_found) >= 10) {
      
      #i returns location name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, ticks_found$Year[1])
      end.year <- c(end.year, ticks_found$Year[nrow(ticks_found)])
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("start year for ", i, ":",  ticks_found$Year[1]))
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("end year for ", i, ":", ticks_found$Year[nrow(ticks_found)]))
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(tryCatch(multiple_breakups(ticks_found), error=function(err) NA))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- tryCatch(stability_time(ticks_found), error=function(err) NA)
      print(st)
      stability.time <- c(stability.time, st)
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("stability time for ", i, ":", st))
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- tryCatch(abs_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("absolute range min for ", i, ":", ar[1]))
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("absolute range max for ", i, ":", ar[2]))
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- tryCatch(relative_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("relative range min for ", i, ":", rr[1]))
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("relative range max for ", i, ":", rr[2]))
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- tryCatch(proportion_significant(ticks_found, significance = 0.05), error=function(err) NA)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("proportion significant for ", i, ":", ps))
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- tryCatch(proportion_wrong(ticks_found, significance = 0.05), error=function(err) NA)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("proportion significantly wrong for ", i, ":", psw))
      
      print("proportion wrong before stability time")
      pwbst <- tryCatch(proportion_wrong_before_stability(ticks_found, significance = 0.05), error=function(err) NA)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("proportion significantly wrong before stability time for ", i, ":", pwbst))
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <-  bestfit(ticks_modeled, "AIC")
      num_phases <- tryCatch(bf$Nfits, error=function(err) NA)
      print(num_phases)
      print(bf)
      if(county.with.10.years.data[county.with.10.years.data==i] == i) {
        num.phases <- c(num.phases, num_phases)
        nymph_tick_count_vect <- c(nymph_tick_count_vect, paste("number of phases for ", i, ":", num_phases))
      }
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/nymph deer ticks found on people in ", i, ", IA ", (Sys.Date()),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("nymph deer ticks found on people in ", i, ", IA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

write.csv(nymph_tick_count_vect, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_nymph_tick_count_text_results_",(Sys.Date()),".csv",sep = ""))

iowa_nymph_tick_count_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
iowa_nymph_tick_count_data_results <- na.omit(iowa_nymph_tick_count_data_results)
write.csv(iowa_nymph_tick_count_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_nymph_tick_count_results_",(Sys.Date()),".csv",sep = ""))


#------------------------------------------------------------------------#
# larval deer ticks found on people in Iowa counties
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from IA
tick_data <- read.xlsx("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/iowa_data/IA_LDS1990-2013_Iscap_02032020.xlsx")

#subset by larvals
tick_larvals <- tick_data[tick_data$Stage == "L",]

#create vector of county in IA sampled
county <- unique(tick_larvals$County)
county <- na.omit(county)
county <- as.character(county)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on larval deer tick data for each county with at least 10 years data
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

larval_tick_count_vect = c()

sink(paste("iowa_tick_larval_count_results_",(Sys.Date()),".txt", sep = ""))
print("larval deer ticks found on people in Iowa counties")
for (i in county) {
  tryCatch({
    #location: i
    tick_loc <- tick_larvals[tick_larvals$County == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "Count")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #aggregate all counts per year in single row
    frame=data.frame(ticks_found)
    ticks_found <- aggregate(frame['Count'], by=frame['Year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(ticks_found) >= 10) {
      
      #i returns location name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      larval_tick_count_vect <- c(larval_tick_count_vect, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, ticks_found$Year[1])
      end.year <- c(end.year, ticks_found$Year[nrow(ticks_found)])
      larval_tick_count_vect <- c(larval_tick_count_vect, paste("start year for ", i, ":",  ticks_found$Year[1]))
      larval_tick_count_vect <- c(larval_tick_count_vect, paste("end year for ", i, ":", ticks_found$Year[nrow(ticks_found)]))
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(tryCatch(multiple_breakups(ticks_found), error=function(err) NA))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- tryCatch(stability_time(ticks_found), error=function(err) NA)
      print(st)
      stability.time <- c(stability.time, st)
      larval_tick_count_vect <- c(larval_tick_count_vect, paste("stability time for ", i, ":", st))
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- tryCatch(abs_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      larval_tick_count_vect <- c(larval_tick_count_vect, paste("absolute range min for ", i, ":", ar[1]))
      larval_tick_count_vect <- c(larval_tick_count_vect, paste("absolute range max for ", i, ":", ar[2]))
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- tryCatch(relative_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      larval_tick_count_vect <- c(larval_tick_count_vect, paste("relative range min for ", i, ":", rr[1]))
      larval_tick_count_vect <- c(larval_tick_count_vect, paste("relative range max for ", i, ":", rr[2]))
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- tryCatch(proportion_significant(ticks_found, significance = 0.05), error=function(err) NA)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      larval_tick_count_vect <- c(larval_tick_count_vect, paste("proportion significant for ", i, ":", ps))
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- tryCatch(proportion_wrong(ticks_found, significance = 0.05), error=function(err) NA)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      larval_tick_count_vect <- c(larval_tick_count_vect, paste("proportion significantly wrong for ", i, ":", psw))
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <-  bestfit(ticks_modeled, "AIC")
      num_phases <- tryCatch(bf$Nfits, error=function(err) NA)
      print(num_phases)
      print(bf)
      if(county.with.10.years.data[county.with.10.years.data==i] == i) {
        num.phases <- c(num.phases, num_phases)
        larval_tick_count_vect <- c(larval_tick_count_vect, paste("number of phases for ", i, ":", num_phases))
      }
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/larval deer ticks found on people in ", i, ", IA ", (Sys.Date()),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("larval deer ticks found on people in ", i, ", IA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

write.csv(larval_tick_count_vect, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_larval_tick_count_text_results_",(Sys.Date()),".csv",sep = ""))

iowa_larval_tick_count_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
iowa_larval_tick_count_data_results <- na.omit(iowa_larval_tick_count_data_results)
write.csv(iowa_larval_tick_count_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_larval_tick_count_results_",(Sys.Date()),".csv",sep = ""))




#------------------------------------------------------------------------#
# Percent of adult deer ticks infected with B. burgdoferi in Iowa counties
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from IA
tick_data <- read.xlsx("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/iowa_data/IA_LDS1990-2013_Iscap_02032020.xlsx")

#subset by adults
tick_adults <- tick_data[tick_data$Stage == "A",]

#aggregate tick count and tick infection count by county and year
agg_tick_data <- aggregate(list(tick_adults$Infection.Count, tick_adults$Count), by = list(tick_adults$Year, tick_adults$County), sum)

#rename colmns
names(agg_tick_data)[1] <- "Year"
names(agg_tick_data)[2] <- "County"
names(agg_tick_data)[3] <- "Infection.Count"
names(agg_tick_data)[4] <- "Count"

#create new column for tick infection percent by dividing infection count by tick count
agg_tick_data$infection_percent <- agg_tick_data$Infection.Count/agg_tick_data$Count

#remove all Infinte and NA values
agg_tick_data$infection_percent[!is.finite(agg_tick_data$infection_percent)] <- NA
agg_tick_data <- na.omit(agg_tick_data)

#create vector of county in IA sampled
county <- unique(agg_tick_data$County)
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
proportion.wrong.before.stability = c()
num.phases = c()

adult_tick_inf_count_vect = c()

sink(paste("iowa_adult_tick_infection_results_",(Sys.Date()),".txt", sep = ""))
print("Percent adult deer ticks infected with B. burgdoferi in Iowa counties")
for (i in county) {
  tryCatch({
    #location: i
    tick_loc <- agg_tick_data[agg_tick_data$County == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "infection_percent")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #data now cleaned
    
    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(ticks_found) >= 10) {
      
      #i returns location name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, ticks_found$Year[1])
      end.year <- c(end.year, ticks_found$Year[nrow(ticks_found)])
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("start year for ", i, ":",  ticks_found$Year[1]))
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("end year for ", i, ":", ticks_found$Year[nrow(ticks_found)]))
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(tryCatch(multiple_breakups(ticks_found), error=function(err) NA))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- tryCatch(stability_time(ticks_found), error=function(err) NA)
      print(st)
      stability.time <- c(stability.time, st)
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("stability time for ", i, ":", st))
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- tryCatch(abs_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("absolute range min for ", i, ":", ar[1]))
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("absolute range max for ", i, ":", ar[2]))
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- tryCatch(relative_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("relative range min for ", i, ":", rr[1]))
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("relative range max for ", i, ":", rr[2]))
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- tryCatch(proportion_significant(ticks_found, significance = 0.05), error=function(err) NA)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("proportion significant for ", i, ":", ps))
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- tryCatch(proportion_wrong(ticks_found, significance = 0.05), error=function(err) NA)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("proportion significantly wrong for ", i, ":", psw))
      
      print("proportion wrong before stability time")
      pwbst <- tryCatch(proportion_wrong_before_stability(ticks_found, significance = 0.05), error=function(err) NA)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("proportion significantly wrong before stability time for ", i, ":", pwbst))
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <-  bestfit(ticks_modeled, "AIC")
      num_phases <- tryCatch(bf$Nfits, error=function(err) NA)
      print(num_phases)
      print(bf)
      if(county.with.10.years.data[county.with.10.years.data==i] == i) {
        num.phases <- c(num.phases, num_phases)
        adult_tick_inf_count_vect <- c(adult_tick_inf_count_vect, paste("number of phases for ", i, ":", num_phases))
      }
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Percent adult deer ticks infected with B. burgdoferi in ", i, ", IA ", (Sys.Date()),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("Percent adult deer ticks infected with B. burgdoferi in ", i, ", IA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

write.csv(adult_tick_inf_count_vect, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_adult_tick_count_infection_text_results_",(Sys.Date()),".csv",sep = ""))

iowa_adult_tick_infection_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, proportion.wrong.before.stability, num.phases)
iowa_adult_tick_infection_data_results <- na.omit(iowa_adult_tick_infection_data_results)
write.csv(iowa_adult_tick_infection_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_adult_tick_infection_results_",(Sys.Date()),".csv",sep = ""))


#------------------------------------------------------------------------#
# Percent of nymph deer ticks infected with B. burgdoferi in Iowa counties
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from IA
tick_data <- read.xlsx("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/iowa_data/IA_LDS1990-2013_Iscap_02032020.xlsx")

#subset by nymphs
tick_nymphs <- tick_data[tick_data$Stage == "N",]

#aggregate tick count and tick infection count by county and year
agg_tick_data <- aggregate(list(tick_nymphs$Infection.Count, tick_nymphs$Count), by = list(tick_nymphs$Year, tick_nymphs$County), sum)

#rename colmns
names(agg_tick_data)[1] <- "Year"
names(agg_tick_data)[2] <- "County"
names(agg_tick_data)[3] <- "Infection.Count"
names(agg_tick_data)[4] <- "Count"

#create new column for tick infection percent by dividing infection count by tick count
agg_tick_data$infection_percent <- agg_tick_data$Infection.Count/agg_tick_data$Count

#remove all Infinte and NA values
agg_tick_data$infection_percent[!is.finite(agg_tick_data$infection_percent)] <- NA
agg_tick_data <- na.omit(agg_tick_data)

#create vector of county in IA sampled
county <- unique(agg_tick_data$County)
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
proportion.wrong.before.stability = c()
num.phases = c()

nymph_tick_inf_count_vect = c()

sink(paste("iowa_nymph_tick_infection_results_",(Sys.Date()),".txt", sep = ""))
print("Percent nymph deer ticks infected with B. burgdoferi in Iowa counties")
for (i in county) {
  tryCatch({
    #location: i
    tick_loc <- agg_tick_data[agg_tick_data$County == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "infection_percent")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #data now cleaned
    
    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(ticks_found) >= 10) {
      
      #i returns location name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, ticks_found$Year[1])
      end.year <- c(end.year, ticks_found$Year[nrow(ticks_found)])
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("start year for ", i, ":",  ticks_found$Year[1]))
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("end year for ", i, ":", ticks_found$Year[nrow(ticks_found)]))
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(tryCatch(multiple_breakups(ticks_found), error=function(err) NA))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- tryCatch(stability_time(ticks_found), error=function(err) NA)
      print(st)
      stability.time <- c(stability.time, st)
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("stability time for ", i, ":", st))
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- tryCatch(abs_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("absolute range min for ", i, ":", ar[1]))
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("absolute range max for ", i, ":", ar[2]))
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- tryCatch(relative_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("relative range min for ", i, ":", rr[1]))
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("relative range max for ", i, ":", rr[2]))
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- tryCatch(proportion_significant(ticks_found, significance = 0.05), error=function(err) NA)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("proportion significant for ", i, ":", ps))
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- tryCatch(proportion_wrong(ticks_found, significance = 0.05), error=function(err) NA)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("proportion significantly wrong for ", i, ":", psw))
      
      print("proportion wrong before stability time")
      pwbst <- tryCatch(proportion_wrong_before_stability(ticks_found, significance = 0.05), error=function(err) NA)
      print(pwbst)
      proportion.wrong.before.stability <- c(proportion.wrong.before.stability, pwbst)
      nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("proportion significantly wrong before stability time for ", i, ":", pwbst))
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <-  bestfit(ticks_modeled, "AIC")
      num_phases <- tryCatch(bf$Nfits, error=function(err) NA)
      print(num_phases)
      print(bf)
      if(county.with.10.years.data[county.with.10.years.data==i] == i) {
        num.phases <- c(num.phases, num_phases)
        nymph_tick_inf_count_vect <- c(nymph_tick_inf_count_vect, paste("number of phases for ", i, ":", num_phases))
      }
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Percent nymph deer ticks infected with B. burgdoferi in ", i, ", IA ", (Sys.Date()),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("Percent nymph deer ticks infected with B. burgdoferi in ", i, ", IA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

write.csv(nymph_tick_inf_count_vect, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_nymph_tick_count_infection_text_results_",(Sys.Date()),".csv",sep = ""))

iowa_nymph_tick_infection_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, proportion.wrong.before.stability, num.phases)
iowa_nymph_tick_infection_data_results <- na.omit(iowa_nymph_tick_infection_data_results)
write.csv(iowa_nymph_tick_infection_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_nymph_tick_infection_results_",(Sys.Date()),".csv",sep = ""))



#------------------------------------------------------------------------#
# Percent of larval deer ticks infected with B. burgdoferi in Iowa counties
# run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from IA
tick_data <- read.xlsx("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/iowa_data/IA_LDS1990-2013_Iscap_02032020.xlsx")

#subset by larvals
tick_larvals <- tick_data[tick_data$Stage == "L",]

#aggregate tick count and tick infection count by county and year
agg_tick_data <- aggregate(list(tick_larvals$Infection.Count, tick_larvals$Count), by = list(tick_larvals$Year, tick_larvals$County), sum)

#rename colmns
names(agg_tick_data)[1] <- "Year"
names(agg_tick_data)[2] <- "County"
names(agg_tick_data)[3] <- "Infection.Count"
names(agg_tick_data)[4] <- "Count"

#create new column for tick infection percent by dividing infection count by tick count
agg_tick_data$infection_percent <- agg_tick_data$Infection.Count/agg_tick_data$Count

#remove all Infinte and NA values
agg_tick_data$infection_percent[!is.finite(agg_tick_data$infection_percent)] <- NA
agg_tick_data <- na.omit(agg_tick_data)

#create vector of county in IA sampled
county <- unique(agg_tick_data$County)
county <- na.omit(county)
county <- as.character(county)

#set up vectors to collect results of stability time, absolute range, relative range 
#run on larval deer tick data for each county with at least 10 years data
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

larval_tick_inf_count_vect = c()

sink(paste("iowa_larval_tick_infection_results_",(Sys.Date()),".txt", sep = ""))
print("Percent larval deer ticks infected with B. burgdoferi in Iowa counties")
for (i in county) {
  tryCatch({
    #location: i
    tick_loc <- agg_tick_data[agg_tick_data$County == i,]
    #get only year and count column
    ticks_found <- tick_loc[,c("Year", "infection_percent")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #data now cleaned
    
    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(ticks_found) >= 10) {
      
      #i returns location name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste(i))
      larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste(i))
      
      #add start and end year to vector to be included in csv
      start.year <- c(start.year, ticks_found$Year[1])
      end.year <- c(end.year, ticks_found$Year[nrow(ticks_found)])
      larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste("start year for ", i, ":",  ticks_found$Year[1]))
      larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste("end year for ", i, ":", ticks_found$Year[nrow(ticks_found)]))
      
      #iterates through our targetted windows and returns summary statistics
      print("Multiple Breakups")
      print(tryCatch(multiple_breakups(ticks_found), error=function(err) NA))
      
      #returns how many years it takes to reach stability
      print("Stability Time")
      st <- tryCatch(stability_time(ticks_found), error=function(err) NA)
      print(st)
      stability.time <- c(stability.time, st)
      larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste("stability time for ", i, ":", st))
      
      #abs_range returns the absolute range of significant findings
      print("Absolute Range")
      ar <- tryCatch(abs_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste("absolute range min for ", i, ":", ar[1]))
      larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste("absolute range max for ", i, ":", ar[2]))
      
      #relative_range returns the absolute over and under estimate compared to the slope of the longest series
      print("Relative Range")
      rr <- tryCatch(relative_range(ticks_found, only_significant = FALSE, significance = 0.05), error=function(err) NA)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste("relative range min for ", i, ":", rr[1]))
      larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste("relative range max for ", i, ":", rr[2]))
      
      #returns the proportion of total windows with statistically significant values
      print("Proportion Significant")
      ps <- tryCatch(proportion_significant(ticks_found, significance = 0.05), error=function(err) NA)
      print(ps)
      proportion.significant <- c(proportion.significant, ps)
      larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste("proportion significant for ", i, ":", ps))
      
      #returns proportion of significant relationships that does not match the direction of the true slope
      print("Proportion Significantly Wrong")
      psw <- tryCatch(proportion_wrong(ticks_found, significance = 0.05), error=function(err) NA)
      print(psw)
      proportion.wrong <- c(proportion.wrong, psw)
      larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste("proportion significantly wrong for ", i, ":", psw))
      
      print("Number of phases")
      #add column for adundance next year to determine best model
      ticks_modeled <- addNt1(ticks_found)
      bf <-  bestfit(ticks_modeled, "AIC")
      num_phases <- tryCatch(bf$Nfits, error=function(err) NA)
      print(num_phases)
      print(bf)
      if(county.with.10.years.data[county.with.10.years.data==i] == i) {
        num.phases <- c(num.phases, num_phases)
        larval_tick_inf_count_vect <- c(larval_tick_inf_count_vect, paste("number of phases for ", i, ":", num_phases))
      }
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Percent larval deer ticks infected with B. burgdoferi in ", i, ", IA ", (Sys.Date()),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("Percent larval deer ticks infected with B. burgdoferi in ", i, ", IA", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
sink()

write.csv(larval_tick_inf_count_vect, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_larval_tick_count_infection_text_results_",(Sys.Date()),".csv",sep = ""))

iowa_larval_tick_infection_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
iowa_larval_tick_infection_data_results <- na.omit(iowa_larval_tick_infection_data_results)
write.csv(iowa_larval_tick_infection_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/iowa_larval_tick_infection_results_",(Sys.Date()),".csv",sep = ""))

