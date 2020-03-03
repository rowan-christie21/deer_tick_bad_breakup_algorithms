#--------------------------------------------------------------------------#

#call needed libaries
library("ggplot2")

#call to bad breakup algorithm developed by Dr. Christie Bahlai of Kent State
source("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/bad_breakup_2/R_model/bad_breakup_script.R")

#call to regime shift detector developed by Dr. Christie Bahlai of Kent State and Dr. Elise Zipkin
source("D:/Ixodes_scapularis_research_2019/regime_shift_tick_algorithms/tick_regime_09192019/monarch_regime-master_07162019/regime_shift_detector.R")

#-----------------------------------------------------------------------------------------------#
sink(paste("harvard_forest_deer_tick_bad_breakup_results",Sys.Date(),".txt", sep = ""))
#Oct 15 2019
# Cleaning up and modeling ticka found on people data from an oppurtunistic study in Harvard forest from here:
# https://portal.lternet.edu/nis/mapbrowse?packageid=knb-lter-hfr.299.2
# then running bad breakup, absolute range, relative range, stabilitiy time on it

##########################################
print("Deer ticks found on people in Harvard forest")
#now run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from Harvard forest
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/harvard_forest/hf299-01-survey.csv", head = T)

#create vector of locations in Harvard forest sampled
locations <- unique(tick_data$location.name)
locations <- na.omit(locations)
locations <- as.character(locations)

for (i in locations) {
  tryCatch({
    #location: i
    tick_loc <- tick_data[tick_data$location.name == i,]
    #look at deer ticks *found*
    #get only year and ticks found column
    ticks_found <- tick_loc[,c("year", "deer.found")]
    #omit nas
    ticks_found <- na.omit(ticks_found)
    #aggregate all counts per year in single row
    frame=data.frame(ticks_found)
    ticks_found <- aggregate(frame['deer.found'], by=frame['year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(ticks_found) >= 10) {
      print(i)
      #Multiple breakups
      print("Multiple Breakups")
      print(multiple_breakups(ticks_found))
      #Stability
      print("Stability Time")
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
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Deer tick found in ", i, ", Harvard Forest ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(ticks_found, title=paste("Deer ticks found in ", i, ", Harvard Forest", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

##########################################
#now run the aggregated count in Harvard forest as a whole for deer ticks found on people
print("Harvard Forest/n")

#read in survey data from Harvard forest
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/harvard_forest/hf299-01-survey.csv", head = T)

#location: Harvard forest
tick_loc <- tick_data
#look at deer ticks *found*
#get only year and ticks found column
ticks_found <- tick_loc[,c("year", "deer.found")]
#omit nas
ticks_found <- na.omit(ticks_found)
#aggregate all counts per year in single row
frame=data.frame(ticks_found)
ticks_found <- aggregate(frame['deer.found'], by=frame['year'], sum)
#data now cleaned for Simes Tract

#Multiple breakups
print("Multiple Breakups")
print(multiple_breakups(ticks_found))
#Stability
print("Stability")
print(stability_time(ticks_found))
#Absolute range
print("Absolute Range")
print(abs_range(ticks_found, only_significant = FALSE, significance = 0.05))
print("Relative Range")
print(relative_range(ticks_found, only_significant = FALSE, significance = 0.05))
print("Proportion Significant")
print(proportion_significant(ticks_found, significance = 0.05))
print("Proportion Significantly Wrong")
print(proportion_wrong(ticks_found, significance = 0.05))
#add column for adundance next year to determine best model
ticks_modeled <- addNt1(ticks_found)
bf <- bestfit(ticks_modeled, "AIC")
print(bf$Nfits)
print(bf)

png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Deer ticks found in Harvard Forest ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
pyramid_plot(ticks_found, title="Deer ticks found in Harvard Forest", plot_insig = TRUE, significance=0.05, rsq_points =TRUE)
dev.off()

#-----------------------------------------------------------------------------------------------#
#Oct 15 2019
# Cleaning up and modeling tick bite data from an oppurtunistic study in Harvard forest from here:
# https://portal.lternet.edu/nis/mapbrowse?packageid=knb-lter-hfr.299.2
# then running bad breakup, absolute range, relative range, stabilitiy time, proportion significant,
# proportion significantly wrong functions on it

print("Deer tick bites on people in Harvard forest")

##########################################
#now run bad breakup, stabiliity time, relative range, absolute range for all locations 

#read in survey data from Harvard forest
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/harvard_forest/hf299-01-survey.csv", head = T)

#create vector of locations in Harvard forest sampled
locations <- unique(tick_data$location.name)
locations <- na.omit(locations)
locations <- as.character(locations)

for (i in locations) {
  tryCatch({
    #location: i
    tick_loc <- tick_data[tick_data$location.name == i,]
    #look at deer ticks *found*
    #get only year and ticks found column
    tick_bites <- tick_loc[,c("year", "deer.bite")]
    #omit nas
    tick_bites <- na.omit(tick_bites)
    #aggregate all counts per year in single row
    frame=data.frame(tick_bites)
    tick_bites <- aggregate(frame['deer.bite'], by=frame['year'], sum)
    #data now cleaned
    
    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(tick_bites) >= 10) {
      print(i)
      #Multiple breakups
      print("Multiple Breakups")
      multiple_breakups(tick_bites)
      #Stability
      print("Stability Time")
      print(stability_time(tick_bites))
      #Absolute range
      print("Absolute Range")
      print(abs_range(tick_bites, only_significant = FALSE, significance = 0.05))
      print("Relative Range")
      print(relative_range(tick_bites, only_significant = FALSE, significance = 0.05))
      print("Proportion Significant")
      print(proportion_significant(ticks_found, significance = 0.05))
      print("Proportion Significantly Wrong")
      print(proportion_wrong(ticks_found, significance = 0.05))
      
      #print pyramid plot
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Deer tick bites in ", i, ", Harvard Forest ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(tick_bites, title=paste("Deer tick bites in ", i, ", Harvard Forest", sep = ''), plot_insig = TRUE, significance=0.05, rsq_points =TRUE))
      dev.off()
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

##########################################
#now run the aggregated count in Harvard forest as a whole for deer tick bites on people
print("Harvard Forest/n")

#read in survey data from Harvard forest
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/harvard_forest/hf299-01-survey.csv", head = T)

#location: Harvard forest
tick_loc <- tick_data
#look at deer ticks *found*
#get only year and ticks found column
tick_bites <- tick_loc[,c("year", "deer.bite")]
#omit nas
tick_bites <- na.omit(tick_bites)
#aggregate all counts per year in single row
frame=data.frame(tick_bites)
tick_bites <- aggregate(frame['deer.bite'], by=frame['year'], sum)

#Multiple breakups
print("Multiple Breakups")
print(multiple_breakups(tick_bites))
#Stability
print("Stability Time")
print(stability_time(tick_bites))
#Absolute range
print("Absolute Range")
print(abs_range(tick_bites, only_significant = FALSE, significance = 0.05))
#successful
print("Relative Range")
print(relative_range(tick_bites, only_significant = FALSE, significance = 0.05))
print("Proportion Significant")
print(proportion_significant(ticks_found, significance = 0.05))
print("Proportion Significantly Wrong")
print(proportion_wrong(ticks_found, significance = 0.05))
#add column for adundance next year to determine best model
ticks_modeled <- addNt1(ticks_found)
bf <- bestfit(ticks_modeled, "AIC")
print(bf$Nfits)
print(bf)

png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Deer tick bites in Harvard Forest ", Sys.Date(),".png", sep = ''), width = 876, height = 604)
pyramid_plot(tick_bites, title="Deer tick bites in Harvard Forest", plot_insig = TRUE, significance=0.05, rsq_points =TRUE)
dev.off()

sink()
