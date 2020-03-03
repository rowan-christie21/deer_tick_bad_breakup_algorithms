#Jan 22 2020
#creating pyramid plots for deer tick data in PA and recording results of multiple breakup, 
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

#call to bad breakup algorithm developed by Dr. Christie Bahlai of Kent State
source("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/bad_breakup_2/R_model/bad_breakup_script.R")

#call to regime shift detector developed by Dr. Christie Bahlai of Kent State and Dr. Elise Zipkin
source("D:/Ixodes_scapularis_research_2019/regime_shift_tick_algorithms/tick_regime_09192019/monarch_regime-master_07162019/regime_shift_detector.R")

#------------------------------------------------------------------------#
#Tick adult density

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

#location: i
county <- tick_data[tick_data$County.where.tick.found == "Delaware",]
#get only year and ticks found column
tick_adults <- county[,c("Year", "Adults")]
#omit nas
tick_adults <- na.omit(tick_adults)
#aggregate all counts per year in single row
frame=data.frame(tick_adults)
tick_adults <- aggregate(frame['Adults'], by=frame['Year'], sum)
#data now cleaned

start.year <- c(start.year, tick_adults$Year[1])
end.year <- c(end.year, tick_adults$Year[nrow(tick_adults)])

print(multiple_breakups(tick_adults))

#sink(paste("pa_adult_tick_density_data_results",Sys.Date(),".txt", sep = ""))
for (i in county) {
  tryCatch({
    #location: i
    county <- tick_data[tick_data$County.where.tick.found == i,]
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
      #i returns grid name
      print(i)
      county.with.10.years.data <- c(county.with.10.years.data, paste("County where tick found", i))
      
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

county_adult_tick_data_results <- data.frame(county.with.10.years.data, start.year, end.year, stability.time)#, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value, proportion.significant, proportion.wrong, num.phases)
write.csv(county_adult_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/pa_adult_tick_data_density_results_",Sys.Date(),".csv",sep = ""))
