#Oct 3rd 2019
#creating pyramid plots for deer tick data in NY counties and recording results of multiple breakup function
#dataset source: New York Deptartment of Health
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
source("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/christie_bahlai_bad_breakup/R_model/bad_breakup_script.R")

#-------------------------------------------------#
#Tick adults

#read in tick adult data
tick_adults_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/data/Deer_Tick_Surveillance__Adults__Oct_to_Dec__excluding_Powassan_virus__Beginning_2008.csv", head = T)

#create vector of counties in NY state in the dataset
counties <- unique(tick_adults_data$County)
counties <- na.omit(counties)
counties <- as.character(counties)

#function_results <- data.frame(county=c(""), stability.time=c(""), Absolute.Range.Min.Value=c(""), Absolute.Range.Max.Value=c(""), Relative.Range.Min.Value=c(""), Relative.Range.Max.Value=c(""), stringsAsFactors=FALSE)
#rbind(function_results,list("tree",43,3,12,43,55))

counties.with.10.years.data = c()
stability.time = c()
absolute.range.min.value = c()
absolute.range.max.value = c()
relative.range.min.value = c()
relative.range.max.value = c()

for (i in counties) {
  tryCatch({
    #location: i
    county <- tick_data[tick_adults_data$County == i,]
    #get only year and ticks found column
    tick_adults <- county[,c("Year", "Tick.Population.Density")]
    #omit nas
    tick_adults <- na.omit(tick_adults)
    #data now cleaned

    #if number of years in location exceeds 9 run bad breakup algorithm
    if(nrow(tick_adults) >= 10) {
      print(i)
      counties.with.10.years.data <- c(counties.with.10.years.data, i)
      print("Multiple Breakups")
      print(multiple_breakups(tick_adults))
      print("Stability Time")
      st <- stability_time(tick_adults, error_multiplyer = 0.5)
      print(st)
      stability.time <- c(stability.time, st)
      print("Absolute Range")
      ar <- abs_range(tick_adults, only_significant = FALSE, significance = 0.5)
      print(ar)
      absolute.range.min.value <- c(absolute.range.min.value, ar[1])
      absolute.range.max.value <- c(absolute.range.max.value, ar[2])
      print("Relative Range")
      rr <- relative_range(tick_adults, only_significant = FALSE, significance = 0.05)
      print(rr)
      relative.range.min.value <- c(relative.range.min.value, rr[1])
      relative.range.max.value <- c(relative.range.max.value, rr[2])
      
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/Adult deer tick density in ", i, ", NY ",Sys.Date(),".png", sep = ''), width = 876, height = 604)
      print(pyramid_plot(tick_adults, title=paste("Adult deer tick density in ", i, ", NY", sep = ''), plot_insig = TRUE, significance=0.1, rsq_points =TRUE, caption_plot = paste("Stability Time:", st)))
      dev.off()

    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

county_adult_tick_data_results <- data.frame(counties.with.10.years.data, stability.time, absolute.range.min.value, absolute.range.max.value, relative.range.min.value, relative.range.max.value)
write.csv(county_adult_tick_data_results, file = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_tick_algorithms/ny_county_adult_tick_data_results_",Sys.Date(),".csv",sep = ""))
