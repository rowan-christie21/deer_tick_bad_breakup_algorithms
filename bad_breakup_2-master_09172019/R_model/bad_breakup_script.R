# How does when you look, and how long, affect the conclusions you reach about your data?
# Are short term studies more likely to yield significant results?
# Are short term studies more likely to find *erroneous* significant trends?
# This script will perform a simple moving window analysis to answer these questions
# for long term data across a variety of domains- essentially, data gets subsetted, 
# we run a simple linear regression on each subset and record summary stats, trends

#assume data is coming in in the form year, response variable
#use this test data set to build stuff
test<-read.csv(file="test.csv", header=TRUE)

# set it so we get our decimal places rather than sci notation in our outputs
options(scipen=10)

# technical things that need to be done:

# data needs to be normalized in some way because we will be operating across domains 
# and likely with dramatically different values. Let's use a standardization approach so
# we don't give undue influence to outliers but also gets responses in the same magnitude
# something like x(standard)=(x-mean)/Stdev, the z score

standardize<-function(data){
  #operate on the second column in data, where our response variable is
  data$stand.response<-(data[,2]-mean(data[,2]))/sd(data[,2])
  #name the columns consistently so we don't get headaches later
  names(data)<-c("year", "response", "stand.response")
  #spit back our new data frame
  return(data)
}

#try it on test data
test1<-standardize(test)
#seems to be functioning

# next we need a function that runs a simple linear model of x=year, y=response variable

linefit<-function (data){
  #fit the model
  model<-lm(stand.response~year, data=data)
  #create a vector of relevant outputs. We want slope, error, P value
  output<-c(min(data$year), #year the analysis started on
            nrow(data), #number of data points the analysis includes
            length(unique(data$year)), #number of years the analysis includes
            summary(model)$coefficients[2,1], # slope
            summary(model)$coefficients[2,2], # se for slope
            summary(model)$coefficients[2,4], #p value
            summary(model)$r.squared, #r-squared
            summary(model)$adj.r.squared) #adjusted r-squared
  return(output)
}

#and try this on test data
linefit(test1)
# functional!

#now we need to think about how to iterate through the dataset. We want a
#function that starts at the first year, counts the number of rows specified
#and then feeds that rsultant data frame to the fittling function. 
#then we want to discard the first row of the data set, and repeat until fewer than
#the number of rows specified remains

breakup<-function(data, window){ #window is the size of the window we want to use
  remaining<-data #create dummy data set to operate on
  output<-data.frame(year=integer(0), #create empty data frame to put our output variables in
                     length=integer(0), 
                     years=integer(0),
                     slope=numeric(0), 
                     slope_se=numeric(0), 
                     p_value=numeric(0),
                     r_square=numeric(0),
                     adj_r_square=numeric(0))
  numyears<-length(unique(data$year))
  while (numyears>(window-1)){ #while there's still more years of data than in the window
    chunk<-subset(remaining, year<(min(year)+window)) #pull out a chunk as big as the window from the top of the data
    out<-linefit(chunk) #fit a linear model and get relevant statistics on chunk
    #add a conditional so that if there's missing data, it's not included in output
    if (window==length(unique(chunk$year))){
      output<-rbind(output, out) #append the stats to the output data frame
    }else{
      output<-output #leave it out if it has missing data
    }
    
    remaining<-subset(remaining, year>min(year)) #cut out the first year of the remaining data + repeat
    numyears<-length(unique(remaining$year))
  }
  names(output)<-c("start_year", "N_data", "N_years", "slope", "slope_se", "p_value",
                   "r_square", "adj_r_square")
  return(output)#output the data frame
}

#and now try this on the test data
breakup(test1, 3)

# now time to write the function that will iterate through our targetted windows
# let's make a decision rule that our test data set must be greater than 10y in length
# let's make this idiot-proof and build the standardize function right into this function
# so we can literally run each properly prepared raw data set with a single line

multiple_breakups<-function(data){
  data1<-standardize(data) #standardize data
  count<-length(data1$year)
  output<-data.frame(year=integer(0), #create empty data frame to put our output variables in
                     length=integer(0), 
                     years=integer(0),
                     slope=numeric(0), 
                     slope_se=numeric(0), 
                     p_value=numeric(0),
                     r_square=numeric(0),
                     adj_r_square=numeric(0))
  for(i in 3:(count-1)){
    outeach<-breakup(data1, i) #fit at each window length
    output<-rbind(output, outeach)#bind it to the frame
  }
  
  outall<-linefit(data1) #fit a line to the complete data set too
  out<-rbind(output, outall)
  return(out)
}

test2<-multiple_breakups(test)
#fan-flipping-tastic! it looks like that works


#let's create a plotting function
library(ggplot2)

pyramid_plot<- function(data, title="", significance=0.05, plot_insig=TRUE, rsq_points=FALSE){
  #suggested change
  colnames(data) <- c("year", "abundance")
  ################
  years<-length(unique(data$year))
  out<-multiple_breakups(data)
  count<-nrow(out)
  #compute mean and sd of longest series for vertical lines
  true_slope<-out[count,4] #find the slope of the longest series
  #remember to convert standard error to standard deviation
  true_error<-(out[count,5])*(sqrt(out[count, 2]))#find the error of the longest series
  max_true<-true_slope+true_error #compute max and min values for slopes we are calling true
  min_true<-true_slope-true_error
  out$significance<-ifelse(out$p_value<significance, "YES", "NO")
  if(rsq_points==TRUE){
    point_scale<-10*out$r_square
    yespt<-1
  }else{
    point_scale<-2
    yespt<-16
  }
  if(plot_insig==FALSE){
    out<-out[which(out$p_value<significance),]
  }
  plot<- ggplot(out) +
    theme_classic() +
    geom_hline(yintercept = true_slope, linetype = 2) +
    geom_hline(yintercept = max_true, linetype = 3, color="grey") +
    geom_hline(yintercept = min_true, linetype = 3, color="grey") +
    aes(y = slope, x = N_years,  ymin = (slope-slope_se), 
        ymax = (slope+slope_se), shape=significance, color=significance) +
    geom_linerange(show.legend = F)+ 
    geom_point(size=point_scale)+ ggtitle(title)+
    scale_shape_manual(values=c("NO"=4,"YES"=yespt))+
    scale_color_manual(values=c("NO"="red","YES"="black"))+
    xlab("Number of years in window")+xlim(3, years)+
    theme(plot.title = element_text(size=22))+
    theme(axis.title.x = element_text(size=17, face = "bold"))+
    theme(axis.title.y = element_text(size=17, face = "bold"))+
    theme(axis.text = element_text(size=14))+
    theme(legend.text = element_text(size = 12))+
    theme(legend.title = element_text(size = 14))+
    coord_flip()
  return(plot)
}

pyramid_plot(test, title="test plot", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)



#now that we have visualization, we need a way to pull relevant metrics out of the computation
#so let's say our longest series is our 'truth', and we want to know how many years it takes 
#to reach 'stability'-so let's define stability as >(some percentage of slopes) occuring within 
#the standard deviation of the slope of the longest series, for a given window length, allow user to change # of SEs

stability_time<-function(data, min_percent=95, error_multiplyer=1){#returns a number 
  test<-multiple_breakups(data)
  count<-nrow(test)
  true_slope<-test[count,4] #find the slope of the longest series
  #remember to convert standard error to standard deviation
  true_error<-(test[count,5])*(sqrt(test[count, 2]))*error_multiplyer #find the error of the longest series
  max_true<-true_slope+true_error #compute max and min values for slopes we are calling true
  min_true<-true_slope-true_error
  windows<-unique(test$N_years)#get a list of unique window lengths
  stability<-max(windows) #start with the assumption that the longest window is the only stable one
  for(i in 1:length(windows)){#for each window length, compute proportion 'correct'
    window_length<-windows[i]
    test_subset<-test[which(test$N_years==window_length),]
    number_of_windows<-nrow(test_subset)#how many windows
    correct_subset<-test_subset[which((test_subset$slope<max_true) & (test_subset$slope>min_true)),]
    number_of_correct<-nrow(correct_subset)#how many windows give the right answer
    percentage_correct<-100*number_of_correct/number_of_windows
    if(percentage_correct > min_percent){
      if(window_length < stability){
        stability<-window_length
      }
    }
  }
  return(stability)
}

#and a test
stability_time(test, error_multiplyer = 1.5)

#now a function that finds the absoloute range of findings, and the absolute 
#range of significant findings

abs_range<- function(data, only_significant=FALSE, significance=0.05){#returns a two unit vector with the max and min slopes
  test<-multiple_breakups(data)
  if(only_significant== TRUE){ #if user specifies only significant values wanted, pull those
    test1<-test[which(test$p_value<significance),]
  }else{
    test1<-test
  }
  max_slope<-max(test1$slope)
  min_slope<-min(test1$slope)
  sloperange<-c(min_slope, max_slope)
  return(sloperange)
  
}

#and try it out
abs_range(test, only_significant = FALSE, significance = 0.5)

#now we want to find the absolute over and under estimate compared to the slope of the 
#longest series

relative_range<- function(data, only_significant=FALSE, significance=0.05){#returns a two unit vector with the max and min slopes
  test<-multiple_breakups(data)
  count<-nrow(test)
  true_slope<-test[count,4] #find the slope of the longest series
  if(only_significant== TRUE){ #if user specifies only significant values wanted, pull those
    test1<-test[which(test$p_value<significance),]
  }else{
    test1<-test
  }
  max_slope<-max(test1$slope)-true_slope
  min_slope<-min(test1$slope)-true_slope
  sloperange<-c(min_slope, max_slope)
  return(sloperange)
  
}

relative_range(test, only_significant = FALSE, significance = 0.05)

#proportion significant- finds the proportion of total windows with statistically significant values

proportion_significant<- function(data, significance=0.05){#returns a single value between 0 and 1
  test<-multiple_breakups(data)
  count<-nrow(test)
  significant_regressions<-test[which(test$p_value<significance),]
  count_sig<-nrow(significant_regressions)
  proportion<-count_sig/count
  return(proportion)
  
}

proportion_significant(test, significance=0.05)

#proportion significantly wrong- we're going to define this as 'directionally wrong'
#where there is a significant relationship that does not match the direction of the true slope

proportion_wrong<- function(data, significance=0.05){#returns a single value between 0 and 1
  test<-multiple_breakups(data)
  count<-nrow(test)
  true_slope<-test[count,4] #find the slope of the longest series
  true_p<-test[count,6]
  #case 1: true slope is not significant
  if (true_p>significance){
    wrong_windows<-test[which(test$p_value<significance),]
  }else{ #true slope is significant
    if(true_slope>0){#true slope is positive
      wrong_windows<test[which(test$slope<0|test$p_value>significance),]#wrong means the slope is the wrong sign or 0
    }else{#true slope is negative
      wrong_windows<test[which(test$slope>0|test$p_value>significance),]#wrong means the slope is the wrong sign or 0
    }
  }
  count_wrong<-nrow(wrong_windows)
  proportion<-count_wrong/count
  return(proportion)
  
}

proportion_wrong(test, significance=0.01)



#proportion wrong by series length- basically the same thing as proportion wrong but looped 
#over all the unique window lengths. Will output a data frame with a window length and proportion
#of outputs are significantly misleading

proportion_wrong_series<- function(data, significance=0.05){#returns a single value between 0 and 1
  test<-multiple_breakups(data)
  count<-nrow(test)
  true_slope<-test[count,4] #find the slope of the longest series
  true_p<-test[count,6]
  windows<-unique(test$N_years)#get a list of unique window lengths
  prop.vec<-c()#create a blank vector to store proportions in
  for(i in 1:length(windows)){#for each window length, compute proportion 'wrong'
    window_length<-windows[i]
    test_subset<-test[which(test$N_years==window_length),]
    number_of_windows<-nrow(test_subset)#how many windows
    #case 1: true slope is not significant
    if (true_p>significance){
      wrong_windows<-test_subset[which(test_subset$p_value<significance),]
    }else{ #true slope is significant
      if(true_slope>0){#true slope is positive
        wrong_windows<test_subset[which(test_subset$slope<0|test_subset$p_value>significance),]#wrong means the slope is the wrong sign or 0
      }else{#true slope is negative
        wrong_windows<test_subset[which(test_subset$slope>0|test_subset$p_value>significance),]#wrong means the slope is the wrong sign or 0
      }
    }
    count_wrong<-nrow(wrong_windows)
    proportion<-count_wrong/number_of_windows
    prop.vec<-c(prop.vec, proportion)
  }
  
  x_name <- "window_length"
  y_name <- "proportion_wrong"
  
  df <- data.frame(windows,prop.vec)
  names(df) <- c(x_name,y_name)
  return(df)
  
}


#test it
proportion_wrong_series(test, significance = 0.1)

#########################################################################################

#now, let's give this a try with some real data

#lets' start with the lampyrid dataset because I'm familiar with it
#and know there is enough data to do it

#bring data in from figshare
lampyrid<-read.csv(file="https://ndownloader.figshare.com/files/3686040",
                   header=T)

#pull in our data cleaning code from https://github.com/cbahlai/lampyrid/blob/master/lampyrid_analysis.R
#details of cleaning in the code in comments found at that link- in summary, get all the typoes
#out and make the date column usable
library(lubridate)
lampyrid$newdate<-mdy(lampyrid$DATE)
lampyrid$year<-year(lampyrid$newdate)
lampyrid$DOY<-yday(lampyrid$newdate)
lampyrid<-na.omit(lampyrid)
lampyrid$TREAT_DESC<-gsub("Early succesional community", "Early successional", lampyrid$TREAT_DESC)
lampyrid$TREAT_DESC<-gsub("Early sucessional community", "Early successional", lampyrid$TREAT_DESC)
lampyrid$TREAT_DESC<-gsub("Succesional", "Successional", lampyrid$TREAT_DESC)
lampyrid$TREAT_DESC<-gsub("Sucessional", "Successional", lampyrid$TREAT_DESC)
lampyrid$TREAT_DESC<-gsub("Biologically based \\(organic\\)", "Organic", lampyrid$TREAT_DESC)
lampyrid$TREAT_DESC<-gsub("Conventional till", "Conventional", lampyrid$TREAT_DESC)
lampyrid$TREAT_DESC<-as.factor(lampyrid$TREAT_DESC)
lampyrid$HABITAT<-as.factor(lampyrid$HABITAT)
lampyrid$REPLICATE<-as.factor(lampyrid$REPLICATE)
lampyrid$STATION<-as.factor(lampyrid$STATION)

# data is currently structured as subsamples, replicates, across treatments. We know from 
# http://biorxiv.org/content/early/2016/09/11/074633 that the organisms were most abundant in
# Alfalfa and no-till treatments, so we should use data from these treatments if we're 
# interested in picking up trends over time. We also don't really care about within
# year dynamics for this experiment- we essentially want a summary measure of what was going on
#within a rep, within a year. So let's reshape the data, drop out irrelevant treatments.

library(reshape2)
#tell R where the data is by melting it, assigning IDs to the columns
lampyrid1<-melt(lampyrid, id=c("DATE","TREAT_DESC","HABITAT","REPLICATE","STATION","newdate", "year", "DOY"))
#cast the data to count up the fireflies
lampyrid2<-dcast(lampyrid1, year+TREAT_DESC+REPLICATE~., sum)
#cast the data to count the traps
lampyrid3<-dcast(lampyrid1, year+TREAT_DESC+REPLICATE~., length)
#let's rename these new vectors within the data frame
names(lampyrid2)[4]<-"ADULTS"
names(lampyrid3)[4]<-"TRAPS"

#rename the data frame and combine the number of traps we counted into it from lampyrid3
lampyrid_summary<-lampyrid2
lampyrid_summary$TRAPS<-lampyrid3$TRAPS

#create a new variable to account for trapping effort in a given year
lampyrid_summary$pertrap<-lampyrid_summary$ADULTS/lampyrid_summary$TRAPS

#get rid of columns we don't need for analysis
lampyrid_summary$REPLICATE<-NULL
lampyrid_summary$ADULTS<-NULL
lampyrid_summary$TRAPS<-NULL

# pull out relevant treatments, create a data frame for each
lampyrid_alfalfa<-subset(lampyrid_summary, TREAT_DESC=="Alfalfa")
lampyrid_notill<-subset(lampyrid_summary, TREAT_DESC=="No till")

#get rid of data that isn't needed for our analysis
lampyrid_alfalfa$TREAT_DESC<-NULL
lampyrid_notill$TREAT_DESC<-NULL

#ok, these data frames should be ready to go

#here goes nothing
multiple_breakups(lampyrid_alfalfa)
# there are some perculiarities because 2007 is missing, but I think we're working now.
# try it with other data too
multiple_breakups(lampyrid_notill)

#-------------------------------------------------------------------
#July 11
#working on excuting bad breakup data with selected tick data from
# https://datadryad.org/handle/10255/dryad.178685

library(ggplot2)

tick_data <- read.csv(file="data/Lyme_Study_Dataset.csv", header=TRUE)

#----------------------------------------------
#Grid GC

#create dataframe for gridcode GC
tick_grid <- tick_data[tick_data$Grid == "GC",]

#subset only tick adult and year
tick_adults <- tick_grid[,c("Year","DOA")]

#omit NAs
tick_adults <- na.omit(tick_adults)
#data successfully cleaned

#get decimal places
options(scipen=10)

#test standarization of data
tick_standarized_data <- standardize(tick_adults)
#successful

#test linear model
linefit(tick_standarized_data)
#successful

#testing breakup code on standarized data
breakup(tick_standarized_data,3)
#successful

#testing mulptiple breakup function on cleaned tick data
multiple_breakups(tick_adults)
#successful

#rename columns of dataframe so that it works with pyramid plot function
colnames(tick_adults) <- c("year", "abundance")

#try pyramid plot
png(filename = paste("Ixodes_scapularis_adults_in_Grid_GC",Sys.Date(),".png"))
pyramid_plot(tick_adults, title="Ixodes scapularis adults in Grid GC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#only shows 3 as the # of years in window, but after adding 
#the colnames line it works as intended

#test stability function
stability_time(tick_adults, error_multiplyer = 1.5)
#successful - # was 11

#test absolute range function
abs_range(tick_adults, only_significant = FALSE, significance = 0.5)

#do the same thing for tick nymphs and larvae data
tick_nymphs <- tick_grid[,c("Year", "DON")]
tick_nymphs <- na.omit(tick_nymphs)
#Multiple breakups
multiple_breakups(tick_nymphs)
png(filename = paste("Ixodes_scapularis_nymphs_in_Grid_GC",Sys.Date(),".png"))
pyramid_plot(tick_nymphs, title="Ixodes scapularis nymphs in Grid GC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

tick_larvae <- tick_grid[,c("Year", "DOL")]
tick_larvae <- na.omit(tick_larvae)
multiple_breakups(tick_larvae)
png(filename = paste("Ixodes_scapularis_larvae_in_Grid_GC",Sys.Date(),".png"))
pyramid_plot(tick_larvae, title="Ixodes scapularis larvae in Grid GC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

#do the same for lyme disease variable:
#Weighted average nymphal infection prevalence
#combining ticks tested using the two different
#methods when necessary
tick_infected <- tick_grid[,c("Year","NIP")]
tick_infected <- na.omit(tick_infected)
#Multiple breakups
multiple_breakups(tick_infected)
png(filename = paste("Fraction_of_nymphal_ticks_infected_with_Borrelia_burgdorferi_in_Grid_GC",Sys.Date(),".png"))
pyramid_plot(tick_infected, title="Fraction of nymphal ticks infected with Borrelia burgdorferi in Grid GC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_infected, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_infected, only_significant = FALSE, significance = 0.5)
#successful

#----------------------------------------------
#Grid GX

#create dataframe for gridcode GC
tick_grid <- tick_data[tick_data$Grid == "GX",]

#subset only tick adult and year
tick_adults <- tick_grid[,c("Year","DOA")]

#omit NAs
tick_adults <- na.omit(tick_adults)
#data successfully cleaned

#get decimal places
options(scipen=10)

#testing mulptiple breakup function on cleaned tick data
multiple_breakups(tick_adults)
#successful

#try pyramid plot
png(filename = paste("Ixodes_scapularis_adults_in_Grid_GX",Sys.Date(),".png"))
pyramid_plot(tick_adults, title="Ixodes scapularis adults in Grid GX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#only shows 3 as the # of years in window, but after adding 
#the colnames line it works as intended

#test stability function
stability_time(tick_adults, error_multiplyer = 1.5)
#successful - # was 11

#test absolute range function
abs_range(tick_adults, only_significant = FALSE, significance = 0.5)

#do the same thing for tick nymphs and larvae data
tick_nymphs <- tick_grid[,c("Year", "DON")]
tick_nymphs <- na.omit(tick_nymphs)
#Multiple breakups
multiple_breakups(tick_nymphs)
png(filename = paste("Ixodes_scapularis_nymphs_in_Grid_GX",Sys.Date(),".png"))
pyramid_plot(tick_nymphs, title="Ixodes scapularis nymphs in Grid GX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

tick_larvae <- tick_grid[,c("Year", "DOL")]
tick_larvae <- na.omit(tick_larvae)
multiple_breakups(tick_larvae)
png(filename = paste("Ixodes_scapularis_larvae_in_Grid_GX",Sys.Date(),".png"))
pyramid_plot(tick_larvae, title="Ixodes scapularis larvae in Grid GX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

#do the same for lyme disease variable:
#Weighted average nymphal infection prevalence
#combining ticks tested using the two different
#methods when necessary
tick_infected <- tick_grid[,c("Year","NIP")]
tick_infected <- na.omit(tick_infected)
#Multiple breakups
multiple_breakups(tick_infected)
png(filename = paste("Fraction_of_nymphal_ticks_infected_with_Borrelia_burgdorferi_in_Grid_GX",Sys.Date(),".png"))
pyramid_plot(tick_infected, title="Fraction of nymphal ticks infected with Borrelia burgdorferi in Grid GX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_infected, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_infected, only_significant = FALSE, significance = 0.5)
#successful

#----------------------------------------------
#Grid HC

#create dataframe for gridcode GC
tick_grid <- tick_data[tick_data$Grid == "HC",]

#subset only tick adult and year
tick_adults <- tick_grid[,c("Year","DOA")]

#omit NAs
tick_adults <- na.omit(tick_adults)
#data successfully cleaned

#get decimal places
options(scipen=10)

#testing mulptiple breakup function on cleaned tick data
multiple_breakups(tick_adults)
#successful

#try pyramid plot
png(filename = paste("Ixodes_scapularis_adults_in_Grid_HC",Sys.Date(),".png"))
pyramid_plot(tick_adults, title="Ixodes scapularis adults in Grid HC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#only shows 3 as the # of years in window, but after adding 
#the colnames line it works as intended

#test stability function
stability_time(tick_adults, error_multiplyer = 1.5)
#successful - # was 11

#test absolute range function
abs_range(tick_adults, only_significant = FALSE, significance = 0.5)

#do the same thing for tick nymphs and larvae data
tick_nymphs <- tick_grid[,c("Year", "DON")]
tick_nymphs <- na.omit(tick_nymphs)
#Multiple breakups
multiple_breakups(tick_nymphs)
png(filename = paste("Ixodes_scapularis_nymphs_in_Grid_HC",Sys.Date(),".png"))
pyramid_plot(tick_nymphs, title="Ixodes scapularis nymphs in Grid HC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

tick_larvae <- tick_grid[,c("Year", "DOL")]
tick_larvae <- na.omit(tick_larvae)
multiple_breakups(tick_larvae)
png(filename = paste("Ixodes_scapularis_larvae_in_Grid_HC",Sys.Date(),".png"))
pyramid_plot(tick_larvae, title="Ixodes scapularis larvae in Grid HC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

#do the same for lyme disease variable:
#Weighted average nymphal infection prevalence
#combining ticks tested using the two different
#methods when necessary
tick_infected <- tick_grid[,c("Year","NIP")]
tick_infected <- na.omit(tick_infected)
#Multiple breakups
multiple_breakups(tick_infected)
png(filename = paste("Fraction_of_nymphal_ticks_infected_with_Borrelia_burgdorferi_in_Grid_HC",Sys.Date(),".png"))
pyramid_plot(tick_infected, title="Fraction of nymphal ticks infected with Borrelia burgdorferi in Grid HC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_infected, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_infected, only_significant = FALSE, significance = 0.5)
#successful

#----------------------------------------------
#Grid HX

#create dataframe for gridcode GC
tick_grid <- tick_data[tick_data$Grid == "HX",]

#subset only tick adult and year
tick_adults <- tick_grid[,c("Year","DOA")]

#omit NAs
tick_adults <- na.omit(tick_adults)
#data successfully cleaned

#get decimal places
options(scipen=10)

#testing mulptiple breakup function on cleaned tick data
multiple_breakups(tick_adults)
#successful

#try pyramid plot
png(filename = paste("Ixodes_scapularis_adults_in_Grid_HX",Sys.Date(),".png"))
pyramid_plot(tick_adults, title="Ixodes scapularis adults in Grid HX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#only shows 3 as the # of years in window, but after adding 
#the colnames line it works as intended

#test stability function
stability_time(tick_adults, error_multiplyer = 1.5)
#successful - # was 11

#test absolute range function
abs_range(tick_adults, only_significant = FALSE, significance = 0.5)

#do the same thing for tick nymphs and larvae data
tick_nymphs <- tick_grid[,c("Year", "DON")]
tick_nymphs <- na.omit(tick_nymphs)
#Multiple breakups
multiple_breakups(tick_nymphs)
png(filename = paste("Ixodes_scapularis_nymphs_in_Grid_HX",Sys.Date(),".png"))
pyramid_plot(tick_nymphs, title="Ixodes scapularis nymphs in Grid HX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

tick_larvae <- tick_grid[,c("Year", "DOL")]
tick_larvae <- na.omit(tick_larvae)
multiple_breakups(tick_larvae)
png(filename = paste("Ixodes_scapularis_larvae_in_Grid_HX",Sys.Date(),".png"))
pyramid_plot(tick_larvae, title="Ixodes scapularis larvae in Grid HX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

#do the same for lyme disease variable:
#Weighted average nymphal infection prevalence
#combining ticks tested using the two different
#methods when necessary
tick_infected <- tick_grid[,c("Year","NIP")]
tick_infected <- na.omit(tick_infected)
#Multiple breakups
multiple_breakups(tick_infected)
png(filename = paste("Fraction_of_nymphal_ticks_infected_with_Borrelia_burgdorferi_in_Grid_HX",Sys.Date(),".png"))
pyramid_plot(tick_infected, title="Fraction of nymphal ticks infected with Borrelia burgdorferi in Grid HX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_infected, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_infected, only_significant = FALSE, significance = 0.5)
#successful

#----------------------------------------------
#Grid TC

#create dataframe for gridcode GC
tick_grid <- tick_data[tick_data$Grid == "TC",]

#subset only tick adult and year
tick_adults <- tick_grid[,c("Year","DOA")]

#omit NAs
tick_adults <- na.omit(tick_adults)
#data successfully cleaned

#get decimal places
options(scipen=10)

#testing mulptiple breakup function on cleaned tick data
multiple_breakups(tick_adults)
#successful

#try pyramid plot
png(filename = paste("Ixodes_scapularis_adults_in_Grid_TC",Sys.Date(),".png"))
pyramid_plot(tick_adults, title="Ixodes scapularis adults in Grid TC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#only shows 3 as the # of years in window, but after adding 
#the colnames line it works as intended

#test stability function
stability_time(tick_adults, error_multiplyer = 1.5)
#successful - # was 11

#test absolute range function
abs_range(tick_adults, only_significant = FALSE, significance = 0.5)

#do the same thing for tick nymphs and larvae data
tick_nymphs <- tick_grid[,c("Year", "DON")]
tick_nymphs <- na.omit(tick_nymphs)
#Multiple breakups
multiple_breakups(tick_nymphs)
png(filename = paste("Ixodes_scapularis_nymphs_in_Grid_TC",Sys.Date(),".png"))
pyramid_plot(tick_nymphs, title="Ixodes scapularis nymphs in Grid TC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

tick_larvae <- tick_grid[,c("Year", "DOL")]
tick_larvae <- na.omit(tick_larvae)
multiple_breakups(tick_larvae)
png(filename = paste("Ixodes_scapularis_larvae_in_Grid_TC",Sys.Date(),".png"))
pyramid_plot(tick_larvae, title="Ixodes scapularis larvae in Grid TC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

#do the same for lyme disease variable:
#Weighted average nymphal infection prevalence
#combining ticks tested using the two different
#methods when necessary
tick_infected <- tick_grid[,c("Year","NIP")]
tick_infected <- na.omit(tick_infected)
#Multiple breakups
multiple_breakups(tick_infected)
png(filename = paste("Fraction_of_nymphal_ticks_infected_with_Borrelia_burgdorferi_in_Grid_TC",Sys.Date(),".png"))
pyramid_plot(tick_infected, title="Fraction of nymphal ticks infected with Borrelia burgdorferi in Grid TC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_infected, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_infected, only_significant = FALSE, significance = 0.5)
#successful

#----------------------------------------------
#Grid TX

#create dataframe for gridcode GC
tick_grid <- tick_data[tick_data$Grid == "TX",]

#subset only tick adult and year
tick_adults <- tick_grid[,c("Year","DOA")]

#omit NAs
tick_adults <- na.omit(tick_adults)
#data successfully cleaned

#get decimal places
options(scipen=10)

#testing mulptiple breakup function on cleaned tick data
multiple_breakups(tick_adults)
#successful

#try pyramid plot
png(filename = paste("Ixodes_scapularis_adults_in_Grid_TX",Sys.Date(),".png"))
pyramid_plot(tick_adults, title="Ixodes scapularis adults in Grid TX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#only shows 3 as the # of years in window, but after adding 
#the colnames line it works as intended

#test stability function
stability_time(tick_adults, error_multiplyer = 1.5)
#successful - # was 11

#test absolute range function
abs_range(tick_adults, only_significant = FALSE, significance = 0.5)

#do the same thing for tick nymphs and larvae data
tick_nymphs <- tick_grid[,c("Year", "DON")]
tick_nymphs <- na.omit(tick_nymphs)
#Multiple breakups
multiple_breakups(tick_nymphs)
png(filename = paste("Ixodes_scapularis_nymphs_in_Grid_TX",Sys.Date(),".png"))
pyramid_plot(tick_nymphs, title="Ixodes scapularis nymphs in Grid TX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

tick_larvae <- tick_grid[,c("Year", "DOL")]
tick_larvae <- na.omit(tick_larvae)
multiple_breakups(tick_larvae)
png(filename = paste("Ixodes_scapularis_larvae_in_Grid_TX",Sys.Date(),".png"))
pyramid_plot(tick_larvae, title="Ixodes scapularis larvae in Grid TX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

#do the same for lyme disease variable:
#Weighted average nymphal infection prevalence
#combining ticks tested using the two different
#methods when necessary
tick_infected <- tick_grid[,c("Year","NIP")]
tick_infected <- na.omit(tick_infected)
#Multiple breakups
multiple_breakups(tick_infected)
png(filename = paste("Fraction_of_nymphal_ticks_infected_with_Borrelia_burgdorferi_in_Grid_TX",Sys.Date(),".png"))
pyramid_plot(tick_infected, title="Fraction of nymphal ticks infected with Borrelia burgdorferi in Grid TX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_infected, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_infected, only_significant = FALSE, significance = 0.5)
#successful

#----------------------------------------------
#All Grids

#subset only tick adult and year
tick_adults <- tick_data[,c("Year","DOA")]

#omit NAs
tick_adults <- na.omit(tick_adults)
#data successfully cleaned

#get decimal places
options(scipen=10)

#testing mulptiple breakup function on cleaned tick data
multiple_breakups(tick_adults)
#successful

#try pyramid plot
png(filename = paste("Ixodes_scapularis_adults_in_all_Grids",Sys.Date(),".png"))
pyramid_plot(tick_adults, title="Adult deer ticks in Cary Forest", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#only shows 3 as the # of years in window, but after adding 
#the colnames line it works as intended

#test stability function
stability_time(tick_adults, error_multiplyer = 1.5)
#successful - # was 11

#test absolute range function
abs_range(tick_adults, only_significant = FALSE, significance = 0.5)

#do the same thing for tick nymphs and larvae data
tick_nymphs <- tick_grid[,c("Year", "DON")]
tick_nymphs <- na.omit(tick_nymphs)
#Multiple breakups
multiple_breakups(tick_nymphs)
png(filename = paste("Ixodes_scapularis_nymphs_in_all_Grids",Sys.Date(),".png"))
pyramid_plot(tick_nymphs, title="Nymphal deer ticks in Cary Forest", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

tick_larvae <- tick_grid[,c("Year", "DOL")]
tick_larvae <- na.omit(tick_larvae)
multiple_breakups(tick_larvae)
png(filename = paste("Ixodes_scapularis_larvae_in_all_Grids",Sys.Date(),".png"))
pyramid_plot(tick_larvae, title="Ixodes scapularis larvae in all Grids", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

#do the same for lyme disease variable:
#Weighted average nymphal infection prevalence
#combining ticks tested using the two different
#methods when necessary
tick_infected <- tick_grid[,c("Year","NIP")]
tick_infected <- na.omit(tick_infected)
#Multiple breakups
multiple_breakups(tick_infected)
png(filename = paste("Fraction_of_nymphal_ticks_infected_with_Borrelia_burgdorferi_in_all_Grids",Sys.Date(),".png"))
pyramid_plot(tick_infected, title="Fraction of nymphal ticks infected with Borrelia burgdorferi in all Grids", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_infected, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_infected, only_significant = FALSE, significance = 0.5)
#successful

#-------------------------------------------------
#July 16 2019
#working on creating pyramid plots for deer data
#variable used: Ave.deer.hr
#Total number of different white-tailed deer seen by bowhunters 
#in the Cary Institute deer management program divided by the total 
#number of hours afield during good light conditions for identifying deer

library(ggplot2)

#---------------------------
#Grid GC

tick_grid <- tick_data[tick_data$Grid == "GC",]
deer <- tick_grid[,c("Year", "Ave.Deer.Hr")]
deer <- na.omit(deer)

#Multiple breakups
multiple_breakups(deer)
png(filename = paste("Average_#_deer_per_hour_of_observation_in_Grid GC",Sys.Date(),".png"))
pyramid_plot(deer, title="Average # deer / hour of observation in Grid GC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(deer, error_multiplyer = 1.5)
#Absolute range
abs_range(deer, only_significant = FALSE, significance = 0.5)
#successful

#---------------------------
#Grid GX

tick_grid <- tick_data[tick_data$Grid == "GX",]
deer <- tick_grid[,c("Year", "Ave.Deer.Hr")]
deer <- na.omit(deer)

#Multiple breakups
multiple_breakups(deer)
png(filename = paste("Average_#_deer_per_hour_of_observation_in_Grid GX",Sys.Date(),".png"))
pyramid_plot(deer, title="Average # deer / hour of observation in Grid GX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(deer, error_multiplyer = 1.5)
#Absolute range
abs_range(deer, only_significant = FALSE, significance = 0.5)
#successful

#---------------------------
#Grid HC

tick_grid <- tick_data[tick_data$Grid == "HC",]
deer <- tick_grid[,c("Year", "Ave.Deer.Hr")]
deer <- na.omit(deer)

#Multiple breakups
multiple_breakups(deer)
png(filename = paste("Average_#_deer_per_hour_of_observation_in_Grid HC",Sys.Date(),".png"))
pyramid_plot(deer, title="Average # deer / hour of observation in Grid HC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(deer, error_multiplyer = 1.5)
#Absolute range
abs_range(deer, only_significant = FALSE, significance = 0.5)
#successful

#---------------------------
#Grid HX

tick_grid <- tick_data[tick_data$Grid == "HX",]
deer <- tick_grid[,c("Year", "Ave.Deer.Hr")]
deer <- na.omit(deer)

#Multiple breakups
multiple_breakups(deer)
png(filename = paste("Average_#_deer_per_hour_of_observation_in_Grid HX",Sys.Date(),".png"))
pyramid_plot(deer, title="Average # deer / hour of observation in Grid HX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(deer, error_multiplyer = 1.5)
#Absolute range
abs_range(deer, only_significant = FALSE, significance = 0.5)
#successful

#---------------------------
#Grid TC

tick_grid <- tick_data[tick_data$Grid == "TC",]
deer <- tick_grid[,c("Year", "Ave.Deer.Hr")]
deer <- na.omit(deer)

#Multiple breakups
multiple_breakups(deer)
png(filename = paste("Average_#_deer_per_hour_of_observation_in_Grid TC",Sys.Date(),".png"))
pyramid_plot(deer, title="Average # deer / hour of observation in Grid TC", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(deer, error_multiplyer = 1.5)
#Absolute range
abs_range(deer, only_significant = FALSE, significance = 0.5)
#successful

#---------------------------
#Grid TX

tick_grid <- tick_data[tick_data$Grid == "TX",]
deer <- tick_grid[,c("Year", "Ave.Deer.Hr")]
deer <- na.omit(deer)

#Multiple breakups
multiple_breakups(deer)
png(filename = paste("Average_#_deer_per_hour_of_observation_in_Grid TX",Sys.Date(),".png"))
pyramid_plot(deer, title="Average # deer / hour of observation in Grid TX", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(deer, error_multiplyer = 1.5)
#Absolute range
abs_range(deer, only_significant = FALSE, significance = 0.5)
#successful

#---------------------------
#All Grids

deer <- tick_data[,c("Year", "Ave.Deer.Hr")]
deer <- na.omit(deer)

#Multiple breakups
multiple_breakups(deer)
png(filename = paste("Average_#_deer_per_hour_of_observation_in_all_Grids",Sys.Date(),".png"))
pyramid_plot(deer, title="Average # deer / hour of observation in all Grids", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(deer, error_multiplyer = 1.5)
#Absolute range
abs_range(deer, only_significant = FALSE, significance = 0.5)
#successful

#--------------------------------------------------------------
#July 22 2019
#Create pyramid plots for NY Dept of Health tick data
#Data for tick adults from here: https://health.data.ny.gov/Health/Deer-Tick-Surveillance-Adults-Oct-to-Dec-excluding/vzbp-i2d4
#Data for tick nymphs from here: https://health.data.ny.gov/Health/Deer-Tick-Surveillance-Nymphs-May-to-Sept-excludin/kibp-u2ip

#call to tick adult data
tick_data <- read.csv(file="D:/Ixodes_scapularis_research_2019/tick_regime_07172019/data/Deer_Tick_Surveillance__Adults__Oct_to_Dec__excluding_Powassan_virus__Beginning_2008.csv", header=TRUE)

#create vector of county names
counties <- c("Albany",
              "Cattaraugus",
              "Chautauqua",
              "Chemung",
              "Clinton",
              "Columbia",
              "Dutchess",
              "Erie",
              "Herkimer",
              "Jefferson",
              "Monroe",
              "Onondaga",
              "Orange",
              "Oswego",
              "Rockland",
              "Saratoga",
              "Schoharie",
              "Schuyler",
              "Seneca",
              "Sullivan",
              "Tompkins",
              "Ulster",
              "Washington",
              "Westchester")

tick_county <- tick_data[tick_data$County == counties[3],]
tick_adults <- tick_county[,c("Year", "Tick.Population.Density")]
tick_adults <- na.omit(tick_adults)

#Multiple breakups
multiple_breakups(tick_adults)
#png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_07222019/R_model/Ixodes scapularis adults in Dutchess county",Sys.Date(),".png"))
pyramid_plot(tick_adults, title="Ixodes scapularis adults in Dutchess county", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
#dev.off()
#Stability
stability_time(tick_adults, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_adults, only_significant = FALSE, significance = 0.5)
#successful

#loop through all counties
for(i in counties) {
  tryCatch({
    tick_county <- tick_data[tick_data$County == i,]
    tick_adults <- tick_county[,c("Year", "Tick.Population.Density")]
    tick_adults <- na.omit(tick_adults)
    #Multiple breakups
    multiple_breakups(tick_adults)
    png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_07222019/R_model/Ixodes scapularis adults in",i,"county",Sys.Date(),".png"))
    print(pyramid_plot(tick_adults, title=paste("Adult deer ticks in",i,"county"), plot_insig = TRUE, significance=0.1, rsq_points =TRUE))
    dev.off()
    #Stability
    stability_time(tick_adults, error_multiplyer = 1.5)
    #Absolute range
    abs_range(tick_adults, only_significant = FALSE, significance = 0.5)
    #successful
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

##########################################
#do the same for nymphs
#call to tick adult data
tick_data <- read.csv(file="D:/Ixodes_scapularis_research_2019/tick_regime_07172019/data/Deer_Tick_Surveillance__Nymphs__May_to_Sept__excluding_Powassan_virus__Beginning_2008.csv", header=TRUE)

#create vector of county names
counties <- c("Albany",
              "Cattaraugus",
              "Chautauqua",
              "Chemung",
              "Clinton",
              "Columbia",
              "Dutchess",
              "Erie",
              "Herkimer",
              "Jefferson",
              "Monroe",
              "Onondaga",
              "Orange",
              "Oswego",
              "Rockland",
              "Saratoga",
              "Schoharie",
              "Schuyler",
              "Seneca",
              "Sullivan",
              "Tompkins",
              "Ulster",
              "Washington",
              "Westchester")

#loop through all counties
for(i in counties) {
  tryCatch({
    tick_county <- tick_data[tick_data$County == i,]
    tick_nymphs <- tick_county[,c("Year", "Tick.Population.Density")]
    tick_nymphs <- na.omit(tick_nymphs)
    #Multiple breakups
    multiple_breakups(tick_nymphs)
    png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_07222019/R_model/Ixodes scapularis nymphs in",i,"county",Sys.Date(),".png"))
    print(pyramid_plot(tick_nymphs, title=paste("Nymphal deer ticks in",i,"county"), plot_insig = TRUE, significance=0.1, rsq_points =TRUE))
    dev.off()
    #Stability
    stability_time(tick_nymphs, error_multiplyer = 1.5)
    #Absolute range
    abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
    #successful
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#------------------------------------------------
#july 25th 2019
#Create pyramid plot of NY with data from all counties aggregated

###################################
#Start with adults
tick_data <- read.csv(file="D:/Ixodes_scapularis_research_2019/tick_regime_07172019/data/Deer_Tick_Surveillance__Adults__Oct_to_Dec__excluding_Powassan_virus__Beginning_2008.csv", header=TRUE)

tick_county <- tick_data
tick_adults <- tick_county[,c("Year", "Tick.Population.Density")]

frame=data.frame(tick_adults)
tick_adults <- aggregate(frame['Tick.Population.Density'], by=frame['Year'], sum)

tick_adults <- na.omit(tick_adults)

#Multiple breakups
multiple_breakups(tick_adults)
png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_07222019/R_model/Ixodes scapularis adults in NY State",Sys.Date(),".png"))
pyramid_plot(tick_adults, title="Ixodes scapularis adults in NY State", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_adults, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_adults, only_significant = FALSE, significance = 0.5)
#successful

###################################
#now do nymphs
tick_data <- read.csv(file="D:/Ixodes_scapularis_research_2019/tick_regime_07172019/data/Deer_Tick_Surveillance__Nymphs__May_to_Sept__excluding_Powassan_virus__Beginning_2008.csv", header=TRUE)

tick_county <- tick_data
tick_nymphs <- tick_county[,c("Year", "Tick.Population.Density")]

frame=data.frame(tick_nymphs)
tick_nymphs <- aggregate(frame['Tick.Population.Density'], by=frame['Year'], sum)

tick_nymphs <- na.omit(tick_nymphs)

#Multiple breakups
multiple_breakups(tick_nymphs)
png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_07222019/R_model/Ixodes scapularis nymphs in NY State",Sys.Date(),".png"))
pyramid_plot(tick_nymphs, title="Ixodes scapularis nymphs in NY State", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(tick_nymphs, error_multiplyer = 1.5)
#Absolute range
abs_range(tick_nymphs, only_significant = FALSE, significance = 0.5)
#successful

#-----------------------------------------------------------------------------------------------#
#Sept 17 2019
#Cleaning up and modeling tick data from an oppurtunistic study in Harvard forest from here:
# https://portal.lternet.edu/nis/mapbrowse?packageid=knb-lter-hfr.299.2

#read in survey data from Harvard forest
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_09172019/R_model/data/harvard_forest/hf299-01-survey.csv", head = T)

#location: Simes Tract
tick_loc <- tick_data[tick_data$location.name == "Simes Tract",]
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
multiple_breakups(ticks_found)
#png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_07222019/R_model/Ixodes scapularis nymphs in NY State",Sys.Date(),".png"))
pyramid_plot(ticks_found, title="Deer ticks in Simes Tract, Harvard Forest", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
#dev.off()
#Stability
stability_time(ticks_found, error_multiplyer = 1.5)
#Absolute range
abs_range(ticks_found, only_significant = FALSE, significance = 0.5)
#successful

##########################################
#now run for all locations

#read in survey data from Harvard forest
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_09172019/R_model/data/harvard_forest/hf299-01-survey.csv", head = T)

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
      #Multiple breakups
      multiple_breakups(ticks_found)
      png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_09172019/R_model/Deer ticks in ", i, ", Harvard Forest ", Sys.Date(),".png", sep = ''))
      print(pyramid_plot(ticks_found, title=paste("Deer ticks in ", i, ", Harvard Forest", sep = ''), plot_insig = TRUE, significance=0.1, rsq_points =TRUE))
      dev.off()
      #Stability
      stability_time(ticks_found, error_multiplyer = 1.5)
      #Absolute range
      abs_range(ticks_found, only_significant = FALSE, significance = 0.5)
      #successful
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

##########################################
#now run the aggregated count in Harvard forest as a whole

#read in survey data from Harvard forest
tick_data <- read.csv("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_09172019/R_model/data/harvard_forest/hf299-01-survey.csv", head = T)

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
multiple_breakups(ticks_found)
png(filename = paste("D:/Ixodes_scapularis_research_2019/bad_breakup_2-master_09172019/R_model/Deer ticks in Harvard Forest ", Sys.Date(),".png", sep = ''))
pyramid_plot(ticks_found, title="Deer ticks in Harvard Forest", plot_insig = TRUE, significance=0.1, rsq_points =TRUE)
dev.off()
#Stability
stability_time(ticks_found, error_multiplyer = 1.5)
#Absolute range
abs_range(ticks_found, only_significant = FALSE, significance = 0.5)
#successful
relative_range(ticks_found, only_significant = FALSE, significance = 0.05)

