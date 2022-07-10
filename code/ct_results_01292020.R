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