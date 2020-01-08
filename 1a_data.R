## data set ####################################################################

# read the file as a data frame
Data = gdata::read.xls("ebola_ng.xlsx")
# delete columns with NA
Data = Data[, colSums(is.na(Data)) != nrow(Data)]
# rename columns
colnames(Data) <- c("i", "vi", "tiE", "tiS", "tiH", "tiD", "tiR")

# convert into date from text
Data$tiE = as.Date(Data$tiE, "%Y-%m-%d")
Data$tiS = as.Date(Data$tiS, "%Y-%m-%d")
Data$tiH = as.Date(Data$tiH, "%Y-%m-%d")
Data$tiD = as.Date(Data$tiD, "%Y-%m-%d")
# convert into numeric from date
Data$tiE = as.numeric(Data$tiE - as.Date("2014-07-17"), units="days")
Data$tiS = as.numeric(Data$tiS - as.Date("2014-07-17"), units="days")
Data$tiH = as.numeric(Data$tiH - as.Date("2014-07-17"), units="days")
Data$tiD = as.numeric(Data$tiD - as.Date("2014-07-17"), units="days")
# the end of the outbreak on 20 Oct 2014 (=Day 95)
Data_tEND = as.numeric(as.Date("2014-10-20") - as.Date("2014-07-17"), units="days")

# existence of time points
Data_obs = matrix(!is.na(c(Data$tiE, Data$tiS, Data$tiH, Data$tiD)), ncol=4 )

# convert into numeric from factor
Data$vi = suppressWarnings( as.integer(as.character(Data$vi)) )
# Note that contacts with Case 20: wi = c(5, 6, 7, 10, 11, 12)
# now: Data$vi[20] = NA; and is selected from the candidate list during MCMC

## observed time lags ##########################################################

# incubation periods
Data$tauiES = Data$tiS - Data$tiE

# time lag from illness onset to hospitalization
Data$tauiSH = Data$tiH - Data$tiS

# serial intervals
Data$tauiSS[Data$i[Data$vi != 0 | is.na(Data$vi)]] = 
  Data$tiS[Data$i[Data$vi != 0 | is.na(Data$vi)]] - Data$tiS[Data$vi]

# take non-NA values
( tauiES = Data$tauiES[!is.na(Data$tauiES)] )
( tauiSH = Data$tauiSH[!is.na(Data$tauiSH)] )
( tauiSS = Data$tauiSS[!is.na(Data$tauiSS)] )