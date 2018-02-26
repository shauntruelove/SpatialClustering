

# library(foreign) 
# library(SpatialEpi)
# library(maptools)
# library(maps)

library(rgdal)
library(ggplot2)
#library(dplyr)
#library(stringi)
library(sp)



####################################################################################
############  ------------  GET DHS DATA     --------------   ###
#Function gets and cleans DHS data
#
#Parameters -
# filename - name of STATA file with the DHS data
# data.year - year of DHS survey
# numvacc.dpt - number of different vaccination variables that can include DPT in that country (i.e. h3, sdhh1,...)
# country - specific country selected
# youngest.age - the youngest age to look at in months...should always be 0
# oldest.age - the oldest age to look at....default to 60 months
#
#Returns -
# data.frame with all data from DHS survey for selected country and year
# Includes: - "censoring" variable (same as dhs.to.vac.cens function), 
#           - recoded DPT and Measles doses and dates

dhs.data.import <- function(filename,
                            youngest.age = 0, oldest.age = 60, country,
                            model.type = "all", data.year=NA) {
  
  #library(foreign)
  `%ni%` = Negate(`%in%`)
  
  ##########################
  #Load and clean the data
  ##########################
  
  dataDHS <- suppressWarnings(foreign::read.dta(filename, convert.dates = TRUE, convert.factors = F, 
                                       missing.type = FALSE, convert.underscore = TRUE, warn.missing.labels = TRUE))
  
  if (country=='Nepal'){
    dataDHS <- convert.nepali.dates(dataDHS, data.year)
    interview.year <- dataDHS$v007
  } else {
    interview.year <- round((dataDHS$v008-1)/12,0)+1900
  }
  
  if (is.na(data.year)){
    data.year.tmp <- min(interview.year)
    if (data.year.tmp>=1950){
      data.year <- data.year.tmp
    } else if (data.year.tmp>50 & data.year.tmp<100){
      data.year <- as.integer(paste0("19", data.year.tmp))
    } else if (data.year.tmp<=20 & data.year.tmp>=0){
      data.year <- as.integer(paste0("20", data.year.tmp))
    }
  }
  
  if (is.null(dataDHS$h9)){
    print('No Measles Vaccination Data')
    return(NULL)
  }
  
  if ("july" %ni% tolower(dataDHS$v006)) {
    dataDHS$v006 <- as.integer(as.character(dataDHS$v006))
  }
  if ("july" %ni% tolower(dataDHS$b1)) {
    dataDHS$b1 <- as.integer(as.character(dataDHS$b1))
  }  
  if ("july" %ni% tolower(dataDHS$h9m)) {
    dataDHS$h9m <- as.integer(as.character(dataDHS$h9m))
  }
  
  dataDHS$b2 <- as.integer(as.character(dataDHS$b2))
  dataDHS$h9y <- as.integer(as.character(dataDHS$h9y))
  dataDHS$h9d <- as.integer(as.character(dataDHS$h9d))
  
  
  #### Recode Years --> Ethiopia years are all off by (-8) years #####
  if (country == "Ethiopia"){
    interview.year  <- as.integer(interview.year + 8)
    dataDHS$b2    <- as.integer(dataDHS$b2 + 8)
    dataDHS$h9y   <- as.integer(dataDHS$h9y + 8)
    dataDHS$h3y   <- as.integer(dataDHS$h3y + 8)
    dataDHS$h5y   <- as.integer(dataDHS$h5y + 8)
    dataDHS$h7y   <- as.integer(dataDHS$h7y + 8)
  }
  
  # Make new variable of measles vacc status #
  if (is.factor(dataDHS$h9)) {
    dataDHS$h9.new[dataDHS$h9=="no"] <- 0                   #unvaccinated get a 0
    dataDHS$h9.new[grep("date on card", dataDHS$h9)] <- 1   #if they are vaccinated and there is a date, make h9 1
    dataDHS$h9.new[grep("reported by mother", dataDHS$h9)] <- 2   #if they are vaccinated and we do not know when,they get a 2
    dataDHS$h9.new[grep("marked on card", dataDHS$h9)] <- 2   #if they are vaccinated and we do not know when,they get a 2
    dataDHS$h9.new[dataDHS$h9=="dk"] <-  NA                 #Put in the NAs
    dataDHS$h9.new[dataDHS$h9=="don't know"] <-  NA                 #Put in the NAs
    dataDHS$h9.new[dataDHS$h9=="dont know"] <-  NA                 #Put in the NAs
    dataDHS$h9.new[is.na(dataDHS$h9)] <- NA
  } else {
    dataDHS$h9.new[dataDHS$h9==0] <- 0                 #assume all unvaccinated
    dataDHS$h9.new[dataDHS$h9==1] <- 1  #1 means the date is on the card...our 1
    dataDHS$h9.new[dataDHS$h9==2] <- 2  #2 means mother said vaccinated only...our 2
    dataDHS$h9.new[dataDHS$h9==3] <- 2  #3 means on card but no date markes....our 2
    dataDHS$h9.new[dataDHS$h9==8] <- NA #8 means don't know
    dataDHS$h9.new[is.na(dataDHS$h9)] <- NA  #Put in the NAs
  }
  
  # Make new variable of measles vacc status #
  if (is.factor(dataDHS$h9)) {
    dataDHS$measlesrec[dataDHS$h9=="no"] <- 0                   #unvaccinated get a 0
    dataDHS$measlesrec[grep("date on card", dataDHS$h9)] <- 1   #if they are vaccinated and there is a date, make h9 1
    dataDHS$measlesrec[grep("reported by mother", dataDHS$h9)] <- 1   #if they are vaccinated and we do not know when,they get a 2
    dataDHS$measlesrec[grep("marked on card", dataDHS$h9)] <- 1   #if they are vaccinated and we do not know when,they get a 2
    dataDHS$measlesrec[dataDHS$h9=="dk"] <-  NA                 #Put in the NAs
    dataDHS$measlesrec[dataDHS$h9=="don't know"] <-  NA                 #Put in the NAs
    dataDHS$measlesrec[dataDHS$h9=="dont know"] <-  NA                 #Put in the NAs
    dataDHS$measlesrec[is.na(dataDHS$h9)] <- NA
  } else {
    dataDHS$measlesrec[dataDHS$h9==0] <- 0 
    dataDHS$measlesrec[dataDHS$h9>=1 & dataDHS$h9<=3 & !is.na(dataDHS$h9)] <- 1
    dataDHS$measlesrec[dataDHS$h9==8] <- NA
    dataDHS$measlesrec[is.na(dataDHS$h9)] <- NA
  }
  
  # Vaccination proof on card
  dataDHS$VaccCard <- NA
  if (is.factor(dataDHS$h9)) {
    dataDHS$VaccCard[dataDHS$h9=="no"] <- 0                   #unvaccinated get a 0
    dataDHS$VaccCard[grep("date on card", dataDHS$h9)] <- 1   #if they are vaccinated and there is a date, make h9 1
    dataDHS$VaccCard[grep("reported by mother", dataDHS$h9)] <- NA   #if they are vaccinated and we do not know when,they get a 2
    dataDHS$VaccCard[grep("marked on card", dataDHS$h9)] <- 1   #if they are vaccinated and we do not know when,they get a 2
    dataDHS$VaccCard[dataDHS$h9=="dk"] <-  NA                 #Put in the NAs
    dataDHS$VaccCard[dataDHS$h9=="don't know"] <-  NA                 #Put in the NAs
    dataDHS$VaccCard[dataDHS$h9=="dont know"] <-  NA                 #Put in the NAs
    dataDHS$VaccCard[is.na(dataDHS$h9)] <- NA
  } else {
    dataDHS$VaccCard[dataDHS$h9==0] <- 0                 #assume all unvaccinated
    dataDHS$VaccCard[dataDHS$h9==1] <- 1  #1 means the date is on the card...our 1
    dataDHS$VaccCard[dataDHS$h9==2] <- NA  #2 means mother said vaccinated only...our NA
    dataDHS$VaccCard[dataDHS$h9==3] <- 1  #3 means on card but no date markes....our 1
    dataDHS$VaccCard[dataDHS$h9==8] <- NA #8 means don't know
    dataDHS$VaccCard[is.na(dataDHS$h9)] <- NA  #Put in the NAs
  }
  
  
  #Make the censoring column
  dataDHS$cens <- dataDHS$h9.new
  
  #survey date
  dataDHS$survey.date <- rep(NA, nrow(dataDHS))
  
  if (data.year>=2000) {
    dataDHS$survey.date  <- as.Date(sprintf("%d-%d-15", interview.year, dataDHS$v006))
  } else if (data.year<2000 & interview.year<100 & interview.year>50 & !is.na(interview.year)) {
    dataDHS$survey.date <- as.Date(sprintf("%d-%d-15", as.integer(paste("19", interview.year, sep="")), dataDHS$v006))
  } else if (data.year<2000 & interview.year>1900 & !is.na(interview.year)) {
    dataDHS$survey.date <- as.Date(sprintf("%d-%d-15", interview.year, dataDHS$v006))
  }
  
  if ("july" %in% tolower(dataDHS$v006)) {
    dataDHS$survey.date<-as.Date(dataDHS$survey.date, format="%Y-%B-%d")
  }
  
  if (country=="Malawi"){
    dataDHS$survey.date[dataDHS$survey.date == "2004-01-15"] <- "2004-10-15"
  }
  
  #age in months at time of survey
  dob.year <- rep(NA, nrow(dataDHS))
  max.yob <- max(dataDHS$b2)
  dob.year <- dataDHS$b2
  if (max.yob>50 & max.yob<100){
    dob.year[!is.na(dataDHS$b2) & dataDHS$b2>50 & dataDHS$b2<100] <- as.integer(paste("19", dataDHS$b2, sep=""))[!is.na(dataDHS$b2) & dataDHS$b2<2000 & dataDHS$b2>50]
  }
  
  dob <- rep(NA, nrow(dataDHS))
  dob <- paste(dob.year, dataDHS$b1, 15, sep="-")
  dob[is.na(dob.year) | is.na(dataDHS$b1)] <- NA
  
  dataDHS$dob <- as.Date(strptime(dob, "%Y-%m-%d"))
  dataDHS$age.in.months.surv <- (as.numeric(dataDHS$survey.date)-as.numeric(dataDHS$dob))/30.41667
  dataDHS$age.in.months <- dataDHS$age.in.months.surv
  
  # Remove individuals without unknown AGE. We cannot know if they were in the age range
  if (sum(is.na(dataDHS$age.in.months)>0)){
    dataDHS <- dataDHS[-which(is.na(dataDHS$age.in.months)),]
  }
  
  #for those with a date on the vacination card get the age in months
  #at the time of the vaccination
  dov.day <- rep(15, nrow(dataDHS))
  dov.day[which(!is.na(dataDHS$h9d) & dataDHS$h9d<=31)] <- dataDHS$h9d[which(!is.na(dataDHS$h9d) & dataDHS$h9d<=31)]
  dov.year <- rep(NA, nrow(dataDHS))
  
  if (sum(!is.na(dataDHS$h9y))>0){
    mean.dov.year <- mean(dataDHS$h9y, na.rm=TRUE)
    if (mean.dov.year>=1950){
      dov.year <- dataDHS$h9y
    } else if (mean.dov.year>50 & mean.dov.year<100){
      dov.year[which(!is.na(dataDHS$h9y))] <- as.integer(paste("19", dataDHS$h9y, sep=""))[which(!is.na(dataDHS$h9y))]
    }
    
    dov <- rep(NA, nrow(dataDHS))
    dov <- paste(dov.year, dataDHS$b1, dov.day, sep="-")
    dov[is.na(dov.year) | is.na(dataDHS$h9m)] <- NA
    dataDHS$dov <- as.Date(strptime(dov, "%Y-%m-%d"))
    
    # if dov were invalid, we will change the day to 15
    invalid.dov <- which(is.na(dataDHS$dov) & !is.na(dov))
    if (length(invalid.dov)>0){
      dov <- as.Date(strptime(paste(dov.year[invalid.dov], dataDHS$b1[invalid.dov], 15, sep="-"), "%Y-%m-%d"))
      dataDHS$dov[invalid.dov] <- dov
    }
    rm(dov)
  }
  
  # Missing Date on Vaccination Card for those that it appeared they had a date
  dataDHS$miss.date.card <- dataDHS$cens==1 & (is.na(dataDHS$h9y) | is.na(dataDHS$h9m) | 
                                                 dataDHS$h9y>data.year | dataDHS$h9m>12)
  dataDHS$cens[dataDHS$miss.date.card] <- 2
  dataDHS$dov[dataDHS$miss.date.card] <- NA 
  
  #print(dov[!is.na(dataDHS$cens) & dataDHS$cens==1])
  dataDHS$age.meas.months <- NA
  dataDHS$age.meas.months[which(dataDHS$cens==1)] <-
    ((as.numeric(dataDHS$dov) - as.numeric(dataDHS$dob)) / 30.41667)[which(dataDHS$cens==1)]
  
  #assume poeple with silly dates were vaccinated but we do not know when
  dataDHS$cens[!is.na(dataDHS$age.meas.months) & dataDHS$age.meas.months>1000 & dataDHS$cens==1] <- 2
  dataDHS$age.meas.months[!is.na(dataDHS$age.meas.months) & dataDHS$age.meas.months>1000] <- NA 
  dataDHS$cens[!is.na(dataDHS$age.meas.months) & dataDHS$age.meas.months<0 & dataDHS$dataDHS$cens==1] <- 2
  dataDHS$age.meas.months[!is.na(dataDHS$age.meas.months) & dataDHS$age.meas.months<0] <- NA   
  dataDHS$age.cens <- ifelse(!is.na(dataDHS$age.meas.months), dataDHS$age.meas.months, dataDHS$age.in.months.surv)
  
  #print(dataDHS$dov[!is.na(dataDHS$age.meas.months) & dataDHS$age.meas.months<0])
  
  # Create a variable of the mother's age at birth of the child #
  dataDHS$mother.age <- round(((dataDHS$b3 - dataDHS$v011)/12), digits = 0)
  
  # Create delivery location variables #
  dataDHS$delivery.home <- NA
  dataDHS$delivery.home[dataDHS$m15=="respondents home" | dataDHS$m15=="other home"] <- 1 
  dataDHS$delivery.home[is.na(dataDHS$delivery.home) & dataDHS$m15!="99" & !is.na(dataDHS$m15)] <- 0
  dataDHS$delivery.home <- factor(dataDHS$delivery.home, levels=c(0,1), labels=c("Non-home delivery", "Home delivery"))
  
  dataDHS$delivery.cat <- NA
  dataDHS$delivery.cat[dataDHS$m15=="respondents home" | dataDHS$m15=="other home"] <- 1
  dataDHS$delivery.cat[dataDHS$m15=="govt. hospital" | dataDHS$m15=="govt. health center" |
                         dataDHS$m15=="govt. health post" | dataDHS$m15=="other public"] <- 2
  dataDHS$delivery.cat[!is.na(dataDHS$m15) & is.na(dataDHS$delivery.cat)] <- 3
  dataDHS$delivery.cat <- factor(dataDHS$delivery.cat, levels=c(1,2,3), labels=c("Home delivery", "Govt. facility", "Other location"))
  
  # Create Education Variables #
  dataDHS$mom.educ <- as.numeric(dataDHS$v106)
  dataDHS$partner.educ[as.numeric(dataDHS$v701)<=4 & !is.na(dataDHS$v701)] <- dataDHS$v701[as.numeric(dataDHS$v701)<=4 & !is.na(dataDHS$v701)]
  dataDHS$parenteduc <-  NA
  educ.tmp <- rep(NA, nrow(dataDHS))
  educ.tmp[dataDHS$partner.educ >= dataDHS$mom.educ] <- 0
  educ.tmp[dataDHS$partner.educ < dataDHS$mom.educ] <- 1
  
  dataDHS$parenteduc[educ.tmp==0 & !is.na(educ.tmp)] <- dataDHS$partner.educ[educ.tmp==0 & !is.na(educ.tmp)]
  dataDHS$parenteduc[educ.tmp==1 & !is.na(educ.tmp)] <- dataDHS$mom.educ    [educ.tmp==1 & !is.na(educ.tmp)]
  rm(educ.tmp)
  
  dataDHS$birth.year <- dataDHS$b2
  
  #drop people outside of our age range
  dataDHS <- dataDHS[which(dataDHS$age.in.months.surv>=youngest.age & dataDHS$age.in.months.surv<=oldest.age & 
                             !is.na(dataDHS$age.in.months.surv)),]
  
  return(dataDHS)
}









#------------------------------------------------------------------------------------
# (2) Import and Clean data

# Function to read and clean data 
# - checks if data exists, then sets it up for use with clustering analysis
read.and.clean.data <- function(country, year, country.dhs, country.gps, age.range=c(13,60), directory=paste0(dirname(getwd()),"/data/")){
  
  if(!is.null(data.info)){
    clean.data <- suppressWarnings(dhs.data.import(filename = paste0(directory,country.dhs),
                                                   youngest.age = age.range[1], oldest.age = age.range[2], 
                                                   country=country, data.year = year, model.type = "all" ))
    
    var.keep <- as.integer(unlist(sapply(c("caseid","hhid2","v001","v002","v005","v024","v101","h9","cens","measlesrec","dov","dob","age.cens",
                                           "birth.year","age.in.months","age.in.months.surv"), grep, colnames(clean.data))))
    clean.data <- clean.data[,var.keep]
    clean.data$id <- paste0(rep(tolower(substr(country,1,4)), nrow(clean.data)),1:nrow(clean.data))
    
    # GPS data
    gps.data <- read.dbf(file=paste0("../data/",data.info$dbfName))
    gps.data <- gps.data[, c("DHSID", "DHSCLUST", "URBAN_RURA", "LATNUM", "LONGNUM", "ADM1FIPSNA")]
    
    # Merge Survey and GPS
    clean.data <- merge(clean.data, gps.data, by.x="v001", by.y="DHSCLUST", all=TRUE)
    rm(gps.data)
    
    # Create new Lat and Long variables in km
    #   * Conversion using latlong2grid()
    #   * Long = x; Lat = y
    coords <- SpatialEpi::latlong2grid(clean.data[,c("LONGNUM","LATNUM")])
    
    # Isolate Clusters and Households
    clean.data$hhid2 <- paste("hh",clean.data$v001,clean.data$v002, sep="-")
    clean.data$weight <- as.numeric(clean.data$v005/1000000)
    
    clean.data <- cbind(clean.data$id, clean.data$hhid2, clean.data$v001, coords, 
                        duplicated(clean.data$v001), duplicated(clean.data$hhid2), 
                        clean.data$URBAN_RURA, clean.data$measlesrec, 
                        clean.data$age.in.months.surv, clean.data$ADM1FIPSNA, clean.data$v024, clean.data$v101,
                        clean.data$LONGNUM, clean.data$LATNUM, clean.data$weight)
    
    colnames(clean.data) <- c("id","hhid2","cluster","x","y","dupl_clust","dupl_hh","urban","measvacc", "age",
                              "district", "region1","region2","long", "lat", "sampwt")
    
    urban <- as.character(clean.data[,"urban"])
    clean.data[,'urban'] <- (urban == 'U')
    clean.data[,'hhid2'] <- as.character(clean.data[,'hhid2'])
    clean.data[,'cluster'] <- as.integer(clean.data[,'cluster'])
    clean.data$measvaccNA0 <- clean.data$measvacc
    clean.data$measvaccNA0[which(is.na(clean.data$measvaccNA0))] <- 0
    
    return(clean.data)
  }
}














#...................................................................................................................
#...................................................................................................................
#...................................................................................................................



# DHS DATA IMPORT ---------------------------------------------------------




####################################################################################
##' 
##' Name: dhs.data.import   
##' Description: Function gets and cleans DHS data (.DTA format) for relevant measles variables
##' 
##' Arguments:
##' - country - specific country selected
##' - data.year - year of DHS survey
##' - filename - name of STATA file with the DHS data. If country and year are entered, filename is not required.
##' - numvacc.dpt - number of different vaccination variables that can include DPT in that country (i.e. h3, sdhh1,...)
##' - youngest.age - the youngest age to look at in months...should always be 0
##' - oldest.age - the oldest age to look at....default to 60 months
##' - directory - local location of DHS data files
##' 
##' Returns:
##' - data.frame with all data from DHS survey for selected country and year
##' 
##' Includes: - "censoring" variable (same as dhs.to.vac.cens function),
##'           - recoded DPT and Measles doses and dates
##'           

dhs.data.import <- cmpfun(function(country, data.year=NA, filename, filename.gps,
                            directory=getwd(), 
                            youngest.age=0, oldest.age=60, 
                            model.type="all", merge.gps=TRUE) {
  
  `%ni%` = Negate(`%in%`)
  
  #---------------------------------------------------------------------------------
  #   Load and clean the data
  #---------------------------------------------------------------------------------
  
  # ------  Import STATA file -------------------------------
 
  # dataDHS <- suppressWarnings(foreign::read.dta(filename, convert.dates = TRUE, convert.factors = TRUE, 
  #                                      missing.type = FALSE, convert.underscore = TRUE, warn.missing.labels = TRUE))
  
  # Use 'haven' package to maintain variable labels
  dataDHS <- haven::read_dta(file.path(directory, filename))
  var.labels <- sjmisc::get_label(dataDHS)
  var.value.labels <- sjmisc::get_labels(dataDHS)
  
  
  # ------  Check for measles vaccination data  -------------
  
  # -- If not present, end function and don't return any data
  if (is.null(dataDHS$h9)){
    print('No Measles Vaccination Data')
    return(NULL)
  }
  
  
  # ------  Fix Measles Vaccination Variables  ----------
  
  # Make new variable of measles vacc status #
  dataDHS$h9.new <- NA
  if (is.factor(dataDHS$h9)) {
    h9.tmp <- tolower(dataDHS$h9)
    dataDHS$h9.new[h9.tmp=="no"] <- 0                   #unvaccinated get a 0
    dataDHS$h9.new[grep("date on card", h9.tmp)] <- 1   #if they are vaccinated and there is a date, make h9 1
    dataDHS$h9.new[grep("reported by mother", h9.tmp)] <- 2   #if they are vaccinated and we do not know when,they get a 2
    dataDHS$h9.new[grep("marked on card", h9.tmp)] <- 2   #if they are vaccinated and we do not know when,they get a 2
    dataDHS$h9.new[h9.tmp=="dk"] <-  NA                 #Put in the NAs
    dataDHS$h9.new[h9.tmp=="don't know"] <-  NA                 #Put in the NAs
    dataDHS$h9.new[h9.tmp=="dont know"] <-  NA                 #Put in the NAs
    dataDHS$h9.new[is.na(h9.tmp)] <- NA
  } else {
    dataDHS$h9.new[dataDHS$h9==0] <- 0                 #assume all unvaccinated
    dataDHS$h9.new[dataDHS$h9==1] <- 1  #1 means the date is on the card...our 1
    dataDHS$h9.new[dataDHS$h9==2] <- 2  #2 means mother said vaccinated only...our 2
    dataDHS$h9.new[dataDHS$h9==3] <- 2  #3 means on card but no date marks....our 2
    dataDHS$h9.new[dataDHS$h9==8] <- NA #8 means don't know
  }
  
  # Make new variable of measles vacc status #
  dataDHS$measlesrec <- NA
  dataDHS$measlesrec[is.na(dataDHS$h9.new)] <- NA
  dataDHS$measlesrec[dataDHS$h9.new==0] <- 0
  dataDHS$measlesrec[dataDHS$h9.new==1] <- 1
  dataDHS$measlesrec[dataDHS$h9.new==2] <- 1
  
  # Vaccination proof on card
  dataDHS$VaccCard <- NA
  if (is.factor(dataDHS$h9)) {
    dataDHS$VaccCard[dataDHS$h9=="no"] <- 0                   #unvaccinated get a 0
    dataDHS$VaccCard[grep("date on card", dataDHS$h9)] <- 1   #if they are vaccinated and there is a date, make h9 1
    dataDHS$VaccCard[grep("reported by mother", dataDHS$h9)] <- NA   #if they are vaccinated and we do not know when,they get a 2
    dataDHS$VaccCard[grep("marked on card", dataDHS$h9)] <- 1   #if they are vaccinated and we do not know when,they get a 2
    dataDHS$VaccCard[dataDHS$h9=="dk"] <-  NA                 #Put in the NAs
    dataDHS$VaccCard[dataDHS$h9=="don't know"] <-  NA                 #Put in the NAs
    dataDHS$VaccCard[dataDHS$h9=="dont know"] <-  NA                 #Put in the NAs
    dataDHS$VaccCard[is.na(dataDHS$h9)] <- NA
  } else {
    dataDHS$VaccCard[dataDHS$h9==0] <- 0                 #assume all unvaccinated
    dataDHS$VaccCard[dataDHS$h9==1] <- 1  #1 means the date is on the card...our 1
    dataDHS$VaccCard[dataDHS$h9==2] <- NA  #2 means mother said vaccinated only...our NA
    dataDHS$VaccCard[dataDHS$h9==3] <- 1  #3 means on card but no date marks....our 1
    dataDHS$VaccCard[dataDHS$h9==8] <- NA #8 means don't know
    dataDHS$VaccCard[is.na(dataDHS$h9)] <- NA  #Put in the NAs
  }
  
  #Make the censoring column
  dataDHS$cens <- dataDHS$h9.new
  
  
  # ------  Fix Dates  --------------------------------------
  
  # Fix Dates for Nepali data
  if (country=='Nepal'){
    dataDHS <- convert.nepali.dates(dataDHS, data.year)
    interview.year <- dataDHS$v007
  } else {
    interview.year <- round((dataDHS$v008-1)/12,0)+1900
  }
  
  if (is.na(data.year)){
    data.year.tmp <- min(interview.year)
    if (data.year.tmp>=1950){
      data.year <- data.year.tmp
    } else if (data.year.tmp>50 & data.year.tmp<100){
      data.year <- as.integer(paste0("19", data.year.tmp))
    } else if (data.year.tmp<=20 & data.year.tmp>=0){
      data.year <- as.integer(paste0("20", data.year.tmp))
    }
  }
  
  # Month of Interview
  if ("july" %ni% tolower(dataDHS$v006)) {
    dataDHS$v006 <- as.integer(as.character(dataDHS$v006))
  }
  # Month of Birth
  if ("july" %ni% tolower(dataDHS$b1)) {
    dataDHS$b1 <- as.integer(as.character(dataDHS$b1))
  }  
  # Month of Measles Vaccination
  if ("july" %ni% tolower(dataDHS$h9m)) {
    dataDHS$h9m <- as.integer(as.character(dataDHS$h9m))
  }
  
  dataDHS$b2 <- as.integer(as.character(dataDHS$b2))      # Year of Birth
  dataDHS$h9y <- as.integer(as.character(dataDHS$h9y))    # Year of Measles Vaccination
  dataDHS$h9d <- as.integer(as.character(dataDHS$h9d))    # Day of Measles Vaccination
  
  
  # Recode Years --> Ethiopia years are all off by (-8) years #####
  if (country == "Ethiopia"){
    interview.year  <- as.integer(interview.year + 8)
    dataDHS$b2    <- as.integer(dataDHS$b2 + 8)
    dataDHS$h9y   <- as.integer(dataDHS$h9y + 8)
    dataDHS$h3y   <- as.integer(dataDHS$h3y + 8)
    dataDHS$h5y   <- as.integer(dataDHS$h5y + 8)
    dataDHS$h7y   <- as.integer(dataDHS$h7y + 8)
  }
  
  #survey date
  dataDHS$survey.date <- NA
  
  if (data.year>=2000) {
    dataDHS$survey.date  <- as.Date(sprintf("%d-%d-15", interview.year, dataDHS$v006))
  } else if (data.year<2000 & interview.year<100 & interview.year>50 & !is.na(interview.year)) {
    dataDHS$survey.date <- as.Date(sprintf("%d-%d-15", as.integer(paste("19", interview.year, sep="")), dataDHS$v006))
  } else if (data.year<2000 & interview.year>1900 & !is.na(interview.year)) {
    dataDHS$survey.date <- as.Date(sprintf("%d-%d-15", interview.year, dataDHS$v006))
  }
  
  if ("july" %in% tolower(dataDHS$v006)) {
    dataDHS$survey.date<-as.Date(dataDHS$survey.date, format="%Y-%B-%d")
  }
  
  if (country=="Malawi"){
    dataDHS$survey.date[dataDHS$survey.date == "2004-01-15"] <- "2004-10-15"
  }
  
  #age in months at time of survey
  dob.year <- rep(NA, nrow(dataDHS))
  max.yob <- max(interview.year)
  dob.year <- dataDHS$b2
  if (max.yob>50 & max.yob<100){
    dob.year[!is.na(dataDHS$b2) & dataDHS$b2>50 & dataDHS$b2<100] <- as.integer(paste("19", dataDHS$b2, sep=""))[!is.na(dataDHS$b2) & dataDHS$b2<2000 & dataDHS$b2>50]
  }
  
  dob <- rep(NA, nrow(dataDHS))
  dob <- paste(dob.year, dataDHS$b1, 15, sep="-")
  dob[is.na(dob.year) | is.na(dataDHS$b1)] <- NA
  
  dataDHS$dob <- as.Date(strptime(dob, "%Y-%m-%d"))  # Date of Birth
  dataDHS$age.in.months.surv <- (as.numeric(dataDHS$survey.date)-as.numeric(dataDHS$dob))/30.41667 # Age in months at time of survey
  dataDHS$age.in.months <- dataDHS$age.in.months.surv 
  
  # Remove individuals without unknown AGE. We cannot know if they were in the age range
  if (sum(is.na(dataDHS$age.in.months)>0)){
    dataDHS <- dataDHS[-which(is.na(dataDHS$age.in.months)),]
  }
  
  #for those with a date on the vacination card get the age in months
  #at the time of the vaccination
  dov.day <- rep(15, nrow(dataDHS))
  dov.day[which(!is.na(dataDHS$h9d) & dataDHS$h9d<=31)] <- dataDHS$h9d[which(!is.na(dataDHS$h9d) & dataDHS$h9d<=31)]
  dov.year <- rep(NA, nrow(dataDHS))
  
  if (sum(!is.na(dataDHS$h9y))>0){
    mean.dov.year <- mean(dataDHS$h9y, na.rm=TRUE)
    if (mean.dov.year>=1950){
      dov.year <- dataDHS$h9y
    } else if (mean.dov.year>50 & mean.dov.year<100){
      dov.year[which(!is.na(dataDHS$h9y))] <- as.integer(paste("19", dataDHS$h9y, sep=""))[which(!is.na(dataDHS$h9y))]
    }
    
    dov <- rep(NA, nrow(dataDHS))
    dov <- paste(dov.year, dataDHS$b1, dov.day, sep="-")
    dov[is.na(dov.year) | is.na(dataDHS$h9m)] <- NA
    dataDHS$dov <- as.Date(strptime(dov, "%Y-%m-%d"))
    
    # if dov were invalid, we will change the day to 15
    invalid.dov <- which(is.na(dataDHS$dov) & !is.na(dov))
    if (length(invalid.dov)>0){
      dov <- as.Date(strptime(paste(dov.year[invalid.dov], dataDHS$b1[invalid.dov], 15, sep="-"), "%Y-%m-%d"))
      dataDHS$dov[invalid.dov] <- dov
    }
    rm(dov)
  }
  
  # Missing Date on Vaccination Card for those that it appeared they had a date
  dataDHS$miss.date.card <- NA
  dataDHS$miss.date.card <- dataDHS$cens==1 & (is.na(dataDHS$h9y) | is.na(dataDHS$h9m) | 
                                                 dataDHS$h9y>(data.year+1) | dataDHS$h9m>12)
  dataDHS$cens[dataDHS$miss.date.card] <- 2
  dataDHS$dov[dataDHS$miss.date.card] <- NA 
  
  #print(dov[!is.na(dataDHS$cens) & dataDHS$cens==1])
  dataDHS$age.meas.months <- NA
  dataDHS$age.meas.months[which(dataDHS$cens==1)] <-
    ((as.numeric(dataDHS$dov) - as.numeric(dataDHS$dob)) / 30.41667)[which(dataDHS$cens==1)]
  
  #assume poeple with silly dates were vaccinated but we do not know when
  dataDHS$cens[(!is.na(dataDHS$age.meas.months) & dataDHS$age.meas.months>1000) & dataDHS$cens==1] <- 2
  dataDHS$cens[(!is.na(dataDHS$age.meas.months) & dataDHS$age.meas.months<0)    & dataDHS$cens==1] <- 2
  dataDHS$age.meas.months[!is.na(dataDHS$age.meas.months) & dataDHS$age.meas.months>1000] <- NA
  dataDHS$age.meas.months[!is.na(dataDHS$age.meas.months) & dataDHS$age.meas.months<0] <- NA   
  dataDHS$age.cens <- ifelse(!is.na(dataDHS$age.meas.months), dataDHS$age.meas.months, dataDHS$age.in.months.surv)
  
  # Create a variable of the mother's age at birth of the child #
  dataDHS$mother.age <- round(((dataDHS$b3 - dataDHS$v011)/12), digits = 0)
  
  # ----  Clean other variables of interest  -------------------------

  # Create Education Variables #
  dataDHS$mom.educ <- dataDHS$partner.educ <- dataDHS$parenteduc <- NA
  dataDHS$mom.educ <- as.numeric(dataDHS$v106)
  dataDHS$partner.educ[as.numeric(dataDHS$v701)<=4 & !is.na(dataDHS$v701)] <- dataDHS$v701[as.numeric(dataDHS$v701)<=4 & !is.na(dataDHS$v701)]
  educ.tmp <- rep(NA, nrow(dataDHS))
  educ.tmp[dataDHS$partner.educ >= dataDHS$mom.educ] <- 0
  educ.tmp[dataDHS$partner.educ < dataDHS$mom.educ] <- 1
  
  dataDHS$parenteduc[educ.tmp==0 & !is.na(educ.tmp)] <- dataDHS$partner.educ[educ.tmp==0 & !is.na(educ.tmp)]
  dataDHS$parenteduc[educ.tmp==1 & !is.na(educ.tmp)] <- dataDHS$mom.educ    [educ.tmp==1 & !is.na(educ.tmp)]
  rm(educ.tmp)
  
  dataDHS$birth.year <- dataDHS$b2 
  
  
  #================================================
  # Merge GPS data with Survey data
  
  if(merge.gps){
    
      # GPS data
      gps.data <- foreign::read.dbf(file=file.path(directory,filename.gps))
      gps.data <- gps.data[, c("DHSID", "DHSYEAR", "DHSCLUST", "URBAN_RURA", "LATNUM", "LONGNUM", 
                               "ADM1FIPSNA", "DHSREGNA", "DATUM")]
      
      gps.data$LATNUM[round(gps.data$LATNUM,8)==0] <- 0
      gps.data$LONGNUM[round(gps.data$LONGNUM,8)==0] <- 0
      
      # Merge Survey and GPS
      # - Keeping all survey data, but not unmatched gps data
      dataDHS <- merge(x=dataDHS, y=gps.data, by.x="v001", by.y="DHSCLUST", all.x=TRUE, all.y=FALSE)
      rm(gps.data)
      
      # Create new Lat and Long variables in km
      #   * Conversion using latlong2grid()
      #   * Long = x; Lat = y
      coords <- SpatialEpi::latlong2grid(dataDHS[,c("LONGNUM","LATNUM")])
      
      # Isolate Clusters and Households
      dataDHS$hhid2 <- paste("hh",dataDHS$v001,dataDHS$v002, sep="-")
      dataDHS$weight <- as.numeric(dataDHS$v005/1000000)
      dataDHS$x.coord.km <- coords$x
      dataDHS$y.coord.km <- coords$y
      
      setnames(dataDHS, 'v001', 'cluster')
      setnames(dataDHS, 'URBAN_RURA', 'urban')
      
      urban <- as.character(dataDHS[,"urban"])
      dataDHS[,'urban'] <- (urban == 'U')  # make Urban/Rural variable TRUE/FALSE 
      dataDHS[,'hhid2'] <- as.character(dataDHS[,'hhid2'])
      dataDHS[,'cluster'] <- as.integer(dataDHS[,'cluster'])
  }
  
  dataDHS$id <- paste0(rep(tolower(substr(country,1,4)), nrow(dataDHS)),1:nrow(dataDHS))
  
  dataDHS$measvaccNA0 <- NA
  dataDHS$measvaccNA0 <- dataDHS$measlesrec
  dataDHS$measvaccNA0[is.na(dataDHS$measvaccNA0)] <- 0
  #================================================
  
  #drop people outside of our age range
  dataDHS <- dataDHS[which(dataDHS$age.in.months.surv>=youngest.age & dataDHS$age.in.months.surv<=oldest.age & 
                             !is.na(dataDHS$age.in.months.surv)),]
  
  return(list(data=dataDHS, var.labels=var.labels, var.value.labels=var.value.labels))
})


#///////////////////////////////////////////////////////////////////////////////////////////
# # EXAMPLE:
#
# Zamb2007 <- dhs.data.import(country='Zambia', data.year=2007, filename=NULL,
#                             directory="C:/Users/Shaun/JHSPH/IDD SVN Repository/Measles/data",
#                             youngest.age=13, oldest.age=60,
#                             model.type="all", merge.gps=TRUE)
#
#///////////////////////////////////////////////////////////////////////////////////////////




# Get rid of unwanted variables, if desired
clean.dhs <- function(dataDHS) {
  dataDHS <- data.frame(dataDHS$id, dataDHS$hhid2, dataDHS$cluster, dataDHS$x.coord.km, dataDHS$y.coord.km,
                   duplicated(dataDHS$cluster), duplicated(dataDHS$hhid2), 
                   dataDHS$urban, dataDHS$measlesrec, 
                   dataDHS$age.in.months.surv, dataDHS$ADM1FIPSNA, dataDHS$v024, dataDHS$v101,
                   dataDHS$LONGNUM, dataDHS$LATNUM, dataDHS$weight)
  
  colnames(clean.data) <- c("id","hhid2","cluster","x.coord.km","y.coord.km","dupl_clust","dupl_hh","urban","measvacc", "age",
                            "district", "region1","region2","long", "lat", "sampwt")
  return(dataDHS)
}





# CLUSTER AGGREGATION -----------------------------------------------------


##' Name: aggr.dhs.clusters
##' Description: Aggregate data to the cluster level
##' 
aggr.dhs.clusters <- function(data, vaccbreaks=seq(0,1,.1), treat_NAs="unvacc"){
  # File: DHS_Import$Clean_v2.R
  
  if (treat_NAs=="remove"){
    data <- data[!is.na(data$measlesrec),]
    meas.vacc <- data$measlesrec
  } else if(treat_NAs=="unvacc"){
    meas.vacc <- data$measvaccNA0
  }
  data$meas.vacc <- meas.vacc
  
  dupl.index <- which(duplicated(data[,"cluster"])==FALSE)
  dat.tmp <- data[dupl.index,]
  
  # If aggregating from mapping, need to include matched region variable
  if (!is.null(dat.tmp$regid)){
    dat.clust <- data.frame(id=dat.tmp$cluster, long=dat.tmp$LONGNUM, lat=dat.tmp$LATNUM, x=dat.tmp$x.coord.km, y=dat.tmp$y.coord.km, 
                            district=dat.tmp$district, region=dat.tmp$region1, region2=dat.tmp$region2,
                            sampwt=dat.tmp$weight, urban=dat.tmp$urban, 
                            num.indiv=NA, num.vacc=NA, vaccrate=NA, vacc.cat=NA, vacc.group=NA, 
                            regid=dat.tmp$regid)
  }else{
    dat.clust <- data.frame(id=dat.tmp$cluster, long=dat.tmp$LONGNUM, lat=dat.tmp$LATNUM, x=dat.tmp$x.coord.km, y=dat.tmp$y.coord.km, 
                            district=dat.tmp$ADM1FIPSNA, region=dat.tmp$v024, region2=dat.tmp$v101,
                            sampwt=dat.tmp$weight, urban=dat.tmp$urban, 
                            num.indiv=NA, num.vacc=NA, vaccrate=NA, vacc.cat=NA, vacc.group=NA)
  }
  dat.clust <- dat.clust[order(dat.clust$id),]
  
  data <- data[, c('cluster', 'meas.vacc')]
  data <- data.table(data)        # Convert data to data.table for faster indexing and tabling
  
  # Tabulate data for clusters
  dat.clust$num.indiv <- data[ ,.N, by=list(cluster)][,N]
  dat.clust$num.vacc <- data[,.(Vacc=sum(meas.vacc)),by=list(cluster)][,Vacc]
  dat.clust$vaccrate <- round(dat.clust$num.vacc / dat.clust$num.indiv, 4)
  vacc.cat.labels <- vaccbreaks[-length(vaccbreaks)]
  vacc.group.labels <- 1:length(vacc.cat.labels)
  
  dat.clust$vacc.cat <- cut(dat.clust$vaccrate, vaccbreaks, include.lowest=T, right=T, labels=vacc.cat.labels)
  dat.clust$vacc.group <- cut(dat.clust$vaccrate, vaccbreaks, include.lowest=T, right=T, labels=vacc.group.labels)
  
  dat.clust$district <- as.character(dat.clust$district)
  dat.clust$region <- as.character(dat.clust$region)
  dat.clust$region2 <- as.character(dat.clust$region2)
  
  return(dat.clust)
}



# Reduce Aggregated Cluster data to the variables of interest and calculate p, alpha, and delta
DHS.to.spatialtau.clust <- function(data=clust.dat){
  data <- data[,c('vaccrate', 'x', 'y', 'num.indiv', 'sampwt','urban')]
  colnames(data) <- c('p','x','y','s','weight','urban')
  
  data[,'p'] <- 1-data[,'p']
  data[,'alpha'] <- round(data[,'s']*data[,'weight'],2)             #  Weighted number of individuals
  data[,'delta'] <- round(data[,'s']*data[,'weight']*data[,'p'],2)  #  weighted number of susceptibles
  
  # Remove clusters without a coordinate
  data <- data[!(round(data$x,5)==0 & round(data$y,5)==0),]
  
  return(data)
}







# MAP DHS -----------------------------------------------------------------


# Map DHS Country and Clusters, ggplot2

# Download Boundaries from web zip file
### http://spatialdata.dhsprogram.com/boundaries/#view=table&countryId=AL


map.country.clust.gg <- cmpfun(function(country='Philippines', year=2008, age.range=c(13,60),
                                        filename, filename.gps, treat_NAs="unvacc",
                                        directory=paste0(dirname(getwd()),"/data/"), shapefile.dir='../data/DHS_boundaries/',
                                        pt.size=2, pt.transpar=.60, scale.dist=50, scale.loc='bottomright'){

  iso.tmp <- maps::iso.alpha(country, 2)
  print(paste0(country,', ', year, ', ', iso.tmp))
  
  
  #///////////////////////////////////////////////////////////////////////////////////////////
  # Vaccination Rates by District and Cluster
  
  # Import and Clean data
  dhs.data.list <- dhs.data.import(country=country, data.year=year, filename=filename, filename.gps=filename.gps,
                                   directory=directory, youngest.age=age.range[1], oldest.age=age.range[2], 
                                   model.type="all", merge.gps=TRUE)
  dhs.data <- dhs.data.list$data
  dhs.data <- dhs.data[!(dhs.data$LATNUM==0 & dhs.data$LONGNUM==0),]
  
  
  if (is.null(dhs.data)){
    return('GPS not available')
  }
  dhs.data <- dhs.data[(dhs.data$LATNUM!=0 & dhs.data$LONGNUM!=0),]
  
  if (sum(is.na(dhs.data$measlesrec)) >= (.5*nrow(dhs.data))){
    return('Measles Vaccination Data Not Available')
  }
  
  if (sum(is.na(dhs.data$LATNUM) | is.na(dhs.data$LONGNUM)) >= (.5*nrow(dhs.data))){
    return('Spatial Data Not Available for DHS Clusters')
  }
  
  
  # if(sum(dhs.data$district=="NULL")>0){
  #   dhs.data$district <- as.character(dhs.data$district)
  #   dhs.data$district[which(dhs.data$district=="NULL")] <- NA
  # }
  # if certain countries, edit the districts to match the correct names that are in the shapefiles
  
  region.values.1 <- as.character(unlist(dhs.data.list$var.value.labels['v024']))
  region.values.2 <- as.character(unlist(dhs.data.list$var.value.labels['v101']))
  
  # Get rid of strange accent problems
  Encoding(region.values.1) <- Encoding(region.values.2) <- 'latin1'
  region.values.1 <- stri_trans_general(region.values.1, 'Latin-ASCII')
  region.values.2 <- stri_trans_general(region.values.2, 'Latin-ASCII')
  
  dhs.data$district <- dhs.data$ADM1FIPSNA
  dhs.data$region1 <- region.values.1[dhs.data$v024]  
  dhs.data$region2 <- region.values.2[dhs.data$v101]  
  
  
  if (country=='Zimbabwe'){
    dhs.data$district <- sub('Matabeland', 'Matebeleland', dhs.data$district)
    dhs.data$district <- sub(' Province', '', dhs.data$district)
  }
  if (country=='Mozambique'){
    dhs.data$district[which(dhs.data$district=='Maputo (city)')] <- "Maputo City"
  }
  if (country=='Tanzania' & year==1999){
    dhs.data$region1[which(dhs.data$region1=='arusha')] <- "arusha/manyara"
    dhs.data$region1[which(dhs.data$region1=='rest of zanzibar')] <- "unguja"
    dhs.data$region1[which(dhs.data$region1=='coast')] <- "pwani"
    dhs.data$region1[which(dhs.data$region1=='shinyinga')] <- "shinyanga"
    dhs.data$region1[which(dhs.data$region1=='siginda')] <- "singida"
  }  
  # if (country=='Benin'){
  #   dhs.data$district <- as.character(dhs.data$region1)
  # }
  if (country=='Jordan'){
    dhs.data$district <- tolower(as.character(dhs.data$district))
    dhs.data$district[grep('at ', dhs.data$district)] <- substring(dhs.data$district[grep('at ', dhs.data$district)], first=4)
    dhs.data$district[grep('az ', dhs.data$district)] <- substring(dhs.data$district[grep('az ', dhs.data$district)], first=4)
    dhs.data$district[grep('al ', dhs.data$district)] <- substring(dhs.data$district[grep('al ', dhs.data$district)], first=4)
    dhs.data$district[which(dhs.data$district=='tafilah')] <- "tafiela"
    dhs.data$district[which(dhs.data$district=='aqabah')] <- "aqaba"
    dhs.data$district[which(dhs.data$district=='maan')] <- "ma'an"
  }
  if (country=='Tanzania' & year==1999){
    dhs.data$region1[which(dhs.data$region1=='arusha')] <- "arusha/manyara"
    dhs.data$region1[which(dhs.data$region1=='rest of zanzibar')] <- "unguja"
    dhs.data$region1[which(dhs.data$region1=='coast')] <- "pwani"
    dhs.data$region1[which(dhs.data$region1=='shinyinga')] <- "shinyanga"
    dhs.data$region1[which(dhs.data$region1=='siginda')] <- "singida"
  } 
  
  
  #///////////////////////////////////////////////////////////////////////////////////////////
  # Get Shapefile and map data
  
  shapefiles <- list.files(shapefile.dir)
  shape_file_name <- paste0(tolower(country), year); shape_file_name <- gsub(" ", "", shape_file_name)
  if (shape_file_name=="tanzania1999") shape_file_name="tanzania1996"
  if (shape_file_name %in% shapefiles){
    country.map <- readOGR(paste0('../data/DHS_boundaries/', shape_file_name, '/shps/sdr_subnational_boundaries.shp'),
                           layer='sdr_subnational_boundaries')
  } else {
    worlddat <- readOGR('../data/DHS_boundaries/sdr_subnational_data_world/shps/sdr_subnational_data.shp', layer='sdr_subnational_data')
    country.map <- subset(worlddat, as.character(worlddat$ISO)==iso.tmp)
    #country.map <- subset(country.map, as.integer(country.map$SVYYEAR)==as.integer(year))
    worlddat <- NULL
  }
  
  method2=FALSE
  
  #///////////////////////////////////////////////////////////////////////////////////////////
  # Identify matching district/region variables
  shp.tmp1 <- sort(unique(tolower(as.character(country.map@data$REGNAME))))
  shp.tmp2 <- sort(unique(tolower(as.character(country.map@data$DHSREGEN))))
  
  if (length(shp.tmp1)==0 & length(shp.tmp2)==0){
    return(paste0('No regions or districts defined in shapefile for ', country, ', ', year))
  }
  
  dhs.tmp1 <- sort(unique(tolower(as.character(dhs.data$district))))
  dhs.tmp2 <- sort(unique(tolower(as.character(dhs.data$region1))))
  dhs.tmp3 <- sort(unique(tolower(as.character(dhs.data$region2))))
  
  shp.list <- list(shp.tmp1, shp.tmp2)
  dhs.list <- list(dhs.tmp1, dhs.tmp2, dhs.tmp3)
  
  match.dists <- c( sum(stringdist(shp.list[[1]], dhs.list[[1]])), 
                    sum(stringdist(shp.list[[1]], dhs.list[[2]])), 
                    sum(stringdist(shp.list[[1]], dhs.list[[3]])),
                    sum(stringdist(shp.list[[2]], dhs.list[[1]])),
                    sum(stringdist(shp.list[[2]], dhs.list[[2]])), 
                    sum(stringdist(shp.list[[2]], dhs.list[[3]])))
  
  match.num <- which.min(match.dists)
  
  # If minimum match distance is still high, stop function
  if (length(match.num)==0){
    #return(paste0('Cannot match districts for ', country, ', ', year))
    method2=TRUE
  }
  
  # Print Best Matches
  print(list(Map.Regions=shp.list[[ifelse(match.num<=3, 1, 2)]],
             DHS.Regions=dhs.list[[match.num-ifelse(match.num<=3, 0, 3)]]))
  
  if (min(match.dists)>10){
    print(min(match.dists))
    #return(paste0('Regions too different to match for ', country, ', ', year))
    method2=TRUE
  }
  
  
  #///////////////////////////////////////////////////////////////////////////////////////////
  # Create variables for IDs based on match
  if (method2==FALSE){
    if (match.num==1){
      country.map@data$id = tolower(as.character(country.map@data$REGNAME))# create ID var in map data
      dhs.reg.var <- tolower(as.character(dhs.data$district))# dhs region variable
    }else if(match.num==2){
      country.map@data$id = tolower(as.character(country.map@data$REGNAME))
      dhs.reg.var <- tolower(as.character(dhs.data$region1)) 
    }else if(match.num==3) {
      country.map@data$id = tolower(as.character(country.map@data$REGNAME))
      dhs.reg.var <- tolower(as.character(dhs.data$region2))
    }else if(match.num==4) {
      country.map@data$id = tolower(as.character(country.map@data$DHSREGEN))
      dhs.reg.var <- tolower(as.character(dhs.data$district))
    }else if(match.num==5) {
      country.map@data$id = tolower(as.character(country.map@data$DHSREGEN))
      dhs.reg.var <- tolower(as.character(dhs.data$region1))
    }else if(match.num==6) {
      country.map@data$id = tolower(as.character(country.map@data$DHSREGEN))
      dhs.reg.var <- tolower(as.character(dhs.data$region2))
    }
    
    map.vals <- sort(unique(country.map@data$id)) # unique map regions
    dhs.reg.index <- match(dhs.reg.var, sort(unique(dhs.reg.var))) # Indicator for match
    dhs.data$regid <- map.vals[dhs.reg.index] # new region id
    
    # Set up cluster data
    clust.dat <- aggr.dhs.clusters(dhs.data, vaccbreaks=seq(0,1,.1), treat_NAs="unvacc")
    clust.dat <- clust.dat[!is.na(clust.dat$sampwt) & !is.na(clust.dat$long) & !is.na(clust.dat$long),]
    
    # Set up Region data
    vaccrates <- as.data.frame.array(t(table(dhs.data$measvaccNA0, dhs.data$regid)))
    vaccrates$vaccrate <- vaccrates[,"1"] / (vaccrates[,"0"] + vaccrates[,"1"])
    vaccrates$id <- tolower(row.names(vaccrates))
    colnames(vaccrates) <- c('unvacc', 'vacc', 'vaccrate', 'id')
    
    distr.vacc <- aggregate((clust.dat$vaccrate*clust.dat$sampwt), by=list(clust.dat$regid), FUN="sum")
    distr.sampwt <- aggregate(clust.dat$sampwt, by=list(clust.dat$regid), FUN="sum")
    vaccrates$vacc.wt <- distr.vacc[,2]/distr.sampwt[,2]
    
    # Set up spatial information for plotting
    country.map.points = fortify(country.map, region="id")
    country.map.df = join(country.map.points, country.map@data, by="id")
    
    # Merge District Vaccination Rate data with polygon dataframe
    country.map.df <- merge(country.map.df, vaccrates, by='id', all=T)
    
    
  }else if(method2==TRUE){
    
    # Set up cluster data
    clust.dat <- aggr.dhs.clusters(dhs.data, vaccbreaks=seq(0,1,.1), treat_NAs="unvacc")
    clust.dat <- clust.dat[!is.na(clust.dat$sampwt) & !is.na(clust.dat$long) & !is.na(clust.dat$long),]
    clust.dat$vacc.wt <- clust.dat$vaccrate*clust.dat$sampwt
    if (length(which(clust.dat$lat==0 & clust.dat$long==0)>0)){
      clust.dat <- clust.dat[-which(clust.dat$lat==0 & clust.dat$long==0),]
    }
    
    country.map <- spTransform(country.map, CRS("+init=epsg:4326"))
    clust.dat2 <- SpatialPointsDataFrame(coords=matrix(c(clust.dat$long, clust.dat$lat), ncol=2), 
                                         data=clust.dat, proj4string= CRS("+init=epsg:4326"))
    countryAg.wt <- aggregate(clust.dat2['sampwt'], country.map, FUN="sum")
    countryAg <- aggregate(clust.dat2['vacc.wt'], by=country.map, FUN="sum")
    countryAg$wt <- countryAg.wt$sampwt
    countryAg$vacc.wt <- countryAg$vacc.wt / countryAg$wt
    
    countryAg@data$id = tolower(as.character(country.map@data$DHSREGEN))
    country.map.df = fortify(countryAg, region="id")
    country.map.df = inner_join(country.map.df, countryAg@data, by="id")
  }
  
  #///////////////////////////////////////////////////////////////////////////////////////////
  # Plot the map and data
  
  # Set up color palatte
  myPalette <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"), space="Lab")
  
  # Plot the map and district-level vaccination
  p1 <- ggplot() + geom_polygon(data=country.map.df, 
                                aes(x=long, y=lat,  group=group, fill=vacc.wt)) + coord_equal() + 
    geom_path(data=country.map.df, aes(x=long, y=lat, group=group), color="white", size=1) +
    scale_fill_gradientn(colours = myPalette(100), name='Vaccination\nRate', limits=c(0,1)) 
  
  # Add Clusters colored by vaccination rate
  p1 <- p1 + geom_point(data=clust.dat, aes(x=long, y=lat, fill=vaccrate), 
                        pch=21, alpha=pt.transpar, size=pt.size, colour='grey20') +
    scale_color_gradientn(colours = myPalette(100), name='Vaccination\nRate', limits=c(0,1)) +
    guides(color=guide_colorbar(barwidth=1, barheight=10)) +
    xlab('Longitude') + ylab('Lattitude') + ggtitle(paste0(country, ', ', year, ' DHS'))
  p1 <- p1 + scalebar(country.map.df, location=scale.loc, st.size=4, dist = scale.dist, dd2km = TRUE, model = 'WGS84')
  return(p1)
})








# CALCULATE TAU -----------------------------------------------------------







#############################################################
#   jitter_dhs function                             ####
#   - Jitters the location of DHS clusters based on their status as
#   - urban or rural. For use with tau estimation to smooth taus
#  * Urban jittered to max 2km
#  * Rural to 5km
#  * 1% of Rural jittered up to 10km

jitter_dhs <- function(data, urb.jitt=2, rur.jitt=5, rur.jitt2=10){
  n <- nrow(data)
  urban <- data[,'urban']==TRUE
  n.urb <- sum(urban)
  n.rur <- sum(!urban)
  n.rur10km <- round(.01*n.rur,0)
  n.rur <- n.rur-n.rur10km
  rur10km_ <- sample(1:n, n.rur10km, replace=FALSE, prob=as.numeric(!urban))
  rur10km <- rep(FALSE, n)
  rur10km[rur10km_] <- TRUE
  
  angle <- runif(n, min = 0, max = 2 * pi)
  
  d.urb <- sqrt(runif(n.urb, min=0, max=urb.jitt^2))
  d.rur <- sqrt(runif(n.rur, min=0, max=rur.jitt^2))
  d.rur10km <- sqrt(runif(n.rur10km, min=0, max=rur.jitt2^2))
  
  data[urban, c('x', 'y')] <- data[urban, c('x', 'y')] + cbind(d.urb * cos(angle[urban]), d.urb * sin(angle[urban]))
  data[!urban & !rur10km, c('x', 'y')] <- data[!urban & !rur10km, c('x', 'y')] + cbind(d.rur * cos(angle[!urban & !rur10km]), d.rur * sin(angle[!urban & !rur10km]))
  data[rur10km, c('x', 'y')] <- data[rur10km, c('x', 'y')] + cbind(d.rur10km * cos(angle[rur10km]), d.rur10km * sin(angle[rur10km]))
  
  return(data)
}



jittered_tau_hh <- function(data, typed.tau=FALSE, r=r.max, r.low=r.min, remove.hh=TRUE, boot.iter=1000, bootstrap=TRUE, moving.window=TRUE, max.dist=10, remove_self=FALSE){
  
  data <- data[!is.na(data[,'weight']),]
  
  # if (moving.window){
  #   r.move.max <- sort(c((r[-1]+r[-length(r)])/2, r))
  #   r.move.min <- sort(c((r.low[-1]+r.low[-length(r.low)])/2, r.low))
  #   r <- r.move.max
  #   r.low <- r.move.min
  # }
  
  taus <- matrix(NA, boot.iter, length(r))
  
  for (i in 1:boot.iter){
    dat.tmp <- jitter_dhs(data)
    if (bootstrap) dat.tmp <- dat.tmp[sample(nrow(dat.tmp), replace=T),]
    
    if (typed.tau) {
      if (remove.hh){
        taus[i,] <- get.tau.typed.clustsurvey.hh.wts(posmat=dat.tmp, typeA=1, typeB=2, r=r, r.low=r.low, remove_self=remove_self)
      } else {
        taus[i,] <- get.tau.typed.clustsurvey.wts(posmat=dat.tmp, typeA=1, typeB=2, r=r, r.low=r.low, remove_self=remove_self)
      }
    } else {
      if (remove.hh){
        taus[i,] <- get.tau.clustsurvey.hh.wts(posmat=dat.tmp, r=r, r.low=r.low, remove_self=remove_self)
      } else {
        taus[i,] <- get.tau.clustsurvey.wts(posmat=dat.tmp, r=r, r.low=r.low, remove_self=remove_self)
      }
    }
    
    taus[i,which(taus[i,]==0)] <- NA
    
  }
  return(list(taus=taus, r=r, r.low=r.low))
}




fits.summary <- function(fits, conf.inter=c(.025,.975)){
  fits <- fits[fits$converged,]
  fits <- fits[,-4]
  rc <- matrix(nrow=2,ncol=ncol(fits)); rownames(rc) <- conf.inter
  for (i in 1:ncol(rc)) {
    rc[,i] <- quantile(fits[,i], probs=c(conf.inter[1], conf.inter[2]), na.rm = TRUE)
  }
  colnames(rc) <- colnames(fits)
  
  return(list(means=colMeans(fits, na.rm = T), sd=apply(fits, 2, sd, na.rm=T), medians=apply(fits, 2, median, na.rm = T),
              conf.limits=rc))
}










# Plot the DHS Cluster Empirical Tau with confidence limits#

plot.tau.dhs.clust.boots <- function(tau.dat, r.mid, ylim=NULL, max.dist, conf.inter=c(0.025, 0.975), 
                                     lwd=c(3), col=c('red'), gridln=FALSE, recenter=FALSE){
  
  #recenter empirical
  if (recenter) {
    means <- rowMeans(tau.dat, na.rm=T)
    tau.dat <- tau.dat/means
  }
  
  rc <- matrix(nrow=2,ncol=ncol(tau.dat))
  rownames(rc) <- c(conf.inter[1],conf.inter[2])
  for (i in 1:ncol(rc)) {
    rc[,i] <- quantile(tau.dat[,i], probs=c(conf.inter[1], conf.inter[2]), na.rm = TRUE)
  }
  
  y <- colMeans(tau.dat, na.rm=TRUE)
  
  if (is.null(ylim)) {
    ylim <- c(min(rc), max(rc))
  }
  
  plot(c(r.mid), y, type='l', lty=1, lwd=lwd[1], col=col[1], ylab=expression(paste(tau*(r))), xlab='distance (km)', 
       ylim=ylim, xlim=c(0, max(r.mid)))
  
  if (gridln){
    axis(2, tck=1, col.ticks="light gray")
  }
  
  polygon(c(r.mid,rev(r.mid)), c(rc[1,], rev(rc[2,])), col='grey', border=NA)
  lines(r.mid, colMeans(tau.dat, na.rm=T), type='l', lty=1, lwd=lwd, col=col)
  # lines(c(0,r.mid), y.lower, lty=2, col='green')
  # lines(c(0,r.mid), y.upper, lty=2, col='green')
  abline(h=1, lty=1, col='darkgrey')
}




# Plot the DHS Cluster Empirical Tau with estimated exponential fit and confidence limits for the fit #

plot.tau.dhs.clust.boots2 <- function(tau.dat, fits, lambda=.5, r.mid, ylim=NULL, max.dist, conf.inter=c(0.025, 0.975), 
                                      lwd=c(3,2), col=c('red','blue'), col2='grey', gridln=FALSE, add.plot=FALSE){
  
  rc <- matrix(nrow=2,ncol=ncol(fits))
  rownames(rc) <- c(conf.inter[1],conf.inter[2])
  for (i in 1:ncol(rc)) {
    rc[,i] <- quantile(fits[,i], probs=c(conf.inter[1], conf.inter[2]), na.rm = TRUE)
  }
  colnames(rc) <- colnames(fits)
  
  #recenter empirical
  means <- rowMeans(tau.dat, na.rm=TRUE)
  tau.dat <- tau.dat/means
  
  y <- calc.tau.from.pars(c(0,r.mid), theta=mean(fits[,'theta'], na.rm=T), lambda=lambda, psi=mean(fits[,'psi'], na.rm=T))
  y.lower <- calc.tau.from.pars(c(0,r.mid), theta=rc[1,'theta'], lambda=lambda, psi=rc[1,'psi'])
  y.upper <- calc.tau.from.pars(c(0,r.mid), theta=rc[2,'theta'], lambda=lambda, psi=rc[2,'psi'])
  
  if (is.null(ylim)) {
    ylim <- c(min(min(rc), min(colMeans(tau.dat, na.rm=T))), max(max(rc), max(colMeans(tau.dat, na.rm=T))))
  }
  
  if (add.plot==FALSE){
    plot(c(0,r.mid), y, type='l', lty=1, lwd=lwd[1], col=col[1], ylab=expression(paste(tau*(r))), xlab='distance (km)', 
         ylim=ylim, xlim=c(0, max(r.mid)))
  }
  
  if (gridln){
    axis(2, tck=1, col.ticks="light gray")
  }
  
  
  polygon(c(0,r.mid,rev(r.mid),0), c(y.upper, rev(y.lower)), col=col2, border=NA)
  lines(c(0,r.mid), y, type='l', lty=1, lwd=lwd[1], col=col[1])
  lines(r.mid, colMeans(tau.dat, na.rm=T), type='l', lty=3, lwd=lwd[2], col=col[2])
  # lines(c(0,r.mid), y.lower, lty=2, col='green')
  # lines(c(0,r.mid), y.upper, lty=2, col='green')
  abline(h=1, lty=1, col='darkgrey')
}



# Plot the DHS Cluster Empirical Tau with estimated exponential fit and confidence limits for the fit #

plot.tau.dhs.clust.boots.both <- function(tau.dat, fits, lambda=.5, r.mid, ylim=NULL, max.dist, conf.inter=c(0.025, 0.975), 
                                          lwd=c(3,2), col=c('red','blue'), col2='grey', gridln=FALSE, add.plot=FALSE){
  
  rc <- matrix(nrow=2,ncol=ncol(fits))
  rownames(rc) <- c(conf.inter[1],conf.inter[2])
  for (i in 1:ncol(rc)) {
    rc[,i] <- quantile(fits[,i], probs=c(conf.inter[1], conf.inter[2]), na.rm = TRUE)
  }
  colnames(rc) <- colnames(fits)
  
  #recenter empirical
  means <- rowMeans(tau.dat, na.rm=TRUE)
  tau.dat <- tau.dat/means
  
  y <- calc.tau.from.pars(c(0,r.mid), theta=mean(fits[,'theta'], na.rm=T), lambda=lambda, psi=mean(fits[,'psi'], na.rm=T))
  y.lower <- calc.tau.from.pars(c(0,r.mid), theta=rc[1,'theta'], lambda=lambda, psi=rc[1,'psi'])
  y.upper <- calc.tau.from.pars(c(0,r.mid), theta=rc[2,'theta'], lambda=lambda, psi=rc[2,'psi'])
  
  if (is.null(ylim)) {
    ylim <- c(min(min(rc), min(colMeans(tau.dat, na.rm=T))), max(max(rc), max(colMeans(tau.dat, na.rm=T))))
  }
  
  if (add.plot==FALSE){
    plot(c(0,r.mid), y, type='l', lty=1, lwd=lwd[1], col=col[1], ylab=expression(paste(tau*(r))), xlab='distance (km)', 
         ylim=ylim, xlim=c(0, max(r.mid)))
  }
  
  if (gridln){
    axis(2, tck=1, col.ticks="light gray")
  }
  
  
  polygon(c(0,r.mid,rev(r.mid),0), c(y.upper, rev(y.lower)), col=col2, border=NA)
  lines(c(0,r.mid), y, type='l', lty=1, lwd=lwd[1], col=col[1])
  lines(r.mid, colMeans(tau.dat, na.rm=T), type='l', lty=3, lwd=lwd[2], col=col[2])
  # lines(c(0,r.mid), y.lower, lty=2, col='green')
  # lines(c(0,r.mid), y.upper, lty=2, col='green')
  abline(h=1, lty=1, col='darkgrey')
}














jittered_tau <- function(data, r=r.max, r.low=r.min, boot.iter=1000, bootstrap=TRUE, moving.window=TRUE, max.dist=10, remove_self=FALSE){
  data <- data[!is.na(data[,'weight']),]
  
  # if (moving.window){
  #   r.move.max <- sort(c((r[-1]+r[-length(r)])/2, r))
  #   r.move.min <- sort(c((r.low[-1]+r.low[-length(r.low)])/2, r.low))
  #   r <- r.move.max
  #   r.low <- r.move.min
  # }
  
  taus <- matrix(NA, boot.iter, length(r))
  
  for (i in 1:boot.iter){
    dat.tmp <- jitter_dhs(data)
    if (bootstrap) dat.tmp <- dat.tmp[sample(nrow(dat.tmp), replace=T),]
    taus[i,] <- get.tau.clustsurvey.wts(posmat=dat.tmp, r=r, r.low=r.low, remove_self=remove_self)
    taus[i,which(taus[i,]==0)] <- NA
  }
  return(list(taus=taus, r=r, r.low=r.low))
}



jittered_tau_hh <- function(data, typed.tau=FALSE, r=r.max, r.low=r.min, remove.hh=TRUE, boot.iter=1000, bootstrap=TRUE, moving.window=TRUE, max.dist=10, remove_self=FALSE){
  
  data <- data[!is.na(data[,'weight']),]
  
  # if (moving.window){
  #   r.move.max <- sort(c((r[-1]+r[-length(r)])/2, r))
  #   r.move.min <- sort(c((r.low[-1]+r.low[-length(r.low)])/2, r.low))
  #   r <- r.move.max
  #   r.low <- r.move.min
  # }
  
  taus <- matrix(NA, boot.iter, length(r))
  
  for (i in 1:boot.iter){
    dat.tmp <- jitter_dhs(data)
    if (bootstrap) dat.tmp <- dat.tmp[sample(nrow(dat.tmp), replace=T),]
    
    if (typed.tau) {
      if (remove.hh){
        taus[i,] <- get.tau.typed.clustsurvey.hh.wts(posmat=dat.tmp, typeA=1, typeB=2, r=r, r.low=r.low, remove_self=remove_self)
      } else {
        taus[i,] <- get.tau.typed.clustsurvey.wts(posmat=dat.tmp, typeA=1, typeB=2, r=r, r.low=r.low, remove_self=remove_self)
      }
    } else {
      if (remove.hh){
        taus[i,] <- get.tau.clustsurvey.hh.wts(posmat=dat.tmp, r=r, r.low=r.low, remove_self=remove_self)
      } else {
        taus[i,] <- get.tau.clustsurvey.wts(posmat=dat.tmp, r=r, r.low=r.low, remove_self=remove_self)
      }
    }
    
    taus[i,which(taus[i,]==0)] <- NA
    
  }
  return(list(taus=taus, r=r, r.low=r.low))
}


#/////////////////////////////////////////////////////////////////////////////////////////////////
# Fit Parameters to Bootstraps                                                     #########
# - Using the full set of bootstrap values 

fit.dhs.cluster.boots_V2 <- function(tau.dat, r_lower, r_upper, theta=.75, lambda=.5, psi=1, fit.psi=TRUE, fit.lambda=FALSE){
  
  log.var.tau.dat <- apply(log(tau.dat), 2, var, na.rm=T)
  var.tau.dat <- exp(log.var.tau.dat)
  
  fits <- data.frame(theta=rep(NA, nrow(tau.dat)), lambda=rep(NA, nrow(tau.dat)), psi=psi, converged=TRUE)
  for (i in 1:nrow(tau.dat)){
    
    errorOccurred <- FALSE
    dat.tmp <- data.frame(tau=tau.dat[i,], wt=(1/var.tau.dat), r_lower=r_lower, r_upper=r_upper)
    
    tryCatch({
      fit <- optim(par=c(theta, lambda), fn=min_funct_mean_wt, data=dat.tmp, method = c("Nelder-Mead"))
      fits[i,'theta'] <- fit$par[1]
      fits[i,'psi'] <- psi
      fits[i,'lambda'] <- fit$par[2]
      
    }, error=function(e){
      cat("ERROR : ", conditionMessage(e), "\n")
      errorOccurred <<- TRUE
    })
    if(errorOccurred) fits$converged[i] <- FALSE
  }
  return(fits)
}  



#/////////////////////////////////////////////////////////////////////////////////////////////////
# - Using the full set of bootstrap values 

fit.dhs.cluster.boots_V3 <- function(tau.dat, r_lower, r_upper, theta=.75, lambda=.5, psi=1, fit.psi=TRUE, fit.lambda=FALSE){
  
  log.var.tau.dat <- apply(log(tau.dat), 2, var, na.rm=T)
  var.tau.dat <- exp(log.var.tau.dat)
  fits <- data.frame(theta=rep(NA, nrow(tau.dat)), lambda=rep(NA, nrow(tau.dat)), psi=psi, converged=TRUE)
  
  for (i in 1:nrow(tau.dat)){
    
    errorOccurred <- FALSE
    dat.tmp <- data.frame(tau=tau.dat[i,], wt=(1/var.tau.dat), r_lower=r_lower, r_upper=r_upper)
    
    tryCatch({
      fit <- nls( tau ~ fit_funct2(r_upper, r_lower, theta, lambda), data=dat.tmp, start=list(theta=theta, lambda=lambda), weight=dat.tmp$wt, algorithm = 'port')
      fits.tmp <- coefficients(fit)
      lambda.fit <- fits.tmp['lambda']
      psi.fit <- psi
      
      fits[i,'theta'] <- fits.tmp['theta']
      fits[i,'psi'] <- psi.fit
      fits[i,'lambda'] <- lambda.fit
      
      
    }, error=function(e){
      cat("ERROR : ", conditionMessage(e), "\n")
      errorOccurred <<- TRUE
    })
    if(errorOccurred) fits$converged[i] <- FALSE
  }
  return(fits)
}  



# Fitting Function for fitting taus as means of r slices
fit_funct <- function(r_upper=r.max.new, r_lower=r.min.new, par=c(theta, lambda)) {
  n.r <- length(r_upper)
  sapply(X=1:n.r, FUN=function(x=X) {integrate(calc.tau.from.pars, lower=r_lower[x], upper=r_upper[x], theta=par[1], lambda=par[2], psi=1)$value / (r_upper[x]-r_lower[x])})
}

fit_funct2 <- function(r_upper=r.max.new, r_lower=r.min.new, theta, lambda) {
  n.r <- length(r_upper)
  sapply(X=1:n.r, FUN=function(x=X) {integrate(calc.tau.from.pars, lower=r_lower[x], upper=r_upper[x], theta=theta, lambda=lambda, psi=1)$value / (r_upper[x]-r_lower[x])})
}

# Minimization functions for use with 'fit_funct' and 'optim'
min_funct_mean <- function(par, data){ 
  data <- data[!is.na(data$tau),]
  sum((data$tau - fit_funct(r_upper=data$r_upper, r_lower=data$r_lower, par))^2)
}

min_funct_mean_wt <- function(par, data){ 
  data <- data[!is.na(data$tau),]
  sum(((data$tau - fit_funct(r_upper=data$r_upper, r_lower=data$r_lower, par))^2)*data$wt)
}






# Calculate Tau from available DHS files
### calc.tau.dhs.data   ####
calc.tau.dhs.data <- cmpfun(function(dat.clust, country, data.year, sampwt=TRUE, remove.hh=FALSE,
                                     r.max=c(2,4,6,9,12,16,20,25,30,35,40,50), fit.lambda=FALSE, fit.psi=FALSE,
                                     theta=0.5, lambda=0.5, psi=1, plot.taus=FALSE, boot.iter=10000,
                                     plot.dir=NA, save.plots=TRUE, refit.nonconv=TRUE, plot.fits=TRUE,
                                     blend.first.r=FALSE
                                     ){
  
  # Set Up "r" vectors
  r.min <- c(0, r.max[-length(r.max)])
  r.mid <- (r.max+r.min)/2
  max.dist <- max(r.max)
  
  # Let the second slice include the first slice
  if (blend.first.r) {
    r.min <- c(0, 0, r.max[-length(r.max)][-1])
    r.mid <- (r.max+r.min)/2
  }

  # Calculate Tau w/ Jittering, Moving Window, Bootstraps
  tau.data <- jittered_tau_hh(data=dat.clust, r=r.max, r.low=r.min, remove.hh=remove.hh, boot.iter=boot.iter, bootstrap=TRUE, moving.window=TRUE, max.dist=max.dist, remove_self=remove_self)

  tau <- tau.data$tau
  colnames(tau) <- tau.data$r
  r.mid.new <- (tau.data$r+tau.data$r.low)/2
  
  # Fit Curves to Empirical Taus
  #fits <- fit.dhs.cluster.boots(tau.dat=tau, r.mid=r.mid.new, theta=theta, lambda=lambda, psi=psi, fit.psi=fit.psi, fit.lambda=fit.lambda)
  fits2 <- fit.dhs.cluster.boots_V2(tau.dat=tau, r_lower=tau.data$r.low, r_upper=tau.data$r, theta=theta, lambda=lambda, psi=psi, fit.psi=fit.psi)
  fits <- fit.dhs.cluster.boots_V3(tau.dat=tau, r_lower=tau.data$r.low, r_upper=tau.data$r, theta=theta, lambda=lambda, psi=psi, fit.psi=fit.psi)
  
  fits_summ <- fits.summary(fits=fits, conf.inter=c(.025,.975))
  fits_summ2 <- fits.summary(fits=fits2, conf.inter=c(.025,.975))
  
  # Plot Fitted Parameters
  if(plot.fits) {
    par(mfrow=c(1,2))
    hist(fits[,1], breaks=100, xlab='theta', main='Fitted Thetas')
    hist(fits[,2], breaks=100, xlab='lambda', main='Fitted Lambdas')
    # hist(fits_v2[,1], breaks=100, xlab='theta', main='Fitted Thetas - Means')
    # hist(fits_v2[,2], breaks=100, xlab='lambda', main='Fitted Lambdas - Means')
  }
  
  # Create Variable indicating the Type of calculation
  if (remove.hh & sampwt)  tau.type <- 'wt.hh'
  if (!remove.hh & sampwt)  tau.type <- 'wt'
  if (remove.hh & !sampwt)  tau.type <- 'hh'
  if (!remove.hh & !sampwt)  tau.type <- ''
  
  # Plot resulting tau & Fits
  if (plot.taus){
    
    if (save.plots) {
      dir.create(plot.dir, recursive=T, showWarnings=F)
    }
    
    if (save.plots) {
      png(file.path(plot.dir, paste0(get.iso(country), '_', as.integer(data.year),'_',tau.type,'_emp&fit.png')))
    }
    par(mfrow=c(2,1), mar=c(1,4.5,3,1))
    plot.tau.dhs.clust.boots(tau.dat=tau, r.mid=r.mid.new, ylim=NULL, max.dist=max.dist, conf.inter=c(0.025, 0.975))
    title(paste0('Empirical Tau(r) - ',country,', ', data.year))
    
    par(mar=c(4,4.5,3,1))
    plot.tau.dhs.clust.boots2(tau.dat=tau, fits=fits, lambda=lambda, r.mid=r.mid.new, ylim=NULL, 
                              max.dist=max.dist, conf.inter=c(0.025, 0.975), lwd=c(3,2), col=c('red','blue'))
    legend((max(r.mid.new, na.rm=T)*.5),(max(tau, na.rm=T)*.9),c(expression(paste(fitted," ",tau*(r))), 
                                                                 expression(paste(mean," ",tau*(r)))),col=c('red','blue'), lty=c(1,3), lwd=c(3,2))
    title(paste0('Fitted Tau(r) - ',country,', ', data.year))

    
    if (save.plots) {
      dev.off()
    }
  }
  
  return(list(tau.type=tau.type, dat.clust=dat.clust, tau.data=tau.data, tau.means=colMeans(tau), fits=fits, fits_summ=fits_summ, fits2=fits2, fits_summ2=fits_summ2,
              tau=tau, r.mid.new=r.mid.new, r.max=r.max, r.mid=r.mid, r.min=r.min))
})









