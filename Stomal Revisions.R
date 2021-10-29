library(readxl)
stomalcomps <- read_excel("~/Desktop/stomalcomps.xlsx")
View(stomalcomps)

##Packages
library(lubridate)
library(tidyverse)

###survival curve of stomal revision based on different facts 
colnames(stomalcomps)

##create tables for all instruments 
blank <- subset(stomalcomps, stomalcomps$redcap_repeat_instrument == "NA" | is.na(redcap_repeat_instrument))
View(blank)

initialsurgery <- subset(stomalcomps, stomalcomps$redcap_repeat_instrument == "initial_reconstruction")
View(initialsurgery)

other_surgeries <- subset(stomalcomps, stomalcomps$redcap_repeat_instrument == "subsequent_surgeries")
View(other_surgeries)

fu <- subset(stomalcomps, stomalcomps$redcap_repeat_instrument == "recent_fu")

#Create table with only variables I want (DOES NOT YET INCLUDE SURGERIES)
temp1.0<- blank[,c(1, 6, 7, 14, 15, 18, 19, 36)]
View(temp1.0)
temp1.1<- initialsurgery[,c(1, 97, 114,117, 153, 154, 155, 156, 157, 159, 161, 166, 167, 168, 172, 179)]
View(temp1.1)
temp1.2<- other_surgeries[,c(1,399,406, 408, 448, 449)]
temp1.3 <- fu[,c(1, 462)]

StomalRev <- merge(temp1.0, temp1.1, by = "mrn")
StomalRev <- merge(StomalRev, temp1.3, by = "mrn")

View(StomalRev)

##view structure of date 
str(StomalRev)
summary(StomalRev)

StomalRev$mrn <- as.character(StomalRev$mrn)
factor_cols <- c("sex","init_prior_surg", "preop_eval_abdsx___3", "preop_eval_abdsx___4", "primary_dx")
StomalRev[factor_cols] <- lapply(StomalRev[factor_cols], as.numeric)

colnames(StomalRev)

###Change to factors and other correct atomic classes

StomalRev <- StomalRev %>% mutate_at(vars(10:23), as.factor)

StomalRev$bmi <- as.numeric(StomalRev$bmi)
StomalRev$add_sx_info___6 <- as.numeric(StomalRev$add_sx_info___6)

StomalRev <- StomalRev %>% mutate_at(vars(12:15, 17:23), as.numeric)

###change dates 
StomalRev <- StomalRev %>% mutate_at(vars(2, 9, 24), ymd)

###clean INF from BMI 
StomalRev$bmi[is.infinite(StomalRev$bmi)] <- NA
StomalRev$bmi

##Boxplots 
boxplot(StomalRev$bmi)
hist(StomalRev$bmi)
plot(StomalRev$bmi, main = "BMI of Augmented Stomal Group")

##Age at surgery 

StomalRev$ageatsurgery <- (difftime(StomalRev$date_surgery, StomalRev$dob, units = "secs"))/31536000
View(StomalRev)
boxplot(StomalRev$ageatsurgery)
summary(StomalRev$ageatsurgery)
StomalRev$ageatsurgery <- as.numeric(StomalRev$ageatsurgery)

##Pearson chi square with multiple variables (had to change all the involved variables for numericals from factors)
colnames(StomalRev)
CHISPREOP <- lapply(StomalRev[,c(3,5,7,8)], function(x) cor.test(StomalRev[,11], x, method = "pearson"))
CHISPREOP

do.call(rbind, CHISPREOP)[,c(1:3)]


##Perform Pearson chi square in operative variables
lapply(StomalRev[,c(10:23)], function(x) levels(x))

CHISOR <- lapply(StomalRev[,c(12:15, 17:23)], function(x) cor.test(StomalRev[,11], x, method = "pearson"))
do.call(rbind, CHISOR)[,c(1:3)]

###Continuous variables 
WILCOXPREOP <- lapply(StomalRev[,c(4, 25)], function(x) wilcox.test(x~StomalRev[,11], na.rm=T))
wilcox.test(StomalRev$bmi ~ StomalRev$add_sx_info___6)
do.call(rbind, WILCOXPREOP)[,c(1:3)]

##Setting up Survival Data with actual revisions
colnames(other_surgeries)
test1 <- other_surgeries[,c(1, 399, 406)] 
class(test1$subs_sx)
test1$subs_sx <- ymd(test1$subs_sx)
View(test1)

test2 <- filter(test1, subs_sx_info___5 == "1")
colnames(test2)

test2 <- test2 %>% group_by(mrn) %>% filter(subs_sx == min(subs_sx)) %>% slice(1)
View(test2)
test2$DOE <- test2$subs_sx

##Setting up Survival Data with nonrevisions by choosing the last time events were recorded. 
##need to combine fu date first
colnames(fu)

fudate <- fu[,c(1,462)]
fudate
test3 <- merge(test1, fudate, by = "mrn")
test3 <- filter(test3, subs_sx_info___5 == "0")
view(test3)

##remove uncensored mrns
test3$mrn %in% test2$mrn
test3 <- test3[!test3$mrn %in% test2$mrn, ]

test3$mrn %in% test2$mrn

##get single unique mrns for censored data. This path will get you 
##the last event row, but its not necessary in this case since 
##censored event is date of follow up.Also it gets rid of events without dates.
class(test3$subs_sx)
test3$subs_sx <- ymd(test3$subs_sx)
View(test3)
View(test3 %>% group_by(mrn) %>% filter(subs_sx == max(subs_sx)) %>% slice(1))

test3 <- test3[!duplicated(test3$mrn),]
colnames(test3)
test3$DOE <- test3$fu_date
test3 <- test3[-c(4)]

test3$DOE <- ymd(test3$DOE)

##combine tables

test4 <- rbind(test2, test3)
View(test4)

##Combine surgery events with original data
RD <- merge(StomalRev, test4, by = "mrn")
View(RD)

##now run survival graph
install.packages("survival")
library(survival)
library(survminer)

##create new column with time diff 
RD$time_diff <- (difftime(RD$DOE, RD$date_surgery, units="days"))
RD$time_diff

class(RD$time_diff)

##survival curve
colnames(RD)
surv_object <- Surv(time=RD$time_diff, event = RD$subs_sx_info___5)
surv_object

fit1<-survfit(surv_object~add_sx_info___6, data = RD)
summary(fit1)

##need to update R base to use survminer

ggsurvplot(fit1, data = RD, pval=TRUE)

##Fit a Cox proportional hazards model
fit.coxph <- coxph(surv_object~ add_sx_info___6, data = RD)
ggforest(fit.coxph, data = RD)
summary(fit.coxph)

##Survival based on channel type
fit2<-survfit(surv_object~type_conduit, data = RD)
summary(fit2)
ggsurvplot(fit2, data = RD, pval=TRUE)
