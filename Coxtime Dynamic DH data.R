library(dplyr)

########## Preparing data for coxtime model ####################################
coxdata <- simdata[,c(2,3,4,5,12:14,6,8,10)]
names(coxdata)[names(coxdata) == "time"] <- "duration"
names(coxdata)[names(coxdata) == "Event"] <- "event"
write.csv(coxdata, file = "coxdata.csv", row.names = FALSE)

#split
ctest_data <- test_data[,c(2,3,4,5,12:14,6,8,10)]
names(ctest_data)[names(ctest_data) == "time"] <- "duration"
names(ctest_data)[names(ctest_data) == "Event"] <- "event"
write.csv(ctest_data, file = "ctest_data.csv", row.names = FALSE)

ctrain_data <- training_data[,c(2,3,4,5,12:14,6,8,10)]
names(ctrain_data)[names(ctrain_data) == "time"] <- "duration"
names(ctrain_data)[names(ctrain_data) == "Event"] <- "event"
write.csv(ctrain_data, file = "ctrain_data.csv", row.names = FALSE)

cvalidation_data <- validation_data[,c(2,3,4,5,12:14,6,8,10)]
names(cvalidation_data)[names(cvalidation_data) == "time"] <- "duration"
names(cvalidation_data)[names(cvalidation_data) == "Event"] <- "event"
write.csv(cvalidation_data, file = "cvalidation_data.csv", row.names = FALSE)


######### Preparing data for dynamic deephit model #############################
label <- aggregate(Event ~ ID, data = simdata,FUN = max)
ddh_data <- merge(simdata,label, by = "ID")
names(ddh_data)[names(ddh_data) == "Event.x"] <- "fstatus"
names(ddh_data)[names(ddh_data) == "Event.y"] <- "label"
names(ddh_data)[names(ddh_data) == "ttp"] <- "tte"
names(ddh_data)[names(ddh_data) == "time"] <- "times"
names(ddh_data)[names(ddh_data) == "ID"] <- "id"
ddh_data <- ddh_data %>% 
select(id,  tte, times, label, X1, X2, X3, X4, X5.1, X5.2, X5.3, X6)
write.csv(ddh_data,file = "ddeephit_data.csv", row.names = FALSE)
