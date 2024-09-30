library(dplyr)
library(caret)

# Load the data
simdata <- read.csv("simdata.csv")

# Aggregate to get unique IDs with their event status
unique_individuals <- simdata %>%
  group_by(ID) %>%
  summarize(Event = max(Event), .groups = 'drop')  # Take max to ensure at least one event is considered

# Separate IDs into those with events and those without events
event_ids <- unique_individuals %>% filter(Event == 1)
non_event_ids <- unique_individuals %>% filter(Event == 0)

# Split event IDs into training and testing (75% training, 25% testing)
set.seed(132)  # Ensure reproducibility
train_event_ids_initial <- createDataPartition(event_ids$ID, p = 0.75, list = FALSE)  # 75% for initial training

# IDs for initial training and testing
train_event_ids_initial <- event_ids$ID[train_event_ids_initial]
test_event_ids <- event_ids$ID[!(event_ids$ID %in% train_event_ids_initial)]

# Combine initial training IDs with all non-event IDs
train_ids_initial <- c(train_event_ids_initial, non_event_ids$ID)
test_ids <- c(test_event_ids, setdiff(unique_individuals$ID, train_ids_initial))

# Split initial training data into final training (75%) and validation (25%)
train_ids_final <- createDataPartition(train_ids_initial, p = 0.75, list = FALSE)

# Final training and validation IDs
train_ids <- train_ids_initial[train_ids_final]
val_ids <- train_ids_initial[!(train_ids_initial %in% train_ids)]

# Create training, validation, and test datasets
training_data <- simdata %>% filter(ID %in% train_ids)
validation_data <- simdata %>% filter(ID %in% val_ids)
test_data <- simdata %>% filter(ID %in% test_ids)

# Check proportions of event in train, validation, and test data
 sum(training_data$Event == 1) / sum(simdata$Event == 1)
 sum(validation_data$Event == 1) / sum(simdata$Event == 1)
 sum(test_data$Event == 1) / sum(simdata$Event == 1)
 
 training_data$Event <- as.numeric(training_data$Event)
 training_data.id <- training_data %>%
   group_by(ID) %>%
   filter(row_number() == n()) %>%
   select(ID, X1, X2, ttp, Event)
 
 test_data.id <- test_data %>%
   group_by(ID) %>%
   filter(row_number() == n()) %>%
   select(ID, X1, X2, ttp, Event)
 
 validation_data.id <- validation_data %>%
   group_by(ID) %>%
   filter(row_number() == n()) %>%
   select(ID, X1, X2, ttp, Event)
 

 write.csv(training_data,file = "training_data.csv", row.names = FALSE)
 write.csv(validation_data,file = "validation_data.csv", row.names = FALSE)
 write.csv(test_data,file = "test_data.csv", row.names = FALSE)
 
