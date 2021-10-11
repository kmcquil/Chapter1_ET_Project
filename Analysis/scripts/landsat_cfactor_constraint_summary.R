###########################################################################
## Look at which scenes ET was not calculated due to c-factor constraints 
## How many were not calculated out of total scenes downloaded 
## Which tiles out of how many unique tiles 
## What time of year 
###########################################################################
library(data.table)

all_scenes <- fread("G:/My Drive/Chapter1_ET_Project/Data/Landsat_ET/downloading/all_scenes.txt", 
                    header = F, 
                    sep = "|")
colnames(all_scenes) <- "Name"

# create a column for the tile, the date, and if it was able or not to generate ET 
all_scenes$Tile <- substr(all_scenes$Name, 11, 16)
all_scenes$Date <- as.Date(substr(all_scenes$Name, 18, 25), format = "%Y%m%d")
all_scenes$Calcualted <- ifelse(grepl("C Factor Constraints", all_scenes$Name), "F", "T")
all_scenes$Month <- month(all_scenes$Date)

notCalculated <- all_scenes[all_scenes$Calcualted == "F",]

# percent of total scenes ET could not be estiamted  = 14.76
perc_not_calculated <- (nrow(notCalculated)/nrow(all_scenes))*100

# list all unique tiles 
all_unique_tiles <-unique(all_scenes$Tile)
noCalc_unique_tiles <- unique(notCalculated$Tile)
notCalculated[,.N, Tile]

# common dates that ET couldn't be generated 
noCalc_byMonth <- notCalculated[,.N, Month]


# how many come from growing season vs non growing season months 
dormant_count <- sum(noCalc_byMonth[Month %in% c(10, 11, 12, 1, 2, 3),]$N)
dormant_perc <- (dormant_count/nrow(notCalculated))*100

gs_count <- sum(noCalc_byMonth[Month %in% c(4,5,6,7,8,9),]$N)
gs_perc <- (gs_count/nrow(notCalculated))*100

# how many of growing season scenes could not be generated 
(gs_count/nrow(all_scenes[Month %in% c(4,5,6,7,8,9),]))*100

# which tiles have a c factor constraint during the growing season 
notCalculated[Month %in% c(4,5,6,7,8,9), .N, Tile]



