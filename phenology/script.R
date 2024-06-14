##### you may need to install some packages first
library(rnaturalearth)
library(lubridate)
library(tidyverse)
library(sf)
library(readr)
library(phenesse)

#set working directory
setwd("C:/Users/robgu/Downloads")

##################### read in data and do some initial filtering
# view the full file below (gistfile2.txt) and then save_as "allDistinctLepObsSubset.csv" to your working directory in R.
# this is occurrence data for a small subset of butterfly species 
leps <- read_csv("allDistinctLepObsSubset.csv")

#lets view this file and take a look at dates
View(leps)

#filter records for bad years
leps <- leps %>%
  filter(year>1867,year<2022)

#remove NA eventDates 
leps <- leps %>%
  drop_na(eventDate)

#Get day of year using yday and add as a column using mutate
#yday comes from the excellent lubridate package
leps <- leps %>% 
  mutate(doy = yday(lubridate::as_date(eventDate)))

#drop any remnant bad dates that didn't generate a DOY
leps <- leps %>% drop_na(doy)

########################################### make a grid across North America
#create sf object of North America using rnaturalearth
test2 <- rnaturalearth::ne_countries(country = c("United States of America", "Mexico", "Canada"),returnclass = "sf")

#set projection
test2<-st_transform(test2, crs=
                      "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#make the grid and deal withh some bookkeeping 
grids = st_make_grid(test2, cellsize = c(250000, 250000))
grids2 = grids[test2]

#number each grid cell
grids2 = mutate(st_sf(geometry = grids2), id_cells = 1:n())

#let's make sure it looks ok
plot(grids2)

############################## Lets assemble obervations to grids
#remove blank decimal lat and long 
#and tell R which columns have lat lons so these are now geometries
leps_sf <- leps %>%
  filter(!is.na(decimalLongitude),
         !is.na(decimalLatitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = "+proj=longlat +datum=WGS84 +no_defs")

#lets_drop_lots_of_duplicates
leps_sf <- distinct(leps_sf, locality, year, doy,recordedBy, .keep_all= TRUE)


#set the same projection for points as for grids
leps_sf <- st_transform(leps_sf,crs =st_crs(grids))

#join the points to the grids - so easy!
leps_sf <- st_join(leps_sf, grids2)

#we can now drop the geometries since we have a grid cell assigned
leps_df <- st_drop_geometry(leps_sf)



#lets just check and see what we have now
View(leps_df)

#write_out_your_final_file_for_later
write.csv(x = leps_df, file = "LepsByGrid.csv", row.names = F)

######################################################
###Lets assemble some phenometrics!
##first we will get a count of observations per each year x sciName x id_cells combination

leps_df_count <- leps_df %>%
  group_by(year, sciName,id_cells) %>%
  summarise(obsnum = n()) 

#lets_next_run_a_fast_phenometric_with_boostrapping_for_CIs
#after_we_do_some_cleanup filtering
#filters include: remvoing december and jan flight recs
#filtering for at least 12 observations per cell
#filter any cases where we have no variation in DOY
leps_flight_10perc_onset <-  leps_df %>% 
  group_by(sciName,year,id_cells) %>% 
  filter(doy > 30 & doy < 330 ) %>%
  filter(n() > 12) %>%
  filter(n_distinct(doy)>1) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$doy, percentile = 0.10, bootstraps=500)))

#Let's combine counts and phenometrics into one output
leps_flight_10perc_onset_wcount <- merge(leps_flight_10perc_onset,leps_df_count, by=c("year", "id_cells", "sciName"))

#How does it look?
View(leps_flight_10perc_onset_wcount)
