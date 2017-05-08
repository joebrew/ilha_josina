# Clusters of 50 children, with 1 km buffers (2 kms between core zones)
# How many clusters can we create?
# agnostic to village boundaries, etc.

# Libraries
library(tidyverse)
library(RColorBrewer)
library(cism)
library(sp)

# Algorithm for creating clusters
source('../../lib/cluster_optimize.R')

# Get data from openhds
if('open_hds_data.RData' %in% dir()){
  load('open_hds_data.RData')
} else {
  membership <- 
    cism::get_data(tab = 'membership',
                   dbname = 'openhds')
  individual <- 
    cism::get_data(tab = 'individual',
                   dbname = 'openhds')
  location <- 
    cism::get_data(tab = 'location',
                   dbname = 'openhds')
  residency <-
    cism::get_data(tab = 'residency',
                   dbname = 'openhds')
  VISIT_REGISTRATION_CORE <-
    cism::get_data(tab = 'VISIT_REGISTRATION_CORE',
                   dbname = 'dssodk')
  
  save(membership,
       individual,
       location,
       residency,
       VISIT_REGISTRATION_CORE,
       file = 'open_hds_data.RData')
}
# Clean up -----------------------------------

# Remove the extra characters in invdividual.extId
individual$extId <- substr(individual$extId,
                           start = 1,
                           stop = 9)

# Remove those with hh in permid, and by permid we mean lastname
individual <- individual %>%
  filter(!grepl('hh', tolower(lastName)))

# Make data objects

residency <- residency %>%
  mutate(startDate = as.Date(startDate),
         endDate = as.Date(endDate, origin = '1970-01-01'))
individual$dob <- as.Date(individual$dob)

# We're going to snapshot on. So, remove
# those observations that come before/after, etc.
snap_shot <- as.Date('2016-05-06')
residency <- residency %>%
  mutate(endDate = ifelse(is.na(endDate), snap_shot, endDate)) %>%
  filter(startDate <= snap_shot,
         endDate >= snap_shot)

# Keep only those people as of the snap_shot date
people <- residency %>%
  dplyr::select(individual_uuid,
                location_uuid) %>%
  left_join(individual %>%
              dplyr::select(#extId,
                            uuid,
                            dob,
                            firstName,
                            gender,
                            lastName,
                            middleName,
                            gender),
            by = c('individual_uuid' = 'uuid')) %>%
  left_join(location %>%
              dplyr::select(extId,
                            uuid,
                            latitude,
                            locationName,
                            longitude),
            by = c('location_uuid' = 'uuid'))

people$longitude <- as.numeric(as.character(people$longitude))
people$latitude <- as.numeric(as.character(people$latitude))

# Join to visit registration core from dssodk
people <- people %>%
left_join(VISIT_REGISTRATION_CORE %>%
          dplyr::select(LOCATION_NAME,
                        COORDINATES_LAT,
                        COORDINATES_LNG) %>%
            filter(!is.na(LOCATION_NAME),
                   !is.na(COORDINATES_LAT),
                   !is.na(COORDINATES_LNG)) %>%
            filter(!duplicated(LOCATION_NAME)),
          by = c('locationName' = 'LOCATION_NAME'))
people$latitude <- 
  ifelse(is.na(people$latitude), people$COORDINATES_LAT, people$latitude)
people$longitude <- 
  ifelse(is.na(people$longitude), people$COORDINATES_LNG, people$longitude)
# Get location code
people$is_ij <- substr(people$locationName, 1, 2) == '33'

# Get ilha josina
ij <- man3
ij <- ij[ij$NAME_3 == 'Ilha Josina Machel',]

# Make spatial loc
people_sp <- people %>% filter(!is.na(longitude))
coordinates(people_sp) <- ~longitude+latitude
proj4string(people_sp) <- proj4string(man3)
people_sp$ij <- !is.na(over(people_sp, polygons(ij)))
ij_points <- people_sp[people_sp$ij,]
people <- people %>%
  left_join(ij_points@data)
people$ij[is.na(people$ij)] <- FALSE
# Keep only ilha josina
people <- people %>% filter(ij)
# people <- people %>% filter(is_ij)
# people <- people %>% filter(!is.na(longitude), !is.na(latitude))

census <- people

# Keep only those under 5 at time of snapshot
people <- people %>% 
  mutate(age_years = as.numeric(snap_shot - dob) / 365.25) %>%
  filter(age_years <= 5)

# Create clusters, starting with the southern most point

# Get children in the format needed for clustering
ij_children <- people %>%
  mutate(id = individual_uuid,
         lat = latitude,
         lon = longitude) %>%
  dplyr::select(id,
                lat,
                lon) %>%
  mutate(x = lon,
         y = lat)
coordinates(ij_children) <- ~x+y
proj4string(ij_children) <- proj4string(ij)


# # Loop
# results_list <- list()
# counter <- 0
# starts <- 'far' #c("far", "close", "random")
# rests <- 'close' # c("far", "close", "random")
# walks <- c(TRUE, FALSE)
# for (ss in starts){
#   for(rr in rests){
#     for(ww in walks){
#       mm <- ifelse(ss == 'random' |
#                      rr == 'random',
#                    10,
#                    1)
#       for (i in 1:mm){
#         counter <- counter + 1
#         message(counter)
#         x = cluster_optimize(cluster_size = 50,
#                              plot_map = FALSE,
#                              locations = ij_children,
#                              shp = ij,
#                              sleep = 0,
#                              start = ss,
#                              rest = rr,
#                              messaging = FALSE,
#                              buffer = 1000,
#                              walk_along = ww)
#         x$start <- ss
#         x$rest <- rr
#         x$walk_along <- ww
#         x$iter <- counter
#         results_list[[counter]] <- x
#       }
#     }
#   }
# }
# results <- bind_rows(results_list)
# results$iter <- rep(1:(nrow(results) / nrow(ij_children)),
#               each = nrow(ij_children))

# Just do one
walk_along <- FALSE
start = 'close'
rest = 'far'
results = cluster_optimize(cluster_size = 50,
                     plot_map = FALSE,
                     locations = ij_children,
                     shp = ij,
                     sleep = 0,
                     start = start,
                     rest = rest,
                     messaging = FALSE,
                     buffer = 1000,
                     walk_along = walk_along) %>%
  mutate(iter = 1,
         start = start,
         rest = rest,
         walk_along = walk_along)
table(results$cluster)
