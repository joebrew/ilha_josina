# Clusters of 50 children, with 1 km buffers (2 kms between core zones)
# How many clusters can we create?
# agnostic to village boundaries, etc.

# Libraries
library(tidyverse)
library(RColorBrewer)
library(cism)
library(sp)
library(rgdal)
library(rgeos)

# Data preparation
if('prepared_data.RData' %in% dir('../data')){
  load('../data/prepared_data.RData')
} else {
  # Algorithm for creating clusters
  source('cluster_optimize.R')
  
  # Get data from openhds
  if('open_hds_data.RData' %in% dir('../data')){
    load('../data/open_hds_data.RData')
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
         file = '../data/open_hds_data.RData')
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
  
  # Make ij wider by 2 kilometers to allow those on the border
  ij_exact <- ij
  ij <- spTransform(ij,
                    CRS("+init=epsg:3347"))
  
  ij <- gBuffer(ij, byid = TRUE, width = 3000)
  # convert back to lat/lon
  ij <- spTransform(ij, CRS(proj4string(man3)))
  
  # Make spatial loc
  people_sp <- people %>% filter(!is.na(longitude))
  coordinates(people_sp) <- ~longitude+latitude
  proj4string(people_sp) <- proj4string(man3)
  people_sp$ij <- !is.na(over(people_sp, polygons(ij_exact)))
  ij_points <- people_sp[people_sp$ij,]
  people <- people %>%
    left_join(ij_points@data)
  people$ij[is.na(people$ij)] <- FALSE
  # Keep only ilha josina
  people <- people %>% filter(ij)
  # people <- people %>% filter(is_ij)
  # people <- people %>% filter(!is.na(longitude), !is.na(latitude))
  
  people <- people %>%
    mutate(age_years = as.numeric(snap_shot - dob) / 365.25) 
  
  # Fix the last name / permid naming issue
  people <- people %>%
    rename(perm_id = lastName)
  
  # Get a household id
  people <- people %>%
    mutate(household_id = substr(x = perm_id, 
                                 start = 1, 
                                 stop = 8))
  
  # Remove duplicates
  people <- people %>%
    filter(!duplicated(perm_id))
  
  # Get a copy of the full census
  census <- people
  
  # Keep only those under 5 at time of snapshot
  people <- people %>% 
    filter(age_years <= 5)
    # filter(age_years < 1)
  
  # Create clusters, starting with the southern most point
  
  # Get children in the format needed for clustering
  ij_children <- people %>%
    mutate(id = individual_uuid,
           lat = latitude,
           lon = longitude) %>%
    dplyr::select(id,
                  perm_id,
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
  start = 'close' # 1 5 change back to close for 5
  rest = 'far'
  results = cluster_optimize(cluster_size = 52,
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
  
  # Create an id in census
  census$id <- census$individual_uuid
  
  # Get results of our simulations
  x <- results %>%
    group_by(iter) %>%
    summarise(start = first(start),
              rest = first(rest),
              walk_along = first(walk_along),
              n = n(),
              ib = length(which(!is.na(cluster)))) %>%
    mutate(pib = ib / n * 100) %>%
    ungroup %>%
    mutate(p_in_core = 100 - pib) %>%
    mutate(n_in_core = n - ib)
  
  # Get that which yields the most number of successes
  best <- x %>% filter(pib == min(pib))
  best <- results[results$iter == best$iter[1],]
  x <- best
  results <- x
  
  # Create centroids  
  centroids <- best %>%
    filter(!is.na(cluster)) %>%
    group_by(cluster) %>%
    summarise(lon = mean(lon),
              lat = mean(lat))
  
  # Assignation to intervention
  
  # Create a dataframe of clusters
  clusters <- results %>%
    group_by(cluster) %>%
    summarise(x_centroid = mean(x),
              y_centroid = mean(y)) %>%
    filter(!is.na(cluster))
  
  # Set a random seed
  set.seed(510)
  
  # Create assignations
  clusters$status <- sample(c('intervention',
                              'control'),
                            nrow(clusters), 
                            replace = TRUE)
  
  # Voronoi tesselation
  
  # Make a spatial version of results
  results_sp <- results %>% filter(!is.na(cluster))
  coordinates(results_sp) <- ~x+y
  
  # Create convex hulls
  ch <- list()
  for (i in 1:nrow(clusters)){
    this_cluster <- clusters$cluster[i]
    sub_results_sp <- results_sp[results_sp@data$cluster == this_cluster,]
    x <- rgeos::gConvexHull(sub_results_sp)
    ch[[i]] <- x
  }
  
  # Create delaunay triangulation / voronoi tiles for entire surface
  voronoi <- function(shp = results_sp){
    
    # Remove any identical ones
    shp <- shp[!duplicated(shp$lon, shp$lat),]
    
    # Helper function to create coronoi polygons (tesselation, not delaunay triangles)
    # http://carsonfarmer.com/2009/09/voronoi-polygons-with-r/
    voronoipolygons = function(layer) {
      require(deldir)
      crds = layer@coords
      z = deldir(crds[,1], crds[,2])
      w = tile.list(z)
      polys = vector(mode='list', length=length(w))
      require(sp)
      for (i in seq(along=polys)) {
        pcrds = cbind(w[[i]]$x, w[[i]]$y)
        pcrds = rbind(pcrds, pcrds[1,])
        polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
      }
      SP = SpatialPolygons(polys)
      voronoi = SpatialPolygonsDataFrame(SP, data=data.frame(x=crds[,1], 
                                                             y=crds[,2], row.names=sapply(slot(SP, 'polygons'), 
                                                                                          function(x) slot(x, 'ID'))))
    }
    # http://gis.stackexchange.com/questions/180682/merge-a-list-of-spatial-polygon-objects-in-r
    appendSpatialPolygons <- function(x) {
      ## loop over list of polygons
      for (i in 2:length(x)) {
        # create initial output polygon
        if (i == 2) {
          out <- maptools::spRbind(x[[i-1]], x[[i]])
          # append all following polygons to output polygon  
        } else {
          out <- maptools::spRbind(out, x[[i]])
        }
      }
      return(out)
    }
    
    tile_polys <- voronoipolygons(shp)
    # Add the cluster numbers
    tile_polys@data$cluster <- the_clusters <- shp$cluster
    cols <- rainbow(as.numeric(factor(tile_polys@data$cluster)))
    
    # Disolve borders
    x = gUnaryUnion(tile_polys, id = tile_polys$cluster)
    
    jdata = SpatialPolygonsDataFrame(Sr=x, 
                                     data=data.frame(cluster = as.numeric(as.character(names(x)))),FALSE)
    
    return(jdata)
  }
  
  # Get voronoi tesselations
  results_spv <- voronoi(shp = results_sp)
  
  # # Plot in current state
  # plot(results_spv)
  # plot(ij, add = T, col = adjustcolor('blue', alpha.f = 0.2))
  
  # Narrow down so as to only keep those areas which are IN Magude
  proj4string(results_spv) <- proj4string(ij)
  out <- gIntersection(results_spv, ij, byid=TRUE)
  
  # Join with data
  row.names(out) <- as.character(1:length(out))
  out <- SpatialPolygonsDataFrame(out, data.frame(clusters), match.ID = FALSE)
  
  # # Plot again
  # plot(ij, col = adjustcolor('red', alpha.f = 0.2))
  # plot(out, add = T, lwd = 0.2)
  
  
  out@data$color <- ifelse(out@data$status == "control", "green",
                           ifelse(out@data$status == "intervention", "red", "grey"))
  # # Plot leaflet
  # library(leaflet)
  # leaflet() %>%
  #   addProviderTiles(provider = 'Stamen.Toner') %>%
  #   addPolygons(data = out,
  #               color = out@data$color,
  #               popup = paste0('Cluster ', as.character(out@data$cluster),
  #                              '. Status: ', as.character(out@data$status)))
  
  # Defining intervention zones.
  
  # Having expanded borders, we can re-calculate buffers based not on cluster number but on intervention status (since there is no need for a buffer between two clusters of identical intervention status).
  # 
  # To do this, we first "dissolve" our polygons by intervention status. This means removing boundaries between areas of the same intervention status. This results in the following intrevention zones (red = intervention; green = control).
  
  # Disolve borders
  collapsed <- gUnaryUnion(out, id = out@data$status)
  
  collapsed <- SpatialPolygonsDataFrame(collapsed, data.frame(zone = 1:length(collapsed@polygons)), 
                                        match.ID = FALSE)
  
  # Get the intervention status
  tester <- clusters %>% filter(!duplicated(status))
  coordinates(tester) <- ~x_centroid+y_centroid
  proj4string(tester) <- proj4string(collapsed)
  x <- over(tester, polygons(collapsed))
  tester$zone <- collapsed@data$zone[x]
  tester <- tester@data %>%
    dplyr::select(zone, status)
  
  # Join the intervention status to collapsed
  collapsed@data <-
    left_join(x = collapsed@data,
              y = tester,
              by = 'zone')
  
  # Add a color
  collapsed@data$color <- ifelse(collapsed@data$status == 'control',
                                 'green',
                                 ifelse(collapsed@data$status == 'intervention',
                                        'red', 'grey'))
  
  # leaflet() %>%
  #   addProviderTiles(provider = 'Stamen.Toner') %>%
  #   addPolygons(data = collapsed,
  #               color = collapsed@data$color)
  
  # We then use spatial buffering to re-define our buffers, going inwards 1,000 meters from the "dissolved" boundaries.
  # 
  
  # Transform to projected form
  collapsed_projected <- spTransform(collapsed, #CRS(proj4string(man3))
                                     CRS( "+init=epsg:3347" ) )
  # Get buffer at 1km
  collapsed_projected_buffered <- rgeos::gBuffer(collapsed_projected, 
                                                 width = -1000,
                                                 byid = TRUE)
  # Un -project collapsed_projected_buffered
  collapsed_buffered <- spTransform(collapsed_projected_buffered, proj4string(ij))
  # leaflet() %>%
  #   addProviderTiles(provider = 'Stamen.Toner') %>%
  #   addPolygons(data = collapsed_buffered,
  #               color = collapsed_buffered@data$color)
  
  
  ## Data products
  
  # Having finished our implementation of cluster delineation, we write software-agnostic shapefiles. The purpose of this step is to remove the data products (cluster and buffer delineations) from this report, so that others can use these clusters, regardless of their GIS software.
  
  # We create 3 data products.
  
  ### 1. `ij_clusters`
  
  # plot(out,
  #      col = rainbow(nrow(out)))
  library(rgdal)
  owd <- getwd()
  setwd('../data_products/')
  writeOGR(obj=out, 
           dsn="ij_clusters", 
           layer = 'ij_clusters',
           driver="ESRI Shapefile",
           overwrite_layer = TRUE)
  setwd(owd)
  ij_clusters <- out
  rm(out)
  
  ### 2. `ij_intervention_zones`
  # `ij_intervention_zones` is a projected shapefile with the exact boundaries of the intervention zones. 
  
  owd <- getwd()
  setwd('../data_products/')
  writeOGR(obj=collapsed, 
           dsn="ij_intervention_zones", 
           layer = 'ij_intervention_zones',
           driver="ESRI Shapefile",
           overwrite_layer = TRUE)
  setwd(owd)
  ij_intervention_zones <- collapsed
  rm(collapsed)
  
  ### 3. `ij_intervention_zones_buffered`
  
  # `ij_intervention_zones_buffered` is a projected shapefile with the exact boundaries of the intervention zones with a -1000 meter inward buffer added to each zone. 
  
  owd <- getwd()
  setwd('../data_products/')
  writeOGR(obj=collapsed_buffered, 
           dsn="ij_intervention_zones_buffered", 
           layer = 'ij_intervention_zones_buffered',
           driver="ESRI Shapefile",
           overwrite_layer = TRUE)
  setwd(owd)
  ij_intervention_zones_buffered <- collapsed_buffered
  rm(collapsed_buffered)
  
  # Increased catchment
  
  # Get zones
  results_sp <- results
  coordinates(results_sp) <- ~x+y
  proj4string(results_sp) <- proj4string(ij)
  x <- over(results_sp, ij_intervention_zones_buffered)
  results_sp@data <- cbind(results_sp@data, x)
  
  # See which ones now get a zone
  results_sp@data <-
    results_sp@data %>%
    mutate(new_zone = is.na(cluster) & !is.na(zone))
  
  results_sp@data$color <- 
    ifelse(is.na(results_sp@data$color),
           'blue',
           results_sp@data$color)
  # plot(ij)
  # points(results_sp, 
  #        col = results_sp@data$color,
  #        pch = '.')
  # 
  new_zone <- results_sp[which(results_sp@data$new_zone),]
  # points(new_zone, 
  #        col = new_zone@data$color,
  #        pch  = 2,
  #        cex = 0.5)
  
  
  # Make a buffered version of clusters
  ij_clusters_buffered <- ij_clusters
  ij_clusters_buffered <- spTransform(ij_clusters_buffered, CRS( "+init=epsg:3347"))
  ij_clusters_buffered <- gBuffer(ij_clusters_buffered, byid = TRUE, width = -1000)
  # convert back to lat/lon
  ij_clusters_buffered <- spTransform(ij_clusters_buffered, CRS(proj4string(man3)))
  
  # Get cluster assignment
  x <- census
  coordinates(x) <- ~longitude+latitude
  proj4string(x) <- proj4string(ij_clusters)
  z <- over(x, polygons(ij_clusters))
  x$cluster <- z
  census_sp <- x
  census_right <- x@data
  
  # Join to census
  census <-
    left_join(x = census,
              y = census_right %>%
                dplyr::select(perm_id, cluster),
              by = 'perm_id')
  
  # Get whether in core or not
  z <- over(x, polygons(ij_intervention_zones_buffered))
  census_sp$in_buffer <- is.na(z)
  census_sp$in_core <- !is.na(z)
  # Join back to census
  census_right <- census_sp@data
  census <-
    left_join(x = census,
              y = census_right %>%
                dplyr::select(perm_id, in_buffer, in_core),
              by = 'perm_id')
  
  save(results,
       results_spv,
       results_sp,
       new_zone,
       census,
       census_sp,
       centroids,
       clusters,
       ij,
       ij_exact,
       ij_children,
       snap_shot,
       ij_clusters,
       ij_clusters_buffered,
       ij_intervention_zones,
       ij_intervention_zones_buffered,
       file = '../data/prepared_data.RData')
}
