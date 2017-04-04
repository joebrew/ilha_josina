#' Optimization of cluster formation
#' 
#' Use the brute force method to assign households to clusters in a way that minimizes inter-cluster distance.
#' @param cluster_size The size of each cluster
#' @param plot_map Whether to plot a map after each cluster formation (slows down operations significantly)
#' @param locations a SpatialPointsDataFrame containing 3 columns: id, lat, and lon
#' @param shp The shapefile to plot (only if plot_map is TRUE)
#' @param sleep How many seconds to sleep between plots
#' @param buffer A buffer size (in meters)
#' @param walk_along Whether to do a walk-based clustering 
#' as opposed to centroid based
#' @return A dataframe of size \code{times} times the number of points in \code{locations}, with columns indicating each \code{simulation_number}
#' @export


# Define our algorithm
cluster_optimize <- function(cluster_size = 10,
                             plot_map = FALSE,
                             sleep = 0,
                             locations, 
                             shp, 
                             start=c("far", "close", "random"),
                             rest=c("far", "close", "random"),
                             messaging = FALSE,
                             buffer = 0,
                             walk_along = FALSE){
  
  # 
  # cluster_size = 50
  # plot_map = FALSE
  # sleep = 0
  # locations = ij_children 
  # shp = ij
  # start='close'
  # rest='close'
  # messaging = TRUE
  # buffer = 1000
  # walk_along = TRUE
  
  require(dplyr)
  require(rgeos)
  
  # Get a distance matrix between all points
  distance_matrix <-
    spDists(x = locations,
            longlat = TRUE)
  
  # # Repeat [times] times the search
  # for (time in 1:times){
  
  ## We are adding a colum with a numeric index for each of the points id
  locations$index<-1:nrow(locations)
  # Create a fresh copy of locations
  locations_fresh <- locations
  # Specify that none of the points have yet been selected
  locations_fresh$selected <- FALSE
  # Create a column for clusters
  locations_fresh$cluster <- NA
  # Place holder for whether in buffer
  locations_fresh$in_buffer <- FALSE
  # # Create column for simulations
  # locations_fresh$simulation_number <- time
  # Create column for indication of whether full sized cluster or not
  locations_fresh$complete_cluster <- TRUE
  
  # Pick a start point
  # ONLY IF THERE ARE MORE THINGS TO BE SELECTED 
  # (we need to add conditionality here)
  # (the point which is furthest from all other points)
  possibles <- distance_matrix[!locations_fresh$selected, !locations_fresh$selected]
  # possibles <- spDists(x = locations_fresh[!locations_fresh$selected,],
  #                      longlat = TRUE)
  
  if(start=="far") {
    start_index <- locations_fresh$index[!locations_fresh$selected][which.max(rowSums(possibles))][1]  
  } else if(start=="close"){
    start_index <- locations_fresh$index[!locations_fresh$selected][which.min(rowSums(possibles))][1]  
  }else if(start=="random") {
    start_index <- sample(which(!locations_fresh$selected),1)
  }else {
    stop("start must be one of (far, close, random)")
  }
  
  # Start the cluster counter 
  cluster <- 1
  
  # Go until all the points are filled
  while(length(which(!locations_fresh$selected)) > 0){
    if(messaging){
      message(paste0('making cluster number ', cluster,'\n',
                     length(which(!locations_fresh$selected)),
                     ' points remaining'))   
    }
    
    # Use the start index to get a start point
    start_point <- locations_fresh[locations_fresh$index == start_index,]
    # Remove that start point from the list of eligibles
    locations_fresh$selected[locations_fresh$index == start_index] <- TRUE
    # Assign the cluster to the start point
    locations_fresh$cluster[locations_fresh$index == start_index] <- cluster

    all_distances <- distance_matrix[start_index,]
    all_distances <- data.frame(index = 1:nrow(locations),
                                distance = all_distances)
    # Remove those rows which are ineligible (already selected/start_point)
    all_distances <- 
      all_distances[! all_distances$index %in% which(locations_fresh$selected),]
    
    # Order by distance
    all_distances <- all_distances[order(all_distances$distance),]
    
    # Define whether this is an incomplete cluster
    incomplete_cluster <- (nrow(all_distances) + 1) < cluster_size
    
    if(incomplete_cluster){
      nearest <- all_distances
      these_indices <- locations_fresh$index[!locations_fresh$selected]
      # Mark if it's a full size cluster or not
      locations_fresh$complete_cluster[these_indices] <- !incomplete_cluster
      locations_fresh$cluster[these_indices] <- cluster
      locations_fresh$selected[these_indices] <- TRUE
      this_cluster_indices <- these_indices
      # If a sequence, just get the next nearest point
    } else if(walk_along){
      # WALK BASED
      so_far_this_cluster <- 1
      walk_starts_here <- start_point$index
      this_cluster_indices <- walk_starts_here
      distance_vector <- all_distances
      # Flag the first point as done in the data
      locations_fresh$cluster[walk_starts_here] <- cluster
      locations_fresh$selected[walk_starts_here] <- TRUE
      while(so_far_this_cluster < cluster_size){
        # Up the counter
        so_far_this_cluster <- so_far_this_cluster + 1
        # Pick the next nearest point
        new_point <- distance_vector$index[1]
        # points(locations_fresh[new_point,], col = 'green', pch = '.')
        # Sys.sleep(0.2)
        # Store the point
        this_cluster_indices <- c(this_cluster_indices,
                                  new_point)
        # Update the data
        locations_fresh$cluster[locations_fresh$index == new_point] <-
          cluster
        # And mark them as selected
        locations_fresh$selected[locations_fresh$index == new_point] <-
          TRUE
        # Create a new vector of distances
        new_distances <- distance_matrix[new_point,]
        new_distances <- data.frame(index = 1:nrow(locations),
                                    distance = new_distances)
        # Remove those rows which are ineligible (already selected/start_point)
        new_distances <- 
          new_distances[! new_distances$index %in% which(locations_fresh$selected),]
        # Reorder
        new_distances <- new_distances %>% arrange(distance)
        # Overwrite the distance vector
        distance_vector <- new_distances
      }
      
    } else {
      # CENTROID BASED
      # Get the cluster_size nearest points
      nearest <- all_distances[1:(cluster_size - 1),]
      # Mark those nearest points as part of the same cluster
      locations_fresh$cluster[nearest$index] <- cluster
      # And mark them as selected
      locations_fresh$selected[nearest$index] <- TRUE
      this_cluster_indices <- nearest$index
    }
    
    # If buffering, define buffer and exclude from future selection
    # those that are in it
    if(buffer > 0){
      lf <- locations_fresh[locations_fresh$selected &
                              !is.na(locations_fresh$cluster) &
                              locations_fresh$cluster == cluster,]
      these_indices <- lf$index
      lf <- rgeos::gConvexHull(lf)
      # Transform to projected form
      lfp <- spTransform(lf, CRS( "+init=epsg:3347" ) )
      # Get buffer at 1km
      lfb <- rgeos::gBuffer(lfp, width = buffer)
      # Get a secondary buffer at double that
      lfb2 <- rgeos::gBuffer(lfp, width = buffer * 2)
      # Un -project lfb
      lfb <- spTransform(lfb, proj4string(locations))
      lfb2 <- spTransform(lfb2, proj4string(locations))
      # Identify which points fall into the buffer
      in_buffer <- over(locations_fresh, polygons(lfb))
      in_buffer2 <- over(locations_fresh, polygons(lfb2))
      in_buffer <- ifelse(is.na(in_buffer), FALSE, TRUE)
      in_buffer2 <- ifelse(is.na(in_buffer2), FALSE, TRUE)
      in_buffer[these_indices] <- FALSE
      in_buffer2[these_indices] <- FALSE

      
      # Mark those in buffer as already if within 2 km
      locations_fresh$selected[in_buffer2] <- TRUE
      
      # Say which cluster it was
      # locations_fresh$cluster[in_buffer &
      #                           is.na(locations_fresh$cluster)] <- cluster
      # locations_fresh$in_buffer <- 
      #   ifelse(in_buffer &) [in_buffer] <- TRUE
    }
    
    # Get the start_point for the next round 
    # possibles <- spDists(x = locations_fresh[!locations_fresh$selected,],
    #                      longlat = TRUE)
    possibles <- distance_matrix[!locations_fresh$selected,
                                 !locations_fresh$selected]
    # Ensure that it remains a matrix regardless of size
    possibles <- as.matrix(possibles)
    if(nrow(possibles) > 0){
      if (rest=="far") {
        start_index <- locations_fresh$index[!locations_fresh$selected][which.max(rowSums(possibles))][1]  
      } else if(rest=="close") {
        start_index <- locations_fresh$index[!locations_fresh$selected][which.min(rowSums(possibles))][1]  
      } else if(rest=="random") {
        remaining <- locations_fresh$index[!locations_fresh$selected]
        if(length(remaining) == 0){
          start_index <- NA
        } else {
          start_index <- sample(remaining, 1)       
        }
      } else {
        stop("rest must be one of (far, close, random)")
      }
      
      # Move the cluster counter up
      cluster <- cluster + 1
      
      # Plot if necessary
      if(plot_map){
        colors <- ifelse(locations_fresh$selected, 'red', 'grey')
        plot(shp)
        maps::map.scale(ratio=FALSE, relwidth=0.2)
        points(locations_fresh, col = colors, pch = 1)
        points(locations_fresh[this_cluster_indices,], col = 'blue', pch = 1)
        plot(lf, add = TRUE)
        plot(lfb,
             add = TRUE)
        legend('topleft',
               legend = c('This cluster',
                          'Already selected',
                          'Not selected yet'),
               pch = c(1, 3, 3),
               col = c('blue', 'red', 'grey'),
               border = FALSE,
               bty = 'n')
        title(main = paste0('Simulation number ', 
                            # time, 
                            '\n',
                            'Cluster number ', cluster))
        Sys.sleep(sleep)
      }
    }
    
  } # all locations have now been selected
  # return the dataframe
  return(data.frame(locations_fresh))
}
