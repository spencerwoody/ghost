
## df contains latitude and longitude; shape is sf object for this area
makeNewPoint <- function(df, shape) {

  require(sf)

  whilecount <- 0
  
  while (TRUE) {

    ## prop <- runif(1)

    prop <- rbeta(1, 0.5, 0.5)

    mypoints <- sample(nrow(df), 2)
    
    newpoint <- c(prop * df$LONGITUDE[mypoints[1]] +
                  (1 - prop) * df$LONGITUDE[mypoints[2]],
                  prop * df$LATITUDE[mypoints[1]] +
                  (1 - prop) * df$LATITUDE[mypoints[2]])

    newpointsf <- st_point(newpoint)

    stopcond <- suppressMessages(st_intersects(newpointsf, shape, sparse = FALSE))
    
    if (stopcond) {
      break
    } else {
      whilecount <- whilecount + 1
      if (whilecount > 1000) {
        stop("Never found a point...")
      }
    }
    
  }

  return(newpoint)

}

## Wrapper for previous function
newPoints <- function (n, df, shape, verbose = TRUE) {
  out <- data.frame(LONGITUDE = rep(NA, n),
                    LATITUDE = rep(NA, n))

  if (verbose) prog <- progress_estimated(n)

  for (i in 1:n) {
    out[i, ] <- makeNewPoint(df, shape)
    if (verbose) prog$tick()$print()
  }
  if (verbose) cat("\n")

  return(out)
}
