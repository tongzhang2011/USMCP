

#Define function to calculate Moran's Index in parallel
moran.Parallel_tz <- function(protein.data, geneNames, nPixels, cores, multiPlate = FALSE, nsim = 999, nrow = 17, ncol = 29,
                              impute = FALSE, lag.max = (nrow*ncol)/2, plot = FALSE,
                              cutoff = 8) {
  registerDoParallel(cores = cores)
  fill.na <- function(x, i=5) {
    if( is.na(x)[i] ) {
      return( round(mean(x, na.rm=TRUE),0) )
    } else {
      #return( round(x[i],0) )
      return(x[i])
    }
  }
  
  a <- foreach(gene = iter(geneNames), .combine = "rbind", .packages = c("raster")) %dopar% {
    errorFlag <- FALSE
    filteredGenes <- dplyr::filter(protein.data, Gene.names == gene)
    vec <- as.numeric(filteredGenes[1,][1:nPixels])
    if (all(is.na(vec))) {
      return(data.frame(gene = gene, moran = NA, p.value = NA, valid.values = 0))
    }
    
    
    r <- raster(matrix(vec, nrow = nrow, ncol = ncol, byrow = F))
    
    
    
    if (impute == TRUE) {
      r <- raster::focal(r, w = matrix(1,3,3), fun = fill.na, pad = TRUE, na.rm = FALSE)
    }
    
    if (plot == TRUE) {
      if (impute == FALSE) {
        plot(as.matrix(r), xlab = '', ylab = '', xaxt = 'n', na.col = 'white',
             col = viridis(24), main = paste(gene))  
      }
      if (impute == TRUE) {
        plot(as.matrix(r), xlab = '', ylab = '', xaxt = 'n', na.col = 'white',
             col = viridis(24), main = paste(gene, "(Imputed)"))
      }
    }
    
    if (sum(!is.na(raster::as.matrix(r))) < cutoff) {
      return(data.frame(gene = gene, moran = NA, p.value = NA, valid.values = sum(!is.na(raster::as.matrix(r)))))
    }
    
    w <- spdep::poly2nb(raster::rasterToPolygons(r, na.rm = F))
    ww <- spdep::nb2listw(w, style = 'B')
    
    tryCatch({
      mi.mc <- spdep::moran.mc(raster::values(r),ww,nsim=nsim,na.action = na.omit, zero.policy = TRUE)
    },
    error = function(e){
      errorFlag <<- TRUE
      message("* Caught an error on gene ", gene)
    })
    if (errorFlag) {
      return(data.frame(gene = gene, moran = NA, p.value = NA, valid.values = NA))
    }
    
    
    return(data.frame(gene = gene, moran = mi.mc$statistic, p.value = mi.mc$p.value, valid.values = sum(!is.na(raster::as.matrix(r)))))
  }
  
  stopImplicitCluster()
  
  rownames(a) <- NULL
  
  return(a)
  
}
