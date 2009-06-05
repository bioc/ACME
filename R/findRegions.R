"findRegions" <-
  function(x,thresh=0.0001) {
    regions <- list()
    for(i in 1:ncol(x)) {
      mat <- list()
      for(chrom in unique(chromosome(x))) {
        chromSub <- which(chromosome(x)==chrom)
        y <- rle(as.vector(vals(x)[chromSub,i]<thresh))
        startind <- min(chromSub)+cumsum(y$lengths)-y$lengths
        endind <- min(chromSub)+cumsum(y$lengths)-1
        mat[[chrom]] <- data.frame(Length=y$lengths,TF=y$values,
                          StartInd=startind,
                          EndInd=endind,
                          Sample=sampleNames(x)[i],
                          Chromosome=chromosome(x)[startind],
                          Start=start(x)[startind],
                          End=start(x)[endind]
                          )
      }
      regions[[sampleNames(x)[i]]] <- do.call(rbind,mat)
    }
    regions.df <- do.call('rbind',regions)
    meds <- apply(regions.df,1,function(tmp) {
      return(median(as.vector(vals(x)[(as.numeric(tmp[3])):(as.numeric(tmp[4]))
                                   ,tmp[5]])))
    }
                  )
    means <- apply(regions.df,1,function(tmp) {
      return(mean(as.vector(vals(x)[(as.numeric(tmp[3])):(as.numeric(tmp[4]))
                                 ,tmp[5]])))
    }
                   )
    regions.df <- data.frame(regions.df,Median=meds,Mean=means)
    return(regions.df)
  }

                        
