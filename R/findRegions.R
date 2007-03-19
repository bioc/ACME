"findRegions" <-
  function(calcobj,thresh=0.0001) {
    vals <- calcobj@vals
    cols <- colnames(calcobj@vals)
    annot <- calcobj@annotation
    regions <- list()
    for(i in cols) {
      for(chrom in unique(as.character(annot$Chromosome))) {
        chromSub <- which(annot$Chromosome==chrom)
        y <- rle(as.vector(vals[chromSub,i]<thresh))
        startind <- min(chromSub)+cumsum(y$lengths)-y$lengths
        endind <- min(chromSub)+cumsum(y$lengths)-1
        mat <- data.frame(Length=y$lengths,TF=y$values,
                          StartInd=startind,
                          EndInd=endind,
                          Sample=i,
                          Chromosome=calcobj@annotation[startind,'Chromosome'],
                          Start=calcobj@annotation[startind,'Location'],
                          End=calcobj@annotation[endind,'Location']
                          )
        if(length(regions[[i]])>0) {
          regions[[i]] <- rbind(regions[[i]],mat)
        } else {
          regions[[i]] <- mat
        }
      }
    }
    regions.df <- do.call('rbind',regions)
    meds <- apply(regions.df,1,function(x) {
      return(median(as.vector(vals[(as.numeric(x[3])):(as.numeric(x[4]))
                                   ,x[5]])))
    }
                  )
    means <- apply(regions.df,1,function(x) {
      return(mean(as.vector(vals[(as.numeric(x[3])):(as.numeric(x[4]))
                                 ,x[5]])))
    }
                   )
    regions.df <- data.frame(regions.df,Median=meds,Mean=means)
    return(regions.df)
  }

                        
