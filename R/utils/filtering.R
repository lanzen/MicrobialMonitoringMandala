# Removes whole taxa if maximum relative abundance is below cutoff
# Good for uneven sampling and uneven distribution, to avoid losing relatively
# rare but discrimnant taxa.
# Anders Lanzen 2018

# Recommended minRelativeAbundance is 10/(the number of reads in the least sequeneced sample)

dropRareByMaxAbundance = function(relativeAbundances, minRelativeAbundance, 
                                  taxaAsRows = FALSE) {
  # Transpose if taxa are rows
  if (taxaAsRows) relativeAbundances = as.data.frame(t(relativeAbundances))
  
  # Filter and return
  maxTaxonAbundances = sapply(relativeAbundances, max)
  filtered = relativeAbundances[,maxTaxonAbundances >= minRelativeAbundance]
  if (taxaAsRows) filtered = (as.data.frame(t(filtered)))
  return(filtered)
}

# Removes whole taxa if average relative abundance is below cutoff
# Good for uneven sampling and similar samples (even distribution)
dropRareByAvgAbundance = function(relativeAbundances, minRelativeAbundance, 
                                  taxaAsRows = FALSE) {
  require(vegan)
  # Transpose if taxa are rows
  if (taxaAsRows) relativeAbundances = as.data.frame(t(relativeAbundances))
  
  # Filter and return
  # avgTaxonAbundances = colMeans(relativeAbundances)
  filtered = relativeAbundances[,avgTaxonAbundances >= minRelativeAbundance]
  if (taxaAsRows) filtered = (as.data.frame(t(filtered)))
  return(filtered)
}

# Only filter individual occurences (bad with uneven sampling -> FNs in samples with many reads)
# Accept and returns either rel. or abs. abundance df (2017-01-24)
dropRareTaxa <- function(df, minRelativeAbundance){
  dfM = as.matrix(df)
  ra = decostand(df,method="total")
  raM = as.matrix(ra)
  
  # Filter below min relative abundances
  fdist=as.data.frame(ifelse(raM>minRelativeAbundance,dfM,0))
  
  names(fdist)=names(df)
  row.names(fdist)=row.names(df)
  nonzero = fdist[,colSums(fdist)>0]
  return(nonzero)
}

# Only filter individual occurences (bad with uneven sampling -> FNs in samples with many reads
# Accept and returns either rel. or abs. abundance df (2017-01-24)
dropRareTaxa <- function(df, maxRelativeAbundance){
  dfM = as.matrix(df)
  require(vegan)
  ra = decostand(df,method="total")
  raM = as.matrix(ra)
  fdist=as.data.frame(ifelse(raM<maxRelativeAbundance,dfM,0))
  names(fdist)=names(df)
  row.names(fdist)=row.names(df)
  nonzero = fdist[,colSums(fdist)>0]
  return(nonzero)
}

dropSingletons <- function(df){
  require(vegan)
  good = df[,colSums(df)>1]
  return(good)
}


# Compares abs. read abundance to total OTU reads/k as lower limit -> 
# FNs in undersampled. Similar to UNCROSS
filterCrossContaminants = function(otus, minAbundanceQuote){
  totals = colSums(otus)
  otus.clean = matrix(nrow=dim(otus)[1],ncol=dim(otus)[2])
  for (i in c(1:dim(otus)[2])) {
    for (j in c(1:dim(otus)[1])) {
      if(otus[j,i] > totals[i] / minAbundanceQuote) {
        otus.clean[j,i] = otus[j,i]
      }
      else {
        otus.clean[j,i] = 0
      }
    }
  }
  otus.clean = as.data.frame(otus.clean)
  row.names(otus.clean) = row.names(otus)
  names(otus.clean) = names(otus)
  return(otus.clean)
}

# Based on relative abundance. Better for uneven sample sizes
filterCrossContaminants2 = function(otus, minRAQuote = 200, taxaAsRows = FALSE){
  if (taxaAsRows) otus = as.data.frame(t(otus))
  otusM = as.matrix(otus)
  require(vegan)
  ra = decostand(otus,method="total")
  raM = as.matrix(ra)
  meanTaxonAbundance = as.vector(colMeans(ra))
  otus.clean = matrix(nrow=dim(otus)[1],ncol=dim(otus)[2])
  #Iterate over each column (OTU)
  for (i in c(1:dim(otus)[2])) {
    otus.clean[,i] = ifelse(raM[,i] >= meanTaxonAbundance[i] / minRAQuote, otusM[,i], 0)
  }
  otus.clean = as.data.frame(otus.clean)
  row.names(otus.clean) = row.names(otus)
  names(otus.clean) = names(otus)
  if (taxaAsRows) otus.clean = as.data.frame(t(otus.clean))
  return(otus.clean)
}