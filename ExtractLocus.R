# Check if a position is within a region
in.region <- function(y, x, window) {
  x1 <- x - window
  x2 <- x + window
  if (x1 < 0){x1<-0}
  return(y >= x1 && y <= x2)
}

# Select the top snp for a given region (region is defined by the window size)
filter.regions <- function(data, window=250000, chr.col="chr", ps.col="pos", p.col="p.value", return.type="names") {
  
  # Split so each chromosome is evaluated individually
  loci <- lapply(unique(data[,chr.col]), function(chr, data, window, return.type){
    snps <- data[data[,chr.col] == chr,]
    snps$loci <- rep(NA, nrow(snps))
    # Loop over and assign the locus id
    for (i in 1:nrow(snps)){
      cur.locus <- sapply(snps[, ps.col], in.region, x=snps[i, ps.col], window=window)
      snps[cur.locus, "loci"] <- i
      if(sum(is.na(snps$loci)) == 0){
        break
      }
    }
    # After all id's have been assigned only select the most
    # significant snp for that locus
    unique.snps <- sapply(unique(snps$loci), function(locus, snps){
      c(rownames(snps[snps$loci == locus,])[snps[snps$loci == locus, p.col] == min(snps[snps$loci == locus, p.col])])[1]
    }, snps=snps)
    # Remove now redundant loci collumn
    snps$loci <- NULL
    
    if (return.type == "names") {
      return(unique.snps)
    } else if (return.type == "df") {
      return(snps[unique.snps,])
    }
  }, data=data, window=window, return.type=return.type)
  
  # Merge the list of dataframes
  if (return.type == "names") {
    return(unlist(loci))
  } else if (return.type == "df") {
    return(do.call("rbind", loci))
  }
}

# Get the SNPs surrounding a significant loci for a given dataset
get.loci <- function(top.snp, query, window=250000, chr.col="chr", ps.col="pos", p.col="p.value"){
  # Select all the SNPs within a certain window of the top SNPs
  loci <- lapply(top.snp, function(snp, snps, window, chr.col, ps.col, p.col){
    chr <- snps[snp, chr.col]
    pos <- snps[snp, ps.col]
    loci <- rownames(snps[snps[,chr.col] == chr,][snps[snps[,chr.col] == chr ,ps.col] > pos - window &
                                                    snps[snps[,chr.col] == chr ,ps.col] < pos + window,])
  }, snps=query, window=window, chr.col=chr.col, ps.col=ps.col, p.col=p.col)
  return(loci)
}