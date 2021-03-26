## ------------------------------------------------------------------------
library(data.table)
library(optparse)

option_list = list(
    make_option(c("-q", "--query"), type="character",
                help="Path to file containing paths to query files. Paths should be separated by a newline char."),
    make_option(c("-r", "--reference"), type="character",
                help="Path to file containing paths to reference files. Paths should be separated by a newline char."),
    make_option(c("-l", "--samplesizesqry"), type="numeric", default=482,
                help="Sample size of the query experiment.\n\t\t[default=%default]"),
    make_option(c("-s", "--samplesizeref"), type="character",
                help="Path to file containing sample sizes of reference files.\n\t\tMust contain the following collumns 'filename',  'N', 'Case', 'Control'. Case & Control or N may be NA not both."),
    make_option(c("-m", "--maf"), type="character",
                help="Path to MAF file. Must contain collumns 'SNP' & 'MAF'"),
    make_option(c("-o", "--out"), type="character", default="ColocResults.RData",
                help="Output file name.\n\t\t[default=./%default]"),
    make_option(c("-w", "--window"), type="numeric", default=250000,
                help="The window size used to extend the loci arround significant SNPs.\n\t\t[default=%default]"),
    make_option(c("-t", "--threshold"), type="numeric", default=5e-8,
                help="The threshold for significant SNPs.\n\t\t[default=%default]"),
    make_option(c("-n", "--colnamesquery"), type="character", default="rs,ps,chr,p_wald",
                help="Collumn names of: <snp id>,<posisiton>,<chromosome>,<pvalue> in query\n\t\t[default=%default]"),
    make_option(c("-p", "--colnamesreference"), type="character", default="rs,P",
                help="Collumn names of: <snp id>,<pvalue> in reference\n\t\t[default=%default]"),
    make_option(c("-f", "--fread"), action="store_true", default=FALSE,
                help="Use fread instead of read table to read the reference files.\n\t\t[default=%default]"),
    make_option(c("-z", "--zipped"), action="store_true", default=FALSE,
                help="Use are the refenence files zipped.\n\t\t[default=%default]"),
    make_option(c("-x", "--minsnps"), type="numeric", default=5,
                help="The minimum amount of SNPs in the query file.\n\t\t[default=%default]")
)

opt_parser    <- OptionParser(option_list=option_list, description="\nScript for running a approximate bayes coloc analysis using multiple input QTLs and compare them to multiple traits. Auto selects loci to compare to traits based on threshold and window size. Returns an .RData file containing the summary's of the coloc.abf function. Selection of quantative or case control is based upon the file given in the option -s. If case & control are NA a quantative analysis is used and vice versa.")
opt           <- parse_args(opt_parser)

tryCatch({
    maf.path      <- opt$maf
    maf           <- fread(maf.path, select=c("SNP", "MAF"), data.table=F, verbose=F, showProgress=F)
    rownames(maf) <- maf$SNP
    maf           <- maf[, 2,drop=F]

    ref.path      <- opt$samplesizeref
    ref.info      <- read.table(ref.path, header=T, stringsAsFactors=F, sep=";", row.names=1)
    rownames(ref.info) <- tolower(rownames(ref.info))
    
    qry.pathfile  <- opt$query
    qry.files     <- read.table(qry.pathfile, stringsAsFactors=F, header=F, sep=";")[,1]

    ref.pathfile  <- opt$reference
    ref.files     <- read.table(ref.pathfile, stringsAsFactors=F, header=F, sep=";")[,1]

    window        <- opt$window
    threshold     <- opt$threshold

    out.file      <- opt$out

    col.q         <- strsplit(opt$colnamesquery, ",")[[1]]
    col.r         <- strsplit(opt$colnamesreference, ",")[[1]]
    maf.col       <- "MAF"
    # These are fixed for the 500FG cytokine QTL mapping for now
    type.q        <- "quant"
    n.q           <- opt$samplesizesqry
    
}, warning=function(w){
    cat(paste0("\n[WARN]\t", w, "\n"))
    #print_help(opt_parser)
    #q(save="no")
},error=function(e){
    cat(paste0("\n[ERROR]\t", e, "\n"))
    print_help(opt_parser)
    q(save="no")
})

library(coloc)

## ------------------------------------------------------------------------
# Check if a position is within a region
in.region <- function(y, x, window) {
    x1 <- x - window
    x2 <- x + window
    if (x1 < 0){x1<-0}
    return(y >= x1 && y <= x2)
}

# Select the top snp for a given region (region is defined by the window size)
filter.regions <- function(data, window=250000, chr.col="chr", ps.col="ps", p.col="p_wald", return.type="names") {
    
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
get.loci <- function(top.snp, query, window=250000, chr.col="chr", ps.col="ps", p.col="p_wald"){
    # Select all the SNPs within a certain window of the top SNPs
    loci <- lapply(top.snp, function(snp, snps, window, chr.col, ps.col, p.col){
        chr <- snps[snp, chr.col]
        pos <- snps[snp, ps.col]
        loci <- rownames(snps[snps[,chr.col] == chr,][snps[snps[,chr.col] == chr ,ps.col] > pos - window &
                                                      snps[snps[,chr.col] == chr ,ps.col] < pos + window,])
    }, snps=query, window=window, chr.col=chr.col, ps.col=ps.col, p.col=p.col)
    return(loci)
}

# Run colocolization analysis using coloc package
# Expects 2 dataframes where the first collumn is the pvalue and the rownames are SNP id.
# In case of quantative data use type.q or type.r 'quant' and specify n.r or n.q = n
# In case of case control style use type.q or type.r 'cc' and specify n.r or n.q = c(n.case, n.control)
coloc.wrapper <- function(query, reference, maf, n.q, n.r, type.q="quant", type.r="quant") {
    tryCatch({
        if (type.q == "quant") {
             # Prepare and harmonize the data
            ds1 <- list(pvalues=query, N=n.q, type=type.q)
        } else if (type.q == "cc") {
            ratio.q <- n.q[1] / sum(n.q)
            ds1 <- list(pvalues=query, N=sum(n.q), s=ratio.q, type=type.q)
        }
        
        if (type.r =="quant") {
            ds2 <- list(pvalues=reference, N=n.r, type=type.r)
        } else if (type.r == "cc") {
            ratio.r <- n.r[1] / sum(n.r)
            ds2 <- list(pvalues=reference, N=sum(n.r), s=ratio.r , type=type.r)
        }
        
        sink("/dev/null")
        res <- coloc.abf(ds1, ds2, MAF=maf)
        sink()
        return(res)
    }, warning=function(w){
        cat(paste0("\n[WARN]\t", w, "\n"))
        return(NA)
    },error=function(e){
        cat(paste0("\n[ERROR]\t", e, "\n"))
        return(NA)
    })
}


## ------------------------------------------------------------------------
results   <- list()
cat("[INFO]\tCalculating coloc for: ", length(qry.files), " query files\n")
cat("[INFO]\tCalculating coloc for: ", length(ref.files), " reference files\n")

# For code readabillity this is done using for loops since it is not time critical
for (qry.file in qry.files) {
    
    cat("[INFO]\tCalculating coloc for query: ", basename(qry.file), "\n")
    # Load & clean the query data
    query           <- fread(qry.file, data.table=F, verbose=F, showProgress=F)
    rownames(query) <- query[,col.q[1]]
    query           <- query[,-which(colnames(query) == col.q[1])]
    query$chr       <- as.numeric(gsub("chr", "", query[,col.q[3]]))
    sig.query       <- query[query[,col.q[4]] < threshold,]
    query.name      <- tail(strsplit(qry.file, "/")[[1]], n=1)

    # If fewer then N SNPs in query skip 
    if (nrow(query) <= opt$minsnps) {
        cat("[WARN] Not enough SNPs in query file, skipping. To overide change --minsps to 0\n") 
        next
    }
    
    if(nrow(sig.query) > 0) {
        # Select which SNPs are the top in their respective region
        if (nrow(sig.query) == 1) {
            top.snp <- rownames(sig.query)
        } else {
            top.snp <- filter.regions(sig.query, window=window, p.col=col.q[4], chr.col=col.q[3], ps.col=col.q[2])
        }

        loci        <- get.loci(top.snp, query=query, window=window, p.col=col.q[4], chr.col=col.q[3], ps.col=col.q[2])
        names(loci) <- top.snp
        # TODO: Perhaps do this using a mclapply for parrallelization
        for (ref.file in ref.files) {
            # Use either fread or read table to read the files. Read table can read .gz
            if (opt$fread) {
                # Read the file as zipped or not
                cat("[INFO]\tReading reference file: ", basename(ref.file), "\n")
                if (opt$zipped){
                    ref       <- fread(paste0("zcat < ", ref.file), data.table=F, header=T, stringsAsFactors=F, verbose=F, showProgress=F)
                } else {
                    ref       <- fread(ref.file, data.table=F, header=T, stringsAsFactors=F, verbose=F, showProgress=F)
                }
            } else {
                ref       <- read.table(ref.file, header=T, stringsAsFactors=F)
                             #colClasses=c(character(), character(), numeric(), numeric()))
            }
            ref           <- ref[!duplicated(ref[,col.r[1]]),]
            rownames(ref) <- ref[,col.r[1]]
            ref.name      <- tail(strsplit(ref.file, "/")[[1]], n=1)
            n.r           <- ref.info[tolower(ref.name), 1]
            type.r        <- "quant"
            if (is.na(n.r)) {
                n.r       <- as.numeric(ref.info[tolower(ref.name), 2:3])
                type.r    <- "cc"
                # If not quant or cc skip the file 
                if (sum(is.na(n.r)) > 0) {
                    cat("[WARN]\tNo sample size found for file:", ref.name,
                        "\n\tSkipping to next file\n")
                    next
                }
            }
            
            for (locus in top.snp) {
                name     <- locus
                lc <- intersect(intersect(loci[[locus]], rownames(ref)), rownames(maf))
                if (length(lc) > 0) {
                    results[[query.name]][[name]][[ref.name]] <- coloc.wrapper(query[lc, col.q[4]],
                                                                           ref[lc, col.r[2]],
                                                                           maf[lc, maf.col],
                                                                           n.q=n.q,
                                                                           n.r=n.r,
                                                                           type.q=type.q,
                                                                           type.r=type.r)[["summary"]]
               } else {
                  cat("[WARN] No overlapping snps for query: ", query.name, " and reference: ", ref.name, "\n") 
                   t        <- as.numeric(rep(NA, 6))
                   names(t) <- c("nsnps", "PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")
                   results[[query.name]][[name]][[ref.name]] <- t
               }
            }
        }
        
    } else if (nrow(sig.query) == 0) {
        cat("[WARN]\tNo significant hits in file: ", basename(qry.file), "\n\tfor threshold: ", threshold, "\n\tSkipping to next file\n")
        next
    } else {
        cat("[ERROR]\tSomething went wrong. Skipping file\n")
        next
    }
}
# Save the results
save(results, file=out.file)

