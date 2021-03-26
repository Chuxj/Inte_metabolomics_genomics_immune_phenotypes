##########MATRIX EQTL POST PROCESSING##########
#
# Name:       Twan Spenkelink
# Date:       October 19, 2017
# Contact:    t.o.spenkelink@umcg.nl
# Purpose:    Post processing for matrixEQTL output.

### LOAD LIBRARIES ###
# Load optparse library to support command line arguments
library(optparse)
# Load qqman to construct Manhattan plots
library(qqman)
# Load data.table to use fread to provide faster data reading 
library(data.table)
# Load RColorBrewer to make use of nice color pallets
library(RColorBrewer)
# Load ggplot2 for nice plotting
library(ggplot2)
# Source functions from ExtractLocus.R script used for finding top SNP's and loci's
source("/groups/umcg-wijmenga/tmp04/umcg-xchu/database/post_impu_500FG_300DM/cellPercent_mapping/post_analysis/ExtractLocus.R")     # <- for cluster
                                                                  
### DEFINE COMMAND LINE OPTIONS ###
option_list = list(
  make_option(c("-r", "--result"), action="store", default=NA, type='character',
              help="Path to the output file which contains the result of the cQTL mapping"),
  
  make_option(c("-p", "--position"), action="store", default=NA, type='character',
              help="Path to the file containing the snp positions"),
  
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="Path to the folder where the plots are going to be stored"),
  
  make_option(c("-g", "--genotype"), action="store", default=NA, type='character',
              help="Path to the file containing the snp information"),
  
  make_option(c("-c", "--covariate"), action="store", default=NA, type='character',
              help="Path to the file that contains the covariate data"),
  
  make_option(c("-t", "--threshold"), action="store", default=5e-7, type='numeric',
              help="Threshold of the pvalue that is being analyzed, default = 5e-7"),
  
  make_option(c("-w", "--regionwindow"), action="store", default=2500000, type='numeric',
              help="Window that has to be used for the region selection, default = 250000")
)

opt_parser    <- OptionParser(option_list=option_list, description="\nScript that provides post processing features for MatrixEQTL output")
opt           <- parse_args(opt_parser)

### READ DATA ###
cat("[INFO]\tReading data\n")

# Read SNP data
genotype.path     <- opt$genotype
genotype.data     <- fread(genotype.path, header = TRUE, data.table = FALSE)

# Read covariate data
covariate.path    <- opt$covariate
covariate.data    <- read.table(covariate.path, header = TRUE, fill=TRUE,sep = "\t")

# Read file that contains cQTL analysis results
result.path       <- opt$result
cqtl.result       <- fread(result.path, header = TRUE, data.table = FALSE, select = c("SNP", "gene", "p-value"))

# Read snp position data
position.path     <- opt$position
snp.position      <- fread(position.path, header = TRUE, data.table = FALSE)

# Define output directory
output.path	      <- opt$output

# Additional information
snp.info.dir      <- opt$snpinfo
region.window     <- opt$regionwindow
threshold         <- opt$threshold
covariate.name    <- tail(rownames(covariate.data), n=1)

# Check if output directory already existst, if not: create the directory
# TODO: add covariate name in manhattan dir name

# Check if output directory already existst, if not: create the directory
if (!dir.exists(paste0(output.path, "/locus_files_",threshold, "_" ,covariate.name,"/"))) {
  dir.create(paste0(output.path, "/locus_files_",threshold, "_" ,covariate.name, "/"), recursive=T)
}

# Check if output directory already existst, if not: create the directory
if (!dir.exists(paste0(output.path, "/snps_per_cytokine/"))) {
  dir.create(paste0(output.path, "/snps_per_cytokine/"), recursive=T)
}

### DEFINE FUNCTIONS ###

# Function used for showing the sample sizes in the boxplots
give.n <- function(x){
  return(c(y = median(x), label = length(x)))
}

### WORK ###
cat("[INFO]\tReading data: done \n")
cat("[INFO]\tProcessing data... \n")
colnames(snp.position)[1] <- "SNP"

snp.pvalue.pos            <- merge(cqtl.result, snp.position, by = "SNP")
colnames(snp.pvalue.pos)  <- c("SNP", "GENE", "P", "CHR", "POS")
col.q                     <- colnames(snp.pvalue.pos)
sig.pvalues               <- subset(snp.pvalue.pos, P < threshold)
stimulations              <- levels(factor(snp.pvalue.pos$GENE))
occurences.count          <- 0
not.sig.count             <- 0
top.snp.info              <- data.frame()

for (stimulus in stimulations){
  
  # Create subset of only significant hits below the threshold
  sig.pvalues.filtered <- sig.pvalues[sig.pvalues[,col.q[2]] == stimulus,]
  pvalues.filtered <- snp.pvalue.pos[snp.pvalue.pos[,col.q[2]] == stimulus,]
  
  # Save each stimuls/cytokine group in a seperate file (used for coloc)
  write.table(pvalues.filtered, file=paste0(output.path, "/snps_per_cytokine/", stimulus, "_cQTL_output.txt"), quote=FALSE, row.names=FALSE)    
  
  # Change rownames to the SNP id's
  rownames(sig.pvalues.filtered)  <- sig.pvalues.filtered[,col.q[1]]
  rownames(pvalues.filtered)  <- pvalues.filtered[,col.q[1]]
  
  # Check if there are significant hits
  if(nrow(sig.pvalues.filtered) > 0) {
    cat("[INFO]\tSignificant hit(s) for threshold: ", threshold, "\n")
    # Check number of significant hits
    if (nrow(sig.pvalues.filtered) == 1) {
      # When there is only one significant SNP, this will automaticaly be the top SNP
      top.snp <- rownames(sig.pvalues.filtered)
    } else {
      
      # Find top SNP(s)
      top.snp <- filter.regions(data = sig.pvalues.filtered, chr.col = col.q[4], ps.col = col.q[5], p.col = col.q[3])
      cat("[INFO]\tFound ",length(top.snp),"significant top SNP(s) for cytokine ",stimulus,"\n")
    }
    
    # Find loci of the top snp(s)
    loci        <- get.loci(top.snp, pvalues.filtered, chr.col = col.q[4], ps.col = col.q[5], p.col = col.q[3])
    names(loci) <- top.snp
    # Save loci information to a file which can be used to create a locuszoom plot
    for (locus in names(loci)) {
      loci.region <- pvalues.filtered[loci[[locus]], col.q, drop=FALSE]
      write.table(loci.region, file=paste0(output.path, "/locus_files_",threshold, "_" ,covariate.name, "/", locus, "_", stimulus,"_", region.window, "_for_locuszoom.txt"), quote=FALSE, row.names=FALSE)    
    }
    
    for (snp.id in top.snp){
      
      # For every snp + cytokine/stimulus a dataframe is made containg all the necesary information
      # needed for the interaction plot
      top.snp.info 		                <- rbind(top.snp.info, subset(sig.pvalues.filtered, sig.pvalues.filtered$SNP == snp.id))
    }
  }
}
write.table(top.snp.info, file=paste0(output.path, "/locus_files_",threshold, "_" ,covariate.name, "/", "all_top_snps.txt"), quote=FALSE, row.names=FALSE)
cat("[INFO]\tNumber of cytokines that do not have significant (",threshold ,") SNPs ", not.sig.count)
cat("[INFO]\tNumber of top snps for threshold ", threshold, " that do not have enough samples per SNP: ",occurences.count)
