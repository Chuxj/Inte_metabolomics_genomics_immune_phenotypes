library(getopt)
parameter=matrix(c("snp","s",1,"character","exp","e",1,"character","cov","c",1,"character","out","o",1,"character"),byrow=T,ncol=4)
opt=getopt(parameter)

library(MatrixEQTL);

useModel = modelLINEAR;
SNP_file_name=paste(opt$snp);
expression_file_name=paste(opt$exp);
covariates_file_name=paste(opt$cov);
output_file_name=paste(opt$out);

pvOutputThreshold = 1;
errorCovariance = numeric();

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;
snps$fileSkipColumns = 1;
snps$fileSliceSize = 2000;
snps$LoadFile(SNP_file_name);

gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;
gene$fileSkipColumns = 1;
gene$fileSliceSize = 2000;
gene$LoadFile(expression_file_name);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

