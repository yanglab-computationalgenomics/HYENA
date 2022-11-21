oldw <- getOption("warn")
options(warn = -1)

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))

option_list <- list( 
  make_option(c("-l", "--filelist"), action="store",
              help="Comma separated list of result files. Make sure the files are listed in an INCREASING order for the number of parameters (e.g. PCs) included in regression model"),
  make_option(c("-p", "--power"), action="store", default= 0.8,
              help="Desired power [default %default]"),
  make_option(c("-f", "--fdrcutoff"), action="store", default= 0.1, 
              help="FDR cut-off for significant genes [default %default]"),
  make_option(c("-w", "--write"), action="store", default="./intermediate/",
              help="Folder for intermediate files [default %default]"),
  make_option(c("--emp"), action="store_true", default=FALSE,
              help="P values to process are empirical p-values (p.emp.onesided) [default %default]"),
  make_option(c("-d", "--dir"), action="store", default= "./results/", 
              help="Output folder [default %default]"),
  make_option(c("-x", "--prefix"), action="store", default="DAT",
              help="Prefix for output files [default %default]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="")
args <- parse_args(parser, positional_arguments=0)
opt <- args$options


filelist <- strsplit(opt$filelist, ",")
filelist <- filelist[[1]]


run <- c()
num <- c()

for (i in 1:length(filelist)) {
	res <- fread(filelist[i])
	
  if (opt$emp) {
    res <- res[res$p.emp.fdr <= opt$fdrcutoff,]
  } else {
    res <- res[res$fdr.onesided <= opt$fdrcutoff,]
  }
	
  sighits <- nrow(res)

	run <- c(run, filelist[i])
	num <- c(num, sighits)
}

tally <- cbind(run, num)
tally <- as.data.frame(tally)
tally$num <- as.numeric(as.character(tally$num))
print(tally)

if (file.exists(opt$write)) {
 write.table(tally, file = paste0(opt$dir, opt$prefix, ".pc_tally.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
} else {
 dir.create(opt$write)
 write.table(tally, file = paste0(opt$dir, opt$prefix, ".pc_tally.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
}


pcmax <- max(tally$num)
powercut <- pcmax*as.numeric(as.character(opt$power))

file.to.copy <- filelist[which(tally$num >= powercut)[1]]
res <- fread(file.to.copy)

if (opt$emp) {
  res <- res[res$p.emp.fdr <= opt$fdrcutoff,]
} else {
  res <- res[res$fdr.onesided <= opt$fdrcutoff,]
}

write.table(res, file = paste0(opt$dir, gsub(".txt", "", str_split(file.to.copy, "\\/")[[1]][length(str_split(file.to.copy, "\\/")[[1]])]), "_sig.txt"), quote = FALSE, sep = "\t", row.names = FALSE)



options(warn = oldw)
