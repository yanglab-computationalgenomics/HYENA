library(optparse)
library(data.table)
library(RMThreshold)

option_list <- list( 
  make_option(c("-e", "--exprs"), action="store", default="exprs.txt",
              help="Expression data file [default %default] to be used as the input"),
  make_option(c("-w", "--write"), action="store", default= getwd(), 
              help="Output directory"),
  make_option(c("-o", "--output"), action="store", default="exprs_noisy.txt",
              help="Output file name [default %default]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="")
args <- parse_args(parser, positional_arguments=0)
opt <- args$options

# load expression data
dat.expression <- fread(opt$exprs, stringsAsFactors = FALSE)
dat.expression <- as.matrix(dat.expression)
print("expression data loaded:")
dat.expression[1:5,1:5]

# annotate rows
rownames(dat.expression) <- dat.expression[,1]
dat.expression <- dat.expression[,-1]
class(dat.expression) <- "numeric"

# remove genes that are 0 for all patients
dat.expression <- dat.expression[rowSums(dat.expression[]) > 0 ,]
dat.expression[1:5,1:5]

# add noise to the entire expression matrix
dat.expression.noisy <- add.Gaussian.noise(as.matrix(dat.expression),
                                           mean = 0.000000001,
                                           stddev = 0.000000001,
                                           symm = FALSE)
dat.expression.noisy <- as.data.frame(dat.expression.noisy)
dat.expression.noisy[1:5,1:5]

# write noisy expression data as output
outputName <- opt$output
write.table(dat.expression.noisy, outputName, sep = "\t", quote = FALSE, row.names = TRUE)
