#!/bin/env
oldw <- getOption("warn")
options(warn = -1)

# load packages
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(optparse))

option_list <- list( 
  make_option(c("-r", "--ref"), action="store", 
              help="Gene annotation with TSS window"),
  make_option(c("-i", "--id"), action="store",
              help="Sample id file"),
  make_option(c("-b", "--bedpe"), action="store", 
              help="Folder for bedpe files"),
  make_option(c("--ingene"), action="store", default=1,
              help="Remove SVs located entirely in a gene [default %default]"),
  make_option(c("--del"), action="store", default=0,
              help="Remove small deletions [default %default]"),
  make_option(c("--dels"), action="store", default=10000,
              help="Small deletion size cutoff [default %default]"),
  make_option(c("--dup"), action="store", default=1,
              help="Remove small duplications [default %default]"),
  make_option(c("--dups"), action="store", default=10000,
              help="Small duplications size cutoff [default %default]"),
  make_option(c("--inv"), action="store", default=0,
              help="Remove small inversion [default %default]"),
  make_option(c("--invs"), action="store", default=10000,
              help="Small inversion size cutoff [default %default]"),
  make_option(c("-w", "--write"), action="store", default="./intermediate/",
              help="Folder for intermediate files [default %default]"),
  make_option(c("-x", "--prefix"), action="store", default="DAT",
              help="Prefix for output files [default %default]"),
  make_option(c("-v","--verbose"), action="store_true", default=FALSE,
              help="Verbose mode for troubleshooting [default %default]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="")
args <- parse_args(parser, positional_arguments=0)
opt <- args$options


########################################################################
################### Map SVs to gene promoter regions ###################
########################################################################

# load sample ids
samples <- fread(opt$id, stringsAsFactors = FALSE)
aliqid_wgs <- samples$aliquot_id_wgs
aliqid_rnaseq <- samples$aliquot_id_rnaseq

# load TSS windows for genes in Reference (GENCODE) and convert into a genomicRanges object
print("Loading TSS windows ...")
tss_windows <- fread(opt$ref, stringsAsFactors = FALSE)

# convert TSS window reference data into GRanges object
tss.ranges <- makeGRangesFromDataFrame(tss_windows,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=FALSE,
                                       # seqinfo=NULL,
                                       seqnames.field=c("seqnames", "seqname",
                                                        "chromosome", "chrom",
                                                        "chr", "chromosome_name",
                                                        "seqid", "CHROM"),
                                       start.field="tss_left", # original:"start"
                                       end.field=c("tss_right"), # original: "c("end, "stop")
                                       strand.field="STRAND",
                                       starts.in.df.are.0based=FALSE)


# map SV breakpoints to nearby genes (within 500Kb TSS window)
## read the "*.somatic.sv.bedpe" files with breakpoint info
print("Identifying SV breakpoint files...")
file_list <- list.files(path = opt$bedpe, pattern = ".bedpe$")


### select for a list of files that are on the sample aliquot information file
### this line allows for running the code only on a list of files that are provided in the "sample_ids_whitelist_tumor" file
file_list <- file_list[(vapply(str_split(file_list,"\\."), `[`, 1, FUN.VALUE=character(1))) %in% aliqid_wgs]

print("Concatenating all breakpoints...")
cat.sv <- data.frame(chr = as.character(),
                     pos = as.integer(),
                     sv_id = as.character(),
                     pe_support = as.integer(),
                     svstrand = as.character(),
                     svclass = as.character(),
                     svmethod = as.character(),
                     pair.chr = as.character(),
                     pair.pos = as.integer(),
                     aliq_id = as.character(),
                     stringsAsFactors = FALSE)

for (i in 1:length(file_list)){
  sv.bedpe <- read.delim(file = paste0(opt$bedpe,file_list[i]), stringsAsFactors = FALSE)
  if (nrow(sv.bedpe) == 0){next}
  
  # filter out the intrachromosomal SVs that are less than user defined size cutoff
  filt.out <- c()
  
  if (opt$del == 1)
  {
		for (j in 1:nrow(sv.bedpe))
		{
    	if (sv.bedpe[j,1] == sv.bedpe[j,4])
    	{
    		if (sv.bedpe[j,9] == "+" && sv.bedpe[j,10] == "-")
    		{
      		size <- sv.bedpe[j,2]-sv.bedpe[j,5]
      		if (abs(size) <= opt$dels)
      		{
        		filt.out <- c(filt.out, j)
     	 		}
    		}
    	}
  	}
  }
  if (opt$dup == 1)
  {
		for (j in 1:nrow(sv.bedpe))
		{
    	if (sv.bedpe[j,1] == sv.bedpe[j,4])
    	{
    		if (sv.bedpe[j,9] == "-" && sv.bedpe[j,10] == "+")
    		{
      		size <- sv.bedpe[j,2]-sv.bedpe[j,5]
      		if (abs(size) <= opt$dups)
      		{
        		filt.out <- c(filt.out, j)
     	 		}
    		}
    	}
  	}
  }
  if (opt$inv == 1)
  {
		for (j in 1:nrow(sv.bedpe))
		{
    	if (sv.bedpe[j,1] == sv.bedpe[j,4])
    	{
    		if (sv.bedpe[j,9] == sv.bedpe[j,10])
    		{
      		size <- sv.bedpe[j,2]-sv.bedpe[j,5]
      		if (abs(size) <= opt$invs)
      		{
        		filt.out <- c(filt.out, j)
     	 		}
    		}
    	}
  	}
  }
  
  if(opt$verbose){print(paste0(file_list[i], " # of small SV filtered: ", length(filt.out)))}
  if (!is.null(filt.out)){sv.bedpe <- sv.bedpe[-filt.out,]}
  
  # linearize breakpoints
  temp1 <- sv.bedpe[,c(1,2,7,8,9,11,12,4,5)]
  colnames(temp1) <- c("chr", "pos", "sv_id", "pe_support" ,"svstrand", "svclass", "svmethod", "pair.chr", "pair.pos")
  temp2 <- sv.bedpe[,c(4,5,7,8,10,11,12,1,2)]
  colnames(temp2) <- c("chr", "pos", "sv_id", "pe_support" ,"svstrand", "svclass", "svmethod", "pair.chr", "pair.pos")
  sv.bedpe <- rbind(temp1, temp2)
  
  aliq_id <- strsplit(file_list, "\\.")[[i]][1]
  sv.bedpe[(ncol(sv.bedpe)+1)] <- rep(aliq_id, nrow(sv.bedpe))
  cat.sv <- rbind(cat.sv, sv.bedpe)
}
colnames(cat.sv) <- c("chr", "pos", "sv_id", "pe_support" ,"svstrand", "svclass", "svmethod", "pair.chr", "pair.pos", "aliq_id")

#if (file.exists(opt$write))
#{
# write.table(cat.sv, file = paste0(opt$write, opt$prefix, ".bedpe_lin.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
#}#else
#{
# dir.create(opt$write)
# write.table(cat.sv, file = paste0(opt$write, opt$prefix, ".bedpe_lin.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
#}


## convert the breakpoint data into a genomicRanges object
print("Converting SV data frame into a GRanges object...")
sv.bedpe.ranges <- makeGRangesFromDataFrame(cat.sv,
                                            keep.extra.columns=TRUE,
                                            ignore.strand=TRUE,
                                            # seqinfo=NULL,
                                            seqnames.field=c("seqnames", "seqname",
                                                             "chromosome", "chrom",
                                                             "chr", "chromosome_name",
                                                             "seqid"),
                                            start.field="pos", # original:"start"
                                            end.field=c("pos"), # original: "c("end, "stop")
                                            #strand.field="strand",
                                            starts.in.df.are.0based=FALSE)

## find the overlap between two genomic range data
print("Mapping SVs to the gene-based TSS windows...")
merged1 <- mergeByOverlaps(tss.ranges, sv.bedpe.ranges)
if(opt$verbose){head(merged1)}
# merged1 <- S4Vectors::DataFrame(x = merged1)
# merged1 <- merged1[,c(1,2,3,5,6,9,15,19)]

chrom <- as.character(merged1$tss.ranges@seqnames)
gene.start <- as.integer(merged1$START)
gene.end <- as.integer(merged1$END)
strand = as.character(merged1$tss.ranges@'strand')
type = as.character(merged1$TYPE)
gene = as.character(merged1$gene_id)
tss = as.integer(merged1$tss)
tts = as.integer(merged1$tts)
sv.pos = as.integer(merged1$sv.bedpe.ranges@ranges@start)
sv.strand = as.character(merged1$svstrand)
sv.id = as.character(merged1$sv_id)
sv.pe.support = as.integer(merged1$pe_support)
sv.class = as.character(merged1$svclass)
sv.method = as.character(merged1$svmethod)
sv.pair.chrom = as.character(merged1$pair.chr)
sv.pair.pos = as.integer(merged1$pair.pos)
aliq_id = as.character(merged1$aliq_id)

if(opt$verbose){
  print(length(chrom))
  print(length(gene.start))
  print(length(gene.end))
  print(length(strand))
  print(length(type))
  print(length(gene))
  print(length(tss))
  print(length(tts))
  print(length(sv.pos))
  print(length(sv.strand))
  print(length(sv.id))
  print(length(sv.pe.support))
  print(length(sv.class))
  print(length(sv.method))
  print(length(sv.pair.chrom))
  print(length(sv.pair.pos))
  print(length(aliq_id))
}



merged1 <- data.frame(chrom = chrom,
                      gene.start = gene.start,
                      gene.end = gene.end,
                      strand = strand,
                      type = type,
                      gene = gene,
                      tss = tss,
                      tts = tts,
                      sv.pos = sv.pos,
                      sv.strand = sv.strand,
                      sv.id = sv.id,
                      sv.pe.support = sv.pe.support,
                      sv.class = sv.class,
                      sv.method = sv.method,
                      sv.pair.chrom = sv.pair.chrom,
                      sv.pair.pos = sv.pair.pos,
                      aliq_id = aliq_id,
                      stringsAsFactors = FALSE)
# colnames(merged1) <- c("chrom", "tss.range.start", "tss.range.end", "strand", "type", "gene", "sv.pos", "aliq_id" )
if(opt$verbose){print(head(merged1))}

#write.table(merged1, file = paste0(opt$write, opt$prefix, ".sv_merged.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


# # subset at the gene level
# # turn this section as an input for the users to choose
merged.gene <- merged1[merged1$type == "gene",]

# select for SV-gene pairs where the gene is kept intact 
x <- merged.gene

keep <- c()
for (i in 1:nrow(x)){
  if(x[i,"sv.strand"] == "+" & (x[i,"sv.pos"] >= x[i,"tss"])) {
    keep <- c(keep,i)
    } else if (x[i,"sv.strand"] == "-" & (x[i,"sv.pos"] <= x[i,"tss"])) {
      keep <- c(keep,i)
      } else {
        next
      }
}
if(opt$verbose){print(head(keep))}
merged1 <- merged.gene[keep,]
if(opt$verbose){print(head(merged1))}

# filter out SV-gene pairs where the entire SV falls within the gene body
if(opt$ingene == 1)
{
	print("Removing gene body SVs...")
	remove <- c()
	for (i in 1:nrow(merged1))
	{
		if (merged1[i,"chrom"] == merged1[i,"sv.pair.chrom"])
		{
 	  	sv.min <- min(merged1[i,"sv.pos"],merged1[i,"sv.pair.pos"])
 	  	sv.max <- max(merged1[i,"sv.pos"],merged1[i,"sv.pair.pos"])
 	  	gene.min <- min(merged1[i,"gene.start"],merged1[i,"gene.end"])
 	 		gene.max <- max(merged1[i,"gene.start"],merged1[i,"gene.end"])
 	  	if(sv.min >= gene.min & sv.max <= gene.max)
 	  	{
 	  		remove <- c(remove,i)
 	  	}
	  }
	}
	if(opt$verbose){print(head(remove))}
	bp.removed <- merged1[remove,]
	merged1 <- merged1[-remove,]
	if(opt$verbose){print(head(merged1))}
	if(opt$verbose){print(head(bp.removed))}
	#write.table(bp.removed, file = paste0(opt$write, opt$prefix, ".genebodySV_removed.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
}

print("Isolating unique SV-gene pairs...")
merged1 <- unique(merged1)
if(opt$verbose){print(head(merged1))}

merged1$sv.stat <- rep("sv", nrow(merged1))
if(opt$verbose){print(head(merged1))}

# fatten by aliquot id
print("Fattening data frame by aliquot id ...")
merged1 <- spread(data = merged1, key = aliq_id, value = sv.stat, fill = "no_sv")
if(opt$verbose){print(head(merged1))}

# add the patients with no SVs at all, as "no_SV" to the data table 
temp <- colnames(merged1)[17:ncol(merged1)]
allwt <- aliqid_wgs[!(aliqid_wgs %in% temp)]

if (!(length(allwt)==0)){
  allwt_mat <- matrix("no_sv", ncol = length(allwt), nrow = nrow(merged1))
  colnames(allwt_mat) <- allwt
  merged1 <- cbind(merged1, allwt_mat)
}


#######################################################################
########## Collapse and Count per ENSG id #############################
#######################################################################

# collapse sv calls per gene
print("Collapsing SV calls per gene (ENSG)...")
data <- as.data.frame(unclass(merged1),stringsAsFactors = TRUE, check.names = FALSE)

# set the levels of the factor to "no_sv" and "sv"
for (i in 17:ncol(data)){
  levels(data[,i]) <- c("no_sv", "sv")
}

is.sv <- function(x){
  ifelse (table(x)[2] > 0, "sv", "no_sv")  
}

data_collapsed <- data %>% group_by(gene) %>% summarise_at(colnames(data)[17:ncol(data)], is.sv)
if(opt$verbose){print(head(data_collapsed))}

num.sv <- rowSums(data_collapsed == "sv")
num.nosv <- rowSums((data_collapsed != 'sv'))-1 # don't include the first columns in the count
data_collapsed$num.sv <- num.sv
data_collapsed$num.nosv <- num.nosv
if(opt$verbose){print(head(data_collapsed))}

if (file.exists(opt$write))
{
	write.table(data_collapsed, file = paste0(opt$write, opt$prefix, ".sv_mapped_filtered_numSV.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
}else
{
	dir.create(opt$write)
	write.table(data_collapsed, file = paste0(opt$write, opt$prefix, ".sv_mapped_filtered_numSV.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
}



options(warn = oldw)