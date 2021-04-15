#!/bin/env Rscript

library(optparse)
library(data.table)
library(stringr)
library(stats)
library(MatrixGenerics)


option_list <- list(
  make_option(c("-q", "--qnormexprs"), action="store", default="qnorm_exprs.txt",
              help="Quantile normalized FPKM-UQ expression data file [default %default] to be used as the input"),
  make_option(c("-s", "--sv"), action="store", default= "sv_mapped.txt", 
              help="File with SV breakpoints mapped to gene loci [default %default]"),
  make_option(c("-i", "--id"), action="store", default= "sample_ids.txt", 
              help="File with donor, wgs and rnaseq sample ids [default %default]"),
  make_option(c("-a", "--annot"), action="store", default= "transctipts_annot.txt", 
              help="Annotation file with transcripts ids and symbols [default %default based on Gencode ]"),
  make_option(c("-p", "--purity"), action="store", default= "purity.txt", 
              help="Purity of the samples [default %default]"),
  make_option(c("-c", "--cna"), action="store", default= "cna.txt", 
              help="Gene level copy number data [default %default] compiled from all samples"),
  make_option(c("-m", "--clinical"), action="store", default= "clindat.txt", 
              help="Clinical metadata for the sample donors [default %default]"),
  make_option(c("-n", "--npc"), action="store", type="integer", default= 5, 
              help="Clinical metadata for the sample donors [default %default]"),
  make_option(c("-d", "--dir"), action="store", default= getwd(), 
              help="Input directory"),
  make_option(c("-w", "--write"), action="store", default= getwd(), 
              help="Output directory"),
  make_option(c("--pur"), action="store_true", default=FALSE,
              help="Include purity in the linear model"),
  make_option(c("--cn"), action="store_true", default=FALSE, 
              help="Include copy number in the linear model"),
  make_option(c("--age"), action="store_true", default=FALSE,
              help="Include age in the linear model"),
  make_option(c("--sex"), action="store_true", default=FALSE,
              help="Include sex in the linear model"),
  make_option(c("--PC"), action="store_true", default=FALSE,
              help="Include principal components in the linear model"),
  make_option(c("-f", "--fcutoff"), action="store", type="integer", default= 5, 
              help="SV Frequncy (%) cut-off [default %default]"),
  make_option(c("-C", "--cncutoff"), action="store", type="integer", default= 8, 
            help="Copy number cut-off for each gene [default %default]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="")
args <- parse_args(parser, positional_arguments=0)
opt <- args$options

setwd(opt$dir)
getwd()

# Load the input expression, SV, purity, CNA, and patient data ####
fpkm.uq.qn <- fread(opt$qnormexprs, stringsAsFactors = FALSE)
fpkm.uq.qn <- as.matrix(fpkm.uq.qn)
print("fpkm-uq-qn loaded:")
fpkm.uq.qn[1:5,1:5]

dat.svmapped <- fread(opt$sv, stringsAsFactors = FALSE)
print("dat.svmapped loaded:")
dat.svmapped[1:5,1:5]

samples <- fread(opt$id)
print("samples loaded:")
head(samples)

transcripts <- fread(opt$annot, stringsAsFactors = FALSE)
print("transcripts loaded:")
head(transcripts)

if(opt$pur) {
  purity <- read.table(opt$purity, stringsAsFactors = FALSE, header = TRUE)
  print("purity loaded:")
  head(purity)
}

if(opt$cn) {
  cna <- fread(opt$cna, stringsAsFactors = FALSE)
  print("cna loaded:")
  cna[1:5,1:5] 
}

if(opt$age | opt$sex) {
  clindat <- fread(opt$clinical, stringsAsFactors = FALSE)
  print("clindat loaded:")
  clindat[1:5,1:10]
}


# Prepare and QC the expression data ####
# annotate rows
rownames(fpkm.uq.qn) <- fpkm.uq.qn[,1]
fpkm.uq.qn <- fpkm.uq.qn[,-1]
class(fpkm.uq.qn) <- "numeric"

# remove genes that are 0 for all patients (this line is potentially obsolete thanks to 'add_noise.R')
fpkm.uq.qn <- fpkm.uq.qn[rowSums(fpkm.uq.qn[]) > 0 ,]

# transpose the expression matrix so that the variable (genes) are in columns
fpkm.uq.qn.t <- t(fpkm.uq.qn)
print("fpkm-uq-qn transposed:")
fpkm.uq.qn.t[1:5,1:5]

# prepare RLE values
fpkm.uq.qn.log <- log(fpkm.uq.qn + 1)
fpkm.uq.qn.log <- as.matrix(fpkm.uq.qn.log)
features_meds <- rowMedians(fpkm.uq.qn.log)
med_devs <- fpkm.uq.qn.log - features_meds
print("RLE values calculated:")
med_devs[1:5,1:5]

# plot RLE
pdf("rle.pdf")
print(boxplot(med_devs, 
              outline=FALSE, 
              xaxt = "n", 
              xlab = paste0("Samples(", ncol(med_devs) ,")"), 
              ylab = "relative log expression", 
              main = "RLE plot"))
dev.off()


# Prepare the SV data ####
print("dat.svmapped is being processed...")
# prep and transpose the data
dat.svmapped <- t(dat.svmapped) # transpose
colnames(dat.svmapped) <- dat.svmapped[1,] # make the ENSG names the column names
dat.svmapped <- dat.svmapped[-c(1, nrow(dat.svmapped)-1, nrow(dat.svmapped)),] # remove first and last two rows

# match the order of rows in dat.svmapped to the order of samples in dat.expression according to the samples table
aliqid_wgs <- samples$aliquot_id_wgs
aliqid_rnaseq <- samples$aliquot_id_rnaseq

fpkm.uq.qn.t <- fpkm.uq.qn.t[match(aliqid_rnaseq, rownames(fpkm.uq.qn.t)),]

dat.svmapped <- dat.svmapped[match(aliqid_wgs, rownames(dat.svmapped)),]
dat.svmapped <- data.frame(dat.svmapped, stringsAsFactors = FALSE)

# only keep the genes that have expression values
genes <- colnames(dat.svmapped)
keep <- genes[genes %in% rownames(fpkm.uq.qn)]
dat.svmapped <- dat.svmapped[,keep]

dim(dat.svmapped)
dat.svmapped[1:5,1:5]

# remove genes that have only single factor level for sv_status 
print("removing single factor genes from dat.svmapped")

check.levels <- function(x) {
  length(levels(as.factor(x))) >= 2
}

keep <- apply(dat.svmapped, 2, check.levels)
dat.svmapped <- dat.svmapped[,keep]

print("dat.svmapped processed:")
print(dat.svmapped[1:5,1:5])

write.table(dat.svmapped, "dat.svmapped.txt", quote = FALSE, sep = "\t")

# Prepare the purity data ####
if(opt$pur) {
  purity <- purity[purity$samplename %in% aliqid_wgs,]
  purity <- purity[match(aliqid_wgs, purity$samplename),]
  print("purity processed:")
  head(purity)
}


# Prepare the CNA data ####
if(opt$cn) {
  cna <- t(cna)
  colnames(cna) <- cna[1,]
  cna <- cna[-c(1,2),]
  cna <- cna[rownames(cna) %in% aliqid_wgs,]
  cna <- cna[match(aliqid_wgs, rownames(cna)),]
  print("cna processed:")
  cna[1:5,1:5]
}

# Prepare the Clinical data ####
if(opt$age | opt$sex) {
  clindat <- subset(clindat, select = c(donor_unique_id, project_code, icgc_donor_id, submitted_donor_id, tcga_donor_uuid,
                                        donor_sex, donor_age_at_diagnosis))
  clindat <- merge(clindat, samples, by.x = "tcga_donor_uuid", by.y = "submitter_donor_id")
  clindat <- clindat[match(aliqid_wgs, clindat$aliquot_id_wgs),]
  print("clindat processed:")
  head(clindat)
}

# Transform expression data to normal scores
ns_transform <- function (y) {
  y <- qqnorm(y, plot.it = FALSE)$x
}

temp <- rownames(fpkm.uq.qn.t)
fpkm.uq.qn.t <- apply(fpkm.uq.qn.t, 2, ns_transform)
rownames(fpkm.uq.qn.t) <- temp
print("PCA input:")
fpkm.uq.qn.t[1:5,1:5]

# Perform  PCA using prcomp() ####
pca.prcomp <- prcomp(fpkm.uq.qn.t, center = FALSE, scale. = FALSE)
var_explained <- (pca.prcomp$sdev^2/sum(pca.prcomp$sdev^2))*100

pdf("var_explained.pdf")
print(barplot(var_explained[1:10], ylab = "% var explained", names = colnames(pca.prcomp$x)[1:10]))
dev.off()

pdf("pca.pdf")
print(plot(pca.prcomp$x[,1:2]))
dev.off()

pca.ind <- as.data.frame(pca.prcomp$x)
pca.ind <- pca.ind[match(aliqid_rnaseq, rownames(pca.ind)),]

write.table(pca.ind, "pca_ind.txt", quote = FALSE, sep = "\t")
write.table(as.data.frame(pca.prcomp$x), "prcomp_centered_scaled_x.txt", quote = FALSE, sep = "\t")
write.table(pca.prcomp$sdev, "prcomp_centered_scaled_sdev.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(pca.prcomp$center, "prcomp_centered_scaled_center.txt", quote = FALSE, sep = "\t", col.names = FALSE)
write.table(pca.prcomp$scale, "prcomp_centered_scaled_scale.txt", quote = FALSE, sep = "\t", col.names = FALSE)
write.table(as.data.frame(pca.prcomp$rotation), "prcomp_centered_scaled_rotation.txt", quote = FALSE, sep = "\t")
write.table(summary(pca.prcomp)$importance, "prcomp_centered_scaled_importance.txt", quote = FALSE, sep = "\t")

print("PCA analysis completed!")


# Test PCs against the SV variable and remove PCs that correlate with SV status based on two-sided ranksum test ####
Gene <- c()
Freq <- c()
Ratio <- c()
Estimate <- c()
StdError <- c()
t <- c()
pvalue <- c()

for (i in 1:ncol(dat.svmapped)) {
  gene_id <- colnames(dat.svmapped)[i]
  print(paste0("Testing gene ", gene_id, ": i = ", i))
  
  sv_status <- dat.svmapped[,gene_id]
  print(paste0("Gene ", gene_id, " levels: ",levels(as.factor(sv_status))))
  print(sv_status)
  
  check.pca <- function(x) {
    test <- cbind(sv_status,x)
    test <- as.data.frame(test)
    test$sv_status <- as.factor(test$sv_status)
    #print(test)
    result <- wilcox.test(as.numeric(as.character(test$x)) ~ test$sv_status, data = test)
    result$p.value
  }
  
  test.result <- apply(pca.ind, 2, check.pca)
  
  keep <- test.result[!(test.result <= 0.1)] # vector with PCs that DO NOT significantly correlate with the SV variable
  keep <- names(keep)
  pca.ind.temp <- pca.ind[,keep]
  
  
  
  # Perform generalized linear model while correcting for principal components ####
  test.mat <- cbind(fpkm.uq.qn.t[,gene_id], sv_status)
  test.mat <- as.data.frame(test.mat)
  colnames(test.mat)[1] <- "expression"
  test.mat$expression <- as.numeric(as.character(test.mat$expression))
  
  
  if(opt$pur) {
    pur <- purity$purity
    if(length(levels(as.factor(pur))) >= 2) {
      test.mat$pur <- as.numeric(as.character(pur))
    }
  }
  
  if(opt$cn) {
    cn <- as.numeric(as.character(cna[,gene_id]))
    if(!is.na(sum(cn))) { # do not include cn in the model if all "NA"s
      if (length(levels(as.factor(cn)) >= 2)) {
        test.mat$cn <- as.numeric(as.character(cn))
      }
    }
  }
  
  if(opt$age) {
    age <- as.numeric(as.character(clindat$donor_age_at_diagnosis))
    if(length(levels(as.factor(age))) >= 2) {
      test.mat$age <- as.numeric(as.character(age))
    }
  }
  
  if(opt$sex) {
    sex <- clindat$donor_sex
    if(length(levels(as.factor(sex))) >= 2) {
      test.mat$sex<- sex
    }
  }
  
  
  if(opt$PC) {
    n <- opt$npc # number of first n PCs to be included in the linear model
    test.mat <- cbind(test.mat, pca.ind.temp[,1:n])
  }
  
  test.mat2 <- na.omit(test.mat) # remove all rows with NA
  test.mat2 <- test.mat2[test.mat2$expression != "Inf",] # remove patients with Inf
  print(test.mat2)
  
  if(table(test.mat2$sv_status)[1] >= (nrow(test.mat2)-2)){
    next
  } # There needs to be at least 3 sv samples
  
  if(table(test.mat2$sv_status)[1] == 0){
    next
  } # if all samples are sv, no comparison can be made
  
  # Apply CN cutoff
  if(opt$cn) {
    if("cn" %in% colnames(test.mat2)){
      if(nrow(test.mat2[test.mat2$cn <= opt$cncutoff,]) >= 10) {
        test.mat2 <- test.mat2[test.mat2$cn <= opt$cncutoff,]
      }
    }
  }

  # If all samples are same sex after sample trimming, remove the sex column from matrix
  if(opt$sex) {
    if(length(levels(as.factor(test.mat2$sex))) < 2) {
      test.mat2 <- subset(test.mat2, select = -c(sex))
    }
  }

  print(test.mat2)

  # Calculate sv ratio and frequency
  if(length(which(test.mat2$sv_status %in% "no_sv")) == 0) {
    no_sv = 0
  } else if(length(which(test.mat2$sv_status %in% "sv")) == 0) {
    sv = 0
  } else {
    sv <- length(which(test.mat2$sv_status %in% "sv"))
    no_sv <- length(which(test.mat2$sv_status %in% "no_sv"))
  }
  
  ratio <- paste0(sv, "/", no_sv)
  freq <- round((sv/(sv+no_sv))*100, 1)
  freq2 <- round((no_sv/(sv+no_sv))*100, 1)
  
  # Apply frequncy cutoff for sv (set by the user)
  if(opt$fcutoff) {
    if(freq < opt$fcutoff){
      next
    }
  }
  
  # Apply frequency cutoff for no_sv (set to 10%)
  if(freq2 < 10){
    next
  }
  

  # Run the linear model (glm() default setting is the same as lm() )
  glm <- glm(expression ~ ., data = test.mat2)
  summary <- summary(glm)
  print(summary)
  
  #extract the pval for sv_status
  Ratio <- c(Ratio, ratio)
  Freq <- c(Freq, freq)
  print(paste0("SV ratio and percent frequency for gene ", gene_id, ": ", ratio, ", ", freq, "%"))
  
  estimate <- summary$coefficients[rownames(summary$coefficients) == "sv_statussv", 1]
  error <- summary$coefficients[rownames(summary$coefficients) == "sv_statussv", 2]
  tval <-  summary$coefficients[rownames(summary$coefficients) == "sv_statussv", 3]
  pval <- summary$coefficients[rownames(summary$coefficients) == "sv_statussv", 4]
  
  print(paste0("Extracted estimate: ", estimate))
  print(paste0("Extracted error: ", error))
  print(paste0("Extracted tval: ", tval))
  print(paste0("Extracted pval: ", pval))
  
  Gene <- c(Gene, gene_id)
  Estimate <- c(Estimate, estimate)
  StdError <- c(StdError, error)
  t <- c(t, tval)
  pvalue <- c(pvalue, pval)
  
  print("Running gene_id:")
  print(tail(Gene,3))
  print("Running estimate:")
  print(tail(Estimate,3))
  print("Running error:")
  print(tail(StdError,3))
  print("Running tval:")
  print(tail(t,3))
  print("Running pval:")
  print(tail(pvalue,3))
  
  print(length(Gene))
  print(length(Estimate))
  print(length(StdError))
  print(length(t))
  print(length(pvalue))
}

print(paste0(i-1, " iterations of the loop have been completed successfully; the loop stalled at i = ", i))


# Compile results into a table ####
# Make a results table 
setwd(opt$write)
getwd()

results <- cbind(Gene, Ratio, Freq, Estimate, StdError, t, pvalue)
results <- as.data.frame(results)
results$Freq <- as.numeric(results$Freq)
results$Estimate <- as.numeric(results$Estimate)
results$StdError <- as.numeric(results$StdError)
results$t <- as.numeric(results$t)
results$pvalue <- as.numeric(results$pvalue)
str(results)
print("Results:")
head(results)
print(paste0("number of genes tested: ", nrow(results)))

outputName1 <- "results_sv"
if(opt$pur) {outputName1 <- paste0(outputName1,"+pur")}
if(opt$cn) {outputName1 <- paste0(outputName1,"+cn")}
if(opt$age) {outputName1 <- paste0(outputName1,"+age")}
if(opt$sex) {outputName1 <- paste0(outputName1,"+sex")}
if(opt$PC){outputName1 <- paste0(outputName1,"+PC",opt$npc)}
outputName1 <- paste0(outputName1, ".txt")

write.table(results, outputName1, quote = FALSE, sep = "\t", row.names = FALSE)


# Annotate and order genes ####
colnames(results)[1] <- "gene_id"
results <- merge(transcripts, results, by = "gene_id", all.y = TRUE)

# Separate Oncogenes from Suppressors and Perform FDR correction ####
# divide p-value by 2 to calculate one-sided p-value
results$p.onesided <- (results$pvalue)/2

# remove genes with negative Estimate (these are putative tumor suppressors)
nrow(results)
results.posEstimate <- results[results$Estimate > 0,] #oncogene candidates
nrow(results.posEstimate)

results.negEstimate <- results[results$Estimate < 0,] #suppressor candidates
nrow(results.negEstimate)

# Calculate FDR and order
results.posEstimate$fdr.onesided <- p.adjust(results.posEstimate$p.onesided, method = "fdr")
results.negEstimate$fdr.onesided <- p.adjust(results.negEstimate$p.onesided, method = "fdr")

results.posEstimate <- results.posEstimate[order(results.posEstimate$fdr.onesided),]
results.negEstimate <- results.negEstimate[order(results.negEstimate$fdr.onesided),]

print("Top oncogene candidates:")
head(results.posEstimate)
print("Top suppressor candidates:")
head(results.negEstimate)

# Output
outputName <- "results_annot_sv"
if(opt$pur) {outputName <- paste0(outputName,"+pur")}
if(opt$cn) {outputName <- paste0(outputName,"+cn")}
if(opt$age) {outputName <- paste0(outputName,"+age")}
if(opt$sex) {outputName <- paste0(outputName,"+sex")}
if(opt$PC){outputName <- paste0(outputName,"+PC",opt$npc)}
outputName.posEstimate <- paste0(outputName, "_posEstimate.txt")
outputName.negEstimate <- paste0(outputName, "_negEstimate.txt")

write.table(results.posEstimate, outputName.posEstimate, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(results.negEstimate, outputName.negEstimate, quote = FALSE, sep = "\t", row.names = FALSE)
