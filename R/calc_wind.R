#!/bin/env

oldw <- getOption("warn")
options(warn = -1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))

option_list <- list( 
  make_option(c("-r", "--ref"), action="store",
              help="Gene annotation"),
  make_option(c("-u", "--up"), action="store", default=500000,
              help="Number of basepairs UPSTREAM of transcription start site [default %default bp]"),
  make_option(c("-d", "--down"), action="store", default=500000,
              help="Number of basepairs DOWNSTREAM of transcription start site [default %default bp]"),
  make_option(c("-w", "--write"), action="store", default= "./tss_windows.txt", 
              help="Output file [default %default]"),
  make_option(c("-v","--verbose"), action="store_true", default=FALSE,
              help="Verbose mode for troubleshooting [default %default]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="")
args <- parse_args(parser, positional_arguments=0)
opt <- args$options

# load the reference genome
ref <- fread(opt$ref, stringsAsFactors = FALSE)

# calculate the TSS window for each gene
tss <- c()
tts <- c()
tss_left <- c()
tss_right <- c()

up.n <- as.integer(opt$up)
down.n <- as.integer(opt$down)

for (i in 1:nrow(ref)) {
  if(ref$STRAND[i] == "+") {
    start <-  as.numeric(ref$START)[i]
    end <- as.numeric(ref$END)[i]
    left <- as.numeric(ref$START)[i] - up.n
    right <- as.numeric(ref$START)[i] + down.n
    tss <- c(tss, start)
    tts <- c(tts, end)
    tss_left <- c(tss_left, left)
    tss_right <- c(tss_right, right)
  } else if (ref$STRAND[i] == "-") {
    start <-  as.numeric(ref$END)[i]
    end <- as.numeric(ref$START)[i]
    left <- as.numeric(ref$END)[i] - down.n
    right <- as.numeric(ref$END)[i] + up.n
    tss <- c(tss, start)
    tts <- c(tts, end)
    tss_left <- c(tss_left, left)
    tss_right <- c(tss_right, right)
  } else {
    start <- NA
    end <- NA
    left <- NA
    right <- NA
    tss <- c(tss, start)
    tts <- c(tts, end)
    tss_left <- c(tss_left, left)
    tss_right <- c(tss_right, right)
    warning("At least one row is missing strand information.")
  }
  
  if(opt$verbose) {
    print(i)
    print(tail(tss_left))
    print(tail(tss_right))
    print(tail(tss))
  }

}

ref$tss <- tss
ref$tts <- tts
ref$tss_left <- tss_left
ref$tss_right <- tss_right
if(opt$verbose){print(head(ref))}


if (file.exists(opt$write)) {
  print(paste0(opt$write," file already exists"))
} else {
  write.table(ref, file = opt$write, quote = FALSE, sep = "\t", row.names = FALSE)
}



options(warn = oldw)