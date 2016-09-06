###############################################################################
## CONTA
###############################################################################
# Copyright (c) 2016 Tobias Meissner, Niels Weinhold

# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
# THE SOFTWARE.

#!/usr/local/bin/r

###############################################################################
# command line options
###############################################################################
packages <- function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x, character.only=TRUE, quietly=TRUE)){
    install.packages(pkgs=x, repos="http://cran.r-project.org", quiet=TRUE)
    require(x, character.only=TRUE, quietly=TRUE)
  }
}
packages(optparse)
option_list = list(
  make_option(c("-t", "--tumor"), type="character", default=NULL, 
              help="tumor bam file(s) [REQUIRED]", metavar="character"),
  make_option(c("-g", "--germline"), type="character", default=NULL, 
              help="germline bam file(s) [REQUIRED]", metavar="character"),
  make_option(c("-p", "--panel"), type="character", default=NULL, 
              help="panel of SNPs for tracking/contaminatin estimation [REQUIRED]", 
              metavar="character"),
  make_option(c("-n", "--sampleName"), type="character", default=NULL, 
              help="sample names(s) [REQUIRED]", metavar="character"),
  make_option(c("-s", "--sampleSheet"), type="character", default=NULL, 
              help="sample sheet that lists sample name(s),
              tumor bam file location(s), germline bam file location(s) 
              [REQUIRED]", metavar="character"),
  make_option(c("-e", "--percentHom"), type="numeric", default=10, 
              help="percent homozygous [default= %default]", 
              metavar="numeric"),
  make_option(c("-r", "--minReads"), type="numeric", default=50, 
              help="minimum read depth at SNP position [default= %default]", 
              metavar="numeric"),
  make_option(c("-c", "--maxContLevelGerm"), type="numeric", default=10, 
              help="max contamination germline sample [default= %default]", 
              metavar="numeric"),
  make_option(c("-q", "--min_base_quality"), type="numeric", default=20, 
              help="minimum base quality [default= %default]", 
              metavar="numeric"),
  make_option(c("-x", "--contPerSNP"), type="numeric", default=0, 
              help="contamination per SNP [default= %default]", 
              metavar="numeric"),
  make_option(c("-a", "--aberrantSNP"), type="numeric", default=5, 
              help="number abberant SNPs [default= %default]", 
              metavar="numeric"),
  make_option(c("-z", "--aberrantSNPPercent"), type="numeric", default=10, 
              help="percent abberant SNPs [default= %default]", 
              metavar="numeric"),
  make_option(c("-m", "--mode"), type="character", default='pair', 
              help="analysis mode [default= %default]", 
              metavar="character"),
  make_option(c("-o", "--out"), type="character", default="conta.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-v", "--verbose"), type="logical", default=F, 
              help="get more detailed output [default= %default]", 
              metavar="logical")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$tumor) & is.null(opt$sampleSheet)){
  print_help(opt_parser)
  stop("(A) tumor sample(s) or sample sheet must be supplied", call.=FALSE)
}
if (is.null(opt$germline)  & is.null(opt$sampleSheet)){
  print_help(opt_parser)
  stop("(A) germline sample(s) or sample sheet must be supplied", call.=FALSE)
}
if (is.null(opt$panel)){
  print_help(opt_parser)
  stop("A SNP panel must be supplied", call.=FALSE)
}
if (is.null(opt$sampleName) & is.null(opt$sampleSheet)){
  print_help(opt_parser)
  stop("(A) sample name(s) or sample sheet must be supplied", call.=FALSE)
}

# if there are multiple tumor/germline input files, split by seperator (comma)
if(is.null(opt$sampleSheet)) {
  opt$tumor <- unlist(strsplit(opt$tumor, split=','))
  opt$germline <- unlist(strsplit(opt$germline, split=','))
  opt$sampleName <- unlist(strsplit(opt$sampleName, split=','))  
} 

###############################################################################
# load and install required packages if needed
###############################################################################

if(opt$verbose) {
  packages(Rsamtools)
  packages(VariantAnnotation)
  packages(GenomicRanges)
  packages(plyr)
} else {
  msg.trap <- capture.output(suppressMessages(packages(Rsamtools)))
  msg.trap <- capture.output(suppressMessages(packages(VariantAnnotation)))
  msg.trap <- capture.output(suppressMessages(packages(GenomicRanges)))
  msg.trap <- capture.output(suppressMessages(packages(plyr)))
}

###############################################################################
# functions
###############################################################################
.pileupFreq <- function(pileupres) {
  nucleotides <- levels(pileupres$nucleotide)
  res <- split(pileupres, pileupres$seqnames)
  res <- lapply(res, function (x) {split(x, x$pos)})
  res <- lapply(res, function (positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.character(unique(each$pos))
      tabcounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
      c(chr,pos, tabcounts)
    })
    nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) <- NULL
    nuctab
  })
  res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) <- NULL
  colnames(res) <- c("seqnames","start",levels(pileupres$nucleotide))
  res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
  res
}

.readBam <- function(bam, vcf.ranges, min_base_quality) {
  bamfile <- as.character(bam)
  bf <- BamFile(bamfile)
  param <- ScanBamParam(which = vcf.ranges)
  flag <- scanBamFlag(isDuplicate = FALSE)
  p_param <- PileupParam(max_depth = 10000,
                         min_base_quality = min_base_quality,
                         ignore_query_Ns = FALSE)
  res <- pileup(bf, 
                scanBamParam=param, 
                scanBamFlag=flag,
                pileupParam=p_param)
  res
}

.qcTab <- function(type, ranges, tab, contPerSNP, maxContLevelGerm, minReads, percentHom, aberrantSNP, aberrantSNPPercent) {
  tab$Ref_A <- as.vector(ranges$Ref)
  tab$Alt_A <- as.vector(ranges$Alt)
  tab$Ref <- NA
  tab$Alt <- NA
  tab$Ref <- ifelse(ranges$Ref == "T", tab$T, tab$Ref)
  tab$Ref <- ifelse(ranges$Ref == "A", tab$A, tab$Ref)
  tab$Ref <- ifelse(ranges$Ref == "C", tab$C, tab$Ref)
  tab$Ref <- ifelse(ranges$Ref == "G", tab$G, tab$Ref)
  tab$Alt <- ifelse(ranges$Alt == "T", tab$T, tab$Alt)
  tab$Alt <- ifelse(ranges$Alt == "A", tab$A, tab$Alt)
  tab$Alt <- ifelse(ranges$Alt == "C", tab$C, tab$Alt)
  tab$Alt <- ifelse(ranges$Alt == "G", tab$G, tab$Alt)
  
  tab$A1 <- ifelse(tab$Ref_A =="A" & tab$Alt_A == "C", tab[,5], 0)
  tab$A1 <- ifelse(tab$Ref_A =="A" & tab$Alt_A == "G", tab[,4], tab$A1)
  tab$A1 <- ifelse(tab$Ref_A =="A" & tab$Alt_A == "T", tab[,4], tab$A1)
  tab$A1 <- ifelse(tab$Ref_A =="C" & tab$Alt_A == "A", tab[,5], tab$A1)
  tab$A1 <- ifelse(tab$Ref_A =="C" & tab$Alt_A == "G", tab[,3], tab$A1)
  tab$A1 <- ifelse(tab$Ref_A =="C" & tab$Alt_A == "T", tab[,3], tab$A1)
  tab$A1 <- ifelse(tab$Ref_A =="G" & tab$Alt_A == "A", tab[,4], tab$A1)
  tab$A1 <- ifelse(tab$Ref_A =="G" & tab$Alt_A == "C", tab[,3], tab$A1)
  tab$A1 <- ifelse(tab$Ref_A =="G" & tab$Alt_A == "T", tab[,3], tab$A1)
  tab$A1 <- ifelse(tab$Ref_A =="T" & tab$Alt_A == "A", tab[,4], tab$A1)
  tab$A1 <- ifelse(tab$Ref_A =="T" & tab$Alt_A == "C", tab[,3], tab$A1)
  tab$A1 <- ifelse(tab$Ref_A =="T" & tab$Alt_A == "G", tab[,3], tab$A1)
  
  tab$A2 <- ifelse(tab$Ref_A == "A" & tab$Alt_A == "C", tab[,6], 0)
  tab$A2 <- ifelse(tab$Ref_A == "A" & tab$Alt_A == "G", tab[,6], tab$A2)
  tab$A2 <- ifelse(tab$Ref_A == "A" & tab$Alt_A == "T", tab[,5], tab$A2)
  tab$A2 <- ifelse(tab$Ref_A == "C" & tab$Alt_A == "A", tab[,6], tab$A2)
  tab$A2 <- ifelse(tab$Ref_A == "C" & tab$Alt_A == "G", tab[,6], tab$A2)
  tab$A2 <- ifelse(tab$Ref_A == "C" & tab$Alt_A == "T", tab[,5], tab$A2)
  tab$A2 <- ifelse(tab$Ref_A == "G" & tab$Alt_A == "A", tab[,6], tab$A2)
  tab$A2 <- ifelse(tab$Ref_A == "G" & tab$Alt_A == "C", tab[,6], tab$A2)
  tab$A2 <- ifelse(tab$Ref_A == "G" & tab$Alt_A == "T", tab[,4], tab$A2)
  tab$A2 <- ifelse(tab$Ref_A == "T" & tab$Alt_A == "A", tab[,5], tab$A2)
  tab$A2 <- ifelse(tab$Ref_A == "T" & tab$Alt_A == "C", tab[,5], tab$A2)
  tab$A2 <- ifelse(tab$Ref_A == "T" & tab$Alt_A == "G", tab[,4], tab$A2)
  tab <- tab[,c(1:2,11:16)]
  
  tab$P_Ref <- round(100/(as.numeric(tab$Ref) + as.numeric(tab$Alt)) * as.numeric(tab$Ref), 0)
  
  tab$minReads <- ifelse(as.numeric(tab$Ref) + as.numeric(tab$Alt) < minReads, "no", "yes")
 
  if (type=='germline') {
    tab$Hom <- ifelse(tab$P_Ref == 100 | tab$P_Ref <= percentHom, "yes", "no")
    tab$Conta <- as.numeric(ifelse(tab$P_Ref >= (100 - maxContLevelGerm), tab$Alt, tab$Ref))
    tab$Geno_A <- ifelse(tab$P_Ref >= (100 - maxContLevelGerm), tab$Ref_A, tab$Alt_A)
    tab$Conta_A <- ifelse(tab$P_Ref >= (100 - maxContLevelGerm), tab$Alt_A, tab$Ref_A)  
    
    tab <- subset(tab, Hom == 'yes' & minReads == 'yes')
    
    tab_QC <- tab
    tab_QC$Norm <- as.numeric(ifelse(tab_QC$P_Ref >= 90, tab_QC$Ref, tab_QC$Alt))
    tab_QC$ProContaAlt <- 100 / (tab_QC$Norm+tab_QC$Conta + tab_QC$A1 + tab_QC$A2) * (tab_QC$Conta)
    tab_QC$ProContaA1 <- 100 / (tab_QC$Norm+tab_QC$Conta + tab_QC$A1 + tab_QC$A2) * (tab_QC$A1)
    tab_QC$ProContaA2 <- 100 / (tab_QC$Norm+tab_QC$Conta + tab_QC$A1 + tab_QC$A2) * (tab_QC$A2)
    
    QC <- list(Type = 'Germline',
               N_Hom = dim(tab_QC)[1],
               NumSNP_Alt = length(which(as.numeric(tab_QC[ ,12]) > contPerSNP)),
               NumSNP_A1 = length(which(as.numeric(tab_QC[ ,7]) > contPerSNP)),
               NumSNP_A2 = length(which(as.numeric(tab_QC[ ,8]) > contPerSNP)),
               Sum_Alt = sum(as.numeric(tab_QC[ ,12])),
               Sum_A1 = sum(as.numeric(tab_QC[ ,7])),
               Sum_A2 = sum(as.numeric(tab_QC[ ,8])),
               Percent_Alt = 100 / dim(tab_QC)[1] * length(which(as.numeric(tab_QC[ ,12]) > contPerSNP)),
               Percent_A1 = 100 / dim(tab_QC)[1] * length(which(as.numeric(tab_QC[ ,7]) > contPerSNP)),
               Percent_A2 = 100 / dim(tab_QC)[1] * length(which(as.numeric(tab_QC[ ,8]) > contPerSNP)),
               PercentCont_Alt = sum(tab_QC$ProContaAlt[which(as.numeric(tab_QC[ ,12]) > contPerSNP)]) / dim(tab)[1],
               PercentCont_A1 = sum(tab_QC$ProContaA1[which(as.numeric(tab_QC[ ,7]) > contPerSNP)]) / dim(tab)[1],
               PercentCont_A2 = sum(tab_QC$ProContaA2[which(as.numeric(tab_QC[ ,8]) > contPerSNP)]) / dim(tab)[1]
    )
    QC$Background <- ifelse(QC$PercentCont_A1 > QC$PercentCont_A2, QC$PercentCont_A1, QC$PercentCont_A2)
    QC$EstCont <- ifelse(QC$N_Hom / 100 * (QC$Percent_Alt >= aberrantSNPPercent) & (QC$NumSNP_Alt >= aberrantSNP), 2 * (QC$PercentCont_Alt - QC$Background), 0)
    QC$EstCont <- round(ifelse(QC$EstCont < 0, 0, QC$EstCont), 2)
  }
  
  if (type=='tumor') {
    tab$Hom <- ifelse(tab$P_Ref >= 90 | tab$P_Ref <= 10 , "yes", "no")
    tab$Conta <- tab$Alt
    tab$Norm=tab$Ref
    tab$ProContaAlt <- 100 / (tab$Norm + tab$Conta + tab$A1 + tab$A2) * (tab$Conta)
    tab$ProContaA1 <-  100 / (tab$Norm + tab$Conta + tab$A1 + tab$A2) * (tab$A1)
    tab$ProContaA2 <- 100 / (tab$Norm + tab$Conta + tab$A1 + tab$A2) * (tab$A2)
    
    QC <- list(Type = 'Tumor',
               N_Hom = dim(tab)[1],
               NumSNP_Alt = length(which(as.numeric(tab[ ,12]) > contPerSNP)),
               NumSNP_A1 = length(which(as.numeric(tab[ ,7]) > contPerSNP)),
               NumSNP_A2 = length(which(as.numeric(tab[ ,8]) > contPerSNP)),
               Sum_Alt = sum(as.numeric(tab[ ,12])),
               Sum_A1 = sum(as.numeric(tab[ ,7])),
               Sum_A2 = sum(as.numeric(tab[ ,8])),
               Percent_Alt = 100 / dim(tab)[1] * length(which(as.numeric(tab[ ,12]) > contPerSNP)),
               Percent_A1 = 100 / dim(tab)[1] * length(which(as.numeric(tab[ ,7]) > contPerSNP)),
               Percent_A2 = 100 / dim(tab)[1] * length(which(as.numeric(tab[ ,8]) > contPerSNP)),
               PercentCont_Alt = sum(tab$ProContaAlt[which(as.numeric(tab[ ,12]) > contPerSNP)]) / dim(tab)[1],
               PercentCont_A1 = sum(tab$ProContaA1[which(as.numeric(tab[ ,7]) > contPerSNP)]) / dim(tab)[1],
               PercentCont_A2 = sum(tab$ProContaA2[which(as.numeric(tab[ ,8]) > contPerSNP)]) / dim(tab)[1]
    )
    QC$Background <- ifelse(QC$PercentCont_A1 > QC$PercentCont_A2, QC$PercentCont_A1, QC$PercentCont_A2)
    QC$EstCont <- ifelse(QC$N_Hom / 100 * (QC$Percent_Alt >= aberrantSNPPercent) & (QC$NumSNP_Alt >= aberrantSNP), 2 * (QC$PercentCont_Alt - QC$Background), 0)
    QC$EstCont <- round(ifelse(QC$EstCont < 0, 0, QC$EstCont), 2)
  }
  return(list(QC=QC, tab=tab))
}

estCont <- function(bamGermline, bamTumor, panel, mode='pair', percentHom=10, minReads=50, maxContLevelGerm=10,
                    min_base_quality=20, contPerSNP=0, aberrantSNP=5, aberrantSNPPercent=10) {
  if(file.exists(bamGermline) & file.exists(bamTumor)) {
    if(mode == 'pair') {
      ## GERMLINE
      vcf.ranges=GRanges(seqnames = panel$contig, 
                         ranges = IRanges(start=panel$position, end=panel$position), 
                         strand = "*", 
                         paramRangeID = NA, 
                         Ref = panel$ref_allele, 
                         Alt = panel$alt_allele, 
                         QUAL = NA, 
                         FILTER = "REJECT", 
                         SNP = panel$SNP)
      tryCatch({
        res <- .readBam(bamGermline, vcf.ranges, min_base_quality)
        tab <- .pileupFreq(res)
        
        if(dim(tab)[1] != length(vcf.ranges)) {
        vcf.ranges <- vcf.ranges[start(vcf.ranges) %in% tab$start]
      }
        
        QCGerm <- .qcTab('germline', vcf.ranges, tab, contPerSNP, maxContLevelGerm, minReads, percentHom, aberrantSNP, aberrantSNPPercent)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      ## TUMOR
      vcf.ranges2 <- GRanges(seqnames = QCGerm$tab$seqnames, 
                             ranges = IRanges(start=as.numeric(QCGerm$tab$start), end=as.numeric(QCGerm$tab$start)), 
                             strand = "*", 
                             paramRangeID = NA, 
                             Ref = QCGerm$tab$Geno_A, 
                             Alt = QCGerm$tab$Conta_A, 
                             QUAL = NA, 
                             FILTER = "REJECT")
      tryCatch({
        resTumor <- .readBam(bamTumor, vcf.ranges2, min_base_quality)
        tabTumor <- .pileupFreq(resTumor)
        
        if(dim(tabTumor)[1] != length(vcf.ranges2)) {
          vcf.ranges2 <- vcf.ranges2[start(vcf.ranges2) %in% tabTumor$start]
        }
        
        QCTumor <- .qcTab('tumor', vcf.ranges2, tabTumor, contPerSNP, maxContLevelGerm, minReads, percentHom, aberrantSNP, aberrantSNPPercent)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    } 
    
    if(mode == 'single') {
      print('No support for this, yet!')
    }
    return(list(QCGerm=QCGerm$QC,
                QCTumor=QCTumor$QC,
                TabGerm=QCGerm$tab,
                TabTumor=QCTumor$tab
                )
    )
  } else {
    return(list(QCGerm=NA,
                QCTumor=NA,
                TabGerm=NA,
                TabTumor=NA
                )
    )
  }
}

###############################################################################
# read in the SNP panel
###############################################################################
panel <- read.csv2(opt$panel, sep=",", stringsAsFactors = F)
panel <- panel[with(panel,order(contig,position)),]
panel$contig <- paste0('chr', panel$contig)

###############################################################################
# estimate contamination 
###############################################################################
if(is.null(opt$sampleSheet)) {
  samples <- data.frame(SAMPLE_ID=opt$sampleName,
                        CONTROL=opt$germline,
                        TUMOR=opt$tumor)  
} else {
  samples <- opt$sampleSheet
  samples <- read.csv2(opt$sampleSheet, sep=',', stringsAsFactors=F)
}

cat('Running analysis with the following settings .. \n')
cat(paste0('  Number of samples to be analyzed: ', dim(samples)[1], '\n'))
cat(paste0('  Mode: ', opt$mode, '\n'))
cat(paste0('  Percent Hom.: ', opt$percentHom, '\n'))
cat(paste0('  Min. Reads: ', opt$minReads, '\n'))
cat(paste0('  Max. Cont. Level Germ.: ', opt$maxContLevelGerm, '\n'))
cat(paste0('  Min. Base Qual: ', opt$min_base_quality, '\n'))
cat(paste0('  Cont. per SNP: ', opt$contPerSNP, '\n'))
cat(paste0('  Abberant SNPs: ', opt$aberrantSNP, '\n'))
cat(paste0('  Abberant SNPs Percent: ', opt$aberrantSNPPercent, '\n'))
cat('\n')

res <- apply(samples, 1, function(x) {
  cat(paste0('Analyzing Sample: '), x['SAMPLE_ID'], '\n')
  estCont(x['CONTROL'], 
          x['TUMOR'], 
          panel,
          opt$mode,
          opt$percentHom,
          opt$minReads,
          opt$maxContLevelGerm,
          opt$min_base_quality,
          opt$contPerSNP,
          opt$abberantSNP,
          opt$abberantSNPPercent
          )
})
names(res) <- samples$SAMPLE_ID

# format the output for tabular display
df <- rbind(do.call(rbind, lapply(res, '[[', 1)),
            do.call(rbind, lapply(res, '[[', 2))
)
# dont throw a warning here...
oldw <- getOption("warn")
options(warn = -1)
df1 <- data.frame(Sample=rownames(df), df, stringsAsFactors = F)
options(warn = oldw)
df1 <- arrange(df1, Sample)
df1[sapply(df1, is.list)] <- apply(df1[sapply(df1, is.list)], 2, function(x) as.vector(unlist(x)))

# writting results
cat(paste0('Writing results to '), opt$out, '\n')
write.table(df1, 
          opt$out, 
          row.names = F, 
          quote = F,
          sep='\t'
)
