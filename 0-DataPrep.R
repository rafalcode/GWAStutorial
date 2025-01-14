#!/usr/bin/env Rscript
# 
if(!require("snpStats")) {
      stop("snpStats is required for this script. Please install it on your system.")
}

# we want to read in a bim/bed/fam file triple, create a 3 element string vector for the path to each one.
pathM <- paste0("Genomics/108Malay_2527458snps", c(".bed", ".bim", ".fam"))
pathI <- paste0("Genomics/105Indian_2527458snps", c(".bed", ".bim", ".fam"))
pathC <- paste0("Genomics/110Chinese_2527458snps", c(".bed", ".bim", ".fam"))

# now use the snpStats package to read them into three objects defined by variable names.
SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])
SNP_I <- read.plink(pathI[1], pathI[2], pathI[3])
SNP_C <- read.plink(pathC[1], pathC[2], pathC[3])


if( ncol(SNP_C) != ncol(SNP_I)) {
    stop("Different number of columns in input files detected. This is not allowed.")
}
if( ncol(SNP_I) != ncol(SNP_M)) {
    stop("Different number of columns in input files detected. This is not allowed.")
}

# append the $genotypes element from each object row-wise (we've checked columns are equal in number)
SNP <- rbind(SNP_M$genotypes, SNP_I$genotypes, SNP_C$genotypes)

# Take one bim map (all 3 maps are based on the same ordered set of SNPs)
map <- SNP_M$map
colnames(map) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")

# Rename SNPs present in the conversion table into rs IDs
load("conversionTable.RData")
mappedSNPs <- intersect(map$SNP, names(conversionTable))
newIDs <- conversionTable[match(map$SNP[map$SNP %in% mappedSNPs], names(conversionTable))]
map$SNP[rownames(map) %in% mappedSNPs] <- newIDs

# Load lipid datasets & match SNP-Lipidomics samples
lipidsMalay <- read.delim("Lipidomic/117Malay_282lipids.txt", row.names = 1)
lipidsIndian <- read.delim("Lipidomic/120Indian_282lipids.txt", row.names = 1)
lipidsChinese <- read.delim("Lipidomic/122Chinese_282lipids.txt", row.names = 1)

all(Reduce(intersect, list(colnames(lipidsMalay),
                           colnames(lipidsIndian),
                           colnames(lipidsChinese))) == colnames(lipidsMalay)) # TRUE
lip <- rbind(lipidsMalay, lipidsIndian, lipidsChinese)

matchingSamples <- intersect(rownames(lip), rownames(SNP))
SNP <- SNP[matchingSamples,]
lip <- lip[matchingSamples,]

genData <- list(SNP = SNP, MAP = map, LIP = lip)
save(genData, file = "PhenoGenoMap.RData")

# Clear memory
rm(list = ls())
