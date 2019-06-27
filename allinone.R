#!/usr/bin/env Rscript
# packages prereqs, check they are all there
source("pkg.prereqs")
# hardcode for time being, need to get rid of this later
prepath <- c("/home/rf73/neiproj/public/")

# we want to read in a bim/bed/fam file triple, create a 3 element string vector for the path to each one.
pathM <- paste0(prepath[1], "Genomics/108Malay_2527458snps", c(".bed", ".bim", ".fam"))
pathI <- paste0(prepath[1], "Genomics/105Indian_2527458snps", c(".bed", ".bim", ".fam"))
pathC <- paste0(prepath[1], "Genomics/110Chinese_2527458snps", c(".bed", ".bim", ".fam"))
if(!file.exists(pathM[1])) {
    stop("some of your data files don't exist")
}

# now use the snpStats package to read them into three objects defined by variable names.
SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])
SNP_I <- read.plink(pathI[1], pathI[2], pathI[3])
SNP_C <- read.plink(pathC[1], pathC[2], pathC[3])


if( ncol(SNP_C$genotypes) != ncol(SNP_I$genotypes)) {
    stop("Different number of columns in input files detected. This is not allowed.")
}
if( ncol(SNP_I$genotypes) != ncol(SNP_M$genotypes)) {
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

# Load lipid datasets
lipidsMalay <- read.delim(paste0(prepath[1], "Lipidomic/117Malay_282lipids.txt"), row.names = 1)
lipidsIndian <- read.delim( paste0(prepath[1], "Lipidomic/120Indian_282lipids.txt"), row.names = 1)
lipidsChinese <- read.delim(paste0(prepath[1], "Lipidomic/122Chinese_282lipids.txt"), row.names = 1)
# match SNP-Lipidomics samples
all(Reduce(intersect, list(colnames(lipidsMalay),
                           colnames(lipidsIndian),
                           colnames(lipidsChinese))) == colnames(lipidsMalay)) # TRUE

lip <- rbind(lipidsMalay, lipidsIndian, lipidsChinese)

matchingSamples <- intersect(rownames(lip), rownames(SNP))
SNP <- SNP[matchingSamples,]
lip <- lip[matchingSamples,]

genData <- list(SNP = SNP, MAP = map, LIP = lip)
save(genData, file = "PhenoGenoMap.RData")

#get GDS part
# Merge all files and write as PLINK again; we need it in order to correct for kinship
SNP_M$genotypes <- rbind(SNP_M$genotypes, SNP_I$genotypes, SNP_C$genotypes)
colnames(map) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
SNP_M$fam<- rbind(SNP_M$fam, SNP_I$fam, SNP_C$fam)

# Rename SNPs present in the conversion table into rs IDs
mappedSNPs <- intersect(SNP_M$map$SNP, names(conversionTable))
newIDs <- conversionTable[match(SNP_M$map$SNP[SNP_M$map$SNP %in% mappedSNPs], names(conversionTable))]
SNP_M$map$SNP[rownames(SNP_M$map) %in% mappedSNPs] <- newIDs

write.plink("convertGDS", snps = SNP_M$genotypes)

# part 3
source("GWASfunction.R")
maf <- 0.1
callRate <- 1
SNPstats <- col.summary(genData$SNP)

maf_call <- with(SNPstats, MAF > maf & Call.rate == callRate)
genData$SNP <- genData$SNP[,maf_call]
genData$MAP <- genData$MAP[maf_call,]
SNPstats <- SNPstats[maf_call,]

# Sample call rate & heterozygosity
callMat <- !is.na(genData$SNP)
Sampstats <- row.summary(genData$SNP)
hetExp <- callMat %*% (2 * SNPstats$MAF * (1 - SNPstats$MAF)) # Hardy-Weinberg heterozygosity (expected)
hetObs <- with(Sampstats, Heterozygosity * (ncol(genData$SNP)) * Call.rate)
Sampstats$hetF <- 1-(hetObs/hetExp)
# Use sample call rate of 100%, het threshold of 0.1 (very stringent)
het <- 0.1 # Set cutoff for inbreeding coefficient;
het_call <- with(Sampstats, abs(hetF) < het & Call.rate == 1)
genData$SNP <- genData$SNP[het_call,]
genData$LIP <- genData$LIP[het_call,]

# LD and kinship coeff
ld <- .2
kin <- .1
cat(paste0("About to run via snpgdsBED2GDS() ...\n"))
snpgdsBED2GDS(bed.fn = "convertGDS.bed", bim.fn = "convertGDS.bim", fam.fn = "convertGDS.fam", out.gdsfn = "myGDS", cvt.chr = "char")
genofile <- snpgdsOpen("myGDS", readonly = F)
gds.ids <- read.gdsn(index.gdsn(genofile,  "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = T)
geno.sample.ids <- rownames(genData$SNP)

# First filter for LD
cat(paste0("About to filter via snpgdsLDpruning() ...\n"))
# quite a demanding function this one, need the threads (not in original script)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld, sample.id = geno.sample.ids, snp.id = colnames(genData$SNP), num.thread = 4, verbose = TRUE)
snpset.ibd <- unlist(snpSUB, use.names = F)
# And now filter for MoM
cat(paste0("About to filter via snpgdsIBDMoM() ...\n"))
ibd <- snpgdsIBDMoM(genofile, kinship = T, sample.id = geno.sample.ids, snp.id = snpset.ibd, num.thread = 4)
cat(paste0("Completed snpgdsIBDMoM(), now snpgdsIBDSelection() ...\n"))
ibdcoef <- snpgdsIBDSelection(ibd)
cat(paste0("Completed snpgdsIBDSelection().\n"))
ibdcoef <- ibdcoef[ibdcoef$kinship >= kin,]

# worked fine up to here.

# Filter samples out
related.samples <- NULL
cat(paste0("About to start snpgdsPCA() ...\n"))
cat(paste0("ibdcoef came out with ", nrow(ibdcoef), "number of rows\n"))
while (nrow(ibdcoef) > 0) {
      # count the number of occurrences of each and take the top one
      sample.counts <- arrange(count(c(ibdcoef$ID1, ibdcoef$ID2)), -freq)
      rm.sample <- sample.counts[1, 'x']
      cat("Removing sample ", as.character(rm.sample), " as it is too closely related to ", sample.counts[1, 'freq'], " other samples.\n")
      # remove from ibdcoef and add to list
      ibdcoef <- ibdcoef[ibdcoef$ID1 != rm.sample & ibdcoef$ID2 != rm.sample,]
      related.samples <- c(as.character(rm.sample), related.samples)
}
genData$SNP <- genData$SNP[!(rownames(genData$SNP) %in% related.samples),]
genData$LIP <- genData$LIP[!(rownames(genData$LIP) %in% related.samples),]

# PCA
cat(paste0("About to start snpgdsPCA() ...\n"))
pcaout <- snpgdsPCA(genofile, sample.id = geno.sample.ids, snp.id = snpset.ibd, num.thread = 4)

cat(paste0("Completed snpgdsPCA(), rendering into data.frame ...\n"))
pctab <- data.frame(sample.id = pcaout$sample.id, PC1 = pcaout$eigenvect[,1], PC2 = pcaout$eigenvect[,2], stringsAsFactors = FALSE)

origin <- read.delim("countryOrigin.txt", sep = "\t")
origin <- origin[match(pcaout$sample.id, origin$sample.id),]

pcaCol <- rep(rgb(0,0,0,.3), length(pcaout$sample.id)) # Set black for chinese
pcaCol[origin$Country == "I"] <- rgb(1,0,0,.3) # red for indian
pcaCol[origin$Country == "M"] <- rgb(0,.7,0,.3) # green for malay

png("PCApopulation.png", width = 500, height = 500)
plot(pctab$PC1, pctab$PC2, xlab = "PC1", ylab = "PC2", col = pcaCol, pch = 16)
abline(h = 0, v = 0, lty = 2, col = "grey")
legend("top", legend = c("Chinese", "Indian", "Malay"), col = 1:3, pch = 16, bty = "n")
dev.off()

# Choose trait for association analysis, use colnames(genData$LIP) for listing
# NOTE: Ignore the first column of genData$LIP (gender)
target <- "Cholesterol"

phenodata <- data.frame("id" = rownames(genData$LIP), "phenotype" = scale(genData$LIP[,target]), stringsAsFactors = FALSE)

# Conduct GWAS (will take a while)
cat(paste0("About to call GWAA() function ...\n"))
start <- Sys.time()
GWAA(genodata = genData$SNP, phenodata = phenodata, filename = paste0(target, ".txt"))
Sys.time() - start # benchmark
cat(paste0("Completed GWAA() function.\n"))

# Manhattan plot
GWASout <- read.table(paste(target, ".txt", sep = ""), header = T, colClasses = c("character", rep("numeric",4)))
GWASout$type <- rep("typed", nrow(GWASout))
GWASout$Neg_logP <- -log10(GWASout$p.value)
GWASout <- merge(GWASout, genData$MAP[,c("SNP", "chr", "position")])
GWASout <- GWASout[order(GWASout$Neg_logP, decreasing = T),]

png(paste(target, ".png", sep = ""), height = 500,width = 1000)
cat(paste0("About to call GWAS_Manhattan() function ...\n"))
GWAS_Manhattan(GWASout)
dev.off()
cat(paste0("Completed GWAS_Manhattan().\n"))

# QQ plot using GenABEL estlambda function
png(paste(target, "_QQplot.png", sep = ""), width = 500, height = 500)
lambda <- estlambda(GWASout$t.value**2, plot = T, method = "median")
dev.off()
