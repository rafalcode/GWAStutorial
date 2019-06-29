#!/usr/bin/env Rscript
# A little workflow (pipeline) for taking a ped file and rendering it to give condition (column 6) as 1 (c0ntrol) if male (col 5 == 1),
# and col6 ==2 if col5 ==2 (i.e. case if female).

args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 1 # expected numargs
if(numargs != enumargs) {
    print("This script takes a ped file (note!) prefix name")
    stop("Stopping right here")
}

# first look for ped files, it's like though, that
inpedf <- paste0(args[1], ".ped")
foundped =0
inmapf <- paste0(args[1], ".map")
if(!file.exists(inpedf)) {
    warning(paste0("Error: cannot find the input file \"", inpedf, "\".\n"))
} else {
    foundped =1
}
if(!file.exists(inmapf)) {
    warning(paste0("Error: \"", inmapf, "\" not found. Although not required, it is a mistake to have a ped file without an accompanying map file.\n"))
}

if(!foundped) {
    print("So, no ped file was found with the that prefix, try the bed/bim/fam trio")
    inbimf <- paste0(args[1], ".bim")
    inbedf <- paste0(args[1], ".bed")
    infamf <- paste0(args[1], ".fam")

    if(!file.exists(inbedf)) {
        stop(paste0("Error: cannot find the input file \"", inbedf, "\".\n"))
    }
    if((!file.exists(inbimf)) | (!file.exists(infamf))) {
        stop(paste0("Error: \"", inbimf, "\" or \"", infamf, "\" not found. Although not required, it is a mistake to have a bed file without accompanying bim and fam files.\n"))
    }

    # so now we're going to use plink to convert the bed/bim/fam to ped
    plinkbin <- Sys.which("plink1.9")
    humopts <- " --silent --allow-no-sex"
    plinkCmd <- paste0(plinkbin, humopts, " --bfile ", args[1], " --recode --out ", args[1])
    rt <- system(plinkCmd, intern=FALSE)
    if (rt!=0) {
        stop("plink bed to ped error.")
    }
}

# OK, even if we had only bed, it's now been convert to ped
outpf <- paste0(args[1], "_femcase.ped")
outmf <- paste0(args[1], "_femcase.map")
# 
awkCall <- paste0("awk '{$6=$5; print $0}' ", inpedf, " > ", outpf)
rt <- system(awkCall, intern=F)
if (rt!=0) {
    stop("awk operation error")
}
file.copy(inmapf, outmf)

# But if it was bed, if's not over, we need to convert back
if(!foundped) {
    print("So, no ped file was found with the that prefix, try the bed/bim/fam trio")
    outfpfx <- paste0(args[1], "_femcase")
    plinkCmd <- paste0(plinkbin, humopts, " --file ", outfpfx, " --recode --out ", outfpfx)
    rt <- system(plinkCmd, intern=FALSE)
    if (rt!=0) {
        stop("plink ped to bed error.")
    }
}
# that's it
