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

inpedf <- paste0(args[1], ".ped")
inmapf <- paste0(args[1], ".map")
if(!file.exists(inpedf)) {
    stop(paste0("Error: cannot find the input file \"", inpedf, "\".\n"))
}
if(!file.exists(inmapf)) {
    stop(paste0("Error: \"", inmapf, "\" not found. Although not required, it is a mistake to have a ped file without an accompanying map file.\n"))
}

outpf <- paste0(args[1], "_femcase.ped")
outmf <- paste0(args[1], "_femcase.map")
# 
awkCall <- paste0("awk '{$6=$5; print $0}' ", inpedf, " > ", outpf)
rt <- system(awkCall, intern=F)

file.copy(inmapf, outmf)
