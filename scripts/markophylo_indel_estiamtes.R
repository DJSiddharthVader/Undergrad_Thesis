#!/usr/bin/env Rscript
#arg 1 is species tree in nexus format
#arg 2 is presenec absence matrix of genefamilies
#arg 3 is column index json
#arg 4 is row    index json

library("ape")
library("markophylo")
suppressMessages(library("jsonlite"))
suppressMessages(library("phytools")) #map warnings

#CLI Args
args <-  commandArgs(trailingOnly=TRUE)
if (length(args) != 4){
    stop("args are species_tree.nexus, pa_matrix.csv, column index, row index .n",call.=FALSE)
}
jsonToList <- function(filepath){
    json <- jsonlite::fromJSON(paste(readLines(filepath),collapse=""))
    lst <- c(sapply(1:length(json)-1,function(x){json[[as.character(x)]]}))
    return(lst)
}
jsonToDF <- function(filepath){
    json <- jsonlite::fromJSON(paste(readLines(filepath),collapse=""))
    df <- do.call("rbind",json)
    df <- t(df)
    return(df)
}
readSpeciesTree <- function(filepath){
    spectree <- ape::read.nexus(file=filepath)
    if (!(ape::is.rooted(spectree))){
        phytools::midpoint.root(spectree)
    }
    return(spectree)
}
readAnnoData <- function(){
    annotation_file <- "/home/sid/thesis_SidReed/pop_annotation_data_frame.json"
    df <- jsonlite::fromJSON(paste(readLines(annotation_file),collapse=""))
    return(df)
}
partionByCRISPR <- function(spectree,annodf){
    tmp <- t(do.call(rbind,annodf))
    taxlist <- spectree$tip.label
    accession_col_idx <- grep("Accession Number", colnames(tmp))
    crispr_col_idx <- grep("isCRISPR", colnames(tmp))
#    partitions <- vector("list", 2)
#    names(partitions) <- c("crispr","not_crispr")
    crispr <- vector()
    not_crispr <- vector()
    for (acc in taxlist){
        rowidx <- match(acc,annodf[[accession_col_idx]])
        iscrispr <- annodf[[crispr_col_idx]][rowidx][[1]]
        if (iscrispr){
            crispr <- c(crispr,acc)
        } else {
            not_crispr <- c(not_crispr,acc)
        }
    }
    crispr <- c(crispr,6:length(taxlist))
    partitions <- list(crispr,not_crispr)
    if (length(partitions[[1]]) == 0){
        return(list(1:length(taxlist),ancestral))
    } else if (length(partitions[[2]]) == 0){
        return(list(1:length(taxlist),ancestral))
    } else {
        return(partitions)
    }
}
main <- function(args){
    #read species tree
    species_tree <- readSpeciesTree(args[1])
    #read P/A matrix, name cols/rows
    pamat <- read.csv(args[2],header=FALSE)
    pamat <- t(pamat)
    row.names(pamat) <- jsonToList(args[3])
    colnames(pamat) <- jsonToList(args[4])
    #create crispr tip partitions
    annodf <- readAnnoData()
    crispr_tip_partition <- partionByCRISPR(species_tree,annodf)
    #markophylo
    estrates <- markophylo::estimaterates(usertree=species_tree,
                                         userphyl=pamat,
                                         alphabet=c(0,1),
                                         bgtype="listofnodes",
                                         bg=crispr_tip_partition,
                                         rootprob="maxlik",
                                         modelmat="BDSYM",
                                         matchtiptodata=TRUE)
    #write to file
    sink("markophylo_results.txt")
    print(estrates)
    sink()
}
main(args)
