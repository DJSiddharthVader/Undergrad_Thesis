#!/usr/bin/env Rscript
#arg 1 is species tree in nexus format
#arg 2 is presenec absence matrix of genefamilies
#arg 3 is column index json
#arg 4 is row    index json
#arg 5 is crispr annotation json

library("ape")
library("rjson")
library("markophylo")
suppressMessages(library("phytools")) #map warnings

#CLI Args
args <-  commandArgs(trailingOnly=TRUE)
if (length(args) != 4){
    stop("args are species_tree.nexus, pa_matrix.txt, column index, row index .n",call.=FALSE)
}
jsonToList <- function(filepath){
    json <- rjson::fromJSON(paste(readLines(filepath),collapse=""))
    lst <- c(sapply(1:length(json)-1,function(x){json[[as.character(x)]]}))
    return(lst)
}
jsonToDF <- function(filepath){
    json <- rjson::fromJSON(paste(readLines(filepath),collapse=""))
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
    df <- rjson::fromJSON(paste(readLines(annotation_file),collapse=""))
    return(df)
}
partionByCRISPR <- function(taxlist,annodf){
    tmp <- t(do.call(rbind,annodf))
    accession_col_idx <- grep("Accession Number", colnames(tmp))
    crispr_col_idx <- grep("isCRISPR", colnames(tmp))
    partitions <- vector("list", 2)
    names(partitions) <- c("crispr","not_crispr")
    for (acc in taxlist){
        rowidx <- match(acc,annodf[[accession_col_idx]])
        iscrispr <- annodf[[crispr_col_idx]][rowidx][[1]]
        if (iscrispr){
            partitions$crispr <- c(partitions$crispr,acc)
        } else {
            partitions$not_crispr <- c(partitions$not_crispr,acc)
        }
    }
    print(partitions)
    return(partitions)
}
main <- function(args){
    #read species tree
    species_tree <- readSpeciesTree(args[1])
    #read P/A matrix, name cols/rows
    pamat <- read.csv(args[2],header=FALSE)
    pamat <- t(pamat)
    colnames(pamat) <- jsonToList(args[3])
    row.names(pamat) <- jsonToList(args[4])
    #create crispr tip partitions
    annodf <- readAnnoData()
    crispr_tip_partition <- partionByCRISPR(species_tree$tip.label,annodf)
    print(length(crispr_tip_partition$crispr))
    print(length(crispr_tip_partition$not_crispr))
    #markophylo
    estrates <- markophylo::estimaterates(usertree=species_tree,
                                         userphyl=pamat,
                                         alphabet=c(0,1),
                                         bgtype="listofnodes",
                                         bg=crispr_tip_partition,
                                         modelmat="BDSYM",
                                         rootprob="maxlik",
                                         matchtiptodata=TRUE)
    #print(estrates)
    #write markophylo results to csv
}
main(args)
