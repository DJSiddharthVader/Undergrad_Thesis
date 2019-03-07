#!/usr/bin/env Rscript
#arg 1 is species tree in nexus format
#arg 2 is presenec absence matrix of genefamilies
#arg 3 is crispr annotation data

library("ape")
library("markophylo")
suppressMessages(library("jsonlite"))
suppressMessages(library("phangorn")) #map warnings

newfp <- function(fp){
    base <- basename(fp)
    newbase <- paste('rooted_',base,sep='')
    dir <- dirname(fp)
    newfp <- paste(dir,newbase,sep='/')
    return(newfp)
}
readSpeciesTree <- function(filepath){
    spt <- ape::read.nexus(file=filepath)
    rspt <- phangorn::midpoint(spt)
    rootfp <- newfp(filepath)
    ape::write.nexus(rspt,file=rootfp)
    rspt <- ape::read.nexus(file=rootfp)
    return(rspt)
}
trimVersion <- function(accnum){
    trimmed <- strsplit(accnum,'.',fixed=TRUE)[[1]][1]
    return(trimmed)
}
partitionByCRISPR <- function(spt,annodf){
    #annotation info
    tmp <- t(do.call(rbind,annodf))
    accession_col_idx <- grep("Accession Number", colnames(tmp))
    crispr_col_idx <- grep("isCRISPR", colnames(tmp))
    #tree nodes
    taxlist <- spt$tip.label
    root <- getRoot(spt)
    #create tip partitions
    crispr <- c()
    not_crispr <- c()
    for (acc in taxlist){
        rowidx <- match(trimVersion(acc),annodf[[accession_col_idx]])
        iscrispr <- annodf[[crispr_col_idx]][rowidx][[1]]
        nodenum <- match(acc,taxlist)
        if (iscrispr){
            crispr <- c(crispr,nodenum)
        } else {
            not_crispr <- c(not_crispr,nodenum)
        }
    }
    #convert tip partitions to brnachlists
    if (length(crispr) == 0 || length(not_crispr) == 0){
        allbranches <- 1:(length(spt$tip.label)+spt$Nnode)
        return(allbranches)
    } else {
        crispr <- which.edge(spt,crispr)
        not_crispr <- which.edge(spt,not_crispr)
        partitions <- list(crispr,not_crispr)
        return(partitions)
    }
}
main <- function(species_tree,pamat,crisprinfo){
    species_tree <- readSpeciesTree(species_tree)
    pamat <- read.csv(pamat,header=TRUE,row.names="gene_family")
    ajson <- fromJSON(crisprinfo)
    ajson['Cas Proteins (CRISPRone)'] <- NULL
    crispr_tip_partition <- partitionByCRISPR(species_tree,ajson)
    print(crispr_tip_partition)
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

if (sys.nframe() == 0){
    #CLI Args
    args <-  commandArgs(trailingOnly=TRUE)
    if (length(args) != 3){
        stop("args are species_tree.nexus, pa_matrix.csv, crispr_annotation_json .n",call.=FALSE)
    }
    main(args[1],args[2],args[3])
}

##DEPRECIATED
#jsonToList <- function(filepath){
#    json <- jsonlite::fromJSON(paste(readLines(filepath),collapse=""))
#    lst <- c(sapply(1:length(json)-1,function(x){json[[as.character(x)]]}))
#    return(lst)
#}
#jsonToDF <- function(filepath){
#    json <- jsonlite::fromJSON(paste(readLines(filepath),collapse=""))
#    df <- do.call("rbind",json)
#    df <- t(df)
#    return(df)
#}
