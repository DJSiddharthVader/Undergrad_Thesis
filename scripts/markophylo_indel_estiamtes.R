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
    newbasenwk <- paste('rooted_',base,'.newick',sep='')
    dir <- dirname(fp)
    newfp <- paste(dir,newbase,sep='/')
    newfpnwk <- paste(dir,newbasenwk,sep='/')
    return(c(newfp,newfpnwk))
}
readSpeciesTree <- function(filepath){
    spt <- ape::read.nexus(file=filepath)
    rspt <- phangorn::midpoint(spt)
    rootfp <- newfp(filepath)
    ape::write.nexus(rspt,file=rootfp[1])
    ape::write.tree(rspt,file=rootfp[2])
    rspt <- ape::read.nexus(file=rootfp[1])
    return(rspt)
}
loadAnnotationInfo <- function(filepath){
    crisprAnnotation <- fromJSON(filepath)
    crisprAnnotation['Cas Proteins (CRISPRone)'] <- NULL
    return(crisprAnnotation)
}
trimVersion <- function(accnum){
    trimmed <- strsplit(accnum,'.',fixed=TRUE)[[1]][1]
    return(trimmed)
}
partitionByCRISPR <- function(spt,crisprAnnotation){
    #annotation info
    tmp <- t(do.call(rbind,crisprAnnotation))
    accession_col_idx <- grep("Accession Number", colnames(tmp))
    crispr_col_idx <- grep("isCRISPR", colnames(tmp))
    #tree nodes
    taxlist <- spt$tip.label
    root <- getRoot(spt)
    #create tip partitions
    crispr <- c()
    not_crispr <- c()
    for (acc in taxlist){
        rowidx <- match(trimVersion(acc),crisprAnnotation[[accession_col_idx]])
        iscrispr <- crisprAnnotation[[crispr_col_idx]][rowidx][[1]]
        nodenum <- match(acc,taxlist)
        pathToRoot <- nodepath(spt,from=nodenum,to=root)
        if (iscrispr){
            crispr <- c(crispr,pathToRoot)
        } else {
            not_crispr <- c(not_crispr,pathToRoot)
        }
    }
    #check partition not empty and dedup
    if (length(crispr) == 0 || length(not_crispr) == 0){
        allnodes <- 1:(length(spt$tip.label)+spt$Nnode)
        partition <- list(allnodes) #all nodes in tree
        return(partitions)
    } else {
        return(list(unique(crispr),unique(not_crispr)))
    }
}
markophyloEstimate <- function(speciesTree,paMatrix,crisprAnnotation){
    partition <- partitionByCRISPR(speciesTree,crisprAnnotation)
    estrates <- markophylo::estimaterates(usertree=speciesTree,
                                         userphyl=paMatrix,
                                         alphabet=c(0,1),
                                         bgtype="listofnodes",
                                         bg=partition,
                                         rootprob="maxlik",
                                         modelmat="BDSYM",
                                         matchtiptodata=TRUE)
    #write to file
    sink("markophylo_results.txt")
    print(estrates)
    sink()
    return(estrates)
}
plotTree <- function(speciesTree,crisprAnnotation){
    pdf(file='labelled_cladogram.pdf',width=25,height=15)
    plot(speciesTree,show.tip.label=FALSE,use.edge.length=TRUE)
#    edf <- as.data.frame(tree$edge)
#    edgelabels(apply(edf[,colnames(edf)],1,paste,collapse='\n'))
    #annotation info
    tmp <- t(do.call(rbind,crisprAnnotation))
    accession_col_idx <- grep("Accession Number", colnames(tmp))
    crispr_col_idx <- grep("isCRISPR", colnames(tmp))
    taxlist <- speciesTree$tip.label
    #create tip partitions
    crispr <- c()
    not_crispr <- c()
    for (acc in taxlist){
        rowidx <- match(trimVersion(acc),crisprAnnotation[[accession_col_idx]])
        iscrispr <- crisprAnnotation[[crispr_col_idx]][rowidx][[1]]
        nodenum <- match(acc,taxlist)
        if (iscrispr){
            crispr <- c(crispr,nodenum)
        } else {
            not_crispr <- c(not_crispr,nodenum)
        }
    }
    tiplabels(paste('CRISPR',crispr,sep=' '),crispr)
    tiplabels(paste('Non-CRISPR',not_crispr,sep=' '),not_crispr)
    nodelabels()
    dev.off()
}
main <- function(speciesTreePath,paMatrixPath,crisprAnnotationPath){
    speciesTree <- readSpeciesTree(speciesTreePath)
    paMatrix <- read.csv(paMatrixPath,header=TRUE,row.names="gene_family")
    crisprAnnotation <- loadAnnotationInfo(crisprAnnotationPath)
    markophyloEstimate(speciesTree,paMatrix,crisprAnnotation)
    plotTree(speciesTree,crisprAnnotation)
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
#partitionByCRISPR <- function(spt,crisprAnnotation){
#    #annotation info
#    tmp <- t(do.call(rbind,crisprAnnotation))
#    accession_col_idx <- grep("Accession Number", colnames(tmp))
#    crispr_col_idx <- grep("isCRISPR", colnames(tmp))
#    #tree nodes
#    taxlist <- spt$tip.label
#    root <- getRoot(spt)
#    #create tip partitions
#    crispr <- c(root)
#    not_crispr <- c(root)
#    for (acc in taxlist){
#        rowidx <- match(trimVersion(acc),crisprAnnotation[[accession_col_idx]])
#        iscrispr <- crisprAnnotation[[crispr_col_idx]][rowidx][[1]]
#        nodenum <- match(acc,taxlist)
#        if (iscrispr){
#            crispr <- c(crispr,nodenum)
#        } else {
#            not_crispr <- c(not_crispr,nodenum)
#        }
#    }
#    #convert tip partitions to branchlists
#    if (length(crispr) == 0 || length(not_crispr) == 0){
#        allbranches <- 1:(length(spt$tip.label)+spt$Nnode)
#        return(allbranches)
#    } else {
#        crispr <- which.edge(spt,crispr)
#        not_crispr <- which.edge(spt,not_crispr)
#        partitions <- list(crispr,not_crispr)
#        return(partitions)
#    }
#}
