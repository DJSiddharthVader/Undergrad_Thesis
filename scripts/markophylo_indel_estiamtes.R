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
    if (!is.binary.tree(rspt)){
        rspt <- multi2di(rspt,random=FALSE)
    }
    rootfp <- newfp(filepath)
    ape::write.nexus(rspt,file=rootfp[1])
 #   ape::write.tree(rspt,file=rootfp[2])
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
    non_crispr <- c()
    for (acc in taxlist){
        rowidx <- match(acc,crisprAnnotation[[accession_col_idx]])
        iscrispr <- crisprAnnotation[[crispr_col_idx]][rowidx][[1]]
        nodenum <- match(acc,taxlist)
        pathToRoot <- nodepath(spt,from=nodenum,to=root)
        tipbranch <- pathToRoot[1:2]
        if (iscrispr){
            crispr <- c(crispr,tipbranch)
        } else {
            non_crispr <- c(non_crispr,tipbranch)
        }
    }
    tipnodes <- 1:length(spt$tip.label)
    allnodes <- 1:(length(spt$tip.label)+spt$Nnode)
    internal <- allnodes[which(!allnodes %in% tipnodes)]
    partitions <- list(crispr=unique(crispr),
                       non_crispr=unique(non_crispr),
                       internal=internal)
    return(partitions)
}
plotTree <- function(speciesTree,partitions){
    with(partitions, {
        pdf(file='labelled_cladogram.pdf',width=25,height=15)
        plot(speciesTree,show.tip.label=FALSE,use.edge.length=TRUE)
        edgelabels(speciesTree$edge.length,frame='none',adj = c(0.5,-0.25))
        taxlist <- speciesTree$tip.label
        #label crispr tips
        if (length(crispr) != 0){
            tc = crispr[which(crispr %in% speciesTree$tip.label)]
            acclist = taxlist[tc]
            acctext = paste('CRISPR',acclist,sep='   ')
            tiplabels(acctext,tc,frame='rect', bg='cadetblue')
        }
        #label non-crispr tips
        if (length(non_crispr) != 0){
            tnc = non_crispr[which(non_crispr %in% speciesTree$tip.label)]
            acclist = taxlist[tnc]
            acctext = paste('Non-CRISPR',acclist,sep='   ')
            tiplabels(acctext,tnc,frame='rect', bg='firebrick1')
        }
        #label internal nodes
        internal <- (length(taxlist)+1):(speciesTree$Nnode+length(taxlist))
        nodelabels(paste('Internal',internal,sep=' '),internal,
                   frame='rect', bg='green')
        dev.off()
    })
}
markophyloEstimate <- function(speciesTree,paMatrix,partitions){
    estrates <- markophylo::estimaterates(usertree=speciesTree,
                                          userphyl=paMatrix,
                                          alphabet=c(0,1),
                                          bgtype="listofnodes",
                                          bg=partitions,
                                          rootprob="maxlik",
                                          modelmat="BDSYM",
                                          numhessian=FALSE,
                                          matchtiptodata=TRUE)
    return(estrates)
}
main <- function(speciesTreePath,paMatrixPath,crisprAnnotationPath){
    speciesTree <- readSpeciesTree(speciesTreePath)
    paMatrix <- read.csv(paMatrixPath,header=TRUE,row.names="gene_family")
    crisprAnnotation <- loadAnnotationInfo(crisprAnnotationPath)
    partitions <- partitionByCRISPR(speciesTree,crisprAnnotation)
    print('plotting tree')
    plotTree(speciesTree,partitions)
    print('running markophylo')
    rates <- markophyloEstimate(speciesTree,paMatrix,partitions)
    sink("markophylo_results.txt")
    print(rates)
    sink()
}

if (sys.nframe() == 0){
    #CLI Args
#    args <-  commandArgs(trailingOnly=TRUE)
#    if (length(args) != 3){
#        stop("args are species_tree.nexus, pa_matrix.csv, crispr_annotation_json .n",call.=FALSE)
#    }
#    main(args[1],args[2],args[3])
    speciesTreePath = './species_tree_files/species_tree/species_tree.con.tre'
    paMatrixPath = 'binary_pa_matrix.csv'
    crisprAnnotationPath = '/home/sid/thesis_SidReed/data/fasta_crispr_annotation_df.json'
    main(speciesTreePath,paMatrixPath,crisprAnnotationPath)
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
        rowidx <- match(acc,crisprAnnotation[[accession_col_idx]])
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
        partitions <- list(non_crispr=allnodes) #all nodes in tree
    } else {
        partitions <- list(crispr=unique(crispr),non_crispr=unique(not_crispr))
    }
    return(partitions)
}
plotTree <- function(speciesTree,crisprAnnotation){
    pdf(file='labelled_cladogram.pdf',width=25,height=15)
    plot(speciesTree,show.tip.label=FALSE,use.edge.length=TRUE)
    edgelabels(speciesTree$edge.length,frame='none')
    #annotation info
    tmp <- t(do.call(rbind,crisprAnnotation))
    accession_col_idx <- grep("Accession Number", colnames(tmp))
    crispr_col_idx <- grep("isCRISPR", colnames(tmp))
    taxlist <- speciesTree$tip.label
    #create tip partitions
    crispr <- c()
    not_crispr <- c()
    for (acc in taxlist){
        rowidx <- match(acc,crisprAnnotation[[accession_col_idx]])
        iscrispr <- crisprAnnotation[[crispr_col_idx]][rowidx][[1]]
        nodenum <- match(acc,taxlist)
        if (iscrispr){
            crispr <- c(crispr,nodenum)
        } else {
            not_crispr <- c(not_crispr,nodenum)
        }
    }
    #label tips
    if (length(crispr) != 0){
        acclist = taxlist[crispr]
        acctext = paste('CRISPR',acclist,sep='   ')
        tiplabels(acctext,crispr,frame='rect', bg='blue')
    }
    if (length(not_crispr) != 0){
        acclist = taxlist[not_crispr]
        acctext = paste('Non-CRISPR',acclist,sep='   ')
        tiplabels(acctext,not_crispr,frame='rect', bg='red')
    }
    internal <- (length(taxlist)+1):(speciesTree$Nnode+length(taxlist))
    nodelabels(paste('Internal',internal,sep=' '),internal,
               frame='circle', bg='green')
    dev.off()
}
