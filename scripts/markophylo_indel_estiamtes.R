#!/usr/bin/env Rscript
#arg 1 is species tree in nexus format
#arg 2 is presenec absence matrix of genefamilies
#arg 3 is crispr annotation data

library("ape")
library("markophylo")
#suppressMessages(library("jsonlite"))
suppressMessages(library("data.table"))
suppressMessages(library("phangorn")) #map warnings

newfp <- function(fp){
    dir <- dirname(fp)
    newfp <- paste(dir,paste('rooted_',basename(fp),sep=''),sep='/')
    newfpnwk <- paste(dir,paste('rooted_',basename(fp),'.newick',sep=''),sep='/')
    return(c(newfp,newfpnwk))
}
readSpeciesTree <- function(filepath){
    #read the species tree from MrBayes, root it and return the rooted tree
    #rooted tree is required by markophylo but the root is arbitrary
    rootedtree <- phangorn::midpoint(ape::read.nexus(file=filepath))
    if (!is.binary.tree(rootedtree)){
        rootedtree <- multi2di(rootedtree,random=FALSE)
    }
    rootedtreefile<- newfp(filepath)
    ape::write.nexus(rootedtree,file=rootedtreefile[1])
 #   ape::write.tree(rootedtree,file=rootfp[2])
    rootedtree <- ape::read.nexus(file=rootedtreefile[1])
    return(rootedtree)
}
loadAnnotationInfo <- function(filepath){
#    crisprAnnotation <- fromJSON(filepath)
    crisprAnnotation <- read.csv(filepath,sep='\t',stringsAsFactors=F)
    colnames(crisprAnnotation) <- gsub('.',' ',colnames(crisprAnnotation),fixed=T)
    crisprAnnotation['Cas Proteins (CRISPRone)'] <- NULL
    crisprAnnotation['isCRISPR'] <- as.logical(crisprAnnotation[['isCRISPR']])
    return(crisprAnnotation)
}
trimVersion <- function(accnum){
    trimmed <- strsplit(accnum,'.',fixed=TRUE)[[1]][1]
    return(trimmed)
}
partitionByCRISPR <- function(spt,crisprAnnotation){
    taxlist <- spt$tip.label
    root <- getRoot(spt)
    crisprAccessions <-  crisprAnnotation[(crisprAnnotation[['Accession Number']] %in% taxlist) &
                                          (crisprAnnotation[['isCRISPR']]),'Accession Number']
    crispr <- lapply(crisprAccessions,function(acc){nodepath(spt,from=match(acc,taxlist),to=root)[1:2]})
    non_crisprAccessions <- setdiff(taxlist,crisprAccessions)
    non_crispr <- lapply(non_crisprAccessions,function(acc){nodepath(spt,from=match(acc,taxlist),to=root)[1:2]})
    tipnodes <- 1:length(spt$tip.label)
    allnodes <- 1:(length(spt$tip.label)+spt$Nnode)
    internal <- allnodes[which(!allnodes %in% tipnodes)]
    partitions <- list(crispr=unlist(unique(crispr)),
                       non_crispr=unlist(unique(non_crispr)),
                       internal=internal)
    return(partitions)
}
plotTree <- function(speciesTree,partitions){
    with(partitions, {
        pdf(file='labelled_cladogram.pdf',width=25,height=15)
        plot(speciesTree,show.tip.label=FALSE,use.edge.length=TRUE)
        edgelabels(speciesTree$edge.length,frame='none',adj = c(0.5,-1.25))
        taxlist <- speciesTree$tip.label
        tipnums <-  1:length(taxlist)
        adj <- 0
        #label crispr tips
        if (length(crispr) != 0){
            tc = crispr[which(crispr %in% tipnums)]
            acclist = taxlist[tc]
            #acctext = paste('CRISPR',acclist,sep='   ')
            acctext = acclist
            tiplabels(acctext,tc,frame='rect',adj=adj,bg='cadetblue')
        }
        #label non-crispr tips
        if (length(non_crispr) != 0){
            tnc = non_crispr[which(non_crispr %in% tipnums)]
            acclist = taxlist[tnc]
            #acctext = paste('Non-CRISPR',acclist,sep='   ')
            acctext = acclist
            tiplabels(acctext,tnc,frame='rect',adj=adj,bg='firebrick1')
        }
        #label internal nodes
        internal <- (length(taxlist)+1):(speciesTree$Nnode+length(taxlist))
        #nodetext <- paste('Internal',internal,sep=' ')
        nodetext <- internal
        nodelabels(nodetext,internal,frame='circle',bg='green')
        dev.off()
    })
}
markophyloEstimate <- function(speciesTree,paMatrix,partitions){
    cleaned <- partitions[lapply(partitions,length)>0]
    estrates <- markophylo::estimaterates(usertree=speciesTree,
                                          userphyl=paMatrix,
                                          alphabet=c('D','B'),
                                          bgtype="listofnodes",
                                          bg=cleaned,
                                          rootprob="maxlik",
                                          modelmat="BDSYM",
                                          #numhessian=FALSE,
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
    args <-  commandArgs(trailingOnly=TRUE)
    if (length(args) != 3){
        stop("args are species tree, PA matrix, crispr annotation data .n",call.=FALSE)
    }
    speciesTreePath = args[1] #'./species_tree_files/species_tree/species_tree.con.tre'
    paMatrixPath = args[2] #'binary_pa_matrix.csv'
    crisprAnnotationPath = args[3] #'/home/sid/thesis_SidReed/data/fasta_crispr_annotation_df.json'
    main(speciesTreePath,paMatrixPath,crisprAnnotationPath)
}

#Estimated parameters on interval bounds
#this error indicates rate estimated just becomes 100.00, not sure why?

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
#partitionByCRISPR <- function(spt,crisprAnnotation){
#    #annotation info
#    tmp <- t(do.call(rbind,crisprAnnotation))
#    accession_col_idx <- grep("Accession Number", colnames(tmp))
#    crispr_col_idx <- grep("isCRISPR", colnames(tmp))
#    #tree nodes
#    taxlist <- spt$tip.label
#    root <- getRoot(spt)
#    #create tip partitions
#    crispr <- c()
#    not_crispr <- c()
#    for (acc in taxlist){
#        rowidx <- match(acc,crisprAnnotation[[accession_col_idx]])
#        iscrispr <- crisprAnnotation[[crispr_col_idx]][rowidx][[1]]
#        nodenum <- match(acc,taxlist)
#        pathToRoot <- nodepath(spt,from=nodenum,to=root)
#        if (iscrispr){
#            crispr <- c(crispr,pathToRoot)
#        } else {
#            not_crispr <- c(not_crispr,pathToRoot)
#        }
#    }
#    #check partition not empty and dedup
#    if (length(crispr) == 0 || length(not_crispr) == 0){
#        allnodes <- 1:(length(spt$tip.label)+spt$Nnode)
#        partitions <- list(non_crispr=allnodes) #all nodes in tree
#    } else {
#        partitions <- list(crispr=unique(crispr),non_crispr=unique(not_crispr))
#    }
#    return(partitions)
#}
#plotTree <- function(speciesTree,crisprAnnotation){
#    pdf(file='labelled_cladogram.pdf',width=25,height=15)
#    plot(speciesTree,show.tip.label=FALSE,use.edge.length=TRUE)
#    edgelabels(speciesTree$edge.length,frame='none')
#    #annotation info
#    tmp <- t(do.call(rbind,crisprAnnotation))
#    accession_col_idx <- grep("Accession Number", colnames(tmp))
#    crispr_col_idx <- grep("isCRISPR", colnames(tmp))
#    taxlist <- speciesTree$tip.label
#    #create tip partitions
#    crispr <- c()
#    not_crispr <- c()
#    for (acc in taxlist){
#        rowidx <- match(acc,crisprAnnotation[[accession_col_idx]])
#        iscrispr <- crisprAnnotation[[crispr_col_idx]][rowidx][[1]]
#        nodenum <- match(acc,taxlist)
#        if (iscrispr){
#            crispr <- c(crispr,nodenum)
#        } else {
#            not_crispr <- c(not_crispr,nodenum)
#        }
#    }
#    #label tips
#    if (length(crispr) != 0){
#        acclist = taxlist[crispr]
#        acctext = paste('CRISPR',acclist,sep='   ')
#        tiplabels(acctext,crispr,frame='rect', bg='blue')
#    }
#    if (length(not_crispr) != 0){
#        acclist = taxlist[not_crispr]
#        acctext = paste('Non-CRISPR',acclist,sep='   ')
#        tiplabels(acctext,not_crispr,frame='rect', bg='red')
#    }
#    internal <- (length(taxlist)+1):(speciesTree$Nnode+length(taxlist))
#    nodelabels(paste('Internal',internal,sep=' '),internal,
#               frame='circle', bg='green')
#    dev.off()
#}

