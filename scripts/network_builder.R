#!/usr/bin/env Rscript
library("ape")
library('progress')
suppressMessages(library('igraph')) #masking
suppressMessages(library("phangorn")) #map warnings

set.seed(9)
#Docs
#STEPS
#1. convert species nexus_tree to newick tree
#2. convert gene nexus_trees to newick trees
#3. split networks into subsets (how yet not sure)
#4. run HiDe on each subset to produce 1 network
#5. parse raw HiDe output to some computer readable format
#6. write to json file for use in networkx/igraph
#produces file direcories as follows
#network_files/
#    all_newick_trees/species.newick
#                     fam1.newick
#                     ...
#    subsets/
#        subset_01/species.newick
#                  fam1.newick
#                  fam2.newick
#                  fam4.newick
#                  ...
#        subset_02/species.newick
#                  fam3.newick
#                  fam2.newick
#                  fam6.newick
#                  ...
#        ...
#    raw_hide/subset_01.hide
#             subset_02.hide
#             subset_03.hide
#             ...
#    ncols/subset_1.ncol
#          subset_2.ncol
#          subset_3.ncol
#          subset_4.ncol
#          ...

makepb <- function(total){
    pb <- progress_bar$new(format = "[:bar] :percent eta: :eta",total = total, clear = FALSE, width= 75)
    return(pb)
}
setUpDirs <- function(rootDir){
    netbase = paste(rootDir,'network_files',sep='/')
    dir.create(netbase)
    dir.create(paste(netbase,'all_newick_trees',sep='/'))
    dir.create(paste(netbase,'subsets',sep='/'))
    dir.create(paste(netbase,'raw_hide',sep='/'))
    dir.create(paste(netbase,'csvs',sep='/'))
}

convertToNwk <- function(pathToNexus,isSpecies=FALSE){
    tree <- ape::read.nexus(file=pathToNexus)
    if (length(grep('copy',tree$tip.label)) != 0){
        return(FALSE)
    }
    if (!is.rooted(tree)){
        tree <- phangorn::midpoint(tree)
    }
    tree$edge.length <- NULL
    if (!is.binary.tree(tree)){
        tree <- multi2di(tree)
    }
    if (isSpecies){
        nwkpath <- 'network_files/all_newick_trees/species.newick'
    } else{
        filename <- strsplit(basename(pathToNexus),'.',fixed=TRUE)[[1]][1]
        nwkpath <- paste('network_files/all_newick_trees/',filename,'.newick',sep='')
    }
    ape::write.tree(tree,nwkpath)
    return(nwkpath)
}

generateSample <- function(diridx,allGeneTrees,size){
    sdir <- paste('network_files/subsets/subset_',diridx,sep='')
    dir.create(sdir)
    file.copy('network_files/all_newick_trees/species.newick',paste(sdir,'species.newick',sep='/'))
    sample <- sample(allGeneTrees,size=size,replace=FALSE)
    for (file in sample){
        file.copy(file,paste(sdir,basename(file),sep='/'))
    }
}
sampleSet <- function(){
    allsubsets <- c()
    subsets <- list.files('network_files/subsets',full.names=TRUE)
    for (set in subsets){
        allsubsets <- c(allsubsets,list.files(set,pattern='*.newick'))
    }
    allsubsets <- unique(allsubsets)
    print(length(allsubsets))
    alltrees <- list.files('network_files/all_newick_trees',pattern='*.newick')
    print(length(alltrees))
    return(setequal(alltrees,allsubsets))
}
sampleUntillAllUsed <- function(size){
    stopsampling <- sampleSet()
    print(stopsampling)
    if (stopsampling){
        return(NULL)
    }
    subsetidx <- length(list.files('network_files/subsets'))+1
    allGeneTrees <- list.files('network_files/all_newick_trees',full.names=TRUE)
    allGeneTrees <- allGeneTrees[!grepl('species.newick',allGeneTrees)]
    while (!stopsampling){
        generateSample(subsetidx,allGeneTrees,size)
        stopsampling <- sampleSet()
        subsetidx <- subsetidxs + 1
    }
}
runHide <- function(treedir){
    lua <- '/usr/bin/lua'
    hideexe <- '/home/sid/thesis_SidReed/hide_program/score.lua'
    edgelist <- system(paste(lua,hideexe,treedir),intern=TRUE,ignore.stderr=TRUE)
    outf <- file(paste('network_files/raw_hide/',basename(treedir),'.hide',sep=''))
    writeLines(edgelist,outf)
    close(outf)
}
normalizeEdge <- function(weight,fname){
    subset <- strsplit(basename(fname),'.',fixed=TRUE)[[1]][1][1]
    subdir <- paste('network_files/subsets/',subset,sep='')
    norm <- length(list.files(subdir))-1
    return(weight/norm)
}
parseEdge <- function(line,hidefile){
    split <- strsplit(line,'\t')
    weight <- as.numeric(split[[1]][1])
    if (weight == 0){
       return(c(-1,-1,-1))
    }
    esplit <- strsplit(split[[1]][2],' ')
    source <- strsplit(esplit[[1]][1],'-')[[1]][2]
    sink <- strsplit(esplit[[1]][3],'-')[[1]][2]
    nweight <- normalizeEdge(weight,hidefile)
    direction <- split[[1]][3]
    edge <- list(source=source,sink=sink,raw_weight=weight,weight=nweight,direction=direction)
    return(edge)
}
hideToEdgeList <- function(hidefile){
    lines <- readLines(hidefile)
    parsed_edge_list <- parseEdge(lines[1],hidefile)
    for (line in lines[2:length(lines)]){
        edge <- parseEdge(line,hidefile)
        if (edge[3][1] == -1){
            next
        }
        parsed_edge_list <- rbind(parsed_edge_list,edge)
    }
    return(parsed_edge_list)
}
main <- function(rootDir,num,size){
    setUpDirs(rootDir)
    print('converting trees to newick')
    #convert species tree to newick
    speciesTreePath <- paste(rootDir,'species_tree_files/species_tree/species_tree.con.tre',sep='')
    speciesTreePath <- convertToNwk(speciesTreePath,isSpecies=TRUE)
    #convert gene tree to newick
    dirlist <- list.files(paste(rootDir,'gene_tree_files/trees',sep=''),pattern='fam*',full.names=TRUE)
    geneTreePaths <- c()
    pb <- makepb(length(dirlist))
    for (dir in dirlist){
        tree <- list.files(dir,patter='*.con.tre',full.names=TRUE)
        gpath <- convertToNwk(tree[1])
        if (typeof(gpath) == 'character'){
            geneTreePaths <- c(geneTreePaths,gpath)
        }
        pb$tick()
    }
    print('generating samples')
    #generate samples of gene trees
    allGeneTrees <- list.files('network_files/all_newick_trees',full.names=TRUE)
    allGeneTrees <- allGeneTrees[!grepl('species.newick',allGeneTrees)]
    pb <- makepb(num)
    for (s in 1:num){
        generateSample(s,allGeneTrees,size)
        pb$tick()
    }
    #sampleUntillAllUsed(size)
    print('running HiDe')
    #build networks
    subsets <- list.files('network_files/subsets',full.names=TRUE)
    pb <- makepb(length(subsets))
    for (sample in subsets){
        runHide(sample)
        pb$tick()
    }
    print('parsing and writing networks')
    #parse and write networks
    hidefiles <- list.files('network_files/raw_hide',full.names=TRUE)
    pb <- makepb(length(hidefiles))
    for (hide in hidefiles){
        edgelist <- hideToEdgeList(hide)
        base <- strsplit(basename(hide),'.',fixed=TRUE)[[1]][1]
        outf <- paste('network_files/csvs/',base,'.csv',sep='')
        write.table(edgelist,file=outf,sep='~',row.names=1:dim(edgelist)[1])
        pb$tick()
    }
}


#if name == main
if (sys.nframe() == 0){
    #CLI Args
    args <-  commandArgs(trailingOnly=TRUE)
    if (length(args) != 1){
        stop("arg is genus dir .n",call.=FALSE)
    }
    main(args[1],1000,50)
}

