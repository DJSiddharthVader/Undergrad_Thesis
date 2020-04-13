#!/usr/bin/env Rscript

suppressMessages(library("ape"))
suppressMessages(library('pryr')) #for partial
suppressMessages(library('doSNOW'))
suppressMessages(library('igraph')) #masking
suppressMessages(library("phangorn")) #map warnings
suppressMessages(library('progress'))

set.seed(9)
#Docs
#GETTING HIDE TO RUN
## Must use the /usr/bin/lua executable lua to run it
## run the script ~/thesis_SidReed/hide2/score.lua on a dir of gene trees and a species tree (1 dir, 1 argument)
## The dir with the hide files (~/thesis_SidReed/hide2/) must be in the LD_ LIBRARY_PATH variable or it wont work
### Give the error "cant find shared object file liblua5.1.so
##Command is
### exprt LD_LIBRARY_PATH=~/thesis_SidReed/hide2/:$LD_LIBRARY_PATH
### /usr/bin/lua ~/thesis_SidReed/hide2/score.lua /dir/with/trees

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
gtdir <- 'gene_tree_files/trees'             #famXXX/famXXX.con.tre

#Utilities
makepb <- function(total,msg=NA){
    pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :tick_rate/s :elapsedfull",
                           total=total,
                           show_after=0,
                           clear=FALSE)
    if (!(is.na(msg))){pb$message(msg)}
    pb$tick(0)
    return(pb)
}
removeFinishedEntries <- function(sourcedir,targetdir,sp,tp){
    # used to prevent redoing every step/file by only looking at the input files without an extant corresponding output
    #source is a list of files used as input
    #target is the list of files that would be the output
    #remove all members in target from source to avoid recomputing files that already were computed i.e exist
    targets <- list.files(targetdir,pattern=tp,full.names=TRUE)
    targets <- lapply(targets,function(tf){strsplit(basename(tf),'.',fixed=T)[[1]][1]})
    sources <- list.files(sourcedir,pattern=sp,full.names=TRUE)
    remaining <- lapply(sources,function(sf){
                            bsf <- strsplit(basename(sf),'.',fixed=T)[[1]][1]
                            if (!(bsf %in% targets)){
                                return(sf) #returns NULL if conditon not met for some reason
                            }})
    return(remaining[lengths(remaining) != 0]) #remove NULL elements
}
runParallel <- function(fnc,iterable,processes,msg,exports=c(),ret=FALSE,...){
    #NOTE: fnc must take an argument for a progress bar to call pb$tick  for progress reporting
    print(msg)
    if (length(iterable) == 0) {
        print('skipped, already computed')
        return(NA)
    }
    pb <- txtProgressBar(max=length(iterable),style=3)
    prog <- function(n) setTxtProgressBar(pb,n)
    #pb <- makepb(length(iterable))
    #prog <- function(n) pb$update(match(n,iterable)/length(iterable))
    registerDoSNOW(makeCluster(processes))
    result <- foreach(i=iterable,
                      .options.snow=list(progress=prog),
                      .export=exports) %dopar% { return(fnc(i,...)) }
    cat('\n')
    if (ret){
        return(result)
    } else {
        return(NULL)
    }
}
setUpDirs <- function(rootDir,speciesTreePath){
    if (length(grep('16S',speciesTreePath)) > 0){ #if WGS in stp
        netbase = paste(rootDir,'network_files_16S',sep='/')
    } else {
        netbase = paste(rootDir,'network_files_WGS',sep='/')
    }
    dir.create(netbase,showWarnings=FALSE)
    dir.create(paste(netbase,'all_newick_trees',sep='/'),showWarnings=FALSE)
    dir.create(paste(netbase,'subsets',sep='/'),showWarnings=FALSE)
    dir.create(paste(netbase,'raw_hide',sep='/'),showWarnings=FALSE)
    dir.create(paste(netbase,'csvs',sep='/'),showWarnings=FALSE)
    return(netbase)
}


#Convert MrBayes output to newick
convertToNwk <- function(treedir,netdir,isSpecies=FALSE){
    treefile <- list.files(treedir,pattern='*.con.tre',full.names=TRUE)[1]
    if (is.na(treefile)){
        return(NULL)
    }
    tree <- ape::read.nexus(file=treefile)
    if (length(grep('copy',tree$tip.label)) != 0){
        return(FALSE)
    }
    if (!(ape::is.rooted(tree))){ tree <- phangorn::midpoint(tree) }
    tree$edge.length <- NULL
    if (!(ape::is.binary.tree(tree))){ tree <- ape::multi2di(tree) }
    if (isSpecies){
        nwkpath <- file.path(netdir,'all_newick_trees','species.newick')
    } else{
        filename <- strsplit(basename(treefile),'.',fixed=TRUE)[[1]][1]
        nwkpath <- file.path(netdir,'all_newick_trees',paste(filename,'.newick',sep=''))
    }
    ape::write.tree(tree,nwkpath)
    return(NULL)
}
convertAllTrees <- function(rootDir,speciesTreePath,processes,netdir){
    #convert species tree to newick
    convertToNwk(file.path(rootDir,speciesTreePath),netdir,isSpecies=TRUE)
    if (!(file.exists(file.path(netdir,'all_newick_trees','species.newick')))){
        print('no species tree, stopping...')
        quit()
    }
    #convert gene tree to newick
    dirlist <- removeFinishedEntries(file.path(rootDir,gtdir),
                                     file.path(rootDir,netdir,'all_newick_trees'),
                                     'fam*','fam*')
    runParallel(convertToNwk, dirlist, processes, 'converting trees...',netdir=netdir)
}

#Create samples of gene trees
generateSample <- function(diridx,allGeneTrees,size,netdir){
    sdir <- paste(netdir,'/subsets/subset_',diridx,sep='')
    dir.create(sdir)
    file.copy(paste(netdir,'/all_newick_trees/species.newick',sep=''),
              paste(sdir,'species.newick',sep='/'))
    sample <- sample(allGeneTrees,size=size,replace=FALSE)
    for (file in sample){
        file.copy(file,paste(sdir,basename(file),sep='/'))
    }
}

whichTreesSampled <- function(netdir,bool=TRUE){
    #check if all gene trees that were computed are represented in everysubset used for a hide network
    subsets <- list.files(paste(netdir,'/subsets',sep=''),full.names=TRUE)
    allsubsets <- unique(sapply(subsets,function(set) list.files(set,pattern='*.newick')))
    alltrees <- list.files(paste(netdir,'/all_newick_trees',sep=''),pattern='*.newick') #all trees generated
    if (bool){
        return(setequal(alltrees,allsubsets)) #boolean value, tree if all trees used
    } else {
        return(setdiff(alltrees,allsubsets)) #return list of unused trees
    }
}

sampleUntillAllUsed <- function(size,netdir){
    stopsampling <- whichTreesSampled(netdir)
    if (stopsampling){#if all trees exists in at least 1 sample
        return(NA) #stop sampling
    }
    subsetidx <- length(list.files(paste(netdir,'/subsets',sep='')))
    allGeneTrees <- whichTreesSampled(bool=F)
    while (!stopsampling){#if not
        subsetidx <- subsetidx + 1
        generateSample(subsetidx,allGeneTrees,size) #make anew sample with unused trees
        stopsampling <- whichTreesSampled() #check if all samples are used, if so stop
    }
}

sampleGeneTreesForNetworks <- function(num,size,processes,netdir){
    allGeneTrees <- list.files(paste(netdir,'/all_newick_trees',sep=''),full.names=TRUE)
    allGeneTrees <- allGeneTrees[!grepl('species.newick',allGeneTrees)]
    subsets <- length(list.files(paste(netdir,'/subsets',sep='')))
    if (subsets >= num){
        print('generating samples...')
        print('skipped, already computed')
        return(NA)
    }
    runParallel(generateSample,
                (subsets:(subsets+num)),
                processes,
                'generating samples...',
                exports=c('generateSample'),
                allGeneTrees=allGeneTrees,size=size,netdir=netdir)
    sampleUntillAllUsed(size,netdir)
    #non-parallel code
#    iters <- subsets:(subsets+num)
#    pb <- makepb(length(iters),msg='generating samples...')
#    lapply(iters,function(i){
#                    generateSample(i,allGeneTrees,size,netdir)
#                    pb$tick()
#                 })
    #scope issues with parallel sampling, not really an issue bc its pretty fast anyways
    #gsfnc <- pryr::partial(generateSample,allGeneTrees=allGeneTrees,size=size,netdir=netdir)
}

#Run HiDe on each samples
runHide <- function(treedir,netdir){
    lua <- '/usr/bin/lua'
    hideexe <- '/home/sid/thesis_SidReed/hide_program/score.lua'
    edgelist <- system(paste(lua,hideexe,treedir),intern=TRUE,ignore.stderr=TRUE)
    outf <- file(paste(netdir,'/raw_hide/',basename(treedir),'.hide',sep=''))
    writeLines(edgelist,outf)
    close(outf)
}
runAllHide <- function(processes,netdir){
    #subsets <- list.files(network_files/subsets,full.names=TRUE)
    subsets <- removeFinishedEntries(paste(netdir,'/subsets',sep=''),
                                     paste(netdir,'/raw_hide',sep=''),
                                     'subset*','subset*')
    runParallel(runHide,subsets,processes,'running HiDe...',netdir=netdir)
}

#Convert HiDe to edge list csvs
parseEdge <- function(hedge,enorm){
    split <- strsplit(hedge,'\t')
    weight <- as.numeric(split[[1]][1])
    if (weight == 0){
       return(NA)
    }
    esplit <- strsplit(split[[1]][2],' ')
    source <- strsplit(esplit[[1]][1],'-')[[1]][2]
    sink <- strsplit(esplit[[1]][3],'-')[[1]][2]
    direction <- split[[1]][3]
    nweight <- weight/enorm
    edge <- list(source=source,
                 sink=sink,
                 raw_weight=weight,
                 weight=nweight,
                 direction=direction)
    return(edge)
}
convertHideToCSV <- function(hidefile,netdir){
    #find a description of the format .hide output here http://acgt.cs.tau.ac.il/hide/
    #converts a .hide file in to a .csv, where each row is a edge i.e.
    #the direction is a probability of the tranfer occuring from (u to v)/(v to u)
    #node u, node v, edge weight, normalized weight, estimated direction
    #NCXX.1, NCYY.1, 38.4432    , 0.872473         , (80%/20%)
    subset <- strsplit(basename(hidefile),'.',fixed=TRUE)[[1]][1][1]
    edge_norm <- length(list.files(paste(netdir,'/subsets/',subset,sep='')))-1
    hedges <- readLines(hidefile)
    parsed_edge_list <- lapply(hedges,
                                function(hedge){
                                edge <- parseEdge(hedge,edge_norm)
                                if (!(is.na(edge))){
                                    return(edge)
                                }})
    parsed_edge_list <- do.call('rbind',parsed_edge_list)
    base <- strsplit(basename(hidefile),'.',fixed=TRUE)[[1]][1]
    outf <- file.path(netdir,'/csvs',paste(base,'.csv',sep=''))
    write.table(parsed_edge_list,file=outf,sep='~',row.names=F)
}
parseAllEdgeLists <- function(processes,netdir){
    hidefiles  <- removeFinishedEntries(paste(netdir,'/raw_hide',sep=''),
                                        paste(netdir,'/csvs',sep=''),
                                        'subset*','subset*')
    runParallel(convertHideToCSV,
                hidefiles,
                processes,
                'formatting networks...',
                exports=c('parseEdge'),
                netdir=netdir)
}

main <- function(rootDir,speciesTreePath,num,size,processes){
    netdir <- setUpDirs(rootDir,speciesTreePath)
    convertAllTrees(rootDir,speciesTreePath,processes,netdir)
    sampleGeneTreesForNetworks(num,size,processes,netdir)
    runAllHide(processes,netdir)
    parseAllEdgeLists(processes,netdir)
}

#if name == main
if (sys.nframe() == 0){
    #CLI Args
    args <-  commandArgs(trailingOnly=TRUE)
    species_tree <- ifelse(length(args) > 0,args[1],'species_tree_WGS') #use species_tree_16S/species_tree for the 16S built species tree
    genus_dir = ifelse(length(args) > 1,args[2],getwd())
    subset_number <- ifelse(length(args) > 2,as.numeric(args[3]),1000)
    subset_size <- ifelse(length(args) > 3,as.numeric(args[4]),50)
    processes <- ifelse(length(args) > 4,as.numeric(args[5]),8)
    main(genus_dir,species_tree,subset_number,subset_size,processes)
}

#DEPRECIATED
#convertAllTrees <- function(rootdir,speciesTreePath){
#    #convert species tree to newick
#    speciesTreePath <- file.path(rootDir,speciesTree)
#    speciesTreePath <- convertToNwk(speciesTreePath,isSpecies=TRUE)
#    #convert gene tree to newick
#    dirlist <- list.files(file.path(rootDir,gtdir),pattern='fam*',full.names=TRUE)
#    pb <- makepb(length(dirlist))
#    geneTreePaths <- c()
#    for (dir in dirlist){
#        tree <- list.files(dir,pattern='*.con.tre',full.names=TRUE)
#        gpath <- convertToNwk(tree[1])
#        if (typeof(gpath) == 'character'){
#            geneTreePaths <- c(geneTreePaths,gpath)
#        }
#        pb$tick()
#    }
#    return(list(sp=speciesTreePath,gt=geneTreePaths))
#}
