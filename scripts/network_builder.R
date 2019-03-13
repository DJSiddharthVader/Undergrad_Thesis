library("ape")
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
#    networks_raw/all_gene_trees.txt
#                 subset_01.txt
#                 ...
#    networks_json/all_gene_trees.json
#                  subset_01.json
#                  ...

setUpDirs <- function(rootDir){
    netbase = paste(rootDir,'networkfiles',sep='/')
    mkdir(netbase)
    mkdir(paste(netbase,'all_newick_trees',sep='/'))
    mkdir(paste(netbase,'subsets',sep='/'))
}

convertToNwk <- function(pathToNexus,isSpecies=FALSE){
    #newbase = .../genusname/
    tree <- ape::read.nexus(file=pathToNexus)
    if (isSpecies){
        nwkpath <- 'networkfiles/all_newick_trees/sepcies.newick'
    } else{
        nwkpath <- paste('networkfiles/all_newick_trees/',strsplit(basename(pathToNexus),'.')[1],'.newicl',sep='')
    }
    ape::write.tree(tree,nwkpath)
}


main <- function(rootDir){
    setUpDirs(rootDir)
    speciesTreePath <- paste(rootDir,'species_tree_files/species_tree/rooted_species_tree.con.tre',sep=''
    convertToNwk(speciesTreePath,isSpecies=TRUE)
    geneTreeFiles <- list.files(paste(rootDir,'gene_tree_files/trees/',sep=''),pattern='**/*.con.tre',full.names=TRUE)
    print(geneTreeFiles)
}

#Markophylo Functions
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
        rspt <- multi2di(rspt)
    }
    rootfp <- newfp(filepath)
    ape::write.nexus(rspt,file=rootfp[1])
    ape::write.tree(rspt,file=rootfp[2])
    rspt <- ape::read.tree(file=rootfp[2])
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
        partitions <- list(non_crispr=allnodes) #all nodes in tree
    } else {
        partitions <- list(crispr=unique(crispr),non_crispr=unique(not_crispr))
    }
    return(partitions)
}
