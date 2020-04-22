#!/usr/bin/Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(grid))
#suppressPackageStartupMessages(library(gridExtra))

network_data_file <- './all_WGS_results.csv'
genus_data_file <- './all_data.tsv'

labelColumns <- function(stat,patterns){
    for (i in patterns){
        if (str_detect(as.character(stat),i)){
            return(i)
        }
    }
    return(NA)
}
splitStatColumn <- function(stat){
    #cnverts stat.crispr_mean to tuple(stat,crispr,mean), for later filtering for ggplot
    s1 <- strsplit(as.character(stat),'.',fixed=TRUE)[[1]]
    stat <- s1[1]
    s2 <- strsplit(paste(s1[2:length(s1)],collapse='.'),'_',fixed=TRUE)[[1]]
    iscrispr <- ifelse(s2[1] %in% c('crispr','non.crispr'),s2[1],NA)
    stype <- s2[2]
    return(c(stat,iscrispr,stype))
}
calculateCrisprFaction <- function(Genus,df){
    crispr.status <- as.logical(df[df$genus == Genus,]$isCRISPR)
    return(sum(crispr.status)/length(crispr.status))
}
loadData <- function(){
    ndf <- as.data.frame(fread(network_data_file),stringsAsFactors=FALSE)
    adf <- as.data.frame(fread(genus_data_file),stringsAsFactors=FALSE)
    df <- merge(ndf,adf,by.x='genus',by.y='Genus',stringsAsFactors=FALSE)
    gdf <- df[,c('genus',colnames(ndf))] %>% distinct(genus,.keep_all=TRUE)
    gdf <- gdf[,which(colnames(gdf)!='genus.1')]
    gdf <- melt(gdf,id.vars='genus',stringsAsFactors=FALSE)
    gdf$value <- as.numeric(gdf$value)
    gdf <- gdf[(!(is.na(gdf$value))),]
    gdf[,c('Stat','CRISPR.Status','StatType')] <- t(rbind(sapply(gdf$variable,splitStatColumn)))
    gdf$CRISPR.Fraction <- sapply(gdf$genus,calculateCrisprFaction,df=df)
    return(gdf)
}

wilcoxtest <- function(stat,df){
    data <- dcast(gdf[(gdf$Stat==stat)&(!(is.na(gdf$CRISPR.Status))),],
                  genus ~ CRISPR.Status + StatType)
    print(dim(data))
    wiltest <- wilcox.test(data$crispr <- mean,data$non.crispr <- mean,paired=TRUE)
    return(wiltest)
}
plotPairedMetric <- function(stat,gdf){
    fdata <- dcast(gdf[(gdf$Stat==stat)&(!(is.na(gdf$CRISPR.Status))),],genus ~ CRISPR.Status + StatType)
    fdata <- fdata[(!(is.na(fdata$non.crispr_mean)))&(!(is.na(fdata$crispr_mean))),]
    print(dim(fdata))
    fit <- lm(non.crispr_mean ~ crispr_mean, data=fdata)
    wiltest <- wilcox.test(fdata$crispr_mean,fdata$non.crispr_mean,paired=TRUE)
    statslabel <- paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                        "\nSlope =",signif(fit$coef[[2]], 5),
                        "\nP value =",signif(summary(fit)$coef[2,4], 5),
                        "\nWilcoxon Rank P value=",signif(wiltest$p.value,4))
    plot <- ggplot(fdata,aes(x=non.crispr_mean,y=crispr_mean)) +
        geom_point(size=1,color='black') +
        geom_smooth(method='lm',col='blue') +
        xlim(min(fdata$non.crispr_mean),max(fdata$non.crispr_mean)) +
        ylim(min(fdata$crispr_mean),max(fdata$crispr_mean)) +
        annotate('text',label=statslabel,color='black',size=2.5,
                 x=max(fdata$non.crispr_mean)*3/4,y=max(fdata$crispr_mean)*3/4) +
        xlab('Non-CRISPR') + ylab('CRISPR') +
        #ggtitle(paste(stat,' for various genera')) +
        theme(legend.position='top',panel.background=element_rect(fill='white'))
    return(plot)
}

