## For drawing different kinds of plot:
## 1. get_colVector()
## 2. draw_densityPlot()
## 3. draw_SigCounts()
## 4. draw_boxplot()
## 5. draw_heatmap()

library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(pheatmap)

# 1. Create a color vector containing 74 colour hexcodes

get_colVector <- function(n){
  col <- brewer.pal.info[brewer.pal.info$category=='qual',] # get max. 74 colours
  col_vector = unlist(mapply(brewer.pal, col$maxcolors, rownames(col)))
  if(n > length(col_vector)){
    vec <- sample(col_vector, n, replace=T)
  } else{
    vec <- sample(col_vector, n, replace=F)
  }
  
  vec
}

# 2. Density plot (x-ais: log-cpm, y-ais: density)

draw_densityPlot <- function(lcpmT, colVector, sampleNum){
  d <- plot(density(lcpmT[,1]), col=colVector[1], 
            lwd=2, ylim=c(0,0.3), las=2,main="",xlab="") +
    abline(v=0, lty=3) + title(xlab="Log-cpm") +
    for(i in 2:sampleNum){
      den <- density(lcpmT[,i])
      lines(den$x, den$y, col=dp_col[i], lwd=2)
    } 
  
  d
}

# 3. Summarize significant DEG numbers from the limma results, threshold: 1)q-val<0.05; 2)q-val<0.1; 3)q-val<0.1&FC > 1.5 

draw_SigCounts <- function(ttable){
  sig_table <- data.frame(
    thres = c("qval<0.05", "qval<0.1", "qval<0.1 & FC > 1.5"),
    counts = c(nrow(ttable[ttable$adj.P.Val < 0.05,]),
               nrow(ttable[ttable$adj.P.Val < 0.1,]),
               nrow(ttable[(ttable$adj.P.Val < 0.1 & abs(ttable$logFC) > log2(1.5)),])))
  
  p <- ggplot(data=sig_table, aes(x=thres, y=counts, fill=thres)) +
    geom_bar(stat="identity", width=0.5) + 
    geom_text(aes(label=counts), vjust=1.6, color="white", 
              position = position_dodge(0.9), size=3.5) +
    scale_fill_brewer(palette="Dark2",direction = -1) +
    guides(fill=FALSE) +
    labs(x="Threshold", y="DEG number") +
    theme_bw() 
  
  p
}

# 4. Boxplot for gene expression level comparison in overall and tissue-specific way.

draw_boxplot <- function(countsMat, tTable, phenoDF, geneid){
  exp_list <- as.data.frame(cpm(countsMat[geneid,], log=T))
  colnames(exp_list) <- c("counts")
  t_pdata <- phenoDF %>% select("sample_barcode", "type", "Metastasis")
  mdf <- merge(exp_list, t_pdata, by=0)
  symbol <- tTable[geneid,]$ID
  qVal <- formatC(tTable[geneid,]$adj.P.Val, format="e", digits=2)
  logFC <- round(tTable[geneid,]$logFC,3)
  
  p <- ggplot(data=mdf,aes(x=Metastasis,y=counts)) + 
    geom_boxplot(aes(fill=Metastasis))  + 
    scale_fill_manual(values = c("lightsteelblue2","wheat1")) + # box fill color
    labs(title=symbol, x="metastasis", y="log-CPM") + guides(fill=F) +
    theme_bw()
  
  b <- ggboxplot(mdf, x="Metastasis", y="counts",
                 color="Metastasis", palette="jco",main=paste0("Overall q-value=",qVal,", logFC=",logFC),
                 xlab="", ylab="log-cpm", add="jitter") +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          strip.background = element_rect(colour = "black", fill = "white")) +
    facet_wrap(~type, ncol=5)  + 
    stat_compare_means() # wilcoxon rank sum test
  
  ggarrange(p, b, ncol=2, nrow=1, widths=c(1,4))
}

# 5. Heat map for DEGs, using normalized counts in voom object, with the tumor stage and type annotations.

## Summarize the stage names first
stage_name <- function(df){
  colnames(df) <- c("type", "stage", "Metastasis")
  df$stage <- gsub("^Stage ", "", df$stage)
  df$stage <- gsub("[ABC]$", "", df$stage)
  df$stage[!(df$stage %in% c("I","II","III","IV"))] <- "NA"
  df
}

draw_heatmap <- function(voomObj,topTable,phenoDF,geneAnnoDF,esList){
  hm_cdr <- phenoDF %>% select(type, ajcc_pathologic_tumor_stage, Metastasis)
  hm_cdr <- stage_name(hm_cdr)
  hm_gene <- geneAnnoDF[esList,]$geneNames
  type <- get_colVector(length(unique(hm_cdr$type)))
  names(type) <- unique(hm_cdr$type)
  stage <- c("#d7191c","#fdae61","#a1d99b","#2b83ba","#bababa")
  names(stage) <- c("I","II","III","IV","NA")
  Metastasis <- c("#5e4fa2", "#fc8d59")
  names(Metastasis) <- unique(hm_cdr$Metastasis)
  anno_colors <- list(type = type,
                      stage = stage,
                      Metastasis = Metastasis) # the name must be consistent
  h <- pheatmap(voomObj$E[esList,],annotation_col=hm_cdr, annotation_colors=anno_colors, 
                labels_row = hm_gene, show_colnames = F)
  h
}





