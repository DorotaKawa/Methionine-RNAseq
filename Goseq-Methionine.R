library("goseq")
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

setwd("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/Methionine_experiment/DEG_v2/")
#background
DEList1 <- read.csv("Methionine_cpmExpression.csv", header = T,stringsAsFactors = F)
head(DEList1)
DEList2 <- DEList1
DEList <- list(DEList1, DEList2)
names(DEList) <- c("Met_up", "Met_down")

Met <- read.csv("DiffExprs/GeneLists/Met_effect.csv", header = T,stringsAsFactors = F)
Met_up <- subset(Met, logFC>=1) 
Met_down <- subset(Met, logFC<=-1)
uniqueList <- list(Met_up, Met_down)
names(uniqueList) <- c("Met_up", "Met_down")
sapply(uniqueList,nrow)

setwd("/Users/Dorota/Dropbox/UC_DAVIS/PROMISE/POP#2/RNA/DEG/4.SQR2/GOseq")
GFFfile = "Sbicolor_454_v3.1.1.gff.txt"

GFF <- import.gff(GFFfile,version="3",feature.type="gene")

grl <- reduce(split(GFF, mcols(GFF)$Name))
reducedGTF <- unlist(grl, use.names=T)
mcols(reducedGTF)$Name <- rep(names(grl), elementNROWS(grl))
reducedGTF
AllLengths <- width(reducedGTF)
names(AllLengths) <- mcols(reducedGTF)$Name
#
head(reducedGTF)
head(AllLengths)

go <- read.delim("Sbicolor_454_v3.1.1.annotation_go.txt",header = F)

colnames(go) <- c("AGI",
                  "GOID"
)

#head(go)
go.goseq <- go[,c("AGI", "GOID")]
head(go.goseq)

GOList <- list()
for (each in names(DEList)){
  print (each) 
  tmp <- DEList[[each]]
  ##
  rownames(tmp) <- tmp$ID
  assayed.genes <- rownames(tmp)
  de.genes <- uniqueList[[each]]$ID
  # de.genes <- rownames(tmp)[which(tmp[,grep("adj",colnames(tmp))] < 0.05)]#rownames(tmp)[tmp$adj.P.Val<0.05]
  length(assayed.genes)
  length(de.genes)
  #
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  head(gene.vector)
  table(gene.vector)
  
  gene.vector2 <- gene.vector
  #names(gene.vector2) <- gsub("..$","",names(gene.vector))
  
  GeneLengths <- AllLengths[names(gene.vector2)]
  ##
  cbind(GeneLengths,gene.vector)
  
  all.genes <- names(gene.vector)
  pwf=nullp(gene.vector, "AGI", id=all.genes, bias.data=GeneLengths)
  head(pwf)
  
  # 
  go.all <- goseq(pwf, "ITAG3.1", gene2cat=go.goseq,method = "Hypergeometric")
  head(go.all)
  dim(go.all)
  table(go.all$over_represented_pvalue < 0.05)
  
  go.sign <- go.all[go.all$over_represented_pvalue < 0.05,]
  
  
  #View(go.sign)
  GOList[[each]] <- go.all
  
}

sapply(GOList, dim)

signGOList <- lapply(GOList, function(x){  x[x$over_represented_pvalue < 0.05,] })

lapply(signGOList, head)

#lapply(signGOList, function(x){x[grep("nitr",x$term),]})
source("metaFunctions_forNetworkAnalysis.R")

list2env(signGOList,envir=.GlobalEnv)

setwd("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/Methionine_experiment/DEG_v2/GO enrichments/")

write.csv(Met_up, "GO_FC1_Met_up.csv")
write.csv(Met_down, "GO_FC1_Met_down.csv")


