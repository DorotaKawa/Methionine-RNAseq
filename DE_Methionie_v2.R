#RNA seq Methionie experiment. Plants grown with or without methionie amendment, with or without Striga.
#2 week post infection, methionie added with Striga conditioning. Soil:sand mix
#Experiment done in November/December 2021 UCDavis

library(edgeR)
library(RColorBrewer)
library(limma)
library(Glimma)
library(reshape)
library(gplots)
library(calibrate)

setwd("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/Methionine_experiment/DEG_v2")
######################## Output folders
outDir = "DiffExprs/"
dir.create(outDir, showWarnings=T)

geneListsDir = "DiffExprs/GeneLists"
dir.create(geneListsDir, showWarnings=T)

######################## Counts file
counts <- read.csv("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/Methionine_experiment/Methionine_Counts.csv",
                   header = T,stringsAsFactors = F,row.names = 1)
dim(counts)

######################## Metadata table
# meta <- data.frame(do.call("rbind",strsplit(colnames(counts),"_")),
#                    row.names = colnames(counts),stringsAsFactors = F)
# colnames(meta) <- c("Infection","Amendment", "Rep")

meta<- read.csv("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/Methionine_experiment/Methionine_meta.csv",
            header = T,stringsAsFactors = F,row.names = 1)
head(meta)
####################### Removing lowly-expressed genes
# Choosing tresholg for CPMs
hist(log(as.matrix(counts)+1,10))
minCPM <- 10^-0.2 # count threshold
sampleMin <- 3 

dge <- DGEList(counts=counts,remove.zeros = T)
#Removing 6437 rows with all zero counts
isexpr <- rowSums(cpm(dge) > minCPM) >= sampleMin 
dge <- dge[isexpr,,keep.lib.size = FALSE]

#######################  Calculate the normalization factors with TMM
dge.norm <- calcNormFactors(dge, method = "TMM")
dge.norm$samples$norm.factors

par(mfrow=c(1,2))
lcpm <- cpm(dge, log=TRUE)
boxplot(lcpm)
title(main="A. Example: Unnormalised data",ylab="Log-cpm")

lcpm2 <- cpm(dge.norm, log=TRUE)
boxplot(lcpm2)
title(main="B. Example: Normalised data",ylab="Log-cpm")

##################### MDS plots
meta$Infection <- as.factor(meta$Infection)
meta$Amendment <- as.factor(meta$Amendment)
meta$Rep <- as.factor(meta$Rep)
glMDSPlot(dge, labels=rownames(dge.norm$samples), 
          groups=meta, launch=TRUE)

par(mfrow=c(2,2))
col.Infection <- meta$Infection
levels(col.Infection) <-  brewer.pal(nlevels(col.Infection), "Set1")
col.Infection <- as.character(col.Infection)
col.met <- meta$Amendment
levels(col.met) <-  brewer.pal(nlevels(col.met), "Set2")
col.met <- as.character(col.met)
col.rep <- meta$Rep
levels(col.rep) <-  brewer.pal(nlevels(col.rep), "Set2")
col.rep <- as.character(col.rep)

plotMDS(lcpm2, labels=meta$Infection,col=col.Infection)
title(main="A. Infection")
plotMDS(lcpm2, labels=meta$Amendment, col=col.met)
title(main="B. Amendment")
plotMDS(lcpm2, labels=meta$Rep, col=col.rep)
title(main="C. Replicate")
plotMDS(lcpm2, labels=rownames(dge.norm$samples))
title(main="D. Sample")

#### Design matrix
#define interaction effect
InfMet <- paste(meta$Infection,meta$Amendment, sep=".")
InfMet <- as.factor(InfMet)
design <- model.matrix(~0+InfMet, data=meta)
colnames(design) <- levels(InfMet)
head(design)

####################### Differential expression - transform with voom (and renormalize with quantile).
v <- voom(dge.norm, design, plot = TRUE,normalize.method = "quantile")

####################### Save Expression Data
cpmExpression <- cpm(dge.norm)
voomExpression <- v$E
colnames(cpmExpression) <- rownames(meta)
colnames(voomExpression) <-  rownames(meta)

save(cpmExpression,file = "Methionine_cpmExpression.RData")
save(voomExpression,file = "Methionine_voomExpression.RData")
write.table(cpmExpression, file = "Methionine_cpmExpression.csv",sep = ",",quote = F,col.names = T )

### Contrasts & DGE

cont.matrix <-makeContrasts(
  Met_effect = (Ctr.Met+Striga.Met)/2 - (Ctr.Ctr+Striga.Ctr)/2,
  Str_effect = (Striga.Ctr+Striga.Met)/2 - (Ctr.Ctr+Ctr.Met)/2,
  Str_effect_woMet = (Striga.Ctr-Ctr.Ctr),
  Str_effect_wMet = (Striga.Met-Ctr.Met),
  Met_effect_woStriga = (Ctr.Met-Ctr.Ctr),
  Met_effect_wStriga = (Striga.Met-Striga.Ctr),
  Interaction = (Striga.Met-Striga.Ctr)-(Ctr.Met-Ctr.Ctr),
  levels=design)

fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))
colSums(abs(decideTests(fit2)))

Contrasts <- colnames(cont.matrix)
Contrasts
Met_effect <- topTable(fit2, coef="Met_effect", adjust="BH",number = Inf,sort.by = "none")
Str_effect <- topTable(fit2, coef="Str_effect", adjust="BH",number = Inf,sort.by = "none")
Interaction <-  topTable(fit2, coef="Interaction", adjust="BH",number = Inf,sort.by = "none")

Str_effect_woMet<- topTable(fit2, coef="Str_effect_woMet", adjust="BH",number = Inf,sort.by = "none")
Str_effect_wMet <- topTable(fit2, coef="Str_effect_wMet", adjust="BH",number = Inf,sort.by = "none")
Met_effect_woStriga <-  topTable(fit2, coef="Met_effect_woStriga", adjust="BH",number = Inf,sort.by = "none")
Met_effect_wStriga <-  topTable(fit2, coef="Met_effect_wStriga", adjust="BH",number = Inf,sort.by = "none")

allDEG <- list(Met_effect,Str_effect,Str_effect_woMet,Str_effect_wMet, Met_effect_woStriga, Met_effect_wStriga,Interaction)
names(allDEG) <-Contrasts

########################## DEG filter Padj < 0.05
significantDE <- lapply(allDEG,function(x){ x[x[,grep("adj.P.Val",colnames(x))] < 0.05,] })
sapply(significantDE,nrow)

sapply(significantDE,nrow)
  shortName <- "Methionine"
  saveDElist <- paste0("DEList_",shortName,".RData")
  save(file = saveDElist,DEList)
  sapply(DEList,colnames)
  
  names(significantDE)
setwd("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/Methionine_experiment/DEG_v2/DiffExprs/GeneLists")  

  write.table( significantDE[["Met_effect"]],file = "Met_effect.csv",sep = ",",quote = F,col.names = NA)
  write.table( significantDE[["Str_effect"]],file = "Str_effect.csv",sep = ",",quote = F,col.names = NA)
  write.table( significantDE[["Interaction"]],file = "Interaction.csv",sep = ",",quote = F,col.names = NA)
  
  write.table( significantDE[["Str_effect_woMet"]],file = "Str_effect_woMet.csv",sep = ",",quote = F,col.names = NA)
  write.table( significantDE[["Str_effect_wMet"]],file = "Str_effect_wMet.csv",sep = ",",quote = F,col.names = NA)
  
  write.table( significantDE[["Met_effect_woStriga "]],file = "Met_effect_woStriga.csv",sep = ",",quote = F,col.names = NA)
  write.table( significantDE[["Met_effect_wStriga"]],file = "Met_effect_wStriga.csv",sep = ",",quote = F,col.names = NA)
  
  
  DESummary <- t(summary(decideTests(fit2)))[,-2]
  colnames(DESummary) = c("Downregulated","Upregulated")
  write.csv(x=DESummary,"DESummary.csv",quote = F,row.names = T)
  