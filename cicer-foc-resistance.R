library(tximport) # importing kalisto transcript counts to geneLevels
library(DESeq2) # rna-seq
library(tidyverse)
library(readr) # Fast readr of files.
library(GenomicFeatures)
library(EnhancedVolcano)
library(dplyr)
library(tidyr)
library(AnnotationHub)
library(clusterProfiler)

#SETUP
ah <- AnnotationHub()
orgdb <- query(ah, c("OrgDb", "maintainer@bioconductor.org"))[[1]]
keytypes(orgdb)
unique(ah$species)
query(ah, "Cicer")
Cicer <- ah[["AH114954"]]
keys(Cicer)
keytypes(Cicer)


metadata <- read.delim("~/samples.txt", header=T, row.names = 2)
metadata
coldata <- metadata
coldata$Treatment <- as.factor(coldata$Treatment)

TxDb <- makeTxDbFromGFF(file = "GCF_000331145.1_ASM33114v1_genomic.gff")
k <- keys(TxDb, keytype = "TXNAME")
tx2gene <- select(TxDb, k, "GENEID", "TXNAME")
tx2gene$TXNAME <- paste("rna-",tx2gene$TXNAME, sep = "")
head(tx2gene)
tx2gene$TXNAME <- str_split_i(tx2gene$TXNAME, pattern = "\\.", 1)

files <- paste("Kallisto_quant" ,  
               list.files(path = "Kallisto_quant", pattern = ".tsv", recursive = TRUE),
               sep = "/")
names(files) <- str_split_i(str_split_i(files, pattern = "/", 3), pattern = "\\.", 1)

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

txi.kallisto.tsv$counts <- round(txi.kallisto.tsv$counts[, rownames(coldata)],digits = 0)
all(rownames(coldata) == colnames(txi.kallisto.tsv$counts))

#ANALYSIS

coldata.sub <- coldata[which(coldata$Feature == "resistant"), ]

cts.sub <- txi.kallisto.tsv$counts[, which(colnames(txi.kallisto.tsv$counts) %in% rownames(coldata.sub))]

dds <- DESeqDataSetFromMatrix(countData = cts.sub,
                              colData = coldata.sub,
                              design = ~ Species + Treatment)
keep <- rowSums(counts(dds)) >= 400 
dds <- dds[keep,]
dds$Treatment <- factor(dds$Treatment, levels = c("control", "Fusarium_treated"))
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
summary(res)


sig.res <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 2.0 ),]

sig.res.desc <- merge(as.data.frame(sig.res), desc, by=0)

as.data.frame(sig.res.desc) %>% 
  write_delim(file="~/Resistant.DEGs.txt", delim = "\t", col_names = T)

goi <- rownames(sig.res)[1:4]
stopifnot(all(goi %in% names(dds)))
goi

#counts plot
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))


ggplot(tcounts, aes(Treatment, expression, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~gene, scales="free_y") +
  labs(x="Fusarium treatment",
       y="Expression (log normalized counts)",
       fill="Treatment",
       title="Top Results")


deg.counts <- log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))
pheatmap::pheatmap(deg.counts, annotation_col = coldata.sub, show_rownames = F)

original_gene_list <- sig.res$log2FoldChange

names(original_gene_list) <- rownames(sig.res)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order 
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800,
             nPermSimple = 10000,
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = Cicer, 
             pAdjustMethod = "none")
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
#emapplot(gse, showCategory = 10)
#cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

as.data.frame(gse) %>% 
  write_delim(file="~/Resistant.GOterms.txt", delim = "\t", col_names = T)

sensdegs <- read.delim("Sensitive.DEGs.txt", header = T, row.names = 1)
head(sensdegs)
resdegs <- read.delim("Resistant.DEGs.txt", header = T, row.names = 1)
head(resdegs)

library(VennDiagram)
venn.diagram(
  x = list(rownames(sensdegs), rownames(resdegs)),
  category.names = c("Sensitive DEGs" , "Resistant DEGs"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)

#for pca plot

table(vsd$Species)
table(vsd$Treatment)
table(vsd$Feature)
table(vsd$PlantId)

plotPCA(vsd, intgroup=c("Species")) +
  geom_point(aes(shape=factor(coldata$Treatment)), size=6) +
  geom_text(aes(label=coldata$PlantId , hjust = 1.2, nudge_x = 0.5))

plotPCA(vsd, intgroup=c("Feature")) +
  geom_point(aes(shape=factor(coldata$Treatment)), size=6) +
  geom_text(aes(label=coldata$PlantId , hjust = 1.2, nudge_x = 0.5))