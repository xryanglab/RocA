setwd("../data")
## read data
RNA <- read.table("RNA.txt",sep="\t",stringsAsFactors = F,header = T,check.names = F)
rownames(RNA) <- RNA$gene_id
RNA <- RNA[,-1]

RPF <- read.table("RPF.txt",sep="\t",stringsAsFactors = F,header = T,check.names = F)
rownames(RPF) <- RPF$gene_id
RPF <- RPF[,-1]

human_info <- read.table("Homo_sapiens.GRCh38.88.transInfo.txt",sep="\t",stringsAsFactors = F,check.names = F,header = T)
longest_info <- read.table("longest.transcripts.info.txt",sep="\t",stringsAsFactors = F,check.names = F,header = T)

common_gene_id <- intersect(rownames(RPF),rownames(RNA))

RNA_common <- RNA[common_gene_id,]
RPF_common <- RPF[common_gene_id,]

common_gene_info <- human_info[which(human_info$gene_id%in%common_gene_id),]

## DE
library(DESeq2)
RocA03_RPF <- RPF[,c(1,2,4)]
indice <- Reduce(intersect,lapply(RocA03_RPF,function(x){which(x>=10)}))
RPF_filtered <- RocA03_RPF[indice,]
RPF_filtered_DESeq2 <- RPF_filtered
# construct sample information
colData <- data.frame(sample=colnames(RPF_filtered_DESeq2),condition=rep(c("DMSO","RocA03"),times=c(2,1)))
rownames(colData) <- as.character(colData$sample)

# construct DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = RPF_filtered_DESeq2,colData = colData,design = ~condition)

# normalization
dds <- DESeq(dds)

# get the results
res <- results(dds)
summary(res)
head(res)

# get differential expression genes
table(res$padj<0.05)
res <- res[order(res$padj),]

DE_genes_table <- subset(res,padj < 0.05 & abs(log2FoldChange)>=1)
write.table(DE_genes_table,"DE_genes_table.txt",sep="\t",quote = F,row.names = TRUE,col.names = TRUE)
DE_genes <- rownames(DE_genes_table)

# get expression matrix of DE genes(normalized) and counts
DE_genes_counts_normalized <- counts(dds,normalize=TRUE)[DE_genes,]
DE_genes_counts <- RPF_filtered_DESeq2[DE_genes,]
up_genes_res <- as.data.frame(subset(res,padj<=0.05&log2FoldChange>=1))
down_genes_res <- as.data.frame(subset(res,padj<=0.05&log2FoldChange<=-1))
Unchanged_genes_res <- as.data.frame(res)[setdiff(res@rownames,union(rownames(up_genes_res),rownames(down_genes_res))),]

up_gene_info <- longest_info[match(rownames(up_genes_res),longest_info$gene_id),]
down_gene_info <- longest_info[match(rownames(down_genes_res),longest_info$gene_id),]
unchanged_gene_info <- longest_info[match(rownames(Unchanged_genes_res),longest_info$gene_id),]

write.table(up_gene_info,"374_up_genes_info.txt",sep="\t",quote=F,row.names = F,col.names = T)
write.table(down_gene_info,"643_down_genes_info.txt",sep="\t",quote=F,row.names = F,col.names = T)
write.table(unchanged_gene_info,"3918_unchanged_genes_info.txt",sep="\t",quote=F,row.names = F,col.names = T)

write.table(up_genes_res,"374_up_genes_data.txt",sep="\t",quote=F,row.names = T,col.names = T)
write.table(down_genes_res,"643_down_genes_data.txt",sep="\t",quote=F,row.names = T,col.names = T)

## visualization
## volcano plot using ggplot2
log2FC <- res$log2FoldChange
padjust <- res$padj
data_for_plot <- data.frame(log2FC=log2FC,padjust=padjust)
data_for_plot <- na.omit(data_for_plot)
data_for_plot$sig[(data_for_plot$padjust>0.05|data_for_plot$padjust=="NA")|(data_for_plot$log2FC>-1|data_for_plot$log2FC<1)] <- "NO"
data_for_plot$sig[(data_for_plot$padjust<=0.05&data_for_plot$log2FC>=1)] <- "up"
data_for_plot$sig[(data_for_plot$padjust<=0.05 & data_for_plot$log2FC<=-1)] <- "down"
summary(data_for_plot$sig=="up")
summary(data_for_plot$sig=="down")
summary(data_for_plot$sig=="NO")
x_lim <- max(log2FC,-log2FC)
library(ggplot2)
library(RColorBrewer)
theme_set(theme_bw())

p <- ggplot(data_for_plot,aes(log2FC,-1*log10(padjust),color = sig))+geom_point()+xlim(-4,4)+labs(x="log2RPF(RocA03/DMSO)",y="-log10(p.adjust)")

p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))+ylim(0,15)
p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p