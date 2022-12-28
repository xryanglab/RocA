setwd("../data")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(org.Hs.eg.db)
GO_analysis <- read.table("GO_analysis.txt",sep="\t",stringsAsFactors = F,header = T,fill = T,check.names = F)
RocA_insensitive <- na.omit(GO_analysis$`RocA-up-regulated`)
RocA_sensitive <- na.omit(GO_analysis$`RocA-down-regulated`)
## construct dataframe
# up
RocA_insensitive_entrezID <- bitr(RocA_insensitive,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
group_up <- rep("RocA-up-regulated",times=length(RocA_insensitive_entrezID$ENTREZID))
data_up <- data.frame(geneID=RocA_insensitive_entrezID$ENTREZID,group=group_up)

# down
RocA_sensitive_entrezID <- bitr(RocA_sensitive,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
group_down <- rep("RocA-down-regulated",times=length(RocA_sensitive_entrezID$ENTREZID))
data_down <- data.frame(geneID=RocA_sensitive_entrezID$ENTREZID,group=group_down)

# all
data_all_GO <- rbind(data_up,data_down)

formula <- compareCluster(geneID~group, data=data_all_GO, fun='enrichGO',OrgDb='org.Hs.eg.db',pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05,ont="BP")
summary(formula)
formula <- simplify(formula)
y <- dotplot(formula,showCategory=20,includeAll=TRUE)
y

formula_dataframe <- as.data.frame(formula)
write.table(formula_dataframe,"GO_analysis_dataframe.txt",quote = F,sep="\t",row.names = F,col.names = T)