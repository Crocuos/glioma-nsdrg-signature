library(pheatmap)

#CGGA#TCGA
datw = read.delim('tcga11.txt', head=T, row.names = 1) 
anno_colw=read.delim('survival11.txt', head=T, row.names = 1)
#anno_rowOA =read.delim(file.choose(), head=T, row.names = 1)
pheatmap(log((datw+1),2),annotation_col=anno_colw,legend=T, 
         col = colorRampPalette(c("blue","white","red"))(120), 
         breaks = seq(-0,6,by=.05),
         show_colnames=F,show_rownames=T, cluster_cols=F, cluster_rows=F)
