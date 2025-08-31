library(clusterProfiler)#??Ҫ?õ?clusterProfiler??
df=read.table("TCGA.txt",header=T,as.is=T)
head(df)

dim(df)

df.id=bitr(df$SYMBOL,
fromType="SYMBOL",
toType="ENTREZID",
OrgDb="org.Hs.eg.db")
head(df.id)
head(df)
easy.df<-merge(df,df.id,by= "SYMBOL",all=F)
sortdf<-easy.df[order(easy.df $logFC, decreasing = T),] 
head(sortdf)
gene.expr = sortdf$logFC
head(gene.expr)
names(gene.expr) <- sortdf $ENTREZID
head(gene.expr)

kk<-gseGO(gene.expr,
          ont = "BP", 
          keyType = "ENTREZID",
          OrgDb = org.Hs.eg.db,
          pvalueCutoff = 0.05, 
          pAdjustMethod = "BH", verbose = TRUE,
          seed = FALSE,
          by = "fgsea")
head(kk)
dim(kk)
sortkk<-kk[order(kk$enrichmentScore, decreasing = T),]
head(sortkk)
dim(sortkk)
write.table(sortkk, "gsea_TCGA_score.txt",sep = "\t",quote = T,col.names = T,row.names = T)
gseaplot(kk, "hsa04510")
gseaplot(kk, "hsa04510", color.line= 'steelblue')
library(enrichplot)
gseaplot2(kk, "hsa04510")
gseaplot2(kk, "GO:0042116", color = "firebrick", rel_heights=c( 1, .2, .6))
gseaplot2( kk, row.names( sortkk) [1:4])
paths<- c( "GO:0042102", "GO:0006959", "GO:0042129", "GO:0019724",'GO:0002456')#GO:0042102 tϸ????ֳ GO:0042129 tϸ?????? 

gseaplot2(kk, paths,color = "firebrick",title = "CGGA_325 dataset", base_size = 18, rel_heights=c( 1, .2, .6))
gseaplot2(kk, paths, subplots=1) 

gseaplot2(kk, paths, subplots=1:2) 

gseaplot2(kk, paths, subplots=c(1,3)) 
gseaplot2(kk, 
          
          row.names(sortkk)[ 39], 
          
          title = "Prion disease", 
          
          base_size = 15, 
          
          color = "green", 
          
          pvalue_table = TRUE, 
          
          ES_geom= "line") 



######up###
#############
############
df=read.table("scorediffup.txt",header=T,as.is=T)
head(df)

dim(df)

df.id=bitr(df$SYMBOL,
           fromType="SYMBOL",
           toType="ENTREZID",
           OrgDb="org.Hs.eg.db")
head(df.id)
head(df)
easy.df<-merge(df,df.id,by= "SYMBOL",all=F)
sortdf<-easy.df[order(easy.df $logFC, decreasing = T),] #????????
head(sortdf)
gene.expr = sortdf$logFC
head(gene.expr)
names(gene.expr) <- sortdf $ENTREZID#????????ȡ??foldchange??Ӧ??ENTREZID
head(gene.expr)
#kk<- gseKEGG(gene.expr, organism = "hsa")  ??KEGG??ѡHSA??ͼ??
kk<-gseGO(gene.expr,
          ont = "BP", 
          keyType = "ENTREZID",
          OrgDb = org.Hs.eg.db,
          pvalueCutoff = 0.05, 
          pAdjustMethod = "BH", verbose = TRUE,
          seed = FALSE,
          by = "fgsea")
head(kk)
dim(kk)
sortkk<-kk[order(kk$enrichmentScore, decreasing = T),]
head(sortkk)
dim(sortkk)
write.table(sortkk, "gsea_CGGA_score-UP.txt",sep = "\t",quote = T,col.names = T,row.names = T)
gseaplot(kk, "hsa04510")
gseaplot(kk, "hsa04510", color.line= 'steelblue')
library(enrichplot)
gseaplot2(kk, "hsa04510")
gseaplot2(kk, "GO:0042116", color = "firebrick", rel_heights=c( 1, .2, .6))
gseaplot2( kk, row.names( sortkk) [1:4])
paths<- c( "GO:0042102", "GO:0006959", "GO:0042129", "GO:0019724",'GO:0002456')#GO:0042102 tϸ????ֳ GO:0042129 tϸ?????? 
#GO:0019724 Bϸ??????
gseaplot2(kk, paths,color = "firebrick",title = "CGGA_325 dataset", base_size = 18, rel_heights=c( 1, .2, .6))
gseaplot2(kk, paths, subplots=1) #ֻҪ??һ??ͼ

gseaplot2(kk, paths, subplots=1:2) #ֻҪ??һ?͵ڶ???ͼ

gseaplot2(kk, paths, subplots=c(1,3)) #ֻҪ??һ?͵?????ͼ
gseaplot2(kk, paths, subplots=c(1,2,3)) #ֻҪ??һ?͵?????ͼ
gseaplot2(kk, #????
          
          row.names(sortkk)[ 39], #????һ?е??ź?ͨ·
          
          title = "Prion disease", #????
          
          base_size = 15, #??????С
          
          color = "green", #????????ɫ
          
          pvalue_table = TRUE, #?Ӳ???pֵ
          
          ES_geom= "line") #?????ߣ???????d??

