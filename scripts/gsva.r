
genefile=read.table("genelist.txt",
                    header=F,
                    row.names = 1,
                    sep="\t")
genefile[1:4,1:4]
genefile=genefile[,-1]     
genefile[1:4,1:4]
genelist<-apply(genefile,1,function(x){
  (x[x!=""])
})


data=read.table("CGGA.txt",
                head=T,
                row.names = 1,
                sep="\t") 
data[1:4,]

library(GSVA)
library(pheatmap)
gsva=gsva(data.matrix(data), genelist,  verbose=FALSE)
pheatmap(gsva,cluster_cols = F)
write.table(gsva,
            "CGGA.gsva_RESULT.txt",
            col.names = NA,
            sep="\t",
            quote=F)
