

install.packages('corrgram')
library(corrgram)
data = read.delim("corgram-signature.txt", row.names = 1)
corrgram(data, order=F, lower.panel=panel.shade, upper.panel=panel.pie, text.panel=panel.txt,main="Correlogram of phenotype correlated markers", col.regions=colorRampPalette(c("green1", "white","firebrick1")))##??ͼ
