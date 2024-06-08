library(CMplot)
df <- read.table("~/Foxtail-millet-data-analysis/Figure1_and_Figure5/data/OTU_106.assoc.txt",header = T)
CMplot(head(df,n = 2000), plot.type = "m", LOG10=TRUE, cex=0.4, threshold=2.2E-5,, threshold.lty=2, threshold.col="black",  amplify=TRUE, signal.col="red", signal.cex=1,width=14, bin.size=1e6, chr.den.col=c("darkgreen", "yellow", "red"), verbose=TRUE, height=5, file="pdf")
data(pig60K) 
head(pig60K)
head(df)
CMplot(pig60K,type="p",plot.type="m",LOG10=TRUE,threshold=NULL,file="jpg",file.name="",dpi=300,
       file.output=TRUE,verbose=TRUE,width=14,height=6,chr.labels.angle=45)
