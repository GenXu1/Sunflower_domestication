library(data.table)
library(scatterplot3d)
a=fread("wild.txt", header=F,data.table=F)
b=fread("modern.txt", header=F,data.table=F)
c=fread("landrace.txt", header=F,data.table=F)
pc=fread("Sunflower_SAM_HanXRQr2_dp1_miss05_maf001.eigenvec", header=F,data.table=F)[,c(2:5)]
eig=fread("Sunflower_SAM_HanXRQr2_dp1_miss05_maf001.eigenval", header=F,data.table=F)[,1]
pc.percent <- 100 * eig[1:10]/sum(eig)
pc.percent=round(pc.percent,2)
plot(1:10,pc.percent,type="b",main="",xlab="PCs")

PC1=pc[,2]
PC2=pc[,3]
PC3=pc[,4]

col=c("#A020F080","orange","chartreuse4")
colnames(pc)=c("ID","PC1","PC2","PC3")
pc$col="gray"
pc[pc[,1]%in%a[,1],5]=col[1]
pc[pc[,1]%in%b[,1],5]=col[2]
pc[pc[,1]%in%c[,1],5]=col[3]


pdf("sunflower.PCA.2D_26wilds_V2.pdf",height = 5,width = 6)
par(mar=c(4,4,2,2),mfrow=c(1,1))
plot(pc[,2], pc[,3],col=pc[,5],pch=16, xlab=qq("PC1 ({pc.percent[1]}%)"), ylab=qq("PC2 ({pc.percent[2]}%)"),las=1,tck=-0.02,bg=col,cex=pc$cex)
legend("topright",c("Wild","Landrace","Modern"),pch=16,col=col[c(1,3,2)],bty="n",ncol=1) ##??????,right:????????????,c("Teosinte","Landrace","Trop")??????????????????,bty???????????????
dev.off()

