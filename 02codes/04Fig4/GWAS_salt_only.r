library("data.table")
library(Ropt)
ch=fread("chr_length_HanXRQr2.txt",data.table=F,head=T)
f=list.files(pattern = "_respo.assoc.txt")
d1=fread("stress_response_GWAS_Overlap_Selection.txt",data.table=F,head=T)
dd=NULL
for (k in 1:17)
{
  sub=subset(d1,d1[,1]==k)
  sub[,5]=sub[,5]+ch[k,3]
  dd=rbind(dd,sub)
}
dd$minP_pos=dd$minP_pos/1e6
dd$pch=25
dd[grep("B2",dd[,8]),9]=16
dd$col=1
col=c("blue","hotpink","#A020F0","orange","chartreuse4")
col=adjustcolor(col,alpha.f=0.6)
dd[dd[,8]=="wild_land.XPCLR",10]=col[1]
dd[dd[,8]=="land_modern.XPCLR",10]=col[2]
dd[dd[,8]=="wild26.B2",10]=col[3]
dd[dd[,8]=="landrace.B2",10]=col[4]
dd[dd[,8]=="modern.B2",10]=col[5]
res=NULL
for(i in f)
{
  d=fread(i, header=T,data.table=F)[,c(1,3,13)]
  d$Trait=gsub(".assoc.txt","",i)
  res=rbind(res,d)
}

d2=NULL
for (k in 1:17)
{
  sub=subset(res,res[,1]==k)
  sub[,2]=sub[,2]+ch[k,3]
  d2=rbind(d2,sub)
}
d2$ps=d2$ps/1e6
d2$col=sapply(d2[,1],function(x){if(x%%2==1){x="black"}else{x="gray"}})
unique(d2$Trait)
tiff("Manhatan_salt_stress2.tiff",res=600,units = "mm",height = 60,width = 210)
par(mar=c(0,4,0,0),mfrow=c(3,1))
tr="Germ_respo"
d3=d2[d2[,4]==tr,]
plot(d3[,2],-log10(d3[,3]),col=d3$col,pch=16,cex=0.5,xlab="",ylab="-log10(p)",axes=F)
axis(2,las=2,tck=-0.02)
#axis(1,at=ch[,4]/1e6,labels=1:17,tck=-0.02)
box()
abline(h=-log10(1/33801),col="red",lty=2)
dd1=subset(dd,dd$Trait==tr)
points(dd1[,5],-log10(dd1[,4]),pch=dd1$pch,col=dd1$col,cex=1,bg=dd1$col)

tr="Radicula_respo"
d3=d2[d2[,4]==tr,]
plot(d3[,2],-log10(d3[,3]),col=d3$col,pch=16,cex=0.5,xlab="",ylab="-log10(p)",axes=F)
axis(2,las=2,tck=-0.02)
axis(1,at=ch[,4]/1e6,labels=1:17,tck=-0.02)
box()
abline(h=-log10(1/33801),col="red",lty=2)
dd1=subset(dd,dd$Trait==tr)
points(dd1[,5],-log10(dd1[,4]),pch=dd1$pch,col=dd1$col,cex=1,bg=dd1$col)
dev.off()

pdf("legend.pdf",height = 4,width = 5)
par(mar=c(4,4,2,2),mfrow=c(1,1))
plot(1:5,col=col,bg=col,pch=c(25,25,24,24,24),cex=2)
dev.off()

###plot effect size
d=fread("All_sweep_effects_with_ancestral_allele_and_effect_salt.txt",header=T,data.table=F)

pdf("Stress_effect_size.pdf",height = 3,width = 8.6)
par(mar=c(3,2,2,1),mfrow=c(1,4))
d1=d[d[,14]=="Germ_respo" & d[,6]=="wild26.B2_99",]
d2=d[d[,14]=="Germ_respo" & d[,6]=="landrace.B2_99",]
d3=d[d[,14]=="Germ_respo" & d[,6]=="modern.B2_99",]
d4=d[d[,14]=="Germ_respo" & d[,6]=="wild_land.XPCLR.99",]
d5=d[d[,14]=="Germ_respo" & d[,6]=="land_modern.XPCLR.99",]
d1=d1[order(d1[,20]),]
d2=d2[order(d2[,20]),]
d3=d3[order(d3[,20]),]
d4=d4[order(d4[,20]),]
d5=d5[order(d5[,20]),]
plot(d1[,20],1:nrow(d1),type="b",pch=21,xlab="Ancestor Allele effect",ylab="Sweep SNPs",main="SNPs under balancing selection",axes=F,col="#A020F080",lwd=1.2,ylim=c(1,35))
axis(1,las=1,tick=-0.02)
axis(2,las=2,tick=-0.02)
abline(v=0,col="red",lty=2,lwd=1.2)
grid()
points(d2[,20],1:nrow(d2),type="b",pch=21,col="orange",lwd=1.2)
points(d3[,20],1:nrow(d3),type="b",pch=21,col="chartreuse4",lwd=1.2)

plot(d4[,20],1:nrow(d4),type="b",pch=25,xlab="Ancestor Allele effect",ylab="Sweep SNPs",main="SNPs under positive selection",axes=F,col="blue",ylim=c(1,140))

axis(1,las=1,tick=-0.02)
axis(2,las=2,tick=-0.02)
abline(v=0,col="red",lty=2)
grid()
points(d5[,20],1:nrow(d5),type="b",pch=25,col="hotpink")

d1=d[d[,14]=="Radicula_respo" & d[,6]=="wild26.B2_99",]
d2=d[d[,14]=="Radicula_respo" & d[,6]=="landrace.B2_99",]
d3=d[d[,14]=="Radicula_respo" & d[,6]=="modern.B2_99",]
d4=d[d[,14]=="Radicula_respo" & d[,6]=="wild_land.XPCLR.99",]
d5=d[d[,14]=="Radicula_respo" & d[,6]=="land_modern.XPCLR.99",]
d1=d1[order(d1[,20]),]
d2=d2[order(d2[,20]),]
d3=d3[order(d3[,20]),]
d4=d4[order(d4[,20]),]
d5=d5[order(d5[,20]),]
plot(d1[,20],1:nrow(d1),type="b",pch=21,xlab="Ancestor Allele effect",ylab="Sweep SNPs",main="SNPs under balancing selection",axes=F,col="#A020F080",lwd=1.2,ylim=c(1,50))
axis(1,las=1,tick=-0.02)
axis(2,las=2,tick=-0.02)
abline(v=0,col="red",lty=2,lwd=1.2)
grid()
points(d2[,20],1:nrow(d2),type="b",pch=21,col="orange",lwd=1.2)
points(d3[,20],1:nrow(d3),type="b",pch=21,col="chartreuse4",lwd=1.2)

plot(d4[,20],1:nrow(d4),type="b",pch=25,xlab="Ancestor Allele effect",ylab="Sweep SNPs",main="SNPs under positive selection",axes=F,col="blue",ylim=c(1,120))

axis(1,las=1,tick=-0.02)
axis(2,las=2,tick=-0.02)
abline(v=0,col="red",lty=2)
grid()
points(d5[,20],1:nrow(d5),type="b",pch=25,col="hotpink")
dev.off()
