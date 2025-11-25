library("data.table")
library(Ropt)
########wild vs landrace
ch=fread("chr_length_HanXRQr2.txt",data.table=F,head=T)
d=fread("Genes_need_highlight.txt",data.table=F,head=T)
sw=fread("wild_land_ov_land_modern.txt",data.table=F,head=F)
a=NULL
for(i in 1:nrow(d))
{
  d1=d[i,]
  d1[,2]=d1[,2]+ch[ch[,1]==d1[1,1],3]
  d1[,3]=d1[,3]+ch[ch[,1]==d1[1,1],3]
  a=rbind(a,d1)
}
d=a

s=NULL
for(i in 1:nrow(sw))
{
  d1=sw[i,]
  d1[,2]=d1[,2]+ch[ch[,1]==d1[1,1],3]
  d1[,3]=d1[,3]+ch[ch[,1]==d1[1,1],3]
  d1[,6]=d1[,6]+ch[ch[,1]==d1[1,1],3]
  d1[,7]=d1[,7]+ch[ch[,1]==d1[1,1],3]
  s=rbind(s,d1)
}
sw=s

d1=d[d[,10]=="WD_vs_LR",]
file="wild_land.xpclr.txt";ylab="XP-CLR (WD_vs_LR)"
out=gsub("txt","tiff",file)
tiff(out,width = 210, height =80,res=600,units="mm" )
re=fread(file,data.table=F)
re=re[,c(2,3,12)]
re=re[order(re[,1],re[,2]),]
colnames(re)=c("Chr","Pos","XP-CLR")
res=NULL
for (i in 1:17)
{
  sub=subset(re,re[,1]==i)
  sub[,2]=sub[,2]+ch[i,3]
  res=rbind(res,sub)
}
res[,2]= res[,2]/1000000
col=sapply(re[,1],function(x){if(x%%2==1){x="black"}else{x="gray"}})
plot(res[,2],res[,3],col=col,pch=16,cex=0.5,ylim=c(0,max(res[,3])),bty="l",xlab="Chromosome",ylab=ylab,axes=F)
abline(h=quantile(re[,3],0.99),lwd=2,lty=2,col="red")
axis(2,las=2)
axis(1,at=ch[,4]/1e6,labels=1:17,tck=-0.02)
box(bty="l")
d2=d1[,c(1,2,4,15)]
d2$col="red"
d2[d2[,4]=="OIL",5]="blue"
d2[d2[,4]=="stress",5]="orange"
d2[d2[,4]=="FT",5]=colours()[257]
points(d2[,2]/1e6,d2[,3],col=d2[,5],pch=16,cex=0.6)
points(sw[,2]/1e6,sw[,4],col=adjustcolor("gold1",alpha.f = .7),pch=17,cex=0.6,type="h",lwd=1.5)
# legend("topright",legend=c("Yield","Oil","Stress","Flowering Time"),col=c("red","blue","orange",colours()[257]),pch=16,bty="n",cex=0.8)
dev.off()




########landrace vs modern
ch=fread("chr_length_HanXRQr2.txt",data.table=F,head=T)
d=fread("Genes_need_highlight.txt",data.table=F,head=T)
sw=fread("wild_land_ov_land_modern.txt",data.table=F,head=F)
a=NULL
for(i in 1:nrow(d))
{
  d1=d[i,]
  d1[,2]=d1[,2]+ch[ch[,1]==d1[1,1],3]
  d1[,3]=d1[,3]+ch[ch[,1]==d1[1,1],3]
  a=rbind(a,d1)
}
d=a

s=NULL
for(i in 1:nrow(sw))
{
  d1=sw[i,]
  d1[,2]=d1[,2]+ch[ch[,1]==d1[1,1],3]
  d1[,3]=d1[,3]+ch[ch[,1]==d1[1,1],3]
  d1[,6]=d1[,6]+ch[ch[,1]==d1[1,1],3]
  d1[,7]=d1[,7]+ch[ch[,1]==d1[1,1],3]
  s=rbind(s,d1)
}
sw=s

d1=d[d[,10]=="LR_vs_MD",]
file="land_modern.xpclr.txt";ylab="XP-CLR (LR vs MD)"
  out=gsub("txt","tiff",file)
  tiff(out,width = 210, height =80,res=600,units="mm" )
  re=fread(file,data.table=F)
  re=re[,c(2,3,12)]
  re=re[order(re[,1],re[,2]),]
  colnames(re)=c("Chr","Pos","XP-CLR")
  res=NULL
  for (i in 1:17)
  {
    sub=subset(re,re[,1]==i)
    sub[,2]=sub[,2]+ch[i,3]
    res=rbind(res,sub)
  }
  res[,2]= res[,2]/1000000
  col=sapply(re[,1],function(x){if(x%%2==1){x="black"}else{x="gray"}})
  plot(res[,2],res[,3],col=col,pch=16,cex=0.5,ylim=c(0,max(res[,3])),bty="l",xlab="Chromosome",ylab=ylab,axes=F)
  abline(h=quantile(re[,3],0.99),lwd=2,lty=2,col="red")
  axis(2,las=2)
  axis(1,at=ch[,4]/1e6,labels=1:17,tck=-0.02)
  box(bty="l")
  d2=d1[,c(1,2,4,15)]
  d2$col="red"
  d2[d2[,4]=="OIL",5]="blue"
  d2[d2[,4]=="stress",5]="orange"
  d2[d2[,4]=="FT",5]=colours()[257]
  points(d2[,2]/1e6,d2[,3],col=d2[,5],pch=16,cex=0.6)
  points(sw[,6]/1e6,sw[,8],col=adjustcolor("gold1",alpha.f = .7),pch=17,cex=0.6,type="h",lwd=1.5)
 # legend("topright",legend=c("Yield","Oil","Stress","Flowering Time"),col=c("red","blue","orange",colours()[257]),pch=16,bty="n",cex=0.8)
  dev.off()