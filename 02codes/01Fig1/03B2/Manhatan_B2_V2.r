library("data.table")
library(Ropt)
ch=fread("chr_length_HanXRQr2.txt",data.table=F,head=T)
d=fread("Gene_need2_highlight.txt",data.table=F,head=T)
a=NULL
for(i in 1:nrow(d))
{
  b=d[i,]
  b[,4]=b[,4]+ch[b[1,1],3]
  a=rbind(a,b)
}
d=a
d[,4]= d[,4]/1e6
d1=d[d[,12]=="Wild" & d[,18]=="Stress",]
d2=d[d[,12]=="Wild" & d[,18]=="Salt_ions",]
  file="Sunflower_SAM_HanXRQr2_wild26.B2_stat.txt"
  ylab="B2(Wild)"
  out=gsub("txt","tiff",file)
  tiff(out,width = 180, height = 100,res=600,units="mm" )
  re=fread(file,data.table=F)
  re=re[,c(1,2,4)]
  re=re[order(re[,1],re[,2]),]
  colnames(re)=c("Chr","Pos","B2")
  res=NULL
  for (i in 1:17)
  {
    sub=subset(re,re[,1]==i)
    sub[,2]=sub[,2]+ch[i,3]
    res=rbind(res,sub)
  }
  res[,2]= res[,2]/1000000
  col=sapply(re[,1],function(x){if(x%%2==1){x="black"}else{x="gray"}})
  plot(res[,2],res[,3],col=col,pch=16,cex=0.3,ylim=c(0,max(res[,3])),bty="l",xlab="Chromosome",ylab=ylab,axes=F)
  abline(h=quantile(re[,3],0.99),lwd=1,lty=2,col="red")
  axis(2,las=2,tck=-0.02)
  axis(1,at=ch[,4]/1e6,labels=1:17,tck=-0.02)
  box(bty="l")
  
  points(d1[,4],d1[,6],col=adjustcolor("orange",alpha.f = 1),pch=16,cex=1)
  points(d2[,4],d2[,6],col=adjustcolor("red",alpha.f = 1),pch=16,cex=1)
  dev.off()
  