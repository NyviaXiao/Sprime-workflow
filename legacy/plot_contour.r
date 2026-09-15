library(MASS)
args = commandArgs(trailingOnly=TRUE)


# da=read.csv(file=args[1],header=T,sep="\t")  
da=read.csv(file=args[1],header=T,sep=" ")  

png(paste(args[2],".png",sep=""),width=1100,height=1100, res=250); 
# ylab="Match to Altai Denisovan"; 
# xlab="Match to Altai Neanderthal"; 

xlab="Match to Denisovan3"; 
ylab="Match to Denisovan25"; 

level1=seq(0.3,0.9,0.1); 
level2=seq(1,10,1)
X<-da$Denisovan3
Y<-da$Denisovan25
# X<-da$AltaiNean
# Y<-da$AltaiDeni
o<-!is.na(X)&!is.na(Y)
X<-X[o]
Y<-Y[o]

contour(kde2d(X,Y,n=100, lims=c(0,1,0,1)),levels=c(level1,level2), xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=args[3],las=1,cex.lab=1.3,cex.axis=1.2,cex.main=1.5,lty=5,labcex=0.8);grid(lty=3,col="gray")
dev.off()
