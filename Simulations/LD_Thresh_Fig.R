###### Use this script to make an LD-Threshold figure

### Timothy Beissinger
### 10-11-2017

setwd("/home/beissinger/Documents/ComplexSelection/Manuscript/Figures/")

Meff_vector <- c()

Decay1000 <-c()
FalsePos1000 <- c()

Decay100 <-c()
FalsePos100 <- c()

Decay50 <-c()
FalsePos50 <- c()

Decay10 <-c()
FalsePos10 <- c()

### Load LD-Thresh workspaces (1000) and save R2, p-value
for (line in 44:58){
  print(line)
  load(paste("../../LD_Thresh_1000/QTL-1000_L",line,".RData",sep=""))
  Decay1000[line-43] <- mean(R2)
  FalsePos1000[line-43] <- length(which(Pvals.drift <= 0.05))/length(Pvals.drift)
  Meff_vector[line-43] <- Meff
}

### Load LD-Thresh workspaces (100) and save R2, p-value
for (line in 44:58){
  print(line)
  load(paste("../../LD_Thresh_100/QTL-100_L",line,".RData",sep=""))
  Decay100[line-43] <- mean(R2)
  FalsePos100[line-43] <- length(which(Pvals.drift <= 0.05))/length(Pvals.drift)
}

### Load LD-Thresh workspaces (50) and save R2, p-value
for (line in 44:58){
  print(line)
  load(paste("../../LD_Thresh_50/QTL-50_L",line,".RData",sep=""))
  Decay50[line-43] <- mean(R2)
  FalsePos50[line-43] <- length(which(Pvals.drift <= 0.05))/length(Pvals.drift)
}

### Load LD-Thresh workspaces (10) and save R2, p-value
for (line in 44:58){
  print(line)
  load(paste("../../LD_Thresh_10/QTL-10_L",line,".RData",sep=""))
  Decay10[line-43] <- mean(R2)
  FalsePos10[line-43] <- length(which(Pvals.drift <= 0.05))/length(Pvals.drift)
}

### Make table
Decay <- cbind(Decay1000,Decay100,Decay50,Decay10)


info <- cbind(Meff_vector,apply(Decay,1,mean),FalsePos1000,FalsePos100,FalsePos50,FalsePos10)


### Make Plot
jpeg("LD_Thresh.jpg",height=600*1.2,width=800*1.2,pointsize=18)
par(mar=c(6,4,4,2)+.1)
plot(Decay1000,FalsePos1000,type="l",lty=2,col=2,ylab="False Positive Rate",xlab="",xaxt="n",lwd=2)
mtext("LD Threshold \n (Effective Marker Number)",side=1,line=4.5)
#lines(Decay1000,FalsePos1000,lty=2,col=1)
lines(Decay100,FalsePos100,lty=3,col=3,lwd=2)
lines(Decay50,FalsePos50,lty=4,col=4,lwd=2)
lines(Decay10,FalsePos10,lty=5,col=5,lwd=2)
abline(h=0.05,col="brown",lwd=2)
axis(1,round(info[c(1,2,3,4,8,14),2],3),at=info[c(1,2,3,4,8,14),2],las=1,line=0)
mtext(paste("(",info[c(1,2,3,4,8,14),1],")",sep=""),1,at=info[c(1,2,3,4,8,14),2],line=2)
legend("topleft","c(x,y)",c("1,000 QTL","100 QTL","50 QTL","10 QTL","0.05 False positive rate"),lwd=2,lty=c(2,3,4,5,1),col=c(2,3,4,5,"brown"),cex=0.75)
dev.off()

### Write False positive table
colnames(info)[2] <- "R2"
write.table(info,file="LD_Thresh_Table.txt",sep="\t",row.names=F,quote=F)

