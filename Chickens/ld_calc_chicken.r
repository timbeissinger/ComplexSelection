####################################################################
#
#        LD-Calculation Chicken Data
#
####################################################################

a <- readRDS("zMat_white.RDS")
amap <- read.table("markerkarte_white.txt",stringsAsFactors=FALSE, header=TRUE)
all(order(amap[,2],amap[,3])==1:nrow(amap))
#[1] TRUE    --> Check whether SNPs are already ordered

### function to calculate correlation pairwise --> more efficient than cor:
cor1 <- function(x1,x2){
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2) 
  sx1 <- t(x1)-colMeans(x1) 
  sx2 <- t(x2)-colMeans(x2)
  cc <- (rowSums(sx1*sx2))/(sqrt(rowSums(sx1^2)*rowSums(sx2^2)))
  return(cc)
}

### package for parallelization
library("parallel")    

### Max distant between SNPs in #SNPs to calculate
mdiff <- 2000     
  
### 32 chromosomes:
alld1 <- matrix(0,nrow=32,ncol=mdiff)

for(rr in 32:1){      ## loop starts with smallest chromosome
       print(paste("Chromosome",rr))
       af <- a[,amap[,2]==rr,drop=FALSE]
       if(ncol(af)==0) next
       print(dim(af))
       ind <- 1:ncol(af)

       fun <- function(i){  # 1:mdiff; function to calculate r^2 for all SNPs with distance i to each other
           if(i>=ncol(af)) next
           ld1 <- mean(cor1(af[,head(ind,-i)],af[,tail(ind,-i)])^2,na.rm=TRUE) 
       return(ld1)
       }
       
       ## calculation of r^2 ; parallelization using 10 cores
       ld1 <- unlist(mclapply(1:2000, FUN=fun,  mc.preschedule = TRUE, mc.cores = 10))

       alld1[rr,] <- ld1
}

save(alld1,file="Res_white_Chrwise.RData")



##########################################################################
#
#        Plotting r^2 results
#
##########################################################################

tag <- "white"       # Brown or white layer

## loading data
a <- readRDS(paste("zMat_",tag, ".RDS",sep=""))
amap <- read.table(paste("markerkarte_",tag,".txt",sep=""),stringsAsFactors=FALSE, header=TRUE)
laeng <- unlist(lapply(by(amap[,3],amap[,2],diff),FUN=median))        ## median distance between snps, chromosomewise 

## load r^2 results
load(file=paste("Res_",tag,"_Chrwise.RData",sep=""))
alld1 <- alld1[-c(29,30),]                ## chromosome 29 and 30 had no SNPs data

## function to get the minimal distance with r^2 < 0.05
f <- function(x){return(min(which(x<0.05)))}
dd <- apply(alld1,1,f)

ltt <- rep(1:5,each=8)   ## to set the line types

## set screens for split.screen
sc <- rbind(c(0,0.8,0,1),c(0.75,1,0,1),c(0.4,0.75,0.5,1))
split.screen(sc)

## Plot LD-curves
screen(1)
plot(alld1[1,],ylab="r^2", xlab="Distance in #SNPs", type="l", ylim=c(0,0.65),col="white",bty="n",main=paste(ifelse(tag=="brown","Brown","White")," Layer;",ncol(a),"SNPs and",nrow(a),"Individuals" ))
for(i in 1:nrow(alld1)){
      lines(alld1[i,], lwd=2, col=i,lty=ltt[i])
}
abline(h=0.05,col="red",lwd=2)

## Plot Legends
screen(2)
par(mai=c(0,0,0,0))
plot(1:1000,1:1000, col="white",axes=FALSE,xlab="", ylab="",bty="n")
legend(1,900, legend=c("Chr: #SNPs : Thres.(LD<0.05)",paste(c(1:32)[-c(29,30)],tab,dd,sep="   :   ")),col=0:nrow(alld1),lty=c(0,ltt[1:nrow(alld1)]),lwd=2,bty="n",cex=0.9)

## Plot mean distance versus minimaö distance
screen(3)
x <- as.numeric(laeng)
y <- dd
fit <-  lm(y~x)$coefficients
plot(x, y,type="p",xlab="Mean Distance (bp) between SNPs", ylab="Threshold",pch=18,col=1:nrow(alld1))
abline(fit,col="red",lwd=2)

close.screen(all=T)




