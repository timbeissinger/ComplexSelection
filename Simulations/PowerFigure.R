########################################################
### This script will create a figure of the power of ###
### Ghat compared to selection mapping.              ###
########################################################

### Tim beissinger
### 8-31-2017

### Make vectors of TP and FP rate for Ghat.
Ghat_TP <- c(0.04,0.54,.94,1)
Ghat_FP <- c(0.03,0.03,0.02,0.03)

Ghat_TPSD <- sqrt(Ghat_TP*(1-Ghat_TP)/100) #standard deviation

### Load info to calculate TP and FP rate for mapping.
load("../../MapSims/Map10.RData")
tp10 <- tp/10
fp10 <- fp/10

load("../../MapSims/Map50.RData")
tp50 <- tp/50
fp50 <- fp/50

load("../../MapSims/Map100.RData")
tp100 <- tp/100
fp100 <- fp/100

load("../../MapSims/Map1000.RData")
tp1000 <- tp/1000
fp1000 <- fp/1000

tpMat <- cbind(tp10,tp50,tp100,tp1000)

Map_tp <- apply(tpMat,2,mean)
Map_tpsd <- apply(tpMat,2,sd)

fpMat <- cbind(fp10,fp50,fp100,fp1000)
Map_fp <- apply(fpMat,2,mean)
Map_fpsd <- apply(fpMat,2,sd)

### Make True positives plot
jpeg("detectionRate.jpg",height=480,width=640)
plot(Ghat_TP,pch=19,ylim=c(0,1.1),xaxt="n",cex=0.65,xlab="Number of QTL Simulated",ylab="Detection rate",col="brown")
lines(Ghat_TP,col="brown")
points(Map_tp,col="blue",pch=19,cex=0.65)
lines(Map_tp,col="blue")
segments(1:4,Map_tp-Map_tpsd,1:4,Map_tp+Map_tpsd,col="blue")
segments(1:4,Ghat_TP-Ghat_TPSD,1:4,Ghat_TP+Ghat_TPSD,col="brown")
#points(Ghat_FP,pch=4,col="brown")
#lines(Ghat_FP,col="brown")
axis(1,at=c(1,2,3,4),lab=c(10,50,100,1000))
legend("topleft","c(x,y)",lty=1,pch=19,col=c("brown","blue"),c("G-hat","Selection Mapping"))
dev.off()

# ### Make false positives plot
# plot(Ghat_FP,pch=19,col="brown",ylim=c(0,8),xlab="Number of QTL Simulated",ylab="False positives/True positve",xaxt="n")
# lines(Ghat_FP,col="brown")
# points(Map_fp,col="blue",pch=19)
# lines(Map_fp,col="blue")
# #segments(1:4,Map_fp-Map_fpsd,1:4,Map_fp+Map_fpsd,col="blue")
# axis(1,at=c(1,2,3,4),lab=c(10,50,100,1000))
# legend("topleft","c(x,y)",lty=1,pch=19,col=c("brown","blue"),c("G-hat","Selection Mapping"))

### Make and export a table of detection rate and number of false positives
table <- rbind(Ghat_TP,Map_tp,Ghat_FP,Map_fp)
colnames(table) = c(10,50,100,1000)
write.table(table,file="TPFP.csv",row.names=T,col.names=T,sep=",")

##################################################################
######################EYE PLOTS - 10 #############################
##################################################################
dir <- "DataForFigs/10."
iter <- "001"
source("DataForFigs/Ghat.R")

pheno <- read.table(paste(dir, "SelectPop_data_",iter,".txt",sep=""),header=T,stringsAsFactors=F)
map <- read.table(paste(dir, "lm_mrk_",iter,".txt",sep=""),header=T,stringsAsFactors=F)
geno <- read.table(paste(dir, "SelectPop_mrk_",iter,".txt-head",sep=""),header=F,stringsAsFactors=F,skip=1,sep="", colClasses=c("numeric","character")) # read genos as characters


#load allele frequencies
freqs0<-read.table(paste(dir, "SelectPop_freq_mrk_",iter,".txt",sep=""),header=T,nrows=nrow(map),fill=T,stringsAsFactors=F)
freqs20<-read.table(paste(dir, "SelectPop_freq_mrk_",iter,".txt",sep=""),header=F,skip={4*nrow(map)+1},fill=T,stringsAsFactors=F)

### Manipulate genotypes by coding to -1,0,1 and so that markers are columns, individuals are rows.
gen <- matrix(NA,nrow=nrow(map),ncol=nrow(geno))
gen <- as.data.frame(gen)
names(gen) <- geno[,1]
for(i in 1:1000){
  print(i)
  tmp <- as.numeric(unlist(strsplit(geno[i,2],split="")))
  tmp[which(tmp == 0)] <- -1
  tmp[which(tmp == 3 | tmp ==4)] <- 0
  tmp[which(tmp==2)] <- 1
  gen[,i] <- tmp
}
gen<-t(gen)
gc()

### Calculate allele frequencies
#Gen 0
names(freqs0)[4]<- "Allele1"
names(freqs0)[5]<- "Allele2"
freqs0$Allele2[which(substr(freqs0$Allele1,1,1)==2)] <- "2:1.00000" # put this in the spot for the second allele
freqs0$Allele1[which(substr(freqs0$Allele1,1,1)==2)] <- "1:0.00000" 
freqs0$Allele2[which(substr(freqs0$Allele1,3,3)==1)] <- "2:0.00000" ##
freqs0$Allele1 <- as.numeric(substr(freqs0$Allele1,3,1000))
freqs0$Allele2 <- as.numeric(substr(freqs0$Allele2,3,1000))

#Gen 20
names(freqs20)[4]<- "Allele1"
names(freqs20)[5]<- "Allele2"
freqs20$Allele2[which(substr(freqs20$Allele1,1,1)==2)] <- "2:1.00000" # put this in the spot for the second allele
freqs20$Allele1[which(substr(freqs20$Allele1,1,1)==2)] <- "1:0.00000" 
freqs20$Allele2[which(substr(freqs20$Allele1,3,3)==1)] <- "2:0.00000"  ##
freqs20$Allele1 <- as.numeric(substr(freqs20$Allele1,3,1000))
freqs20$Allele2 <- as.numeric(substr(freqs20$Allele2,3,1000))

#Calculate change
change2<-freqs20$Allele2-freqs0$Allele2

### Setup for function
geno <- gen
phen <- as.matrix(pheno[1:1000,10])
rownames(phen)<-pheno[1:1000,1]
change=change2
perms <- 10000
blockSize=100000/100

# Run function
test<- Ghat_func(geno=geno,phen=phen,change=change2,method = "scale", num_eff = 6667,  perms=1000,plot="Both")
#Pvals.selection[as.numeric(iter)] <- test$p.val

### Save data for plotting
change10 <- change2
effects10 <- test$effects



##################################################################
######################EYE PLOTS - 50 #############################
##################################################################
dir <- "DataForFigs/50."
iter <- "001"
source("DataForFigs/Ghat.R")

pheno <- read.table(paste(dir, "SelectPop_data_",iter,".txt",sep=""),header=T,stringsAsFactors=F)
map <- read.table(paste(dir, "lm_mrk_",iter,".txt",sep=""),header=T,stringsAsFactors=F)
geno <- read.table(paste(dir, "SelectPop_mrk_",iter,".txt-head",sep=""),header=F,stringsAsFactors=F,skip=1,sep="", colClasses=c("numeric","character")) # read genos as characters


#load allele frequencies
freqs0<-read.table(paste(dir, "SelectPop_freq_mrk_",iter,".txt",sep=""),header=T,nrows=nrow(map),fill=T,stringsAsFactors=F)
freqs20<-read.table(paste(dir, "SelectPop_freq_mrk_",iter,".txt",sep=""),header=F,skip={4*nrow(map)+1},fill=T,stringsAsFactors=F)

### Manipulate genotypes by coding to -1,0,1 and so that markers are columns, individuals are rows.
gen <- matrix(NA,nrow=nrow(map),ncol=nrow(geno))
gen <- as.data.frame(gen)
names(gen) <- geno[,1]
for(i in 1:1000){
  print(i)
  tmp <- as.numeric(unlist(strsplit(geno[i,2],split="")))
  tmp[which(tmp == 0)] <- -1
  tmp[which(tmp == 3 | tmp ==4)] <- 0
  tmp[which(tmp==2)] <- 1
  gen[,i] <- tmp
}
gen<-t(gen)
gc()

### Calculate allele frequencies
#Gen 0
names(freqs0)[4]<- "Allele1"
names(freqs0)[5]<- "Allele2"
freqs0$Allele2[which(substr(freqs0$Allele1,1,1)==2)] <- "2:1.00000" # put this in the spot for the second allele
freqs0$Allele1[which(substr(freqs0$Allele1,1,1)==2)] <- "1:0.00000" 
freqs0$Allele2[which(substr(freqs0$Allele1,3,3)==1)] <- "2:0.00000" ##
freqs0$Allele1 <- as.numeric(substr(freqs0$Allele1,3,1000))
freqs0$Allele2 <- as.numeric(substr(freqs0$Allele2,3,1000))

#Gen 20
names(freqs20)[4]<- "Allele1"
names(freqs20)[5]<- "Allele2"
freqs20$Allele2[which(substr(freqs20$Allele1,1,1)==2)] <- "2:1.00000" # put this in the spot for the second allele
freqs20$Allele1[which(substr(freqs20$Allele1,1,1)==2)] <- "1:0.00000" 
freqs20$Allele2[which(substr(freqs20$Allele1,3,3)==1)] <- "2:0.00000"  ##
freqs20$Allele1 <- as.numeric(substr(freqs20$Allele1,3,1000))
freqs20$Allele2 <- as.numeric(substr(freqs20$Allele2,3,1000))

#Calculate change
change2<-freqs20$Allele2-freqs0$Allele2

### Setup for function
geno <- gen
phen <- as.matrix(pheno[1:1000,10])
rownames(phen)<-pheno[1:1000,1]
change=change2
perms <- 10000
blockSize=100000/100

# Run function
test<- Ghat_func(geno=geno,phen=phen,change=change2,method = "scale", num_eff = 6667,  perms=1000,plot="Both")
#Pvals.selection[as.numeric(iter)] <- test$p.val

### Save data for plotting
change50 <- change2
effects50 <- test$effects



##################################################################
######################EYE PLOTS - 100 #############################
##################################################################
dir <- "DataForFigs/100."
iter <- "001"
source("DataForFigs/Ghat.R")

pheno <- read.table(paste(dir, "SelectPop_data_",iter,".txt",sep=""),header=T,stringsAsFactors=F)
map <- read.table(paste(dir, "lm_mrk_",iter,".txt",sep=""),header=T,stringsAsFactors=F)
geno <- read.table(paste(dir, "SelectPop_mrk_",iter,".txt-head",sep=""),header=F,stringsAsFactors=F,skip=1,sep="", colClasses=c("numeric","character")) # read genos as characters


#load allele frequencies
freqs0<-read.table(paste(dir, "SelectPop_freq_mrk_",iter,".txt",sep=""),header=T,nrows=nrow(map),fill=T,stringsAsFactors=F)
freqs20<-read.table(paste(dir, "SelectPop_freq_mrk_",iter,".txt",sep=""),header=F,skip={4*nrow(map)+1},fill=T,stringsAsFactors=F)

### Manipulate genotypes by coding to -1,0,1 and so that markers are columns, individuals are rows.
gen <- matrix(NA,nrow=nrow(map),ncol=nrow(geno))
gen <- as.data.frame(gen)
names(gen) <- geno[,1]
for(i in 1:1000){
  print(i)
  tmp <- as.numeric(unlist(strsplit(geno[i,2],split="")))
  tmp[which(tmp == 0)] <- -1
  tmp[which(tmp == 3 | tmp ==4)] <- 0
  tmp[which(tmp==2)] <- 1
  gen[,i] <- tmp
}
gen<-t(gen)
gc()

### Calculate allele frequencies
#Gen 0
names(freqs0)[4]<- "Allele1"
names(freqs0)[5]<- "Allele2"
freqs0$Allele2[which(substr(freqs0$Allele1,1,1)==2)] <- "2:1.00000" # put this in the spot for the second allele
freqs0$Allele1[which(substr(freqs0$Allele1,1,1)==2)] <- "1:0.00000" 
freqs0$Allele2[which(substr(freqs0$Allele1,3,3)==1)] <- "2:0.00000" ##
freqs0$Allele1 <- as.numeric(substr(freqs0$Allele1,3,1000))
freqs0$Allele2 <- as.numeric(substr(freqs0$Allele2,3,1000))

#Gen 20
names(freqs20)[4]<- "Allele1"
names(freqs20)[5]<- "Allele2"
freqs20$Allele2[which(substr(freqs20$Allele1,1,1)==2)] <- "2:1.00000" # put this in the spot for the second allele
freqs20$Allele1[which(substr(freqs20$Allele1,1,1)==2)] <- "1:0.00000" 
freqs20$Allele2[which(substr(freqs20$Allele1,3,3)==1)] <- "2:0.00000"  ##
freqs20$Allele1 <- as.numeric(substr(freqs20$Allele1,3,1000))
freqs20$Allele2 <- as.numeric(substr(freqs20$Allele2,3,1000))

#Calculate change
change2<-freqs20$Allele2-freqs0$Allele2

### Setup for function
geno <- gen
phen <- as.matrix(pheno[1:1000,10])
rownames(phen)<-pheno[1:1000,1]
change=change2
perms <- 10000
blockSize=100000/100

# Run function
test<- Ghat_func(geno=geno,phen=phen,change=change2,method = "scale", num_eff = 6667,  perms=1000,plot="Both")
#Pvals.selection[as.numeric(iter)] <- test$p.val

### Save data for plotting
change100 <- change2
effects100 <- test$effects









##################################################################
######################EYE PLOTS - 1000 #############################
##################################################################
dir <- "DataForFigs/1000."
iter <- "001"
source("DataForFigs/Ghat.R")

pheno <- read.table(paste(dir, "SelectPop_data_",iter,".txt",sep=""),header=T,stringsAsFactors=F)
map <- read.table(paste(dir, "lm_mrk_",iter,".txt",sep=""),header=T,stringsAsFactors=F)
geno <- read.table(paste(dir, "SelectPop_mrk_",iter,".txt-head",sep=""),header=F,stringsAsFactors=F,skip=1,sep="", colClasses=c("numeric","character")) # read genos as characters


#load allele frequencies
freqs0<-read.table(paste(dir, "SelectPop_freq_mrk_",iter,".txt",sep=""),header=T,nrows=nrow(map),fill=T,stringsAsFactors=F)
freqs20<-read.table(paste(dir, "SelectPop_freq_mrk_",iter,".txt",sep=""),header=F,skip={4*nrow(map)+1},fill=T,stringsAsFactors=F)

### Manipulate genotypes by coding to -1,0,1 and so that markers are columns, individuals are rows.
gen <- matrix(NA,nrow=nrow(map),ncol=nrow(geno))
gen <- as.data.frame(gen)
names(gen) <- geno[,1]
for(i in 1:1000){
  print(i)
  tmp <- as.numeric(unlist(strsplit(geno[i,2],split="")))
  tmp[which(tmp == 0)] <- -1
  tmp[which(tmp == 3 | tmp ==4)] <- 0
  tmp[which(tmp==2)] <- 1
  gen[,i] <- tmp
}
gen<-t(gen)
gc()

### Calculate allele frequencies
#Gen 0
names(freqs0)[4]<- "Allele1"
names(freqs0)[5]<- "Allele2"
freqs0$Allele2[which(substr(freqs0$Allele1,1,1)==2)] <- "2:1.00000" # put this in the spot for the second allele
freqs0$Allele1[which(substr(freqs0$Allele1,1,1)==2)] <- "1:0.00000" 
freqs0$Allele2[which(substr(freqs0$Allele1,3,3)==1)] <- "2:0.00000" ##
freqs0$Allele1 <- as.numeric(substr(freqs0$Allele1,3,1000))
freqs0$Allele2 <- as.numeric(substr(freqs0$Allele2,3,1000))

#Gen 20
names(freqs20)[4]<- "Allele1"
names(freqs20)[5]<- "Allele2"
freqs20$Allele2[which(substr(freqs20$Allele1,1,1)==2)] <- "2:1.00000" # put this in the spot for the second allele
freqs20$Allele1[which(substr(freqs20$Allele1,1,1)==2)] <- "1:0.00000" 
freqs20$Allele2[which(substr(freqs20$Allele1,3,3)==1)] <- "2:0.00000"  ##
freqs20$Allele1 <- as.numeric(substr(freqs20$Allele1,3,1000))
freqs20$Allele2 <- as.numeric(substr(freqs20$Allele2,3,1000))

#Calculate change
change2<-freqs20$Allele2-freqs0$Allele2

### Setup for function
geno <- gen
phen <- as.matrix(pheno[1:1000,10])
rownames(phen)<-pheno[1:1000,1]
change=change2
perms <- 10000
blockSize=100000/100

# Run function
test<- Ghat_func(geno=geno,phen=phen,change=change2,method = "scale", num_eff = 6667,  perms=1000,plot="Both")
#Pvals.selection[as.numeric(iter)] <- test$p.val

### Save data for plotting
change1000 <- change2
effects1000 <- test$effects



#######################
#######################
### Remove large objects ###
rm(gen)
rm(geno)
gc()
#######################
#######################


#################################################################
######## Plot correlation and test for selection data (10) ###########
#################################################################

## some pretty colors
colfunc <- colorRampPalette(c("lightblue", "darkred","red","yellow"))
#colfunc(n)

### Replace NA with 0 for plotting
change10[which(is.na(change10))]<-0

## compute 2D kernel density, see MASS book, pp. 130-131
library(MASS)
z_10 <- kde2d(change10, effects10, n=200)

plot(change10,effects10,pch=19,ylab="Effect size", xlab="Frequency change",col="#25282E",cex.lab=2,mgp=c(2.5,1,0))
contour(z_10,drawlabels=F,nlevels=200,col=colfunc(200),add=T)
abline(lm(effects10~change10),col="darkred",lwd=3)
#text(-.2,1e-07,paste("r =", round(cor(test10$effects,change3_10,use="pairwise.complete.obs"),3)),cex=1.4)


#################################################################
######## Plot correlation and test for selection data (50) ###########
#################################################################

## some pretty colors
colfunc <- colorRampPalette(c("lightblue", "darkred","red","yellow"))
#colfunc(n)

### Replace NA with 0 for plotting
change50[which(is.na(change50))]<-0
gc()

## compute 2D kernel density, see MASS book, pp. 130-131
library(MASS)
z_50 <- kde2d(change50, effects50, n=200)

plot(change50,effects50,pch=19,ylab="Effect size", xlab="Frequency change",col="#25282E",cex.lab=2,mgp=c(2.5,1,0))
contour(z_50,drawlabels=F,nlevels=200,col=colfunc(200),add=T)
abline(lm(effects50~change50),col="darkred",lwd=3)
#text(-.2,1e-07,paste("r =", round(cor(test10$effects,change3_10,use="pairwise.complete.obs"),3)),cex=1.4)


#################################################################
######## Plot correlation and test for selection data (100) ###########
#################################################################

## some pretty colors
colfunc <- colorRampPalette(c("lightblue", "darkred","red","yellow"))
#colfunc(n)

### Replace NA with 0 for plotting
change100[which(is.na(change100))]<-0

## compute 2D kernel density, see MASS book, pp. 130-131
library(MASS)
z_100 <- kde2d(change100, effects100, n=200)
gc()

plot(change100,effects100,pch=19,ylab="Effect size", xlab="Frequency change",col="#25282E",cex.lab=2,mgp=c(2.5,1,0))
contour(z_100,drawlabels=F,nlevels=200,col=colfunc(200),add=T)
abline(lm(effects100~change100),col="darkred",lwd=3)
#text(-.2,1e-07,paste("r =", round(cor(test10$effects,change3_10,use="pairwise.complete.obs"),3)),cex=1.4)




#################################################################
######## Plot correlation and test for selection data (1000) ###########
#################################################################

## some pretty colors
colfunc <- colorRampPalette(c("lightblue", "darkred","red","yellow"))
#colfunc(n)

### Replace NA with 0 for plotting
change1000[which(is.na(change1000))]<-0

## compute 2D kernel density, see MASS book, pp. 130-131
library(MASS)
z_1000 <- kde2d(change1000, effects1000, n=200)
gc()

plot(change1000,effects1000,pch=19,ylab="Effect size", xlab="Frequency change",col="#25282E",cex.lab=2,mgp=c(2.5,1,0))
contour(z_1000,drawlabels=F,nlevels=200,col=colfunc(200),add=T)
abline(lm(effects1000~change1000),col="darkred",lwd=3)
#text(-.2,1e-07,paste("r =", round(cor(test10$effects,change3_10,use="pairwise.complete.obs"),3)),cex=1.4)


#################################################################################
#################################################################################
#################    Make Total Plot with all 5 panels   ########################
#################################################################################

jpeg("PowerPlot.jpg",height=600,width=800,pointsize=18)
layout(matrix(c(1,1,1,1,2,3,4,5),2,4,byrow=T),heights=c(1.5,1))
par(mar=c(4,4,3,2))

### Power summary
plot(Ghat_TP,pch=19,ylim=c(0,1.1),xaxt="n",cex=0.75,xlab="Number of QTL simulated",ylab="Detection rate",col="brown")
axis(1,c("10 QTL","50 QTL","100 QTL","1000 QTL"),at=c(1,2,3,4))
lines(Ghat_TP,col="brown",lwd=3)
points(Map_tp,col="blue",pch=19,cex=0.75)
lines(Map_tp,col="blue",lwd=3)
segments(1:4,Map_tp-Map_tpsd,1:4,Map_tp+Map_tpsd,col="blue",lwd=3)
segments(1:4,Ghat_TP-Ghat_TPSD,1:4,Ghat_TP+Ghat_TPSD,col="brown",lwd=3)
legend("topleft","c(x,y)",lty=1,lwd=2,pt.cex=0.75,pch=19,col=c("brown","blue"),c("G-hat","Selection Mapping"))

### Eye 10
plot(change10,effects10,pch=19,ylab="Effect size", xlab="Frequency change",col="#25282E",cex.lab=1,mgp=c(2.5,1,0),main="10 QTL")
contour(z_10,drawlabels=F,nlevels=400,col=colfunc(400),add=T)
abline(lm(effects10~change10),col="darkred",lwd=3)

### Eye 50
plot(change50,effects50,pch=19,ylab="Effect size", xlab="Frequency change",col="#25282E",cex.lab=1,mgp=c(2.5,1,0),main="50 QTL")
contour(z_50,drawlabels=F,nlevels=400,col=colfunc(400),add=T)
abline(lm(effects50~change50),col="darkred",lwd=3)

### Eye 100
plot(change100,effects100,pch=19,ylab="Effect size", xlab="Frequency change",col="#25282E",cex.lab=1,mgp=c(2.5,1,0),main="100 QTL")
contour(z_100,drawlabels=F,nlevels=400,col=colfunc(400),add=T)
abline(lm(effects100~change100),col="darkred",lwd=3)

### Eye 1000
plot(change1000,effects1000,pch=19,ylab="Effect size", xlab="Frequency change",col="#25282E",cex.lab=1,mgp=c(2.5,1,0),main="1000 QTL")
contour(z_1000,drawlabels=F,nlevels=400,col=colfunc(400),add=T)
abline(lm(effects1000~change1000),col="darkred",lwd=3)

dev.off()
?layout


###################### Save and Image!! ###########################
save.image("PowerPlot.RData")
