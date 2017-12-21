###########################################################################
### This is a function to test for selection on a trait under polygenic ###
### control.                                                            ### 
###########################################################################


### Timothy M Beissinger
### 5/1/2017
###########################################################################
### This version of the function applies to raw genotype/phenotype data ###
###########################################################################
# Inputs are listed below
#geno:   Matrix of genotype data. Individuals in rows, genotypes (-1, 0, 1) in columns.
#phen:   Matrix of phenotype data. Individuals in rows, phenotype in single column.
#change: Change in allele frequency. Must record the change at the allele corresponding to "1" in the genotype file.
#method: "trim", "scale", or "vanilla". How should LD be accounted for?
#perms:  Number of permutations to run
#plot:   "Ghat", "Cor", or "Both", Should a plot of the Ghat or correlation test be returned?
#blockSize: How large should blocks for trimming be? Only required if method = "trim"
#num_eff: What is the effective number of independent markers?

Ghat_func <- function(geno=geno,phen=phen,change=change2, method = "scale", perms=1000 ,plot="Both",blockSize=1000, num_eff = NULL){
  require(rrBLUP) 
  ##################
  ##################
  ### Some empty variables###
  sd_eff <- NA
  ###########################
    if(method == "vanilla"){
    blupModel <- mixed.solve(phen, Z=geno,K=NULL,SE=F,return.Hinv=FALSE,method="ML")  
    effects <- blupModel$u
    effects.mat <- as.matrix(effects)
    
    # Compute test statistics
    Ghat <- sum(change*effects,na.rm=T)
    Cor <- cor(change,effects,use="pairwise.complete.obs")
    
    # Do shuffling permutation test
    Ghat_perm <- c()
    Cor_perm <- c()
    print("Doing Permutation Test")
    for(test in 1:perms){
      Ghat_perm[test] <- sum(change*effects[sample(length(effects),replace=F)],na.rm=T)
      Cor_perm[test] <- cor(change,effects[sample(length(effects),replace=F)],use="pairwise.complete.obs")
    }
    
    
    # Make plot
    if(plot=="Ghat"){
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution",freq=F)
      abline(v=Ghat,lwd=3,col="darkblue")
      #x<-seq(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm)),length.out=1000)  #these two lines for plotting only
      #lines(x,dnorm(x,mean=mean(Ghat_perm),sd=sd(Ghat_perm)),col="red",lwd=3)
    }
    if(plot=="Cor"){
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    if(plot=="Both"){
      par(mfrow=c(1,2))
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue")
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    
    # calculate p-value
    if(Ghat > 0)  p.val <- length(which(Ghat_perm >= Ghat | Ghat_perm <= -1*Ghat))/length(Ghat_perm)
    if(Ghat < 0)  p.val <- length(which(Ghat_perm <= Ghat | Ghat_perm >= -1*Ghat))/length(Ghat_perm)
  }
  ######################
  ######################
  if(method == "trim"){
    
    ########Trim###########
    start <- sample(length(change)/blockSize,1)
    index<-seq(start,length(change),by=blockSize)
    change<-change[index]
    ########END Trim########
    
    blupModel <- mixed.solve(phen, Z=geno[,index],K=NULL,SE=F,return.Hinv=FALSE,method="ML")
    effects <- blupModel$u
    effects.mat <- as.matrix(effects)
    
    # Compute test statistics
    Ghat <- sum(change*effects,na.rm=T)
    Cor <- cor(change,effects,use="pairwise.complete.obs")
    
    # Do shuffling permutation test
    Ghat_perm <- c()
    Cor_perm <- c()
    print("Doing Permutation Test")
    for(test in 1:perms){
      Ghat_perm[test] <- sum(change*effects[sample(length(effects),replace=F)],na.rm=T)
      Cor_perm[test] <- cor(change,effects[sample(length(effects),replace=F)],use="pairwise.complete.obs")
    }
    
    # Make plot
    if(plot=="Ghat"){
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue")
    }
    if(plot=="Cor"){
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    if(plot=="Both"){
      par(mfrow=c(1,2))
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue")
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    
    # calculate p-value
    if(Ghat > 0)  p.val <- length(which(Ghat_perm >= Ghat | Ghat_perm <= -1*Ghat))/length(Ghat_perm)
    if(Ghat < 0)  p.val <- length(which(Ghat_perm <= Ghat | Ghat_perm >= -1*Ghat))/length(Ghat_perm)
  }
  ########################
  ########################
  if(method == "scale"){
    blupModel <- mixed.solve(phen, Z=geno,K=NULL,SE=F,return.Hinv=FALSE,method="ML")  
    effects <- blupModel$u
    effects.mat <- as.matrix(effects)
    
    # Compute test statistics
    Ghat <- sum(change*effects,na.rm=T)
    Cor <- cor(change,effects,use="pairwise.complete.obs")
    
    # Do shuffling permutation test
    Ghat_perm <- c()
    Cor_perm <- c()
    print("Doing Permutation Test")
    for(test in 1:perms){
      Ghat_perm[test] <- sum(change*effects[sample(length(effects),replace=F)],na.rm=T)
      Cor_perm[test] <- cor(change,effects[sample(length(effects),replace=F)],use="pairwise.complete.obs")
    }
    
    # Do scaling (Ghat)
    scale = sqrt(length(effects))/sqrt(num_eff)
    sd_eff = sd(Ghat_perm)*scale
    
    # Do scaling (Cor)
    sd_eff_cor = sd(Cor_perm)*scale
    
    # Make plot
    if(plot=="Ghat"){
      left <- min(mean(Ghat_perm)-4*sd_eff,Ghat)
      right <- max(mean(Ghat_perm)+4*sd_eff,Ghat)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Ghat_perm),sd=sd_eff),col="red",lwd=3,type="l",ylab="Density",xlab="Null Distribution",main="Scaled Ghat test for selection")
      abline(v=Ghat,lwd=3,col="darkblue")
    }
    
    if(plot=="Cor"){
      left <- min(mean(Cor_perm)-4*sd_eff_cor,Ghat)
      right <- max(mean(Cor_perm)+4*sd_eff_cor,Ghat)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Cor_perm),sd=sd_eff_cor),col="red",lwd=3,type="l",ylab="Density",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    
    if(plot=="Both"){
      par(mfrow=c(1,2))
      left <- min(mean(Ghat_perm)-4*sd_eff,Ghat)
      right <- max(mean(Ghat_perm)+4*sd_eff,Ghat)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Ghat_perm),sd=sd_eff),col="red",lwd=3,type="l",ylab="Density",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue") 
      left <- min(mean(Cor_perm)-4*sd_eff_cor,Cor)
      right <- max(mean(Cor_perm)+4*sd_eff_cor,Cor)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Cor_perm),sd=sd_eff_cor),col="red",lwd=3,type="l",ylab="Density",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    
    # calculate p-value
    if(Ghat > 0)  p.val <- pnorm(Ghat,mean=mean(Ghat_perm),sd=sd_eff,lower.tail=F)*2
    if(Ghat < 0)  p.val <- pnorm(Ghat,mean=mean(Ghat_perm),sd=sd_eff,lower.tail=T)*2
  }
  
  # return values
  return(list(p.val=p.val, Ghat=Ghat, Ghat_perm=Ghat_perm, sd_eff=sd_eff, Cor=Cor, Cor_perm=Cor_perm, effects=effects))
}



###############################################################################
### This version of the function applies with externally estimated effects. ###
###############################################################################
###########################################################################
### This version of the function applies to raw genotype/phenotype data ###
###########################################################################
# Inputs are listed below
#effects: Vector of effect sizes
#change: Change in allele frequency. Must record the change at the allele corresponding to "1" in the genotype file.
#method: "trim", "scale", or "vanilla". How should LD be accounted for?
#perms:  Number of permutations to run
#plot:   "Ghat", "Cor", or "Both", Should a plot of the Ghat or correlation test be returned?
#blockSize: How large should blocks for trimming be? Only required if method = "trim"
#num_eff: What is the effective number of independent markers?

Ghat_func_effectsKnown <- function(effects=effects,change=change2, method = "scale", perms=1000 ,plot="Both",blockSize=1000, num_eff = NULL){
  ### Some empty variables###
  sd_eff <- NA
  ###########################
  
  ##################
  ##################
  if(method == "vanilla"){
    effects.mat <- as.matrix(effects)
    
    # Compute test statistics
    Ghat <- sum(change*effects,na.rm=T)
    Cor <- cor(change,effects,use="pairwise.complete.obs")
    
    # Do shuffling permutation test
    Ghat_perm <- c()
    Cor_perm <- c()
    print("Doing Permutation Test")
    for(test in 1:perms){
      Ghat_perm[test] <- sum(change*effects[sample(length(effects),replace=F)],na.rm=T)
      Cor_perm[test] <- cor(change,effects[sample(length(effects),replace=F)],use="pairwise.complete.obs")
    }
    
    
    # Make plot
    if(plot=="Ghat"){
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution",freq=F)
      abline(v=Ghat,lwd=3,col="darkblue")
      #x<-seq(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm)),length.out=1000)  #these two lines for plotting only
      #lines(x,dnorm(x,mean=mean(Ghat_perm),sd=sd(Ghat_perm)),col="red",lwd=3)
    }
    if(plot=="Cor"){
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    if(plot=="Both"){
      par(mfrow=c(1,2))
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue")
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    
    # calculate p-value
    if(Ghat > 0)  p.val <- length(which(Ghat_perm >= Ghat | Ghat_perm <= -1*Ghat))/length(Ghat_perm)
    if(Ghat < 0)  p.val <- length(which(Ghat_perm <= Ghat | Ghat_perm >= -1*Ghat))/length(Ghat_perm)
  }
  ######################
  ######################
  if(method == "trim"){
    
    ########Trim###########
    start <- sample(length(change)/blockSize,1)
    index<-seq(start,length(change),by=blockSize)
    change<-change[index]
    ########END Trim########
    
    effects.mat <- as.matrix(effects)
    
    # Compute test statistics
    Ghat <- sum(change*effects,na.rm=T)
    Cor <- cor(change,effects,use="pairwise.complete.obs")
    
    # Do shuffling permutation test
    Ghat_perm <- c()
    Cor_perm <- c()
    print("Doing Permutation Test")
    for(test in 1:perms){
      Ghat_perm[test] <- sum(change*effects[sample(length(effects),replace=F)],na.rm=T)
      Cor_perm[test] <- cor(change,effects[sample(length(effects),replace=F)],use="pairwise.complete.obs")
    }
    
    # Make plot
    if(plot=="Ghat"){
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue")
    }
    if(plot=="Cor"){
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    if(plot=="Both"){
      par(mfrow=c(1,2))
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue")
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    
    # calculate p-value
    if(Ghat > 0)  p.val <- length(which(Ghat_perm >= Ghat | Ghat_perm <= -1*Ghat))/length(Ghat_perm)
    if(Ghat < 0)  p.val <- length(which(Ghat_perm <= Ghat | Ghat_perm >= -1*Ghat))/length(Ghat_perm)
  }
  ########################
  ########################
  if(method == "scale"){
    effects.mat <- as.matrix(effects)
    
    # Compute test statistics
    Ghat <- sum(change*effects,na.rm=T)
    Cor <- cor(change,effects,use="pairwise.complete.obs")
    
    # Do shuffling permutation test
    Ghat_perm <- c()
    Cor_perm <- c()
    print("Doing Permutation Test")
    for(test in 1:perms){
      Ghat_perm[test] <- sum(change*effects[sample(length(effects),replace=F)],na.rm=T)
      Cor_perm[test] <- cor(change,effects[sample(length(effects),replace=F)],use="pairwise.complete.obs")
    }
    
    # Do scaling (Ghat)
    scale = sqrt(length(effects))/sqrt(num_eff)
    sd_eff = sd(Ghat_perm)*scale
    
    # Do scaling (Cor)
    sd_eff_cor = sd(Cor_perm)*scale
    
    # Make plot
    if(plot=="Ghat"){
      left <- min(mean(Ghat_perm)-4*sd_eff,Ghat)
      right <- max(mean(Ghat_perm)+4*sd_eff,Ghat)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Ghat_perm),sd=sd_eff),col="red",lwd=3,type="l",ylab="Density",xlab="Null Distribution",main="Scaled Ghat test for selection")
      abline(v=Ghat,lwd=3,col="darkblue")
    }
    
    if(plot=="Cor"){
      left <- min(mean(Cor_perm)-4*sd_eff_cor,Ghat)
      right <- max(mean(Cor_perm)+4*sd_eff_cor,Ghat)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Cor_perm),sd=sd_eff_cor),col="red",lwd=3,type="l",ylab="Density",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    
    if(plot=="Both"){
      par(mfrow=c(1,2))
      left <- min(mean(Ghat_perm)-4*sd_eff,Ghat)
      right <- max(mean(Ghat_perm)+4*sd_eff,Ghat)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Ghat_perm),sd=sd_eff),col="red",lwd=3,type="l",ylab="Density",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue") 
      left <- min(mean(Cor_perm)-4*sd_eff_cor,Cor)
      right <- max(mean(Cor_perm)+4*sd_eff_cor,Cor)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Cor_perm),sd=sd_eff_cor),col="red",lwd=3,type="l",ylab="Density",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    
    # calculate p-value
    if(Ghat > mean(Ghat_perm))  p.val <- pnorm(Ghat,mean=mean(Ghat_perm),sd=sd_eff,lower.tail=F)*2
    if(Ghat < mean(Ghat_perm))  p.val <- pnorm(Ghat,mean=mean(Ghat_perm),sd=sd_eff,lower.tail=T)*2
  }
  
  # return values
  return(list(p.val=p.val, Ghat=Ghat, Ghat_perm=Ghat_perm, sd_eff=sd_eff, Cor=Cor, Cor_perm=Cor_perm, effects=effects))
}