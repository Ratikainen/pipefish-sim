######################################################################
##          Simulation of selection index with random mating        ##
######################################################################

#package needed for standard error
library(plotrix)

#Make sure sampling function behaves the way it should
resample <- function(x, ...) x[sample.int(length(x), ...)]

#Import data
Matings<-read.table("mating.csv", dec=".", sep=";", header=T)
Fecundity<-read.table("Matingsystemdata1.csv", dec=".", sep=";", header=T)

Female_fecundity<-Fecundity[Fecundity$sex==1&Fecundity$osr==2,]$embryo
Male_fecundity<-Fecundity[Fecundity$sex==2&Fecundity$full==100,]$embryo
Male_fecundity<-Male_fecundity[1:155]
egg_vector<-Matings[Matings$OSR==2,]$totalt                       #NEEDS to be changed to the correct sex-ratio male-biased:2 female-biased:1
egg_vector<-sort(egg_vector)

#For male-biased
egg_vector<-egg_vector[c(65:126,166:length(egg_vector))]          #We remove some data with few eggs transferred, assuming that at least some of these are due to the male being full or the female going empty
#egg_vector<-egg_vector[c(16:74,100:length(egg_vector))] 

#For female-biased 
egg_vector<-egg_vector[c(29:89,119:length(egg_vector))]  

#A quick check of the data
summary(tapply(Matings$totalt, Matings$hunn, sum, na.rm=T))
summary(tapply(Matings$totalt, list(Matings$hann,Matings$treatment), sum, na.rm=T))
summary(egg_vector)
plot(Fecundity[Fecundity$sex==1&Fecundity$osr==2,]$lengde,Fecundity[Fecundity$sex==1&Fecundity$osr==2,]$embryo)
plot(Fecundity[Fecundity$sex==2&Fecundity$full==100,]$lengde,Fecundity[Fecundity$sex==2&Fecundity$full==100,]$embryo)
hist(Fecundity[Fecundity$sex==1&Fecundity$osr==2,]$embryo)

#Set parameters for the model
Nf<-10 #number of females
#mean.eggs.females <-84  #Mean number of eggs in a female 
#mean.eggs.trans<-70 #Mean number of eggs transferred in one mating

Nm<-30-Nf #number of males
#mean.space.males <-72  #Mean space in a male
#max.eggs.trans<-80  # Maximum number of eggs that can be transferred during one mating

X<-100 #number of times to repeat simulation


#These are just matrices for storing data, needs to be reset before running a simulation
females<-rep(1:Nf,3)
dim(females)<-c(Nf,3)
males<-rep(1:Nm,3)
dim(males)<-c(Nm,3)
matings<-rep(0,Nm*Nf*X)
dim(matings)<-c(X,Nm,Nf)
repr.success.f<-rep(0,X*Nf)
dim(repr.success.f)<-c(X,Nf)
repr.success.m<-rep(0,X*Nm)
dim(repr.success.m)<-c(X,Nm)
tot.matings.m<-rep(0,X*Nm)
dim(tot.matings.m)<-c(X,Nm)
tot.matings.f<-rep(0,X*Nf)
dim(tot.matings.f)<-c(X,Nf)
OSS<-rep(0,X*2)
dim(OSS)<-c(X,2)
OS<-rep(0,X*2)
dim(OS)<-c(X,2)

egg_vector2<-NA


for (x in 1:X) {        #We do the whole simulation X times
  females[,2]<-sample(Female_fecundity,size=Nf)           #pick random number of eggs 
  females[,3]<-females[,2]                                #copy so we remember
  males[,2]<-sample(Male_fecundity,size=Nm)               #pick random space for eggs 
  males[,3]<-males[,2]                                    #copy so we remember
  while (sum(females[,2])>0&&sum(males[,2])>0){           #While there are still females with eggs and males with space
    females.with.eggs<-resample(females[females[,2]>0,1]) #randomise females
    fem<-females.with.eggs[1]                             #pick a random female
    males.with.space<-resample(males[males[,2]>0,1])      #randomise males
    male<-males.with.space[1]                             #pick random male
    max.eggs<-min(females[fem,2],males[male,2])           #maximum number of eggs that can be transferred is this (lower if female have few eggs or male has little space
    if (max.eggs<4){                                      #if max number of eggs that can be transferred is below 5 then that number will be transferred
      eggs.trans<-max.eggs
      }else {
      eggs.trans<-min(max.eggs,sample(egg_vector,size=1)) # pick a random number of eggs to transfer 
      }
    tranferred.eggs<-eggs.trans[1]                        #pick random number of eggs to transfer
    egg_vector2<-c(egg_vector2,tranferred.eggs)
    if(tranferred.eggs==0){break}
    matings[x,male,fem]<-1                                #register that these individuals have mated
    females[fem,2]<-females[fem,2]-tranferred.eggs        #calculate and store eggs left in female
    males[male,2]<-males[male,2]-tranferred.eggs          #calculate and store space left in male
  }
  repr.success.f[x,]<-females[,3]-females[,2]             #reproductive success=starting number of eggs - number of eggs left
  repr.success.m[x,]<-males[,3]-males[,2]                 #reproductive success=starting space - space left
  for (m in 1:Nm){                                        #for all males
    tot.matings.m[x,m]<-sum(matings[x,m,])                #calculate total number of partners
    }
  for (f in 1:Nf){                                        #for all females
    tot.matings.f[x,f]<-sum(matings[x,,f])                #calculate total number of partners
    }
  OSS[x,1]<-var(tot.matings.m[x,])/(mean(tot.matings.m[x,]))^2   #Calculate opportunity for sexual selection males  
  OSS[x,2]<-var(tot.matings.f[x,])/(mean(tot.matings.f[x,]))^2   #Calculate opportunity for sexual selection females 
  OS[x,1]<-var(repr.success.m[x,])/(mean(repr.success.m[x,]))^2  #Calculate opportunity for selection males 
  OS[x,2]<-var(repr.success.f[x,])/(mean(repr.success.f[x,]))^2  #Calculate opportunity for selection females
  }


#############################################################################################################################
## RESULTS

#Create storage matrix for final results 
 results<-rep(0,8*2)
 dim(results)<-c(2,8)
 colnames(results)<-c('mean_repr','var_repr','mean_mat','var_mat','I','I_SE','Is','Is_SE')
 rownames(results)<-c('females','males')

## Calculate and insert results

# Reproductive success
results[1,1]<-mean(repr.success.f)  
results[1,2]<-mean((var(repr.success.f))^2)    
results[2,1]<-mean(repr.success.m)  
results[2,2]<-mean((var(repr.success.m))^2)

# Matings
results[1,3]<-mean(tot.matings.f)  
results[1,4]<-mean((var(tot.matings.f))^2)             
results[2,3]<-mean(tot.matings.m)  
results[2,4]<-mean((var(tot.matings.m))^2) 
    
#Mean opportunity for sexual selection males   
results[2,7]<-mean(OSS[,1])  
results[2,8]<-std.error(OSS[,1])  

#Mean opportunity for sexual selection females   
results[1,7]<-mean(OSS[,2])  
results[1,8]<-std.error(OSS[,2])
    
#Mean opportunity for selection males   
results[2,5]<-mean(OS[,1])  
results[2,6]<-std.error(OS[,1])  

#Mean opportunity for selection females   
results[1,5]<-mean(OS[,2])  
results[1,6]<-std.error(OS[,2])
    
results #Print table with results
