## R script to extract the length of the longest mapped read from each ind
## to the wolbachia genome (ccnice -- AAalbopictus_sequence.fasta)

# vshastry -- Dec 2019

library(foreach)
library(doParallel)

library(stringr)
library(iterators)

#setwd("/Volumes/evolgen/vshastry/wolbachia/analyses")

files<-system("ls ../matches/samfiles/Scaf1260/ | grep 'hitFiles_*'", intern=T)

ids<-str_match(files,"hitFiles_(.+).sam")[,2]

pop<-unique(substr(ids,1,3))
num.pop<-length(unique(substr(ids,1,3)))

system("module load swset gcc samtools")
# running in parallel over 16 threads
cl<-parallel::makeCluster(32)
doParallel::registerDoParallel(cl)
longest.read.melan<-foreach(i=1:length(files), .combine='rbind') %dopar% {
	as.numeric(system(sprintf("samtools view -F 2432 ../matches/samfiles/Super/%s | cut -f 10 | perl -ne 'chomp;print length($_) . \"\n\"' | sort -n | uniq -c | tail -n 1 | awk '{print $1}{print $2}'",files[i]),intern=T))[1:2]
}
# this command gets all the reads sorted by length and number of matches 
scaf.reads<-foreach(i=1:length(files)) %dopar% {
	test<-system(sprintf("samtools view -F 2432 ../matches/samfiles/Scaf1260/%s | cut -f 10 | perl -ne 'chomp;print length($_) . \"\n\"' | sort -n | uniq -c",files[i]),intern=T)
	if(length(test)==0)test<-rep(0,2)
	matrix(unlist(strsplit(trimws(test,which="left")," ")),ncol=2,byrow=T)
}
parallel::stopCluster(cl)

# sum all reads across all inds for one reference genome
tot.melan.sum<-0
for(i in 1:2388){
    tot.melan.sum = tot.melan.sum+sum(as.numeric(melan.reads[i][[1]][,1])*as.numeric(melan.reads[i][[1]][,2]))
}

### repeating the same for Aalbo and Lectu
#which.na<-sapply(aalbo.reads,function(x){sum(is.na(x[[1]]))>0},simplify=T)
#idx.na<-which(which.na==T,arr.ind=T)
#new.aalbo.reads<-replace(aalbo.reads,idx.na,list(cbind(0,0)))
tot.aalbo.sum<-0
for(i in 1:2388){
        tot.aalbo.sum = tot.aalbo.sum+sum(as.numeric(aalbo.reads[i][[1]][,1])*as.numeric(aalbo.reads[i][[1]][,2]))
}
tot.lectu.sum<-0
for(i in 1:2388){
        tot.lectu.sum = tot.lectu.sum+sum(as.numeric(super.reads[i][[1]][,1])*as.numeric(super.reads[i][[1]][,2]))
}

## get number of reads mapping to each ind
tot.reads<-rep(0,2388)
for(i in 1:2388){
    temp<-matrix(as.numeric(melan.reads[[i]]),ncol=2)
    tot.reads[i]<-sum(temp[which(temp[,2]>80),1])
}

## unlist reads to put elements into a matrix (contains reads from all inds)
melan.reads<-do.call("rbind",sapply(melan.reads,function(x){matrix(as.numeric(x),ncol=2)},simplify=T))
scaf.reads<-do.call("rbind",sapply(scaf.reads,function(x){matrix(as.numeric(x),ncol=2)},simplify=T))
super.reads<-do.call("rbind",sapply(super.reads,function(x){matrix(as.numeric(x),ncol=2)},simplify=T))
aalbo.reads<-do.call("rbind",sapply(aalbo.reads,function(x){matrix(as.numeric(x),ncol=2)},simplify=T))
lectu.reads<-do.call("rbind",sapply(lectu.reads,function(x){matrix(as.numeric(x),ncol=2)},simplify=T))
table.melan<-rep(NA,max(melan.reads[,2]))
for(i in 1:max(melan.reads[,2])){
    table.melan[i]<-sum(melan.reads[which(melan.reads[,2]==i),1])
}
table.scaf<-rep(NA,max(scaf.reads[,2]))
for(i in 1:max(scaf.reads[,2])){
        table.scaf[i]<-sum(scaf.reads[which(scaf.reads[,2]==i),1])
}
table.super<-rep(NA,max(super.reads[,2]))
for(i in 1:max(super.reads[,2])){
        table.super[i]<-sum(super.reads[which(super.reads[,2]==i),1])
}
table.aalbo<-rep(NA,max(aalbo.reads[,2]))
for(i in 1:max(aalbo.reads[,2])){
        table.aalbo[i]<-sum(aalbo.reads[which(aalbo.reads[,2]==i),1])
}
table.lectu<-rep(NA,max(lectu.reads[,2]))
for(i in 1:max(lectu.reads[,2])){
        table.lectu[i]<-sum(lectu.reads[which(lectu.reads[,2]==i),1])
}

#get population-level infection statistics
longest.read.melan[is.na(longest.read.melan)]<-0
# set 80bp threshold to count as an infection match
infected<-rep(F,length(ids))
# array of whether an individual is infected or not (across all inds)
infected<-longest.read.melan[,2]>86 & longest.read.melan[,1]>1

# comparing results with PCR for different numbers of mapped reads
pcr.res<-read.table("pcr_results.txt",header=T)
metrics<-matrix(NA,nrow=3,ncol=25)
rownames(metrics)<-c("accuracy","TPR","FPR") # comparing bioinformatics with PCR (truth)
for(i in 1:25){
    infected<-longest.read.melan[,2]>=85 & longest.read.melan[,1]>=i
    metrics[1,i]<-(table(infected,pcr.res$pcr)[1,1]+table(infected,pcr.res$pcr)[2,2])*100/129
    metrics[2,i]<-table(infected,pcr.res$pcr)[1,1]*100/(table(infected,pcr.res$pcr)[1,1]+table(infected,pcr.res$pcr)[1,2])
    metrics[3,i]<-table(infected,pcr.res$pcr)[2,1]*100/(table(infected,pcr.res$pcr)[2,2]+table(infected,pcr.res$pcr)[2,1])
}
write.table(metrics,file="metrics_PCR85bp.txt",quote=F)

metrics<-matrix(NA,nrow=3,ncol=6)
rownames(metrics)<-c("accuracy","TPR","FPR") # comparing bioinformatics with PCR (truth)
for(i in 1:ncol(metrics)){
    infected<-tot.reads>=5^(i-1)
    metrics[1,i]<-(table(infected,pcr.res$pcr)[1,1]+table(infected,pcr.res$pcr)[2,2])*100/129
    metrics[2,i]<-table(infected,pcr.res$pcr)[1,1]*100/(table(infected,pcr.res$pcr)[1,1]+table(infected,pcr.res$pcr)[1,2])
    metrics[3,i]<-table(infected,pcr.res$pcr)[2,1]*100/(table(infected,pcr.res$pcr)[2,2]+table(infected,pcr.res$pcr)[2,1])
}
write.csv(metrics,file="metrics_PCRtot.txt",quote=F)

# contains MOST RECENT infection status -- May 27, 2021
infected.10g80<-tot.reads>=10
infected.20g80<-tot.reads>=20

# count number of infected inds for length and number of reads
num.inf<-seq(1,1500,1)
len.inf<-85:87
infected.mat<-matrix(0,nrow=length(num.inf),ncol=4)
#row.names(infected.mat)<-paste0("X",num.inf)
#col.names(infected.mat)<-paste0("X",len.inf)
for(i in 1:nrow(infected.mat)){
    for(j in 1:ncol(infected.mat)){
        infected.mat[i,j]<-sum(longest.read.melan[,1]>=num.inf[i] & longest.read.melan[,2]>=len.inf[j])
    }
}
## plotting a single line for each length of reads
#pdf("infpcbylength.pdf")
plot(num.inf[1:100],infected.mat[1:100,1]/2388,type="l",col="darkgoldenrod1",
     xlab="Number of reads",ylab="Percent infected",lwd=2,ylim=c(0,1.1))
lines(num.inf[1:100],infected.mat[1:100,2]/2388,col="darkgoldenrod2",lwd=2)
lines(num.inf[1:100],infected.mat[1:100,3]/2388,col="darkgoldenrod3",lwd=2)
legend("topright",col=c(paste0("darkgoldenrod",1:3)),legend=len.inf,lwd=2)
par(new=T, fig=c(0.25,0.80,0.45,1.0))
plot(num.inf,infected.mat[,1]/2388,type="l",col="darkgoldenrod1",lwd=2,
     xlab='',ylab='',ylim=c(0,1))
lines(num.inf,infected.mat[,2]/2388,col="darkgoldenrod2",lwd=2)
lines(num.inf,infected.mat[,3]/2388,col="darkgoldenrod3",lwd=2)
#heatmap(infected.mat,Rowv=NA,Colv=NA,labRow=paste0(num.inf),
#        labCol=paste0(len.inf),xlab="Length of reads", ylab="Number of reads",
#        revC=T,col=rev(heat.colors(n=50)))

count.infected.pop<-rep(0,length(num.pop))
count.total.pop<-table(as.factor(substr(ids,1,3)))
for(i in 1:num.pop){
	count.infected.pop[i]<-sum(infected[substr(ids,1,3)==pop[i]])
}

write.table(cbind(names(count.total.pop),count.infected.pop,count.total.pop),"infectedmelan.csv",
           quote=F,row.names=F,col.names=c("code","inf","total"))

hist(count.infected.pop/count.total.pop,breaks=seq(0,1,0.05))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# get infection status by ind instead of pop
max.reads<-matrix(0,nrow=length(files),ncol=2)
for(i in 1:length(files)){
	max.reads[i,]<-as.numeric(tail(melan.reads[i][[1]],n=1))
}

max.reads<-longest.read.melan

#no.matches<-which(is.na(max.reads[,2]))

infection<-data.frame(matrix(ncol=8,nrow=length(files)))
# scaf1260 numbers
#colnames(infection)<-c("Ind","NoMatch","1of81","594of81","1of87","6of86")
# super numbers
colnames(infection)<-c("Ind","NoMatch","1of81","524of81","1of87","7of87","10of>80","20of>80")


# aalbo numbers
#colnames(infection)<-c("Ind","NoMatch","1of80","238of80","1of86","27of86")

infection$Ind<-ids

infected<-rep(F, length(files))
infected[which(longest.read.melan[,2]==0)]<-T
infection$NoMatch<-infected

infected<-rep(F, length(files))
infected[max.reads[,1]>0 & max.reads[,2]>80]<-T
infection[,3]<-infected

infected<-rep(F, length(files))
infected[max.reads[,1]>237 & max.reads[,2]>80]<-T
infection[,4]<-infected

infected<-rep(F, length(files))
infected[max.reads[,1]>0 & max.reads[,2]>86]<-T
infection[,5]<-infected

infected<-rep(F, length(files))
infected[max.reads[,1]>27 & max.reads[,2]>86]<-T
infection[,6]<-infected

infection[,7]<-infected.10g80
infection[,8]<-infected.20g80

write.csv(rbind(infection,c("total",apply(infection[-1,2:8],2,sum))), file="wolinfectedbyindsuper.csv", quote=F)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# count the number of reads that map to Wolbachia over total number of reads
# average and sd coverage
cl<-parallel::makeCluster(32)
doParallel::registerDoParallel(cl)
mean.sd.cov<-foreach(i=1:length(files)) %dopar% {
    as.numeric(system(sprintf("samtools depth ../matches/bamfiles/Super/sortedbamfiles/%s | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR; print sqrt(sumsq/NR - (sum/NR)**2)}' | head -n 2",gsub('.sam','.sorted.bam',files[i])),intern=T))
}
parallel::stopCluster(cl)

#summary(unlist(mean.sd.cov)[seq(1,2388,2)])
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#  0.9767   3.2208   8.5648  22.2329  22.1446 307.6790 
#summary(unlist(mean.sd.cov)[seq(2,2388,2)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   5.095  19.452  52.371  61.840 854.252 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# obtaining infection status for each pop and plotting it 

gbs.list<-read.csv("GBS_list_2019.csv",stringsAsFactors=F)

# extracting pop codes from csv file
inf.ind<-read.csv("wolinfectedbyindmelan.csv")
inf.ind<-inf.ind[-nrow(inf.ind),]

# create a column of code names
inf.ind$code<-gsub('[0-9]+|_','',tolower(substr(inf.ind$Ind,1,4)))

# get unique pop names (in the same format as GBS list)
inf.pop<-gsub('[0-9]+|_','',tolower(unique(substr(inf.ind$Ind,1,4))))

inf.counts<-rep(0,length(inf.pop))
for(idx in 1:length(inf.pop)){
    inf.counts[idx]<-table(inf.ind$X1of80[inf.ind$code==inf.pop[idx]])[3]
}

