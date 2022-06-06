## R script to analyze haplotypes from fasta file 
## vshastry -- Sep 29, 2020

library(ape)
library(pegas)

naso<-read.dna("filt75pc.super.fasta",format="fasta")
## need to remove individuals that are uninfected
inf.ind<-read.csv("../wolinfectedbyindsuper.csv",stringsAsFactors=F)
inf.ind<-inf.ind[-nrow(inf.ind),]

naso<-naso[-which(inf.ind[,4]=="FALSE"),]
# 20x of >80bp
naso<-naso[-which(inf.ind[,9]=="FALSE"),]

miss.pc<-apply(as.character(naso),1,function(x){sum(x=="n")/ncol(naso)})

## finding nucleotide diversity (Nei 1987) for each pop separately with variance
pop2specs<-read.table("../../../lyc2019analysis/vcffiles/pop2specsfull.txt",sep=",")
popnames<-pop2specs[-which(inf.ind[,4]=="FALSE"),][miss.pc==0,1]
start<-rep(NA,length(unique(popnames)))
end<-rep(NA,length(unique(popnames)))
nucdiv<-data.frame("code"=rep(NA,length(unique(popnames))),"mean"=rep(NA,length(unique(popnames))),"var"=rep(NA,length(unique(popnames))))
for(i in 1:length(unique(popnames))){
    start[i]<-min(which(popnames==unique(popnames)[i]))
    end[i]<-max(which(popnames==unique(popnames)[i]))

    nucdiv[i,1]<-as.character(unique(popnames)[i])
    nucdiv[i,2]<-nuc.div(naso[start[i]:end[i],],var=T)[1]
    nucdiv[i,3]<-nuc.div(naso[start[i]:end[i],],var=T)[2]
}

write.csv(file="nucdivhaplopop.txt",x=nucdiv,quote=F,row.names=F)

nasohaps<-haplotype(naso)

# atleast 10 inds with maximum N pc of 25 (=17 missing bp out of 68)
summary(subset(nasohaps,minfreq=10,maxna=0.25))
red.nasohaps<-subset(nasohaps,minfreq=10,maxna=0.25)

# II    IV  LIII   LIV LXVII  XCII   CXI 
# 30  1127    11    31    39    23   207 

## strict=F (default)
# II     IV      V  XVIII LXXIII  LXXIX    CXV    CXX  CXCIV 
# 34    440     13    739     59     23    216     15     12

## strict=T
# VI        XX      CIII     CVIII     CXIII      CLIV   CLXXVII    CLXXXV 
# 261       596        59        23        39        29       157        15 
# CCCXXVIII 
#       12

#red.nasonet<-haploNet(red.nasohaps)
#plot(red.nasonet,size=attr(red.nasonet,"freq")/5,fast=F)

## include species nominal labels and plot the pie charts
pop2specs<-read.table("../../../lyc2019analysis/vcffiles/pop2specsfull.txt",sep=",")
rownames(naso)<-pop2specs$V1[-which(inf.ind[,4]=="FALSE")]
ind.hap.pop<-with(stack(setNames(attr(nasohaps,"index"),rownames(nasohaps))),table(hap=ind,inds=rownames(naso)[values]))
ind.hap.pop.all<-with(stack(setNames(attr(nasohaps,"index"),rownames(nasohaps))),table(hap=ind,inds=inf.ind$Ind[-which(inf.ind[,4]=="FALSE")][values]))

rownames(naso)<-pop2specs$V2[-which(inf.ind[,4]=="FALSE")]
nasospecshaps<-nasohaps #haplotype(naso)
spec.hap.all<-with(stack(setNames(attr(nasospecshaps,"index"),rownames(nasospecshaps))),table(hap=ind,specs=rownames(naso)[values]))

red.nasospechaps<-subset(nasospecshaps,minfreq=10,maxna=0.25)
spec.hap<-with(stack(setNames(attr(red.nasospechaps,"index"),rownames(red.nasospechaps))),table(hap=ind,specs=rownames(naso)[values]))

## get number of differences in haplotype sequences
diff.haplo<-matrix(0,9,9)
for(i in 1:9){
    for(j in i:9){
        diff.haplo[i,j]<-length(diffHaplo(red.nasohaps,i,j)$pos)
        diff.haplo[j,i]<-diff.haplo[i,j]
    }
}
rownames(diff.haplo)<-dimnames(red.nasohaps)[[1]]
colnames(diff.haplo)<-dimnames(red.nasohaps)[[1]]

#diff.super<-matrix(-1,3,511)
#rownames(diff.super)<-c("IV","CXI","LIII")
#colnames(diff.super)<-rownames(nasohaps)
#for(i in c("IV","CXI","LIII")){
#     for(j in rownames(nasohaps)){
#         diff.super[i,j]<-length(diffHaplo(nasohaps,i,j)$pos)
#     }
#}

diff.super<-matrix(-1,9,dim(nasohaps)[1])
rownames(diff.super)<- rownames(diff.haplo)
colnames(diff.super)<-rownames(nasohaps)
for(i in rownames(diff.super)){
    for(j in rownames(nasohaps)){
         diff.super[i,j]<-length(diffHaplo(nasohaps,i,j)$pos)
     }
}

dist.dna.super<-matrix(-1,9,dim(nasohaps)[1])
rownames(dist.dna.super)<- rownames(diff.haplo)
colnames(dist.dna.super)<-rownames(nasohaps)
dist.dna.super<-dist.dna(nasohaps,"JC69",p=T,as.matrix=T)[rownames(diff.haplo),]

# get unique haps based on integer nucleotide differences
dist.dna.N<-dist.dna(nasohaps,"N",p=T,as.matrix=T)
rownames(nasohaps)[(dist.dna.N["I",]==0 | dist.dna.N["III",]==0 | dist.dna.N["XVIII",]==0 | dist.dna.N["CCCLXXIII",]==0 | dist.dna.N["CCLXII",]==0 | dist.dna.N["XXVII",]==0 | dist.dna.N["LXX",]==0 | dist.dna.N["XLVIII",]==0)] ## 378 out of 465 haps


## concatenating all the genomes into the top 3 using differences
#hap1<-c("II","IV","V","XVIII","CXCIV")
hap1<-rownames(diff.haplo)[c(1,2,3,4,9)]
#hap2<-c("LXXIII","LXXIX","CXX")
hap2<-rownames(diff.haplo)[c(5,6,8)]
#hap3<-"CXV"
hap3<-rownames(diff.haplo)[7]

#sg1<-unique(append(names(which(diff.super[1,]<4)),c(names(which(diff.super[2,]<4)),names(which(diff.super[3,]<3)),names(which(diff.super[4,]<4)),names(which(diff.super[9,]<4)))))
#sg2<-unique(append(names(which(diff.super[5,]<4)),c(names(which(diff.super[6,]<4)),names(which(diff.super[8,]<4)))))
#sg3<-names(which(diff.super[7,]<3))
mat<-1*round(sqrt(attr(dist.dna(red.nasohaps,as.matrix=T,p=T,variance=T),"variance")),3)
mat2<-matrix(0,9,9)
mat2[lower.tri(mat2,diag=F)]<-mat

## extend these sd to obtain fewer overlaps: trade-off tho...
sg1<-names(which(dist.dna.super[1,]<0.025))
sg2<-names(which(dist.dna.super[5,]<0.04))
sg3<-names(which(dist.dna.super[7,]<0.02))

overlap.12<-sg1[sg1%in%sg2]
overlap.13<-sg1[sg1%in%sg3]
overlap.23<-sg2[sg2%in%sg3]

super.hapbypop<-matrix(-1,5,109)
rownames(super.hapbypop)<-c("melissa","hybrid","third","unaccounted","overlap")
super.hapbypop[1,]<-colSums(ind.hap.pop[sg1[!(sg1%in%append(overlap.13,overlap.23))],])
super.hapbypop[2,]<-colSums(ind.hap.pop[sg2[!(sg2%in%append(overlap.13,overlap.23))],])
super.hapbypop[3,]<-colSums(ind.hap.pop[sg3[!(sg3%in%append(overlap.13,overlap.23))],])
super.hapbypop[4,]<-colSums(ind.hap.pop[!(rownames(nasohaps)%in%c(sg1,sg2,sg3)),])
super.hapbypop[5,]<-colSums(ind.hap.pop[c(overlap.13,overlap.12,overlap.23),])
colnames(super.hapbypop)<-colnames(ind.hap.pop)
write.table(super.hapbypop,file="super.hapbypop.csv")

## note: LIII and LIV are 1 step away from each other but form a group on their own
#diff.super["LIII","LIV"]<-0
#
#overlap.12<-names(which(diff.super["IV",]==0))[names(which(diff.super["IV",]==0)) %in% names(which(diff.super["CXI",]==0))]
## these next two contain haps with >90% missing data
#overlap.13<-names(which(diff.super["IV",]==0))[names(which(diff.super["IV",]==0)) %in% names(which(diff.super["LIII",]==0))]
#overlap.23<-names(which(diff.super["CXI",]==0))[names(which(diff.super["CXI",]==0)) %in% names(which(diff.super["LIII",]==0))]
#
### remove the above haps from the following calculations
#super.hapbypop<-matrix(-1,3,109)
#rownames(super.hapbypop)<-c("IV","CXI","LIII")
#for(i in c("IV","CXI","LIII")){
#    super.hapbypop[i,]<-colSums(ind.hap.pop.all[names(which(diff.super[i,]==0))[!(names(which(diff.super[i,]==0)) %in% append(overlap.12,overlap.13))],])
#}
#colnames(super.hapbypop)<-colnames(ind.hap.pop.all)
#rownames(super.hapbypop)<-c("melissa","hybrid","third")
#write.table(super.hapbypop,file="super.hapbypop.csv")
#
### haplotype indicator variable for each ind (for linear model analysis)
#super.hapbyind<-matrix(NA,3,sum(inf.ind[,4]=="TRUE"))
#rownames(super.hapbyind)<-c("IV","CXI","LIII")
#for(i in c("IV","CXI","LIII")){
#        super.hapbyind[i,]<-apply(ind.hap.pop.all[names(which(diff.super[i,]==0))[!(names(which(diff.super[i,]==0)) %in% append(overlap.12,overlap.13))],],2,sum)
#}
## contains TLA pop names
#colnames(super.hapbyind)<-tolower(substr(inf.ind$Ind[-which(inf.ind[,4]=="FALSE")],1,3))
#rownames(super.hapbyind)<-c("melissa","hybrid","third")
#write.table(super.hapbyind,file="super.hapbyind.csv")

## getting haplotype sequences into fasta format
unlink("allhaps.fasta")
for(i in rownames(nasohaps)){
    inds <- which(attr(nasohaps, "dimnames")[[1]] == i)
    seq <- attr(nasohaps, "index")[[inds]][1]

    new.hap<-naso[seq,]
    rownames(new.hap)<-rownames(nasohaps)[which(i==rownames(nasohaps))]

    write.FASTA(new.hap,file="allhaps.super.90pc.fasta",append=T)
}

unlink("redhaps.fasta")
for(i in rownames(red.nasohaps)){
    inds <- which(attr(red.nasohaps, "dimnames")[[1]] == i)
    seq <- attr(red.nasohaps, "index")[[inds]][5]

    new.hap<-naso[seq,]
    rownames(new.hap)<-rownames(red.nasohaps)[which(i==rownames(red.nasohaps))]

    write.FASTA(new.hap,file="redhaps.super.fasta",append=T)
}

## constructing haplotype tree
library(ggtree)

grp<-list(warsam="IV",melidas=c("LIII","LIV","LXVII"),annwhi="CXI",warsam2="II",warsam3="XCII")
redtree<-read.tree("redtrees_file")
redtree_plot<-ggtree(redtree)+geom_treescale()+geom_tiplab(size=4)
groupOTU(redtree_plot,grp,'Species')+aes(color=Species)+theme(legend.position="bottom")


alltree<-read.tree("alltrees_file")
grp<-list(reduced=redtree$tip.label,full=alltree$tip.label[!(alltree$tip.label %in% redtree$tip.label)])
alltree_plot<-ggtree(alltree)+geom_tiplab(size=2)
groupOTU(alltree_plot,grp,'Representation')+aes(color=Representation)+scale_color_manual(values=c("black","red"))+theme(legend.position="bottom")

### using the 'haplotypes' package
library(haplotypes)
xall<-read.fas("allhaps.fasta")

phap<-haplotype(as.dna(naso[miss.pc==0,]))
phap@haplist$haplotype33

pnet<-parsimnet(as.dna(naso[miss.pc==0,]),prob=NULL)

pnet.simp<-parsimnet(as.dna(naso[miss.pc==0,]))
rownames(pnet.simp@d$net1)

pdf("supernetA.pdf")
plot(pnet.simp,net=1,label=c(paste0("A",c(1:9)),rep(NA,6)),main="Net A",inter.labels=F)
dev.off()
pdf("supernetB.pdf")
plot(pnet.simp,net=2,label=c(paste0("B",c(1:30))),main="Net B",inter.labels=F)
dev.off()
pdf("supernetC.pdf")
plot(pnet.simp,net=3,label=c(paste0("C",c(1:4)),rep(NA,5)),main="Net C",inter.labels=F)
dev.off()
pdf("supernetALL.pdf")
plot(parsimnet(as.dna(naso[miss.pc==0,]),prob=NULL),main="Complete network",label=rep(NA,50))
dev.off()

lab<-rep(NA,50)
for(hap in 1:50){
    temp.net<-substr(names(which(unlist(pnet.simp@rowindex)==hap)),4,4)

    temp.hap<-eval(parse(text=paste0("which(pnet.simp@rowindex$net",temp.net,"==",hap,")")))

    if(temp.net=="1")
        lab[hap]<-paste0("A",temp.hap,"(",phap@freq[hap],")")
    else if(temp.net=="2")
        lab[hap]<-paste0("B",temp.hap,"(",phap@freq[hap],")")
    else if(temp.net=="3")
        lab[hap]<-paste0("C",temp.hap,"(",phap@freq[hap],")")
}

pdf("supernetAwinds.pdf")
plot(pnet.simp,net=1,label=c(lab[pnet.simp@rowindex$net1],rep(NA,6)),main="Net A",inter.labels=F)
dev.off()
pdf("supernetBwinds.pdf")
plot(pnet.simp,net=2,label=lab[pnet.simp@rowindex$net2],main="Net B",inter.labels=F)
dev.off()
pdf("supernetCwinds.pdf")
plot(pnet.simp,net=3,label=c(lab[pnet.simp@rowindex$net3],rep(NA,5)),main="Net C",inter.labels=F)
dev.off()
pdf("supernetALLwinds.pdf")
plot(pnet,label=c(lab,rep(NA,75)),main="Complete network",inter.labels=F,label.cex=0.6)
dev.off()

plot(pnet)

## capturing inds from haplotypes of each major network 
net1.inds<-unlist(sapply(pnet.simp@rowindex$net1,function(x){eval(parse(text=paste0("phap@haplist$haplotype",x)))}))
net2.inds<-unlist(sapply(pnet.simp@rowindex$net2,function(x){eval(parse(text=paste0("phap@haplist$haplotype",x)))}))
net3.inds<-unlist(sapply(pnet.simp@rowindex$net3,function(x){eval(parse(text=paste0("phap@haplist$haplotype",x)))}))
write.csv(net1.inds,"supernet1inds.txt")
write.csv(net2.inds,"supernet2inds.txt")
write.csv(net3.inds,"supernet3inds.txt")

xred<-read.fas("redhapsII.fasta")
pred<-parsimnet(xred,prob=NULL)
plot(pred)


