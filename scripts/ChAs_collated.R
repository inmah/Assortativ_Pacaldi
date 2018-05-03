#ChAs orginal scripts from VP
#Modified and collated by KC
#Calculates ChAs of chromatin features, prepares 3D data and overlaps with features
#This script first processes binarised ChIP-seq and then PCHiC, HiCap and ChiA-PET
#Once processed, can generate plots
#########
########from process_chromatinFeatures.R - prepares chromatin feature object

library(igraph)
library(GenomicRanges)

load('/home/rico/Repos/Assortativity/chromfeatures_mESC.Rdata') ##CrPs
CrPs[,1] <- paste(rep('chr', nrow(CrPs)),CrPs[,1], sep='')

##this bit needed to calculate ChAs to begin with
##puts CrPs into a GRanges object
bedepi <- with(CrPs, GRanges(chr, ranges=IRanges(start, end),strand=Rle(strand(rep('*', nrow(CrPs)))))   )
mcols(bedepi) <- CrPs[,-c(1,2,3)]
names(bedepi) <- paste(CrPs[,1], CrPs[,2], sep='_')

save(bedepi, file='bedepi.Rdata')

#########from Process_PCHiC_may16.R - prepares PCHiC data

source('assort_script.R') #contains function for calculating ChAs

#Read PCHiC data
#Table contains maps of bait regions and other ends. 'baitnames' and 'chicme' is in PCHiC_fragtypes.Rdata
#chicmore=read.table('mESC_PCHiC/mESC_wt_and_KO.txt', sep='\t', header=T) #not available here
#load('PCHiC_fragtypes.Rdata')
load('/data/Projects/kat/Projects/Assortativity/commsandmods.Rdata') ##chicmore table in here, contains PChIC interaction map
baitnames <- unique(chicmore$baitStart) #unique bait starts

length(which(chicmore$oeStart %in% baitnames)) ##OE same as bait
chicme <- unique(chicmore$oeStart[-which(chicmore$oeStart %in% baitnames)]) #unique other end starts

save(baitnames, chicme,file= 'PCHiC_fragtypes.Rdata')

#####map chromatin fragments to bedfile
#bait ends
v1 <- chicmore[,c(1,2,3)]
colnames(v1) <- c('chr', 'start', 'end')
#other ends
v2 <- chicmore[,c(6,7,8)]
colnames(v2) <- c('chr', 'start', 'end')
#combine and remove duplicates, concatenate chr with chr numbers
chicmore_all <- rbind(v1, v2)
chicmore_allnodup <- chicmore_all[!duplicated(chicmore_all),]
chicmore_allnodup[,1] <- paste(rep('chr', nrow(chicmore_allnodup)),chicmore_allnodup[,1], sep='')

#create GRanges object of PHiC interaction data
bedchicmore <- with(chicmore_allnodup, GRanges(chr, IRanges(start, end)))
bedchicmore$ID <- chicmore_allnodup[,2]

#######create feature table
#bedepi from process_chromatinFeatures.R and is the processed ChIP-seq data
#finding overlaps between PCHiC fragments and processed ChIP-seq features
#findOverlaps returns a 'Hits' object with number of matches of row number of bedchicmore and row number of bedepi
overlaps <- findOverlaps(bedchicmore,bedepi)

##matching overlapping bedchicmore IDs (start) with chip-seq features
match_hit <- data.frame(bedchicmore$ID[queryHits(overlaps)],as.data.frame(mcols(bedepi)[subjectHits(overlaps),] ),stringsAsFactors=T)
colnames(match_hit)[1] <- 'fragment'

##aggregate windows in fragments, collapses fragments
#install.packages("data.table")
library(data.table)
data.dt <- data.table(match_hit)
setkey(data.dt, fragment) #sorts ascending by fragment
b <- data.frame(data.dt[, lapply(.SD, mean), by = fragment]) #mean of ChIP-seq features by fragment
rownames(b) <- b[,1]
agchic <- b
agchic <- agchic[,-1]

save(agchic, file='ag_CHiC.Rdata')

#start making igraphs PCHiC
#Gchic <- graph.data.frame(chicmore[,c(2,7, 12,13,14)], vertices=NULL, directed=F)
#Gchic <- simplify(Gchic, remove.multiple=TRUE) #removes multiple edges going from and to the same place
#prom <- which(V(Gchic)$name %in% baitnames)

###filter interactions with CHiCAGO score>5 in mESC wt
chicwt <- chicmore[which(chicmore$mESC_wt>5),]

Gchicwt <- graph.data.frame(chicwt[,c(2,7, 12,13,14)], vertices=NULL, directed=F)  #IGRAPH UN-- 55845 72231

#remove edges AB BA
Gchicwt <- simplify(Gchicwt, remove.multiple=TRUE) #IGRAPH UN-- 55845 69987 

promwt <- which(V(Gchicwt)$name %in% baitnames)
oewt <- V(Gchicwt)$name[which(V(Gchicwt)$name %in% chicme)]

#remove nodes that are OE starts, i.e. get promoters only
Gchicwtprom <- delete.vertices(Gchicwt, V(Gchicwt)[which(V(Gchicwt)$name %in% chicme)])
#Select OE starts only and put into igraph
chicpromoe <- chicwt[which(chicwt[,7] %in% chicme),] 
Gchicwtpromoe <- graph.data.frame(chicpromoe[,c(2,7, 12,13,14)], vertices=NULL, directed=F)


save(Gchicwt, Gchicwtprom, Gchicwtpromoe, file='PCHiC_graphs.Rdata')

#calculate assortativity in all PChiC, PP only and PO only
ass_wtallepi <- calc_assort(Gchicwt,agchic)
ass_wtepiprom <- calc_assort(Gchicwtprom,agchic)
ass_wtepipromoe <- calc_assort(Gchicwtpromoe,agchic)

#calculates abundance of features
pchicpromab <- colMeans(agchic[which(rownames(agchic) %in% V(Gchicwtprom)$name),])
pchicpromoeab <- colMeans(agchic[which(rownames(agchic) %in% V(Gchicwtpromoe)$name),])
pchicab <- colMeans(agchic[which(rownames(agchic) %in% V(Gchicwt)$name),])

save(ass_wtallepi, ass_wtepiprom, ass_wtepipromoe, pchicab, pchicpromab, pchicpromoeab, file='PCHiCdata.Rdata')

############from Process_HiCap.R

#Making networks of HiCap
#Data from Sahlen et al
cappp <- read.table('Sahlen_interactions_promprom.csv', sep='\t', skip=1,header=T)
cappromoe <- read.table('Sahlen_interactions_promoe.csv', sep='\t', skip=4,header=T)
capoeoe <- read.table('Sahlen_interactions_oeoe.csv', sep='\t', skip=1,header=T)

##make networks, bait chr_start plus OE chr_start columns
capnetpp <- cbind(paste(cappp[,3], cappp[,4], sep='_'), paste(cappp[,11], cappp[,12], sep='_'))
capnetoeoe <- cbind(paste(capoeoe[,1], capoeoe[,2], sep='_'), paste(capoeoe[,5], capoeoe[,6], sep='_'))
capnetpoe<- cbind(paste(cappromoe[,6], cappromoe[,7], sep='_'), paste(cappromoe[,9], cappromoe[,10], sep='_'))
capnet <- rbind(capnetpp, capnetpoe, capnetoeoe) ##all interactions combined

##check if pp net has self edges
length(which(capnetpp[,1]==capnetpp[,2]))

Gcappp <- graph.data.frame(capnetpp, directed=F)
Gcapoe <- graph.data.frame(capnetoeoe, directed=F)
Gcappoe <- graph.data.frame(capnetpoe, directed=F)

save(Gchicwt, Gchicwtprom, Gchicwtpromoe, file='PCHiC_graphs.Rdata')

##read gff of fragments
samp1 <- read.table('GSE60494_HiCap_PromoterEnhancer_Interactions.gff', sep='\t', header=F)

##make bed file
colnames(samp1)[c(1,4,5)]=c('chr', 'start', 'end')
samp1bed <- with(samp1[,c(1,4,5)], GRanges(chr, ranges=IRanges(start, end)))
names(samp1bed) <- paste(samp1[,1], samp1[,4], sep='_')
mcols(samp1bed) <- samp1[,3]

##overlap epigenetic data with chromatin fragments
sampoverlaps <- findOverlaps(samp1bed,bedepi)
#create the match
samp_hit <- data.frame(names(samp1bed)[queryHits(sampoverlaps)],as.data.frame(mcols(bedepi)[subjectHits(sampoverlaps),] ),stringsAsFactors=T)
colnames(samp_hit)[1] <- 'fragment'

##aggregate windows in fragments
data.dt <- data.table(samp_hit)
setkey(data.dt, fragment)
b=data.frame(data.dt[, lapply(.SD, mean), by =fragment])
rownames(b) <- b[,1]
newsampag <- b

##make list of HiCap promoters
proms1 <- cappp[,c(3,4,4,1)]
colnames(proms1) <- c('chr', 'start', 'end','tss')
proms2 <- cappp[,c(11,12, 12,9)]
colnames(proms2) <- c('chr', 'start', 'end', 'tss')
proms3 <- cappromoe[,c(6,7,7,1)]
colnames(proms3) <- c('chr', 'start', 'end', 'tss')
proms <- rbind(proms1, proms2, proms3)
proms <- proms[!duplicated(proms),]
rownames(proms) <- paste(proms[,1], paste(proms[,2],proms[,4], sep='_'), sep='_')
namesproms <- rownames(proms)

save(namesproms, file='HiCap_proms.Rdata')

#make list of HiCap oes
oes <- cappromoe[,c(9,10,11)]

#give names
colnames(oes) <- c('chr', 'start', 'end')
oes <- oes[!duplicated(oes),]
namesoes <- paste(oes[,1], oes[,2], sep='_')

#make bedfiles
bedoes <- with(oes, GRanges(chr, ranges=IRanges(start, end)))
names(bedoes) <- namesoes

#Annotating chromatin states in fragments
####overlaps of proms fragment coordinates with tss coordinate for network
bedproms <- with(proms,GRanges(chr, ranges=IRanges(start, end) ))
names(bedproms) <- paste('TSS',rownames(proms), sep='')

oeoverlap <- findOverlaps(samp1bed,bedoes)
oehit <- data.frame(names(samp1bed)[queryHits(oeoverlap)],names(bedoes)[subjectHits(oeoverlap)] ,stringsAsFactors=T)
colnames(oehit)[1] <- 'fragment'
colnames(oehit)[2] <- 'node'
oehit <- oehit[!duplicated(oehit),]

promoverlap <- findOverlaps(samp1bed,bedproms)
promhit <- data.frame(names(samp1bed)[queryHits(promoverlap)],names(bedproms)[subjectHits(promoverlap)] ,stringsAsFactors=T)
colnames(promhit)[1] <- 'fragment'
colnames(promhit)[2] <- 'node'
promhit <- promhit[!duplicated(promhit),]

save(promhit,file='HiCap_proms.Rdata')

#renaming
promhit[,3] <- sapply(as.vector(promhit$node), function(x){
  res=unlist(strsplit(x, split='_'))
  res=paste(res[1], res[2],sep='_')
  res=gsub('TSS', '', res)
})

colnames(promhit) <- c('fragment', 'TSS', 'node')

promnear <- nearest(bedproms, samp1bed, select=c("arbitrary"),
                 ignore.strand=TRUE)

newpromhit <- cbind(names(bedproms), names(samp1bed)[promnear])

newpromhit[,1] <- sapply((newpromhit[,1]), function(x){
  res=unlist(strsplit(x, split='_'))
  res=paste(res[1], res[2],sep='_')
  res=gsub('TSS', '', res)
})

####making a dictionary of TSS coordinates and fragment starts
dicl <- list()
for (i in 1:nrow(newpromhit)){
  dicl[[newpromhit[i,1]]]=as.vector(newpromhit[i,2])
}

fraghit <- rbind(oehit,promhit[,c(1,3)])

promhit2 <- promhit[,c(1,2)]
colnames(promhit2) <- c('fragment', 'node')
fraghittss <- rbind(oehit,promhit2)

####see when tss is the same
sel <- grep('TSS',fraghittss[,2])
a <- sapply(as.vector(fraghittss[sel,2]),  function(x){
  vec=unlist(strsplit(x,split='_'))
  print(paste(vec[1], vec[3], collapse='_'))
  return(paste(vec[1], vec[3], collapse='_'))
})
fraghit_tssnew <- as.matrix(fraghittss)
fraghit_tssnew <- (cbind(fraghit_tssnew, rep('NA', nrow(fraghit_tssnew))))
fraghit_tssnew[sel,3] <- as.vector(a)

fraghit_nodup <- fraghit[!duplicated(fraghit_tssnew[,c(1,3)]),]

fraghit_fin <- fraghit_nodup[!duplicated(fraghit_nodup[,1]),]
rownames(fraghit_fin) <- fraghit_fin[,1]

fraghit_fin[as.vector(promhit[,1]),2] <- as.vector(promhit[,3])
frag_hitunique <- fraghit_fin[!duplicated(fraghit_fin),]

###rename rownames of newsampag with names of things in the network
common <- intersect(rownames(newsampag), frag_hitunique[,1])
newsampag2 <- cbind(frag_hitunique[common,2],newsampag[common,])

###adjust to make feature table like for pchic
ag_HiCap <- newsampag2[,-c(1,2)]

colnames(ag_HiCap) <- gsub('X5', '5', colnames(ag_HiCap))

save(ag_HiCap, file='ag_HiCap.Rdata')

############################################################################
###Calculate Chas and abs for HiCap

capnetppfin <- apply(capnetpp, 1, function(x){
  return(c(dicl[[x[1]]], dicl[[x[2]]] ))
})
capnetppfin <- t(capnetppfin)

capnetpoefin <- apply(capnetpoe, 1, function(x){
  return(c(dicl[[x[1]]], x[2] ))
})
capnetpoefin <- t(capnetpoefin)

captotal <- rbind(capnetppfin,capnetpoefin)

Gcapnetppfin <- graph.data.frame(capnetppfin[which(capnetppfin[,1]!=capnetppfin[,2]),], directed=F)
Gcapnetpoefin <- graph.data.frame(capnetpoefin, directed=F)
Gcaptotal <- graph.union(Gcapnetppfin, Gcapnetpoefin, byname=TRUE)

save(Gcapnetppfin, Gcapnetpoefin, Gcaptotal, file='HiCap_graphs.Rdata')

#don't think we need to save tables below?
#write.table(cbind(V(Gcapnetppfin)$name, rep('p',length(V(Gcapnetppfin)))), 'HiCap_ppfinproms.txt', quote=F,sep='\t', row.names=F)
#write.table(cbind(get.edgelist(Gcapnetppfin), rep('pp', length(E(Gcapnetppfin)))), 'HiCap_capnetppedges.txt', quote=F,sep='\t', row.names=F)

Nass_cappp <- calc_assort(Gcapnetppfin,ag_HiCap)
Nass_cappoe <- calc_assort(Gcapnetpoefin,ag_HiCap)
Nass_captotal <- calc_assort(Gcaptotal,ag_HiCap)

abcapnew <- colMeans(ag_HiCap[which(rownames(ag_HiCap) %in% V(Gcaptotal)$name),])
abcapppnew <- colMeans(ag_HiCap[which(rownames(ag_HiCap) %in% V(Gcapnetppfin)$name),])
abcappoenew <- colMeans(ag_HiCap[which(rownames(ag_HiCap) %in% V(Gcapnetpoefin)$name),])

save(Nass_cappp, Nass_cappoe, Nass_captotal, abcapnew, abcapppnew, abcappoenew, file='HiCapdata.Rdata')

##########from 3579114_process_chiapets.R
#process SMC1 mESC ChIA-PET

net <- read.table('CHiAPET_smc1_mESC.txt', sep='\t', header=T, skip=2)

colnames(net)[1:6] <- c('chr', 'start', 'end','chr', 'start', 'end')

Gchia <- graph.data.frame(net[,c(2,5)], directed=F)

bedchia <- with(rbind(net[,c(1:3)], net[,4:6]), GRanges(chr, IRanges(start, end)))
bedchia$ID <- c(net[,2], net[,5])
names(bedchia) <- bedchia$ID

####check overlaps of smc1 fragments with themselves
chiachia <- findOverlaps(bedchia,bedchia)
hit <- data.frame(names(bedchia)[queryHits(chiachia)],names(bedchia)[subjectHits(chiachia)] ,stringsAsFactors=T)
length(which(hit[,1]!=hit[,2]))


###overlap of features with fragments
chiaoverlaps <- findOverlaps(bedchia,bedepi)

#create the match
cmatch_hit <- data.frame(names(bedchia)[queryHits(chiaoverlaps)],as.data.frame(mcols(bedepi)[subjectHits(chiaoverlaps),] ),stringsAsFactors=T)
colnames(cmatch_hit)[1]='fragment'

##aggregate windows in fragments
data.dt <- data.table(cmatch_hit)
setkey(data.dt, fragment)
b <- data.frame(data.dt[, lapply(.SD, mean), by =fragment])
rownames(b) <- b[,1]
agchia <- b

agchia <- agchia[,-1]

colnames(agchia) <- gsub('X5', '5', colnames(agchia))
ab_SMC1=colMeans(agchia[which(rownames(agchia) %in% V(Gchia)$name),])

ass_chia <- calc_assort(Gchia,agchia)
ass_SMC1 <- ass_chia

#########from 5499299_process_chiapets.R

##RNA POL2
##from the paper: http://www.nature.com/nature/journal/v504/n7479/full/nature12716.html
#pol2net=read.table('Zhang2013_pol2chiaPET/ESC_pol2chia', sep='\t', header=T, skip=3) #not given
pol2net <- read.table('mESC_Pol2ChiA.csv', sep=",", header=T, skip=2) # supp table 3 from paper
colnames(pol2net)[1:6]=c('chr', 'start', 'end','chr', 'start', 'end')
##make bed file
#pol2net[,1]=gsub('chr', '', pol2net[,1]) #remove 'chr' from chr1 etc
#pol2net[,4]=gsub('chr', '', pol2net[,4])
Gpol2 <- graph.data.frame(pol2net[,c(2,5)], directed=F)

#save(Gpol2, file='ChiaPol2_graph.Rdata')
#load('ChiaPol2_graph.Rdata') #get Gpol2 from here, ESC_pol2chia table not given

bedpol2 <- with(rbind(pol2net[,c(1:3)], pol2net[,4:6]), GRanges(chr, IRanges(start, end)))
bedpol2$ID <- c(pol2net[,2], pol2net[,5])
names(bedpol2) <- bedpol2$ID

pol2pol2 <- findOverlaps(bedpol2,bedpol2)
hitpol2 <- data.frame(names(bedpol2)[queryHits(pol2pol2)],names(bedpol2)[subjectHits(pol2pol2)] ,stringsAsFactors=T)

pol2overlaps <- findOverlaps(bedpol2,bedepi)

#create the match
pol2match_hit <- data.frame(names(bedpol2)[queryHits(pol2overlaps)],as.data.frame(mcols(bedepi)[subjectHits(pol2overlaps),] ),stringsAsFactors=T)
colnames(pol2match_hit)[1] <- 'fragment'

##aggregate windows in fragments
data.dt <- data.table(pol2match_hit)
setkey(data.dt, fragment)
b <- data.frame(data.dt[, lapply(.SD, mean), by=fragment])
rownames(b) <- b[,1]
agpol2=b

ass_pol2 <- calc_assort(Gpol2,agpol2[,-c(1)])

ab_pol2 <- colMeans(agpol2[which(rownames(agpol2) %in% V(Gpol2)$name),-c(1)])

#contains RNAPII and SMC1 ChIA-PET data
save(agchia, agpol2, ass_pol2, ass_SMC1, ab_SMC1, ab_pol2, file='chiapets.Rdata')

###################################################################################################
#Now ready for plots!

###########from forpaper_allcombined.R

Chas_HiCap <- cbind(unlist(Nass_captotal),unlist(Nass_cappp), unlist(Nass_cappoe))
Ab_HiCap <- cbind(abcapnew, abcapppnew, abcappoenew)

Chas_PCHiC <- cbind(unlist(ass_wtallepi),unlist(ass_wtepiprom), unlist(ass_wtepipromoe))
Ab_PCHiC <- cbind(pchicab, pchicpromab, pchicpromoeab)

Chas_ChiaPets <- cbind(unlist(ass_SMC1),unlist(ass_pol2))
Ab_ChiaPets <- cbind(ab_SMC1, ab_pol2)

###table of values, need to check here they are all same length (79)
Chas_table <- cbind(Chas_PCHiC, Chas_HiCap, Chas_ChiaPets, Ab_PCHiC, Ab_HiCap, Ab_ChiaPets)
colnames(Chas_table) <- c('Chas_PCHiC','Chas_PCHiC_PP', 'Chas_PCHiC_PO','Chas_HiCap','Chas_HiCap_PP', 'Chas_HiCap_PO','Chas_ChiaSMC1','Chas_ChiaPol2'
                       ,'Ab_PCHiC','Ab_PCHiC_PP', 'Ab_PCHiC_PO','Ab_HiCap','Ab_HiCap_PP','Ab_HiCap_PO','Ab_ChiaSMC1','Ab_ChiaPol2')

#write.table(Chas_table, 'Chas_ab_table.txt', quote=F, sep='\t')

#####Now do correlation of everything
##Function to perform the correlation test for each pair of feature

#load('Chas_ab_table.Rdata')

corlist <- list()

findcor <-function(x){ 
  s=cor.test(mat[,x[1]],mat[,x[2]]) 
  return(c(s$est, s$p.value)) 
} 


mat <- Chas_table[,grep('Chas', colnames(Chas_table))]
names <- expand.grid(colnames(mat), colnames(mat)) ##transpose to make it in columns

corChas <- apply(names,1, findcor)

cors <- corChas[1,]
corps <- corChas[2,]

cmat <- cbind(names, cors, corps)

cormat <- matrix(0,nrow=ncol(mat), ncol=ncol(mat))
rownames(cormat) <- colnames(mat)
colnames(cormat) <- colnames(mat)
for(i in 1:ncol(mat)){
  for (j in 1:ncol(mat)){
    cormat[i,j]=cmat[which(cmat[,1]==colnames(mat)[i] & cmat[,2]==colnames(mat)[j]),3]
  }
}

pmat <- matrix(0,nrow=ncol(mat), ncol=ncol(mat))
rownames(pmat) <- colnames(mat)
colnames(pmat) <- colnames(mat)
for(i in 1:ncol(mat)){
  for (j in 1:ncol(mat)){
    pmat[i,j]=cmat[which(cmat[,1]==colnames(mat)[i] & cmat[,2]==colnames(mat)[j]),4]
  }
}
#############################repeat all for abundance


findcorab <-function(x){ 
  s=cor.test(abmat[,x[1]],abmat[,x[2]]) 
  return(c(s$est, s$p.value)) 
} 

abmat <- Chas_table[,grep('Ab', colnames(Chas_table))]
abnames <- expand.grid(colnames(abmat), colnames(abmat)) ##transpose to make it in columns

corab <- apply(abnames,1, findcorab)

abcors <- corab[1,]
abcorps <- corab[2,]

abcmat <- cbind(abnames, abcors, abcorps)

abcormat <- matrix(0,nrow=ncol(abmat), ncol=ncol(abmat))
rownames(abcormat)=colnames(abmat)
colnames(abcormat)=colnames(abmat)
for(i in 1:ncol(abmat)){
  for (j in 1:ncol(abmat)){
    abcormat[i,j]=abcmat[which(abcmat[,1]==colnames(abmat)[i] & abcmat[,2]==colnames(abmat)[j]),3]
  }
}

abpmat <- matrix(0,nrow=ncol(abmat), ncol=ncol(abmat))
rownames(abpmat)=colnames(abmat)
colnames(abpmat)=colnames(abmat)
for(i in 1:ncol(abmat)){
  for (j in 1:ncol(abmat)){
    abpmat[i,j]=abcmat[which(abcmat[,1]==colnames(abmat)[i] & abcmat[,2]==colnames(abmat)[j]),4]
  }
}

corlist[['pearson']]=list(cormat, pmat,abcormat, abpmat)
names(corlist[['pearson']])=c('cormat', 'pmat', 'abcormat', 'abpmat')

save(cormat, pmat,abcormat, abpmat, file='corandp_chasabcomparisons.Rdata')

#####repeat with spearman

findcorsp <-function(x){ 
  s=cor.test(mat[,x[1]],mat[,x[2]], method='sp') 
  return(c(s$est, s$p.value)) 
} 



mat <- Chas_table[,grep('Chas', colnames(Chas_table))]
names <- expand.grid(colnames(mat), colnames(mat)) ##transpose to make it in columns

corChas <- apply(names,1, findcorsp)

cors <- corChas[1,]
corps <- corChas[2,]

cmat <- cbind(names, cors, corps)

cormat <- matrix(0,nrow=ncol(mat), ncol=ncol(mat))
rownames(cormat)=colnames(mat)
colnames(cormat)=colnames(mat)
for(i in 1:ncol(mat)){
  for (j in 1:ncol(mat)){
    cormat[i,j]=cmat[which(cmat[,1]==colnames(mat)[i] & cmat[,2]==colnames(mat)[j]),3]
  }
}

pmat <- matrix(0,nrow=ncol(mat), ncol=ncol(mat))
rownames(pmat)=colnames(mat)
colnames(pmat)=colnames(mat)
for(i in 1:ncol(mat)){
  for (j in 1:ncol(mat)){
    pmat[i,j]=cmat[which(cmat[,1]==colnames(mat)[i] & cmat[,2]==colnames(mat)[j]),4]
  }
}
#############################repeat all for abundance

findcorabsp <-function(x){ 
  s=cor.test(abmat[,x[1]],abmat[,x[2]], method='sp') 
  return(c(s$est, s$p.value)) 
} 

abmat <- Chas_table[,grep('Ab', colnames(Chas_table))]
abnames <- expand.grid(colnames(abmat), colnames(abmat)) ##transpose to make it in columns

corab <- apply(abnames,1, findcorabsp)

abcors <- corab[1,]
abcorps=corab[2,]

abcmat <- cbind(abnames, abcors, abcorps)

abcormat <- matrix(0,nrow=ncol(abmat), ncol=ncol(abmat))
rownames(abcormat)=colnames(abmat)
colnames(abcormat)=colnames(abmat)
for(i in 1:ncol(abmat)){
  for (j in 1:ncol(abmat)){
    abcormat[i,j]=abcmat[which(abcmat[,1]==colnames(abmat)[i] & abcmat[,2]==colnames(abmat)[j]),3]
  }
}

abpmat <- matrix(0,nrow=ncol(abmat), ncol=ncol(abmat))
rownames(abpmat)=colnames(abmat)
colnames(abpmat)=colnames(abmat)
for(i in 1:ncol(abmat)){
  for (j in 1:ncol(abmat)){
    abpmat[i,j]=abcmat[which(abcmat[,1]==colnames(abmat)[i] & abcmat[,2]==colnames(abmat)[j]),4]
  }
}

save(cormat, pmat,abcormat, abpmat, file='corandp_chasabcomparisons_spearman_corrected.Rdata')

corlist[['spearman']]=list(cormat, pmat,abcormat, abpmat)
names(corlist[['spearman']])=c('cormat', 'pmat', 'abcormat', 'abpmat')

save(corlist, file='cor_chasandab_all_corrected.Rdata')

library(gplots)

#heatmap of ChA correlations
svg('figures/heatmap_corchas.svg')
heatmap.2(cormat, trace='none', revC=T, cexRow=0.7, cexCol=0.7)
dev.off()

svg('figures/heatmap_pvalchas.svg')
heatmap.2((pmat), trace='none', revC=T, cexRow=0.7, cexCol=0.7)
dev.off()

#color code features

namesvec <- colnames(agchic)
namesvec <- gsub('\\.', '_', namesvec)
nodecats <- (read.table('epinet_nodescats.txt', sep='\t'))

rownames(nodecats) <- nodecats[,1]
cats <- as.vector(unique(nodecats[,2]))
cols10 <- rainbow(11)
names(cols10) <- cats

cols10['Other'] <- 'grey'
cols4plot <- rep('grey', 78)
names(cols4plot) <- namesvec
for (c in cats){
  rel=which(nodecats[,2]==c)
  cols4plot[rownames(nodecats)[rel]]=cols10[c]
}

##################################  Plots of ChAs vs Abundance

svg('figures/AA_PCHiC_proms.svg')
labs <- c(1:ncol(agchic))
labs2 <- rep('', ncol(agchic))
sel <- which(ass_wtepiprom>0.13)
labs2[sel] <- namesvec[sel]
plot(pchicpromab, unlist(ass_wtepiprom),pch=20, col=cols4plot, main='PCHi-C P-P subnetwork', xlab='Abundance', 
     ylab='ChAS', xlim=c(0, 0.3))
text((pchicpromab),jitter(unlist(ass_wtepiprom), 2),pos=2,offset=0.2, labels=labs, cex=0.7)
text((pchicpromab),jitter(unlist(ass_wtepiprom), 2),pos=4,offset=0.2, labels=labs2, cex=0.5)

legend('bottomright', legend=as.vector(unique(nodecats[,2])), pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.7)
dev.off()


svg('figures/AA_PCHiC_promotherend.svg')
labs <- c(1:ncol(agchic))
labs2 <- rep('', ncol(agchic))
sel <- which(ass_wtepipromoe>=0.08 | ass_wtepipromoe<(-0.03) | pchicpromoeab>0.05)
labs2[sel]=namesvec[sel]
plot(pchicpromoeab, unlist(ass_wtepipromoe),pch=20, col=cols4plot, main='PCHi-C P-O subnetwork',xlim=c(0,0.2), xlab='Abundance', ylab='ChAS')
text((pchicpromoeab),jitter(unlist(ass_wtepipromoe), 2),pos=2,offset=0.2, labels=labs, cex=0.7)
text((pchicpromoeab),jitter(unlist(ass_wtepipromoe), 2),pos=4,offset=0.2, labels=labs2, cex=0.5)
abline(h=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])),bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.7)
dev.off()


svg('figures/AA_PCHiC_general.svg')
labs <- c(1:ncol(agchic))
labs2 <- rep('', ncol(agchic))
sel <- which(ass_wtallepi>0.11 | ass_wtallepi<0 | pchicab>0.05)
labs2[sel] <- namesvec[sel]
plot(pchicab, unlist(ass_wtallepi),pch=20,xlim=c(0,0.21), col=cols4plot, main='PCHi-C network', xlab='Abundance', ylab='ChAS')
text((pchicab),jitter(unlist(ass_wtallepi), 2),pos=2,offset=0.2, labels=labs, cex=0.7)
text((pchicab),jitter(unlist(ass_wtallepi), 2),pos=4,offset=0.2, labels=labs2, cex=0.7)
abline(h=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])),bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.7)
dev.off()


svg('figures/AA_HiCap_prom.svg')
labs <- c(1:length(Nass_captotal))
labs2 <- rep('', length(Nass_cappp))
sel=which(Nass_cappp>=0.06 | Nass_cappp<(0))
labs2[sel] <- namesvec[sel]
plot(abcapppnew, unlist(Nass_cappp),pch=20,xlim=c(0,0.5), col=cols4plot, main='HiCap P-P subnetwork', xlab='Abundance', ylab='ChAS')
text(abcapppnew, unlist(Nass_cappp),pos=2,offset=0.2, labels=labs, cex=0.7)
text(abcapppnew, unlist(Nass_cappp),pos=4,offset=0.2, labels=labs2, cex=0.7)
abline(h=0)
legend('topright', legend=as.vector(unique(nodecats[,2])), pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.7)
dev.off()

svg('figures/AA_HiCap_promoe.svg')
labs <- c(1:length(Nass_cappoe))
labs2 <- rep('', length(Nass_cappoe))
sel <- which(Nass_cappoe>=0.06 | Nass_cappoe<(0) | abcappoenew>0.1)
labs2[sel] <- namesvec[sel]
plot(abcappoenew, unlist(Nass_cappoe),pch=20,xlim=c(0,0.4), col=cols4plot, main='HiCap P-O subnetwork', xlab='Abundance', ylab='ChAS')
text(abcappoenew, unlist(Nass_cappoe),pos=2,offset=0.2, labels=labs, cex=0.7)
text(abcappoenew, unlist(Nass_cappoe),pos=4,offset=0.2, labels=labs2, cex=0.7)
abline(h=0)
legend('topright', legend=as.vector(unique(nodecats[,2])), pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.7)
dev.off()


svg('figures/AA_HiCap_general.svg')
labs <- c(1:length(Nass_captotal))
labs2 <- rep('', length(Nass_captotal))
sel <- which(Nass_captotal>=0.06 | Nass_captotal<(0)| abcapnew>0.1)
labs2[sel] <- namesvec[sel]
plot(abcapnew, unlist(Nass_captotal),pch=20, col=cols4plot,xlim=c(0,0.3), main='HiCap network', xlab='Abundance', ylab='ChAS')
text(abcapnew, unlist(Nass_captotal),pos=2,offset=0.2, labels=labs, cex=0.7)
text(abcapnew, unlist(Nass_captotal),pos=4,offset=0.2, labels=labs2, cex=0.7)
abline(h=0)
legend('topright', legend=as.vector(unique(nodecats[,2])), pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.7)
dev.off()

svg('figures/AA_Pol2chia_general.svg')
labs <- c(1:length(ass_pol2))
labs2 <- rep('', length(ass_pol2))
sel <- which(ass_pol2>=0.12 | ass_pol2<(0)| abpol2>0.2)
labs2[sel] <- namesvec[sel]
plot(abpol2, unlist(ass_pol2),pch=20, col=cols4plot,xlim=c(0,0.6), main='RNAPII ChIA-PET network', xlab='Abundance', ylab='ChAS')
text(abpol2, unlist(ass_pol2),pos=2,offset=0.2, labels=labs, cex=0.7)
text(abpol2, unlist(ass_pol2),pos=4,offset=0.2, labels=labs2, cex=0.7)
abline(h=0)
legend('topright', legend=as.vector(unique(nodecats[,2])), pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.7)
dev.off()

svg('figures/AA_SMC1chia_general.svg')
labs <- c(1:length(ass_SMC1))
labs2 <- rep('', length(ass_SMC1))
sel <- which(ass_SMC1>=0.35 | ab_SMC1>0.2)
labs2[sel] <- namesvec[sel]
plot(ab_SMC1, unlist(ass_SMC1),pch=20, col=cols4plot,xlim=c(0,0.4), main='SMC1 ChIA-PET network', xlab='Abundance', ylab='ChAS')
text(ab_SMC1, unlist(ass_SMC1),pos=2,offset=0.2, labels=labs, cex=0.7)
text(ab_SMC1, unlist(ass_SMC1),pos=4,offset=0.2, labels=labs2, cex=0.7)
abline(h=0)
legend('topright', legend=as.vector(unique(nodecats[,2])), pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.7)
dev.off()

#################################

##HiCap abundance PP vs PO
svg('figures/Ab_HiCap_pppo.svg')
labs <- c(1:length(Nass_cappp))
labs2 <- rep('', length(Nass_cappp))
sel <- which(abcapppnew>=0.1 | abcappoenew>0.05)
labs2[sel] <- names(Nass_cappp)[sel]
plot(abcapppnew, abcappoenew,pch=20,xlim=c(0,0.5), col=cols4plot, main='HiCap subnetworks Abundance', xlab='Abundance P-P', ylab='Abundance P-O')
text(abcapppnew, abcappoenew,pos=2,offset=0.2, labels=labs, cex=0.4)
text(abcapppnew, abcappoenew,pos=4,offset=0.2, labels=labs2, cex=0.4)
abline(h=0)
legend('topright', legend=as.vector(unique(nodecats[,2])), pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
dev.off()


##HiCap ChAs PP vs PO
svg('figures/CA_HiCap_pppo.svg')
labs <- c(1:length(Nass_cappp))
sel <- which(Nass_cappoe>0.09 |  Nass_cappoe<(-0.2)| Nass_cappp>0.05)
labs2 <- rep('', length(Nass_cappp))
labs2[sel] <- names(Nass_cappp)[sel]
plot(unlist(Nass_cappp), unlist(Nass_cappoe), pch=20, col=cols4plot,xlim=c(0,0.3), main='HiCap subnetworks ChAS', xlab='ChAS in P-P subnetwork', ylab='ChAS in P-O subnetwork')
text(unlist(Nass_cappp), unlist(Nass_cappoe),pos=2,offset=0.2, labels=labs, cex=0.4)
text(unlist(Nass_cappp), unlist(Nass_cappoe),pos=4,offset=0.2, labels=labs2, cex=0.4)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
dev.off()

###########PCHiC vs HiCap plots

svg('figures/Pchicvshicap.svg', height=8, width=5)
par(mfrow=c(3,2))
#ab pch ab hic NHiCap_vs_PCHIC_ab.svg')
labs <- c(1:ncol(agchic))
sel <- which(abcapnew>0.05 | abcapnew>0.04)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(pchicab, abcapnew,xlab='PCHi-C', ylab='HiCap',ylim=c(0,0.25), main='Abundance\n(HiCap vs PCHi-C)', pch=20,col=cols4plot)
text(pchicab, abcapnew, labels=labs,pos=2, offset=0.2, cex=0.4)
text(pchicab, jitter(abcapnew,10), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

#ass pchi ass hic NHiCap_vs_PCHIC_assort.svg')
labs <- c(1:ncol(agchic))
sel <- which(Nass_captotal>0.05 | ass_wtallepi>0.1)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(unlist(ass_wtallepi), unlist(Nass_captotal),pch=20, col=cols4plot,xlim=c(0,0.4), main='ChAS\n(HiCap vs PCHi-C)', xlab='PCHiC', ylab='HiCap')
text(unlist(ass_wtallepi),unlist(Nass_captotal), labels=labs,  pos=2, offset=0.2, cex=0.4)
text(unlist(ass_wtallepi),unlist(Nass_captotal), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

#ab PP pchh vs hica NHiCap_vs_PCHIC_promab.svg')
labs <- c(1:ncol(agchic))
sel <- which(abcapppnew>0.05 | abcapppnew>0.04)
labs2=rep('', ncol(agchic))
labs2[sel]=colnames(agchic)[sel]
plot(pchicpromab, abcapppnew,xlab='PCHi-C', ylab='HiCap',ylim=c(0,0.25), main='Abundance PP subnetwork\n(HiCap vs PCHi-C)', pch=20,col=cols4plot)
text(pchicpromab, abcapppnew, labels=labs,pos=2, offset=0.2, cex=0.4)
text(pchicpromab, jitter(abcapppnew,10), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

#ass PP vs hicap NewHiCap_vs_PCHIC_promass.svg')
labs <- c(1:ncol(agchic))
sel <- which( ass_wtepiprom>(0.1) | Nass_cappp>(0.1))
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(ass_wtepiprom, Nass_cappp,xlab='PCHiC ChAS in PP subnetwork', xlim=c(0,0.4), ylim=c(0, 0.4),ylab='HiCap ChAS in PP subnetwork', main='ChAS PP subnetwork\n(HiCap vs PCHi-C)', pch=20,col=cols4plot)
text(ass_wtepiprom, Nass_cappp, labels=labs,pos=2, offset=0.2, cex=0.4)
text(ass_wtepiprom, Nass_cappp,pos=4,offset=0.2, labels=labs2, cex=0.4)
legend('topleft', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
abline(a=0, b=1)

#ab PO pchic vs hicap nHiCap_vs_PCHIC_promoeab.svg')
labs <- c(1:ncol(agchic))
sel <- which(abcappoenew>0.05 | abcappoenew>0.04)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(pchicpromoeab, abcappoenew,xlab='PCHi-C', ylab='HiCap',ylim=c(0,0.25), main='Abundance PO subnetwork\n(HiCap vs PCHi-C)', pch=20,col=cols4plot)
text(pchicpromoeab, abcappoenew, labels=labs,pos=2, offset=0.2, cex=0.4)
text(pchicpromoeab, jitter(abcappoenew,10), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

#ass PO pchic vs hicap
labs <- c(1:ncol(agchic))
sel <- which(Nass_cappoe>0.05 | ass_wtepipromoe>0.1 | Nass_cappoe<(-0.1))
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(unlist(ass_wtepipromoe), unlist(Nass_cappoe),pch=20, col=cols4plot, main='ChAS PO subnetwork\n(HiCap vs PCHi-C)', xlab='PCHi-C', ylab='HiCap')
text(unlist(ass_wtepipromoe), unlist(Nass_cappoe), labels=labs,  pos=2, offset=0.2, cex=0.4)
text(unlist(ass_wtepipromoe), unlist(Nass_cappoe), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

dev.off()

#ChAs vs Chas plots

#PCHiC ChAs PP vs PO
svg('figures/CAplot_PCHiC.svg')
labs <- c(1:ncol(agchic))
sel <- which(ass_wtepiprom>0.1 |  ass_wtepipromoe<(-0.05)| ass_wtepipromoe>0.046)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(unlist(ass_wtepiprom), unlist(ass_wtepipromoe), pch=20, col=cols4plot,xlim=c(0,0.4), main='PCHi-C subnetworks ', xlab='ChAS in PP subnetwork', ylab='ChAS in PO subnetwork')
text(unlist(ass_wtepiprom), unlist(ass_wtepipromoe),pos=2,offset=0.2, labels=labs, cex=0.7)
text(unlist(ass_wtepiprom), unlist(ass_wtepipromoe),pos=4,offset=0.2, labels=labs2, cex=0.7)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('topleft', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.7)
dev.off()

#HiCap ChAs PP vs PO
svg('figures/HiCap_CAplot.svg')
labs <- c(1:length(Nass_cappp))
sel <- which(Nass_cappoe>0.09 |  Nass_cappoe<(-0.2)| Nass_cappp>0.05)
labs2 <- rep('', length(Nass_cappp))
labs2[sel] <- names(Nass_cappp)[sel]
plot(unlist(Nass_cappp), unlist(Nass_cappoe), pch=20, col=cols4plot,xlim=c(0,0.3), main='HiCap subnetworks', xlab='ChAS in PP subnetwork', ylab='ChAS in PO subnetwork')
text(unlist(Nass_cappp), unlist(Nass_cappoe),pos=2,offset=0.2, labels=labs, cex=0.7)
text(unlist(Nass_cappp), unlist(Nass_cappoe),pos=4,offset=0.2, labels=labs2, cex=0.7)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.7)
dev.off()

#############Figure of comparison pchic to hicap

#ChAs PCHiC PP vs HiCap PP
svg('figures/NewHiCap_vs_PCHIC_promass.svg')
labs <- c(1:ncol(agchic))
sel <- which( ass_wtepiprom>(0.1) | Nass_cappp>(0.1))
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(ass_wtepiprom, Nass_cappp,xlab='PCHiC ChAS in PP subnetwork', xlim=c(0,0.4), ylim=c(0, 0.4),ylab='HiCap ChAS in PP subnetwork', main='ChAS of PP subnetwork\n(HiCap vs PCHi-C)', pch=20,col=cols4plot)
text(ass_wtepiprom, Nass_cappp, labels=labs,pos=2, offset=0.2, cex=0.4)
text(ass_wtepiprom, Nass_cappp,pos=4,offset=0.2, labels=labs2, cex=0.4)
legend('topleft', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
#abline(a=0, b=1)
dev.off()

#Difference in ChAs PCHiC PO-PP and HiCap PO-PP
Nass_diff <- unlist(Nass_cappoe)-unlist(Nass_cappp)
Nass_diff_norm <- unlist(Nass_cappoe)/abcappoenew-unlist(Nass_cappp)/abcapppnew
diff_chic <- unlist(ass_wtepipromoe)-unlist(ass_wtepiprom)
diff_chic_norm <- unlist(ass_wtepipromoe)/pchicpromoeab-unlist(ass_wtepiprom)/pchicpromab

#pdf or svg?
#pdf('NewHiCap_vs_PCHIC_assdiff.pdf')
svg('figures/NHiCap_vs_PCHIC_assdiff.svg')
diffcor <- cor.test(diff_chic, (Nass_diff), method='pe')
labs <- c(1:ncol(agchic))
sel <- which( Nass_diff>(0.05) | diff_chic<(-0.1))
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(diff_chic, Nass_diff,xlab='ChAs difference in PCHi-C (PO-PP)',xlim=c(-0.31,0.1), ylim=c(-0.35, 0.2), ylab='ChAs difference in HiCap (PO-PP)',
     main=paste('ChAs difference PCHi-C vs. HiCap\n Pearson r =',round(diffcor$est,2),'p=0'), pch=20,col=cols4plot)
#text(diff_chic, Nass_diff, labels=labs,pos=2, offset=0.2, cex=0.4)
text(diff_chic, Nass_diff,pos=4,offset=0.2, labels=labs2, cex=0.4)
#legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
abline(h=0)
abline(v=0)
dev.off()

#Abundance PCHiC vs HiCap
svg('figures/NHiCap_vs_PCHIC_ab.svg')
labs <- c(1:ncol(agchic))
sel <- which(abcapnew>0.05 | abcapnew>0.04)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(pchicab, abcapnew,xlab='PCHi-C', ylab='HiCap',ylim=c(0,0.25), main='Abundance of features in whole network\n(HiCap vs PCHi-C)', pch=20,col=cols4plot)
text(pchicab, abcapnew, labels=labs,pos=2, offset=0.2, cex=0.4)
text(pchicab, jitter(abcapnew,10), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
dev.off()

#Abundance PCHiC PP vs HiCap PP
svg('figures/NHiCap_vs_PCHIC_promab.svg')
labs <- c(1:ncol(agchic))
sel <- which(abcapppnew>0.05 | abcapppnew>0.04)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(pchicpromab, abcapppnew,xlab='PCHi-C', ylab='HiCap',ylim=c(0,0.25), main='Abundance of features PP subnetwork\n(HiCap vs PCHi-C)', pch=20,col=cols4plot)
text(pchicpromab, abcapppnew, labels=labs,pos=2, offset=0.2, cex=0.4)
text(pchicpromab, jitter(abcapppnew,10), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
dev.off()

#Abundance PCHiC PO vs HiCap PO
svg('figures/NHiCap_vs_PCHIC_promoeab.svg')
labs <- c(1:ncol(agchic))
sel <- which(abcappoenew>0.05 | abcappoenew>0.04)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(pchicpromoeab, abcappoenew,xlab='PCHi-C', ylab='HiCap',ylim=c(0,0.25), main='Abundance of features PO subnetwork\n(HiCap vs PCHi-C)', pch=20,col=cols4plot)
text(pchicpromoeab, abcappoenew, labels=labs,pos=2, offset=0.2, cex=0.4)
text(pchicpromoeab, jitter(abcappoenew,10), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
dev.off()

#PCHiC ChAs vs HiCap ChAs all
svg('figures/NHiCap_vs_PCHIC_assort.svg')
labs <- c(1:ncol(agchic))
sel <- which(Nass_captotal>0.05 | ass_wtallepi>0.1)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(unlist(ass_wtallepi), unlist(Nass_captotal),pch=20, col=cols4plot,xlim=c(0,0.4), main='ChAs of whole network\n(HiCap vs PCHi-C)', xlab='PCHiC', ylab='HiCap')
text(unlist(ass_wtallepi),unlist(Nass_captotal), labels=labs,  pos=2, offset=0.2, cex=0.4)
text(unlist(ass_wtallepi),unlist(Nass_captotal), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
dev.off()

#PCHiC PP ChAs vs HiCap PP ChAs
svg('figures/NHiCap_vs_PCHIC_promassort.svg')
labs <- c(1:ncol(agchic))
sel <- which(Nass_cappp>0.05 | ass_wtepiprom>0.1)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(unlist(ass_wtepiprom), unlist(Nass_cappp),pch=20, col=cols4plot,xlim=c(0,0.4), main='ChAs of PP subnetwork\n(HiCap vs PCHi-C)', xlab='PCHi-C', ylab='HiCap')
text(unlist(ass_wtepiprom), unlist(Nass_cappp), labels=labs,  pos=2, offset=0.2, cex=0.4)
text(unlist(ass_wtepiprom), unlist(Nass_cappp), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
dev.off()

#PCHiC PO ChAs vs HiCap PO ChAs
svg('figures/NHiCap_vs_PCHIC_promoeassort.svg')
labs <- c(1:ncol(agchic))
sel <- which(Nass_cappoe>0.05 | ass_wtepipromoe>0.1 | Nass_cappoe<(-0.1))
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(unlist(ass_wtepipromoe), unlist(Nass_cappoe),pch=20, col=cols4plot, main='ChAs of PO subnetwork\n(HiCap vs PCHi-C)', xlab='PCHi-C', ylab='HiCap')
text(unlist(ass_wtepipromoe), unlist(Nass_cappoe), labels=labs,  pos=2, offset=0.2, cex=0.4)
text(unlist(ass_wtepipromoe), unlist(Nass_cappoe), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)
dev.off()

##Figures CA
#figures of PCHiC ChAs vs abundance and PO-PP HiCap vs PChIC

svg('figures/FigureCA.svg')
par(mfrow=c(2,2)) #4 panel plot

####AApromsPCHIC
labs <- c(1:ncol(agchic))
labs2 <- rep('', ncol(agchic))
sel <- which(ass_wtepiprom>0.13)
labs2[sel] <- colnames(agchic)[sel]
plot(pchicpromab, unlist(ass_wtepiprom),pch=20, col=cols4plot, main='PCHi-C PP subnetwork', xlab='Abundance', 
     ylab='ChAs', xlim=c(0, 0.3))
text((pchicpromab),jitter(unlist(ass_wtepiprom), 2),pos=2,offset=0.2, labels=labs, cex=0.4)
text((pchicpromab),jitter(unlist(ass_wtepiprom), 2),pos=4,offset=0.2, labels=labs2, cex=0.4)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

##AAPCHiC poe
labs <- c(1:ncol(agchic))
labs2 <- rep('', ncol(agchic))
sel <- which(ass_wtepipromoe>=0.08 | ass_wtepipromoe<(-0.03) | pchicpromoeab>0.05)
labs2[sel] <- colnames(agchic)[sel]
plot(pchicpromoeab, unlist(ass_wtepipromoe),pch=20, col=cols4plot, main='PCHi-C PO subnetwork',xlim=c(0,0.2), xlab='Abundance', ylab='ChAs')
text((pchicpromoeab),jitter(unlist(ass_wtepipromoe), 2),pos=2,offset=0.2, labels=labs, cex=0.4)
text((pchicpromoeab),jitter(unlist(ass_wtepipromoe), 2),pos=4,offset=0.2, labels=labs2, cex=0.4)
abline(h=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])),bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

##CApCHIC
labs <- c(1:ncol(agchic))
sel <- which(ass_wtepiprom>0.1 |  ass_wtepipromoe<(-0.05)| ass_wtepipromoe>0.046)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(unlist(ass_wtepiprom), unlist(ass_wtepipromoe), pch=20, col=cols4plot,xlim=c(0,0.4), main='PCHi-C subnetworks ', xlab='ChAS in PP subnetwork', ylab='ChAs in PO subnetwork')
text(unlist(ass_wtepiprom), unlist(ass_wtepipromoe),pos=2,offset=0.2, labels=labs, cex=0.4)
text(unlist(ass_wtepiprom), unlist(ass_wtepipromoe),pos=4,offset=0.2, labels=labs2, cex=0.4)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('topleft', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

####assdiff
labs <- c(1:ncol(agchic))
sel <- which( Nass_diff>(0.05) | diff_chic<(-0.2))
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(diff_chic, Nass_diff,xlab='PCHiC ChAS difference (PO-PP)',xlim=c(-0.4,0.2), ylim=c(-0.4, 0.2), ylab='HiCap ChAS difference (PO-PP)', main='Differential ChAs comparison\n(HiCap vs PCHi-C)', pch=20,col=cols4plot)
text(diff_chic, Nass_diff, labels=labs,pos=2, offset=0.2, cex=0.4)
text(diff_chic, Nass_diff,pos=4,offset=0.2, labels=labs2, cex=0.4)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

dev.off()

#######supfigs compare pchic hicap

svg('figures/Comparisons.svg', width=7, height=5)
par(mfrow=c(2,3))

#NHiCap_vs_PCHIC_ab.svg')
#Abundance HiCap vs PCHiC
labs <- c(1:ncol(agchic))
sel <- which(abcapnew>0.05 | abcapnew>0.04)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(pchicab, abcapnew,xlab='PCHi-C', ylab='HiCap',ylim=c(0,0.25), main='Abundance of features in whole network\n(HiCap vs PCHi-C)', cex.main=0.8, pch=20,col=cols4plot)
text(pchicab, abcapnew, labels=labs,pos=2, offset=0.2, cex=0.4)
text(pchicab, jitter(abcapnew,10), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

#Abundance HiCap PP vs PCHiC PP
#NHiCap_vs_PCHIC_promab.svg')
labs <- c(1:ncol(agchic))
sel <- which(abcapppnew>0.05 | abcapppnew>0.04)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(pchicpromab, abcapppnew,xlab='PCHi-C', ylab='HiCap',ylim=c(0,0.25), main='Abundance of features PP subnetwork\n(HiCap vs PCHi-C)', cex.main=0.8, pch=20,col=cols4plot)
text(pchicpromab, abcapppnew, labels=labs,pos=2, offset=0.2, cex=0.4)
text(pchicpromab, jitter(abcapppnew,10), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

#Abundance HiCap PO vs PCHiC PO
#NHiCap_vs_PCHIC_promoeab.svg')
labs <- c(1:ncol(agchic))
sel <- which(abcappoenew>0.05 | abcappoenew>0.04)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(pchicpromoeab, abcappoenew,xlab='PCHi-C', ylab='HiCap',ylim=c(0,0.25), main='Abundance of features PO subnetwork\n(HiCap vs PCHi-C)', cex.main=0.8, pch=20,col=cols4plot, cex=0.4)
text(pchicpromoeab, abcappoenew, labels=labs,pos=2, offset=0.2, cex=0.4)
text(pchicpromoeab, jitter(abcappoenew,10), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

#HiCap ChAs vs PCHiC ChAs all
#NHiCap_vs_PCHIC_assort.svg')
labs <- c(1:ncol(agchic))
sel <- which(Nass_captotal>0.05 | ass_wtallepi>0.1)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(unlist(ass_wtallepi), unlist(Nass_captotal),pch=20, col=cols4plot,xlim=c(0,0.4), main='ChAs of whole network\n(HiCap vs PCHi-C)', xlab='PCHiC', ylab='HiCap', cex.main=0.8)
text(unlist(ass_wtallepi),unlist(Nass_captotal), labels=labs,  pos=2, offset=0.2, cex=0.4)
text(unlist(ass_wtallepi),unlist(Nass_captotal), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

#HiCap ChAs PP vs PCHiC ChAs PP
#NHiCap_vs_PCHIC_promassort.svg')
labs <- c(1:ncol(agchic))
sel <- which(Nass_cappp>0.05 | ass_wtepiprom>0.1)
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(unlist(ass_wtepiprom), unlist(Nass_cappp),pch=20, col=cols4plot,xlim=c(0,0.4), main='ChAs of PP subnetwork\n(HiCap vs PCHi-C)', xlab='PCHi-C', ylab='HiCap', cex.main=0.8)
text(unlist(ass_wtepiprom), unlist(Nass_cappp), labels=labs,  pos=2, offset=0.2, cex=0.4)
text(unlist(ass_wtepiprom), unlist(Nass_cappp), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

#HiCap ChAs PO vs PCHiC ChAs PO
#NHiCap_vs_PCHIC_promoeassort.svg')
labs <- c(1:ncol(agchic))
sel <- which(Nass_cappoe>0.05 | ass_wtepipromoe>0.1 | Nass_cappoe<(-0.1))
labs2 <- rep('', ncol(agchic))
labs2[sel] <- colnames(agchic)[sel]
plot(unlist(ass_wtepipromoe), unlist(Nass_cappoe),pch=20, col=cols4plot, main='ChAs of PO subnetwork\n(HiCap vs PCHi-C)', xlab='PCHi-C', ylab='HiCap', cex.main=0.8)
text(unlist(ass_wtepipromoe), unlist(Nass_cappoe), labels=labs,  pos=2, offset=0.2, cex=0.4)
text(unlist(ass_wtepipromoe), unlist(Nass_cappoe), labels=labs2, pos=4, offset=0.2, cex=0.4)
abline(a=0, b=1)
abline(h=0)
abline(v=0)
legend('bottomright', legend=as.vector(unique(nodecats[,2])), bg='white', pch=20,col=cols10[as.vector(unique(nodecats[,2]))], cex=0.4)

dev.off()

#########RNA pol2 figures

polsel=names(ass_wtepiprom)[c(37,41, 43, 44, 42)]
polcol=c('#d0d1e6','#a6bddb','#67a9cf','#1c9099','#016c59')

#suppfigs plo2 abundance boxplot

load('HiCAP_proms.Rdata')

hicapproms=promhit$fragment #promoter bait regions

##Abundance of RNAPII features in PCHiC P and O, HiCap P and O
svg('figures/Pol2_frag_abundance.svg', width=5, height=8)
par(mfrow=c(2,2))
boxplot((agchic[which(rownames(agchic) %in% baitnames),polsel]), col=polcol, ylab='Fragment % covered by peak', names=polsel, outline=F, main='PCHi-C P', cex.axis=0.7, las=2)
boxplot((agchic[which(rownames(agchic) %in% chicme),polsel]), col=polcol, ylab='Fragment % covered by peak', names=polsel, outline=F, main='PCHi-C O', cex.axis=0.7, las=2)
boxplot((ag_HiCap[which(rownames(ag_HiCap) %in% hicapproms),polsel]), col=polcol, ylab='Fragment % covered by peak', names=polsel, outline=F, main='HiCap P', cex.axis=0.7, las=2)
boxplot((ag_HiCap[-which(rownames(ag_HiCap) %in% hicapproms),polsel]), col=polcol, ylab='Fragment % covered by peak', names=polsel, outline=F, main='HiCap O', cex.axis=0.7, las=2)

dev.off()


load('PCHiC_rewiredf.Rdata')
load('PCHiC_rewirenodupdf.Rdata')

#ChAs of RNAPII features in PChIC PP and PO
svg('figures/PCHiC_POL2Chas.svg', height=5, width=4)
barplot(cbind(unlist(ass_wtepiprom)[polsel], unlist(ass_wtepipromoe)[polsel]), col=polcol,ylim=c(0,0.45),cex.axis=0.7, cex.names=0.7, cex.lab=0.7,
        beside=T,ylab='ChAS', names=c('PCHiC PP', 'PCHiC PO'))
legend("topleft", 
       legend = polsel, cex=0.8, fill =polcol)
dev.off()

#Abundance of RNAPII features in PCHiC PP and PO
svg('figures/PCHiC_POL2Ab.svg', height=5, width=4)
barplot(cbind(pchicpromab[polsel],pchicpromoeab[polsel]), col=polcol,ylim=c(0,0.25),
        beside=T,main='Abundance', names=c('PCHiC PP', 'PCHiC PO'))
legend("topleft", 
       legend = polsel, cex=0.8, fill =polcol)
dev.off()

#Abundance of RNAPII features in HiCap PP and PO
svg('figures/HiCap_POL2Ab.svg', height=5, width=4)
barplot(cbind(abcapppnew[polsel],abcappoenew[polsel]), col=polcol,ylim=c(0,0.6),
        beside=T,main='Abundance', names=c('HiCap PP', 'HiCap PO'))
legend("topleft", 
       legend = polsel, cex=0.8, fill =polcol)
dev.off()

#ChAs of RNAPII features in HiCap PP and PO
svg('figures/Hicap_POL2Chas.svg', height=5, width=4)
barplot(cbind(unlist(Nass_cappp)[polsel], unlist(Nass_cappoe)[polsel]), col=polcol,cex.axis=0.7, cex.names=0.7, cex.lab=0.7,
        beside=T,ylab='ChAS', names=c('HiCap PP','HiCap PO') )
legend("topleft", 
       legend = polsel, cex=0.8, fill =polcol)
dev.off()

#ChAs of RNAPII features in PCHiC, HiCap, RNAPII and SMC1 CHIA-PET networks
svg('figures/chia_POL2Chas.svg', height=5, width=7)
barplot(cbind(unlist(ass_wtallepi)[polsel], unlist(Nass_captotal)[polsel], unlist(ass_pol2)[polsel], unlist(ass_SMC1)[polsel]), col=polcol,
        beside=T,ylab='ChAS', names=c('PCHiC', 'HiCap', 'RNAPII CHiA-PET', 'SMC1 ChIA-PET'))
legend("topleft", 
       legend = polsel, cex=0.8, fill =polcol)
dev.off()

#Abundance of RNAPII features in PCHiC, HiCap, RNAPII and SMC1 CHIA-PET networks
svg('figures/chia_POL2Ab.svg', height=5, width=7)
barplot(cbind(pchicab[polsel],abcapnew[polsel], abpol2[polsel], ab_SMC1[polsel]), col=polcol,
        beside=T,ylab='Abundance', names=c('PCHiC', 'HiCap', 'RNAPII CHiA-PET', 'SMC1 ChIA-PET'))
legend("topleft", 
       legend = polsel, cex=0.8, fill =polcol)
dev.off()

#############more RNAPII plots

#Abundance of RNAPII features in PCHiC and HiCap
svg('figures/POL2Ab_PCHICHicap.svg')
barplot(cbind(pchicab[polsel],abcapnew[polsel]), col=polcol,ylim=c(0,0.25),
        beside=T,ylab='Abundance', names=c('PCHi-C', 'HiCap'))
legend("topleft", 
       legend = polsel, cex=0.8, fill =polcol)
dev.off()

#ChAs of RNAPII features in PCHiC and HiCap
svg('figures/POL2ChAS_PCHICHicap.svg')
barplot(cbind(unlist(ass_wtallepi)[polsel],unlist(Nass_captotal)[polsel]), col=polcol,ylim=c(0,0.25),
        beside=T,ylab='ChAs', names=c('PCHi-C', 'HiCap'))
legend("topleft", 
       legend = polsel, cex=0.8, fill =polcol)
dev.off()


#figure pol2 chas  and ab
svg('figures/POL2_wholePCHICHicap.svg', height=5, width=10)
par(mfrow=c(1,2))

#svg('figures/POL2ChAS_PCHICHicap.svg')
barplot(cbind(unlist(ass_wtallepi)[polsel],unlist(Nass_captotal)[polsel]), col=polcol,ylim=c(-0.2,0.25),
        beside=T,ylab='ChAs', names=c('PCHi-C', 'HiCap'))
legend("topleft", 
       legend = polsel, cex=0.8, fill =polcol)

#svg('figures/POL2Ab_PCHICHicap.svg')
barplot(cbind(pchicab[polsel],abcapnew[polsel]), col=polcol,ylim=c(0,0.25),
        beside=T,ylab='Abundance', names=c('PCHi-C', 'HiCap'))
legend("topleft", 
       legend = polsel, cex=0.8, fill =polcol)

dev.off()

#### pp po comp ass and ab
##ChAs and abundance of RNAPII in PCHiC and HiCap PP and PO

svg('figures/suppfigs/Pol2_pppo_chasab.svg')
par(mfrow=c(2,2))
#svg('figures/PCHiC_POL2Chas.svg')
barplot(cbind(unlist(ass_wtepiprom)[polsel], unlist(ass_wtepipromoe)[polsel]), col=polcol,ylim=c(0,0.4),
        beside=T,ylab='ChAs', names=c('PCHiC PP', 'PCHiC PO'))
legend("topleft", 
       legend = polsel, cex=0.6, fill =polcol)

#svg('figures/PCHiC_POL2Ab.svg')
barplot(cbind(pchicpromab[polsel],pchicpromoeab[polsel]), col=polcol,ylim=c(0,0.25),
        beside=T,ylab='Abundance', names=c('PCHiC PP', 'PCHiC PO'))
legend("topleft", 
       legend = polsel, cex=0.6, fill =polcol)

#svg('figures/suppfigs/Hicap_POL2Chas.svg')
barplot(cbind(unlist(Nass_cappp)[polsel], unlist(Nass_cappoe)[polsel]), col=polcol,
        beside=T,ylab='ChAs', names=c('HiCap PP','HiCap PO') )
legend("topleft", 
       legend = polsel, cex=0.6, fill =polcol)

barplot(cbind(abcapppnew[polsel],abcappoenew[polsel]), col=polcol,ylim=c(0,0.6),
        beside=T,ylab='Abundance', names=c('HiCap PP', 'HiCap PO'))
legend("topleft", 
       legend = polsel, cex=0.6, fill =polcol)

dev.off()


##############supp fig pol2 4 comparison

#ChAs and Abundance of RNAPII features in all networks
svg('figures/Pol2_4comparisons.svg', height=5, width=10)
par(mfrow=c(1,2))

barplot(cbind(unlist(ass_wtallepi)[polsel], unlist(Nass_captotal)[polsel], unlist(ass_pol2)[polsel], unlist(ass_SMC1)[polsel]), col=polcol,
        beside=T,ylab='ChAs', names=c('PCHiC', 'HiCap', 'RNAPII\nCHiA-PET', 'SMC1\nChIA-PET'), cex.names=0.8)
legend("topleft", legend = polsel, cex=0.8, fill =polcol)

barplot(cbind(pchicab[polsel],abcapnew[polsel], abpol2[polsel], ab_SMC1[polsel]), col=polcol,
        beside=T,ylab='Abundance', names=c('PCHiC', 'HiCap', 'RNAPII\nCHiA-PET', 'SMC1\nChIA-PET'),cex.names=0.8)
legend("topleft", legend = polsel, cex=0.8, fill=polcol)

dev.off()

##try to identify active enhancers vs inactive enhancers
#load('HiCapG_nov15.Rdata')

########### Separation of enhancers in different types

Nhic_proms <- unique(promhit[,1])

#PCHiC types
enhancers <- rownames(agchic)[which(agchic[,'H3K4me1']>0)]
genicenh <- rownames(agchic)[which(agchic[,'H3K4me1']>0 & agchic[,'H3K36me3']>0)]
genebodynoenh <- rownames(agchic)[which(agchic[,'H3K4me1']==0 & agchic[,'H3K36me3']>0)]
enhancers_active <- rownames(agchic)[which(agchic[,'H3K4me1']>0 & agchic[,'H3K27ac']>0 )]
enhancers_poised <- rownames(agchic)[which(agchic[,'H3K4me1']>0 & agchic[,'H3K27me3']>0 )]
oenonenh <- chicme[-which(chicme %in% enhancers)]

#HiCap types
newsampag3 <- ag_HiCap
Nhic_enhancers <- rownames(newsampag3)[which(newsampag3[,'H3K4me1']>0)]
Nhic_enhancers_active <- rownames(newsampag3)[which(newsampag3[,'H3K4me1']>0 & newsampag3[,'H3K27ac']>0 )]
Nhic_enhancers_poised <- rownames(newsampag3)[which(newsampag3[,'H3K4me1']>0 & newsampag3[,'H3K27me3']>0 )]
Nhic_genebodynoenh <- rownames(newsampag3)[which(newsampag3[,'H3K4me1']==0 & newsampag3[,'H3K36me3']>0)]


##abundance in PCHiC enhancer types
ab_enh <- colMeans(agchic[which(rownames(agchic) %in% enhancers),])
ab_enhac <- colMeans(agchic[which(rownames(agchic) %in% enhancers_active),])
ab_enhpoi <- colMeans(agchic[which(rownames(agchic) %in% enhancers_poised),])
ab_proms <- colMeans(agchic[which(rownames(agchic) %in% baitnames),])
ab_noenh <- colMeans(agchic[which(rownames(agchic) %in% genebodynoenh),])

#Abundance of ChIP-seq features in PChIC non-enhancers
svg('figures/Nonenh_ab.svg', height=5, width=8)
barplot(sort(ab_noenh, dec=T), las=2, cex.names=0.5, main='Abundance of features in Other end fragments with no H3K4me1 enhancer mark')
dev.off()

##hic abundance in HiCap enhancer types
Nhicab_enh <- colMeans(newsampag3[which(rownames(newsampag3) %in% Nhic_enhancers),-c(1,2)])
Nhicab_enhac <- colMeans(newsampag3[which(rownames(newsampag3) %in% Nhic_enhancers_active),-c(1,2)])
Nhicab_enhpoi <- colMeans(newsampag3[which(rownames(newsampag3) %in% Nhic_enhancers_poised),-c(1,2)])
Nhicab_proms <- colMeans(newsampag3[which(rownames(newsampag3) %in% Nhic_proms),-c(1,2)])
Nhicab_noenh <- colMeans(newsampag3[which(rownames(newsampag3) %in% Nhic_genebodynoenh),-c(1,2)])


###sum of pol2 in PCHiC enhancers types
sum_enh <- colSums(agchic[which(rownames(agchic) %in% enhancers),])
sum_enhac <- colSums(agchic[which(rownames(agchic) %in% enhancers_active),])
sum_enhpoi <- colSums(agchic[which(rownames(agchic) %in% enhancers_poised),])
sum_proms <- colSums(agchic[which(rownames(agchic) %in% baitnames),])

###RNAPII feature abundance in PChIC enhancer types
#png('Pol2abundance_barplotenhtypes.png')
svg('figures/NPol2abundance_barplotenhtypes.svg')
#pdf('Pol2abundance_barplotenhtypes.pdf')
barplot(cbind(ab_proms[polsel], ab_enhac[polsel], ab_enhpoi[polsel]),col=polcol,
        beside=T, ylim=c(0, 0.2),ylab='Abundance', main='RNA Pol 2 abundance PCHiC',names=c('P', 'Active enhancer','Poised enhancer'), cex.names=1)

dev.off()

#pdf('Pol2abundance_NHICAPbarplotenhtypes.pdf')
svg('Pol2abundance_NHICAPbarplotenhtypes.svg')
barplot(cbind(Nhicab_proms[polsel],Nhicab_enhac[polsel], Nhicab_enhpoi[polsel]),col=polcol,
        beside=T, ylab='Abundance',ylim=c(0,0.6), main='RNA Pol 2 abundance HiCap',names=c('P', 'Active enhancer','Poised enhancer'), cex.names=0.7)
legend("topleft", legend = polsel, cex=0.6, fill =polcol)
dev.off()


pdf('Pol2sum_barplotenhtypes.pdf')
barplot(cbind(sum_proms[polsel],sum_enh[polsel],sum_enh[polsel], sum_enhac[polsel]),col=polcol,
        beside=T, ylab='Sum', main='Pol 2 Sum',names=c('P','Enhancer', 'Active enhancer','Poised enhancer'), cex.names=0.7)
legend("topleft", legend = polsel, cex=0.6, fill =polcol)
dev.off()

#write.table(cbind(oenonenh,rep('nonenh', length(oenonenh))), 'nonenh_forcyto.txt', quote=F, sep='\t', row.names=F, col.names=F)

#####calculate assortativity of RNAPII features in enhancer types

#enhancer subnets of PCHIC
Gchicwtpromoe_enh=induced.subgraph(Gchicwtpromoe, V(Gchicwtpromoe)[which(V(Gchicwtpromoe)$name %in% c(baitnames,enhancers))])
Gchicwtpromoe_enhac=induced.subgraph(Gchicwtpromoe, V(Gchicwtpromoe)[which(V(Gchicwtpromoe)$name %in% c(baitnames,enhancers_active))])
Gchicwtpromoe_enhpoi=induced.subgraph(Gchicwtpromoe, V(Gchicwtpromoe)[which(V(Gchicwtpromoe)$name %in% c(baitnames,enhancers_poised))])
Gchicwtpromoe_promnonenh=induced.subgraph(Gchicwtpromoe, V(Gchicwtpromoe)[which(V(Gchicwtpromoe)$name %in% c(baitnames,oenonenh))])
Gchicwtpromoe_nonenh=induced.subgraph(Gchicwtpromoe, V(Gchicwtpromoe)[which(V(Gchicwtpromoe)$name %in% c(baitnames,oenonenh))])

#enhancer subnets of HICAP
NGcappoe_enh=induced.subgraph(Gcapnetppfin, V(Gcapnetppfin)[which(V(Gcapnetppfin)$name %in% c(Nhic_proms,Nhic_enhancers))])
NGcappoe_enhac=induced.subgraph(Gcapnetppfin, V(Gcapnetppfin)[which(V(Gcapnetppfin)$name %in% c(Nhic_proms,Nhic_enhancers_active))])
NGcappoe_enhpoi=induced.subgraph(Gcapnetppfin, V(Gcapnetppfin)[which(V(Gcapnetppfin)$name %in% c(Nhic_proms,Nhic_enhancers_poised))])
NGcappoe_noenh=induced.subgraph(Gcapnetppfin, V(Gcapnetppfin)[which(V(Gcapnetppfin)$name %in% c(Nhic_proms,Nhic_genebodynoenh))])

##RNAPII ChAs in PCHiC enhancer types
ass_promoeenh=calc_assort(Gchicwtpromoe_enh,agchic)
ass_promoeenhac=calc_assort(Gchicwtpromoe_enhac,agchic)
ass_promoeenhpoi=calc_assort(Gchicwtpromoe_enhpoi,agchic)
ass_promoenonenh=calc_assort(Gchicwtpromoe_promnonenh,agchic)

#oNGcap_enhac=induced.subgraph(Gcaptotal, V(Gcaptotal)[which(V(Gcaptotal)$name %in% c(Nhic_proms,Nhic_enhancers_active))])

##RNAPII ChAs in HiCap enhancer types
Nass_promoeenh=calc_assort(NGcappoe_enh,newsampag3)
Nass_promoeenhac=calc_assort(NGcappoe_enhac,newsampag3)
Nass_promoeenhpoi=calc_assort(NGcappoe_enhpoi,newsampag3)
Nass_promoenoenh=calc_assort(NGcappoe_noenh,newsampag3)

####comparison of enhancer types combined PO network

#RNAPII ChAs in PCHiC enhancer types
svg('figures/enhancertypes_PCHIC_onlyPOnet.svg', height=8, width=11)
barplot(cbind(unlist(ass_wtepiprom[polsel]), unlist(ass_promoeenhac[polsel]),unlist(ass_promoeenhpoi[polsel]),unlist(ass_promoenonenh[polsel])),main='PCHi-C ChAS', 
        col=polcol, beside=T,ylim=c(-0.18,0.32),cex.names=1.5,cex.axis=1.5 ,cex.lab=1.5, ylab='ChAS', names=c('PP', 'P-\nActive enhancer', 'P-\nPoised Enhancer',' \nP-\nNon-enhancer'))
legend("topleft", legend = names(ass_wtepiprom[polsel]), cex=1,fill =polcol, horiz=F)
dev.off()

#RNAPII ChAs in HiCap enhancer types
svg('figures/enhancertypes_HiCap_onlyPOnet.svg', height=8, width=11)
barplot(cbind(unlist(Nass_cappp[polsel]), unlist(Nass_promoeenhac[polsel]),unlist(Nass_promoeenhpoi[polsel]),unlist(Nass_promoenoenh[polsel])),main='HiCap ChAs',
        col=polcol, beside=T,ylim=c(-0.18,0.32),cex.names=1.5,cex.axis=1.5 ,cex.lab=1.5, ylab='ChAS', names=c('PP', 'P-\nActive enhancer', 'P-\nPoised Enhancer','P-\nNon-enhancer'))
legend("topleft", 
       legend = names(ass_wtepiprom[polsel]), cex=1,
       fill =polcol, horiz=F)
dev.off()

##Abundance of RNAPII features in enhancer types in PCHIC and HICAP
svg('figures/enhancertypes_ab.svg') #added
par(mfrow=c(1,2)) #added

##pchic ab
barplot(cbind(ab_proms[polsel], ab_enhac[polsel], ab_enhpoi[polsel]),col=polcol,
        beside=T, ylim=c(0, 0.2),ylab='Abundance', main='PCHi-C Abundance',names=c('P', 'Active enhancer','Poised enhancer'), cex.names=0.7)
legend("topleft", 
       legend = polsel, cex=0.6, fill =polcol)

#Hicap ab
barplot(cbind(Nhicab_proms[polsel],Nhicab_enhac[polsel], Nhicab_enhpoi[polsel]),col=polcol,
        beside=T, ylab='Abundance',ylim=c(0,0.6), main='HiCAP Abundance',names=c('P', 'Active enhancer','Poised enhancer'), cex.names=0.7)
legend("topleft", 
       legend = polsel, cex=0.6, fill =polcol)

dev.off()

#############next network rewiring, go to netrewirings.R to generate data. This part bit unclear



########Now calculate beanplots of different types of enhancers using partially rewired networks

###data from netrewirings.R
#run loops in netrewirings from command line, takes a while. Once run, save RData and load here
load('PCHiC_prerewired.Rdata')
load('HiCap_prerewired.Rdata')
load('HiCap_rewireenhtypesdf.Rdata')
load('PCHiC_rewirenodupdf.Rdata')
load('HiCap_rewiredf.Rdata')

##rewired pchic
enhlist=list()
for (i in 1:5){
  enhlist[[i]]=as.data.frame(cbind(as.numeric(capassprerewdf_enh[polsel[i],]), rep(polsel[i], nrow(capassprerewdf_enh))))
  colnames(enhlist[[i]])=c('val', 'var')
  enhlist[[i]]$var=as.factor(enhlist[[i]]$var)
}

enh=rbind(enhlist[[1]],enhlist[[2]], enhlist[[3]], enhlist[[4]], enhlist[[5]])

enhaclist=list()
for (i in 1:5){
  enhaclist[[i]]=as.data.frame(cbind(as.numeric(capassprerewdf_enhac[polsel[i],]), rep(polsel[i], nrow(capassprerewdf_enhac))))
  colnames(enhaclist[[i]])=c('val', 'var')
  enhaclist[[i]]$var=as.factor(enhaclist[[i]]$var)
}

enhac=rbind(enhaclist[[1]],enhaclist[[2]], enhaclist[[3]], enhaclist[[4]], enhaclist[[5]])

enhpoilist=list()
for (i in 1:5){
  enhpoilist[[i]]=as.data.frame(cbind(as.numeric(capassprerewdf_enhpoi[polsel[i],]), rep(polsel[i], nrow(capassprerewdf_enhpoi))))
  colnames(enhpoilist[[i]])=c('val', 'var')
  enhpoilist[[i]]$var=as.factor(enhpoilist[[i]]$var)
}

enhpoi=rbind(enhpoilist[[1]],enhpoilist[[2]], enhpoilist[[3]], enhaclist[[4]], enhpoilist[[5]])

noenhlist=list()
for (i in 1:5){
  noenhlist[[i]]=as.data.frame(cbind(as.numeric(capassprerewdf_noenh[polsel[i],]), rep(polsel[i], nrow(capassprerewdf_noenh))))
  colnames(noenhlist[[i]])=c('val', 'var')
  noenhlist[[i]]$var=as.factor(noenhlist[[i]]$var)
}

noenh=rbind(noenhlist[[1]],noenhlist[[2]], noenhlist[[3]], noenhlist[[4]], noenhlist[[5]])


svg('figures/rewired_PCHiC_enhancertypes.svg')
par(mfrow=c(1,3))

beanplot(as.numeric(as.vector(enhac$val))~enhac$var, what=c(1,1,1,0), col=as.list(polcol), las=2, ylim=c(-0.1,0.25), main='Rewired PCHi-C\nP-Active Enhancer', ylab='ChAs', names=rep('',5))
abline(h=0)
legend("topleft", legend = polsel, cex=0.8, fill =polcol)

beanplot(as.numeric(as.vector(enhpoi$val))~enhpoi$var, what=c(1,1,1,0), col=as.list(polcol), las=2, ylim=c(-0.1,0.25), main='Rewired PCHi-C\nP-Poised Enhancer', ylab='ChAs', names=rep('',5))
abline(h=0)

beanplot(as.numeric(as.vector(noenh$val))~noenh$var, what=c(1,1,1,0), col=as.list(polcol), las=2, ylim=c(-0.1,0.25), main='Rewired PCHi-C\nP-Non-Enhancer', ylab='ChAs', names=rep('',5))
abline(h=0)
dev.off()


###make lists of different enhancer types HiCap
capenhlist=list()
for (i in 1:5){
  capenhlist[[i]]=as.data.frame(cbind(as.numeric(capassprerewdf_enh[polsel[i],]), rep(polsel[i], nrow(capassprerewdf_enh))))
  colnames(capenhlist[[i]])=c('val', 'var')
  capenhlist[[i]]$var=as.factor(capenhlist[[i]]$var)
}

capenh=rbind(capenhlist[[1]],capenhlist[[2]], capenhlist[[3]], capenhlist[[4]], capenhlist[[5]])

capenhaclist=list()
for (i in 1:5){
  capenhaclist[[i]]=as.data.frame(cbind(as.numeric(capassprerewdf_enhac[polsel[i],]), rep(polsel[i], nrow(capassprerewdf_enhac))))
  colnames(capenhaclist[[i]])=c('val', 'var')
  capenhaclist[[i]]$var=as.factor(capenhaclist[[i]]$var)
}

capenhac=rbind(capenhaclist[[1]],capenhaclist[[2]], capenhaclist[[3]], capenhaclist[[4]], capenhaclist[[5]])

capenhpoilist=list()
for (i in 1:5){
  capenhpoilist[[i]]=as.data.frame(cbind(as.numeric(capassprerewdf_enhpoi[polsel[i],]), rep(polsel[i], nrow(capassprerewdf_enhpoi))))
  colnames(capenhpoilist[[i]])=c('val', 'var')
  capenhpoilist[[i]]$var=as.factor(capenhpoilist[[i]]$var)
}

capenhpoi=rbind(capenhpoilist[[1]],capenhpoilist[[2]], capenhpoilist[[3]], capenhaclist[[4]], capenhpoilist[[5]])

nocapenhlist=list()
for (i in 1:5){
  nocapenhlist[[i]]=as.data.frame(cbind(as.numeric(capassprerewdf_noenh[polsel[i],]), rep(polsel[i], nrow(capassprerewdf_enh))))
  colnames(nocapenhlist[[i]])=c('val', 'var')
  nocapenhlist[[i]]$var=as.factor(nocapenhlist[[i]]$var)
}

nocapenh=rbind(nocapenhlist[[1]],nocapenhlist[[2]], nocapenhlist[[3]], nocapenhlist[[4]], nocapenhlist[[5]])


###beanplot
#install.packages("beanplot")
library(beanplot)

svg('figures/rewired_HiCap_enhancertypes.svg')
par(mfrow=c(1,3))

beanplot(as.numeric(as.vector(capenhac$val))~capenhac$var, what=c(1,1,1,0), col=as.list(polcol), las=2, ylim=c(-0.1,0.25), main='Rewired HiCap\nP-Active Enhancer', ylab='ChAs', names=rep('',5))
abline(h=0)
legend("topleft", legend = polsel, cex=0.8, fill =polcol)

beanplot(as.numeric(as.vector(capenhpoi$val))~capenhpoi$var, what=c(1,1,1,0), col=as.list(polcol), las=2, ylim=c(-0.1,0.25), main='Rewired HiCap\nP-Poised Enhancer', ylab='ChAs', names=rep('',5))
abline(h=0)

beanplot(as.numeric(as.vector(nocapenh$val))~nocapenh$var, what=c(1,1,1,0), col=as.list(polcol), las=2, ylim=c(-0.1,0.25), main='Rewired HiCap\nP-Non-Enhancer', ylab='ChAs', names=rep('',5))
abline(h=0)
dev.off()


