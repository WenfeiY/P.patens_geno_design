setwd("./")
library(GenomicFeatures)
library(Biostrings)
library(GenomicRanges)
library(dplyr)
library(rtracklayer)
library(BSgenome)
library(IRanges)
library(DECIPHER)
library(data.table)
library(tidyverse)

# 1. Import and process sequence information
fasta_genome<-readDNAStringSet('P.patens.genome.fasta')# wildtype genome sequence 
gff_genome<-import('P.patens.genome.gff3')# wildtype genome annotation
annotation_genome<- makeTxDbFromGFF('P.patens.genome.gff3')
genome_gene<-genes(annotation_genome)
seqlengths(genome_gene)<-as.vector(width(fasta_genome))

# 2. Introduction of PCRmarks
# 2.1 Distinguish between intron and CDS, and standardize the design range according to pcrmarks requirements
inByTr_new<-intronsByTranscript(annotation_genome,use.names = TRUE)# introns grouped by transcript
inByTr_new<-inByTr_new[lengths(inByTr_new)!=0]
for(i in 1:length(inByTr_new)){
  inByTr_new[[i]]$trc<-names(inByTr_new[i])
}
inByTr_new_ul<-unlist(inByTr_new)
inByTr_new_ul$gene<-sapply(inByTr_new_ul$trc,function(x){strsplit(x,'\\.')[[1]][1]},USE.NAMES = F)
inByTr_new_ul$type<-'intron'
inByTr_new_unique<-unique(inByTr_new_ul)

cdsByTr<-cdsBy(annotation_genome,by='tx',use.names=T)# CDS grouped by transcript
for (i in seq_len(length(cdsByTr))) {
  if(unique(as.vector(strand(cdsByTr[[i]])))=='+'){
    phase <- (3-cumsum(width(cdsByTr[[i]]))%%3)%%3
    cdsByTr[[i]]$phase <- c(0,phase[-length(phase)])}
  else{
    cdsByTr[[i]]<-rev(cdsByTr[[i]])
    phase <- (3-cumsum(width(cdsByTr[[i]]))%%3)%%3
    cdsByTr[[i]]$phase <- c(0,phase[-length(phase)])
  }
}

cdsByTr_unique<-unique(unlist(cdsByTr))
cdsByTr_unique$wid<-width(cdsByTr_unique)
cdsByTr_unique$trc<-names(cdsByTr_unique)
cdsByTr_unique$gene<-sapply(names(cdsByTr_unique),function(x){strsplit(x,'\\.')[[1]][1]},USE.NAMES = F)
cdsByTr_unique$type<-'cds'
cdsByTr_unique<-cdsByTr_unique[,-c(1:3)]

cdsByTr_1<-unlist(cdsByTr)
cdsByTr_1$wid<-width(cdsByTr_1)
cdsByTr_1$trc<-names(cdsByTr_1)
cdsByTr_1$gene<-sapply(names(cdsByTr_1),function(x){strsplit(x,'\\.')[[1]][1]},USE.NAMES = F)
cdsByTr_1$type<-'cds'
cdsByTr_1<-cdsByTr_1[,-c(1:3)]
cdsByTr_1<-split(cdsByTr_1,cdsByTr_1$trc)

cds_o_intr<-(subjectHits(findOverlaps(cdsByTr_unique,inByTr_new_unique)))
inByTr_new_unique<-inByTr_new_unique[-unique(cds_o_intr)]
disjoin_intron<-function(x){
  gr2<-x
  gr2<-split(gr2,gr2$gene)
  for(i in 1:length(gr2)){
    gr2[[i]]<-disjoin(gr2[[i]])
  }
  gr2<-unlist(gr2)
  gr2$wid<-width(gr2)
  gr2$gene<-names(gr2)
  gr2$type<-'intron'
  return(gr2)
}# cut introns with overlapping sequences in different transcripts
inByTr_new_unique<-disjoin_intron(inByTr_new_unique)
inByTr_new_unique<-inByTr_new_unique[inByTr_new_unique$wid>84]# remove introns whose length is less than 10+50+24 = 84 nt
mod_range<-function(x,st,ed){
  for(i in 1:length(x)){
    if(as.vector(strand(x[i]))=='+'){
      start(x[i])<-start(x[i])+st
      end(x[i])<-end(x[i])-ed
    }else{
      start(x[i])<-start(x[i])+ed
      end(x[i])<-end(x[i])-st
    }
  }
  return(x)
}# narrow down the range of intron/cds
inByTr_new_unique<-mod_range(inByTr_new_unique,10,50)
inByTr_new_unique_ls<-inByTr_new_unique[order(start(inByTr_new_unique),decreasing=F)]
inByTr_new_unique_ls<-split(inByTr_new_unique,inByTr_new_unique$gene)
cdsByTr_unique <- GenomicRanges::shift(cdsByTr_unique,shift=cdsByTr_unique$phase)# Adjust the length of each cds region so that the first base is in frame of the gene
cdsByTr_unique<-cdsByTr_unique[width(cdsByTr_unique)>=36]
cdsByTr_unique<-mod_range(cdsByTr_unique,6,6)
cdsByTr_unique_ls<-split(cdsByTr_unique,cdsByTr_unique$trc)
cdsByTr_unique_gene_ls<-split(cdsByTr_unique,cdsByTr_unique$gene)
largeintron<-unique(inByTr_new_unique[width(inByTr_new_unique)>=200]$gene)
two_intron<-inByTr_new_unique_ls[elementNROWS(inByTr_new_unique_ls)>=2]
near_intron<-data.frame(gene=unique(unlist(two_intron)$gene),neardis=1)
for(i in 1:length(two_intron)){
  two_intron[[i]]<-two_intron[[i]][order(start(two_intron[[i]]),decreasing=F)]
  index<-which.min(width(gaps(two_intron[[i]]))[-1])
  two_intron[[i]]<-two_intron[[i]][c(index,index+1),]
  two_intron[[i]]<-two_intron[[i]][order(start(two_intron[[i]]),decreasing=F)]
  near_intron[i,2]<-start(two_intron[[i]][2])-end(two_intron[[i]][1])
}
intron_cds_intron_gene<-filter(near_intron,neardis<=800)
include_select<-function(x,y){
  s<-names(x)%in%y
  ss<-grep(F,s)
  x<-x[-ss]
  return(x)
}
nearest_intron<-function(x){
  for(i in 1:length(x)){
    x[[i]]<-x[[i]][order(start(x[[i]]),decreasing=F)]
    index<-which.min(width(gaps(x[[i]]))[-1])
    x[[i]]<-x[[i]][c(index,index+1),]
    x[[i]]<-x[[i]][order(start(x[[i]]),decreasing=F)]
  }
  return(x)
}
intron_cds_intron_gene_ls<-include_select(inByTr_new_unique_ls,intron_cds_intron_gene$gene)
intron_cds_intron_gene_ls<-nearest_intron(intron_cds_intron_gene_ls)


# 2.2 codon_score &codon replace
aa_codon_score<-data.frame(GENETIC_CODE)
aa_codon_score<-rownames_to_column(aa_codon_score, var = "rowname")
names(aa_codon_score)[1]<-'codon'
names(aa_codon_score)[2]<-'AA'
aa_codon_score$score <-1
for(i in 1:64){
  if(aa_codon_score[i,2]=='*'|aa_codon_score[i,2]=='M'|aa_codon_score[i,2]=='W')
    aa_codon_score[i,3]<-0
  else if(aa_codon_score[i,2]=='R'|aa_codon_score[i,2]=='L')
    aa_codon_score[i,3]<-2
  else if(aa_codon_score[i,2]=='S')
    aa_codon_score[i,3]<-3
  else 
    aa_codon_score[i,3]<-1
}
aa_codon_replace<-aa_codon_score
aa_codon_replace$AA <- as.character(aa_codon_replace$AA)
aa_codon_replace <- mapply(function(x, y){
  incs <- aa_codon_replace[aa_codon_replace$AA == y,]$codon
  incs <- incs[!grepl(x, incs)]
  return(incs)
},aa_codon_replace$codon, aa_codon_replace$AA)
aa_codon_replace['AGG']<-list(c('CGA','CGT','CGC')) # find the synonymous codon with the largest difference in each codon sequence
aa_codon_replace['AGA']<-list(c('CGG','CGT','CGC'))
aa_codon_replace['CGG']<-list('AGA')
aa_codon_replace['CGA']<-list('AGG')
aa_codon_replace['CGT']<-list(c('AGA','AGG'))
aa_codon_replace['CGC']<-list(c('AGA','AGG'))
aa_codon_replace['TTG']<-list(c('CTA','CTC','CTT'))
aa_codon_replace['TTA']<-list(c('CTG','CTC','CTT'))
aa_codon_replace['CTA']<-list('TTG')
aa_codon_replace['CTG']<-list('TTA')
aa_codon_replace['CTC']<-list(c('TTG','TTA'))
aa_codon_replace['CTT']<-list(c('TTG','TTA'))
aa_codon_replace['AGT']<-list(c('TCG','TCA','TCC'))
aa_codon_replace['AGC']<-list(c('TCG','TCA','TCT'))
aa_codon_replace['TCA']<-list(c('AGC','AGT'))
aa_codon_replace['TCG']<-list(c('AGC','AGT'))
aa_codon_replace['TCT']<-list('AGC')
aa_codon_replace['TCC']<-list('AGT')
aa_codon_replace['ATG']<-list('ATG')
aa_codon_replace['TGG']<-list('TGG')

replace_sy<-function(x){
  sy<-sample(aa_codon_replace[x][[1]],1)
  return(sy)}
replace_synonym<-function(x){
  seq <- substring(x,seq(1,nchar(x)-2,3),seq(3,nchar(x),3))
  seq<-sapply(seq,replace_sy,USE.NAMES = F)
  return(paste(seq,collapse = ''))
}

makebed<-function(x,address){
  for(i in 1:length(x)){
    name<-names(x[i])
    export.bed(x[[i]],paste0(address,'/',name,'.bed'))
  }
}

importbed<-function(gr,address){
  tm <- function(x){
    as.numeric(system(paste0("./primer3/src/oligotm ",x),intern = T))
  }
  shuffle<-function(x,length){
    x<-as.character(x)
    seq<-substring(x,seq(1,nchar(x),1),seq(1,nchar(x),1))
    sf<-paste(sample(seq,size=length,replace=F),collapse = '')
    sff<-unlist(sf)
    return(sff)
  }
  reversseq<-function(x,lian,str){
    revers<-x
    for(j in seq_len(length(x))){
      strand<-as.character(lian)
      if(strand[j]!=str){
        revers[j]<-reverseComplement(revers[j])
      }
      else{revers[j]<-revers[j]}
    }
    return(revers)
  }
  
  for(i in 1:length(gr)){
    gr[[i]]<-gr[[i]][1]
  }
  chr<-unique(seqnames(unlist(gr)))
  gr<-as.vector(strand(unlist(gr)))
  file<-grep("win.bed",dir(address), value = T)
  print(file)
  genename<-file
  for(i in 1:length(file)){
    genename[i]<-strsplit(file[i],'\\_')[[1]][1]
  }
  SS_5 <- "AGGT[AG]AGT"
  BP<-'CT[AG]A[CT]'
  SS_3 <-'[CT]AGG'
  out<-GRangesList()
  test<-function(i){
    x<-import(paste(address, file[i], sep = "/"))
    x$gene<-genename[i]
    print(genename[i])
    strand(x)<-gr[i]
    x<-x[width(x)==24]
    x$seq<-getSeq(fasta_genome,x)
    x$BP<-grepl(BP, as.character(x$seq))
    x<-x[x$BP==F]
    x$tm<-sapply(x$seq,tm)
    x<-x[x$tm>=58]
    x<-x[x$tm<=65]
    if(length(x)<=2){
      x_3<-GRanges()
    }
    else{
      primer<-as.data.frame(matrix(1:(dim(combn(1:length(x),2))[2]*5),ncol = 5))
      names(primer)[1]<-'i_start'
      names(primer)[2]<-'i_end'
      names(primer)[3]<-'j_start'
      names(primer)[4]<-'j_end'
      names(primer)[5]<-'product_length'
      for(z in 1:(length(x)-1)){
        for(j in (z+1):length(x)){
          order<-length(x)*(z-1)-(z*(z+1)/2-z)+(j-z)
          print(order)
          primer[order,1]<-start(x[z])
          primer[order,2]<-end(x[z])
          primer[order,3]<-start(x[j])
          primer[order,4]<-end(x[j])
          sp<-x[c(z,j),]
          primer[order,5]<-length<-max(end(sp))-min(start(sp))
        }
      }
      primer$product_length<-as.numeric(primer$product_length)
      x_1<-primer[primer$product_length<=1000,]
      x_1<-x_1[x_1$product_length>=100,]
      if(nrow(x_1)==0){
        x_3<-GRanges()
      }
      else{
        x_11<-GRanges(seqnames = rep(chr,nrow(x_1)),
                      ranges = IRanges(start=x_1$i_start,end = x_1$i_end),
                      strand = '+',
                      num = c(1:nrow(x_1)))
        x_12<-GRanges(seqnames = rep(chr,nrow(x_1)),
                      ranges = IRanges(x_1$j_start,end = x_1$j_end),
                      strand = '-',
                      num = c(1:nrow(x_1)))
        
        x_2<-c(x_11,x_12)
        x_2$wt<-getSeq(fasta_genome,x_2)
        x_2$sf<-DNAStringSet(sapply(x_2$wt,function(x)shuffle(x,length = 24),USE.NAMES=F))
        x_2$rev<-reversseq(x_2$sf,strand(x_2),gr[i])
        x_2$SS5<-grepl(SS_5, as.character(x_2$rev))
        x_2$BP<-grepl(BP, as.character(x_2$rev))
        x_2$SS3<-grepl(SS_3, as.character(x_2$rev))
        x_2<-x_2[x_2$BP==F]
        x_2<-x_2[x_2$SS5==F]
        x_2<-x_2[x_2$SS3==F]
        if(length(x_2)==0){
          x_3<-GRanges()
        }
        else{
          x_2$sf_Tm<-sapply(x_2$sf,tm)
          x_2<-x_2[x_2$sf_Tm>=58.0]
          x_2<-x_2[x_2$sf_Tm<=65.0]
          if(length(x_2)==0){
            x_3<-GRanges()
          }
          else{
            x_2list<-split(x_2,x_2$num)
            x_2list<-x_2list[lengths(x_2list,use.names = F)==2]
            x_3<-unlist(x_2list)
            if(length(x_3)==0){
              x_3<-GRanges()
            }
            else{
              x_3$gene<-genename[i]
            }
          }
        }
      }
    }
    return(x_3)
  }
  for(i in 1:length(file)){
    x_3<-try(test(i),silent = F)
    if('try-error' %in% class(x_3))           
    {
      x_3<-GRanges()                            
    }
    out[[i]]<-x_3
  }
  names(out)<-genename
  return(out)
}
makedf<-function(x){
  len<-elementNROWS(x)
  print(len)
  total<-list()
  for(i in 1:length(x)){
    if(len[i]>10){
      select<-sample(1:(len[i]/2),size=5,replace = F)
      select_1<-select*2
      select_2<-select*2-1
      select<-c(select_1,select_2)
      x[[i]]<-x[[i]][select,]
    }
    else{x[[i]]<-x[[i]]}
    df_wt<-data.frame(num=as.vector(x[[i]]$num),seq=as.vector(x[[i]]$wt))
    df_wt$seq<-as.character(df_wt$seq)
    df_wt<-df_wt%>%
      group_by(num)%>%
      summarize(count=n(),pair=list(seq))
    df_wt$gene<-names(x[i])
    df_wt$type<-'wt'    
    df_syn<-data.frame(num=as.vector(x[[i]]$num),seq=as.vector(x[[i]]$sf))
    df_syn$seq<-as.character(df_syn$seq)
    df_syn<-df_syn%>%
      group_by(num)%>%
      summarize(count=n(),pair=list(seq))
    df_syn$gene<-names(x[i])
    df_syn$type<-'syn'    
    df<-rbind(df_wt,df_syn)
    total[[i]]<-df
  }
  total<-do.call(rbind,total)
  return(total)  
}

fun <- function(x){
  tmp <- AmplifyDNA(x,fasta_genome,annealingTemp=58, P=4e-7,maxProductSize=3000)
  return(tmp)}#in silico pcr

makebed(intron_cds_intron_gene_ls,'./intron_exon_intron/')
# 2.3 Cut the sequence of the design range into small fragments of ~24bp
# operate in Linux using Shell
# 2.4 Design PCRmarks and verify them in-silico (in introns)
in_in_L_1<-importbed(intron_cds_intron_gene_ls,'./intron_exon_intron/')
in_in_L_2<-makedf(in_in_L_1)
in_in_L_3<- do.call('list',parallel::mclapply(in_in_L_2$pair,function(x){fun(x)},mc.cores = 90))

select_pcrmarker<-function(df,ls){
  df$product<-lengths(ls)
  for(i in 1:nrow(df)){
    if(df$product[i]==1){
      df$width[i]<-c(width(ls[[i]]))
      df$names[i]<-c(names(ls[[i]]))
      df$possible[i]<-strsplit(df$names[i],'\\%')[[1]][1]
      df$primer_F[i]<-substr(strsplit(df$names[i],'\\ ')[[1]][2],2,2)
      df$primer_R[i]<-substr(strsplit(df$names[i],'\\ ')[[1]][4],1,1)
    }else{df$width[i]<-0
    df$names[i]<-0
    df$possible[i]<-0
    df$primer_F[i]<-0
    df$primer_R[i]<-0}
  }
  df$possible<-as.numeric(df$possible)
  df_1<-filter(df %>% filter(type=='wt'),product==1,width<1000,width>100,primer_F!=primer_R,possible>95)
  df_2<-filter(df %>% filter(type=='syn'),product==0)
  df_3<-rbind(df_1,df_2)
  
  df_3<-df_3%>%
    group_by(num,gene)%>%
    summarize(count=n(),pair=list(pair))
  df_3<-filter(df_3,count==2)
  df_3<-df_3%>%
    group_by(gene)%>%
    summarize(number=list(num))
  return(df_3)
}
in_in_L_4<-select_pcrmarker(in_in_L_2,in_in_L_3)

in_in_rest_gene<-setdiff(names(intron_cds_intron_gene_ls),in_in_L_4$gene)

select_pcrtar_from_grangelist<-function(df,grls){
  df_1<-GRangesList()
  for(i in 1:nrow(df)){
    df_1[[i]]<-grls[[df$gene[i]]][which(grls[[df$gene[i]]]$num %in% df$number[[i]]==T)]
  }
  return(df_1)
}
in_in_L_5<-select_pcrtar_from_grangelist(in_in_L_4,in_in_L_1)

# 2.5 Design PCRmarks and verify them in-silico (in intron-CDS-intron)
largeintron<-setdiff(largeintron,genome_gene$gene)
largeintron<-c(largeintron,in_in_rest_gene)
intron_intron_gene_ls<-include_select(inByTr_new_unique_ls,largeintron)
gene_in_arm_inin<-unlist(intron_intron_gene_ls)

for(i in 1:length(intron_intron_gene_ls)){
  intron_intron_gene_ls[[i]]<-intron_intron_gene_ls[[i]][which.max(width(intron_intron_gene_ls[[i]]))]
}
intron_intron_gene_ls
makebed(intron_intron_gene_ls,'./intron_intron/')
# Cut the sequence of the design range into small fragments of ~24bp
# operate in Linux using Shell
in_L_1<-importbed(intron_intron_gene_ls,'./intron_intron/')
in_L_2<-makedf(in_L_1)
in_L_3<- do.call('list',parallel::mclapply(in_L_2$pair,function(x){fun(x)},mc.cores = 5))
in_L_4<-select_pcrmarker(in_L_2,in_L_3)
in_rest_gene<-setdiff(names(intron_intron_gene_ls),in_L_4$gene)
in_L_5<-select_pcrtar_from_grangelist(in_L_4,in_L_1)


# 2.6 Design PCRmarks and verify them in-silico (in CDS)
cds_tag<-setdiff(names(cdsByTr_unique_gene_ls),union(genome_gene$gene,gene_in_arm_inin$gene))###设计cds pcrmarker的gene
cds_tag<-c(cds_tag,in_rest_gene)
cds_cds_gene_ls<-include_select(cdsByTr_unique_gene_ls,cds_tag)
cds_cds_gene_ls<-split(unlist(cds_cds_gene_ls),unlist(cds_cds_gene_ls)$trc)
importbed_exon<-function(gr,address,scores,tm_l,tm_h){
  tm<- function(x){
    as.numeric(system(paste0("/Dell/Dell6/zhangyl/software/primer3/src/oligotm ",x),intern = T))
  }
  
  reversseq<-function(x,lian,str){
    revers<-x
    for(j in seq_len(length(x))){
      strand<-as.character(lian)
      if(strand[j]!=str){
        revers[j]<-reverseComplement(revers[j])
      }
      else{revers[j]<-revers[j]}
    }
    return(revers)
  }
  getscore<-function(x){
    seq <- substring(x,seq(1,nchar(x)-2,3),seq(3,nchar(x),3))
    s<-data.frame(seq)
    names(s)[1]<-'codon'
    s<-merge(s,aa_codon_score,by='codon',all.x=T)
    ss<-sum(s$score)
    return(ss)
  }
  start_stop_codon<-function(x){
    seq <- substring(x,seq(1,nchar(x)-2,3),seq(3,nchar(x),3))
    n<-length(which(seq=='ATG'|seq=='TGA'|seq=='TAA'|seq=='TAG'|seq=='TGG'))
    return(n)
  }
  
  for(i in 1:length(gr)){
    gr[[i]]<-gr[[i]][1]
  }
  chr<-unique(seqnames(unlist(gr)))
  gr<-as.vector(strand(unlist(gr)))
  file<-grep("win.bed",dir(address), value = T)
  genename<-file
  for(i in 1:length(file)){
    genename[i]<-strsplit(file[i],'\\_')[[1]][1]
  }
  out<-GRangesList()
  test<-function(i){
    x<-import(paste(address, file[i], sep = "/"))
    x$gene<-genename[i]
    print(genename[i])
    strand(x)<-gr[i]
    x<-x[width(x)==24]
    x$seq<-getSeq(fasta_genome,x)
    x$score<-sapply(as.character(x$seq),getscore)
    x<-x[x$score>=scores]
    print(length(x))
    if(length(x)>=30){
      x<-x[sample(1:length(x),size=30,replace = F),]
    }else if(length(x)<2){
      x_3<-GRanges()
    }else{x<-x}
    x$tm<-sapply(x$seq,tm)
    x<-x[x$tm>=tm_l]
    x<-x[x$tm<=tm_h]
    if(length(x)<=2){
      x_3<-GRanges()
    }
    else{
      x<-x[order(start(x),decreasing=F)]
      print(x)
      primer<-as.data.frame(matrix(1:(dim(combn(1:length(x),2))[2]*5),ncol = 5))
      names(primer)[1]<-'i_start'
      names(primer)[2]<-'i_end'
      names(primer)[3]<-'j_start'
      names(primer)[4]<-'j_end'
      names(primer)[5]<-'product_length'
      for(z in 1:(length(x)-1)){
        for(j in (z+1):length(x)){
          order<-length(x)*(z-1)-(z*(z+1)/2-z)+(j-z)
          print(order)
          primer[order,1]<-start(x[z])
          primer[order,2]<-end(x[z])
          primer[order,3]<-start(x[j])
          primer[order,4]<-end(x[j])
          sp<-x[c(z,j),]
          primer[order,5]<-length<-max(end(sp))-min(start(sp))
        }
      }
      primer$product_length<-as.numeric(primer$product_length)
      x_1<-primer[primer$product_length<=1000,]
      x_1<-x_1[x_1$product_length>=100,]
      if(nrow(x_1)==0){
        x_3<-GRanges()
      }
      else{
        x_11<-GRanges(seqnames = rep(chr,nrow(x_1)),
                      ranges = IRanges(start=x_1$i_start,end = x_1$i_end),
                      strand = '+',
                      num = c(1:nrow(x_1)))
        x_12<-GRanges(seqnames = rep(chr,nrow(x_1)),
                      ranges = IRanges(x_1$j_start,end = x_1$j_end),
                      strand = '-',
                      num = c(1:nrow(x_1)))
        
        x_2<-c(x_11,x_12)
        x_2$wt<-getSeq(fasta_genome,x_2)
        x_2$sf<-x_2$wt
        x_2$sf<-reversseq(x_2$sf,strand(x_2),gr[i])
        x_2$sf<-as.vector(x_2$sf)
        x_2$sf<-DNAStringSet(sapply(x_2$sf,replace_synonym))
        x_2$sf<-reversseq(x_2$sf,strand(x_2),gr[i])
        x_2$sf_Tm<-sapply(x_2$sf,tm)
        x_2<-x_2[x_2$sf_Tm>=tm_l]
        x_2<-x_2[x_2$sf_Tm<=tm_h]
        if(length(x_2==0)){
          x_3<-GRanges()
        }
        else{
          x_2list<-split(x_2,x_2$num)
          x_2list<-x_2list[lengths(x_2list,use.names = F)==2]
          if(length(x_2list)==0){
            x_3<-GRanges()
          }else{
            x_3<-unlist(x_2list)
            x_3$gene<-genename[i]
          }
        }
      }
    }
    return(x_3)
  }
  for(i in 1:length(file)){
    x_3<-try(test(i),silent = F)
    if('try-error' %in% class(x_3))            
    {
      x_3<-GRanges()                           
    }
    out[[i]]<-x_3
  }
  names(out)<-genename
  return(out)
}

makebed(cds_cds_gene_ls,'./cds_cds/')
# Cut the sequence of the design range into small fragments of ~24bp
# operate in Linux using Shell
cds_L_1<-importbed_exon(cds_cds_gene_ls,'./cds_cds/',
                           scores=11,tm_l=58,tm_h=65)
cds_L_2<-makedf(cds_L_1)
cds_L_3<- do.call('list',parallel::mclapply(cds_L_2$pair,function(x){fun(x)},mc.cores = 70))
cds_L_4<-select_pcrmarker(cds_L_2,cds_L_3)
cds_L_5<-select_pcrtar_from_grangelist(cds_L_4,cds_L_1)

# 2.7 Integrate all PCRmarks and replace sequences
PCRmarks<-c(unlist(in_in_L_5)[,c(1,2,3,9)],unlist(in_L_5)[,c(1,2,3,9)],unlist(cds_L_5)[,c(1,2,3,5)])
for(i in 1:length(PCRmarks)){
  PCRmarks$gene[i]<-strsplit(PCRmarks$gene[i],'\\.')[[1]][1]
}
PCRmarks<-split(PCRmarks,PCRmarks$gene)
manage_total_pcrmarker<-function(x){
  x<-unlist(x)
  x$wt_tm<-sapply(x$wt,tm)
  x$sf_tm<-sapply(x$sf,tm)
  x<-split(x,x$gene)
  return(x)
}# add Tm
PCRmarks<-manage_total_pcrmarker(PCRmarks)

manage_total_pcrmarker_2<-function(x){
  len<-lengths(x)
  for(i in 1:length(x)){
    if(len[i]>2){
      select<-sample(1:(len[i]/2),size=1,replace = F)
      select_1<-select*2
      select_2<-select*2-1
      select<-c(select_1,select_2)
      x[[i]]<-x[[i]][select,]
    }else{x[[i]]<-x[[i]]}
  }
  return(x)
}# randomly select one PCRmark for each gene
PCRmarks_2<-manage_total_pcrmarker_2(PCRmarks)
PCRmarks_2<-unlist(PCRmarks_2)

SYN_sequence<-fasta_genome
PCRmarks_2$sf<-as.character(PCRmarks_2$sf)
pcrmarker_swith<-function(pcrmarker_range,chromosome_fasta,wid){
  start_index <- start(pcrmarker_range)
  strand_index <- as.character(strand(pcrmarker_range))
  for(i in seq_len(length(pcrmarker_range))){
    tag<-as.character(DNAString(pcrmarker_range[i]$sf))
    anti_tag<-as.character(reverseComplement(DNAString(pcrmarker_range[i]$sf)))
    if(strand_index[i] == '+'){
      subseq(chromosome_fasta, start = start_index[i], width = wid) <-tag
    } 
    else{
      subseq(chromosome_fasta, start = start_index[i], width = wid) <-anti_tag
    }
  }
  return(chromosome_fasta)
}
SYN_sequence<-pcrmarker_swith(PCRmarks_2,SYN_sequence,24)

PCRmarks_2$type<-'pcrmarker_pair'
PCRmarks_2$wt<-as.character(PCRmarks_2$wt)
PCRmarks_2$sf<-as.character(PCRmarks_2$sf)
names(PCRmarks_2)<-as.character(PCRmarks_2$wt)
writeXStringSet(SYN_sequence,filepath="./SYN_sequence.fasta",format="fasta")
