library(optparse)
option_list<-list(
  make_option(c('-g','--genome_file'), type='character', default = FALSE, help='genome file in fasta format'),
  make_option(c('-f','--gff_file'), type='character', default = FALSE, help='genome annotation file in gff format'),
  make_option(c('-c','--cen_file'), type='character', default = FALSE, help='centromere file in csv format'),
  make_option(c('-n','--chr'), type='integer', default = FALSE, help='chromosome number (1, 2...)'),
  make_option(c('-a','--arm'), type='character', default = FALSE, help='specify chromosome arm to design (L or R)'),
  make_option(c('-t','--thread'), type='integer', default = 1, help='thread count'),
  make_option(c('-m','--tm_location'), type='character', default = FALSE, help='location of the oligotm function in Primer3'),
  make_option(c('-o','--output_path'), type='character', default = FALSE, help='the output path for the generated PCRmark information')
)
opt_parser = OptionParser(usage='uasge: %prog -g <Genome_file> -f <gff_file> -c <Cen_file> -n <Chr_number(INT)> -a <L/R> -t <thread> -m <tm_location> -o <output_path>',
                          option_list=option_list) 
opt = parse_args(opt_parser) 
geno<-opt$genome_file
gff<-opt$gff_file
centromere<-opt$cen_file
chr_num<-opt$chr
arm<-opt$arm
thr<-opt$thread
path<-opt$output_path
oligotm<-opt$tm_location

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

chr_name<-c(paste(0,1:9,sep=''),10:26)[chr_num]
chr<-paste('Chr',chr_name,sep='')

# 1. Import and process sequence information
fasta_genome<-readDNAStringSet(geno)# wildtype genome sequence 
annotation_genome<- makeTxDbFromGFF(gff)# wildtype genome annotation
genome_gene<-genes(annotation_genome)
seqlengths(genome_gene)<-as.vector(width(fasta_genome))
cen<-fread(centromere)
# 2. Introduction of PCRmarks
# 2.1 Distinguish between intron and CDS, and standardize the design range according to pcrmarks requirements
inByTr_new<-intronsByTranscript(annotation_genome,use.names = TRUE)# introns grouped by transcript
inByTr_new<-inByTr_new[lengths(inByTr_new)!=0]
add_trc_column <- function(gr_list) {
  trc_values <- names(gr_list)
  gr_list_with_trc <- lapply(seq_along(gr_list), function(i) {
    mcols(gr_list[[i]])$trc <- trc_values[i]
    return(gr_list[[i]])
  })
  return(gr_list_with_trc)
}
inByTr_new<-add_trc_column(inByTr_new)
inByTr_new_ul<-do.call(c,inByTr_new)
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
}#narrow down the range of intron/cds
inByTr_new_unique<-mod_range(inByTr_new_unique,10,50)
inByTr_new_unique_ls<-inByTr_new_unique[order(start(inByTr_new_unique),decreasing=F)]
inByTr_new_unique_ls<-split(inByTr_new_unique,inByTr_new_unique$gene)
cdsByTr_unique <- GenomicRanges::shift(cdsByTr_unique,shift=cdsByTr_unique$phase)# Adjust the length of each cds region so that the first base is in frame of the gene
cdsByTr_unique<-cdsByTr_unique[width(cdsByTr_unique)>=36]
cdsByTr_unique<-mod_range(cdsByTr_unique,6,6)

#Filter out introns that can design PCRMarks
largeintron<-inByTr_new_unique[width(inByTr_new_unique)>=200]
largeintron<-split(largeintron,largeintron$gene)
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

# 2.2 Define codon_score & codon replace
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
aa_codon_replace['TAA']<-list('TAA')
aa_codon_replace['TGA']<-list('TGA')
aa_codon_replace['TAG']<-list('TAG')


# 2.3 Design PCRmarks and verify them in-silico (on two introns across a cds, i.e. intron-cds-intron PCRmarks)
#Determine the scope of the design, i.e. chromosomes and arms
intron_cds_intron_gene_ls_chr<-unlist(intron_cds_intron_gene_ls)[seqnames(unlist(intron_cds_intron_gene_ls))==chr]
intron_cds_intron_gene_ls_chr<-split(intron_cds_intron_gene_ls_chr,intron_cds_intron_gene_ls_chr$gene)
#if L or R arm
if(arm=='L'){
  gene_inin_arm<-restrict(unlist(intron_cds_intron_gene_ls_chr),1,cen$Start[chr_num])
}else{gene_inin_arm<-restrict(unlist(intron_cds_intron_gene_ls_chr),cen$End[chr_num],cen$`Chromosome length (bp)`[chr_num])}
intron_cds_intron_gene_ls_chr_arm<-split(gene_inin_arm,gene_inin_arm$gene)
#Cut the sequence of the design range into small fragments of ~24bp (step=12bp)
make_windows_grl<-function(grl, window_size, step_size){
  new_grl<-lapply(grl,function(gr){
    gr_split <- split(gr, seq_along(gr))
    for(i in 1:length(gr)){
      new_starts <- start(gr[i]) + seq(0, width(gr[i]) - window_size, by = step_size)
      gr_split[[i]] <- GRanges(
        seqnames = seqnames(gr[i]),
        ranges = IRanges(start = new_starts, end = new_starts + window_size - 1)
      )}
    new_gr<-unlist(gr_split)
    return(new_gr)
  })
  new_grl<-do.call(c,new_grl)
  return(new_grl)
}
in_in_cut<-make_windows_grl(intron_cds_intron_gene_ls_chr_arm,24,12)
#List all optional PCRmarks for subsequent in-silico PCR
pcrmark_intron<-function(gr,gr_cut){
  tm <- function(x){
    as.numeric(system(paste0(oligotm,' ',x),intern = T))
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
  genename<-names(gr_cut)
  SS_5 <- "AGGT[AG]AGT"
  BP<-'CT[AG]A[CT]'
  SS_3 <-'[CT]AGG'
  out<-GRangesList()
  test<-function(i){
    x<-gr_cut[[i]]
    x$gene<-genename[i]
    strand(x)<-gr[i]
    x<-x[width(x)==24]
    if(length(x)==0){
      x_3<-GRanges()
    }
    else{
      x$seq<-getSeq(fasta_genome,x)
      x$BP<-grepl(BP, as.character(x$seq))
      x<-x[x$BP==F]
      x$tm<-sapply(x$seq,tm)
      x<-x[x$tm>=57]
      x<-x[x$tm<=68]
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
            x_2<-x_2[x_2$sf_Tm>=57.0]
            x_2<-x_2[x_2$sf_Tm<=68.0]
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
    }
    return(x_3)
  }
  for(i in 1:length(genename)){
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
in_in_pcrmark_1<-pcrmark_intron(intron_cds_intron_gene_ls_chr_arm,in_in_cut)
makedf<-function(x){
  len<-elementNROWS(x)
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
insilicopcr <- function(x){
  tmp <- AmplifyDNA(x,fasta_genome,annealingTemp=58, P=4e-7,maxProductSize=1000)
  return(tmp)}
select_pcrmark<-function(df,ls){
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
select_pcrmark_from_grangelist<-function(df,grls){
  df_1<-GRangesList()
  for(i in 1:nrow(df)){
    df_1[[i]]<-grls[[df$gene[i]]][which(grls[[df$gene[i]]]$num %in% df$number[[i]]==T)]
  }
  return(df_1)
}
if(max(lengths(in_in_pcrmark_1))==0){
  in_in_rest_gene<-names(intron_cds_intron_gene_ls_chr_arm)
  in_in_pcrmark_5<-in_in_pcrmark_1[1]
  in_in_pcrmark_5[[1]] <-in_in_pcrmark_5[[1]][-c(1:length(in_in_pcrmark_5[[1]]))]
}else{
  in_in_pcrmark_2<-makedf(in_in_pcrmark_1)
  #in-silico PCR
  in_in_pcrmark_3<- do.call('list',parallel::mclapply(in_in_pcrmark_2$pair,function(x){insilicopcr(x)},mc.cores = thr))
  in_in_pcrmark_4<-select_pcrmark(in_in_pcrmark_2,in_in_pcrmark_3)
  if(nrow(in_in_pcrmark_4)==0){
    in_in_rest_gene<-names(intron_cds_intron_gene_ls_chr_arm)
    in_in_pcrmark_5<-in_in_pcrmark_1[1]
    in_in_pcrmark_5[[1]] <-in_in_pcrmark_5[[1]][-c(1:length(in_in_pcrmark_5[[1]]))]
  }else{in_in_pcrmark_5<-select_pcrmark_from_grangelist(in_in_pcrmark_4,in_in_pcrmark_1)
  in_in_rest_gene<-setdiff(names(intron_cds_intron_gene_ls_chr_arm),in_in_pcrmark_4$gene)}
}

# 2.4 Design PCRmarks and verify them in-silico (on one intron, i.e. intron-intron PCRmarks)
#Determine the scope of the design, i.e. chromosomes and arms
largeintron_chr<-unlist(largeintron)[seqnames(unlist(largeintron))==chr]
largeintron_chr<-split(largeintron_chr,largeintron_chr$gene)
largeintron_chr<-setdiff(names(largeintron_chr),gene_inin_arm$gene)
largeintron_chr<-c(largeintron_chr,in_in_rest_gene)
intron_intron_gene_ls_chr<-include_select(inByTr_new_unique_ls,largeintron_chr)
#if L or R arm
if(arm=='L'){
  gene_in_arm<-restrict(unlist(intron_intron_gene_ls_chr),1,cen$Start[chr_num])
}else{gene_in_arm<-restrict(unlist(intron_intron_gene_ls_chr),cen$End[chr_num],cen$`Chromosome length (bp)`[chr_num])}
intron_intron_gene_ls_chr_arm<-split(gene_in_arm,gene_in_arm$gene)
for(i in 1:length(intron_intron_gene_ls_chr_arm)){
  intron_intron_gene_ls_chr_arm[[i]]<-intron_intron_gene_ls_chr_arm[[i]][which.max(width(intron_intron_gene_ls_chr_arm[[i]]))]
}
#Cut the sequence of the design range into small fragments of ~24bp (step=12bp)
in_cut<-make_windows_grl(intron_intron_gene_ls_chr_arm,24,12)
#List all optional PCRmarks for subsequent in-silico PCR
in_pcrmark_1<-pcrmark_intron(intron_intron_gene_ls_chr_arm,in_cut)
if(max(lengths(in_pcrmark_1))==0){
  in_rest_gene<-names(intron_intron_gene_ls_chr_arm)
  in_pcrmark_5<-in_pcrmark_1[1]
  in_pcrmark_5[[1]] <-in_pcrmark_5[[1]][-c(1:length(in_pcrmark_5[[1]]))]
}else{
  in_pcrmark_2<-makedf(in_pcrmark_1)
  #in-silico PCR
  in_pcrmark_3<- do.call('list',parallel::mclapply(in_pcrmark_2$pair,function(x){insilicopcr(x)},mc.cores = thr))
  in_pcrmark_4<-select_pcrmark(in_pcrmark_2,in_pcrmark_3)
  if(nrow(in_pcrmark_4)==0){
    in_rest_gene<-names(intron_intron_gene_ls_chr_arm)
    in_pcrmark_5<-in_pcrmark_1[1]
    in_pcrmark_5[[1]] <-in_pcrmark_5[[1]][-c(1:length(in_pcrmark_5[[1]]))]
  }else{in_pcrmark_5<-select_pcrmark_from_grangelist(in_pcrmark_4,in_pcrmark_1)
  in_rest_gene<-setdiff(names(intron_intron_gene_ls_chr_arm),in_pcrmark_4$gene)}
}


# 2.5 Design PCRmarks and verify them in-silico (on cds, i.e. cds PCRmarks)
#Determine the scope of the design, i.e. chromosomes and arms
cdsByTr_unique_chr<-cdsByTr_unique[seqnames(cdsByTr_unique)==chr]
cdsByTr_unique_ls_chr<-split(cdsByTr_unique_chr,cdsByTr_unique_chr$gene)
cds_tag<-setdiff(names(cdsByTr_unique_ls_chr),union(gene_inin_arm$gene,gene_in_arm$gene))###设计cds pcrmarker的gene
cds_tag<-c(cds_tag,in_rest_gene)
cds_gene_ls<-include_select(cdsByTr_unique_ls_chr,cds_tag)
#if L or R arm
if(arm=='L'){
  gene_cds_arm<-restrict(unlist(cds_gene_ls),1,cen$Start[chr_num])
}else{gene_cds_arm<-restrict(unlist(cds_gene_ls),cen$End[chr_num],cen$`Chromosome length (bp)`[chr_num])}
cds_gene_ls_chr_arm<-split(gene_cds_arm,gene_cds_arm$gene)
#Cut the sequence of the design range into small fragments of ~24bp (step=12bp)
cds_cut<-make_windows_grl(cds_gene_ls_chr_arm,24,12)
#List all optional PCRmarks for subsequent in-silico PCR
pcrmark_cds<-function(gr,gr_cut,scores){
  tm <- function(x){
    as.numeric(system(paste0(oligotm,' ',x),intern = T))
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
  replace_sy<-function(x){
    sy<-sample(aa_codon_replace[x][[1]],1)
    return(sy)}
  replace_synonym<-function(x){
    seq <- substring(x,seq(1,nchar(x)-2,3),seq(3,nchar(x),3))
    seq<-sapply(seq,replace_sy,USE.NAMES = F)
    return(paste(seq,collapse = ''))}
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
  genename<-names(gr_cut)
  out<-GRangesList()
  test<-function(i){
    x<-gr_cut[[i]]
    x$gene<-genename[i]
    strand(x)<-gr[i]
    x<-x[width(x)==24]
    x$seq<-getSeq(fasta_genome,x)
    x$score<-sapply(as.character(x$seq),getscore)
    x<-x[x$score>=scores]
    if(length(x)>=30){
      x<-x[sample(1:length(x),size=30,replace = F),]
    }else if(length(x)<2){
      x_3<-GRanges()
    }else{x<-x}
    x$tm<-sapply(x$seq,tm)
    x<-x[x$tm>=57]
    x<-x[x$tm<=68]
    if(length(x)<=2){
      x_3<-GRanges()
    }
    else{
      x<-x[order(start(x),decreasing=F)]
      primer<-as.data.frame(matrix(1:(dim(combn(1:length(x),2))[2]*5),ncol = 5))
      names(primer)[1]<-'i_start'
      names(primer)[2]<-'i_end'
      names(primer)[3]<-'j_start'
      names(primer)[4]<-'j_end'
      names(primer)[5]<-'product_length'
      for(z in 1:(length(x)-1)){
        for(j in (z+1):length(x)){
          order<-length(x)*(z-1)-(z*(z+1)/2-z)+(j-z)
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
        x_2<-x_2[x_2$sf_Tm>=57]
        x_2<-x_2[x_2$sf_Tm<=68]
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
  for(i in 1:length(genename)){
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
cds_pcrmark_1<-pcrmark_cds(cds_gene_ls_chr_arm,cds_cut,scores=11)
if(max(lengths(cds_pcrmark_1))==0){
  cds_rest_gene<-names(cds_gene_ls_chr_arm)
  cds_pcrmark_5<-cds_pcrmark_1[1]
  cds_pcrmark_5[[1]] <-cds_pcrmark_5[[1]][-c(1:length(cds_pcrmark_5[[1]]))]
}else{
  cds_pcrmark_2<-makedf(cds_pcrmark_1)
  #in-silico PCR
  cds_pcrmark_3<- do.call('list',parallel::mclapply(cds_pcrmark_2$pair,function(x){insilicopcr(x)},mc.cores =thr))
  cds_pcrmark_4<-select_pcrmark(cds_pcrmark_2,cds_pcrmark_3)
  if(nrow(cds_pcrmark_4)==0){
    cds_rest_gene<-names(cds_gene_ls_chr_arm)
    cds_pcrmark_5<-cds_pcrmark_1[1]
    cds_pcrmark_5[[1]] <-cds_pcrmark_5[[1]][-c(1:length(cds_pcrmark_5[[1]]))]
  }else{cds_pcrmark_5<-select_pcrmark_from_grangelist(cds_pcrmark_4,cds_pcrmark_1)
  cds_rest_gene<-setdiff(names(cds_gene_ls_chr_arm),cds_pcrmark_4$gene)}
}
#Lowering standards for genes for which PCRmarks have not been designed
if(length(cds_rest_gene)==0){
  cds_rest_pcrmark_5<-cds_pcrmark_1[1]
  cds_rest_pcrmark_5[[1]] <-cds_rest_pcrmark_5[[1]][-c(1:length(cds_rest_pcrmark_5[[1]]))]}else{
    cds_rest_gene_ls<-include_select(cds_gene_ls_chr_arm,cds_rest_gene)
    cds_rest_cut<-make_windows_grl(cds_rest_gene_ls,24,3)
    cds_rest_pcrmark_1<-pcrmark_cds(cds_rest_gene_ls,cds_rest_cut,scores=10)
    if(max(lengths(cds_rest_pcrmark_1))==0){
      cds_rest_gene<-cds_rest_gene
      cds_rest_pcrmark_5<-cds_pcrmark_1[1]
      cds_rest_pcrmark_5[[1]] <-cds_rest_pcrmark_5[[1]][-c(1:length(cds_rest_pcrmark_5[[1]]))]
    }else{
      cds_rest_pcrmark_2<-makedf(cds_rest_pcrmark_1)
      cds_rest_pcrmark_3<- do.call('list',parallel::mclapply(cds_rest_pcrmark_2$pair,function(x){insilicopcr(x)},mc.cores = thr))
      cds_rest_pcrmark_4<-select_pcrmark(cds_rest_pcrmark_2,cds_rest_pcrmark_3)
      if(nrow(cds_rest_pcrmark_4)==0){
        cds_rest_gene<-cds_rest_gene
        cds_rest_pcrmark_5<-cds_pcrmark_1[1]
        cds_rest_pcrmark_5[[1]] <-cds_rest_pcrmark_5[[1]][-c(1:length(cds_rest_pcrmark_5[[1]]))]
      }else{cds_rest_pcrmark_5<-select_pcrmark_from_grangelist(cds_rest_pcrmark_4,cds_rest_pcrmark_1)
      cds_rest_gene<-setdiff(names(cds_rest_gene_ls),cds_rest_pcrmark_4$gene)}
    }
  }

if(length(cds_rest_gene)==0){
  cds_rest_pcrmark_5<-cds_rest_pcrmark_5}else{
    makedf_1<-function(x,num){
      len<-elementNROWS(x)
      total<-list()
      for(i in 1:length(x)){
        if(len[i]>num*2){
          select<-sample(1:(len[i]/2),size=num,replace = F)
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
    rest<-cds_rest_pcrmark_1[names(cds_rest_pcrmark_1)==cds_rest_gene]
    if(max(lengths(rest))==0){
      cds_rest_gene<-cds_rest_gene
      cds_rest_pcrmark_5<-cds_rest_pcrmark_5
    }else{
      rest_2<-makedf_1(rest,50)
      rest_3<- do.call('list',parallel::mclapply(rest_2$pair,function(x){insilicopcr(x)},mc.cores =thr ))
      rest_4<-select_pcrmark(rest_2,rest_3)
      if(nrow(rest_4)==0){
        cds_rest_gene<-cds_rest_gene
        cds_rest_pcrmark_5<-cds_rest_pcrmark_5
      }else{rest_5<-select_pcrmark_from_grangelist(rest_4,rest)
      cds_rest_pcrmark_5<-c(unlist(cds_rest_pcrmark_5),unlist(rest_5))
      cds_rest_pcrmark_5<-split(cds_rest_pcrmark_5,cds_rest_pcrmark_5$gene)
      cds_rest_gene<-setdiff(names(rest),rest_4$gene)}
    }
  }
# 2.6 Integrate all PCRmarks and replace sequences
PCRmarks<-c(unlist(in_in_pcrmark_5)[,c(1,2,3,9)],unlist(in_pcrmark_5)[,c(1,2,3,9)],unlist(cds_pcrmark_5)[,c(1,2,3,5)],unlist(cds_rest_pcrmark_5)[,c(1,2,3,5)])
PCRmarks<-split(PCRmarks,PCRmarks$gene)
manage_total_pcrmarker<-function(x){
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
PCRmarks<-manage_total_pcrmarker(PCRmarks)
manage_total_pcrmarker_2<-function(x){
  tm <- function(x){
    as.numeric(system(paste0(oligotm,' ',x),intern = T))
  }
  x<-unlist(x)
  x$wt_tm<-sapply(x$wt,tm)
  x$sf_tm<-sapply(x$sf,tm)
  #x<-split(x,x$gene)
  return(x)
}# add Tm
PCRmarks_2<-manage_total_pcrmarker_2(PCRmarks)
#output the information of PCRmarks
PCRmarks_df<- data.frame(
  chromosomee = as.character(seqnames(PCRmarks_2)),
  gene=PCRmarks_2$gene,
  start = start(PCRmarks_2),
  end = end(PCRmarks_2),
  wt_tag=as.character(PCRmarks_2$wt),
  syn_tag=as.character(PCRmarks_2$sf),
  strand=as.character(strand(PCRmarks_2)),
  wt_tm=PCRmarks_2$wt_tm,
  syn_tm=PCRmarks_2$sf_tm
)
PCRmarks_df2<- data.frame(
  Chr = as.character(seqnames(PCRmarks_2)),
  start = start(PCRmarks_2),
  end = end(PCRmarks_2),
  Edit_type="REPLACE",
  Seq=as.character(PCRmarks_2$sf),
  Feature_name=as.character(PCRmarks_2$gene),
  Feature_type="primer_bind",
  strand=as.character(strand(PCRmarks_2)),
  Ori_seq=as.character(PCRmarks_2$wt)
)
PCRmarks_df2$strand<-as.character(PCRmarks_df2$strand)
PCRmarks_df2$Feature_name<-as.vector(PCRmarks_df2$Feature_name)
for(i in 1:nrow(PCRmarks_df2)){
  if(PCRmarks_df2$strand[i]=='+'){
    PCRmarks_df2$strand[i] <- '1'
    PCRmarks_df2$Feature_name[i]=paste(PCRmarks_df2$Feature_name[i],'.syn.F',sep='')
  }
  else{
    PCRmarks_df2$strand[i] <- '-1'
    PCRmarks_df2$Feature_name[i]=paste(PCRmarks_df2$Feature_name[i],'.syn.R',sep='')}
}

write.table(PCRmarks_df2, paste(path,'Pp.Chr',chr_num,arm,"_PCRmark.txt",sep=''), sep = "\t", quote = FALSE, row.names = FALSE)

