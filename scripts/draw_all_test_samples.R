#!/usr/bin/Rscript
var_log<-snakemake@log[['log']]
var_ref<-snakemake@params[['ref']]
var_samples<-snakemake@params[['samples']]
var_maximum_alt_percent_in_control<-as.numeric(snakemake@params[['maximum_alt_percent_in_control']])
var_minimum_editing_site_coverage_in_test<-as.integer(snakemake@params[['minimum_editing_site_coverage_in_test']])
var_sequencelogo_span_len<-as.integer(snakemake@params[['sequencelogo_span_len']])


#var_log<-"log/draw_all_test_sample.stdout.txt"
#var_ref<-"/home/fengchr/staff/RNA_edit_pipeline/test_E_colli/reference/GCF_000005845.2_ASM584v2_genomic.fna"
#var_samples<-"P8E_1 P8E_2 P1723_1 P1723_2"
#var_maximum_alt_percent_in_control<-as.numeric(0.01)
#var_minimum_editing_site_coverage_in_test<-as.integer(10)
#var_sequencelogo_span_len<-as.integer(5)


log <- file(var_log, open="wt")
sink(log)

library(ggplot2)
library(Biostrings)
library(ggseqlogo)
library(dplyr)
library(qqman)
library(viridis)
fa<-readDNAStringSet(var_ref)
strip<-function(x){return(unlist(lapply(x, function(x){return(unlist(strsplit(x,split=" "))[1])})))}
names(fa)<-strip(names(fa))

sample<-unlist(strsplit(var_samples,split=" "))
res<-data.frame(sample=character(),fraction=numeric(),p_value=numeric())
logosets<-list()

print("log:")
print(snakemake@log[['log']])
print("ref:")
print(snakemake@params[['ref']])
print("samples:")
print(snakemake@params[['samples']])
print("maximum_alt_percent_in_control:")
print(snakemake@params[['maximum_alt_percent_in_control']])
print("minimum_editing_site_coverage_in_test:")
print(snakemake@params[['minimum_editing_site_coverage_in_test']])
print("sequencelogo_span_len:")
print(snakemake@params[['sequencelogo_span_len']])

for (base in sample) {
  data<-read.csv(paste0(base,".res"),sep = "\t",header = F,quote = '"',stringsAsFactors = F)
  num_controls<-(ncol(data)-7)/2
  data0<-data[,1:7]
  names(data0)<-c("chrom","pos","ref","alt","filter","test_ref_cov","test_alt_cov")
  data1<-data[8:ncol(data)]
  cons_ref_cov<-rowSums(data1[,seq(from=1,by=2,length.out=num_controls)])
  cons_alt_cov<-rowSums(data1[,seq(from=2,by=2,length.out=num_controls)])
  data<-cbind(data0,cons_ref_cov=cons_ref_cov,cons_alt_cov=cons_alt_cov)
  pass<-data[data$filter=="PASS",]
  bed<-read.csv(paste0(base,".annotation"),header = F,sep = "\t",stringsAsFactors = F,quote = '"')
  names(bed)<-c("chrom","pos","start","end","gene_id","gene_name","transcript_id","transcript_name","gene_type","transcript_type","strand")
  bed[is.na(bed)]<-""
  mbed<-bed %>% group_by(chrom,pos) %>% summarise(n_plus=sum(strand=="+"),n_minus=sum(strand=="-"))
  pass1<-merge(pass,mbed,by=c("chrom","pos"),all.x=T)
  tmp<-pass1[pass1$alt %in% c("G","C"),]
  tmp<-na.omit(tmp)
  if(nrow(tmp)==0){next}
  index<-tmp$n_plus-tmp$n_minus
  strand<-unlist(lapply(index, function(x){if(x>0){return("F")}else if(x<0){return("R")}else{return(NA)}}))
  tmp<-cbind(tmp,strand)
  tmp<-na.omit(tmp)
  used<-tmp[(tmp$ref=="A" & tmp$alt=="G" & tmp$strand=="F")|(tmp$ref=="T" & tmp$alt=="C" & tmp$strand=="R"),]
  #去除对照组alt%>maximum_alt_percent_in_control | 测试组覆盖度<minimum_editing_site_coverage_in_test的突变
  used<-used[used$cons_alt_cov/(used$cons_ref_cov+used$cons_alt_cov) < var_maximum_alt_percent_in_control & (used$test_ref_cov+used$test_alt_cov) >= var_minimum_editing_site_coverage_in_test,]
  #计算显著性
  fisher_ind<-unlist(apply(used[,6:9], 1, function(x){return(fisher.test(matrix(c(x[1],x[2],x[3],x[4]),nrow=2,byrow = F),alternative = "less")$p.value)}))
  used<-cbind(used,filsher_test=fisher_ind)
  #输出最终编辑位点-基因表
  gene<-merge(used,bed,by=c("chrom","pos"),all.x=T)
  iii<-1:length(names(fa))
  names(iii)<-names(fa)
  gene<-gene[order(iii[gene$chrom]),]
  write.table(gene,paste0(base,"_gene.tsv"),sep = "\t",row.names = F,quote = F)
  #输出bed文件
  outbed0<-data.frame(chrom=gene$chrom,start=gene$pos-1,end=gene$pos,name=gene$gene_id,strand=gene$strand.x)
  outbed0$strand[outbed0$strand=="F"]<-"+"
  outbed0$strand[outbed0$strand=="R"]<-"-"
  outbed<-outbed0 %>% group_by(chrom,start,end) %>% summarise(name=paste(name,collapse=","),score=100,strand=strand[1])
  write.table(outbed,paste0(base,".RNA_edits.bed"),sep = "\t",row.names = F,col.names = F, quote = F)
  
  vec<-used$test_alt_cov/(used$test_alt_cov+used$test_ref_cov)
  res<-rbind(res,data.frame(sample=rep(base,length(vec)),fraction=vec,p_value=fisher_ind))
  
  #sequence logo
  span<-var_sequencelogo_span_len
  faa<-DNAStringSet()
  logoset<-DNAStringSet(apply(used[,c("chrom","pos","strand")], 1, function(x){
    seq<-fa[[x[1]]][(as.integer(x[2])-span):(as.integer(x[2])+span)]
    if (x[3]=='F') {
      return(seq)
    }else{return(reverseComplement(seq))}
  }))
  logosets[[base]]<-as.character(logoset)
  
  #draw chrom plot
  mandata<-data.frame(chrom=used$chrom,pos=used$pos,freq=used$test_alt_cov/(used$test_alt_cov+used$test_ref_cov))
  mandata$chrom<-unlist(lapply(mandata$chrom,function(x){return(which(names(fa)==x))}))
  graphics.off()
  pdf(paste0(base,"_manhattan.pdf"),width = 7,height = 5.5)
  print(manhattan(mandata,chr = "chrom",bp = "pos",p = "freq",snp = "pos",chrlabs = names(fa)[sort(unique(mandata$chrom))],logp = F,col = rainbow(length(unique(mandata$chrom))),ylab="RNA A-to-G editing(%)"))
  dev.off()
}

ind<-1:length(sample)
names(ind)<-sample

graphics.off()
pdf("Number of RNA A-G edits.pdf",width = 6,height = 4)
numdata<-as.data.frame(table(res$sample))
numdata$Var1<-factor(numdata$Var1,levels = unique(numdata$Var1)[order(ind[unique(numdata$Var1)])])
ggplot(data=numdata,aes(x=Var1,y=Freq))+geom_bar(fill="steel blue",stat = "identity")+theme_bw()+labs(x="Sample",y="Number of A-G edits")+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
dev.off()

res$sample<-factor(res$sample,levels = unique(res$sample)[order(ind[unique(res$sample)])])
graphics.off()
pdf("A-G edit distribution2.pdf",width = 6,height = 4)
ggplot(data=res,aes(x=sample,y=fraction,color=-log10(p_value)))+geom_jitter(shape=16,size=1)+scale_color_viridis(name="-log(p)")+theme_bw()+labs(x="Sample",y="Fraction")+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
dev.off()

graphics.off()
pdf("A-G edit distribution.pdf",width = 6,height = 4)
ggplot(data=res,aes(x=sample,y=fraction))+geom_jitter(color="steel blue",shape=16,size=1)+theme_bw()+labs(x="Sample",y="Fraction")+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
dev.off()

graphics.off()
pdf("Sequence logo of A-G edits.pdf",width = (length(logosets)+0.2)*2.4,height = 2)
ggseqlogo(logosets, nrow=1)+theme(axis.text.x = element_text(size = 7))
dev.off()

write("success","success")