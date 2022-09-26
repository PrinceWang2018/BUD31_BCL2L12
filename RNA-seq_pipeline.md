# Trim adapter

## Install

``` shell
source activate rna
conda install -c bioconda trim-galore
trim_galore --help
```

## Single end fastq files

``` shell
conda activate rna
cd /home/wzx/workdir
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
            /home/wzx/workdir/Samples.fastq.gz  \
            --gzip -o /home/wzx/outputdir/cleanfolder
```

Batch processing

``` shell
cd /home/wzx/project4t/PRJNA730205_Cell-Cell_EphB1_LYW_20220920/raw
ls *.fastq.gz | while read id ; do (   \
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
            ./$id  \
            --gzip -o ../clean \
); done
```

## Pair ends fastq files

``` shell
cd /home/wzx/workdir
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
            --paired /home/wzx/workdir/Sample1_1.fastq /home/wzx/workdir/Sample1_2.fastq  \
            --gzip -o /home/wzx/outputdir/cleanfolder
```

Batch processing

Need to prepare file_list like this:

```css
SRR7696207_1.fastq.gz   SRR7696207_2.fastq.gz
SRR8517853_1.fastq.gz   SRR8517853_2.fastq.gz
SRR8517854_1.fastq.gz   SRR8517854_2.fastq.gz
SRR8517855_1.fastq.gz   SRR8517855_2.fastq.gz
```

``` shell 
dir=/home/wzx/outputdir/cleanfolder/
cat file_list |while read id
do
      arr=${id}
      fq1=${arr[0]}
      fq2=${arr[1]}
      nohup trim_galore -q 20 --phred33 --length 20 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 &
done
```

# Files quality control

Fastqc

``` SHELL 
cd /home/wzx/outputdir/cleanfolder
ls /home/wzx/outputdir/cleanfolder/*gz | xargs fastqc -t 4 -o /home/wzx/outputdir/qc
```

Muitiqc

``` shell
multiqc /home/wzx/outputdir/cleanfolder
```

# Hisat2 align

## Index build

``` shell
hisat2-build /home/wzx/BIG_REF/miniconda3/envs/wes/bin/index/hisat2/hg38/hg38.fa /home/wzx/BIG_REF/miniconda3/envs/wes/bin/index/hisat2/hg38
```

## Single end

``` shell
#To sam
ls *.fq.gz | while read id ; do (hisat2 -p 4 --dta -x /home/wzx/BIG_REF/miniconda3/envs/wes/bin/index/hisat2/hg38/hg38 -U /home/wzx/outputdir/cleanfolder/$id -S /home/wzx/outputdir/alignfolder/$id.sam); done
#To bam
ls *.gz | while read id ; do (hisat2 -p 8 --dta -x /home/wzx/BIG_REF/miniconda3/envs/wes/bin/index/hisat2/hg38/hg38 -U /home/wzx/outputdir/cleanfolder/$id |  samtools sort -O bam -@ 8 -o - > /home/wzx/outputdir/alignfolder/$id.sort.bam); done
```

## Pair end

``` shell
hisat2 -p 8 --dta -x /home/wzx/BIG_REF/miniconda3/envs/wes/bin/index/hisat2/hg38/hg38 -1 /home/wzx/outputdir/cleanfolder/Sample1_1_val_1.fq.gz  -2 /home/wzx/outputdir/cleanfolder/Sample1_2_val_2.fq.gz  |  samtools sort -O bam -@ 8 -o - > /home/wzx/outputdir/alignfolder/Sample1.sort.bam
```

# Samtools

``` shell
#Samtools_sort
samtools sort -@ 4 -o /home/wzx/outputdir/alignfolder/Sample1.sort.bam /home/wzx/outputdir/cleanfolder/Sample1.sam

cd /home/wzx/outputdir/alignfolder/
ls *.fq.gz.sam | while read id ; do (samtools sort -@ 4 -o /home/wzx/outputdir/alignfolder/$id.bam  /home/wzx/outputdir/alignfolder/$id ); done


#Samtools_index
samtools index /home/wzx/outputdir/alignfolder/Sample1.bam -@ 8
ls *.bam | xargs -i samtools index {} -@ 8
```

# Bam files QC

``` SHELL
#bam qc
ls *.bam | while read id ;do (samtools flagstat $id > $(basename $id ".bam" ).stat ); done
```



# Deduplicated (Recommonded)

## Samtools_deduplicate

``` shell
#samtools_deduplicate
ls *.bam | while read id ;do (samtools markdup -r $id $(basename $id ".bam").rmdup.bam & );done
ls *.bam | while read id ;do (samtools fixmate -r $id $(basename $id ".bam").rmdup.bam -@ 6 & );done
```

## Picard

``` shell
#under java 1.8 (wes environment)
java -jar /home/wzx/biosoft/picard/picard.jar MarkDuplicates VERBOSITY=ERROR QUIET=true \
CREATE_INDEX=true REMOVE_SEQUENCING_DUPLICATES=true REMOVE_DUPLICATES=false \
INPUT=/home/wzx/outputdir/alignfolder/Sample1.bam  \
OUTPUT=/home/wzx/outputdir/dedup/Sample1_dedup.bam  \
M=/home/wzx/outputdir/dedup/Sample1_dedup.log  \
VALIDATION_STRINGENCY=SILENT MAX_OPTICAL_DUPLICATE_SET_SIZE=30000

# Batch processing
ls *bam  | while read id ;do (  
java -jar /home/wzx/biosoft/picard/picard.jar MarkDuplicates VERBOSITY=ERROR QUIET=true \
CREATE_INDEX=true REMOVE_DUPLICATES=true \
INPUT=/home/wzx/outputdir/alignfolder/$id  \
OUTPUT=/home/wzx/outputdir/dedup/"$id"_dedup.bam  \
M=/home/wzx/outputdir/dedup/"$id"_dedup.log
 );done
 
 # REMOVE_DUPLICATES=true and REMOVE_SEQUENCING_DUPLICATES=true have the same consequence
```

Note: Samtools_deduplicate will remove more sequence.



# FeatureCounts quantification

``` shell
#Single end sequencing
featureCounts -T 4 -t exon -g gene_id   \
-a /home/wzx/BIG_REF/reference_gencode_genome_rna_gtf/gencode.v33.annotation_hg38.gtf  \
-o  /home/wzx/outputdir/featurecounts/Samples_hg38_3v3.txt  \
/home/wzx/outputdir/alignfolder/*.bam

#-p pair ends sequencing
featureCounts -T 8 -t exon -g gene_id  -p \
-a /home/wzx/BIG_REF/reference_gencode_genome_rna_gtf/gencode.v33.annotation_hg38.gtf \
-o  /home/wzx/outputdir/featurecounts/Samples_hg38_3v3.txt  \
/home/wzx/outputdir/alignfolder/*.bam
```

# Differential expression in R

DESeq2 pipeline

``` R
##Dataframe preparation
rm(list = ls())
options(stringsAsFactors= F)
setwd('/home/wzx/outputdir/featurecounts/')
a=read.table('Samples_hg38_3v3.txt', head=T)
meta=a[,1:6]
exprSet=a[,7:ncol(a)]
colnames(exprSet)
colnames(exprSet)<-c("Treatment1","Treatment2","Treatment3","Control1","Control2","Control3")
gsub("\\.*$"," ",meta[,1])
rownames(exprSet)=meta[,1]
write.table(exprSet, file = "exprSet_counts.txt", sep = "\t")
group_list=c('Treatment','Treatment','Treatment','Control','Control','Control')  #OR colname of specific table

##Heatmap corelation
library(corrplot)
library(pheatmap)
pheatmap(scale(cor(log2(exprSet+1))))

## hcluster
# colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                cex = 0.7, col = "blue")
hc=hclust(dist(t(log2(exprSet+1))))
par(mar=c(5,5,5,10))
#png('hclust.png',res=120)
plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
dev.off()

## Analysis
##################################################################################
####################################DESEQ2########################################
##################################################################################
library("DESeq2")
suppressMessages(library(DESeq2))
(colData <- data.frame(row.names = colnames(exprSet),
                       group_list=group_list))
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)           #CALCULATE
dds <- DESeq(dds)

##diff_gene
res <-results(dds,
              contrast = c("group_list","Treatment","Control"))  ##TREAT / CONTROL  #READ RESULT
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG_result=as.data.frame(resOrdered)

#DEG_result = na.omit(DEG_result)  #TO DELETE NA ROW
nrDEG=DEG_result

##HEATMAP
library(pheatmap)
choose_gene=head(rownames(nrDEG),100)  ###50 maybe better
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix,filename='DEG_TOP100_HEATMAP.png')

logFC_cutoff <- with(DEG_result, mean(abs(log2FoldChange))+ 2*sd(abs(log2FoldChange)))
logFC_cutoff=1      # cutoff value should be set as required
DEG_result$change = as.factor(ifelse(DEG_result$pvalue < 0.05 &abs(DEG_result$log2FoldChange) > logFC_cutoff,
                              ifelse(DEG_result$log2FoldChange > logFC_cutoff, "UP","DOWN"),'NOT')
)
write.table(DEG_result,"DEG_RESULT_new",sep = "\t")
this_tile <- paste0('Cutoff for log FC is ', logFC_cutoff)

##vocalno plot
library(ggplot2)
g = ggplot(data=DEG_result,aes(x=log2FoldChange,y=-log10(pvalue),
                               color=change))+
  geom_point(alpha=0.4, size=1.75)+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change")+ylab("-log10 p-value")+
  ggtitle(this_tile)+theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) 
print(g)
ggsave(g,filename = 'volcano.png')

png("qc_dispersions.png",1000,1000,pointsize = 20)
plotDispEsts(dds,main="Dispersion plot")
dev.off()

##normaliation
rld <- rlogTransformation(dds)
exprMatrix_rlog=assay(rld)
exprMatrix_rlog_frame<- as.data.frame(exprMatrix_rlog)
write.table(exprMatrix_rlog_frame,'exprMatrix.rlog.txt', sep = "\t")

#check_normalization
png("DEseq_RAWvsNORM.png",height=800,width=800)
par(cex=0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex=0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(exprSet,col=cols, main="expression value",las=2)
boxplot(exprMatrix_rlog,col=cols,main="expression value",las=2)
hist(as.matrix(exprSet))
hist(exprMatrix_rlog)
dev.off()
```

# Annotation and FPKM transformation

``` R
###############counts_rpkm##################
mycounts<-read.table("/home/wzx/outputdir/featurecounts/Samples_hg38_3v3.txt", header = T)
n = grep(pattern = "PAR_Y",mycounts$Geneid)
mycounts <- mycounts[-n,]
name1 <- mycounts[,1]
name2 <- gsub("\\..*$","",name1)
rownames(mycounts) <- name2
mycounts <- mycounts[,-1:-5]
kb <- mycounts$Length /1000
countdata <- mycounts[,2:7]
rpk<-countdata / kb
tpm <- t(t(rpk)/colSums(rpk)*1000000)
fpkm <- t(t(rpk)/colSums(countdata)*1000000)
#ID_NAME
ENSG_NAME <- row.names(fpkm)
library("clusterProfiler")
library(org.Hs.eg.db)
gene.df <- bitr(ENSG_NAME, fromType = "ENSEMBL", #fromType indicates which type of data
                toType = c("SYMBOL"), #toType output data type (allow more than one type)
                OrgDb = org.Hs.eg.db)#Orgdb indicates which annotation packages
head(gene.df)
gene.df <- data.frame (gene.df)
mycounts <- cbind(rownames(mycounts),mycounts)
merge<- merge(mycounts,gene.df, by.y="ENSEMBL", by.x="rownames(mycounts)",all.x= TRUE)
fpkm <- data.frame (fpkm)
fpkm <- cbind(rownames(fpkm),fpkm)
merge<- merge(merge,fpkm, by.y="rownames(fpkm)", by.x="rownames(mycounts)",all.x= TRUE)
##diff
degseq2<- read.table("/home/wzx/outputdir/featurecounts/DEG_RESULT_new")
degseq2 <- cbind(rownames(degseq2),degseq2)
n = grep(pattern = "PAR_Y",degseq2$`rownames(degseq2)`)
degseq2 <- degseq2[-n,]
name1 <- degseq2[,1]
name2 <- gsub("\\..*$","",name1)
rownames(degseq2) <- name2
degseq2 <- degseq2[,-1]
degseq2 <- cbind(rownames(degseq2),degseq2)
merge<- merge(merge,degseq2, by.y="rownames(degseq2)", by.x="rownames(mycounts)",all.x= TRUE)
write.table(merge,"/home/wzx/outputdir/featurecounts/merge_fpkm_samples_date.txt", sep = "\t", row.names = FALSE)
```

# Acknowledgement

Thanks for the generous help and source code from Biotrainee group and github maintainers. Give respect to all the specialists in bioinformatics. If you have any questions, please contant wangzixiang@live.com without hesitaion.





