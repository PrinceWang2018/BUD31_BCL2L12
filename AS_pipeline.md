# rMATS

``` SHELL
conda activate sf
cd /home/wzx/workdir/rMATS
```

## From bam files

``` shell
rmats.py  \
--b1 /home/wzx/workdir/rMATS/b1.txt --b2 /home/wzx/workdir/rMATS/b2.txt \
--gtf /home/wzx/BIG_REF/reference_gencode_genome_rna_gtf/gencode.v33.annotation_hg38.gtf \
--od /home/wzx/workdir/rMATS/MATS_result \
--tmp /home/wzx/workdir/rMATS/tmp  \
-t paired \
--readLength 150 \
--cstat 0.0001 \
--nthread 14 \
--libType fr-unstranded

#b1.txt
/home/wzx/workdir/alignfile/Treatment_1_Clean.sam.bam,/home/wzx/workdir/alignfile/Treatment_2_Clean.sam.bam,/home/wzx/workdir/alignfile/Treatment_3_Clean.sam.bam

#b2.txt
/home/wzx/workdir/alignfile/NC_1_Clean.sam.bam,/home/wzx/workdir/alignfile/NC_2_Clean.sam.bam,/home/wzx/workdir/alignfile/NC_3_Clean.sam.bam
```



## From fastq

```shell
#from fastq
python rmats.py --s1 /path/to/s1.txt \          
                --s2 /path/to/s2.txt \             
                --gtf /path/to/the.gtf \           
                --bi /path/to/STAR_binary_index \   
                -t paired \                         
                --readLength 150 \                 
                --nthread 4 \                       
                --od /path/to/output \              
                --tmp /path/to/tmp_output           
```

# MISO 

> downloaded from: hollywood.mit.edu/burgelab/miso/annotations/gene-models/Mus_musculus.NCBIM37.65.gff.zip

Note: MISO analysis with hg19 annotation was strongly recommended.

## Build index

``` shell
#version1
index_gff --index ./hg19/A3SS.hg19.gff3 ./hg19/A3SS
index_gff --index ./hg19/A5SS.hg19.gff3 ./hg19/A5SS
index_gff --index ./hg19/MXE.hg19.gff3 ./hg19/MXE
index_gff --index ./hg19/RI.hg19.gff3 ./hg19/RI
index_gff --index ./hg19/SE.hg19.gff3 ./hg19/SE
index_gff --index ./hg19/TandemUTR.hg19.gff3 ./hg19/TandemUTR
index_gff --index ./hg19/AFE.hg19.gff3 ./hg19/AFE
index_gff --index ./hg19/ALE.hg19.gff3 ./hg19/ALE

#version2_hg19
index_gff --index ./hg19/A3SS.hg19.gff3 ./hg19/A3SS
index_gff --index ./hg19/A5SS.hg19.gff3 ./hg19/A5SS
index_gff --index ./hg19/MXE.hg19.gff3 ./hg19/MXE
index_gff --index ./hg19/RI.hg19.gff3 ./hg19/RI
index_gff --index ./hg19/SE.hg19.gff3 ./hg19/SE

#version2_hg38 (liftover from hg19)
> Note: This version just liftover the AS site but not the events name.
index_gff --index ./hg38/A3SS.hg38.gff3 ./hg38/A3SS
index_gff --index ./hg38/A5SS.hg38.gff3 ./hg38/A5SS
index_gff --index ./hg38/MXE.hg38.gff3 ./hg38/MXE
index_gff --index ./hg38/RI.hg38.gff3 ./hg38/RI
index_gff --index ./hg38/SE.hg38.gff3 ./hg38/SE

```

##  Calculate paired end insert distribution-bbmap

``` SHELL
#use bbmap
/home/wzx/miniconda3/envs/sf/opt/bbmap-38.79-0/bbmap.sh in1=/home/wzx/workdir/clean/Sample1_Clean_Data1.fq \
 in2=/home/wzx/workdir/clean/Sample1_Clean_Data2.fq \
 ref=/home/wzx/BIG_REF/miniconda3/envs/wes/bin/index/hg19/hg19.fa ihist=/home/wzx/reference_miso_AS/ihist.txt reads=2m pairlen=1000 threads=2 -Xmx28g
```

## Fragment Insertion Size Distribution-samtools

``` shell
samtools view -f Ox40 Repl .final .bam lperl -ane 'print abs($F[B]);print "\n'';' lsort -n >insertSize .txt

Rscript plotinsertSize.R insertSize .txt insertSize .png
```

plotinsertSize.R

```R
cmd=commandArgs(trailingOnly=TRUE);
input=cmd[l];
output=cmd[2];
d=read.table(input);
png(file=output);
hist(d$Vl ,main="Insertion Size distribution" ,ylab="Read Count" ,xlab="Insert Size" ,xaxt="n" ,breaks=seq(0,max (d$Vl ,by=l )));
axis(side=l ,at=seq (0,max (d$Vl) ,by=l00) ,labels=seq (0,max (d$Vl) ,by=l00 ));
dev.off()
```



## Run MISO

> BAM files created by hisat2 should include --dta parameter.

``` shell
miso --run /home/wzx/BIG_REF/reference_miso_AS/hg19/SE \
/home/wzx/outdir/align/Sample1_val_hg19.sam.bam \
--output-dir /home/wzx/outdir/MISO/BUD31_3_val_hg19_SE_re \
--read-len 141 \
--paired-end 240 117 \
--settings-filename /home/wzx/miniconda3/envs/sf/lib/python2.7/site-packages/misopy/settings/miso_settings.txt
```

Batch processing

``` shell
cd /home/wzx/outdir/align/
ls *.bam | while read id; do (
miso --run /home/wzx/BIG_REF/reference_miso_AS/hg19/SE \
/home/wzx/outdir/align/$id \
--output-dir /home/wzx/outdir/MISO/$id-miso-SE \
--read-len 141 \
--paired-end 240 117 \
--settings-filename /home/wzx/miniconda3/envs/sf/lib/python2.7/site-packages/misopy/settings/miso_settings.txt \
)done;
```

## Summary results (option)

``` shell
summarize_miso --summarize-samples /home/wzx/outdir/MISO/Sample1_hg19_SE_re /home/wzx/outdir/MISO/Sample1_hg19_SE_re 

cd /home/wzx/outdir/MISO/
ls -d * | while read id ;do (summarize_miso --summarize-samples /home/wzx/outdir/MISO/$id /home/wzx/outdir/MISO/$id/) done;
```

## Compare results

``` shell
compare_miso --compare-samples /home/wzx/outdir/MISO/Sample1_val_hg19_SE_re /home/wzx/outdir/MISO/Sample2_val_hg19_SE_re /home/wzx/outdir/MISO/comparsions_hg19
```

## Filter results

``` SHELL
filter_events \
--filter  /home/wzx/outdir/MISO/comparsions_hg19/NC_2_val_hg19_SE_re_vs_BUD31_3_val_hg19_SE_re/bayes-factors/NC_2_val_hg19_SE_re_vs_BUD31_3_val_hg19_SE_re.miso_bf \
--num-inc 1 \
--num-exc 1 \
--num-sum-inc-exc 10 \
--delta-psi 0.2 \
--bayes-factor 10 \
--output-dir /home/wzx/outdir/MISO/filter

#defalt bayes-factor is 10
```

## Gene annotation (Personal)

```	r
#MISO_MERGE_COMPARE_RESULT_ANNOATION
#MERGE_TABLES
setwd("/home/wzx/Project/BUD31_MATS_MISO_20200506/MISO_RESULT/filter")
table1<- read.delim("NC_1_Clean.sam.bam-miso_vs_IR_BUD31_1.miso_bf.filtered")
table2<- read.delim("NC_2_Clean.sam.bam-miso_vs_IR_BUD31_2.miso_bf.filtered")
MISO_Compare_merge_IR<- merge(table1,table2, by.x="event_name", by.y="event_name",all.x= FALSE)
write.table(MISO_Compare_merge_IR,file = "./MISO_Compare_merge_IR.txt", sep = "\t", row.names = FALSE)

#transfrom to bed
library(tidyr)
chr_site <- as.data.frame(MISO_Compare_merge_IR$event_name)
chr_site$`MISO_Compare_merge_IR$event_name` <- gsub("@.*$","",chr_site$`MISO_Compare_merge_IR$event_name`)
chr_site <- separate(chr_site,`MISO_Compare_merge_IR$event_name`,into = c("V1","V2","V4"),sep = ":")
chr_site <- separate(chr_site,`V2`,into = c("V2","V3"),sep = "-")
chr_site$V3 <- as.numeric(chr_site$V2)+1
write.table(chr_site,"./MISO_Compare_merge_IR_LOC.bed",sep = "\t", row.names = FALSE,col.names  = FALSE)

MISO_Compare_merge_IR$first_start <- chr_site$V3  #inorder to merge gene name
#annotation
library(devtools)
library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(clusterProfiler)
peak <- readPeakFile( "./MISO_Compare_merge_IR_LOC.bed" )

##remove_useless_chr
keepChr= !grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
##report peak number
cat(paste0('there are ',length(peak),' peaks for this data'))
#annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno = annotatePeak(peak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")  
peakAnno_df <- as.data.frame(anno)
geneinfor <- peakAnno_df[,16:18]
geneinfor$first_start <- peakAnno_df[,2]
MISO_Compare_merge_IR_GENE_Anno<- merge(geneinfor,MISO_Compare_merge_IR, by.x="first_start", by.y="first_start",all.x= TRUE)

write.table (MISO_Compare_merge_IR_GENE_Anno, "MISO_Compare_merge_IR_GENE_Anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```



## Sashimi plot

``` shell
sashimi_plot \
--plot-event "chr22:18165980:18166087:+@chr22:18171752:18171908:+@chr22:18178907:18178976:+@chr22:18185009:18185152:+" \
/home/wzx/BIG_REF/reference_miso_AS/hg19/SE/ \
/home/wzx/outdir/MISO/sashimi_plot_settings.txt  \
--output-dir /home/wzx/outdir/MISO/sashimi_plot
```

Example file of /home/wzx/outdir/MISO/sashimi_plot_settings.txt

```	css
[data]

# directory where BAM files are

bam_prefix = ./test-data/bam-data/

# directory where MISO output is

miso_prefix = ./test-data/miso-data/

bam_files = [
    "heartWT1.sorted.bam",
    "heartWT2.sorted.bam",
    "heartKOa.sorted.bam",
    "heartKOb.sorted.bam"]

miso_files = [
    "heartWT1",
    "heartWT2",
    "heartKOa",
    "heartKOb"]

[plotting]

# Dimensions of figure to be plotted (in inches)

fig_width = 7
fig_height = 5

# Factor to scale down introns and exons by

intron_scale = 30
exon_scale = 4

# Whether to use a log scale or not when plotting

logged = False
font_size = 6

# Max y-axis

ymax = 150

# Whether to plot posterior distributions inferred by MISO

show_posteriors = True

# Whether to show posterior distributions as bar summaries

bar_posteriors = False

# Whether to plot the number of reads in each junction

number_junctions = True

resolution = .5
posterior_bins = 40
gene_posterior_ratio = 5

# List of colors for read denisites of each sample

colors = [
    "#CC0011",
    "#CC0011",
    "#FF8800",
    "#FF8800"]

# Number of mapped reads in each sample

# (Used to normalize the read density for RPKM calculation)

coverages = [
    6830944,
    14039751,
    4449737,
    6720151]

# Bar color for Bayes factor distribution

# plots (--plot-bf-dist)

# Paint them blue

bar_color = "b"

# Bayes factors thresholds to use for --plot-bf-dist

bf_thresholds = [0, 1, 2, 5, 10, 20]
```

# Spliceseq

> https://bioinformatics.mdanderson.org/public-software/spliceseq/

Avaliable version for windows:

bowtie 1.0.0

jdk1.8.0_29164 (8u291 windows x64)

mysql 5.5.60-winx64

``` shell
mysql -u root -p
```

``` mysql
Welcome to the MySQL monitor. Commands end with ; or \g.
Server version: 5.0.51b-community-nt MySQL Community Edition (GPL)
mysql> CREATE DATABASE SpliceGraph CHARACTER SET latin1 COLLATE latin1_general_ci;
Query OK, 1 row affected (0.05 sec)

mysql>GRANT SELECT ON SpliceGraph.* TO 'sguser'@'localhost' IDENTIFIED BY 'sgpass';
mysql>GRANT SELECT ON SpliceGraph.* TO 'sguser'@'%' IDENTIFIED BY 'sgpass';

mysql>GRANT ALL ON SpliceGraph.* TO 'sgload'@'localhost' IDENTIFIED BY 'sg4ld!';
mysql>GRANT ALL ON SpliceGraph.* TO 'sgload'@'%' IDENTIFIED BY 'sg4ld!';
mysql>exit;
```

``` she	
$ unzip SpliceGraphDB.zip
$ mysql -u sgload -p SpliceGraph < SpliceGraphDB.sql
Enter password: ******
```

# SpliceR

## Hisat2 align and cufflink assemble

This classical R package required cufflink aligned files.

```shell
###cufflink_stringtie_gtf
#hisat2 bam prepare
hisat2 --fr -p 6 --dta-cufflinks -x /home/wzx/miniconda3/envs/wes/bin/index/hisat2/grch37_snp_tran/genome_snp_tran -1 /home/wzx/workdir/BUD31_2_Clean_Data1.fq.gz  -2 /home/wzx/workdir/BUD31_2_Clean_Data2.fq.gz  |  samtools sort -O bam -@ 8 -o - > /home/wzx/workdir/hisat2_hg19_tran_snp/hisat2_hg19_tranBUD31_2_Clean_hg19.bam

###cufflink assamble
#install under r3
conda activate r3
cufflinks  -p 12 -G /home/wzx/BIG_REF/reference_NCBI_GTF_REFERENCE/Homo_sapiens.GRCh37.75.gtf -o /home/wzx/outputdir/cufflinks/nc2  /home/wzx/workdir/hisat2_hg19_tran_snp/hisat2_hg19_tranNC_2_Clean_hg19.bam

###cuffmerge
cuffmerge -o ./merged_asm_new -p 12 -g /home/wzx/reference_NCBI_GTF_REFERENCE/Homo_sapiens.GRCh37.75.gtf ./cuffmerge_list.txt 
#-g /home/wzx/reference_gencode_genome_rna_gtf/gencode.v19.annotation.gtf -s /home/wzx/miniconda3/envs/wes/bin/index/hisat2/hg19/hg19.fa \
```

###./cuffmerge_list.txt 
./nc/transcripts.gtf
./bud/transcripts.gtf
./nc3/transcripts.gtf
./bud3/transcripts.gtf

## Differential expression with cuffdiff

```shell
###cuffdiff
ctrl_bam=/home/wzx/workdir/hisat2_hg19_tran_snp/hisat2_hg19_tranNC_2_Clean_hg19.bam,/home/wzx/workdir/hisat2_hg19_tran_snp/hisat2_hg19_tranNC_3_Clean_hg19.bam
treat_bam=/home/wzx/workdir/hisat2_hg19_tran_snp/hisat2_hg19_tranBUD31_2_Clean_hg19.bam,/home/wzx/workdir/hisat2_hg19_tran_snp/hisat2_hg19_tranBUD31_3_Clean_hg19.bam
label=hey_NC,hey_shBUD
cuffdiff -o ./diff_out1 -p 12 --labels $label --min-reps-for-js-test 2 /home/wzx/workdir/merged_asm_new/merged.gtf $ctrl_bam $treat_bam

#-min-reps-for-js-test how many repeats in treat group (>=2)
```

## Stringtie

#install under rna

```shell
conda activate rna
stringtie -p 8 -G /home/wzx/BIG_REF/reference_gencode_genome_rna_gtf/gencode.v19.annotation.gtf -o /home/wzx/workdir/NC_2_val_hg19.gtf -l NC_2_val_hg19 -A gene_abund_nc2.out /home/wzx/workdir/NC_2_val_hg19.sam.bam
```

## SingleR in R

```R
##conda envrionment r3 (python 2.7)
##R version 3.4.1 (<3.5)
##R code from vignette source 'vignettes/spliceR/inst/doc/spliceR.Rnw'

library("spliceR")
source("http://bioconductor.org/biocLite.R")
biocLite("spliceR")

###################################################
### code chunk number 1: Cufflinks_workflow
###################################################
library("spliceR")
setwd( "/home/wzx/workdir/cufflinks/diff_out1/")
dir <- "/home/wzx/workdir/cufflinks/diff_out1"

cuffDB <- readCufflinks(dir=dir,
                          gtf="./merged.gtf",genome="hg19" ,rebuild=TRUE)
cuffDB_spliceR <- prepareCuff(cuffDB)

###################################################
### code chunk number 2: Helper_functions
###################################################
myTranscripts <- transcripts(cuffDB_spliceR)
write.table(myTranscripts,file = "./myTranscripts.txt", sep = "\t")

myExons <- exons(cuffDB_spliceR)
conditions(cuffDB_spliceR)

###################################################
### code chunk number 8: preSpliceRFilter
###################################################
cuffDB_spliceR_filtered <- preSpliceRFilter(
  cuffDB_spliceR,
  filters=c("expressedIso", "isoOK", "expressedGenes", "geneOK")
)


###################################################
### code chunk number 9: spliceR
###################################################
#Commented out due to problems with vignettes and progress bars.
mySpliceRList <- spliceR(
  cuffDB_spliceR,
  compareTo="preTranscript",
  filters=c("expressedGenes","geneOK", "isoOK", "expressedIso", "isoClass"),
  useProgressBar=F
)

###################################################
### code chunk number 10: annotatePTC
###################################################
ucscCDS <- getCDS(selectedGenome="hg19", repoName="UCSC")

source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
#biocLite("BSgenome.Hsapiens.UCSC.hg38")

require("BSgenome.Hsapiens.UCSC.hg19",character.only = TRUE)
#Commented out due to problems with vignettes and progress bars.
PTCSpliceRList <- annotatePTC(cuffDB_spliceR, cds=ucscCDS, Hsapiens, 
	PTCDistance=50)

PTCSpliceRList2 <- spliceR(
  PTCSpliceRList,
  compareTo="preTranscript",
  filters=c("expressedGenes","geneOK", "isoOK", "expressedIso", "isoClass"),
  useProgressBar=F
)

save.image("./PTC_SpliceR_anno.RData")
load("./PTC_SpliceR_anno.RData")

topisoshift_ptc_result <- topIsoShift(PTCSpliceRList,n=1858180)   ##n = infinity
write.table(topisoshift_ptc_result, file = "./topisoshift_ptc_result_1858180.txt", sep = "\t")

###################################################
### code chunk number 11: generateGTF
###################################################
generateGTF(mySpliceRList2, filters=
              c("geneOK", "isoOK", "expressedGenes", "expressedIso"), 
            scoreMethod="local",
            useProgressBar=F)

generateGTF(PTCSpliceRList2)

###################################################
### code chunk number 12: plot1
###################################################
#Plot the average number of transcripts pr gene
mySpliceRList <- spliceRPlot(PTCSpliceRList2, 
                             evaluate="nr_transcript_pr_gene")


###################################################
### code chunk number 13: plot2
###################################################
#Plot the average number of Exon skipping/inclucion evens pr gene
mySpliceRList <- spliceRPlot(PTCSpliceRList2, evaluate="nr_AS", 
                             asType="All")
###evaluate :   nr_transcript  nr_gene   nr_AS  nr_transcript_pr_gene  
#Upon inital usage of spliceRPlot, theSpliceRListis initiated with internal data, allowing for fasterreplotting.  If theSpliceRListchanges because of filtering or other manipulation, rerun spliceR-Plot withreset=T. For the evaulate parameter, the following are valid:  ’nr_transcript’,’nr_genes’,’nr_transcript_pr_gene’, ’nr_AS’, ’mean_AS_gene’, ’mean_AS_transcript’, ’mean_transcript_exp’,’mean_gene_exp’.  ’nr_transcript’ outputs number of transcripts, ’nr_AS’ outputs number of alter-native splicing events, ’mean_as’ outputs the average number of AS events per gene, ’mean_transcript_exp’outputs the mean transcript expression and ’mean_gene_exp’ output the mean gene expression. Fora detailed description of filters, seesplice
###asType:  ESI MEE MESI ISI A5 A3 ATSS ATTS All

###################################################
### code chunk number 14: plot3
###################################################
#Plot the average number of Exon skipping/inclucion evens pr gene, 
#but only using transcripts that are significantly differntially expressed
mySpliceRList <- spliceRPlot(PTCSpliceRList2, evaluate="mean_AS_transcript", 
                             asType="ESI", filters="sigIso",reset=TRUE)

####topisoshift########
topisoshift_result <- topIsoShift(mySpliceRList,n=20)
write.table(topisoshift_result, file = "./topisoshift_result.txt", sep = "\t")

#####topisoshift_annotatePTC##################
topisoshift_ptc_result <- topIsoShift(PTCSpliceRList2,n=20)
write.table(topisoshift_ptc_result, file = "./topisoshift_ptc_result.txt", sep = "\t")

#####totalNumberofAS##############################
mySpliceRList <- spliceRPlot(mySpliceRList, evaluate= "nr_AS", asType="ATTS")

TotalNumberOfAS_result <- totalNumberOfAS(mySpliceRList)
write.table(TotalNumberOfAS_result, file = "./TotalNumberOfAS_result.txt", sep = "\t")

final_all_spliceR_output <- GenomicRanges::as.data.frame(PTCSpliceRList2[[1]] )    ###by prince wang 20201013
write.table(final_all_spliceR_output, file = "./try.txt", sep = "\t")
########Transcripts/Exon_information####################
library("GenomicRanges")
#these two are GRanges Objects
myTranscripts <- transcripts(cuffDB_spliceR)
myExons <- exons(cuffDB_spliceR)
```

# Acknowledgement

Thanks for the generous help and source code from Biotrainee group and github maintainers. Give respect to all the specialists in bioinformatics. If you have any questions, please contant wangzixiang@live.com without hesitaion.
