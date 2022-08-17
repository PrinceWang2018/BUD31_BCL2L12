# Trim adapter

## Install

``` shell
source activate rna
conda install -c bioconda trim-galore
trim_galore --help
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

# Remove rRNA sequence

## Build bowtie2 rRNA index

```shell
#BOWTIE2_INDEX
bowtie2-build /home/wzx/ncRNA_fasta/RNA_AND_TAXONOMY9606_AND_so_rna_type_nameRRNA.fasta /home/wzx/ncRNA_fasta/index/Cen_rRNA
```

## Remove rRNA

``` SHELL
1.single-end
bowtie2 -x ~/RNASEQ/index/rRNA_index/hg38_rRNA --un-gz ${i}_IP_rmrRNA.fastq.gz -U ${i}_IP.read1_Clean.fastq.gz -p 8 -S ${i}_rRNA.sam; rm ${i}_rRNA.sam
2.paired-end
bowtie2 -x ~/RNASEQ/index/rRNA_index/hg38_rRNA --un-conc-gz ${i}_input_rmrRNA.fastq.gz  -1 ${i}_in.read1_Clean.fastq.gz -2 ${i}_in.read2_Clean.fastq.gz -p 8 -S ${i}_rRNA.sam; rm ${i}_rRNA.sam; 

Note: This will write 2 files ending with .1 and .2, which were the fastq files removed rRNA.
```

Examples:

``` shell
bowtie2 -x /home/wzx/ncRNA_fasta/index/Cen_rRNA --un-conc-gz /home/wzx/project4t/CLIP_BUD_DATAREDUCTION_20201201/Input5_2_rmrRNA.fastq.gz -1 /home/wzx/project4t/bud_original_clip_20201009/X101SC20091110-Z01-J001/CleanFq/Input5_2_clean_1.fq.gz -2 /home/wzx/project4t/bud_original_clip_20201009/X101SC20091110-Z01-J001/CleanFq/Input5_2_clean_2.fq.gz -p 6 -S /home/wzx/project4t/CLIP_BUD_DATAREDUCTION_20201201/Input5_2_rRNA.sam; rm /home/wzx/project4t/CLIP_BUD_DATAREDUCTION_20201201/Input5_2_rRNA.sam; 

bowtie2 -x /home/wzx/ncRNA_fasta/index/Cen_rRNA --un-conc-gz /home/wzx/project4t/CLIP_BUD_DATAREDUCTION_20201201/Spy5_rmrRNA.fastq.gz -1 /home/wzx/project4t/bud_original_clip_20201009/X101SC20091110-Z01-J001/CleanFq/Spy5_clean_1.fq.gz -2 /home/wzx/project4t/bud_original_clip_20201009/X101SC20091110-Z01-J001/CleanFq/Spy5_clean_2.fq.gz -p 8 -S /home/wzx/project4t/CLIP_BUD_DATAREDUCTION_20201201/Spy5_rRNA.sam; rm /home/wzx/project4t/CLIP_BUD_DATAREDUCTION_20201201/Spy5_rRNA.sam; 
```

# Alignment with STAR

``` SHELL
#Index files were downloaded from Ensmble website.

STAR \
    --genomeDir /home/wzx/index_STAR_hg38/ \
     --runThreadN 8 \
    --readFilesIn /home/wzx/workdir/CleanFq/Spy5_clean_1.fq.gz /home/wzx/workdir/CleanFq/Spy5_clean_2.fq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix /home/wzx/workdir/bam_star_hg19/Spy5_star_hg38 \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 4  \
    --outFilterMultimapScoreRange 1 \
    --outFilterScoreMin 10  \
    --outFilterMultimapNmax 1
```

# Remove duplicates

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

# Crosslink sites identification-PURECLIP

``` SHELL
pureclip -i /home/wzx/workdir/Spy5_star_hg19Aligned.sortedByCoord.out.bam  \
      -bai /home/wzx/workdir/Spy5_star_hg19Aligned.sortedByCoord.out.bam.bai  \
      -g /home/wzx/index_STAR_hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa  \
      -ld -nt 16 -dm 80  \
      -o ./spy_bud_single.crosslink_site.bed  \
      -or ./spy_bud_single.crosslink_region.bed
```

# Region annotation and find motif-HOMER

## Install HOMER

```SHELL
Option1:
perl /home/wzx/miniconda3/envs/wes/share/homer-4.10-0/.//configureHomer.pl -install hg38
ls -lh  /home/wzx/miniconda3/envs/wes/share/homer-4.10-0/data

Option2:
cd ~/biosoft
mkdir homer &&  cd homer
wget http://homer.salk.edu/homer/configureHomer.pl 
perl configureHomer.pl -install
perl configureHomer.pl -install hg38
```

## Find motif

``` shell
##homer_find_motif
cd  /home/wzx/workdir/motif
for id in /home/wzx/workdir/motif/*.narrowPeak;
do
echo $id
file=$(basename $id )
sample=${file%%.*} 
echo $sample  
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $id >homer_peaks.tmp  
findMotifsGenome.pl homer_peaks.tmp hg38 ${sample}_motifDir -len 6,8,10,12 -S 10 -rna -p 4
annotatePeaks.pl  homer_peaks.tmp hg38  1>${sample}.peakAnn.txt 2>${sample}.annLog.txt 
done 

#####note!!!
#STAR created .bed must add "chr" in line1
```

# Find sequence belongs to motif

```shell
Example:
scanMotifGenomeWide.pl <motif file> <genome> [options]
scanMotifGenomeWide.pl .motif hg38 -bed > site.hg38.bed

conda activate wes
cd  /home/wzx//home/wzx/workdir/motif_find_seq_bed/
for id in /home/wzx//home/wzx/workdir/motif_find_seq_bed/*.motif;
do
echo $id
file=$(basename $id )
sample=${file%%.*} 
echo $sample 
scanMotifGenomeWide.pl $id hg38 -bed -p 8 > ${sample}.site.hg38.bed
done 
```

# Deeptools visualization

``` shell
bamCoverage -b spy_bud_single.crosslink_add80_region_up3_sort.bam -o spy_bud_single.crosslink_add80_region_up3_sort.bw
computeMatrix reference-point \
--referencePoint center \
-b 200 -a 200 -p 8  \
-R ./spy_bud_single.crosslink_80_site.bed \
-S ./spy_bud_single.crosslink_add80_region_up3_sort.bw \
--skipZeros \
-o ./spy_bud_single.crosslink_add80_region_up3_sort_-200-200.gz   

plotHeatmap -m spy_bud_input_single.crosslink_region_up3_newinput-200-200.gz   -out spy_bud_input_single.crosslink_region_up3_newinput-200-200.pdf --plotFileFormat pdf --dpi 720 --colorMap Blues
```



# Upset plot

``` Upset plot
###MOTIF1234
setwd("/home/wzx/workdir/MOTIF_UPSET/")
MOTIF1<-read.table("./crosslink_regionup3_motif1.bed",header = FALSE)
MOTIF2<-read.table("./crosslink_regionup3_motif2.bed",header = FALSE)
MOTIF3<-read.table("./crosslink_regionup3_motif3.bed",header = FALSE)
MOTIF4<-read.table("./crosslink_regionup3_motif4.bed",header = FALSE)
MOTIF12<- merge(MOTIF1,MOTIF2, by.x = "V4",  by.y = "V4",all = TRUE)
MOTIF123<- merge(MOTIF12,MOTIF3, by.x = "V4",  by.y = "V4",all = TRUE)
MOTIF1234<- merge(MOTIF123,MOTIF4, by.x = "V4",  by.y = "V4",all = TRUE)

upset_data <- MOTIF1234[,c(1,7,13,19,25)]
colnames(upset_data)<- c("clusterid","motif1","motif2","motif3","motif4")

df=as.data.frame(upset_data)
df[,-1]=(df[,-1]!=0)*1

upset(df, order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE),mb.ratio = c(0.75, 0.25))
      , nsets = 4, nintersects = 15, mb.ratio = c(0.7, 0.3),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))

```

# Correlation plot

``` R
##Corelation_dot_cluster_corelation
setwd("/home/wzx/workdir")
raw_counts_data <- read.table("./input_spy_22_corelation.count",header = FALSE)
#calculate_corelation
cor_num_count<-raw_counts_data[,7:10]
colnames(cor_num_count)<- c("Input_1","Input_2","SpyCLIP_1","SpyCLIP_2")
##log10&deal with infinity
cor_num_count_log10<- log10(cor_num_count)
cor_num_count_log10[cor_num_count_log10=="-Inf"]<- 0

cor(cor_num_count_log10)
install.packages("corrgram")
library(corrgram)
corrgram(cor_num_count_log10,order = TRUE, lower.panel = panel.pts,upper.panel = panel.conf, text.panel = panel.txt, diag.panel = panel.density)

#another task
cor_num_count<-cor_num_count[,-2]
cor_num_count<-as.data.frame(cor_num_count)
cor_num_count[cor_num_count=="-Inf"]<- 0

###cor
cor(as.data.frame(cor_num_count))
cor.test(cor_num_count$log2score,cor_num_count$log2spyc)
cor.test(cor_num_count$log2score,cor_num_count$log2inputc)
###log2
raw_counts_data$log2score <- log2(raw_counts_data$V5)
raw_counts_data$log2inputc <- log2(raw_counts_data$V7)
raw_counts_data$log2spyc <- log2(raw_counts_data$V8)
####smoothScatter
with(raw_counts_data,
     smoothScatter(log2inputc,log2spyc,nbin = 1024))
  abline(abline(a = 0, b = 1), col = 'gray')

with(raw_counts_data,
     smoothScatter(log2score,log2inputc,nbin = 1024))
####hexin_plot
install.packages("hexbin")
library("hexbin")
with(raw_counts_data,{
  bin <- hexbin(log2inputc,log2spyc,xbins = 50)
  plot(bin)
})    
```



# Acknowledgement

Thanks for the generous help and source code from github maintainers. Give respect to all the specialists in bioinformatics. If you have any questions, please contant wangzixiang@live.com without hesitaion.