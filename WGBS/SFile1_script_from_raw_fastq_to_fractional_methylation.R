######################################################################################
#0-run FastQC on raw fastq files for diagnostic purposes
######################################################################################
module load fastqc
fastqc file.fastq.gz

######################################################################################
#1- quality and adapter trimming. Versions used:
#Trim Galore version: 0.4.1
#Cutadapt version: 1.8.1
#Trimming mode: paired-end
#Quality Phred score cutoff: 20
#Quality encoding type selected: ASCII+33
#Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
#Maximum trimming error rate: 0.1 (default)
#Minimum required adapter overlap (stringency): 1 bp
#Minimum required sequence length for both reads before a sequence pair gets removed: 20 bp
#All sequences will be trimmed by 1 bp on their 3' end to avoid problems with invalid paired-end alignments with Bowtie 1
#Running FastQC on the data once trimming has completed
######################################################################################
module load python/2.7
module load cutadapt
module load fastqc
trim_galore --paired --trim1 R1.fastq.gz R2.fastq.gz --fastqc

######################################################################################
#2-mapping with Bismarkv0.14.5 and bowtie-1.1.2
#genome GRCh37 version from  ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq/
######################################################################################
module load samtools
bismark_v0.14.5/bismark --bowtie1 --path_to_bowtie bowtie-1.1.2 --multicore 7 -n 1 /path_to_genome/GRCh37/ -1 fastq.1.gz -2 fastq.2.gz -o /outdir/ --gzip --bam

######################################################################################
#3-deduplication step
######################################################################################
bismark_v0.14.5/deduplicate_bismark -p file.bam --bam

######################################################################################
#4-methylation extraction step
######################################################################################
module load samtools
bismark_v0.14.5/bismark_methylation_extractor --multicore 5 -p --no_overlap --comprehensive deduplicated.bam --gzip -o /outdir/ --buffer_size 50%  --ample_memory --report --bedGraph --cytosine_report --genome_folder /path_to_genome/GRCh37/  

######################################################################################
#read each individual's CpG methylation results from Bismark cytosine report to R using bsseq R library
#https://bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq.html
######################################################################################
library(bsseq)
bs <- read.bismark("sample1_CpG.txt",sep=''),strandCollapse=T,rmZeroCov=T,fileType="cytosineReport");
sampleNames(bs) <- "sample1"
#combine all samples (all controls for cell-type analyses or sz-control per cell-type for disease analyses) to run DSS later
cum.bs <- combine(sample1.bs,sample2.bs);
saveRDS(cum.bs,paste("bsseq.Rds",sep=''))
