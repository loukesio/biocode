# samtools 

## Installing samtools 

```r
cd ~
# optional. you may already have a src directory
mkdir src
cd ~/src
git clone https://github.com/samtools/htslib
git clone https://github.com/samtools/samtools
cd samtools
make
cp samtools ~/bin
```

now put in your folder the a sam file. My sam file is_ **barcodes_alignment.sam**_


## 1. Check the header of your sam file 

```
head barcodes_alignment.sam

>
@SQ	SN:Reference_barcodes	LN:172
@PG	ID:minimap2	PN:minimap2	VN:2.15-r905	CL:minimap2 -a Reference_barcodes.fasta subset_0-1-batch_S1_L00
```
your header should start with the sign "@",  which is an indicator for a header line. If you don't see lines starting with the "@" sign, the header information is most likely missing.

If the header is absent from the SAM file use the command below, where reference.fa is the reference fasta file used to map the reads:

samtools view -bT Reference_barcodes.fasta barcodes_alignment.sam > barcodes_alignment.bam

## 2. Convert your sam file to bam file 

```
samtools view -S -b barcodes_alignment.sam > barcodes_alignment.bam
```
S indicates that give it a sam file 
b indicates that you want bam out of it

## 3. Sort the bam file 

The sequences have the same order as in the fasta file that you used to map the reference. Always sort your BAM files; many downstream programs only take sorted BAM files.

```
samtools sort barcodes_alignment.bam -o barcodes_alignment_sorted.bam
```
## 4. Index the sorted bam file 
A BAM index file is usually needed when visualising a BAM file.

```
samtools index barcodes_alignment_sorted.bam > index_barcodes_alignment_sorted.bam
```
## 5. Count the mapped reads 

* count **all reads** in the bam file 

```
samtools idxstats index_barcodes_alignment_sorted.bam | awk '{s+=$3+$4} END {print s}'
```

* count only the **mapped** reads in the bam file
```
samtools idxstats index_barcodes_alignment_sorted.bam | awk '{s+=$3} END {print s}'

samtools view -c -F 4 index_barcodes_alignment_sorted.bam
```

* count only the **unmapped** reads in the bam file 

samtools view -c -f 4 index_barcodes_alignment_sorted.bam

## 6. Plot a depth of coverage across a bam file 

[1] https://labs.epi2me.io/notebooks/Introduction_to_SAM_and_BAM_files.html

Other important references: 
[1] http://quinlanlab.org/tutorials/samtools/samtools.html  : Notes from Quinlan lab <br>
[2] https://davetang.org/wiki/tiki-index.php?page=SAMTools  : Notes from wiki cms <br>
[3] http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/1040-2/ : extract unmapped reads from novo align <br>
[4] https://www.metagenomics.wiki/tools/samtools/number-of-reads-in-bam-file : count unmapped reads, metagenomics notes <br>
[5] http://www.bioinf.uni-leipzig.de/publications/supplements/13-008 : reuse unmapped reads with segemehl <br>
[6] https://davetang.github.io/learning_bam_file/#filtering-unmapped-reads : Notes from dave tang <br>
[7] https://qnot.org/2012/04/14/counting-the-number-of-reads-in-a-bam-file/ : Count number of reads thats very interesting


# Minimap 

## 1. Map to a masked reference 
Minimap has worked perfectly for me to map millions of sequences into a masked referece which looks like this 
```
> Reference_barcodes.fasta
GCTTGGGCCGATGTCCACGAAGCTCTCCTACGNNNNNNNNNNNNNNNNNNNNNNNNNCAGTCCAGCGCCAACCAGATA
```
the nice think of minimap is that you do not have to create an index of the reference. 
* Map to one sample

```
minimap2 -a Reference_barcodes.fasta subset_0-1-batch_S1_L001_R1_001.fasta > alignment.sam 
```
* Map to multiple samples using a for loop (snakemake might a better approach)

https://github.com/lh3/minimap2/blob/master/cookbook.md#map-sr

# Screenshots MAC
A script to save screenshots in a directory other than Desktop

# Weather MAC
A script to see weather in your terminal

#  Open Terminal type and press enter
defaults write com.apple.screencapture type jpg  <Drag the folder you want to save to into the Terminal command box> 

# Note!
you can save your screenshots in a different format than jpg e.g., pdf 
defaults write com.apple.screencapture type pdf




# lectures 

https://stat447.com/lectures/01-shell/ # these lectures are great

# everyday_commands 
material for the course Bash for Beginners - February 2022

```r
# make file  and move files inside them 
mkdir 1_assembled
mv *.assembled.fastq 1_assembled

mkdir 1_discarded
mv *.discarded.fastq 1_discarded

mkdir 1_unassembled
mv *.unassembled.* 1_unassembled
```
\mkdir 7_assembled
mv *.assembled.fastq 7_assembled

mkdir 7_discarded
mv *.discarded.fastq 7_discarded

mkdir 7_unassembled
mv *.unassembled.* 7_unassembled


/home/theodosiou/PROJECTS/BARCODES_PROJECT/BARCODES_DATA/Networks/181121_MiSeq/merged/0_assembled/extract_barcodes/cluster_barcodes/

/home/theodosiou/PROJECTS/BARCODES_PROJECT/BARCODES_DATA/Networks/combine_cluster017

# if you want to copy-paste the content then you use the command cp -R

 cp -R /home/theodosiou/PROJECTS/BARCODES_PROJECT/BARCODES_DATA/Networks/251121_MiSeq/merged/7_assembled/extract_barcodes/cluster_barcodes/* /home/theodosiou/PROJECTS/BARCODES_PROJECT/BARCODES_DATA/Networks/combine_cluster017

 #To determine how many files there are in the current directory, put in 

 ls -1 | wc -l
 
 # How to make fastq to fasta 
 ```
 sed -n '1~4s/^@/>/p;2~4p' INFILE.fastq > OUTFILE.fasta
 awk 'NR==1, NR==20' 0-1-batch_S1_L001_R1_001.fasta > subset_0-1-batch_S1_L001_R1_001.fasta
 blastn -query subset_0-1-batch_S1_L001_R1_001.fasta -subject Reference_barcodes.fasta -out test_blast.txt
```
 
 https://bioinformaticsworkbook.org/dataWrangling/fastaq-manipulations/converting-fastq-format-to-fasta.html#gsc.tab=0
 
 # Blast in your MAC 
 
 https://www.youtube.com/watch?v=sU44ZtZlzo4&list=LL&index=24&t=90s
 
 ### Get coverage with SAMtools 
 https://www.metagenomics.wiki/tools/samtools
 
 
 ### minimap2 
 to map a masked reference 
 minimap2 -a Reference_barcodes.fasta subset_0-1-batch_S1_L001_R1_001.fasta > alignment.sam 
 
 # References 
 
 If you want to learn about alignment 
https://www.youtube.com/watch?v=XU8atPxM0VQ
