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
