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
