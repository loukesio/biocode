![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)
[![forthebadge](https://forthebadge.com/images/badges/works-on-my-machine.svg)](https://forthebadge.com)

# MINIMAP2 <br>
---

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



References: <br>
[1] https://github.com/lh3/minimap2/blob/master/cookbook.md#map-sr : Cookbook on minimap2 <br>
[2] https://forthebadge.com/ : How to add a badge <br>
[3] https://github.com/adam-p/markdown-here/wiki/Markdown-Here-Cheatsheet : Markdown Cheatseat <br>

# SAMTOOLS <br>
---

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

**Important** : Use biosyntaxt to beautify the outcome of your sam file https://github.com/bioSyntax/bioSyntax
samtools view -h NA19238.bam | sam-less -

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

## 7. Remove sequences that map to a specific area in the genome 

bedtools intersect -abam file.bam -b filter.bed -v > filtered.bam
filter.bed should contain

chr    start     end

and maybe 
samtools view input.bam -b -h -o output_inRegions.bam -U output_outRegions.bam -L Regions.bed



Other important references: <br>
[1] http://quinlanlab.org/tutorials/samtools/samtools.html  : Notes from Quinlan lab <br>
[2] https://davetang.org/wiki/tiki-index.php?page=SAMTools  : Notes from wiki cms <br>
[3] http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/1040-2/ : extract unmapped reads <br>
[4] https://www.metagenomics.wiki/tools/samtools/number-of-reads-in-bam-file : count unmapped reads, metagenomics notes <br>
[5] http://www.bioinf.uni-leipzig.de/publications/supplements/13-008 : reuse unmapped reads with segemehl <br>
[6] https://davetang.github.io/learning_bam_file/#filtering-unmapped-reads : Notes from dave tang <br>
[7] https://qnot.org/2012/04/14/counting-the-number-of-reads-in-a-bam-file/ : Count number of reads thats very interesting <br>
[8] https://comppopgenworkshop2019.readthedocs.io/en/latest/contents/02_bam_files/bam_files.html : Impressive Converting, Filtering SAM and BAM Files


# Comments: 

## Use taxonomic profile softwares in your raw data 
You can use 
* [Kraken2](https://github.com/DerrickWood/kraken2/wiki)
* [Kaiju](https://github.com/bioinformatics-centre/kaiju) <br>
and you can visualise everything using <br>
* [Kronatools](https://github.com/marbl/Krona)

### Kraken2 - Krona tools
1. Make a database 
```
kraken2-build --standard --threads 24 --db kraken_db
```
2. run your data over the db
```
kraken2 --paired --threads {threads} -db {params.db} --confidence 0.5 --output output/ cleandata/sample1_clean_R1.fq cleandata/sample1_clean_R2.fq ... cleandata/samplen_clean_R1.fq cleandata/samplen_clean_R2.fq
```
3. use kraken output to Krona 
```
ktImportTaxonomy -q 2 -t 3 Sample1.txt Sample2.txt -o krona.html 
```
References: <br>
[1] https://genomics.sschmeier.com/ngs-taxonomic-investigation/index.html : thats the best kraken turorial I ve found so far

### Merge files from different lanes 
If your files come into multiple lanes you can concatenate them using the following code: 

```
for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

    do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > "$i"_ME_L001_R1_001.fastq.gz

       echo "Merging R2"

cat "$i"_L00*_R2_001.fastq.gz > "$i"_ME_L001_R2_001.fastq.gz

done
```

Also check the following commands which might be very useful to you <br>
1. https://unix.stackexchange.com/questions/394479/concatenating-multiple-fastq-files <br>
2. https://unix.stackexchange.com/questions/436771/how-to-concatenate-rna-seq-files-generated-in-differnt-lanes <br>
3. https://unix.stackexchange.com/questions/615815/for-loop-to-catenate-files-with-two-variables <br>


### Copy-paste the contents of folder to another folder 

```
cp -a /source/. /dest/
```

The -a option is an improved recursive option, that preserve all file attributes, and also preserve symlinks.

The . at end of the source path is a specific cp syntax that allow to copy all files and folders, included hidden ones.



## Add your favourite websites that contain bioinfromatic software 

https://astrobiomike.github.io/unix/installing_tools#my-bioinformatics-tools-bit <br>
https://astrobiomike.github.io/unix/modifying_your_path <br>
http://maasha.github.io/biopieces/ <br>


## How to find tRNAs and rRNAs in your genome? 

+ To find tRNAs you can use <br>
    - A BLAT search in the USSC browser https://genome.ucsc.edu/cgi-bin/hgBlat <br>
    - tRNA scan is useful for sequences less than a 1Million bp http://trna.ucsc.edu/tRNAscan-SE/ <br>
    - tRNA finder seems crap https://ei4web.yz.yamagata-u.ac.jp/~kinouchi/tRNAfinder/

+ To find rRNAs you can use <br>
    - barrnap https://github.com/tseemann/barrnap
    - https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04316-z this paper contains a lot of info about rRNA prediction
    - BLAST using an rRNA database



