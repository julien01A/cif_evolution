# The diversity and spread of cytoplasmic incompatibility genes among maternally inherited symbionts: *Wolbachia* do not walk alone - GitHub page

Amoros J, Buysse M, Floriano AM, Moumen B, Vavre F, Bouchon D, Duron D

This GitHub page shares the command lines and additional data used in this study. It is organized according to the typical workflow used to screen *cifA-cifB* genes from `.fasta` bacterial genome files retrieved from new assemblies or public databases, verify the presence of *cif* genes, determine their flanking genetic environment, visualize the results, and perform phylogenetic analyses to place each element in an evolutionary context.

The page is divided as follows: 

#### 1. Assembling of new genomes and collection of publicly available genomes
This section provides details about the methods used to assemble the eight *Rickettsia* genomes from metagenomes of tick organs (assembly, correction, and alignment command lines).

**Ajout de Marie**  
**Long-reads**  
*De novo* assembly of `XXX` dataset was performed from long-reads (MinION, Oxford Nanopore) using `Flye` (https://github.com/fenderglass/Flye, M.K., D.M.B., B.B., A.G., M.R., S.B.S., K.K., J.Y., E.P., T.P.L.S. and P.A.P., (2020) metaFlye: scalable long-read metagenome assembly using repeat graphs, Nature Methods. doi:s41592-020-00971-x), as follows:
```
module load bioinfo/Flye/2.4.1
flye --nano-raw $ech_porechopped_all.fq.gz --out-dir ./$ech-Flye --threads 12 --iterations 5 --meta --min-overlap 8000 --debug --genome-size 200000
```
The resulting file `assembly.fasta` is then corrected using only the long-reads' dataset using `Medaka` Oxford Nanopore tool (https://github.com/nanoporetech/medaka):
```
mkdir ./Flye-Medaka
samtools faidx assembly.fasta
module load bioinfo/medaka/1.5
medaka_consensus -i $ech_porechopped_all.fq.gz -d assembly.fasta -m r941_min_fast_g303 -o Flye-Medaka -t 4
```
The identification of the XXX MAG is based on the taxonomic assignation of the 16S rDNA sequence of a circular contig visualized with `Bandage` (https://rrwick.github.io/Bandage/, Wick R.R, Schultz M.B., Zobel J., Holt K.E. (2015) Bandage: Interactive visualization of de novo genome assemblies. Bioinformatics. doi: 10.1093/bioinformatics/btv383) (assignation based on the online NCBI BLAST tool). 
**After this step, the contig is considered as XXX MAG, designated as `XXX`.**

**Short-reads**  
## 1.1. Assembly
Reads were assembled using `MEGAHIT` (https://github.com/voutcn/megahit, Li D., Liu C-M., Luo R., Sadakane K., and Lam T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics. doi: 10.1093/bioinformatics/btv033): 
```
megahit -1 $ech-R1.fastq.gz -2 $ech-R2.fastq.gz --k-list 59,77,99 -t 6 -o $ech-metaMEGAHIT
```
## 1.2. Binning and retrieving *R. lusitaniae* MAGs
*R. lusitaniae* MAGs were retrieved from assemblies using `CONCOCT` (https://github.com/BinPro/CONCOCT, Alneberg J., Smári Bjarnason B., de Bruijn I., Schirmer M., Quick J., Ijaz U.Z., Lahti L., Loman N.J., Andersson A.F., & Quince C., (2014) Binning metagenomic contigs by coverage and composition, Nature Methods. doi: 10.1038/nmeth.3103) and the `anvi'o` pipeline (https://anvio.org/, Eren A.M., Kiefl E., Shaiber A. et al., (2021) Community-led, integrated, reproducible multi-omics with anvi’o. Nature Microbiology. doi: 10.1038/s41564-020-00834-3). 
First, the contigs were renamed to match the requirements of `anvi'o`. 
```
sed '/^>/s/ .*//' $ech-final.contigs.fa > $ech-final.contigs-rename.fa
rm $ech-final.contigs.fa
mv $ech-final.contigs-rename.fa $ech-final.contigs.fa
```
The contigs were binned using `CONCOCT`: 
```
bwa index $ech-final.contigs.fa
bwa mem -t 4 $ech-final.contigs.fa $ech-R1.fastq.gz $ech-R2.fastq.gz | samtools sort -@ 4 -T mapped -O BAM -o $ech-reads-mapped.bam
samtools index $ech-reads-mapped.bam

cut_up_fasta.py $ech-final.contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
concoct_coverage_table.py contigs_10K.bed $ech-reads-mapped.bam > coverage_table.tsv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b $ech-concoctMEGAHIT_output/ -t 4
merge_cutup_clustering.py $ech-concoctMEGAHIT_output/clustering_gt1000.csv > $ech-concoctMEGAHIT_output/clustering_merged.csv
mkdir $ech-concoctMEGAHIT_output/fasta_bins
extract_fasta_bins.py $ech-final.contigs.fa $ech-concoctMEGAHIT_output/clustering_merged.csv --output_path $ech-concoctMEGAHIT_output/fasta_bins
```
The `clustering_merged` CSV file produced by CONCOCT needed to be exported in a tabular-delimited TEXT file. 
The bins' names had to be renamed to match the requirements of `anvi'o` : 
```
awk '$2="bin"$2 {print}' bins-to-format.txt > bins-renamed.txt
```
The `bins-renamed.txt` file had to be transformed again to correspond to a tabular-delimited TEXT file, called `bins.txt` hereafter.  
The bins were taxonomically assigned using the `anvi'o` pipeline below: 
```
## To create the contigs database
anvi-gen-contigs-database -f $ech-final.contigs.fa -o $ech-metaMEGAHIT.db --ignore-internal-stop-codons -n Binning -T 4
anvi-run-hmms -c $ech-metaMEGAHIT.db

## To create the profile database
anvi-profile -i $ech-reads-mapped.bam -c $ech-metaMEGAHIT.db --min-contig-length 250 -T 4 -o $ech-PROFILE --cluster-contigs

## To import the bins as a collection
anvi-import-collection bins.txt -c $ech-metaMEGAHIT.db -p $ech-PROFILE/PROFILE.db -C bins --contigs-mode

## To assign the bins and visualize the results
anvi-run-scg-taxonomy -c $ech-metaMEGAHIT.db -T 2
anvi-estimate-scg-taxonomy -c $ech-metaMEGAHIT.db --output-file $ech-TAXONOMY.txt -p $ech-PROFILE/PROFILE.db -C bins --compute-scg-coverages -T 2
anvi-summarize -p $ech-PROFILE/PROFILE.db -c $ech-metaMEGAHIT.db -C bins -o $ech-SUMMARY 
```
Let's see the `$ech-SUMMARY` file (html format) to check the taxnonomy results and some stats about the bins. 
For Oerra sample: 
![Oerra-binning](https://user-images.githubusercontent.com/58982033/185303309-8845797a-9c10-4fc0-9d1f-adfdafc526ba.png)
Then, the `bin22` of Oerra sample was considered as raw *R. lusitaniae* MAGs. The bin was aligned using MAUVE (https://darlinglab.org/mauve/mauve.html, Darling A.C., Mau B., Blattner F.R., Perna N.T. (2004) Mauve: multiple alignment of conserved genomic sequence with rearrangements. Genome Research. doi: 10.1101/gr.2289704). Contigs that showed no homology with *Rickettsia*reference genomes were referenced. These contigs were then BLAST against the NCBI database and those that didn't show sequence similarity with any *Rickettsia* representative were removed, as follows:
```
FastaToTbl bin.fa | grep -wf contigsID_to_remove.txt | TblToFasta > bin_checked.fasta
```
**After this step, the bin was considered as final *R. lusitaniae* MAG.**  
**Fin ajout**  

#### 2. Genome annotation
This step analyzes gene boundaries and provides structural context for each genome using the Prokka annotation tool.

#### 3. Screening of *cifA* and *cifB* and their predicted protein domains
This part presents all the tools and command lines used to detect *cifA* and *cifB* genes and their predicted protein domains in bacterial genomes. The final *cifA* and *cifB* sequences, along with their predicted domains, used in this study are shared and can be used as queries to screen for *cif* presence in future genome analyses. 

#### 4. Determining *cif* flanking region
Here are the guidelines to identify *cif* flanking regions based on Prokka `.faa` protein files, RAGE and WO prophage databases, and BlastP analyses.

#### 5. Drawing of *cif* operons and genetic environment
This section proposes a detailed method to represent *cifA* and *cifB* structures and their flanking regions using `.gbk` files and homemade R scripts. Example files are provided to test the pipeline on a reduced dataset.

#### 6. Phylogenies
All the scripts used for phylogenetic analyses, from single-gene phylogenies to whole-genome phylogenies. Example files are provided to test the workflow on a reduced dataset.
