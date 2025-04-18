# Assembling of new genomes and collection of publicly available genomes

## New sequenced and assembled *Rickettsia* genomes

Organ metagenomes were sequenced for eight individual or pool of ticks. Eight Rickettsia genomes were de novo assembled (one by organ metagenome): (i) *Rickettsia lusitaniae* strain R-Oe isolated from *Ornithodoros erraticus*, (ii) *R. vini* strain IarboMS4 isolated from *Ixodes arboricola*, (iii) *R. raoultii* strains Dreti_100P and Dreti_100F isolated from *Dermacentor reticulatus*, and (iv) *Rickettsia sp.* strains AdisF19, AdisM1, AdisP2, and AdisP3Sn isolated from *Amblyomma dissimile*.

Three different sequencing and assembling methods were used depending on metagenomes. Below, we present the three methods:
#### Method 1. *Rickettsia lusitaniae* strain R-Oe
#### Method 2. *R. vini* strain IarboMS4, *R. raoultii* strains Dreti_100P and Dreti_100F, *Rickettsia sp.* strains AdisF19, AdisM1 and AdisP3Sn
#### Method 3. *Rickettsia sp.* strain AdisP2

### Method 1

The organ metagenome of *O. erraticus* was sequenced using Illumina short-read technology on a HiSeq 2500 platform.

Reads were assembled using `MEGAHIT (v1.2.9)` (<https://github.com/voutcn/megahit>, doi: <https://10.1093/bioinformatics/btv033>): 
```
megahit -1 $ech-R1.fastq.gz -2 $ech-R2.fastq.gz --k-list 59,77,99 -t 6 -o $ech-metaMEGAHIT
```
where `$ech` is the sequencing output name.

To retrieve the *R. lusitaniae* genome from the MEGAHIT meta-assembly contigs, we used `CONCOCT (v1.1.0)` (<https://github.com/BinPro/CONCOCT>, doi: <https://10.1038/nmeth.3103>) and the `anvi'o` pipeline (<https://anvio.org/>, doi: <https://10.1038/s41564-020-00834-3>) as follows: 
First, the contigs were renamed to match the requirements of `anvi'o`. 
```
sed '/^>/s/ .*//' $ech-final.contigs.fa > $ech-final.contigs-rename.fa
rm $ech-final.contigs.fa
mv $ech-final.contigs-rename.fa $ech-final.contigs.fa
```

Then, the contigs were binned using `CONCOCT (v1.1.0)`: 
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

The `clustering_merged` CSV file produced by `CONCOCT (v1.1.0)` needed to be exported in a tabular-delimited TEXT file. 
The bins' names had to be renamed to match the requirements of `anvi'o`: 
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
Let's see the `$ech-SUMMARY` file (html format) to check the taxnonomy results and some stats about the bins. Here is the `$ech-SUMMARY` file for the bin22 of *O.erraticus* which matched with one *Rickettsia* (next confirmed as *R. lusitaniae*): <https://user-images.githubusercontent.com/58982033/185303309-8845797a-9c10-4fc0-9d1f-adfdafc526ba.png>
The contigs constituing the `bin22` were aligned using MAUVE (<https://darlinglab.org/mauve/mauve.html>, doi: <https://10.1101/gr.2289704>). Contigs that showed no homology with *Rickettsia* reference genomes were referenced. These contigs were then `BLAST` (<https://blast.ncbi.nlm.nih.gov/Blast.cgi>) against the NCBI database and those that didn't show sequence similarity with any *Rickettsia* representative were removed, as follows:
```
FastaToTbl bin.fa | grep -wf contigsID_to_remove.txt | TblToFasta > bin_checked.fasta
```

After this step, the curated bin was considered the final *R. lusitaniae*, called strain R-Oe, and was used for downstream analyses in the rest of the study.  

### Method 2

**POUR MARIE**
Sequencing: Illumina HiSeq 2500 platform for *R. vini* strain IarboMS4, Illumina NovaSeq 6000 platform for *R. raoultii* strains Dreti_100P and Dreti_100F, *Rickettsia sp.* strains AdisF19, AdisM1 and AdisP3Sn
A modifief version of the original Blobology pipeline (<https://github.com/blaxterlab/blobology>, doi: <https://10.3389/fgene.2013.00237>) was used. See the details in the following GitHub page (<https://github.com/annamariafloriano/RickettsiellaComparative>). 
This pipeline includes:
SPAdes (doi: <https://doi.org/10.1089/cmb.2012.0021>) -> assembling
bowtie2 (doi: <https://doi.org/10.1038/nmeth.1923>)
samtools (doi: <https://doi.org/10.1093/bioinformatics/btp352>).
Bandage (v0.8.1) (<https://github.com/rrwick/Bandage>, doi: <https://10.1093/bioinformatics/btv383>)
QUAST (v4.6.3)and miComplete (v1.1.1, -hmms Bact105) -> quality check

### Method 3
**POUR MARIE**
Sequencing: Oxford Nanopore MinION device (R9.4.1 flow cell)
Flye (v2.9) 
Medaka (v1.2.2) 
BUSCO (v5.3.2) 

*De novo* assembly was performed from long-reads using `Flye` (<https://github.com/fenderglass/Flye>, doi:<https://10.1038/s41592-020-00971-x>), as follows:
```
module load bioinfo/Flye/2.4.1
flye --nano-raw $ech_porechopped_all.fq.gz --out-dir ./$ech-Flye --threads 12 --iterations 5 --meta --min-overlap 8000 --debug --genome-size 200000
```
The resulting file `assembly.fasta` was then corrected using `Medaka` Oxford Nanopore tool (<https://github.com/nanoporetech/medaka>):
```
mkdir ./Flye-Medaka
samtools faidx assembly.fasta
module load bioinfo/medaka/1.5
medaka_consensus -i $ech_porechopped_all.fq.gz -d assembly.fasta -m r941_min_fast_g303 -o Flye-Medaka -t 4
```
The identification of the *Rickettsia* AdisP2 genome is based on the taxonomic assignation of the 16S rDNA sequence of a circular contig visualized with `Bandage(v0.8.1)` (<https://github.com/rrwick/Bandage>, doi: <https://10.1093/bioinformatics/btv383>) (assignation based on the online NCBI BLAST tool). 
**After this step, the contig is considered as XXX MAG, designated as `XXX`.**

## Plasmid identification
In *Rickettsia sp.* strains AdisF19 and AdisM1, plasmids were identified through visual inspection of the assembly graph using `Bandage (v0.8.1)` (<https://github.com/rrwick/Bandage>, doi: <https://10.1093/bioinformatics/btv383>). Contigs forming connected subgraphs with a total length compatible with known *Rickettsia* plasmids, were flagged as candidate plasmids. A `BLAST` search was performed against the NCBI database. Hits specific to *Rickettsia* plasmid sequences confirmed the plasmidic nature of the contigs. Based on distinct coverage values, lengths, and gene content, the plasmid candidates were retained as separate from chromosomal scaffolds.

## Genome visualization
The eight *Rickettsia* genome representations were performed using the `Proksee` online tool (<https://proksee.ca/>, doi: <https://doi.org/10.1093/nar/gkad326>)

## Genomes from public database

All the accession numbers of publicly available genomes analyzed in this study are listed in **Table S1** and retrievable on Genbank (<https://www.ncbi.nlm.nih.gov/datasets/genome/>), except few of them only available on publication supplementaries (<https://10.1093/sysbio/syv084>, <https://10.7717/peerj.5486>). `.fasta` files were retrieved for each genome and annotated with Prokka.
