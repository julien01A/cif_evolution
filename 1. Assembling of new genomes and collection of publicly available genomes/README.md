# Assembling of new genomes and collection of publicly available genomes

## New sequenced and assembled *Rickettsia* genomes

Organ metagenomes were sequenced for eight individual or pool of ticks. Eight Rickettsia genomes were de novo assembled (one by organ metagenome): (i) *Rickettsia lusitaniae* strain R-Oe isolated from *Ornithodoros erraticus*, (ii) *R. vini* strain IarboMS4 isolated from *Ixodes arboricola*, (iii) *R. raoultii* strains Dreti_100P and Dreti_100F isolated from *Dermacentor reticulatus*, and (iv) *Rickettsia sp.* strains AdisF19, AdisM1, AdisP2, and AdisP3Sn isolated from *Amblyomma dissimile*.

Three different sequencing and assembling methods were used depending on metagenomes. Below, we present the three methods:
Method 1. *Rickettsia lusitaniae* strain R-Oe
Method 2. *R. vini* strain IarboMS4, *R. raoultii* strains Dreti_100P and Dreti_100F, *Rickettsia sp.* strains AdisF19, AdisM1 and AdisP3Sn
Method 3. *Rickettsia sp.* strain AdisP2

### Method 1

**POUR MARIE**
Sequencing: Illumina HiSeq 2500 platform 
MEGAHIT (v1.2.9) -> assembling
CONCOCT (v1.1.0) -> binning
Anvi’o pipeline (v7.1) -> taxonomic assignment
anvi-refine tool with the Anvi’o interactive interface
BLAST of all contigs

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

## Plasmid identification
In *Rickettsia sp.* strains AdisF19 and AdisM1, plasmids were identified through visual inspection of the assembly graph using `Bandage (v0.8.1)` (<https://github.com/rrwick/Bandage>, doi: <https://10.1093/bioinformatics/btv383>). Contigs forming connected subgraphs with a total length compatible with known *Rickettsia* plasmids, were flagged as candidate plasmids. A `BLAST` search was performed against the NCBI database. Hits specific to *Rickettsia* plasmid sequences confirmed the plasmidic nature of the contigs. Based on distinct coverage values, lengths, and gene content, the plasmid candidates were retained as separate from chromosomal scaffolds.

## Genome visualization
The eight *Rickettsia* genome representations were performed using the `Proksee` online tool (<https://proksee.ca/>, doi: <https://doi.org/10.1093/nar/gkad326>)

## Genomes from public database

All the accession numbers of publicly available genomes analyzed in this study are listed in **Table S1** and retrievable on Genbank (<https://www.ncbi.nlm.nih.gov/datasets/genome/>), except few of them only available on publication supplementaries (<https://10.1093/sysbio/syv084>, <https://10.7717/peerj.5486>). `.fasta` files were retrieved for each genome and annotated with Prokka.
