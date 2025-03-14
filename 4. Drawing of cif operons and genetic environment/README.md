# Genome annotation

We get the `.fasta` of the 786 genomes (786 $ech) on a repertory. Then, each genome were annotaded using Prokka (Seemann, 2014, <doi:10.1093/bioinformatics/btu153>) with the following command:

```
prokka $ech-genome.fasta --locustag $ech --prefix $ech --outdir Prokka-$ech --rfam --compliant --cpus 6
```


# Drawing of *cif* operons and genetic environment

## Drawing of CifA and CifB protein domains

First, an information matrix is a prerequisited for cifA and cifB. This matrix must include the total gene length, the start and end positions of each domain, and positions of ORF-disturbing mutations. The cifA and cifB matrices constructed for this study are attached to this GitHub repository in `cifA.txt` and `cifB.txt`. In `cifA.txt`, 'Apopt' is the abbreviation for the Apoptosis regulator-like domain and 'RNA' stands for the RNA-binding-like domain. Same for the `cifB.txt` where 'OTU' stands for the OTU-like cysteine protease, 'AAA' for the AAA-ATPases-like, 'IPPDE' for the PD-(D/E)XK nuclease N-terminal, 'IIPDDE' for the PD-(D/E)XK nuclease C-terminal, 'TPR' for the TPR repeats, 'DUB' for the Deubiquitinase DUB domain, 'TOXIN' for the Pore forming toxin TcdA/TcdB domain, 'DUF' for the DUF3491 domain, 'RTX' for the RTX toxin, 'SAL' for the salivary-gland toxin, 'ANK' for the Ankyrin repeat domains. and 'LATRO' for the Latrotoxin domain.

Then, these two matrices are imported into RStudio v4.3.1 (http://www.rstudio.com/):




for the attached R code
