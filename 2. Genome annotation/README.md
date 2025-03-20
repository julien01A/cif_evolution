# Genome annotation

We get the `.fasta` of the 786 genomes (407 $ech) on a repertory. Then, each genome were annotaded using `Prokka v1.14.6` (Seemann, 2014, <doi:10.1093/bioinformatics/btu153>) with the following command:

```
prokka $ech-genome.fasta --locustag $ech --prefix $ech --outdir Prokka-$ech --rfam --compliant --cpus 6
```

For each genome, Prokka give the list of all the protein annotated in a `.faa` file, which are necessary for Orthofinder analysis (Protein screening against specific databases, see part.4, or whole genome phylogenies, see part.5). The `.gbk` file will also be used to detect the position of *cif* genes in the bacterial genome and trace their neighboring genetic environments.
