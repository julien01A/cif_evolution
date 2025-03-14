# Genome annotation

We get the `.fasta` of the 786 genomes (407 $ech) on a repertory. Then, each genome were annotaded using Prokka (Seemann, 2014, doi:10.1093/bioinformatics/btu153) with the following command:

```
prokka $ech-genome.fasta --locustag $ech --prefix $ech --outdir Prokka-$ech --rfam --compliant --cpus 6
```
